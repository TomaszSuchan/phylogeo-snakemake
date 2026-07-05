#!/usr/bin/env Rscript
# Align ancestry clusters (Q-matrix columns) across K within a method, and across
# methods, so that a given biological cluster keeps the same column - and therefore
# the same colour - in every barplot/map. Downstream plotting colours columns purely
# by position (structure_colors[1:n_clusters]), so permuting columns here is all that
# is needed for consistent colours.
#
# Strategy
#   * The matching itself is done by clue::solve_LSAP (an optimal linear-sum
#     assignment / Hungarian solver). We only build the column-similarity cost
#     matrix and the K->K+1 loop around it - the algorithm is the package.
#   * One method is the reference (ADMIXTURE by default, else the first enabled
#     ancestry method). The reference is aligned across K by progressive chaining
#     (CLUMPAK-style): persistent clusters keep their slot, a splitting cluster gets
#     the next new slot appended.
#   * Every other method is matched, per K, to the reference at the same K. Because
#     the reference is itself consistent across K, this makes each method consistent
#     both across K and with the reference (hence with every other method).
#
# Reading/writing formats:
#   * Most methods write a headerless, tab/space-delimited numeric Q matrix.
#   * DAPC writes membership_probs.txt with a header and an `Ind` column.
#   * All aligned outputs are written as uniform headerless numeric Q matrices.

suppressPackageStartupMessages(library(clue))

## ---- IO helpers -------------------------------------------------------------

# Headerless numeric Q matrix (structure/snmf/tess3/alstructure/construct Qmatrix.txt,
# ADMIXTURE .Q, fastStructure .meanQ, and our own aligned outputs).
read_q_matrix <- function(path) {
  as.matrix(read.table(path, header = FALSE, sep = ""))
}

# DAPC membership_probs.txt: header row + an `Ind` column, remaining columns are the
# per-cluster posterior membership probabilities.
read_dapc_membership <- function(path) {
  df <- read.table(path, header = TRUE, sep = "\t", check.names = FALSE)
  ind_col <- which(tolower(colnames(df)) == "ind")
  if (length(ind_col) == 0) {
    stop(sprintf("DAPC membership file is missing an 'Ind' column: %s", path))
  }
  as.matrix(df[, -ind_col[1], drop = FALSE])
}

read_method_q <- function(path, method) {
  if (identical(method, "dapc")) read_dapc_membership(path) else read_q_matrix(path)
}

write_q_matrix <- function(mat, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(mat, file = path, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
}

## ---- Core alignment ---------------------------------------------------------

# Dissimilarity between the columns of A (individuals x a) and B (individuals x b),
# returned as an a-by-b matrix. Default metric is 1 - Pearson correlation, which is
# scale-free and works well on ancestry columns. Zero-variance columns (correlation
# undefined) fall back to 1 - shared mass so that identical constant columns still
# match. All entries are non-negative, as required by clue::solve_LSAP.
column_cost <- function(A, B) {
  n <- nrow(A)
  a <- ncol(A)
  b <- ncol(B)
  cost <- matrix(1, nrow = a, ncol = b)
  for (i in seq_len(a)) {
    for (j in seq_len(b)) {
      if (stats::sd(A[, i]) > 0 && stats::sd(B[, j]) > 0) {
        cost[i, j] <- 1 - stats::cor(A[, i], B[, j])
      } else {
        cost[i, j] <- 1 - sum(pmin(A[, i], B[, j])) / n
      }
    }
  }
  cost
}

# Permutation that orders the columns of `target` to correspond to the columns of
# `reference` (same number of clusters, same individuals in the same row order).
# target[, match_to_reference(target, reference)] lines up with reference columns 1..K.
match_to_reference <- function(target, reference) {
  stopifnot(ncol(target) == ncol(reference))
  cost <- column_cost(reference, target)   # reference rows x target columns
  as.integer(clue::solve_LSAP(cost))       # for reference slot r, the target column
}

# Progressive K->K+1 chaining for a method with no external reference. `q_by_k` is a
# list of matrices ordered by ascending K. Returns a list (same order) of column
# permutations: persistent clusters keep their slot, split clusters append new slots.
chain_across_k <- function(q_by_k) {
  n <- length(q_by_k)
  perms <- vector("list", n)
  # Smallest K: deterministic seed order (largest cluster first).
  perms[[1]] <- order(colSums(q_by_k[[1]]), decreasing = TRUE)
  if (n == 1) return(perms)
  prev_aligned <- q_by_k[[1]][, perms[[1]], drop = FALSE]
  for (i in 2:n) {
    cur <- q_by_k[[i]]
    # Existing slots (rows, fewer) matched to their best current column (cols, more).
    cost <- column_cost(prev_aligned, cur)
    slot_to_col <- as.integer(clue::solve_LSAP(cost))
    leftover <- setdiff(seq_len(ncol(cur)), slot_to_col)  # the split cluster(s)
    perm <- c(slot_to_col, leftover)                      # keep slots, append new
    perms[[i]] <- perm
    prev_aligned <- cur[, perm, drop = FALSE]
  }
  perms
}

# Column permutations for every K of a method. When reference matrices are supplied
# (keyed by K), each K is matched to the reference at the same K; otherwise the method
# is the reference and is aligned across K by chaining.
alignment_perms <- function(mats_by_k, ref_by_k = NULL) {
  k_names <- names(mats_by_k)
  perms <- vector("list", length(k_names))
  names(perms) <- k_names
  if (is.null(ref_by_k) || length(ref_by_k) == 0) {
    chained <- chain_across_k(mats_by_k)
    for (i in seq_along(k_names)) perms[[i]] <- chained[[i]]
  } else {
    for (kk in k_names) {
      ref <- ref_by_k[[kk]]
      if (!is.null(ref) && ncol(ref) == ncol(mats_by_k[[kk]])) {
        perms[[kk]] <- match_to_reference(mats_by_k[[kk]], ref)
      } else {
        perms[[kk]] <- seq_len(ncol(mats_by_k[[kk]]))  # no comparable reference: identity
      }
    }
  }
  perms
}

## ---- Snakemake entry point --------------------------------------------------

run_from_snakemake <- function(snakemake) {
  # Redirect output to the rule log, if one is configured.
  log_files <- tryCatch(as.character(unlist(snakemake@log)), error = function(e) character(0))
  if (length(log_files) > 0 && nzchar(log_files[1])) {
    dir.create(dirname(log_files[1]), recursive = TRUE, showWarnings = FALSE)
    log_con <- file(log_files[1], open = "wt")
    sink(log_con, type = "output")
    sink(log_con, type = "message")
    on.exit({
      while (sink.number(type = "message") > 0) sink(type = "message")
      while (sink.number(type = "output") > 0) sink(type = "output")
      close(log_con)
    }, add = TRUE)
  }

  method    <- snakemake@params[["method"]]
  k_values  <- as.character(unlist(snakemake@params[["k_values"]]))
  raw_files <- as.character(unlist(snakemake@input[["raw"]]))
  target_k  <- as.character(snakemake@params[["target_k"]])
  out_file  <- snakemake@output[["aligned"]]

  stopifnot(length(raw_files) == length(k_values))
  # Sort by ascending K so chaining sees increasing cluster counts.
  ord <- order(as.integer(k_values))
  k_values  <- k_values[ord]
  raw_files <- raw_files[ord]

  mats_by_k <- lapply(raw_files, read_method_q, method = method)
  names(mats_by_k) <- k_values

  ref_files <- tryCatch(as.character(unlist(snakemake@input[["reference"]])),
                        error = function(e) character(0))
  ref_kvals <- tryCatch(as.character(unlist(snakemake@params[["ref_k_values"]])),
                        error = function(e) character(0))
  ref_by_k <- NULL
  if (length(ref_files) > 0) {
    stopifnot(length(ref_files) == length(ref_kvals))
    ref_by_k <- lapply(ref_files, read_q_matrix)  # aligned reference is headerless numeric
    names(ref_by_k) <- ref_kvals
  }

  perms <- alignment_perms(mats_by_k, ref_by_k)
  perm <- perms[[target_k]]
  aligned <- mats_by_k[[target_k]][, perm, drop = FALSE]
  write_q_matrix(aligned, out_file)

  cat(sprintf("Aligned %s K=%s (%s reference) -> %s\n",
              method, target_k,
              if (is.null(ref_by_k)) "no" else "matched to", out_file))
  cat(sprintf("  column order applied: %s\n", paste(perm, collapse = ", ")))
  flush.console()
  quit(save = "no", status = 0)
}

# Only run when Snakemake injects the `snakemake` object; sourcing the file for tests
# just loads the functions above.
if (exists("snakemake")) {
  run_from_snakemake(snakemake)
}
