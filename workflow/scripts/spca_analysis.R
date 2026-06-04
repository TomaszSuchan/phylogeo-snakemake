#!/usr/bin/env Rscript
# Spatial PCA (sPCA) analysis with adegenet (2.1.x)

Sys.setenv(R_INSTALL_PKG = "")
options(device = function(file = NULL, ...) {
  if (is.null(file)) {
    pdf(NULL)
  } else {
    pdf(file, ...)
  }
})

log_file <- snakemake@log[[1]]
dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")
on.exit({
  while (sink.number(type = "message") > 0) sink(type = "message")
  while (sink.number(type = "output") > 0) sink(type = "output")
  close(log_con)
}, add = TRUE)

suppressPackageStartupMessages({
  library(adegenet)
  library(vcfR)
})

parse_optional_int <- function(x) {
  if (is.null(x) || length(x) == 0) return(NULL)
  if (is.list(x) && length(x) == 1) x <- x[[1]]
  chr <- as.character(x)
  if (is.na(chr) || chr == "" || toupper(chr) %in% c("NULL", "NONE", "NA")) return(NULL)
  as.integer(chr)
}

parse_optional_num <- function(x) {
  if (is.null(x) || length(x) == 0) return(NULL)
  if (is.list(x) && length(x) == 1) x <- x[[1]]
  chr <- as.character(x)
  if (is.na(chr) || chr == "" || toupper(chr) %in% c("NULL", "NONE", "NA")) return(NULL)
  as.numeric(chr)
}

parse_optional_bool <- function(x, default = FALSE) {
  if (is.null(x) || length(x) == 0) return(default)
  if (is.list(x) && length(x) == 1) x <- x[[1]]
  if (is.logical(x)) return(isTRUE(x))
  chr <- toupper(as.character(x))
  if (is.na(chr) || chr == "" || chr %in% c("NULL", "NONE", "NA")) return(default)
  chr %in% c("TRUE", "T", "1", "YES", "Y")
}

add_choosecn_param <- function(args, name, value) {
  if (!is.null(value)) args[[name]] <- value
  args
}

has_duplicate_xy <- function(xy, tol = .Machine$double.eps^0.5) {
  coords <- as.data.frame(xy)
  colnames(coords) <- c("Lon", "Lat")
  any(duplicated(coords) | duplicated(coords, fromLast = TRUE))
}

spread_duplicate_xy <- function(xy, radius, seed = NULL) {
  if (is.null(radius) || radius <= 0) {
    stop("jitter_sd must be > 0 when spreading duplicate coordinates.")
  }
  xy_out <- xy
  coords <- as.data.frame(xy_out)
  colnames(coords) <- c("Lon", "Lat")
  coord_key <- paste(coords$Lon, coords$Lat, sep = "|")
  groups <- split(seq_len(nrow(xy_out)), coord_key)
  if (!is.null(seed)) set.seed(seed)

  n_adjusted <- 0L
  for (idx in groups) {
    if (length(idx) == 1L) next
    center <- colMeans(xy_out[idx, , drop = FALSE])
    angles <- if (length(idx) == 2L) {
      c(0, pi)
    } else {
      seq(0, 2 * pi, length.out = length(idx) + 1L)[seq_len(length(idx))]
    }
    offsets <- cbind(
      radius * cos(angles),
      radius * sin(angles)
    )
    xy_out[idx, ] <- center + offsets
    n_adjusted <- n_adjusted + length(idx)
  }

  if (has_duplicate_xy(xy_out)) {
    stop(
      "Duplicate coordinates remain after spreading; increase parameters.spca.jitter_sd ",
      "or use graph type 5, 6, or 7."
    )
  }

  list(xy = xy_out, n_adjusted = n_adjusted)
}

spca_eigenvalues <- function(spca_obj) {
  if (!is.null(spca_obj$eig)) return(as.numeric(spca_obj$eig))
  if (!is.null(spca_obj$val)) return(as.numeric(spca_obj$val))
  stop("sPCA object has no $eig or $val component.")
}

write_randtest <- function(test_obj, path) {
  capture.output(print(test_obj), file = path)
}

cat("=== sPCA Analysis ===\n")
cat("adegenet version:", as.character(packageVersion("adegenet")), "\n")
cat("VCF:", snakemake@input[["vcf"]], "\n")
cat("indpopdata:", snakemake@input[["indpopdata"]], "\n\n")

nfposi <- as.integer(snakemake@params[["nfposi"]])
nfnega <- as.integer(snakemake@params[["nfnega"]])
cn_type <- as.integer(snakemake@params[["type"]])
cn_d1 <- parse_optional_num(snakemake@params[["d1"]])
cn_d2 <- parse_optional_num(snakemake@params[["d2"]])
cn_k <- parse_optional_int(snakemake@params[["k"]])
cn_a <- parse_optional_num(snakemake@params[["a"]])
cn_dmin <- parse_optional_num(snakemake@params[["dmin"]])
jitter_duplicate_coords <- parse_optional_bool(snakemake@params[["jitter_duplicate_coords"]])
jitter_sd <- parse_optional_num(snakemake@params[["jitter_sd"]])
jitter_seed <- parse_optional_int(snakemake@params[["jitter_seed"]])
n_pca <- parse_optional_int(snakemake@params[["n_pca"]])
nperm <- as.integer(snakemake@params[["nperm"]])
rtest_k <- as.integer(snakemake@params[["rtest_k"]])

cat("Parameters:\n")
cat("  nfposi:", nfposi, "\n")
cat("  nfnega:", nfnega, "\n")
cat("  type (connection network):", cn_type, "\n")
cat("  d1:", cn_d1, " d2:", cn_d2, " k:", cn_k, " a:", cn_a, " dmin:", cn_dmin, "\n")
cat("  jitter_duplicate_coords:", jitter_duplicate_coords, " jitter_sd:", jitter_sd, " jitter_seed:", jitter_seed, "\n")
cat("  n_pca:", if (is.null(n_pca)) "all loci (no PCA reduction)" else n_pca, "\n")
cat("  nperm (Monte Carlo tests):", nperm, "\n")
cat("  rtest_k:", rtest_k, "\n\n")

if (cn_type == 5L && (is.null(cn_d1) || is.null(cn_d2))) {
  stop("spca type 5 (distance band) requires both parameters.spca.d1 and .d2 (same units as Lon/Lat).")
}
if (cn_type == 6L && is.null(cn_k)) {
  stop("spca type 6 (k nearest neighbours) requires parameters.spca.k.")
}
if (cn_type == 7L && (is.null(cn_a) || is.null(cn_dmin))) {
  stop("spca type 7 (inverse distance) requires parameters.spca.a and .dmin.")
}

indpopdata <- read.table(
  snakemake@input[["indpopdata"]],
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)
required_cols <- c("Ind", "Lat", "Lon")
if (!all(required_cols %in% colnames(indpopdata))) {
  stop(
    "indpopdata must contain columns: Ind, Lat, Lon. Found: ",
    paste(colnames(indpopdata), collapse = ", ")
  )
}
if (any(is.na(indpopdata$Lat) | is.na(indpopdata$Lon))) {
  stop("Missing Lat/Lon coordinates in indpopdata.")
}

cat("Loading VCF...\n")
vcf <- read.vcfR(snakemake@input[["vcf"]], verbose = FALSE)
gen <- vcfR2genind(vcf, sep = "/")
cat("genind:", nInd(gen), "individuals,", nLoc(gen), "loci\n")

common <- intersect(indNames(gen), indpopdata$Ind)
if (length(common) < 3) {
  stop("sPCA requires at least 3 individuals with coordinates in both VCF and indpopdata.")
}
if (length(setdiff(indNames(gen), common)) > 0) {
  cat("WARNING: dropping VCF samples without indpopdata/coordinates\n")
}
if (length(setdiff(indpopdata$Ind, common)) > 0) {
  cat("WARNING: dropping indpopdata rows not in VCF\n")
}

gen <- gen[common, ]
indpopdata <- indpopdata[match(common, indpopdata$Ind), , drop = FALSE]
xy <- as.matrix(indpopdata[, c("Lon", "Lat")])
rownames(xy) <- indpopdata$Ind

dup_rows <- duplicated(as.data.frame(xy)) | duplicated(as.data.frame(xy), fromLast = TRUE)
if (any(dup_rows) && cn_type %in% 1:4) {
  if (!isTRUE(jitter_duplicate_coords)) {
    stop(
      "Duplicate Lat/Lon coordinates detected and spca type ", cn_type,
      " cannot use duplicate points. Set parameters.spca.jitter_duplicate_coords: true ",
      "or choose graph type 5, 6, or 7."
    )
  }
  if (is.null(jitter_sd) || jitter_sd <= 0) {
    stop("parameters.spca.jitter_sd must be > 0 when jitter_duplicate_coords is true.")
  }
  spread <- spread_duplicate_xy(xy, radius = jitter_sd, seed = jitter_seed)
  xy <- spread$xy
  cat(
    "Spread", spread$n_adjusted,
    "samples with duplicate Lon/Lat on a circle of radius", jitter_sd,
    "coordinate units.\n"
  )
}

cat("Matched samples:", length(common), "\n\n")

cat("Running spca()...\n")
spca_args <- list(
  obj = gen,
  xy = xy,
  nfposi = nfposi,
  nfnega = nfnega,
  type = cn_type,
  ask = FALSE,
  scannf = FALSE,
  plot.nb = FALSE
)
spca_args <- add_choosecn_param(spca_args, "d1", cn_d1)
spca_args <- add_choosecn_param(spca_args, "d2", cn_d2)
spca_args <- add_choosecn_param(spca_args, "k", cn_k)
spca_args <- add_choosecn_param(spca_args, "a", cn_a)
spca_args <- add_choosecn_param(spca_args, "dmin", cn_dmin)
if (!is.null(n_pca)) {
  if ("n.pca" %in% names(formals(spca))) {
    spca_args$n.pca <- n_pca
  } else if ("retain" %in% names(formals(spca))) {
    spca_args$retain <- n_pca
  } else {
    warning("Installed adegenet::spca does not accept n.pca/retain; ignoring n_pca.")
  }
}

spca_obj <- tryCatch(
  do.call(spca, spca_args),
  error = function(e) stop("spca() failed: ", conditionMessage(e))
)

gen_tab <- tab(gen, freq = FALSE)
if (!is.matrix(gen_tab)) gen_tab <- as.matrix(gen_tab)

nfnega_fit <- as.integer(spca_obj$nfnega)
eig_vals <- spca_eigenvalues(spca_obj)

dir.create(dirname(snakemake@output[["results_rds"]]), recursive = TRUE, showWarnings = FALSE)
saveRDS(
  list(
    spca = spca_obj,
    gen_tab = gen_tab,
    genind = gen,
    nfposi = nfposi,
    nfnega = nfnega_fit
  ),
  snakemake@output[["results_rds"]]
)

scores_global <- as.data.frame(spca_obj$ls)
scores_global$Ind <- rownames(scores_global)
scores_global <- scores_global[, c("Ind", setdiff(names(scores_global), "Ind")), drop = FALSE]
write.table(
  scores_global,
  snakemake@output[["scores_global"]],
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

scores_local <- as.data.frame(spca_obj$li)
scores_local$Ind <- rownames(scores_local)
scores_local <- scores_local[, c("Ind", setdiff(names(scores_local), "Ind")), drop = FALSE]
write.table(
  scores_local,
  snakemake@output[["scores_local"]],
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

eig_df <- data.frame(
  axis = seq_along(eig_vals),
  eigenvalue = eig_vals
)
write.table(
  eig_df,
  snakemake@output[["eigenvalues"]],
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

capture.output(summary(spca_obj), file = snakemake@output[["summary_txt"]])

cat("Running global/local Monte Carlo tests (tutorial: global.rtest / local.rtest)...\n")
if (is.null(spca_obj$lw)) {
  stop("sPCA object has no $lw component; cannot run global/local tests.")
}
global_test <- global.rtest(gen_tab, spca_obj$lw, k = rtest_k, nperm = nperm)
local_test <- local.rtest(gen_tab, spca_obj$lw, k = rtest_k, nperm = nperm)
write_randtest(global_test, snakemake@output[["global_rtest_txt"]])
write_randtest(local_test, snakemake@output[["local_rtest_txt"]])
saveRDS(
  list(global = global_test, local = local_test),
  snakemake@output[["rtest_rds"]]
)

cat("\nGlobal axes:", ncol(spca_obj$ls), "\n")
cat("Local axes:", ncol(spca_obj$li), "\n")
cat("Global test:\n")
print(global_test)
cat("Local test:\n")
print(local_test)
cat("=== sPCA Analysis Complete ===\n")

if (file.exists("Rplots.pdf")) {
  tryCatch(file.remove("Rplots.pdf"), error = function(e) invisible(NULL))
}
