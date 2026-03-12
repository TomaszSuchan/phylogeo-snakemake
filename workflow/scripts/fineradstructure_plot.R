#!/usr/bin/env Rscript
# Plot fineRADstructure results using the upstream heatmap workflow.
# Population labels are derived from indpopdata metadata instead of sample-name
# parsing, so the plot works with arbitrary individual IDs.

library(ape)
library(XML)

log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")
on.exit({
  sink(type = "message")
  sink(type = "output")
  close(log_file)
}, add = TRUE)

snakemake@source("FinestructureLibrary.R")

NameLessSummary <- function(x) {
  paste(x, collapse = ";")
}

safe_read_chunks <- function(path) {
  as.matrix(read.table(
    path,
    row.names = 1,
    header = TRUE,
    skip = 1,
    check.names = FALSE,
    comment.char = "",
    quote = ""
  ))
}

normalize_pop_list <- function(pop_list) {
  lapply(pop_list, function(x) {
    vals <- as.character(x)
    vals[nzchar(vals)]
  })
}

summarize_cluster_label <- function(cluster_inds, indpopdata, label_by, max_values = 3) {
  if (tolower(label_by) == "none") {
    return(sprintf("n=%d", length(cluster_inds)))
  }

  idx <- match(cluster_inds, indpopdata$Ind)
  if (any(is.na(idx))) {
    missing_inds <- cluster_inds[is.na(idx)]
    stop(
      "Individuals from fineRADstructure are missing in indpopdata: ",
      paste(head(missing_inds, 10), collapse = ", ")
    )
  }

  values <- indpopdata[[label_by]][idx]
  values <- trimws(as.character(values))
  values[is.na(values) | values == ""] <- "NA"

  counts <- sort(table(values), decreasing = TRUE)
  top_counts <- head(counts, max_values)
  parts <- paste0(names(top_counts), " (", as.integer(top_counts), ")")
  label <- paste(parts, collapse = ", ")
  if (length(counts) > max_values) {
    label <- paste0(label, ", ...")
  }
  label
}

make_unique_labels <- function(labels) {
  ave(
    seq_along(labels),
    labels,
    FUN = function(x) {
      if (length(x) == 1) {
        ""
      } else {
        paste0("_", LETTERS[seq_along(x)])
      }
    }
  ) -> suffix
  paste0(labels, suffix)
}

validate_inputs <- function(tree, dataraw, indpopdata, label_by) {
  if (!all(tree$tip.label %in% rownames(dataraw))) {
    missing_tips <- setdiff(tree$tip.label, rownames(dataraw))
    stop(
      "Tree tips are missing from chunks.out matrix: ",
      paste(head(missing_tips, 10), collapse = ", ")
    )
  }
  if (!"Ind" %in% names(indpopdata)) {
    stop("indpopdata must contain an 'Ind' column.")
  }
  if (tolower(label_by) != "none" && !label_by %in% names(indpopdata)) {
    stop(
      "Configured fineradstructure.plot.label_by column '", label_by,
      "' not found in indpopdata. Available columns: ",
      paste(names(indpopdata), collapse = ", ")
    )
  }
}

mcmc_tree_xml <- snakemake@input[["mcmcTree"]]
mcmc_xml <- snakemake@input[["mcmc"]]
chunks_file <- snakemake@input[["chunks"]]
indpopdata_file <- snakemake@input[["indpopdata"]]

simple_pdf <- snakemake@output[["simple_pdf"]]
popavg_pdf <- snakemake@output[["popavg_pdf"]]
labeled_pdf <- snakemake@output[["labeled_pdf"]]
output_rds <- snakemake@output[["rds"]]

label_by <- snakemake@params[["label_by"]]
max_indv <- as.numeric(snakemake@params[["max_indv"]])
max_pop <- as.numeric(snakemake@params[["max_pop"]])
max_label_values <- as.integer(snakemake@params[["max_label_values"]])

cat("=== fineRADstructure Plotting ===\n")
cat("MCMC Tree XML:", mcmc_tree_xml, "\n")
cat("MCMC XML:", mcmc_xml, "\n")
cat("Chunks file:", chunks_file, "\n")
cat("indpopdata:", indpopdata_file, "\n")
cat("Label populations by:", label_by, "\n")
cat("Output simple PDF:", simple_pdf, "\n")
cat("Output popavg PDF:", popavg_pdf, "\n")
cat("Output labeled PDF:", labeled_pdf, "\n")
cat("Output RDS:", output_rds, "\n")
cat("===============================\n\n")

tryCatch({
  cat("Reading fineRADstructure inputs...\n")
  dataraw <- safe_read_chunks(chunks_file)
  treexml <- xmlTreeParse(mcmc_tree_xml)
  mcmcxml <- xmlTreeParse(mcmc_xml)
  indpopdata <- read.table(
    indpopdata_file,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    comment.char = "",
    quote = "",
    stringsAsFactors = FALSE
  )

  cat("Chunk matrix dimensions:", nrow(dataraw), "x", ncol(dataraw), "\n")
  cat("indpopdata dimensions:", nrow(indpopdata), "x", ncol(indpopdata), "\n")
  cat("indpopdata columns:", paste(names(indpopdata), collapse = ", "), "\n")

  # This mirrors the upstream script's XML parsing and tree-to-dendrogram workflow.
  ttree <- extractTree(treexml)
  tdend <- myapetodend(ttree, factor = 1)
  validate_inputs(ttree, dataraw, indpopdata, label_by)

  fullorder <- labels(tdend)
  datamatrix <- dataraw[fullorder, fullorder, drop = FALSE]
  cat("Ordered matrix dimensions:", nrow(datamatrix), "x", ncol(datamatrix), "\n")

  mcmc_rows <- length(extractValue(mcmcxml, "Number", getNames = FALSE))
  cat("Read", mcmc_rows, "iterations from MCMC XML\n")

  map_state <- extractValue(treexml, "Pop")
  if (length(map_state) == 0) {
    stop("No Pop state found in mcmcTree XML.")
  }
  mapstatelist <- normalize_pop_list(popAsList(map_state[[1]]))
  names(mapstatelist) <- vapply(mapstatelist, NameLessSummary, character(1))
  cat("MAP state contains", length(mapstatelist), "clusters\n")

  popdend <- makemydend(tdend, mapstatelist, summary = "NameLessSummary")
  cluster_keys_in_order <- labels(popdend)
  ordered_clusters <- lapply(cluster_keys_in_order, function(x) strsplit(x, ";", fixed = TRUE)[[1]])

  cluster_labels <- vapply(
    ordered_clusters,
    summarize_cluster_label,
    character(1),
    indpopdata = indpopdata,
    label_by = label_by,
    max_values = max_label_values
  )
  cluster_labels <- make_unique_labels(cluster_labels)

  cat("Cluster labels:\n")
  for (i in seq_along(cluster_labels)) {
    cat("  ", i, ":", cluster_labels[i], "\n")
  }

  some.colors_end <- MakeColorYRP(final = c(0.2, 0.2, 0.2))

  tmpmat <- datamatrix
  tmpmat[tmpmat > max_indv] <- max_indv
  pdf(file = simple_pdf, height = 25, width = 25)
  plotFinestructure(
    tmpmat,
    dimnames(tmpmat)[[1]],
    dend = tdend,
    cols = some.colors_end,
    cex.axis = 0.45,
    edgePar = list(p.lwd = 0, t.srt = 90, t.off = -0.1, t.cex = 1.2),
    main = "Simple coancestry"
  )
  dev.off()
  cat("Saved simple coancestry heatmap\n")

  popmeanmatrix <- getPopMeanMatrix(datamatrix, mapstatelist)
  tmpmat <- popmeanmatrix
  tmpmat[tmpmat > max_pop] <- max_pop
  pdf(file = popavg_pdf, height = 20, width = 20)
  plotFinestructure(
    tmpmat,
    dimnames(tmpmat)[[1]],
    dend = tdend,
    cols = some.colors_end,
    cex.axis = 0.45,
    edgePar = list(p.lwd = 0, t.srt = 90, t.off = -0.1, t.cex = 1.2),
    main = "Population-averaged coancestry"
  )
  dev.off()
  cat("Saved population-averaged coancestry heatmap\n")

  mappopsizes <- vapply(ordered_clusters, length, integer(1))
  labellocs <- PopCenters(mappopsizes)
  pdf(file = labeled_pdf, height = 25, width = 25)
  plotFinestructure(
    tmpmat,
    labelsx = cluster_labels,
    labelsatx = labellocs,
    labelsaty = labellocs,
    dend = tdend,
    cols = some.colors_end,
    xcrt = 90,
    ycrt = 0,
    cex.axis = 0.8,
    edgePar = list(p.lwd = 0, t.srt = 90, t.off = -0.1, t.cex = 1.2),
    hmmar = c(3, 0, 0, 1),
    main = paste("Population labels by", label_by)
  )
  dev.off()
  cat("Saved labeled population-averaged coancestry heatmap\n")

  saveRDS(
    list(
      tree = ttree,
      tree_dendrogram = tdend,
      population_dendrogram = popdend,
      individual_order = fullorder,
      map_clusters = ordered_clusters,
      cluster_labels = cluster_labels,
      label_by = label_by,
      chunk_matrix = datamatrix,
      popmean_matrix = popmeanmatrix
    ),
    output_rds
  )
  cat("Saved plot objects to RDS\n")

  cat("\n=== fineRADstructure Plotting Complete ===\n")
  cat("Simple PDF:", simple_pdf, "\n")
  cat("Population-average PDF:", popavg_pdf, "\n")
  cat("Labeled PDF:", labeled_pdf, "\n")
  cat("RDS:", output_rds, "\n")
  cat("=========================================\n")
}, error = function(e) {
  cat("ERROR in fineRADstructure plotting:", conditionMessage(e), "\n")
  stop("fineRADstructure plotting failed: ", conditionMessage(e))
})

