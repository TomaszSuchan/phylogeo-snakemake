#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(phangorn)
})

# Prevent accidental Rplots.pdf creation.
pdf(NULL)

# Redirect output and messages to Snakemake log.
log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

dist_file <- snakemake@input[["dist"]]
output_net <- snakemake@output[["net"]]

cat("Reading Euclidean distance matrix:", dist_file, "\n")
dist_df <- read.table(
  dist_file,
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (nrow(dist_df) < 3) {
  stop("NeighborNet requires at least 3 samples; found ", nrow(dist_df), ".")
}

if (!identical(rownames(dist_df), colnames(dist_df))) {
  stop("Distance matrix row names and column names must match and be in the same order.")
}

dist_mat <- as.matrix(dist_df)
diag(dist_mat) <- 0

if (any(!is.finite(dist_mat))) {
  stop("Distance matrix contains non-finite values.")
}

if (any(dist_mat < 0, na.rm = TRUE)) {
  stop("Distance matrix contains negative values.")
}

dist_obj <- as.dist(dist_mat)
net <- phangorn::neighborNet(dist_obj)

dir.create(dirname(output_net), recursive = TRUE, showWarnings = FALSE)
saveRDS(net, output_net)

cat("NeighborNet inferred successfully.\n")
cat("Samples:", nrow(dist_df), "\n")
cat("Output:", output_net, "\n")
