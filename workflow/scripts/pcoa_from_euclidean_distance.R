#!/usr/bin/env Rscript

# PCoA on Euclidean genetic distance matrix produced by calculate_euclidean_distance.R
# Output format mirrors PCA eigvec/eigval files so it can be plotted with plot_pca_single.R:
#   - First row: header line that will be skipped (skip = 1)
#   - First column: sample names (also used as FID/IID)
#   - Second column: duplicate of sample names (IID), then PC axes (PC1, PC2, ...)

# Basic safety check
if (!exists("snakemake")) {
  stop("This script is intended to be run via Snakemake (snakemake@input / snakemake@output).")
}

dist_file <- snakemake@input[["dist"]]
eigvecs_out <- snakemake@output[["eigvecs"]]
eigvals_out <- snakemake@output[["eigvals"]]

message("Reading Euclidean distance matrix from: ", dist_file)

# Read distance matrix: first column are row names (sample IDs), remaining columns are distances
dist_mat <- as.matrix(
  read.table(
    dist_file,
    header = TRUE,
    row.names = 1,
    check.names = FALSE,
    sep = "\t"
  )
)

# Ensure matrix is numeric
storage.mode(dist_mat) <- "double"

message("Distance matrix dimensions: ", paste(dim(dist_mat), collapse = " x "))

# Convert to dist object and run classical multidimensional scaling (PCoA)
message("Running cmdscale (PCoA)...")
dist_obj <- as.dist(dist_mat)

# Number of axes = min(n-1, number of columns)
n_ind <- nrow(dist_mat)
k_max <- n_ind - 1

pcoa_res <- cmdscale(dist_obj, eig = TRUE, k = k_max)

coords <- as.data.frame(pcoa_res$points)

# Name axes as PC1, PC2, ...
pc_names <- paste0("PC", seq_len(ncol(coords)))
colnames(coords) <- pc_names

samples <- rownames(coords)

# Build eigenvector-like table:
# First column: FID (use sample IDs)
# Second column: IID (use sample IDs)
eigvecs_df <- data.frame(
  FID = samples,
  IID = samples,
  coords,
  check.names = FALSE
)

message("Writing eigenvectors to: ", eigvecs_out)

# Header line that will be skipped by plot_pca_single.R (skip = 1)
header_line <- paste(c("#FID", "IID", pc_names), collapse = "\t")

dir.create(dirname(eigvecs_out), recursive = TRUE, showWarnings = FALSE)

con <- file(eigvecs_out, open = "wt")
writeLines(header_line, con = con)
write.table(
  eigvecs_df,
  file = con,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)
close(con)

message("Writing eigenvalues to: ", eigvals_out)

dir.create(dirname(eigvals_out), recursive = TRUE, showWarnings = FALSE)

eigvals <- pcoa_res$eig
write(eigvals, file = eigvals_out, ncolumns = 1, sep = "\n")

message("PCoA completed successfully.")

