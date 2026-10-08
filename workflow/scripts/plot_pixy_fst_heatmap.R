#!/usr/bin/env Rscript
# Create FST heatmap with dendrogram based on pixy FST results

library(tidyverse)
library(pheatmap)
library(gridExtra)
library(grid)

# Prevent creation of Rplots.pdf
pdf(NULL)

# Redirect all output to log file
log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

# Snakemake inputs/outputs
fst_summary_file <- snakemake@input[["fst_summary"]]
output_pdf <- snakemake@output[["pdf"]]
output_rds <- snakemake@output[["rds"]]
plot_width <- as.numeric(snakemake@params[["width"]])
plot_height <- as.numeric(snakemake@params[["height"]])
axis_title_size <- as.numeric(snakemake@params[["axis_title_size"]])
axis_text_size <- as.numeric(snakemake@params[["axis_text_size"]])
if (is.na(plot_width)) plot_width <- 12
if (is.na(plot_height)) plot_height <- 10
if (is.na(axis_title_size)) axis_title_size <- 10
if (is.na(axis_text_size)) axis_text_size <- 8

message("\n=== READING FST SUMMARY ===\n")
fst_df <- read.table(fst_summary_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
message(sprintf("Loaded %d FST comparisons\n", nrow(fst_df)))

# Create symmetric FST matrix
message("\n=== CREATING FST MATRIX ===\n")
populations <- unique(c(fst_df$pop1, fst_df$pop2))
n_pops <- length(populations)
message(sprintf("Found %d populations: %s\n", n_pops, paste(populations, collapse = ", ")))

# Check if we have at least 2 populations
if (n_pops < 2) {
  stop("FST heatmap requires at least 2 populations. Found only ", n_pops, " population(s).")
}

# Initialize matrix with zeros (diagonal)
fst_matrix <- matrix(0, nrow = n_pops, ncol = n_pops, dimnames = list(populations, populations))

# Fill matrix with FST values
for (i in 1:nrow(fst_df)) {
  pop1 <- fst_df$pop1[i]
  pop2 <- fst_df$pop2[i]
  fst_val <- fst_df$mean_fst[i]
  
  # Set both directions (symmetric)
  fst_matrix[pop1, pop2] <- fst_val
  fst_matrix[pop2, pop1] <- fst_val
}

message("FST matrix created: ", nrow(fst_matrix), " x ", ncol(fst_matrix), "\n")

# Negative FST estimates no differentiation, not a negative distance.
# Floor at 0 only for the tree; the heatmap still prints the raw estimates.
fst_dist <- fst_matrix
n_neg <- sum(fst_dist[lower.tri(fst_dist)] < 0, na.rm = TRUE)
fst_dist[fst_dist < 0] <- 0
message(sprintf("Floored %d negative pairwise FST values to 0 for clustering\n", n_neg))
dist_matrix <- as.dist(fst_dist)

# Hierarchical clustering
message("\n=== CALCULATING DENDROGRAM ===\n")
hc <- hclust(dist_matrix, method = "average")
pop_order <- hc$labels[hc$order]
message("Populations ordered by dendrogram: ", paste(pop_order, collapse = " -> "), "\n")

# Create heatmap with dendrogram using pheatmap.
# Pass the unordered matrix: pheatmap reorders by hc$order internally.
# Pre-reordering here would double-permute labels vs dendrogram tips.
#
# Color by the observed pairwise range (not from 0). Diagonal self-comparisons
# would otherwise pin the scale at 0 and compress contrast among real FST values.
message("\n=== CREATING HEATMAP ===\n")

pairwise <- fst_matrix[lower.tri(fst_matrix)]
z_min <- min(pairwise, na.rm = TRUE)
z_max <- max(pairwise, na.rm = TRUE)
if (!is.finite(z_min) || !is.finite(z_max)) {
  stop("No finite pairwise FST values available for coloring.")
}
if (isTRUE(all.equal(z_min, z_max))) {
  pad <- max(abs(z_min) * 0.01, 1e-8)
  z_min <- z_min - pad
  z_max <- z_max + pad
}
message(sprintf("Color scale from pairwise FST range: [%.6g, %.6g]\n", z_min, z_max))

fst_plot <- fst_matrix
diag(fst_plot) <- NA
color_breaks <- seq(z_min, z_max, length.out = 101)

p <- pheatmap::pheatmap(
  fst_plot,
  cluster_rows = hc,
  cluster_cols = hc,
  display_numbers = TRUE,
  number_format = "%.5f",
  fontsize = axis_title_size,
  fontsize_number = axis_text_size,
  color = colorRampPalette(c("white", "yellow", "orange", "red"))(100),
  breaks = color_breaks,
  na_col = "grey90",
  main = "",
  silent = FALSE
)

# Save plots
message("\n=== SAVING OUTPUT ===\n")
message(sprintf("Output PDF: %s\n", output_pdf))
message(sprintf("Output RDS: %s\n", output_rds))

dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)
pdf(output_pdf, width = plot_width, height = plot_height)
grid::grid.draw(p$gtable)
dev.off()

dir.create(dirname(output_rds), recursive = TRUE, showWarnings = FALSE)
saveRDS(p, output_rds)

message("\n=== COMPLETED SUCCESSFULLY ===\n")

# Close log file sinks
sink(type = "message")
sink(type = "output")
close(log_file)
