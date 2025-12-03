#!/usr/bin/env Rscript
# Script to create structure barplots using mapmixture

library(tidyverse)
library(mapmixture)

# Prevent creation of Rplots.pdf
pdf(NULL)

# Snakemake inputs/outputs
qmatrix_file <- snakemake@input[["qmatrix"]]
popmap_file <- snakemake@input[["popmap"]]
output_barplot <- snakemake@output[["barplot"]]
output_barplot_rds <- snakemake@output[["barplot_rds"]]
output_prefix <- snakemake@params[["output_prefix"]]

# Plot dimension parameters
width <- as.numeric(snakemake@params[["width"]])
height <- as.numeric(snakemake@params[["height"]])
dpi <- as.numeric(snakemake@params[["dpi"]])

# Barplot-specific parameters
cluster_cols <- snakemake@params[["cluster_cols"]]
site_dividers <- as.logical(snakemake@params[["site_dividers"]])
divider_width <- as.numeric(snakemake@params[["divider_width"]])
site_order <- snakemake@params[["site_order"]]
flip_axis <- as.logical(snakemake@params[["flip_axis"]])

# Fixed parameters optimized for 1.2 inch height
site_ticks_size <- -0.05
site_labels_y <- -0.35
site_labels_size <- 2.2

# Read Q matrix (space or tab delimited, no headers)
qmatrix <- read.table(qmatrix_file, header = FALSE, sep = "")
n_individuals <- nrow(qmatrix)
n_clusters <- ncol(qmatrix)
cat(sprintf("Q matrix dimensions: %d individuals x %d clusters\n", n_individuals, n_clusters))

# Read popmap file, first column is individual, second population
popmap <- read.table(popmap_file, header = FALSE, sep = "\t")
colnames(popmap) <- c("Ind", "Site")
cat(sprintf("Original popmap contains %d individuals with %d columns\n", nrow(popmap), ncol(popmap)))

# Remove duplicate individuals (keep first occurrence)
popmap <- popmap %>%
  distinct(Ind, .keep_all = TRUE)

cat(sprintf("After removing duplicates, popmap contains %d individuals\n", nrow(popmap)))

# Filter popmap to match the number of individuals in Q matrix
if (nrow(popmap) > nrow(qmatrix)) {
  cat(sprintf("Filtering popmap from %d to %d individuals to match Q matrix\n", nrow(popmap), nrow(qmatrix)))
  popmap <- popmap[1:nrow(qmatrix), ]
}

# Check dimensions match
if (nrow(qmatrix) != nrow(popmap)) {
  stop(sprintf("BŁĄD: Q matrix has %d individuals but popmap has %d individuals",
               nrow(qmatrix), nrow(popmap)))
}

# Combine Q matrix with individual-level data, rename colnames
qmatrix_with_data <- cbind(popmap$Site, popmap$Ind, qmatrix)
colnames(qmatrix_with_data)[1] <- "Site"
colnames(qmatrix_with_data)[2] <- "Ind"
colnames(qmatrix_with_data) <- gsub("V", "Cluster", colnames(qmatrix_with_data))

# For barplots, we don't need popdata (geographic coordinates)
# We only use the site information from the popmap file

# Parse cluster_cols parameter
# Snakemake automatically converts Python lists to R vectors
if (is.null(cluster_cols)) {
  cluster_cols_val <- NULL
} else if (length(cluster_cols) == 1 && (cluster_cols == "NULL" || cluster_cols == "")) {
  cluster_cols_val <- NULL
} else if (is.vector(cluster_cols) || is.list(cluster_cols)) {
  # Already a vector/list from Python list in config
  cluster_cols_val <- unlist(cluster_cols)
} else {
  # Try to parse as R expression (backward compatibility)
  cluster_cols_val <- tryCatch(
    eval(parse(text = cluster_cols)),
    error = function(e) NULL
  )
}

# Parse site_order parameter
# Snakemake automatically converts Python lists to R vectors
if (is.null(site_order)) {
  site_order_val <- NULL
} else if (length(site_order) == 1 && (site_order == "NULL" || site_order == "")) {
  site_order_val <- NULL
} else if (is.vector(site_order) || is.list(site_order)) {
  # Already a vector/list from Python list in config
  site_order_val <- unlist(site_order)
} else {
  # Try to parse as R expression (backward compatibility)
  site_order_val <- tryCatch(
    eval(parse(text = site_order)),
    error = function(e) NULL
  )
}

# Create traditional structure barplot
structure_barplot <- structure_plot(
  admixture_df = qmatrix_with_data,
  type = "structure",
  cluster_cols = cluster_cols_val,
  site_dividers = site_dividers,
  divider_width = divider_width,
  site_order = site_order_val,
  labels = "site",
  flip_axis = flip_axis,
  site_ticks_size = site_ticks_size,
  site_labels_y = site_labels_y,
  site_labels_size = site_labels_size
) +
  theme(
    axis.title.y = element_text(size = 8, hjust = 1),
    axis.text.y = element_text(size = 5),
  )

# Calculate plot dimensions
n_individuals <- nrow(qmatrix_with_data)
cat(sprintf("Number of individuals: %d\n", n_individuals))

# Fixed height, width scales with number of individuals
height <- 1.2  # fixed at 1.2 inches
width <- max(8, 1 + 0.01 * n_individuals)  # at least 8 inches, scaling with individuals
cat(sprintf("Plot dimensions: %.2f x %.2f inches\n", width, height))

# Save barplot
ggsave(
  filename = output_barplot,
  plot = structure_barplot,
  width = width,
  height = height,
  dpi = dpi,
  device = "pdf"
)

# Save barplot ggplot object as RDS
saveRDS(structure_barplot, file = output_barplot_rds)

cat("Successfully created structure barplot\n")
