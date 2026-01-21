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
site_dividers <- as.logical(snakemake@params[["site_dividers"]])
divider_width <- as.numeric(snakemake@params[["divider_width"]])
site_order <- snakemake@params[["site_order"]]
flip_axis <- as.logical(snakemake@params[["flip_axis"]])
site_labels_angle <- as.numeric(snakemake@params[["site_labels_angle"]])

# Fixed parameters optimized for 1.2 inch height
site_ticks_size <- -0.05
# Adjust site_labels_y based on angle to prevent clipping
# For rotated labels (45 or 90 degrees), need more space
if (site_labels_angle == 90) {
  site_labels_y <- -0.5  # More space for perpendicular labels
} else if (site_labels_angle == 45) {
  site_labels_y <- -0.45  # More space for diagonal labels
} else {
  site_labels_y <- -0.35  # Default for horizontal labels
}
site_labels_size <- 2.2

# Read Q matrix (space or tab delimited, no headers)
qmatrix <- read.table(qmatrix_file, header = FALSE, sep = "")
n_individuals <- nrow(qmatrix)
n_clusters <- ncol(qmatrix)
cat(sprintf("Q matrix dimensions: %d individuals x %d clusters\n", n_individuals, n_clusters))

# Read popmap file, first column is individual, second population
cat(sprintf("Reading popmap file: %s\n", popmap_file))
tryCatch({
  popmap <- read.table(popmap_file, header = FALSE, sep = "\t", comment.char = "", fill = TRUE, blank.lines.skip = TRUE)
  colnames(popmap) <- c("Ind", "Site")
  cat(sprintf("Successfully read popmap: %d individuals with %d columns\n", nrow(popmap), ncol(popmap)))

  # Check for lines with missing data
  incomplete_lines <- which(is.na(popmap$Site) | popmap$Site == "")
  if (length(incomplete_lines) > 0) {
    cat(sprintf("WARNING: Found %d line(s) with missing Site column:\n", length(incomplete_lines)))
    for (i in incomplete_lines[1:min(5, length(incomplete_lines))]) {
      cat(sprintf("  Line %d: Individual='%s', Site='%s'\n", i, popmap$Ind[i], popmap$Site[i]))
    }
    if (length(incomplete_lines) > 5) {
      cat(sprintf("  ... and %d more lines\n", length(incomplete_lines) - 5))
    }
    stop(sprintf("ERROR: Popmap file has %d incomplete line(s). Please check the file format.", length(incomplete_lines)))
  }
}, error = function(e) {
  cat(sprintf("ERROR reading popmap file '%s'\n", popmap_file))
  cat(sprintf("Error message: %s\n", e$message))
  # Try to read line by line to find the problematic line
  lines <- readLines(popmap_file)
  cat(sprintf("Total lines in file: %d\n", length(lines)))
  for (i in seq_along(lines)) {
    fields <- strsplit(lines[i], "\t")[[1]]
    if (length(fields) != 2) {
      cat(sprintf("Line %d has %d field(s) instead of 2: '%s'\n", i, length(fields), lines[i]))
      if (i <= 195 && i >= 185) {
        cat(sprintf("  Context (lines %d-%d):\n", max(1, i-2), min(length(lines), i+2)))
        for (j in max(1, i-2):min(length(lines), i+2)) {
          cat(sprintf("    %d: '%s'\n", j, lines[j]))
        }
      }
      stop(sprintf("Popmap file has malformed line %d", i))
    }
  }
  stop(e$message)
})

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

# Get color palette from config - subset first n_clusters colors
cat(sprintf("\n=== SETTING UP COLOR PALETTE ===\n"))
structure_colors_full <- unlist(snakemake@params[["structure_colors"]])
cluster_cols_val <- structure_colors_full[1:n_clusters]
cat(sprintf("Using first %d colors from structure_colors palette: %s\n",
            n_clusters, paste(cluster_cols_val, collapse = ", ")))

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
  site_labels_size = site_labels_size,
  site_labels_angle = site_labels_angle
) +
  coord_cartesian(clip = "off") +  # Allow labels to extend beyond plot area
  theme(
    axis.title.y = element_text(size = 8, hjust = 1),
    axis.text.y = element_text(size = 5),
    # Remove ALL facet strip "boxes" around site names (handles different ggplot2 versions)
    # Based on mapmixture source: site labels use facets, boxes come from strip backgrounds
    strip.background = element_blank(),
    strip.background.x = element_blank(),
    strip.background.y = element_blank(),
    strip.text = element_text(),
    strip.text.x = element_text(),  # Explicitly set strip text (no background)
    strip.text.y = element_text(),  # Explicitly set strip text (no background)
    panel.background = element_blank(),
    plot.background = element_blank(),
    # Remove any annotation backgrounds
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.key = element_blank(),
    # Add extra margin at bottom for rotated labels to prevent clipping
    plot.margin = margin(
      t = 5, 
      r = 5, 
      b = ifelse(site_labels_angle >= 45, max(30, site_labels_size * 5), 5), 
      l = 5, 
      unit = "pt"
    )
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
