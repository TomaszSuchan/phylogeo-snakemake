#!/usr/bin/env Rscript

# Plot phylogenetic tree from IQTree output
# This script reads a tree file, plots it with ggtree, and saves as PDF and ggplot object

library(ggtree)
library(treeio)
library(ggplot2)

# Get parameters from Snakemake object
treefile <- snakemake@input[["treefile"]]
output_pdf <- snakemake@output[["pdf"]]
output_rds <- snakemake@output[["rds"]]
support_threshold <- snakemake@params[["support_threshold"]]

# Read the tree file
tree <- read.iqtree(treefile)

# Count number of tips to determine plot size
n_tips <- length(tree@phylo$tip.label)

# Calculate appropriate plot dimensions
# Base height: 0.3 inches per tip, minimum 6 inches, maximum 50 inches
plot_height <- max(6, min(50, n_tips * 0.3))
# Width: proportional to height with reasonable bounds
plot_width <- max(8, min(20, plot_height * 0.8))

# Create the base plot
p <- ggtree(tree, layout = "rectangular") +
  geom_tiplab(size = 2.5, hjust = -0.05) +
  theme_tree2()

# Add support values only if they exceed the threshold
# IQTree stores bootstrap support in node labels
if (!is.null(tree@data$UFboot)) {
  # Ultrafast bootstrap support
  p <- p + geom_nodelab(
    aes(
      label = ifelse(!is.na(UFboot) & UFboot >= support_threshold, UFboot, ""),
      subset = !isTip
    ),
    size = 2,
    hjust = -0.2,
    vjust = -0.5,
    color = "blue"
  )
} else if (!is.null(tree@data$label)) {
  # Standard bootstrap or other support values
  # Convert to numeric, handling any non-numeric values
  tree@data$support_numeric <- suppressWarnings(as.numeric(as.character(tree@data$label)))

  p <- p + geom_nodelab(
    aes(
      label = ifelse(!is.na(support_numeric) & support_numeric >= support_threshold,
                     support_numeric, ""),
      subset = !isTip
    ),
    size = 2,
    hjust = -0.2,
    vjust = -0.5,
    color = "blue"
  )
}

# Add scale bar
p <- p + geom_treescale(x = 0, y = 0, width = NULL, offset = 1)

# Adjust x-axis limits to accommodate tip labels
# Get the maximum branch length to adjust spacing
max_x <- max(p$data$x, na.rm = TRUE)
p <- p + xlim(NA, max_x * 1.3)  # Add 30% extra space for labels

# Save as PDF with calculated dimensions
ggsave(
  output_pdf,
  plot = p,
  width = plot_width,
  height = plot_height,
  units = "in",
  limitsize = FALSE
)

# Save as RDS (ggplot object) for further customization
saveRDS(p, output_rds)

# Print summary information
cat("Tree plotted successfully!\n")
cat(sprintf("  Number of tips: %d\n", n_tips))
cat(sprintf("  Plot dimensions: %.1f x %.1f inches\n", plot_width, plot_height))
cat(sprintf("  Support threshold: %.0f\n", support_threshold))
cat(sprintf("  Output PDF: %s\n", output_pdf))
cat(sprintf("  Output RDS: %s\n", output_rds))
