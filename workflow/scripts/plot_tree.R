#!/usr/bin/env Rscript

# Plot phylogenetic tree from IQ-TREE output.
# Produces both unrooted and rooted variants; rooted uses outgroup if provided,
# otherwise midpoint rooting.

library(ggtree)
library(treeio)
library(ggplot2)
library(ape)
library(phangorn)

# Get parameters from Snakemake object
treefile <- snakemake@input[["treefile"]]
support_threshold <- snakemake@params[["support_threshold"]]
outgroup_raw <- snakemake@params[["outgroup"]]

parse_outgroup <- function(outgroup_text) {
  if (is.null(outgroup_text) || outgroup_text == "") {
    return(character(0))
  }
  values <- trimws(strsplit(outgroup_text, ",", fixed = TRUE)[[1]])
  values[nzchar(values)]
}

add_support_labels <- function(plot_obj, tree_phylo, threshold) {
  if (is.null(tree_phylo$node.label)) {
    return(plot_obj)
  }

  support_numeric <- suppressWarnings(as.numeric(tree_phylo$node.label))
  node_ids <- seq.int(from = length(tree_phylo$tip.label) + 1,
                      length.out = tree_phylo$Nnode)
  support_df <- data.frame(node = node_ids, support = support_numeric)

  plot_obj %<+% support_df +
    geom_nodelab(
      aes(
        label = ifelse(!is.na(support) & support >= threshold, support, ""),
        subset = !isTip
      ),
      size = 2,
      hjust = -0.2,
      vjust = -0.5,
      color = "blue"
    )
}

build_tree_plot <- function(tree_phylo, threshold, layout = "rectangular") {
  p <- ggtree(tree_phylo, layout = layout)

  if (layout == "daylight") {
    p <- p + geom_tiplab(size = 2.5) + theme_tree()
  } else {
    p <- p + geom_tiplab(size = 2.5, hjust = -0.05) + theme_tree()
  }

  p <- add_support_labels(p, tree_phylo, threshold)
  p + geom_treescale(x = 0, y = 0)
}

plot_dimensions <- function(tree_phylo) {
  n_tips <- length(tree_phylo$tip.label)
  plot_height <- max(6, min(50, n_tips * 0.3))
  plot_width <- max(8, min(20, plot_height * 0.8))
  list(n_tips = n_tips, width = plot_width, height = plot_height)
}

save_plot_outputs <- function(tree_phylo, pdf_path, rds_path, threshold, layout = "rectangular") {
  dims <- plot_dimensions(tree_phylo)
  p <- build_tree_plot(tree_phylo, threshold, layout = layout)

  ggsave(
    pdf_path,
    plot = p,
    width = dims$width,
    height = dims$height,
    units = "in",
    limitsize = FALSE
  )
  saveRDS(p, rds_path)

  list(plot = p, n_tips = dims$n_tips, width = dims$width, height = dims$height)
}

# Read the tree file
tree <- read.iqtree(treefile)
base_tree <- unroot(tree@phylo)
outgroup_tips <- parse_outgroup(outgroup_raw)

valid_outgroup <- outgroup_tips[outgroup_tips %in% base_tree$tip.label]
invalid_outgroup <- setdiff(outgroup_tips, base_tree$tip.label)
if (length(invalid_outgroup) > 0) {
  warning(sprintf(
    "Outgroup taxa not found in tree and will be ignored: %s",
    paste(invalid_outgroup, collapse = ",")
  ))
}

if (length(valid_outgroup) > 0) {
  rooted_tree <- root(base_tree, outgroup = valid_outgroup, resolve.root = TRUE)
  rooting_method <- sprintf("outgroup (%s)", paste(valid_outgroup, collapse = ","))
} else {
  rooted_tree <- midpoint(base_tree)
  rooting_method <- "midpoint"
}

unrooted_result <- save_plot_outputs(
  base_tree,
  snakemake@output[["unrooted_pdf"]],
  snakemake@output[["unrooted_rds"]],
  support_threshold,
  layout = "daylight"
)
rooted_result <- save_plot_outputs(
  rooted_tree,
  snakemake@output[["rooted_pdf"]],
  snakemake@output[["rooted_rds"]],
  support_threshold,
  layout = "rectangular"
)

write.tree(base_tree, file = snakemake@output[["unrooted_tree"]])
write.tree(rooted_tree, file = snakemake@output[["rooted_tree"]])

# Print summary information
cat("Tree plotting complete.\n")
cat(sprintf("  Rooting method: %s\n", rooting_method))
cat(sprintf("  Number of tips: %d\n", unrooted_result$n_tips))
cat(sprintf("  Unrooted plot dimensions: %.1f x %.1f inches\n", unrooted_result$width, unrooted_result$height))
cat(sprintf("  Rooted plot dimensions: %.1f x %.1f inches\n", rooted_result$width, rooted_result$height))
cat(sprintf("  Support threshold: %.0f\n", support_threshold))
cat(sprintf("  Unrooted PDF: %s\n", snakemake@output[["unrooted_pdf"]]))
cat(sprintf("  Rooted PDF: %s\n", snakemake@output[["rooted_pdf"]]))
cat(sprintf("  Unrooted Newick: %s\n", snakemake@output[["unrooted_tree"]]))
cat(sprintf("  Rooted Newick: %s\n", snakemake@output[["rooted_tree"]]))
