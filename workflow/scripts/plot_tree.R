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

compute_scale_width <- function(tree_phylo, fraction = 0.1) {
  root_to_tip <- node.depth.edgelength(tree_phylo)[seq_along(tree_phylo$tip.label)]
  total_depth <- max(root_to_tip, na.rm = TRUE)
  if (!is.finite(total_depth) || total_depth <= 0) {
    return(1)
  }

  target <- total_depth * fraction
  exponent <- floor(log10(target))
  candidates <- c(1, 2, 5) * (10 ^ exponent)
  candidates <- c(candidates, c(1, 2, 5) * (10 ^ (exponent - 1)), c(1, 2, 5) * (10 ^ (exponent + 1)))
  candidates <- sort(unique(candidates[candidates > 0]))
  candidates[which.min(abs(candidates - target))]
}

build_tree_plot <- function(tree_phylo, threshold, layout = "rectangular", scale_width = 1) {
  p <- ggtree(tree_phylo, layout = layout)

  if (layout == "daylight") {
    p <- p + geom_tiplab(size = 2.5) + theme_tree() + coord_equal()
  } else {
    p <- p + geom_tiplab(size = 2.5, hjust = -0.05) + theme_tree()
  }

  p <- add_support_labels(p, tree_phylo, threshold)
  if (layout == "daylight") {
    x_range <- range(p$data$x, na.rm = TRUE)
    y_range <- range(p$data$y, na.rm = TRUE)
    x_span <- diff(x_range)
    y_span <- diff(y_range)
    x_start <- x_range[1] + 0.06 * x_span
    y_pos <- y_range[1] + 0.06 * y_span
    label_y <- y_pos + 0.03 * y_span

    p +
      annotate(
        "segment",
        x = x_start, xend = x_start + scale_width,
        y = y_pos, yend = y_pos,
        linewidth = 0.5, color = "black"
      ) +
      annotate(
        "text",
        x = x_start + scale_width / 2,
        y = label_y,
        label = format(scale_width, scientific = FALSE, trim = TRUE),
        size = 3
      )
  } else {
    p + geom_treescale(x = 0, y = 0, width = scale_width)
  }
}

plot_dimensions <- function(tree_phylo, layout = "rectangular") {
  n_tips <- length(tree_phylo$tip.label)
  if (layout == "daylight") {
    # Unrooted daylight trees render best on square canvases.
    side <- max(8, min(20, 4 + sqrt(n_tips) * 0.8))
    plot_height <- side
    plot_width <- side
  } else {
    plot_height <- max(6, min(50, n_tips * 0.3))
    plot_width <- max(8, min(20, plot_height * 0.8))
  }
  list(n_tips = n_tips, width = plot_width, height = plot_height)
}

save_plot_outputs <- function(tree_phylo, pdf_path, rds_path, threshold, layout = "rectangular", scale_width = 1) {
  dims <- plot_dimensions(tree_phylo, layout = layout)
  p <- build_tree_plot(tree_phylo, threshold, layout = layout, scale_width = scale_width)

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
shared_scale_width <- compute_scale_width(rooted_tree, fraction = 0.1)

unrooted_result <- save_plot_outputs(
  base_tree,
  snakemake@output[["unrooted_pdf"]],
  snakemake@output[["unrooted_rds"]],
  support_threshold,
  layout = "daylight",
  scale_width = shared_scale_width
)
rooted_result <- save_plot_outputs(
  rooted_tree,
  snakemake@output[["rooted_pdf"]],
  snakemake@output[["rooted_rds"]],
  support_threshold,
  layout = "rectangular",
  scale_width = shared_scale_width
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
cat(sprintf("  Shared scale width: %s\n", format(shared_scale_width, scientific = FALSE, trim = TRUE)))
cat(sprintf("  Unrooted PDF: %s\n", snakemake@output[["unrooted_pdf"]]))
cat(sprintf("  Rooted PDF: %s\n", snakemake@output[["rooted_pdf"]]))
cat(sprintf("  Unrooted Newick: %s\n", snakemake@output[["unrooted_tree"]]))
cat(sprintf("  Rooted Newick: %s\n", snakemake@output[["rooted_tree"]]))
