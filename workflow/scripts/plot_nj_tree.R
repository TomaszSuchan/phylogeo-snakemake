#!/usr/bin/env Rscript

# Plot unrooted neighbour-joining tree from rapidNJ Newick output.

library(ggtree)
library(ggplot2)
library(ape)

ggsave_utils <- tryCatch(
  file.path(dirname(normalizePath(snakemake@script)), "plot_ggsave_utils.R"),
  error = function(e) "workflow/scripts/plot_ggsave_utils.R"
)
if (file.exists(ggsave_utils)) {
  source(ggsave_utils)
} else {
  source("workflow/scripts/plot_ggsave_utils.R")
}

treefile <- snakemake@input[["treefile"]]
support_threshold <- snakemake@params[["support_threshold"]]

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

plot_dimensions <- function(tree_phylo) {
  n_tips <- length(tree_phylo$tip.label)
  side <- max(8, min(20, 4 + sqrt(n_tips) * 0.8))
  list(n_tips = n_tips, width = side, height = side)
}

tree <- unroot(read.tree(treefile))

p <- ggtree(tree, layout = "daylight") +
  geom_tiplab(size = 2.5) +
  theme_tree() +
  coord_equal(clip = "off") +
  theme(plot.margin = margin(20, 20, 20, 20))
p <- add_support_labels(p, tree, support_threshold)

dims <- plot_dimensions(tree)
ggsave_pdf(
  snakemake@output[["unrooted_pdf"]],
  plot = p,
  width = dims$width,
  height = dims$height,
  units = "in",
  limitsize = FALSE
)
saveRDS(p, snakemake@output[["unrooted_rds"]])
write.tree(tree, file = snakemake@output[["unrooted_tree"]])

cat("NJ tree plotting complete.\n")
cat(sprintf("  Number of tips: %d\n", dims$n_tips))
cat(sprintf("  Plot dimensions: %.1f x %.1f inches\n", dims$width, dims$height))
cat(sprintf("  Support threshold: %.0f\n", support_threshold))
cat(sprintf("  Unrooted PDF: %s\n", snakemake@output[["unrooted_pdf"]]))
cat(sprintf("  Unrooted Newick: %s\n", snakemake@output[["unrooted_tree"]]))
