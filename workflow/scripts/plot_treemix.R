#!/usr/bin/env Rscript

# Wrapper around TreeMix's original base-R plotting functions:
#   plot_tree(stem)
#   plot_resid(stem, pop_order)

pdf(NULL)

log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")
on.exit({
  while (sink.number(type = "message") > 0) sink(type = "message")
  while (sink.number(type = "output") > 0) sink(type = "output")
  close(log_file)
}, add = TRUE)

source("workflow/scripts/treemix_plotting_funcs.R")

stem <- sub("\\.treeout\\.gz$", "", snakemake@input[["treeout"]])
populations_file <- snakemake@input[["populations"]]
m_value <- as.character(snakemake@params[["m"]])
plot_width <- as.numeric(snakemake@params[["width"]])
plot_height <- as.numeric(snakemake@params[["height"]])
migration_arrow_length <- as.numeric(snakemake@params[["migration_arrow_length"]])

cat("Using original TreeMix plotting functions.\n")
cat("Stem:", stem, "\n")
cat("Migration edges:", m_value, "\n")

populations <- read.table(
  populations_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)
pop_order <- tempfile(pattern = "treemix_pop_order_", fileext = ".txt")
writeLines(populations$treemix_population, pop_order)

dir.create(dirname(snakemake@output[["graph_pdf"]]), recursive = TRUE, showWarnings = FALSE)

grDevices::pdf(snakemake@output[["graph_pdf"]], width = plot_width, height = plot_height)
tree_result <- plot_tree(stem, cex = 0.8, arrow = migration_arrow_length)
title(sprintf("TreeMix graph, m = %s", m_value))
grDevices::dev.off()

grDevices::pdf(snakemake@output[["residual_pdf"]], width = plot_width, height = plot_height)
residual_result <- plot_resid(stem = stem, pop_order = pop_order, cex = 0.8)
title(sprintf("TreeMix residuals, m = %s", m_value))
grDevices::dev.off()

saveRDS(
  list(
    stem = stem,
    m = m_value,
    pop_order = populations$treemix_population,
    plotting_function = "plot_tree",
    result = tree_result
  ),
  snakemake@output[["graph_rds"]]
)

saveRDS(
  list(
    stem = stem,
    m = m_value,
    pop_order = populations$treemix_population,
    plotting_function = "plot_resid",
    result = residual_result
  ),
  snakemake@output[["residual_rds"]]
)

cat("Saved graph plot:", snakemake@output[["graph_pdf"]], "\n")
cat("Saved residual plot:", snakemake@output[["residual_pdf"]], "\n")
