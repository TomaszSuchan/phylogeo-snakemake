#!/usr/bin/env Rscript
# Plot conStruct MAP log-posterior and layer contributions across K.
# Uses plot_choose_k_utils.R for Evanno-matched styling and config-driven size.

suppressPackageStartupMessages({
  library(ggplot2)
})

pdf(NULL)

script_dir <- tryCatch(
  dirname(normalizePath(snakemake@script)),
  error = function(e) "workflow/scripts"
)
source(file.path(script_dir, "plot_choose_k_utils.R"))
plot_dims <- read_choose_k_plot_dims(snakemake@params)

log_file <- snakemake@log[[1]]
dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")
on.exit({
  while (sink.number(type = "message") > 0) sink(type = "message")
  while (sink.number(type = "output") > 0) sink(type = "output")
  close(log_con)
}, add = TRUE)

lpd_summary <- read.table(
  snakemake@input[["lpd_summary"]],
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)
layer_summary <- read.table(
  snakemake@input[["layer_contribution_summary"]],
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

lpd_plot <- plot_choose_k_line(
  data = lpd_summary,
  x = "K",
  y = "median",
  ymin = "min",
  ymax = "max",
  ylab = "MAP log-posterior"
)

layer_summary$Layer <- factor(layer_summary$Layer, levels = sort(unique(layer_summary$Layer)))
layer_plot <- ggplot(layer_summary, aes(x = factor(K), y = Contribution, fill = Layer)) +
  geom_col(position = "stack", width = 0.7) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  xlab("K") +
  ylab("Layer contribution") +
  choose_k_plot_theme()

dir.create(dirname(snakemake@output[["lpd_pdf"]]), recursive = TRUE, showWarnings = FALSE)

choose_k_ggsave(
  filename = snakemake@output[["lpd_pdf"]],
  plot = lpd_plot,
  width = plot_dims$width,
  height = plot_dims$height,
  dpi = plot_dims$dpi
)
choose_k_ggsave(
  filename = snakemake@output[["layer_pdf"]],
  plot = layer_plot,
  width = plot_dims$width,
  height = plot_dims$height,
  dpi = plot_dims$dpi
)

saveRDS(lpd_plot, file = snakemake@output[["lpd_rds"]])
saveRDS(layer_plot, file = snakemake@output[["layer_rds"]])

cat("Saved conStruct MAP lpd plot:", snakemake@output[["lpd_pdf"]], "\n")
cat("Saved conStruct layer contribution plot:", snakemake@output[["layer_pdf"]], "\n")
