#!/usr/bin/env Rscript
# Line plot of model-comparison scores across K (Evanno-style ggplot aesthetics).
# Shared theme/size helpers live in plot_choose_k_utils.R so tess3r, ADMIXTURE,
# conStruct, and DAPC plots stay consistent with plot_evanno.R.

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

cv_summary <- read.table(
  snakemake@input[["cv_summary"]],
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)
crossvalid <- isTRUE(snakemake@params[["crossvalid"]])
crossentropy <- isTRUE(snakemake@params[["crossentropy"]])
ylab_override <- snakemake@params[["ylab"]]
plot_df <- data.frame(
  K = cv_summary$K,
  Value = cv_summary$median,
  Min = cv_summary$min,
  Max = cv_summary$max
)
ylab <- if (!is.null(ylab_override) && nzchar(as.character(ylab_override))) {
  as.character(ylab_override)
} else if (crossentropy) {
  if (crossvalid) "Cross-validation cross-entropy" else "Cross-entropy"
} else if (crossvalid) {
  "Cross-validation score"
} else {
  "RMSE"
}

cv_plot <- plot_choose_k_line(
  data = plot_df,
  x = "K",
  y = "Value",
  ymin = "Min",
  ymax = "Max",
  ylab = ylab
)

dir.create(dirname(snakemake@output[["pdf"]]), recursive = TRUE, showWarnings = FALSE)
choose_k_ggsave(
  filename = snakemake@output[["pdf"]],
  plot = cv_plot,
  width = plot_dims$width,
  height = plot_dims$height,
  dpi = plot_dims$dpi
)
saveRDS(cv_plot, file = snakemake@output[["rds"]])

cat("Saved choose-K score plot:", snakemake@output[["pdf"]], "\n")
