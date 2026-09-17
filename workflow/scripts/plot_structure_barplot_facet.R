#!/usr/bin/env Rscript
# Stack per-K barplots (K >= 2) in one column with side K labels.

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
  library(gtable)
})

pdf(NULL)

script_dir <- tryCatch(
  dirname(normalizePath(snakemake@script)),
  error = function(e) "workflow/scripts"
)
source(file.path(script_dir, "plot_ggsave_utils.R"))
source(file.path(script_dir, "plot_structure_barplot_facet_utils.R"))

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

params <- snakemake@params
if (length(params) == 1L && is.list(params[[1L]]) && length(names(params)) == 0L) {
  params <- params[[1L]]
}

output_pdf <- snakemake@output[["pdf"]]
output_rds <- snakemake@output[["rds"]]
method_label <- params[["method_label"]]
flip_axis <- isTRUE(as.logical(params[["flip_axis"]]))
label_width <- as.numeric(params[["label_width"]])
panel_gap_pt <- as.numeric(params[["panel_gap"]])
legend_pad_in <- as.numeric(params[["legend_pad"]])
if (!is.finite(legend_pad_in) || legend_pad_in < 0) {
  legend_pad_in <- 0.4
}

barplot_input_names <- grep("^barplot_k", names(snakemake@input), value = TRUE)
if (length(barplot_input_names) == 0L) {
  stop("No per-K barplot RDS inputs found (expected names like barplot_k2, barplot_k3, ...).")
}

extract_k <- function(name) as.integer(sub("^barplot_k", "", name))
ord <- order(vapply(barplot_input_names, extract_k, integer(1)))
barplot_input_names <- barplot_input_names[ord]
barplot_ks <- vapply(barplot_input_names, extract_k, integer(1))

message(sprintf(
  "Building %s barplot facet from K = %s",
  method_label,
  paste(barplot_ks, collapse = ", ")
))

barplot_paths <- snakemake@input[barplot_input_names]
panels_raw <- lapply(barplot_paths, function(path) {
  p <- readRDS(path)
  if (!inherits(p, "ggplot")) {
    stop("Expected ggplot in ", path)
  }
  p
})

facet <- build_barplot_facet_gtable(
  panels_raw = panels_raw,
  barplot_ks = barplot_ks,
  flip_axis = flip_axis,
  label_width_in = label_width,
  panel_gap_pt = panel_gap_pt,
  legend_pad_in = legend_pad_in
)

message(sprintf(
  "Saving barplot facet PDF: %s (%.2f x %.2f in)",
  output_pdf,
  facet$width,
  facet$height
))
ggsave_pdf(output_pdf, facet$gt, width = facet$width, height = facet$height)

message(sprintf("Saving barplot facet RDS: %s", output_rds))
saveRDS(facet$gt, file = output_rds)

message("Barplot facet plot completed successfully.")
