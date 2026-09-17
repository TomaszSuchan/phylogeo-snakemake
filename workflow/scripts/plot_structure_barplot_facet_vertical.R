#!/usr/bin/env Rscript
# Vertical multi-K ancestry barplot facet: per-K panels tiled left-to-right with
# the lowest K on the left and the highest K on the right.
#
# Implementation: build the same single-column facet as the column layout but
# with K in DESCENDING order (highest K at the top, lowest at the bottom), then
# rotate the whole gtable 90 degrees clockwise. That maps the top row to the
# right and the bottom row to the left, so K ends up ascending left-to-right,
# the shared site-label strip lands on the left, and individuals run top-to-
# bottom. Reusing the (well-tested) column builder keeps the two facets
# consistent and handles multi-site dividers/labels correctly. There is no
# legend (the column layout does not draw one either).

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
label_width <- as.numeric(params[["label_width"]])
panel_gap_pt <- as.numeric(params[["panel_gap"]])

barplot_input_names <- grep("^barplot_k", names(snakemake@input), value = TRUE)
if (length(barplot_input_names) == 0L) {
  stop("No per-K barplot RDS inputs found (expected names like barplot_k2, barplot_k3, ...).")
}

extract_k <- function(name) as.integer(sub("^barplot_k", "", name))
# Descending K order: after the 90-degree clockwise rotation this becomes
# ascending left-to-right (lowest K on the left, highest K on the right).
ord <- order(vapply(barplot_input_names, extract_k, integer(1)), decreasing = TRUE)
barplot_input_names <- barplot_input_names[ord]
barplot_ks <- vapply(barplot_input_names, extract_k, integer(1))

message(sprintf(
  "Building vertical %s barplot facet from K = %s (rendered top-to-bottom, rotated to left-to-right)",
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

# Build the column-layout gtable (no legend padding; the vertical facet has no
# legend), then rotate it 90 degrees clockwise.
facet <- build_barplot_facet_gtable(
  panels_raw = panels_raw,
  barplot_ks = barplot_ks,
  flip_axis = FALSE,
  label_width_in = label_width,
  panel_gap_pt = panel_gap_pt,
  legend_pad_in = 0
)

rotated <- rotate_grob_cw(facet$gt, facet$width, facet$height)

message(sprintf(
  "Saving vertical barplot facet PDF: %s (%.2f x %.2f in)",
  output_pdf,
  rotated$width,
  rotated$height
))
ggsave_pdf(output_pdf, rotated$grob, width = rotated$width, height = rotated$height)

message(sprintf("Saving vertical barplot facet RDS: %s", output_rds))
saveRDS(rotated$grob, file = output_rds)

message("Vertical barplot facet plot completed successfully.")
