#!/usr/bin/env Rscript
# Interpolated ancestry map using tess3r::plot.tess3Q (tess3r vignette workflow).

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

if (!requireNamespace("tess3r", quietly = TRUE)) {
  stop("tess3r is not installed in this rule environment.")
}
suppressPackageStartupMessages({
  library(tess3r)
})

results_rds <- snakemake@input[["results_rds"]]
output_pdf <- snakemake@output[["pdf"]]

k <- as.integer(snakemake@params[["k"]])
map_method <- as.character(snakemake@params[["map_method"]])
map_resolution <- as.integer(unlist(snakemake@params[["map_resolution"]]))
interpolation_knots <- as.integer(snakemake@params[["interpolation_knots"]])
structure_colors <- unlist(snakemake@params[["structure_colors"]])
use_custom_palette <- length(structure_colors) > 0

results <- readRDS(results_rds)
tess3_obj <- results$tess3
coords <- results$coordinates
qmat <- qmatrix(tess3_obj, K = k)

window <- c(
  min(coords[, 1], na.rm = TRUE),
  max(coords[, 1], na.rm = TRUE),
  min(coords[, 2], na.rm = TRUE),
  max(coords[, 2], na.rm = TRUE)
)

palette <- if (use_custom_palette) {
  colors <- structure_colors[seq_len(min(k, length(structure_colors)))]
  if (length(colors) < k) {
    colors <- rep(colors, length.out = k)
  }
  CreatePalette(colors, 9)
} else {
  CreatePalette()
}
dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)

cat("Drawing tess3r interpolated map for K =", k, "...\n")
tryCatch({
  pdf(output_pdf, width = 10, height = 8)
  plot_args <- list(
    x = qmat,
    coord = coords,
    method = map_method,
    resolution = map_resolution,
    window = window,
    col.palette = palette,
    xlab = "Longitude",
    ylab = "Latitude",
    main = "",
    cex = 0.4
  )
  plot_tess3Q_fn <- getS3method("plot", "tess3Q")
  if ("interpol" %in% names(formals(plot_tess3Q_fn))) {
    plot_args$interpol <- FieldsKrigModel(interpolation_knots)
  } else {
    plot_args$interpolation.model <- FieldsKrigModel(interpolation_knots)
  }
  do.call(plot_tess3Q_fn, plot_args)
  dev.off()
}, error = function(e) {
  cat("WARNING: tess3r interpolation plot failed:", conditionMessage(e), "\n")
  if (dev.cur() != 1) dev.off()
  pdf(output_pdf, width = 10, height = 8)
  plot.new()
  text(0.5, 0.5, paste("Interpolation plot failed:", conditionMessage(e)))
  dev.off()
})

cat("Interpolated map saved to:", output_pdf, "\n")
