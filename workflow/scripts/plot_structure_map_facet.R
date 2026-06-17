#!/usr/bin/env Rscript
# Align pre-rendered per-K map grobs (K >= 2) with patchwork.

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

pdf(NULL)

script_dir <- tryCatch(
  dirname(normalizePath(snakemake@script)),
  error = function(e) "workflow/scripts"
)
source(file.path(script_dir, "plot_panel_utils.R"))
source(file.path(script_dir, "plot_ggsave_utils.R"))

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
ncol_maps <- as.integer(params[["ncol"]])
panel_width <- as.numeric(params[["panel_width"]])
panel_height <- as.numeric(params[["panel_height"]])

map_input_names <- grep("^map_k", names(snakemake@input), value = TRUE)
if (length(map_input_names) == 0L) {
  stop("No per-K map RDS inputs found (expected names like map_k2, map_k3, ...).")
}

extract_k <- function(name) as.integer(sub("^map_k", "", name))
map_ks <- vapply(map_input_names, extract_k, integer(1))
ord <- order(map_ks)
map_input_names <- map_input_names[ord]
map_ks <- map_ks[ord]

message(sprintf("Building %s map facet from K = %s", method_label, paste(map_ks, collapse = ", ")))

map_panels <- lapply(snakemake@input[map_input_names], read_panel_grob)
ncol_maps <- max(1L, min(ncol_maps, length(map_panels)))
nrow_maps <- max(1L, ceiling(length(map_panels) / ncol_maps))

combined_plot <- wrap_plots(map_panels, ncol = ncol_maps) +
  plot_layout(widths = 1, heights = 1) +
  plot_annotation(
    tag_levels = list(paste0("K = ", map_ks)),
    theme = theme(
      plot.tag = element_text(face = "bold", size = 11),
      plot.tag.position = c(0.5, 1)
    )
  )

page_width <- ncol_maps * panel_width
page_height <- nrow_maps * panel_height

message(sprintf(
  "Saving facet PDF: %s (%d columns, %.1f x %.1f in)",
  output_pdf, ncol_maps, page_width, page_height
))
ggsave_pdf(
  filename = output_pdf,
  plot = combined_plot,
  width = page_width,
  height = page_height
)

message(sprintf("Saving facet RDS: %s", output_rds))
saveRDS(combined_plot, file = output_rds)

message("Map facet plot completed successfully.")
