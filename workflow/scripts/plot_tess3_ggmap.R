#!/usr/bin/env Rscript
# tess3r ggtess3Q ancestry surface on the shared mapmixture basemap.
#
# ggtess3Q returns a geom_tile layer in WGS84 lon/lat. The mapmixture basemap is
# built first; tess tiles and scale_fill_identity() are appended on top (not below).

pdf(NULL)

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

suppressPackageStartupMessages({
  library(ggplot2)
  library(sf)
  library(tidyverse)
  library(mapmixture)
  library(terra)
})

if (!requireNamespace("tess3r", quietly = TRUE)) {
  stop("tess3r is not installed in this rule environment.")
}
suppressPackageStartupMessages({
  library(tess3r)
})

common_functions <- tryCatch({
  script_dir <- dirname(normalizePath(snakemake@script))
  file.path(script_dir, "common_map_functions.R")
}, error = function(e) "workflow/scripts/common_map_functions.R")
if (file.exists(common_functions)) {
  source(common_functions)
} else {
  source("workflow/scripts/common_map_functions.R")
}

append_ggplot_layers <- function(base_plt, overlay_plt) {
  if (is.null(overlay_plt$layers) || length(overlay_plt$layers) == 0) {
    return(base_plt)
  }
  base_plt$layers <- c(base_plt$layers, overlay_plt$layers)
  if (!is.null(overlay_plt$scales) && length(overlay_plt$scales$scales) > 0) {
    base_plt$scales$scales <- c(base_plt$scales$scales, overlay_plt$scales$scales)
  }
  base_plt
}

tess3_cluster_palette <- function(structure_colors, n_clusters, palette_length = 9) {
  if (length(structure_colors) == 0) {
    return(CreatePalette())
  }
  colors <- structure_colors[seq_len(min(n_clusters, length(structure_colors)))]
  if (length(colors) < n_clusters) {
    colors <- rep(colors, length.out = n_clusters)
  }
  CreatePalette(colors, palette_length)
}

tess3_map_polygon <- function() {
  if (!requireNamespace("rnaturalearth", quietly = TRUE)) {
    return(NULL)
  }
  tryCatch(
    rnaturalearth::ne_countries(scale = "medium", returnclass = "sp"),
    error = function(e) NULL
  )
}

results_rds <- snakemake@input[["results_rds"]]
indpopdata_file <- snakemake@input[["indpopdata"]]
output_pdf <- snakemake@output[["pdf"]]
output_rds <- snakemake@output[["rds"]]

k <- as.integer(snakemake@params[["k"]])
map_resolution <- as.integer(unlist(snakemake@params[["map_resolution"]]))
interpolation_knots <- as.integer(snakemake@params[["interpolation_knots"]])
structure_colors <- unlist(snakemake@params[["structure_colors"]])
point_size <- as.numeric(snakemake@params[["point_size"]])
use_elevation_bg <- isTRUE(snakemake@params[["use_elevation_bg"]])
basemap <- snakemake@params[["basemap"]]
crs <- as.numeric(snakemake@params[["crs"]])

if (crs != 4326) {
  stop(
    "plot_tess3_ggmap requires map_background.crs = 4326 because ggtess3Q ",
    "interpolates in geographic longitude/latitude."
  )
}

plot_title <- snakemake@params[["plot_title"]]
if (is.null(plot_title)) {
  plot_title <- ""
}

map_params <- list(
  boundary = snakemake@params[["boundary"]],
  crs = crs,
  basemap = basemap,
  land_colour = snakemake@params[["land_colour"]],
  sea_colour = snakemake@params[["sea_colour"]],
  expand = as.logical(snakemake@params[["expand"]]),
  arrow = as.logical(snakemake@params[["arrow"]]),
  arrow_size = as.numeric(snakemake@params[["arrow_size"]]),
  arrow_position = snakemake@params[["arrow_position"]],
  scalebar = as.logical(snakemake@params[["scalebar"]]),
  scalebar_size = as.numeric(snakemake@params[["scalebar_size"]]),
  scalebar_position = snakemake@params[["scalebar_position"]],
  plot_title = plot_title,
  axis_title_size = as.numeric(snakemake@params[["axis_title_size"]]),
  axis_text_size = as.numeric(snakemake@params[["axis_text_size"]]),
  basemap_border = as.logical(snakemake@params[["basemap_border"]]),
  basemap_border_col = snakemake@params[["basemap_border_col"]],
  basemap_border_lwd = as.numeric(snakemake@params[["basemap_border_lwd"]])
)
map_params$basemap <- resolve_map_basemap(use_elevation_bg, snakemake@input, basemap)
map_params$raster_is_elevation_dem <- isTRUE(use_elevation_bg)

results <- readRDS(results_rds)
tess3_obj <- results$tess3
coords <- results$coordinates
qmat <- qmatrix(tess3_obj, K = k)
n_clusters <- ncol(as.matrix(qmat))

indpopdata <- read.table(
  indpopdata_file,
  header = TRUE,
  sep = "\t",
  comment.char = "",
  fill = TRUE,
  blank.lines.skip = TRUE
)
indpopdata <- coerce_indpopdata_lat_lon(indpopdata)
coords_df <- indpopdata %>%
  dplyr::select(Site, Lat, Lon) %>%
  dplyr::distinct(Site, .keep_all = TRUE)

window <- c(
  min(coords[, 1], na.rm = TRUE),
  max(coords[, 1], na.rm = TRUE),
  min(coords[, 2], na.rm = TRUE),
  max(coords[, 2], na.rm = TRUE)
)
boundary <- map_params$boundary
if (!(length(boundary) == 0 || is.null(boundary) || boundary == "NULL")) {
  boundary_parsed <- tryCatch(
    eval(parse(text = boundary)),
    error = function(e) NULL
  )
  if (!is.null(boundary_parsed) && length(boundary_parsed) == 4L) {
    window <- boundary_parsed
  }
}

map_polygon <- tess3_map_polygon()
palette <- tess3_cluster_palette(structure_colors, n_clusters)

cat("Building mapmixture basemap...\n")
base_plt <- create_basemap(coords_df, map_params)

cat("Adding tess3r ggtess3Q ancestry surface (cropped to land)...\n")
tess_plt <- ggtess3Q(
  Q = qmat,
  coord = coords,
  resolution = map_resolution,
  window = window,
  background = TRUE,
  map.polygon = map_polygon,
  interpolation.model = FieldsKrigModel(interpolation_knots),
  col.palette = palette
)

# Basemap first, tess tiles on top; include scale_fill_identity() from ggtess3Q.
combined_plt <- append_ggplot_layers(base_plt, tess_plt)

coord_df <- as.data.frame(coords)
colnames(coord_df) <- c("Lon", "Lat")

combined_plt <- combined_plt +
  ggplot2::geom_point(
    data = coord_df,
    ggplot2::aes(x = Lon, y = Lat),
    size = point_size,
    colour = "black",
    alpha = 0.7
  )

dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)
ggplot2::ggsave(
  filename = output_pdf,
  plot = combined_plt,
  width = as.numeric(snakemake@params[["width"]]),
  height = as.numeric(snakemake@params[["height"]]),
  dpi = as.numeric(snakemake@params[["dpi"]])
)
saveRDS(combined_plt, file = output_rds)

cat("Saved tess3r ggplot map:", output_pdf, "\n")
