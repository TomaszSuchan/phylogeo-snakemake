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
if (!requireNamespace("ggnewscale", quietly = TRUE)) {
  stop(
    "plot_tess3_ggmap requires ggnewscale (elevation DEM fill vs ggtess3Q tiles). ",
    "Re-run install_tess3 or rebuild the tess3 conda env after adding r-ggnewscale."
  )
}
suppressPackageStartupMessages({
  library(ggnewscale)
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

params <- snakemake_rule_params()

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

set_ggplot_layer_alpha <- function(plt, alpha) {
  if (is.null(plt$layers) || length(plt$layers) == 0) {
    return(plt)
  }
  for (i in seq_along(plt$layers)) {
    plt$layers[[i]]$aes_params$alpha <- alpha
  }
  plt
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

k <- as.integer(params[["k"]])
map_resolution <- as.integer(unlist(params[["map_resolution"]]))
interpolation_knots <- as.integer(params[["interpolation_knots"]])
structure_colors <- unlist(params[["structure_colors"]])
point_size <- as.numeric(params[["point_size"]])
tile_alpha <- as.numeric(params[["tile_alpha"]])
use_elevation_bg <- isTRUE(params[["use_elevation_bg"]])
basemap <- params[["basemap"]]
config_crs <- as.numeric(params[["crs"]])
# ggtess3Q tiles are always WGS84 lon/lat; plot in 4326 even when map_background.crs differs.
crs <- 4326L
if (config_crs != 4326L) {
  message(
    "plot_tess3_ggmap uses EPSG:4326 (ggtess3Q is geographic); ",
    "map_background.crs=", config_crs, " applies to other map rules only.\n"
  )
}

plot_title <- params[["plot_title"]]
if (is.null(plot_title)) {
  plot_title <- ""
}

map_params <- list(
  boundary = params[["boundary"]],
  crs = crs,
  basemap = basemap,
  land_colour = params[["land_colour"]],
  sea_colour = params[["sea_colour"]],
  expand = as.logical(params[["expand"]]),
  arrow = as.logical(params[["arrow"]]),
  arrow_size = as.numeric(params[["arrow_size"]]),
  arrow_position = params[["arrow_position"]],
  scalebar = as.logical(params[["scalebar"]]),
  scalebar_size = as.numeric(params[["scalebar_size"]]),
  scalebar_position = params[["scalebar_position"]],
  plot_title = plot_title,
  axis_title_size = as.numeric(params[["axis_title_size"]]),
  axis_text_size = as.numeric(params[["axis_text_size"]]),
  basemap_border = as.logical(params[["basemap_border"]]),
  basemap_border_col = params[["basemap_border_col"]],
  basemap_border_lwd = as.numeric(params[["basemap_border_lwd"]])
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

cat(
  "Adding tess3r ggtess3Q ancestry surface (cropped to land), resolution c(",
  paste(map_resolution, collapse = ", "),
  ") from map_background width/height/dpi...\n"
)
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
if (!is.na(tile_alpha) && tile_alpha < 1) {
  tess_plt <- set_ggplot_layer_alpha(tess_plt, tile_alpha)
}

# Separate fill scales: DEM uses continuous fill; ggtess3Q uses scale_fill_identity().
combined_plt <- base_plt + ggnewscale::new_scale_fill()
combined_plt <- append_ggplot_layers(combined_plt, tess_plt)

coord_df <- as.data.frame(coords)
colnames(coord_df) <- c("Lon", "Lat")

combined_plt <- combined_plt +
  ggplot2::geom_point(
    data = coord_df,
    ggplot2::aes(x = Lon, y = Lat),
    size = point_size,
    colour = "black"
  )

dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)
ggsave_pdf(
  filename = output_pdf,
  plot = combined_plt,
  width = as.numeric(params[["width"]]),
  height = as.numeric(params[["height"]]),
  dpi = as.numeric(params[["dpi"]])
)
saveRDS(combined_plt, file = output_rds)

cat("Saved tess3r ggplot map:", output_pdf, "\n")
