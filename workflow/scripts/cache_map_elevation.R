#!/usr/bin/env Rscript
# Writes results/{project}/maps/elevation_basemap.tif via elevatr.
# Output CRS matches map_background.crs: bbox stays WGS84 for z heuristics, then the
# polygon is transformed before get_elev_raster() so elevatr warps once (same idea as NE cache).

log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

out_tif <- snakemake@output[[1]]

dir.create(dirname(out_tif), recursive = TRUE, showWarnings = FALSE)

message("\n=== CACHE MAP ELEVATION (elevatr) ===\n")

library(sf)
library(terra)
library(elevatr)

common_functions <- tryCatch({
  script_dir <- dirname(normalizePath(snakemake@script))
  file.path(script_dir, "common_map_functions.R")
}, error = function(e) "workflow/scripts/common_map_functions.R")
if (file.exists(common_functions)) {
  source(common_functions)
} else {
  source("workflow/scripts/common_map_functions.R")
}

indpopdata_file <- snakemake@input[["indpopdata"]]
boundary <- snakemake@params[["boundary"]]
plot_crs <- as.numeric(snakemake@params[["crs"]])
width_in <- as.numeric(snakemake@params[["width"]])
height_in <- as.numeric(snakemake@params[["height"]])
dpi <- as.numeric(snakemake@params[["dpi"]])
elevatr_z_param <- snakemake@params[["elevatr_z"]]

message(sprintf("indpopdata: %s\n", indpopdata_file))
message(sprintf("output: %s\n", out_tif))

indpopdata <- read.table(indpopdata_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
coords_df <- indpopdata %>%
  dplyr::select(Site, Lat, Lon) %>%
  dplyr::distinct(Site, .keep_all = TRUE) %>%
  dplyr::mutate(Lat = as.numeric(Lat), Lon = as.numeric(Lon))

boundary_parsed <- NULL
if (!(length(boundary) == 0 || is.null(boundary) || identical(boundary, "NULL"))) {
  boundary_parsed <- eval(parse(text = boundary))
}

bb <- wgs84_bbox_for_elevation(coords_df, boundary_parsed, plot_crs, expand_frac = 0.10)
message(sprintf("WGS84 bbox: xmin=%.5f xmax=%.5f ymin=%.5f ymax=%.5f\n",
                bb["xmin"], bb["xmax"], bb["ymin"], bb["ymax"]))

if (is.null(elevatr_z_param) || length(elevatr_z_param) == 0 ||
    (is.character(elevatr_z_param) && identical(elevatr_z_param, "NULL"))) {
  z <- elevatr_z_auto(bb, width_in, height_in, dpi)
  message(sprintf("elevatr z (auto): %d\n", z))
} else {
  z <- as.integer(elevatr_z_param)
  message(sprintf("elevatr z (config): %d\n", z))
}

poly_wgs84 <- sf::st_as_sf(sf::st_as_sfc(sf::st_bbox(bb, crs = sf::st_crs(4326))))
plot_epsg <- as.integer(round(plot_crs))
message("Requesting elevatr download in EPSG:4326 (WGS84)\n")
r_elev <- elevatr::get_elev_raster(poly_wgs84, z = z)
r_terra <- terra::rast(r_elev)

# elevatr returns a raster in lon/lat; project explicitly if the plot CRS differs.
if (plot_epsg != 4326L) {
  message(sprintf("Projecting elevation raster to EPSG:%d\n", plot_epsg))
  r_terra <- terra::project(r_terra, paste0("EPSG:", plot_epsg), method = "bilinear")
}

terra::writeRaster(r_terra, out_tif, overwrite = TRUE, datatype = "FLT4S")
message(sprintf("Wrote %s\n", out_tif))

sink(type = "message")
sink(type = "output")
close(log_file)
