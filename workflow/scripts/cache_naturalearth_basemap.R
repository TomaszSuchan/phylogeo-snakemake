#!/usr/bin/env Rscript
# Clip + warp Natural Earth HR global GeoTIFF (from fetch_naturalearth_hr_global) to study bbox
# and map_background.crs → naturalearth_basemap.tif. Does not download — CRS/bbox changes only rerun this rule.

log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

out_tif <- snakemake@output[[1]]
global_tif <- snakemake@input[["hr_global"]]

message("\n=== CACHE NATURAL EARTH RASTER (crop + warp) ===\n")
message(sprintf("global source: %s\n", global_tif))
message(sprintf("output: %s\n", out_tif))

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

suppressPackageStartupMessages(library(dplyr))
library(terra)

if (!file.exists(global_tif) || !nzchar(global_tif)) {
  stop("hr_global input missing: ", global_tif)
}

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
message(sprintf(
  "WGS84 clip bbox: xmin=%.5f xmax=%.5f ymin=%.5f ymax=%.5f\n",
  bb["xmin"], bb["xmax"], bb["ymin"], bb["ymax"]
))

rst <- terra::rast(global_tif)
e <- terra::ext(
  as.numeric(bb["xmin"]),
  as.numeric(bb["xmax"]),
  as.numeric(bb["ymin"]),
  as.numeric(bb["ymax"])
)
rst_c <- terra::crop(rst, e, snap = "out")
if (terra::ncell(rst_c) == 0L) {
  stop("Natural Earth crop returned empty extent; check coordinates and map_boundary.")
}
message(sprintf(
  "Cropped to study extent (WGS84): %d cells, nlyr=%d\n",
  terra::ncell(rst_c),
  terra::nlyr(rst_c)
))

target_epsg <- as.integer(round(plot_crs))
target_txt <- paste0("EPSG:", target_epsg)
rst_out <- rst_c
if (target_epsg != 4326L) {
  pm <- basemap_terra_project_method(rst_c)
  message(sprintf("Reprojecting to %s (method=%s)...\n", target_txt, pm))
  rst_out <- terra::project(rst_c, target_txt, method = pm)
  corners <- matrix(
    c(
      bb["xmin"], bb["ymin"],
      bb["xmin"], bb["ymax"],
      bb["xmax"], bb["ymin"],
      bb["xmax"], bb["ymax"]
    ),
    ncol = 2L,
    byrow = TRUE
  )
  cv <- terra::vect(corners, crs = "EPSG:4326")
  cv_t <- terra::project(cv, target_txt)
  ex_t <- terra::ext(cv_t)
  rst_out <- terra::crop(rst_out, ex_t, snap = "out")
  message(sprintf(
    "After warp: %d cells, nlyr=%d, crs=%s\n",
    terra::ncell(rst_out),
    terra::nlyr(rst_out),
    target_txt
  ))
} else {
  message("Plot CRS is 4326; basemap left in WGS84.\n")
}

dir.create(dirname(out_tif), recursive = TRUE, showWarnings = FALSE)
terra::writeRaster(
  rst_out,
  out_tif,
  overwrite = TRUE,
  gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2", "TILED=YES")
)
message(sprintf("Wrote %s\n", out_tif))

sink(type = "message")
sink(type = "output")
close(log_file)
