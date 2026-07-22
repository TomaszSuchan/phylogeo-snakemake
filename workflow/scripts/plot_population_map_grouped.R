#!/usr/bin/env Rscript
# Plot population locations on a map with points colored by an indpopdata column.

common_functions <- tryCatch({
  script_dir <- dirname(normalizePath(snakemake@script))
  file.path(script_dir, "common_map_functions.R")
}, error = function(e) "workflow/scripts/common_map_functions.R")
if (file.exists(common_functions)) {
  source(common_functions)
} else {
  source("workflow/scripts/common_map_functions.R")
}

group_utils <- tryCatch({
  script_dir <- dirname(normalizePath(snakemake@script))
  file.path(script_dir, "plot_group_utils.R")
}, error = function(e) "workflow/scripts/plot_group_utils.R")
if (file.exists(group_utils)) {
  source(group_utils)
} else {
  source("workflow/scripts/plot_group_utils.R")
}

params <- snakemake_rule_params()

pdf(NULL)

log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

indpopdata_file <- snakemake@input[["indpopdata"]]
output_plot <- snakemake@output[["plot"]]
output_plot_rds <- snakemake@output[["plot_rds"]]

group_col <- params[["group_col"]]
group_colors_param <- params[["group_colors"]]
group_sort_by_param <- params[["group_sort_by"]]

if (is.null(indpopdata_file) || indpopdata_file == "NULL" || indpopdata_file == "" || !file.exists(indpopdata_file)) {
  stop("ERROR: indpopdata file not found or not specified. Map plotting requires a popdata file with geographic coordinates.")
}
if (is.null(group_col) || !nzchar(as.character(group_col))) {
  stop("ERROR: group_col parameter is required for grouped population maps.")
}

width <- as.numeric(params[["width"]])
height <- as.numeric(params[["height"]])
dpi <- as.numeric(params[["dpi"]])
boundary <- params[["boundary"]]
crs <- as.numeric(params[["crs"]])
basemap <- params[["basemap"]]
land_colour <- params[["land_colour"]]
sea_colour <- params[["sea_colour"]]
expand <- as.logical(params[["expand"]])
arrow <- as.logical(params[["arrow"]])
arrow_size <- as.numeric(params[["arrow_size"]])
arrow_position <- params[["arrow_position"]]
scalebar <- as.logical(params[["scalebar"]])
scalebar_size <- as.numeric(params[["scalebar_size"]])
scalebar_position <- params[["scalebar_position"]]
plot_title <- params[["plot_title"]]
axis_title_size <- as.numeric(params[["axis_title_size"]])
axis_text_size <- as.numeric(params[["axis_text_size"]])
basemap_border <- as.logical(params[["basemap_border"]])
basemap_border_col <- params[["basemap_border_col"]]
basemap_border_lwd <- as.numeric(params[["basemap_border_lwd"]])
use_elevation_bg <- isTRUE(params[["use_elevation_bg"]])

color_point_size <- as.numeric(params[["color_point_size"]])
point_shape <- as.numeric(params[["point_shape"]])
map_outline <- as.logical(params[["map_outline"]])
if (is.na(map_outline)) {
  map_outline <- TRUE
}

message(sprintf("\n=== GROUPED POPULATION MAP (%s) ===\n", group_col))

message("\n=== READING INDPOPDATA ===\n")
indpopdata <- read.table(indpopdata_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
indpopdata <- coerce_indpopdata_lat_lon(indpopdata)

if (!group_col %in% colnames(indpopdata)) {
  stop(sprintf(
    "ERROR: column '%s' not found in indpopdata. Available columns: %s",
    group_col,
    paste(colnames(indpopdata), collapse = ", ")
  ))
}
if (!all(c("Site", "Lat", "Lon") %in% colnames(indpopdata))) {
  stop("ERROR: indpopdata must contain Site, Lat, and Lon columns for map plotting.")
}

plot_data <- indpopdata %>%
  dplyr::select(Site, Lat, Lon, dplyr::all_of(group_col)) %>%
  dplyr::distinct(Site, .keep_all = TRUE) %>%
  dplyr::mutate(Lat = as.numeric(Lat), Lon = as.numeric(Lon)) %>%
  dplyr::filter(!is.na(Lat), !is.na(Lon), !is.na(.data[[group_col]]))

if (nrow(plot_data) == 0) {
  stop(sprintf("ERROR: No sites with valid coordinates and '%s' values to plot.", group_col))
}

group_levels_vals <- group_levels(plot_data, group_col, group_sort_by_param)
plot_data[[group_col]] <- factor(plot_data[[group_col]], levels = group_levels_vals)
palette_vals <- group_fill_values(group_colors_param)

message(sprintf("Plotting %d sites colored by %s\n", nrow(plot_data), group_col))
message(sprintf("Groups: %s\n", paste(levels(plot_data[[group_col]]), collapse = ", ")))

map_params <- list(
  boundary = boundary,
  crs = crs,
  basemap = basemap,
  land_colour = land_colour,
  sea_colour = sea_colour,
  expand = expand,
  arrow = arrow,
  arrow_size = arrow_size,
  arrow_position = arrow_position,
  scalebar = scalebar,
  scalebar_size = scalebar_size,
  scalebar_position = scalebar_position,
  plot_title = plot_title,
  axis_title_size = axis_title_size,
  axis_text_size = axis_text_size,
  basemap_border = basemap_border,
  basemap_border_col = basemap_border_col,
  basemap_border_lwd = basemap_border_lwd
)
map_params$basemap <- resolve_map_basemap(use_elevation_bg, snakemake@input, basemap)
map_params$raster_is_elevation_dem <- isTRUE(use_elevation_bg)
map_params$elevation_style <- params[["elevation_style"]]
map_params$width <- width
map_params$dpi <- dpi

coords_df <- plot_data %>% dplyr::select(Site, Lat, Lon)

message("\n=== TRANSFORMING COORDINATES ===\n")
message(sprintf("Target CRS: %d\n", crs))
plot_data_transformed <- transform_coordinates(plot_data, crs)

message("\n=== CREATING MAPMIXTURE BASEMAP ===\n")
p <- create_basemap(coords_df, map_params)

message("\n=== ADDING COLORED POINTS ===\n")
p <- p +
  geom_point(
    data = plot_data_transformed,
    aes(x = Lon, y = Lat, color = .data[[group_col]]),
    size = color_point_size,
    shape = point_shape
  )

if (map_outline) {
  p <- p +
    geom_point(
      data = plot_data_transformed,
      aes(x = Lon, y = Lat),
      color = "black",
      size = color_point_size,
      shape = 1
    )
}

if (!is.null(palette_vals) && length(palette_vals) > 0) {
  p <- p + scale_color_manual(values = palette_vals, drop = FALSE, guide = "none")
} else {
  p <- p + scale_color_discrete(drop = FALSE, guide = "none")
}

p <- p + theme(legend.position = "none")

message("\n=== SAVING OUTPUT ===\n")
message(sprintf("Output PDF: %s\n", output_plot))
message(sprintf("Output RDS: %s\n", output_plot_rds))
message(sprintf("Plot dimensions: %.1f x %.1f inches @ %d dpi\n", width, height, dpi))

ggsave_pdf(
  filename = output_plot,
  plot = p,
  width = width,
  height = height,
  dpi = dpi
)
saveRDS(p, file = output_plot_rds)

message("\n=== COMPLETED SUCCESSFULLY ===\n")

sink(type = "message")
sink(type = "output")
close(log_file)
