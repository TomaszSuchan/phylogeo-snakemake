#!/usr/bin/env Rscript
# Plot population locations on map with labels and leader lines
# Uses mapmixture for basemap (same as structure plots), then adds labels

# Source common map functions (loads mapmixture and shared helpers)
# Try to get script directory, fallback to relative path
common_functions <- tryCatch({
  script_dir <- dirname(normalizePath(snakemake@script))
  file.path(script_dir, "common_map_functions.R")
}, error = function(e) "workflow/scripts/common_map_functions.R")
if (file.exists(common_functions)) {
  source(common_functions)
} else {
  # Final fallback
  source("workflow/scripts/common_map_functions.R")
}

params <- snakemake_rule_params()

# Prevent creation of Rplots.pdf
pdf(NULL)

# Redirect all output to log file
log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

# Snakemake inputs/outputs
indpopdata_file <- snakemake@input[["indpopdata"]]
output_plot <- snakemake@output[["plot"]]
output_plot_rds <- snakemake@output[["plot_rds"]]

# Check if indpopdata file exists and is not NULL
if (is.null(indpopdata_file) || indpopdata_file == "NULL" || indpopdata_file == "" || !file.exists(indpopdata_file)) {
  stop("ERROR: indpopdata file not found or not specified. Map plotting requires a popdata file with geographic coordinates. Please specify 'popdata' in your config.yaml under parameters.")
}

# Reuse mapmixture parameters (same as structure plots)
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

# Label-specific parameters
point_size <- as.numeric(params[["point_size"]])
point_color <- params[["point_color"]]
point_shape <- as.numeric(params[["point_shape"]])
label_size <- as.numeric(params[["label_size"]])
label_color <- params[["label_color"]]
label_fontface <- params[["label_fontface"]]
show_points <- as.logical(params[["show_points"]])
show_labels <- as.logical(params[["show_labels"]])

# ggrepel parameters
force <- as.numeric(params[["force"]])
force_pull <- as.numeric(params[["force_pull"]])
max.overlaps <- params[["max_overlaps"]]
if (is.null(max.overlaps)) {
  max.overlaps <- Inf
} else if (is.character(max.overlaps) && length(max.overlaps) == 1L &&
           toupper(max.overlaps) %in% c("INF", "INFINITY", ".INF")) {
  max.overlaps <- Inf
} else {
  max.overlaps <- as.numeric(max.overlaps)
}
min.segment.length <- as.numeric(params[["min_segment_length"]])
segment.color <- params[["segment_color"]]
segment.size <- as.numeric(params[["segment_size"]])

message("\n=== READING INDPOPDATA ===\n")
indpopdata <- read.table(indpopdata_file, header = TRUE, sep = "\t")
indpopdata <- coerce_indpopdata_lat_lon(indpopdata)
message(sprintf("indpopdata file: %s\n", indpopdata_file))
message(sprintf("indpopdata dimensions: %d rows x %d columns\n", nrow(indpopdata), ncol(indpopdata)))
message(sprintf("Columns: %s\n", paste(colnames(indpopdata), collapse = ", ")))

# Create coords_df with unique sites
message("\n=== CREATING COORDINATES DATAFRAME ===\n")
coords_df <- indpopdata %>%
  select(Site, Lat, Lon) %>%
  distinct(Site, .keep_all = TRUE)

message(sprintf("Coords dataframe has %d unique sites\n", nrow(coords_df)))
message(sprintf("Coordinate range: Lat [%.2f, %.2f], Lon [%.2f, %.2f]\n",
            min(coords_df$Lat), max(coords_df$Lat),
            min(coords_df$Lon), max(coords_df$Lon)))

# Prepare map parameters
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

# Transform coordinates to match map CRS if needed
message("\n=== TRANSFORMING COORDINATES ===\n")
message(sprintf("Target CRS: %d\n", crs))
coords_df_transformed <- transform_coordinates(coords_df, crs)
message(sprintf("Coordinates transformed. Range: Lon [%.2f, %.2f], Lat [%.2f, %.2f]\n",
            min(coords_df_transformed$Lon), max(coords_df_transformed$Lon),
            min(coords_df_transformed$Lat), max(coords_df_transformed$Lat)))

# Create basemap using mapmixture (with invisible pies)
message("\n=== CREATING MAPMIXTURE BASEMAP ===\n")
p <- create_basemap(coords_df, map_params)

# Add population points if desired
if (show_points) {
  message("\n=== ADDING POPULATION POINTS ===\n")
  p <- p + geom_point(
    data = coords_df_transformed,
    aes(x = Lon, y = Lat),
    size = point_size,
    color = point_color,
    shape = point_shape
  )
}

# Add labels with leader lines using ggrepel
if (show_labels) {
  message("\n=== ADDING LABELS WITH GGREPEL ===\n")
  p <- p + geom_text_repel(
    data = coords_df_transformed,
    aes(x = Lon, y = Lat, label = Site),
    size = label_size,
    color = label_color,
    fontface = label_fontface,
    force = force,
    force_pull = force_pull,
    max.overlaps = max.overlaps,
    min.segment.length = min.segment.length,
    segment.color = segment.color,
    segment.size = segment.size,
    point.padding = 0.2,
    box.padding = 0.5
  )
}

message("\n=== SAVING OUTPUT ===\n")
message(sprintf("Output PDF: %s\n", output_plot))
message(sprintf("Output RDS: %s\n", output_plot_rds))
message(sprintf("Plot dimensions: %.1f x %.1f inches @ %d dpi\n", width, height, dpi))

# Save plot
ggsave_pdf(
  filename = output_plot,
  plot = p,
  width = width,
  height = height,
  dpi = dpi
)
saveRDS(p, file = output_plot_rds)

message("\n=== COMPLETED SUCCESSFULLY ===\n")
message(sprintf("Final map includes:\n"))
message(sprintf("  - %d sites\n", nrow(coords_df)))
if (show_points) message(sprintf("  - Population points displayed\n"))
if (show_labels) message(sprintf("  - Population labels with leader lines\n"))

# Close log file sinks
sink(type = "message")
sink(type = "output")
close(log_file)

