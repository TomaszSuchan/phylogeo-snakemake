#!/usr/bin/env Rscript
# Plot population locations on map with labels and leader lines
# Uses mapmixture for basemap (same as structure plots), then adds labels

# Install mapmixture if not available (same as structure plots)
if (!require("mapmixture", quietly = TRUE)) {
  install.packages("mapmixture", repos = "https://cloud.r-project.org/")
}

library(tidyverse)
library(mapmixture)
library(ggrepel)

# Source common map functions
# Get the script directory relative to workflow root
script_dir <- dirname(normalizePath(snakemake@script))
common_functions <- file.path(script_dir, "common_map_functions.R")
if (file.exists(common_functions)) {
  source(common_functions)
} else {
  # Fallback: try relative to current directory
  source("workflow/scripts/common_map_functions.R")
}

# Prevent creation of Rplots.pdf
pdf(NULL)

# Redirect all output to log file
log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

# Snakemake inputs/outputs
indpopdata_file <- snakemake@input[["indpopdata"]]
output_plot <- snakemake@output[["plot"]]

# Check if indpopdata file exists and is not NULL
if (is.null(indpopdata_file) || indpopdata_file == "NULL" || indpopdata_file == "" || !file.exists(indpopdata_file)) {
  stop("ERROR: indpopdata file not found or not specified. Map plotting requires a popdata file with geographic coordinates. Please specify 'popdata' in your config.yaml under parameters.")
}

# Reuse mapmixture parameters (same as structure plots)
width <- as.numeric(snakemake@params[["width"]])
height <- as.numeric(snakemake@params[["height"]])
dpi <- as.numeric(snakemake@params[["dpi"]])
boundary <- snakemake@params[["boundary"]]
crs <- as.numeric(snakemake@params[["crs"]])
basemap <- snakemake@params[["basemap"]]
land_colour <- snakemake@params[["land_colour"]]
sea_colour <- snakemake@params[["sea_colour"]]
expand <- as.logical(snakemake@params[["expand"]])
arrow <- as.logical(snakemake@params[["arrow"]])
arrow_size <- as.numeric(snakemake@params[["arrow_size"]])
arrow_position <- snakemake@params[["arrow_position"]]
scalebar <- as.logical(snakemake@params[["scalebar"]])
scalebar_size <- as.numeric(snakemake@params[["scalebar_size"]])
scalebar_position <- snakemake@params[["scalebar_position"]]
plot_title <- snakemake@params[["plot_title"]]
axis_title_size <- as.numeric(snakemake@params[["axis_title_size"]])
axis_text_size <- as.numeric(snakemake@params[["axis_text_size"]])
basemap_border <- as.logical(snakemake@params[["basemap_border"]])
basemap_border_col <- snakemake@params[["basemap_border_col"]]
basemap_border_lwd <- as.numeric(snakemake@params[["basemap_border_lwd"]])

# Label-specific parameters
point_size <- as.numeric(snakemake@params[["point_size"]])
point_color <- snakemake@params[["point_color"]]
point_shape <- as.numeric(snakemake@params[["point_shape"]])
label_size <- as.numeric(snakemake@params[["label_size"]])
label_color <- snakemake@params[["label_color"]]
label_fontface <- snakemake@params[["label_fontface"]]
show_points <- as.logical(snakemake@params[["show_points"]])
show_labels <- as.logical(snakemake@params[["show_labels"]])

# ggrepel parameters
force <- as.numeric(snakemake@params[["force"]])
force_pull <- as.numeric(snakemake@params[["force_pull"]])
# Handle Inf for max.overlaps (show all labels)
max_overlaps_val <- snakemake@params[["max_overlaps"]]
if (is.character(max_overlaps_val) && (tolower(max_overlaps_val) == "inf" || tolower(max_overlaps_val) == "infinity")) {
  max.overlaps <- Inf
} else if (is.infinite(as.numeric(max_overlaps_val))) {
  max.overlaps <- Inf
} else {
  max.overlaps <- as.numeric(max_overlaps_val)
}
min.segment.length <- as.numeric(snakemake@params[["min_segment_length"]])
segment.color <- snakemake@params[["segment_color"]]
segment.size <- as.numeric(snakemake@params[["segment_size"]])

message("\n=== READING INDPOPDATA ===\n")
indpopdata <- read.table(indpopdata_file, header = TRUE, sep = "\t")
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

# Create basemap using mapmixture (with invisible pies)
message("\n=== CREATING MAPMIXTURE BASEMAP ===\n")
p <- create_basemap(coords_df, map_params)

# Add population points if desired
if (show_points) {
  message("\n=== ADDING POPULATION POINTS ===\n")
  p <- p + geom_point(
    data = coords_df,
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
    data = coords_df,
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
message(sprintf("Plot dimensions: %.1f x %.1f inches @ %d dpi\n", width, height, dpi))

# Save plot
ggsave(
  filename = output_plot,
  plot = p,
  width = width,
  height = height,
  dpi = dpi,
  device = "pdf"
)

message("\n=== COMPLETED SUCCESSFULLY ===\n")
message(sprintf("Final map includes:\n"))
message(sprintf("  - %d sites\n", nrow(coords_df)))
if (show_points) message(sprintf("  - Population points displayed\n"))
if (show_labels) message(sprintf("  - Population labels with leader lines\n"))

# Close log file sinks
sink(type = "message")
sink(type = "output")
close(log_file)

