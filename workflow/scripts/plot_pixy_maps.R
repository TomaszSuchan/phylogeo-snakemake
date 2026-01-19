#!/usr/bin/env Rscript
# Plot pixy PI values on maps using the same approach as plot_population_map.R
# Note: Only PI makes sense for maps (FST and DXY are pairwise statistics)

library(tidyverse)
library(mapmixture)
library(RColorBrewer)
library(grid) # for unit()

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
stat_type <- snakemake@params[["stat_type"]]  # Only "pi" makes sense for maps
summary_file <- snakemake@input[["summary"]]
output_pdf <- snakemake@output[["pdf"]]
output_rds <- snakemake@output[["rds"]]

# Check if indpopdata file exists
if (is.null(indpopdata_file) || indpopdata_file == "NULL" || indpopdata_file == "" || !file.exists(indpopdata_file)) {
  stop("ERROR: indpopdata file not found or not specified. Map plotting requires an indpopdata file with geographic coordinates.")
}

# Extract map parameters
map_params <- extract_map_params(snakemake@params)

# Prepare coordinates from indpopdata
message("\n=== READING INDPOPDATA ===\n")
indpopdata <- read.table(indpopdata_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
coords_df <- indpopdata %>%
  dplyr::select(Site, Lat, Lon) %>%
  dplyr::distinct(Site, .keep_all = TRUE) %>%
  dplyr::mutate(Lat = as.numeric(Lat), Lon = as.numeric(Lon)) %>%
  dplyr::filter(!is.na(Lat) & !is.na(Lon))
message(sprintf("Coords dataframe: %d unique sites\n", nrow(coords_df)))

# Read summary file
summary_df <- read.table(summary_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Merge data with coordinates - only pi makes sense for maps (FST and DXY are pairwise statistics)
if (stat_type != "pi") {
  stop(sprintf("ERROR: stat_type '%s' is not valid for maps. Only 'pi' makes sense for population-level maps. FST and DXY are pairwise statistics and should be visualized as heatmaps instead.", stat_type))
}

# For pi: one value per population
# All data is already in popdata, so simple merge by Site = population
value_col <- "mean_pi"
plot_data <- merge(coords_df, summary_df, by.x = "Site", by.y = "population", all.x = TRUE)
legend_title <- "π"

# Remove sites without PI values and ensure values are numeric
plot_data <- plot_data[!is.na(plot_data[[value_col]]), ]
plot_data[[value_col]] <- as.numeric(plot_data[[value_col]])  # Explicitly ensure numeric
message(sprintf("Merged data: %d populations with PI values\n", nrow(plot_data)))
message(sprintf("Value range: [%.6f, %.6f]\n", 
                min(plot_data[[value_col]], na.rm = TRUE),
                max(plot_data[[value_col]], na.rm = TRUE)))
message(sprintf("Values are numeric: %s\n", is.numeric(plot_data[[value_col]])))

# Create basemap using mapmixture (with invisible pies)
message("\n=== CREATING MAPMIXTURE BASEMAP ===\n")
p <- create_basemap(coords_df, map_params)

# Add points colored by value on top of basemap
# Use shape=19 (solid) with a continuous color scale
message("\n=== ADDING COLORED POINTS ===\n")
p <- p +
  geom_point(
    data = plot_data,
    aes(x = Lon, y = Lat, color = !!sym(value_col)),
    size = map_params$point_size,
    shape = 19
  ) +
  scale_color_gradientn(
    name = legend_title,
    colours = RColorBrewer::brewer.pal(9, "YlOrBr"),
    na.value = "gray90",
    guide = guide_colorbar(
      title.position = "top",
      barwidth = 1,
      barheight = 10
    )
  )

# Labels are not shown for pixy maps (only points with values)

# No title for pixy maps
p <- p + theme(legend.position = "right")

# Save plots
message("\n=== SAVING OUTPUT ===\n")
message(sprintf("Output PDF: %s\n", output_pdf))
message(sprintf("Output RDS: %s\n", output_rds))
message(sprintf("Plot dimensions: %.1f x %.1f inches @ %d dpi\n", 
                map_params$width, map_params$height, map_params$dpi))

dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)
ggsave(
  filename = output_pdf,
  plot = p,
  width = map_params$width,
  height = map_params$height,
  dpi = map_params$dpi,
  device = "pdf"
)

dir.create(dirname(output_rds), recursive = TRUE, showWarnings = FALSE)
saveRDS(p, output_rds)

message("\n=== COMPLETED SUCCESSFULLY ===\n")

# Close log file sinks
sink(type = "message")
sink(type = "output")
close(log_file)

