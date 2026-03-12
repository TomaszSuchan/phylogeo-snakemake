#!/usr/bin/env Rscript
# Create continuous ancestry probability surfaces using PopMaps approach
# Based on: https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.13902
# PopMaps creates continuous surfaces from discrete sampling locations

library(tidyverse)
library(sf)
library(terra)
library(fields)
library(akima)
library(viridis)
library(gridExtra)

# Try to load rnaturalearth, but make it optional
tryCatch({
  library(rnaturalearth)
  library(rnaturalearthdata)
  rnaturalearth_available <- TRUE
}, error = function(e) {
  rnaturalearth_available <- FALSE
  message("Warning: rnaturalearth not available, using simplified map\n")
})

# Prevent creation of Rplots.pdf
pdf(NULL)

# Redirect all output to log file
log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

# Snakemake inputs/outputs
qmatrix_file <- snakemake@input[["qmatrix"]]
indpopdata_file <- snakemake@input[["indpopdata"]]
output_pdf <- snakemake@output[["pdf"]]
output_rds <- snakemake@output[["rds"]]

# Parameters
method_name <- snakemake@params[["method"]]
k_value <- snakemake@params[["k"]]
resolution <- as.numeric(snakemake@params[["resolution"]])  # Resolution in degrees
smooth_factor <- as.numeric(snakemake@params[["smooth_factor"]])  # Smoothing parameter
crs_geographic <- as.integer(snakemake@params[["crs_geographic"]])
width <- as.numeric(snakemake@params[["width"]])
height <- as.numeric(snakemake@params[["height"]])
dpi <- as.numeric(snakemake@params[["dpi"]])

message("\n=== READING INPUT DATA ===\n")

# Read Q matrix
qmatrix <- read.table(qmatrix_file, header = FALSE, sep = "")
n_individuals <- nrow(qmatrix)
n_clusters <- ncol(qmatrix)
message(sprintf("Q matrix: %d individuals x %d clusters\n", n_individuals, n_clusters))

# Read individual population data (coordinates)
indpopdata <- read.table(indpopdata_file, header = TRUE, sep = "\t")
message(sprintf("Individual population data: %d individuals\n", nrow(indpopdata)))

# Check required columns
if (!all(c("Ind", "Lat", "Lon") %in% colnames(indpopdata))) {
  stop("indpopdata must contain columns: Ind, Lat, Lon")
}

# Match individuals between Q matrix and coordinates using indpopdata order
indpopdata_unique <- indpopdata %>%
  distinct(Ind, .keep_all = TRUE)

if (nrow(indpopdata_unique) > nrow(qmatrix)) {
  indpopdata_unique <- indpopdata_unique[1:nrow(qmatrix), ]
}

if (nrow(qmatrix) != nrow(indpopdata_unique)) {
  stop(sprintf("Q matrix has %d individuals but indpopdata has %d usable rows", nrow(qmatrix), nrow(indpopdata_unique)))
}

# Merge data
plot_data <- data.frame(
  Ind = indpopdata_unique$Ind,
  Site = indpopdata_unique$Site,
  Q = qmatrix
)
colnames(plot_data)[3:(2+n_clusters)] <- paste0("Cluster", 1:n_clusters)

# Add coordinates
plot_data <- merge(plot_data, indpopdata[, c("Ind", "Lat", "Lon")], by = "Ind", all.x = TRUE)

# Check for missing coordinates
missing_coords <- is.na(plot_data$Lat) | is.na(plot_data$Lon)
if (any(missing_coords)) {
  n_missing <- sum(missing_coords)
  message(sprintf("WARNING: %d individuals missing coordinates, removing from analysis\n", n_missing))
  plot_data <- plot_data[!missing_coords, ]
}

message(sprintf("Final data: %d individuals with coordinates\n", nrow(plot_data)))

# Get coordinate bounds
lon_range <- range(plot_data$Lon, na.rm = TRUE)
lat_range <- range(plot_data$Lat, na.rm = TRUE)
message(sprintf("Coordinate range: Lon [%.4f, %.4f], Lat [%.4f, %.4f]\n",
                lon_range[1], lon_range[2], lat_range[1], lat_range[2]))

# Create grid for interpolation
# Add buffer around data points
lon_buffer <- diff(lon_range) * 0.1
lat_buffer <- diff(lat_range) * 0.1

lon_seq <- seq(lon_range[1] - lon_buffer, lon_range[2] + lon_buffer, by = resolution)
lat_seq <- seq(lat_range[1] - lat_buffer, lat_range[2] + lat_buffer, by = resolution)

grid_data <- expand.grid(Lon = lon_seq, Lat = lat_seq)
message(sprintf("Grid created: %d x %d = %d points\n",
                length(lon_seq), length(lat_seq), nrow(grid_data)))

# Interpolate each cluster separately
message("\n=== INTERPOLATING ANCESTRY PROBABILITIES ===\n")
cluster_surfaces <- list()

for (cluster_idx in 1:n_clusters) {
  cluster_name <- paste0("Cluster", cluster_idx)
  message(sprintf("Interpolating %s...\n", cluster_name))
  
  # Get values at sample locations
  z_values <- plot_data[[cluster_name]]
  
  # Use akima::interp for interpolation (bicubic interpolation)
  interp_result <- tryCatch({
    akima::interp(
      x = plot_data$Lon,
      y = plot_data$Lat,
      z = z_values,
      xo = lon_seq,
      yo = lat_seq,
      linear = FALSE,  # Use bicubic interpolation
      extrap = TRUE,    # Extrapolate beyond data points
      duplicate = "mean"  # Handle duplicate coordinates
    )
  }, error = function(e) {
    message(sprintf("Warning: akima interpolation failed for %s, using linear interpolation\n", cluster_name))
    akima::interp(
      x = plot_data$Lon,
      y = plot_data$Lat,
      z = z_values,
      xo = lon_seq,
      yo = lat_seq,
      linear = TRUE,
      extrap = TRUE,
      duplicate = "mean"
    )
  })
  
  # Convert to data frame
  interp_df <- expand.grid(Lon = interp_result$x, Lat = interp_result$y)
  interp_df$Value <- as.vector(interp_result$z)
  
  # Ensure values are between 0 and 1 (probabilities)
  interp_df$Value <- pmax(0, pmin(1, interp_df$Value))
  
  cluster_surfaces[[cluster_name]] <- interp_df
}

# Create plots for each cluster
message("\n=== CREATING PLOTS ===\n")
plots_list <- list()

for (cluster_idx in 1:n_clusters) {
  cluster_name <- paste0("Cluster", cluster_idx)
  surface_data <- cluster_surfaces[[cluster_name]]
  
  # Create base map using sf
  # Create simple world map outline
  if (rnaturalearth_available) {
    world_map <- tryCatch({
      rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
    }, error = function(e) {
      # Fallback: create simple bounding box
      bbox <- st_bbox(c(xmin = lon_range[1] - lon_buffer, 
                        xmax = lon_range[2] + lon_buffer,
                        ymin = lat_range[1] - lat_buffer,
                        ymax = lat_range[2] + lat_buffer),
                      crs = st_crs(crs_geographic))
      st_as_sfc(bbox)
    })
  } else {
    # Create simple bounding box
    bbox <- st_bbox(c(xmin = lon_range[1] - lon_buffer, 
                      xmax = lon_range[2] + lon_buffer,
                      ymin = lat_range[1] - lat_buffer,
                      ymax = lat_range[2] + lat_buffer),
                    crs = st_crs(crs_geographic))
    world_map <- st_as_sfc(bbox)
  }
  
  # Create plot
  p <- ggplot() +
    # Add world map
    geom_sf(data = world_map, fill = "gray90", color = "gray70", linewidth = 0.3) +
    # Add interpolated surface
    geom_raster(data = surface_data, aes(x = Lon, y = Lat, fill = Value), interpolate = TRUE) +
    scale_fill_viridis_c(
      name = "Ancestry\nProbability",
      limits = c(0, 1),
      option = "plasma",
      na.value = "transparent"
    ) +
    # Add sample points
    geom_point(data = plot_data, aes(x = Lon, y = Lat), 
               color = "black", size = 0.5, alpha = 0.5) +
    coord_sf(crs = crs_geographic, expand = FALSE) +
    labs(
      title = sprintf("PopMaps: %s Ancestry Surface\n%s, K=%s", cluster_name, method_name, k_value),
      x = "Longitude",
      y = "Latitude"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      legend.position = "right"
    )
  
  plots_list[[cluster_name]] <- p
}

# Combine all cluster plots
if (length(plots_list) > 1) {
  library(gridExtra)
  combined_plot <- do.call(grid.arrange, c(plots_list, ncol = min(2, n_clusters)))
} else {
  combined_plot <- plots_list[[1]]
}

# Save plots
message("\n=== SAVING OUTPUT ===\n")
message(sprintf("Output PDF: %s\n", output_pdf))
message(sprintf("Output RDS: %s\n", output_rds))

dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)
ggsave(
  filename = output_pdf,
  plot = combined_plot,
  width = width,
  height = height,
  dpi = dpi,
  device = "pdf"
)

dir.create(dirname(output_rds), recursive = TRUE, showWarnings = FALSE)
saveRDS(combined_plot, output_rds)

message("\n=== COMPLETED SUCCESSFULLY ===\n")

# Close log file sinks
sink(type = "message")
sink(type = "output")
close(log_file)

