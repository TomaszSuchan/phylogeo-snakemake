#!/usr/bin/env Rscript
# conStruct Plotting Script
# Generates spatial maps and barplots from conStruct analysis results

library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(rcolorbrewer)
library(scales)
library(tidyr)
library(dplyr)

# Prevent creation of Rplots.pdf
pdf(NULL)

# Redirect all output to log file
log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

# ==============================================================================
# Read inputs from Snakemake
# ==============================================================================

# Input files
results_rds <- snakemake@input[['results_rds']]
indpopdata_file <- snakemake@input[['indpopdata']]

# Output files
map_plot <- snakemake@output[['map_plot']]
map_plot_rds <- snakemake@output[['map_plot_rds']]
barplot <- snakemake@output[['barplot']]
barplot_rds <- snakemake@output[['barplot_rds']]

# Parameters
crs_param <- snakemake@params[['crs']]
boundary_param <- snakemake@params[['boundary']]
land_colour <- snakemake@params[['land_colour']]
sea_colour <- snakemake@params[['sea_colour']]

cat("=== conStruct Plotting Script ===\n")
cat("Results RDS:", results_rds, "\n")
cat("Individual population data:", indpopdata_file, "\n")
cat("================================\n\n")

# ==============================================================================
# 1. Load data
# ==============================================================================

cat("Loading conStruct results...\n")
construct_results <- readRDS(results_rds)

cat("Loading individual population data...\n")
indpopdata <- read.table(indpopdata_file, header = TRUE, sep = "\t")

# Check required columns
required_cols <- c("Ind", "Lat", "Lon")
if (!all(required_cols %in% colnames(indpopdata))) {
  stop("indpopdata must contain columns: Ind, Lat, Lon")
}

# ==============================================================================
# 2. Extract layer proportions
# ==============================================================================

cat("Extracting layer proportions...\n")

# Extract layer proportions from results
# conStruct results structure may vary, try different access methods
layer_props <- NULL

if (is.list(construct_results)) {
  # Try different possible locations for layer proportions
  if ("layer.proportions" %in% names(construct_results)) {
    layer_props <- construct_results$layer.proportions
  } else if ("results" %in% names(construct_results) && "layer.proportions" %in% names(construct_results$results)) {
    layer_props <- construct_results$results$layer.proportions
  } else if (length(construct_results) > 0 && is.list(construct_results[[1]])) {
    layer_props <- construct_results[[1]]$layer.proportions
  } else {
    # Try to find any matrix in the results
    for (i in seq_along(construct_results)) {
      if (is.matrix(construct_results[[i]]) || is.data.frame(construct_results[[i]])) {
        layer_props <- construct_results[[i]]
        break
      }
    }
  }
}

# If still not found, try reading from the layer_proportions file if it exists
if (is.null(layer_props)) {
  layer_props_file <- gsub("\\.results\\.rds$", ".layer_proportions.txt", results_rds)
  if (file.exists(layer_props_file)) {
    cat("Reading layer proportions from file:", layer_props_file, "\n")
    layer_props <- read.table(layer_props_file, header = TRUE, row.names = 1, sep = "\t")
  } else {
    stop("Could not extract layer proportions from conStruct results")
  }
}

# Ensure it's a matrix or data frame
if (!is.matrix(layer_props) && !is.data.frame(layer_props)) {
  stop("Layer proportions is not in expected format (matrix or data.frame)")
}

# Convert to data frame if needed
if (is.matrix(layer_props)) {
  layer_props <- as.data.frame(layer_props)
}

# Get individual names
if (is.null(rownames(layer_props))) {
  # Try to get from indpopdata
  ind_names <- indpopdata$Ind[1:nrow(layer_props)]
  rownames(layer_props) <- ind_names
} else {
  ind_names <- rownames(layer_props)
}

# Get number of layers
n_layers <- ncol(layer_props)
cat("Number of layers:", n_layers, "\n")
cat("Number of individuals:", nrow(layer_props), "\n")

# Rename columns if needed
if (is.null(colnames(layer_props))) {
  colnames(layer_props) <- paste0("Layer", 1:n_layers)
}

# ==============================================================================
# 3. Merge with geographic data
# ==============================================================================

cat("Merging with geographic data...\n")

# Match individuals
layer_props$Ind <- ind_names
plot_data <- merge(layer_props, indpopdata, by = "Ind", all.x = TRUE)

# Check for missing coordinates
if (any(is.na(plot_data$Lat) | is.na(plot_data$Lon))) {
  n_missing <- sum(is.na(plot_data$Lat) | is.na(plot_data$Lon))
  cat("WARNING:", n_missing, "individuals have missing coordinates\n")
  plot_data <- plot_data[!is.na(plot_data$Lat) & !is.na(plot_data$Lon), ]
}

cat("Plot data prepared:", nrow(plot_data), "individuals\n")

# ==============================================================================
# 4. Create color palette
# ==============================================================================

cat("Setting up color palette...\n")

# Use RColorBrewer palette, or create custom colors
if (n_layers <= 11) {
  colors <- brewer.pal(max(3, n_layers), "Set3")[1:n_layers]
} else {
  colors <- colorRampPalette(brewer.pal(11, "Set3"))(n_layers)
}

cat("Using", n_layers, "colors for", n_layers, "layers\n")

# ==============================================================================
# 5. Create spatial map plot
# ==============================================================================

cat("Creating spatial map plot...\n")

# Create output directory
output_dir <- dirname(map_plot)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Convert to sf object
plot_sf <- st_as_sf(plot_data, coords = c("Lon", "Lat"), crs = crs_param)

# Get world map for background
world <- ne_countries(scale = "medium", returnclass = "sf")

# Determine plot boundaries
if (!is.null(boundary_param) && boundary_param != "NULL") {
  # Parse boundary if provided as string (e.g., "c(xmin=-15, xmax=16, ymin=40, ymax=62)")
  # For now, use data extent
  xmin <- min(plot_data$Lon, na.rm = TRUE)
  xmax <- max(plot_data$Lon, na.rm = TRUE)
  ymin <- min(plot_data$Lat, na.rm = TRUE)
  ymax <- max(plot_data$Lat, na.rm = TRUE)
} else {
  xmin <- min(plot_data$Lon, na.rm = TRUE) - 1
  xmax <- max(plot_data$Lon, na.rm = TRUE) + 1
  ymin <- min(plot_data$Lat, na.rm = TRUE) - 1
  ymax <- max(plot_data$Lat, na.rm = TRUE) + 1
}

# Create map plot with pie charts for each individual
# For simplicity, we'll create a scatter plot colored by dominant layer
plot_data$DominantLayer <- apply(plot_data[, colnames(layer_props)], 1, which.max)

# Create map plot
p_map <- ggplot() +
  geom_sf(data = world, fill = land_colour, color = "white", size = 0.1) +
  geom_sf(data = plot_sf, aes(color = factor(DominantLayer)), size = 2, alpha = 0.7) +
  scale_color_manual(values = colors, name = "Dominant\nLayer") +
  coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = sea_colour, color = NA),
    panel.grid = element_blank(),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.position = "right"
  ) +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = paste0("conStruct Spatial Structure (K=", n_layers, ")")
  )

# Save map plot
ggsave(map_plot, plot = p_map, width = 10, height = 8, dpi = 300)
saveRDS(p_map, file = map_plot_rds)

cat("Map plot saved to:", map_plot, "\n")

# ==============================================================================
# 6. Create barplot
# ==============================================================================

cat("Creating barplot...\n")

# Prepare data for barplot (long format)
barplot_data <- plot_data[, c("Ind", "Site", colnames(layer_props))]
barplot_long <- tidyr::pivot_longer(
  barplot_data,
  cols = all_of(colnames(layer_props)),
  names_to = "Layer",
  values_to = "Proportion"
)

# Order individuals by site (if Site column exists)
if ("Site" %in% colnames(barplot_long)) {
  barplot_long$Ind <- factor(barplot_long$Ind, 
                             levels = unique(barplot_long$Ind[order(barplot_long$Site)]))
} else {
  barplot_long$Ind <- factor(barplot_long$Ind, levels = unique(barplot_long$Ind))
}

# Create barplot
p_barplot <- ggplot(barplot_long, aes(x = Ind, y = Proportion, fill = Layer)) +
  geom_bar(stat = "identity", position = "stack", width = 1) +
  scale_fill_manual(values = colors, name = "Layer") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.position = "right",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Individual",
    y = "Layer Proportion",
    title = paste0("conStruct Layer Proportions (K=", n_layers, ")")
  )

# Add site dividers if Site column exists
if ("Site" %in% colnames(barplot_long)) {
  # Get unique sites and their positions
  site_positions <- barplot_long %>%
    dplyr::group_by(Site) %>%
    dplyr::summarise(
      start = min(as.numeric(Ind)),
      end = max(as.numeric(Ind)),
      .groups = "drop"
    )
  
  # Add vertical lines between sites
  for (i in 1:(nrow(site_positions) - 1)) {
    x_pos <- site_positions$end[i] + 0.5
    p_barplot <- p_barplot + 
      geom_vline(xintercept = x_pos, linetype = "dashed", color = "black", alpha = 0.5)
  }
}

# Save barplot
ggsave(barplot, plot = p_barplot, width = max(12, nrow(plot_data) * 0.1), height = 6, dpi = 300)
saveRDS(p_barplot, file = barplot_rds)

cat("Barplot saved to:", barplot, "\n")

# ==============================================================================
# 7. Summary
# ==============================================================================

cat("\n=== conStruct Plotting Complete ===\n")
cat("Map plot:", map_plot, "\n")
cat("Barplot:", barplot, "\n")

# Close log file
sink(type = "output")
sink(type = "message")
close(log_file)

