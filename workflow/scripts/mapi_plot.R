#!/usr/bin/env Rscript
# MAPI Plotting Script
# Modern ggplot2-compatible plotting for MAPI results

library(sf)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(RColorBrewer)

# ==============================================================================
# Read inputs from Snakemake
# ==============================================================================

# Input files
mapi_gpkg <- snakemake@input[['mapi_gpkg']]
upper_tails_gpkg <- snakemake@input[['upper_tails_gpkg']]
lower_tails_gpkg <- snakemake@input[['lower_tails_gpkg']]
indpopdata_file <- snakemake@input[['indpopdata']]

# Output files
mapi_plot <- snakemake@output[['mapi_plot']]

# Parameters
fill_var <- snakemake@params[['fill_var']]

cat("=== MAPI Plotting Script ===\n")
cat("Input MAPI results:", mapi_gpkg, "\n")
cat("Input upper tails:", upper_tails_gpkg, "\n")
cat("Input lower tails:", lower_tails_gpkg, "\n")
cat("Input indpopdata:", indpopdata_file, "\n")
cat("Output plot:", mapi_plot, "\n")
cat("Fill variable:", fill_var, "\n")
cat("==============================\n\n")

# ==============================================================================
# Modern MAPI plotting function (compatible with modern ggplot2)
# ==============================================================================

plot_mapi <- function(
  mapi_results,
  upper_tails = NULL,
  lower_tails = NULL,
  indpopdata_sf = NULL,
  fill_var = "avg_value",
  crs_plot = 4326,
  limits = NULL
) {

  # mapi_results, upper_tails, and lower_tails should be sf objects

  # --- Get MAPI bbox for restricting country boundaries ---
  mapi_bbox <- st_bbox(mapi_results)

  # --- world boundaries ---
  world <- ne_countries(scale = "medium", returnclass = "sf")

  # --- CRS harmonization ---
  world <- st_transform(world, st_crs(mapi_results))

  # --- Crop world boundaries to MAPI bbox ---
  world_cropped <- st_crop(world, mapi_bbox)

  if (!is.null(upper_tails)) upper_tails <- st_transform(upper_tails, st_crs(mapi_results))
  if (!is.null(lower_tails)) lower_tails <- st_transform(lower_tails, st_crs(mapi_results))

  # --- base plot ---
  p <- ggplot() +
    geom_sf(data = mapi_results,
            aes(fill = .data[[fill_var]]),
            color = NA) +
    scale_fill_gradientn(
      colors = rev(RColorBrewer::brewer.pal(11, "RdYlGn")),
      name = fill_var,
      na.value = "transparent"
    ) +
    # --- Country boundaries overlay (above MAPI, outline only) ---
    geom_sf(data = world_cropped,
            fill = NA,
            color = "grey30",
            linewidth = 0.4) +
    coord_sf(
      crs = st_crs(mapi_results),
      xlim = limits$x,
      ylim = limits$y,
      expand = FALSE
    ) +
    theme_void() +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8)
    )

  # --- overlays for tails (merged polygons) ---
  if (!is.null(upper_tails) && nrow(upper_tails) > 0) {
    # Merge adjacent upper tail cells
    upper_merged <- st_union(upper_tails)
    p <- p +
      geom_sf(
        data = upper_merged,
        fill = NA,
        color = "black",
        linewidth = 0.5
      )
  }

  if (!is.null(lower_tails) && nrow(lower_tails) > 0) {
    # Merge adjacent lower tail cells
    lower_merged <- st_union(lower_tails)
    p <- p +
      geom_sf(
        data = lower_merged,
        fill = NA,
        color = "grey50",
        linewidth = 0.5
      )
  }

  # --- Add individual locations if provided ---
  if (!is.null(indpopdata_sf) && nrow(indpopdata_sf) > 0) {
    # Transform to match mapi_results CRS if needed
    indpopdata_transformed <- st_transform(indpopdata_sf, st_crs(mapi_results))
    p <- p +
      geom_sf(
        data = indpopdata_transformed,
        color = "black",
        size = 0.8,
        alpha = 0.6,
        shape = 16
      )
  }

  return(p)
}

# ==============================================================================
# Load data
# ==============================================================================

cat("Loading MAPI results...\n")
mapi_results <- st_read(mapi_gpkg, layer = "euclidean_results", quiet = TRUE)
cat("Loaded", nrow(mapi_results), "grid cells\n")

cat("Loading upper tails...\n")
upper_tails <- st_read(upper_tails_gpkg, layer = "euclidean_upper", quiet = TRUE)
cat("Loaded", nrow(upper_tails), "upper tail cells\n")

cat("Loading lower tails...\n")
lower_tails <- st_read(lower_tails_gpkg, layer = "euclidean_lower", quiet = TRUE)
cat("Loaded", nrow(lower_tails), "lower tail cells\n")

# ==============================================================================
# Load individual location data
# ==============================================================================

cat("Loading individual location data...\n")
indpopdata <- read.table(indpopdata_file, header = TRUE, sep = "\t")
cat("Loaded", nrow(indpopdata), "individuals\n")

# Check required columns
required_cols <- c("Ind", "Lat", "Lon")
if (!all(required_cols %in% colnames(indpopdata))) {
  warning("indpopdata missing required columns (Ind, Lat, Lon). Individual locations will not be plotted.")
  indpopdata_sf <- NULL
} else if (any(is.na(indpopdata$Lat) | is.na(indpopdata$Lon))) {
  warning("Some individuals have missing coordinates. Individual locations will not be plotted.")
  indpopdata_sf <- NULL
} else {
  # Create sf object with geographic coordinates (assuming WGS84)
  indpopdata_sf <- st_as_sf(indpopdata, coords = c("Lon", "Lat"), crs = 4326)
  cat("Created sf object with", nrow(indpopdata_sf), "individual locations\n")
}

# ==============================================================================
# Create plot
# ==============================================================================

cat("\nGenerating plot...\n")

if (nrow(mapi_results) > 0) {
  tryCatch({
    # Use the modern plotting function
    pl <- plot_mapi(
      mapi_results = mapi_results,
      upper_tails = if(nrow(upper_tails) > 0) upper_tails else NULL,
      lower_tails = if(nrow(lower_tails) > 0) lower_tails else NULL,
      indpopdata_sf = indpopdata_sf,
      fill_var = fill_var
    )

    ggsave(mapi_plot, plot = pl, width = 10, height = 8, dpi = 300)
    cat("Plot saved to:", mapi_plot, "\n")
  }, error = function(e) {
    warning("MAPI plotting failed: ", e$message)
    # Create minimal empty plot as placeholder
    pl_empty <- ggplot() +
      theme_void()
    ggsave(mapi_plot, plot = pl_empty, width = 10, height = 8, dpi = 300)
    cat("Empty plot created (plotting function failed)\n")
  })
} else {
  warning("No MAPI results to plot")
  pl_empty <- ggplot() +
    theme_void()
  ggsave(mapi_plot, plot = pl_empty, width = 10, height = 8, dpi = 300)
}

cat("\nPlotting complete!\n")
