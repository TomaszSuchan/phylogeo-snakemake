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
mapi_plot_rds <- snakemake@output[['mapi_plot_rds']]

# Parameters
fill_var <- snakemake@params[['fill_var']]
width <- as.numeric(snakemake@params[['width']])
height <- as.numeric(snakemake@params[['height']])
dpi <- as.numeric(snakemake@params[['dpi']])
boundary <- snakemake@params[['boundary']]
crs_plot <- as.numeric(snakemake@params[['crs_plot']])
land_colour <- snakemake@params[['land_colour']]
sea_colour <- snakemake@params[['sea_colour']]
expand <- as.logical(snakemake@params[['expand']])
plot_title <- snakemake@params[['plot_title']]
axis_title_size <- as.numeric(snakemake@params[['axis_title_size']])
axis_text_size <- as.numeric(snakemake@params[['axis_text_size']])
basemap_border <- as.logical(snakemake@params[['basemap_border']])
basemap_border_col <- snakemake@params[['basemap_border_col']]
basemap_border_lwd <- as.numeric(snakemake@params[['basemap_border_lwd']])
point_size <- as.numeric(snakemake@params[['point_size']])
point_color <- snakemake@params[['point_color']]
point_alpha <- as.numeric(snakemake@params[['point_alpha']])
tail_linewidth <- as.numeric(snakemake@params[['tail_linewidth']])
upper_tail_color <- snakemake@params[['upper_tail_color']]
lower_tail_color <- snakemake@params[['lower_tail_color']]

cat("=== MAPI Plotting Script ===\n")
cat("Input MAPI results:", mapi_gpkg, "\n")
cat("Input upper tails:", upper_tails_gpkg, "\n")
cat("Input lower tails:", lower_tails_gpkg, "\n")
cat("Input indpopdata:", indpopdata_file, "\n")
cat("Output plot:", mapi_plot, "\n")
cat("Output plot RDS:", mapi_plot_rds, "\n")
cat("Fill variable:", fill_var, "\n")
cat("CRS (plot):", crs_plot, "\n")
cat("==============================\n\n")

# ==============================================================================
# Modern MAPI plotting function (compatible with modern ggplot2)
# ==============================================================================

parse_boundary_limits <- function(boundary_param, target_crs) {
  if (is.null(boundary_param) || length(boundary_param) == 0) return(NULL)
  if (is.character(boundary_param) && length(boundary_param) == 1 &&
      (!nzchar(boundary_param) || boundary_param == "NULL")) {
    return(NULL)
  }
  b <- boundary_param
  if (is.character(b) && length(b) == 1) {
    b <- tryCatch(eval(parse(text = b)), error = function(e) NULL)
  }
  if (is.null(b)) return(NULL)
  if (is.list(b)) b <- unlist(b, use.names = TRUE)
  if (is.numeric(b) && is.null(names(b)) && length(b) == 4) {
    names(b) <- c("xmin", "xmax", "ymin", "ymax")
  }
  if (!is.numeric(b) || !all(c("xmin", "xmax", "ymin", "ymax") %in% names(b))) {
    return(NULL)
  }

  bbox_wgs84 <- st_bbox(c(
    xmin = as.numeric(b[["xmin"]]),
    xmax = as.numeric(b[["xmax"]]),
    ymin = as.numeric(b[["ymin"]]),
    ymax = as.numeric(b[["ymax"]])
  ), crs = st_crs(4326))

  bbox_target <- if (as.integer(target_crs) == 4326L) {
    bbox_wgs84
  } else {
    st_bbox(st_transform(st_as_sfc(bbox_wgs84), st_crs(target_crs)))
  }

  list(
    x = c(unname(bbox_target[["xmin"]]), unname(bbox_target[["xmax"]])),
    y = c(unname(bbox_target[["ymin"]]), unname(bbox_target[["ymax"]]))
  )
}

default_boundary_limits <- function(source_sf, target_crs, expand_frac = 0.10) {
  # Match mapmixture::calc_default_bbox pattern:
  # compute bbox in WGS84, expand proportionally, then transform to plot CRS.
  source_wgs84 <- st_transform(source_sf, 4326)
  bb <- st_bbox(source_wgs84)
  b <- c(
    xmin = unname(bb[["xmin"]]),
    xmax = unname(bb[["xmax"]]),
    ymin = unname(bb[["ymin"]]),
    ymax = unname(bb[["ymax"]])
  )

  if (!is.null(expand_frac) && is.finite(expand_frac) && expand_frac > 0) {
    b["xmin"] <- ifelse(b["xmin"] < 0, b["xmin"] + b["xmin"] * expand_frac, b["xmin"] - b["xmin"] * expand_frac)
    b["xmax"] <- ifelse(b["xmax"] < 0, b["xmax"] + abs(b["xmax"] * expand_frac), b["xmax"] + b["xmax"] * expand_frac)
    b["ymin"] <- ifelse(b["ymin"] < 0, b["ymin"] + b["ymin"] * expand_frac, b["ymin"] - b["ymin"] * expand_frac)
    b["ymax"] <- ifelse(b["ymax"] < 0, b["ymax"] + abs(b["ymax"] * expand_frac), b["ymax"] + b["ymax"] * expand_frac)
  }

  bbox_wgs84 <- st_bbox(c(
    xmin = as.numeric(b["xmin"]),
    xmax = as.numeric(b["xmax"]),
    ymin = as.numeric(b["ymin"]),
    ymax = as.numeric(b["ymax"])
  ), crs = st_crs(4326))

  bbox_target <- if (as.integer(target_crs) == 4326L) {
    bbox_wgs84
  } else {
    st_bbox(st_transform(st_as_sfc(bbox_wgs84), st_crs(target_crs)))
  }

  list(
    x = c(unname(bbox_target[["xmin"]]), unname(bbox_target[["xmax"]])),
    y = c(unname(bbox_target[["ymin"]]), unname(bbox_target[["ymax"]]))
  )
}

merge_tail_cells <- function(tails_sf, tol_frac = 0.03) {
  if (is.null(tails_sf) || nrow(tails_sf) == 0) return(NULL)

  tails_sf <- suppressWarnings(st_make_valid(tails_sf))
  tails_sf <- st_collection_extract(tails_sf, "POLYGON", warn = FALSE)
  if (nrow(tails_sf) == 0) return(NULL)

  med_area <- suppressWarnings(stats::median(as.numeric(st_area(tails_sf)), na.rm = TRUE))
  tol <- sqrt(med_area) * tol_frac
  if (!is.finite(tol) || tol <= 0) tol <- 1

  merged_geom <- tails_sf |>
    st_buffer(tol) |>
    st_union() |>
    st_buffer(-tol) |>
    st_make_valid()

  merged_geom <- st_collection_extract(merged_geom, "POLYGON", warn = FALSE)
  if (length(merged_geom) == 0) return(NULL)

  st_as_sf(data.frame(id = seq_along(merged_geom)), geometry = merged_geom)
}

plot_mapi <- function(
  mapi_results,
  upper_tails = NULL,
  lower_tails = NULL,
  indpopdata_sf = NULL,
  fill_var = "avg_value",
  crs_plot = 4326,
  boundary_limits = NULL,
  land_colour = "#d9d9d9",
  sea_colour = "#deebf7",
  expand = FALSE,
  plot_title = "",
  axis_title_size = 10,
  axis_text_size = 8,
  basemap_border = TRUE,
  basemap_border_col = "black",
  basemap_border_lwd = 0.1,
  point_size = 0.8,
  point_color = "black",
  point_alpha = 0.6,
  tail_linewidth = 0.4,
  upper_tail_color = "#B2182B",
  lower_tail_color = "#1B7837"
) {

  # mapi_results, upper_tails, and lower_tails should be sf objects

  mapi_results <- st_transform(mapi_results, crs_plot)
  if (!is.null(upper_tails) && nrow(upper_tails) > 0) upper_tails <- st_transform(upper_tails, crs_plot)
  if (!is.null(lower_tails) && nrow(lower_tails) > 0) lower_tails <- st_transform(lower_tails, crs_plot)

  # --- world boundaries ---
  world <- ne_countries(scale = "medium", returnclass = "sf")

  # --- CRS harmonization ---
  world <- st_transform(world, st_crs(mapi_results))

  mapi_bbox <- st_bbox(mapi_results)
  crop_bbox <- if (!is.null(boundary_limits)) {
    st_bbox(c(
      xmin = boundary_limits$x[1], xmax = boundary_limits$x[2],
      ymin = boundary_limits$y[1], ymax = boundary_limits$y[2]
    ), crs = st_crs(mapi_results))
  } else {
    mapi_bbox
  }
  world_cropped <- suppressWarnings(st_crop(world, crop_bbox))
  if (nrow(world_cropped) == 0) world_cropped <- world

  # --- base plot ---
  p <- ggplot() +
    # Land background uses configured land_colour and sits below MAPI cells.
    geom_sf(
      data = world_cropped,
      fill = land_colour,
      color = NA
    ) +
    geom_sf(data = mapi_results,
            aes(fill = .data[[fill_var]]),
            color = NA) +
    scale_fill_gradientn(
      colors = rev(RColorBrewer::brewer.pal(11, "RdYlGn")),
      name = fill_var,
      na.value = "transparent"
    ) +
    coord_sf(
      xlim = if (!is.null(boundary_limits)) boundary_limits$x else NULL,
      ylim = if (!is.null(boundary_limits)) boundary_limits$y else NULL,
      expand = expand,
      crs = st_crs(mapi_results)
    ) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = sea_colour, color = NA),
      panel.grid = element_line(color = "white", linewidth = 0.1),
      axis.text = element_text(size = axis_text_size),
      axis.title = element_text(size = axis_title_size),
      legend.position = "right",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.3),
      plot.margin = margin(0, 0, 0, 0)
    ) +
    labs(x = "Longitude", y = "Latitude", title = plot_title)

  # --- overlays for tails (merged polygons) ---
  if (!is.null(upper_tails) && nrow(upper_tails) > 0) {
    upper_merged <- merge_tail_cells(upper_tails)
    if (!is.null(upper_merged)) {
      p <- p +
        geom_sf(
          data = upper_merged,
          fill = NA,
          color = upper_tail_color,
          linewidth = tail_linewidth
        )
    }
  }

  if (!is.null(lower_tails) && nrow(lower_tails) > 0) {
    lower_merged <- merge_tail_cells(lower_tails)
    if (!is.null(lower_merged)) {
      p <- p +
        geom_sf(
          data = lower_merged,
          fill = NA,
          color = lower_tail_color,
          linewidth = tail_linewidth
        )
    }
  }

  # --- country boundaries overlay as transparent outlines on top ---
  p <- p +
    geom_sf(
      data = world_cropped,
      fill = NA,
      color = ifelse(isTRUE(basemap_border), basemap_border_col, NA),
      linewidth = basemap_border_lwd
    )

  # --- Add individual locations if provided ---
  if (!is.null(indpopdata_sf) && nrow(indpopdata_sf) > 0) {
    # Transform to match mapi_results CRS if needed
    indpopdata_transformed <- st_transform(indpopdata_sf, st_crs(mapi_results))
    p <- p +
      geom_sf(
        data = indpopdata_transformed,
        color = point_color,
        size = point_size,
        alpha = point_alpha,
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

output_dir <- dirname(mapi_plot)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

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

boundary_limits <- parse_boundary_limits(boundary, crs_plot)
if (is.null(boundary_limits)) {
  # Prefer unique sampling-site coordinates for map framing.
  boundary_source <- if (!is.null(indpopdata_sf) && nrow(indpopdata_sf) > 0) {
    unique(indpopdata_sf["Site"])
  } else {
    mapi_results
  }
  boundary_limits <- default_boundary_limits(boundary_source, crs_plot, expand_frac = 0.10)
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
      fill_var = fill_var,
      crs_plot = crs_plot,
      boundary_limits = boundary_limits,
      land_colour = land_colour,
      sea_colour = sea_colour,
      expand = expand,
      plot_title = plot_title,
      axis_title_size = axis_title_size,
      axis_text_size = axis_text_size,
      basemap_border = basemap_border,
      basemap_border_col = basemap_border_col,
      basemap_border_lwd = basemap_border_lwd,
      point_size = point_size,
      point_color = point_color,
      point_alpha = point_alpha,
      tail_linewidth = tail_linewidth,
      upper_tail_color = upper_tail_color,
      lower_tail_color = lower_tail_color
    )

    ggsave(mapi_plot, plot = pl, width = width, height = height, dpi = dpi)
    saveRDS(pl, file = mapi_plot_rds)
    cat("Plot saved to:", mapi_plot, "\n")
  }, error = function(e) {
    warning("MAPI plotting failed: ", e$message)
    # Create minimal empty plot as placeholder
    pl_empty <- ggplot() +
      theme_void()
    ggsave(mapi_plot, plot = pl_empty, width = width, height = height, dpi = dpi)
    saveRDS(pl_empty, file = mapi_plot_rds)
    cat("Empty plot created (plotting function failed)\n")
  })
} else {
  warning("No MAPI results to plot")
  pl_empty <- ggplot() +
    theme_void()
  ggsave(mapi_plot, plot = pl_empty, width = width, height = height, dpi = dpi)
  saveRDS(pl_empty, file = mapi_plot_rds)
}

cat("\nPlotting complete!\n")
