#!/usr/bin/env Rscript
# Common mapping functions for population genetics plots
# Used by plot_population_map.R and pixy map plots

library(tidyverse)
library(mapmixture)
library(ggrepel)
library(sf)

#' Create basemap using mapmixture (without pie charts)
#' This function creates a styled basemap with all the same parameters as structure plots
create_basemap <- function(coords_df, params) {
  # Create dummy Q matrix (single cluster, all 1s) for mapmixture
  dummy_qmatrix <- data.frame(
    Site = coords_df$Site,
    Ind = coords_df$Site,  # Use site name as individual ID
    Cluster1 = rep(1, nrow(coords_df))
  )
  
  # Parse boundary if provided
  boundary <- if (length(params$boundary) == 0 || is.null(params$boundary) || params$boundary == "NULL") {
    NULL
  } else {
    eval(parse(text = params$boundary))
  }
  
  # Parse basemap if provided
  basemap <- if (length(params$basemap) == 0 || is.null(params$basemap) || params$basemap == "NULL") {
    NULL
  } else {
    params$basemap
  }
  
  # Create mapmixture plot with pie_size = 0 (invisible pies)
  p <- mapmixture::mapmixture(
    admixture_df = dummy_qmatrix,
    coords_df = coords_df,
    cluster_cols = c("#000000"),  # Dummy color (won't be visible)
    boundary = boundary,
    crs = params$crs,
    basemap = basemap,
    pie_size = 0,  # Make pies invisible
    pie_border = 0,
    pie_border_col = "transparent",
    pie_opacity = 0,
    land_colour = params$land_colour,
    sea_colour = params$sea_colour,
    expand = params$expand,
    arrow = params$arrow,
    arrow_size = params$arrow_size,
    arrow_position = params$arrow_position,
    scalebar = params$scalebar,
    scalebar_size = params$scalebar_size,
    scalebar_position = params$scalebar_position,
    plot_title = params$plot_title,
    axis_title_size = params$axis_title_size,
    axis_text_size = params$axis_text_size,
    basemap_border = params$basemap_border,
    basemap_border_col = params$basemap_border_col,
    basemap_border_lwd = params$basemap_border_lwd
  )
  
  # Remove legend
  p <- p + ggplot2::theme(legend.position = "none")
  
  # Remove layers and scales that use fill (for invisible pies)
  # This allows pixy maps to add their own fill scale without conflicts
  # Remove layers that use fill aesthetic
  if (!is.null(p$layers) && length(p$layers) > 0) {
    p$layers <- p$layers[sapply(p$layers, function(l) {
      if (is.null(l$mapping)) return(TRUE)
      !("fill" %in% names(l$mapping))
    })]
  }
  
  # Remove fill scales
  if (!is.null(p$scales) && length(p$scales$scales) > 0) {
    fill_scale_indices <- sapply(p$scales$scales, function(s) {
      if (is.null(s$aesthetics)) return(FALSE)
      "fill" %in% s$aesthetics
    })
    if (any(fill_scale_indices)) {
      p$scales$scales <- p$scales$scales[!fill_scale_indices]
    }
  }
  
  return(p)
}

#' Prepare coordinates dataframe from popmap and popdata
prepare_coordinates <- function(popmap_file, popdata_file) {
  # Read popmap to get unique populations
  popmap <- read.table(popmap_file, header = FALSE, sep = "\t")
  colnames(popmap) <- c("Ind", "Site")
  unique_sites <- unique(popmap$Site)
  
  # Read popdata
  popdata <- read.table(popdata_file, header = TRUE, sep = "\t")
  
  # Filter popdata to only include populations from popmap
  coords_df <- popdata %>%
    dplyr::filter(Site %in% unique_sites) %>%
    dplyr::select(Site, Lat, Lon) %>%
    dplyr::distinct(Site, .keep_all = TRUE) %>%
    dplyr::mutate(Lat = as.numeric(Lat), Lon = as.numeric(Lon))
  
  # Check for any NA values after conversion
  if (any(is.na(coords_df$Lat)) || any(is.na(coords_df$Lon))) {
    warning("WARNING: Some coordinates could not be converted to numeric and will be removed.")
    coords_df <- coords_df %>%
      dplyr::filter(!is.na(Lat) & !is.na(Lon))
  }
  
  # Check for populations in popmap but not in popdata
  missing_sites <- setdiff(unique_sites, coords_df$Site)
  if (length(missing_sites) > 0) {
    warning(sprintf("WARNING: %d populations in popmap are missing from popdata:\n", length(missing_sites)))
    warning(paste(missing_sites, collapse = ", "))
    warning("\nThese populations will not be plotted.\n")
  }
  
  return(coords_df)
}

#' Transform coordinates from WGS84 (4326) to target CRS
#' This ensures points align with the map when using projected CRS
transform_coordinates <- function(coords_df, target_crs) {
  # If CRS is 4326 (WGS84), no transformation needed
  if (target_crs == 4326) {
    return(coords_df)
  }
  
  # Create sf object with geographic coordinates (WGS84)
  coords_sf <- st_as_sf(coords_df, coords = c("Lon", "Lat"), crs = 4326)
  
  # Transform to target CRS
  coords_transformed <- st_transform(coords_sf, crs = target_crs)
  
  # Extract transformed coordinates
  coords_matrix <- st_coordinates(coords_transformed)
  
  # Update coordinates in dataframe
  coords_df$Lon <- coords_matrix[, "X"]
  coords_df$Lat <- coords_matrix[, "Y"]
  
  return(coords_df)
}

#' Extract map background parameters from snakemake params
#' Also extracts analysis-specific parameters (point_size, labels, etc.) if present
extract_map_params <- function(snakemake_params) {
  list(
    width = as.numeric(snakemake_params[["width"]]),
    height = as.numeric(snakemake_params[["height"]]),
    dpi = as.numeric(snakemake_params[["dpi"]]),
    boundary = snakemake_params[["boundary"]],
    crs = as.numeric(snakemake_params[["crs"]]),
    basemap = snakemake_params[["basemap"]],
    land_colour = snakemake_params[["land_colour"]],
    sea_colour = snakemake_params[["sea_colour"]],
    expand = as.logical(snakemake_params[["expand"]]),
    arrow = as.logical(snakemake_params[["arrow"]]),
    arrow_size = as.numeric(snakemake_params[["arrow_size"]]),
    arrow_position = snakemake_params[["arrow_position"]],
    scalebar = as.logical(snakemake_params[["scalebar"]]),
    scalebar_size = as.numeric(snakemake_params[["scalebar_size"]]),
    scalebar_position = snakemake_params[["scalebar_position"]],
    plot_title = snakemake_params[["plot_title"]],
    axis_title_size = as.numeric(snakemake_params[["axis_title_size"]]),
    axis_text_size = as.numeric(snakemake_params[["axis_text_size"]]),
    basemap_border = as.logical(snakemake_params[["basemap_border"]]),
    basemap_border_col = snakemake_params[["basemap_border_col"]],
    basemap_border_lwd = as.numeric(snakemake_params[["basemap_border_lwd"]]),
    point_size = as.numeric(snakemake_params[["point_size"]]),
    point_color = snakemake_params[["point_color"]],
    point_shape = as.numeric(snakemake_params[["point_shape"]]),
    label_size = as.numeric(snakemake_params[["label_size"]]),
    label_color = snakemake_params[["label_color"]],
    label_fontface = snakemake_params[["label_fontface"]],
    show_points = as.logical(snakemake_params[["show_points"]]),
    show_labels = as.logical(snakemake_params[["show_labels"]]),
    force = as.numeric(snakemake_params[["force"]]),
    force_pull = as.numeric(snakemake_params[["force_pull"]]),
    max_overlaps = snakemake_params[["max_overlaps"]],
    min_segment_length = as.numeric(snakemake_params[["min_segment_length"]]),
    segment_color = snakemake_params[["segment_color"]],
    segment_size = as.numeric(snakemake_params[["segment_size"]])
  )
}

