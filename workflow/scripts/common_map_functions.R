#!/usr/bin/env Rscript
# Common mapping functions for population genetics plots
# Used by plot_population_map.R and pixy map plots

library(tidyverse)
library(mapmixture)
library(ggrepel)
library(sf)
library(terra)

#' Build mapmixture-style ggplot with a terra raster basemap but without the dummy
#' legend geom (discrete fill). mapmixture::mapmixture adds geom_point + scale_fill_manual
#' for the cluster legend; that clashes with ggspatial::layer_spatial raster (continuous fill)
#' and triggers vctrs errors in recent ggplot2.
#' @noRd
.basemap_raster_without_discrete_legend <- function(
    coords_df,
    params,
    basemap_raster,
    boundary_parsed,
    dummy_qmatrix,
    cluster_cols
) {
  crs <- params$crs
  admixture_df <- mapmixture::standardise_data(dummy_qmatrix, type = "admixture")
  coords_std <- mapmixture::standardise_data(coords_df, type = "coordinates")
  cs <- as.character(coords_std[[1]])
  uq <- as.character(unique(admixture_df[[1]]))
  if (!identical(cs, uq)) {
    stop(
      "Site names in coordinates do not match admixture dummy matrix (mapmixture check)."
    )
  }
  admixture_tf <- mapmixture::transform_admix_data(admixture_df)
  admix_coords <- mapmixture::merge_coords_data(coords_std, admixture_tf)
  world_boundaries <- rnaturalearthdata::countries50[, c("geometry")]
  world_boundaries <- sf::st_transform(world_boundaries, crs = crs)
  admix_coords <- mapmixture::transform_df_coords(admix_coords, crs = crs)
  if (is.null(boundary_parsed)) {
    boundary <- mapmixture::calc_default_bbox(admix_coords, expand = 0.10)
  } else {
    boundary <- mapmixture::transform_bbox(boundary_parsed, crs)
    if (is.null(names(boundary)) && length(boundary) == 4L) {
      names(boundary) <- c("xmin", "xmax", "ymin", "ymax")
    }
  }
  # Crop in the raster's native CRS, then project once only if CRS differs from plot CRS
  # (naturalearth_basemap.tif is written already in map_background.crs to avoid warp gaps).
  crop_ext <- basemap_crop_extent_in_raster_crs(boundary, crs, basemap_raster)
  r <- terra::crop(basemap_raster, crop_ext, snap = "out")
  plot_crs_txt <- paste0("EPSG:", as.integer(crs))
  if (!isTRUE(tryCatch(terra::same.crs(r, plot_crs_txt), error = function(e) FALSE))) {
    r <- terra::project(r, plot_crs_txt, method = basemap_terra_project_method(r))
  }
  # ggspatial renders multi-band as RGB when SpatRaster has RGB() set (e.g. NE1 HR land cover).
  if (!isTRUE(params$raster_is_elevation_dem) && terra::nlyr(r) >= 3L) {
    r <- terra::subset(r, 1:3)
    suppressWarnings(tryCatch(terra::RGB(r) <- 1:3, error = function(e) NULL))
  }
  # Greys only for elevatr single-band DEM; never for Natural Earth / other imagery.
  plt <- ggplot2::ggplot() + ggspatial::layer_spatial(r)
  use_grey_dem_scale <- isTRUE(params$raster_is_elevation_dem) && terra::nlyr(r) == 1L
  if (isTRUE(use_grey_dem_scale)) {
    plt <- plt +
      ggplot2::scale_fill_distiller(
        palette = "Greys",
        na.value = params$sea_colour,
        direction = 1,
        guide = "none"
      )
  }
  plt <- plt +
    ggplot2::coord_sf(
      xlim = c(boundary[["xmin"]], boundary[["xmax"]]),
      ylim = c(boundary[["ymin"]], boundary[["ymax"]]),
      expand = params$expand,
      crs = crs
    ) +
    mapmixture::add_pie_charts(
      df = admix_coords,
      admix_columns = 4:ncol(admix_coords),
      lat_column = "lat",
      lon_column = "lon",
      pie_size = 0,
      pie_colours = cluster_cols,
      border = 0,
      border_col = "transparent",
      opacity = 0
    ) +
    ggplot2::ggtitle(params$plot_title) +
    ggplot2::xlab("Longitude") +
    ggplot2::ylab("Latitude") +
    ggplot2::theme(
      axis.text = ggplot2::element_text(colour = "black", size = params$axis_text_size),
      axis.title = ggplot2::element_text(colour = "black", size = params$axis_title_size),
      panel.grid = ggplot2::element_line(colour = "white", linewidth = 0.1),
      panel.background = ggplot2::element_rect(fill = params$sea_colour),
      panel.border = ggplot2::element_rect(fill = NA, colour = "black", linewidth = 0.3),
      plot.title = ggplot2::element_text(
        size = 12,
        face = "bold",
        margin = ggplot2::margin(0, 0, 10, 0)
      )
    )
  if (isTRUE(params$scalebar)) {
    height_size <- params$scalebar_size * 0.15
    width_size <- params$scalebar_size * 0.10
    text_size <- params$scalebar_size * 0.5
    plt <- plt +
      ggspatial::annotation_scale(
        data = world_boundaries,
        location = params$scalebar_position,
        width_hint = width_size,
        bar_cols = c("black", "white"),
        line_width = 0.5,
        height = grid::unit(height_size, "cm"),
        text_cex = text_size
      )
  }
  if (isTRUE(params$arrow)) {
    height_size <- params$arrow_size * 0.3
    width_size <- params$arrow_size * 0.3
    text_size <- params$arrow_size * 2
    pad_size <- ifelse(
      isTRUE(params$scalebar) && params$scalebar_position == params$arrow_position,
      params$arrow_size * 0.5,
      0.25
    )
    plt <- plt +
      ggspatial::annotation_north_arrow(
        data = world_boundaries,
        which_north = "true",
        location = params$arrow_position,
        height = grid::unit(height_size, "cm"),
        width = grid::unit(width_size, "cm"),
        pad_y = grid::unit(pad_size, "cm"),
        style = ggspatial::north_arrow_orienteering(
          text_size = text_size,
          line_width = 0.5
        )
      )
  }
  plt
}

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
  
  # Parse basemap if provided (never use == on SpatRaster vs "NULL" — terra dispatches [==])
  bm <- params$basemap
  basemap <- if (length(bm) == 0L || is.null(bm)) {
    NULL
  } else if (is.character(bm) && length(bm) == 1L && bm == "NULL") {
    NULL
  } else {
    bm
  }
  
  cluster_cols <- c("#000000")
  if (inherits(basemap, "SpatRaster")) {
    p <- .basemap_raster_without_discrete_legend(
      coords_df,
      params,
      basemap,
      boundary,
      dummy_qmatrix,
      cluster_cols
    )
  } else {
    p <- mapmixture::mapmixture(
      admixture_df = dummy_qmatrix,
      coords_df = coords_df,
      cluster_cols = cluster_cols,
      boundary = boundary,
      crs = params$crs,
      basemap = basemap,
      pie_size = 0,
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
  }
  
  p <- p + ggplot2::theme(legend.position = "none")
  
  if (!inherits(basemap, "SpatRaster")) {
    if (!is.null(p$layers) && length(p$layers) > 0) {
      p$layers <- p$layers[sapply(p$layers, function(l) {
        if (is.null(l$mapping) || !("fill" %in% names(l$mapping))) return(TRUE)
        !inherits(l$geom, "GeomPoint")
      })]
    }
    if (!is.null(p$scales) && length(p$scales$scales) > 0) {
      drop_fill <- sapply(p$scales$scales, function(s) {
        aesv <- s$aesthetics
        if (is.null(aesv) || !("fill" %in% aesv)) return(FALSE)
        grepl("^ScaleDiscrete", class(s)[1])
      })
      if (any(drop_fill)) {
        p$scales$scales <- p$scales$scales[!drop_fill]
      }
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

#' Fractional bbox padding (same sign logic as mapmixture::calc_default_bbox)
expand_bbox_mapmixture_style <- function(bbox, expand_frac = 0.10) {
  if (is.null(expand_frac) || expand_frac <= 0) {
    return(bbox)
  }
  b <- bbox
  b["xmin"] <- ifelse(b["xmin"] < 0, b["xmin"] + b["xmin"] * expand_frac, b["xmin"] - b["xmin"] * expand_frac)
  b["xmax"] <- ifelse(b["xmax"] < 0, b["xmax"] + abs(b["xmax"] * expand_frac), b["xmax"] + b["xmax"] * expand_frac)
  b["ymin"] <- ifelse(b["ymin"] < 0, b["ymin"] + b["ymin"] * expand_frac, b["ymin"] - b["ymin"] * expand_frac)
  b["ymax"] <- ifelse(b["ymax"] < 0, b["ymax"] + abs(b["ymax"] * expand_frac), b["ymax"] + b["ymax"] * expand_frac)
  b
}

#' WGS84 bounding box (named xmin,xmax,ymin,ymax) for elevatr from sites + optional map boundary
wgs84_bbox_for_elevation <- function(coords_df, boundary_parsed, plot_crs, expand_frac = 0.10) {
  coords_df <- coords_df %>%
    dplyr::filter(!is.na(.data$Lat), !is.na(.data$Lon))
  if (nrow(coords_df) == 0) {
    stop("No valid coordinates for elevation bbox")
  }
  if (is.null(boundary_parsed)) {
    # Robust bbox from numeric ranges (avoids edge cases in sf::st_bbox() when
    # coordinates are read/converted oddly in some environments).
    xmin <- min(coords_df$Lon, na.rm = TRUE)
    xmax <- max(coords_df$Lon, na.rm = TRUE)
    ymin <- min(coords_df$Lat, na.rm = TRUE)
    ymax <- max(coords_df$Lat, na.rm = TRUE)
    bb <- c(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
    bb <- expand_bbox_mapmixture_style(bb, expand_frac)
  } else {
    # Accept boundary in multiple shapes:
    # - list(xmin=..., xmax=..., ymin=..., ymax=...)  (e.g. parsed from YAML dict)
    # - numeric vector c(xmin=..., xmax=..., ymin=..., ymax=...) (common in config as R code)
    # - unnamed numeric length-4 vector (assume xmin,xmax,ymin,ymax)
    b <- boundary_parsed
    if (is.list(b)) {
      b <- unlist(b, use.names = TRUE)
    }
    if (is.numeric(b) && is.null(names(b)) && length(b) == 4L) {
      names(b) <- c("xmin", "xmax", "ymin", "ymax")
    }
    if (!is.numeric(b) || length(b) < 4L) {
      stop("map_boundary must be numeric with xmin,xmax,ymin,ymax")
    }
    if (!all(c("xmin", "xmax", "ymin", "ymax") %in% names(b))) {
      stop("map_boundary must have names xmin,xmax,ymin,ymax")
    }

    # Convention: map_boundary is always provided in EPSG:4326 (lon/lat degrees),
    # regardless of the target plot CRS.
    bb <- sf::st_bbox(
      c(xmin = b[["xmin"]], xmax = b[["xmax"]], ymin = b[["ymin"]], ymax = b[["ymax"]]),
      crs = sf::st_crs(4326)
    )
  }
  structure(
    c(xmin = unname(bb["xmin"]), xmax = unname(bb["xmax"]),
      ymin = unname(bb["ymin"]), ymax = unname(bb["ymax"])),
    names = c("xmin", "xmax", "ymin", "ymax")
  )
}

#' Choose elevatr zoom z from figure size and geographic span (Web Mercator tile resolution heuristic)
elevatr_z_auto <- function(bbox_wgs84, width_in, height_in, dpi, z_min = 6L, z_max = 14L) {
  lat_mid <- (unname(bbox_wgs84["ymin"]) + unname(bbox_wgs84["ymax"])) / 2 * pi / 180
  dx_deg <- abs(unname(bbox_wgs84["xmax"]) - unname(bbox_wgs84["xmin"]))
  dy_deg <- abs(unname(bbox_wgs84["ymax"]) - unname(bbox_wgs84["ymin"]))
  ground_m <- max(dx_deg * 111320 * cos(lat_mid), dy_deg * 110574)
  px <- max(width_in * dpi, height_in * dpi)
  if (!is.finite(ground_m) || ground_m <= 0 || px <= 0) {
    return(as.integer(z_min))
  }
  target_m_per_px <- ground_m / px
  for (z in z_min:z_max) {
    res <- 20037508.34 * cos(lat_mid) / (2^z * 256)
    if (res <= target_m_per_px * 1.2) {
      return(as.integer(z))
    }
  }
  as.integer(z_max)
}

#' Resampling for terra::project: "near" keeps colour tables; "bilinear" for RGB / float DEM.
#' @noRd
basemap_terra_project_method <- function(r) {
  if (terra::nlyr(r) >= 3L) {
    return("bilinear")
  }
  tab <- suppressWarnings(tryCatch(
    terra::coltab(r, layer = 1L),
    error = function(e1) tryCatch(terra::coltab(r[[1L]]), error = function(e2) NULL)
  ))
  if (is.data.frame(tab) && nrow(tab) > 0L) {
    return("near")
  }
  "bilinear"
}

#' Map plot boundary (in plot CRS) to a SpatExtent in the CRS of `ras` for terra::crop.
#' @noRd
basemap_crop_extent_in_raster_crs <- function(boundary, plot_epsg_int, ras) {
  plot_crs_txt <- paste0("EPSG:", as.integer(plot_epsg_int))
  rx <- terra::crs(ras)
  plot_ext <- terra::ext(
    as.numeric(boundary[["xmin"]]), as.numeric(boundary[["xmax"]]),
    as.numeric(boundary[["ymin"]]), as.numeric(boundary[["ymax"]])
  )
  if (is.null(rx) || !nzchar(as.character(rx))) {
    return(plot_ext)
  }
  same <- tryCatch(
    terra::same.crs(ras, plot_crs_txt),
    error = function(e) identical(trimws(as.character(rx)), trimws(plot_crs_txt))
  )
  if (isTRUE(same)) {
    return(plot_ext)
  }
  p <- sf::st_as_sfc(sf::st_bbox(c(
    xmin = as.numeric(boundary[["xmin"]]),
    xmax = as.numeric(boundary[["xmax"]]),
    ymin = as.numeric(boundary[["ymin"]]),
    ymax = as.numeric(boundary[["ymax"]])
  ), crs = sf::st_crs(plot_crs_txt)))
  p2 <- sf::st_transform(p, sf::st_crs(rx))
  exo <- sf::st_bbox(p2)
  terra::ext(exo[["xmin"]], exo[["xmax"]], exo[["ymin"]], exo[["ymax"]])
}

#' Resolve basemap SpatRaster from optional cached elevation tif, optional rnaturalearth
#' cache, optional GeoTIFF path, or NULL.
resolve_map_basemap <- function(use_elevation_bg, inputs_named_list, basemap_param) {
  ueb <- isTRUE(use_elevation_bg)
  if (ueb) {
    tif <- inputs_named_list[["elevation_basemap"]]
    if (is.null(tif) || length(tif) == 0 || !nzchar(tif) || !file.exists(tif)) {
      stop("raster basemap elevatr selected but elevation_basemap input is missing or not found")
    }
    r <- terra::rast(tif)
    if (terra::nlyr(r) > 1L) {
      r <- r[[1L]]
    }
    return(r)
  }
  ne <- inputs_named_list[["naturalearth_basemap"]]
  if (!is.null(ne) && length(ne) > 0L) {
    nep <- as.character(ne)[1L]
    if (nzchar(nep) && file.exists(nep)) {
      return(terra::rast(nep))
    }
    stop("raster basemap ne1/gray_earth selected but naturalearth_basemap input is missing or not found")
  }
  gt <- inputs_named_list[["basemap_geotiff"]]
  if (!is.null(gt) && length(gt) > 0L) {
    path <- as.character(gt)[1L]
    if (nzchar(path) && file.exists(path)) {
      return(terra::rast(path))
    }
  }
  if (is.character(basemap_param) && length(basemap_param) == 1L && nzchar(basemap_param) &&
      basemap_param != "NULL" &&
      grepl("\\.(tif|tiff|vrt)$", basemap_param, ignore.case = TRUE) &&
      file.exists(basemap_param)) {
    return(terra::rast(basemap_param))
  }
  if (length(basemap_param) == 0 || is.null(basemap_param) || identical(basemap_param, "NULL")) {
    return(NULL)
  }
  basemap_param
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

