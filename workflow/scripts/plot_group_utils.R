# Shared helpers for metadata-grouped plots (ROH boxplots, pixy barplots, etc.).

# Named level -> color map from parameters.{section}.{group_col}.colors, or NULL.
group_fill_values <- function(group_colors_param) {
  if (is.null(group_colors_param)) {
    return(NULL)
  }
  colors <- unlist(group_colors_param, use.names = TRUE)
  if (length(colors) == 0 || all(names(colors) == "")) {
    return(NULL)
  }
  colors
}

# Parse parameters.{section}.{group_col}.sort_by; NULL means alphabetical order.
group_sort_by <- function(group_sort_by_param) {
  if (is.null(group_sort_by_param)) {
    return(NULL)
  }
  sort_by <- as.character(unlist(group_sort_by_param, use.names = FALSE))
  sort_by <- sort_by[!is.na(sort_by) & sort_by != "" & sort_by != "NULL"]
  if (length(sort_by) == 0) {
    return(NULL)
  }
  sort_by
}

# Factor levels for a group column: explicit list, column-based sort, or alphabetical.
group_levels <- function(data, col, sort_by = NULL) {
  group_levels_vals <- unique(data[[col]][!is.na(data[[col]])])
  sort_by <- group_sort_by(sort_by)
  if (is.null(sort_by)) {
    return(sort(group_levels_vals))
  }

  if (length(sort_by) == 1 && sort_by %in% colnames(data) && sort_by != col) {
    order_df <- data %>%
      dplyr::filter(!is.na(.data[[col]])) %>%
      dplyr::distinct(.data[[col]], .keep_all = TRUE) %>%
      dplyr::arrange(.data[[sort_by]], .data[[col]])
    return(order_df[[col]])
  }

  c(sort_by[sort_by %in% group_levels_vals], sort(setdiff(group_levels_vals, sort_by)))
}

# Build one metadata row per population for ordering pixy barplot x-axis labels.
population_metadata <- function(populations, popdata, site_col = "Site") {
  populations <- unique(as.character(populations))
  if (length(populations) == 0) {
    return(data.frame())
  }

  if (site_col %in% colnames(popdata) && all(populations %in% popdata[[site_col]])) {
    idx <- match(populations, popdata[[site_col]])
    meta <- popdata[idx, , drop = FALSE]
    meta$population <- meta[[site_col]]
    return(meta)
  }

  meta <- data.frame(population = populations, stringsAsFactors = FALSE)
  for (col in colnames(popdata)) {
    if (col %in% c(site_col, "Lat", "Lon", "Ind", "Sample")) {
      next
    }
    vals <- unique(as.character(popdata[[col]]))
    if (all(populations %in% vals)) {
      meta[[col]] <- meta$population
    }
  }
  meta
}

population_levels <- function(populations, popdata, sort_by = NULL, site_col = "Site") {
  meta <- population_metadata(populations, popdata, site_col)
  if (nrow(meta) == 0) {
    return(sort(unique(as.character(populations))))
  }
  group_levels(meta, "population", sort_by)
}
