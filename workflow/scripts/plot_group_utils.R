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

# Parse a sort_by / sort_order parameter; NULL or empty means unset.
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

column_setting <- function(column_settings, col, field) {
  if (is.null(column_settings) || is.null(column_settings[[col]])) {
    return(NULL)
  }
  block <- column_settings[[col]]
  if (is.null(block) || !(field %in% names(block))) {
    return(NULL)
  }
  block[[field]]
}

# Explicit value order, then any remaining levels alphabetically.
explicit_level_order <- function(present, sort_order) {
  sort_order <- group_sort_by(sort_order)
  present <- unique(as.character(present[!is.na(present) & present != ""]))
  if (is.null(sort_order)) {
    return(sort(present))
  }
  c(sort_order[sort_order %in% present], sort(setdiff(present, sort_order)))
}

# Factor levels for a group column.
# sort_order: explicit values of `col` (preferred).
# sort_by: one metadata column name, or a legacy explicit value vector.
# Neither: alphabetical.
group_levels <- function(data, col, sort_by = NULL, sort_order = NULL) {
  group_levels_vals <- unique(as.character(data[[col]][!is.na(data[[col]]) & data[[col]] != ""]))
  if (!is.null(group_sort_by(sort_order))) {
    return(explicit_level_order(group_levels_vals, sort_order))
  }

  sort_by <- group_sort_by(sort_by)
  if (is.null(sort_by)) {
    return(sort(group_levels_vals))
  }

  if (length(sort_by) == 1 && sort_by %in% colnames(data) && sort_by != col) {
    keep <- !is.na(data[[col]]) & as.character(data[[col]]) != ""
    order_df <- data[keep, , drop = FALSE]
    order_df <- order_df[!duplicated(as.character(order_df[[col]])), , drop = FALSE]
    order_df <- order_df[order(order_df[[sort_by]], as.character(order_df[[col]])), , drop = FALSE]
    return(as.character(order_df[[col]]))
  }

  explicit_level_order(group_levels_vals, sort_by)
}

# Order barplot group labels.
# site_order, when set, is the full label list and wins.
# Otherwise the first label column with sort_order is used (matching the label
# itself, or ordering labels by that column's values). sort_by on that column,
# then fallback_sort_by, then alphabetical order, are the later options.
order_group_labels <- function(label_df, label_columns, column_settings = NULL,
                               fallback_sort_by = NULL, label_col = "Site") {
  sites <- unique(as.character(label_df[[label_col]]))
  for (col in label_columns) {
    sort_order <- group_sort_by(column_setting(column_settings, col, "sort_order"))
    if (is.null(sort_order)) {
      next
    }
    if (any(sort_order %in% sites)) {
      return(explicit_level_order(sites, sort_order))
    }
    if (col %in% colnames(label_df) && any(sort_order %in% as.character(label_df[[col]]))) {
      value_order <- explicit_level_order(label_df[[col]], sort_order)
      label_df$.key <- match(as.character(label_df[[col]]), value_order)
      label_df <- label_df[order(label_df$.key, label_df[[label_col]]), , drop = FALSE]
      return(as.character(label_df[[label_col]]))
    }
  }

  for (col in label_columns) {
    sort_by <- group_sort_by(column_setting(column_settings, col, "sort_by"))
    if (!is.null(sort_by)) {
      return(group_levels(label_df, label_col, sort_by = sort_by))
    }
  }

  group_levels(label_df, label_col, sort_by = fallback_sort_by)
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

population_levels <- function(populations, popdata, sort_by = NULL, site_col = "Site",
                              sort_order = NULL) {
  meta <- population_metadata(populations, popdata, site_col)
  if (nrow(meta) == 0) {
    return(explicit_level_order(populations, sort_order))
  }
  group_levels(meta, "population", sort_by = sort_by, sort_order = sort_order)
}
