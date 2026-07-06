#!/usr/bin/env Rscript
# Script to create structure barplots using mapmixture

# Redirect all output to log file before loading packages.
log_file <- snakemake@log[[1]]
dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")
on.exit({
  while (sink.number(type = "message") > 0) sink(type = "message")
  while (sink.number(type = "output") > 0) sink(type = "output")
  close(log_con)
}, add = TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(mapmixture)
})

common_functions <- tryCatch({
  script_dir <- dirname(normalizePath(snakemake@script))
  file.path(script_dir, "common_map_functions.R")
}, error = function(e) "workflow/scripts/common_map_functions.R")
if (file.exists(common_functions)) {
  source(common_functions)
} else {
  source("workflow/scripts/common_map_functions.R")
}

params <- snakemake_rule_params()

# Prevent creation of Rplots.pdf
pdf(NULL)

# Snakemake inputs/outputs
qmatrix_file <- snakemake@input[["qmatrix"]]
indpopdata_file <- snakemake@input[["indpopdata"]]
output_barplot <- snakemake@output[["barplot"]]
output_barplot_rds <- snakemake@output[["barplot_rds"]]
output_prefix <- params[["output_prefix"]]

# Plot dimension parameters
width <- as.numeric(params[["width"]])
height <- as.numeric(params[["height"]])
dpi <- as.numeric(params[["dpi"]])

# Barplot-specific parameters
site_dividers <- as.logical(params[["site_dividers"]])
divider_width <- as.numeric(params[["divider_width"]])
site_order <- params[["site_order"]]
population_sort_by <- params[["population_sort_by"]]
flip_axis <- as.logical(params[["flip_axis"]])
site_labels_angle <- as.numeric(params[["site_labels_angle"]])
population_labels <- params[["population_labels"]]

# Fixed parameters optimized for 1.2 inch height
site_ticks_size <- -0.05
# Adjust site_labels_y based on angle to prevent clipping
# For rotated labels (45 or 90 degrees), need more space
if (site_labels_angle == 90) {
  site_labels_y <- -0.5  # More space for perpendicular labels
} else if (site_labels_angle == 45) {
  site_labels_y <- -0.45  # More space for diagonal labels
} else {
  site_labels_y <- -0.35  # Default for horizontal labels
}
site_labels_size <- 2.2

# Read Q matrix (space or tab delimited, no headers)
qmatrix <- read.table(qmatrix_file, header = FALSE, sep = "")
n_individuals <- nrow(qmatrix)
n_clusters <- ncol(qmatrix)
cat(sprintf("Q matrix dimensions: %d individuals x %d clusters\n", n_individuals, n_clusters))

# Read indpopdata and use its row order to match the Q matrix order
cat(sprintf("Reading indpopdata file: %s\n", indpopdata_file))
indpopdata <- read.table(indpopdata_file, header = TRUE, sep = "\t", comment.char = "", fill = TRUE, blank.lines.skip = TRUE)
indpopdata <- indpopdata %>%
  distinct(Ind, .keep_all = TRUE)

if (nrow(indpopdata) > nrow(qmatrix)) {
  cat(sprintf("Filtering indpopdata from %d to %d individuals to match Q matrix\n", nrow(indpopdata), nrow(qmatrix)))
  indpopdata <- indpopdata[1:nrow(qmatrix), ]
}

if (nrow(qmatrix) != nrow(indpopdata)) {
  stop(sprintf("ERROR: Q matrix has %d individuals but indpopdata has %d individuals",
               nrow(qmatrix), nrow(indpopdata)))
}

# Parse population_labels:
# - NULL => hide labels
# - character vector/list => columns from indpopdata used to build labels/grouping
if (is.null(population_labels)) {
  display_population_labels <- FALSE
  population_label_columns_val <- "Site"
} else if (length(population_labels) == 1 &&
           (is.na(population_labels) || population_labels == "NULL" || population_labels == "")) {
  display_population_labels <- FALSE
  population_label_columns_val <- "Site"
} else if (is.logical(population_labels) && length(population_labels) == 1) {
  # Backward compatibility: TRUE => Site labels, FALSE => no labels.
  display_population_labels <- isTRUE(population_labels)
  population_label_columns_val <- "Site"
} else if (is.vector(population_labels) || is.list(population_labels)) {
  display_population_labels <- TRUE
  population_label_columns_val <- unlist(population_labels)
} else {
  display_population_labels <- TRUE
  population_label_columns_val <- tryCatch(
    eval(parse(text = population_labels)),
    error = function(e) "Site"
  )
}
population_label_columns_val <- unique(as.character(population_label_columns_val))
if (length(population_label_columns_val) == 0) {
  population_label_columns_val <- "Site"
}

missing_label_columns <- setdiff(population_label_columns_val, colnames(indpopdata))
if (length(missing_label_columns) > 0) {
  stop(sprintf(
    "ERROR: population_label_columns not found in indpopdata: %s",
    paste(missing_label_columns, collapse = ", ")
  ))
}

cat(sprintf(
  "Using population label columns: %s (display labels: %s)\n",
  paste(population_label_columns_val, collapse = ", "),
  ifelse(display_population_labels, "TRUE", "FALSE")
))

population_label_df <- indpopdata[, population_label_columns_val, drop = FALSE] %>%
  mutate(across(everything(), ~ {
    x <- as.character(.x)
    x[is.na(x) | x == ""] <- "NA"
    x
  }))

population_label_values <- apply(population_label_df, 1, paste, collapse = " | ")

# Combine Q matrix with individual-level data, rename colnames
qmatrix_with_data <- cbind(population_label_values, indpopdata$Ind, qmatrix)
colnames(qmatrix_with_data)[1] <- "Site"
colnames(qmatrix_with_data)[2] <- "Ind"
colnames(qmatrix_with_data) <- gsub("V", "Cluster", colnames(qmatrix_with_data))

# For barplots we only need Ind and selected population label columns from indpopdata.

# Get color palette from config - subset first n_clusters colors
cat(sprintf("\n=== SETTING UP COLOR PALETTE ===\n"))
structure_colors_full <- unlist(params[["structure_colors"]])
cluster_cols_val <- structure_colors_full[1:n_clusters]
cat(sprintf("Using first %d colors from structure_colors palette: %s\n",
            n_clusters, paste(cluster_cols_val, collapse = ", ")))

# Parse site_order parameter
# Snakemake automatically converts Python lists to R vectors
if (is.null(site_order)) {
  site_order_val <- NULL
} else if (length(site_order) == 1 && (site_order == "NULL" || site_order == "")) {
  site_order_val <- NULL
} else if (is.vector(site_order) || is.list(site_order)) {
  # Already a vector/list from Python list in config
  site_order_val <- unlist(site_order)
} else {
  # Try to parse as R expression (backward compatibility)
  site_order_val <- tryCatch(
    eval(parse(text = site_order)),
    error = function(e) NULL
  )
}

# Parse population_sort_by: indpopdata column used to order populations when site_order is unset
if (is.null(population_sort_by)) {
  population_sort_by_val <- NULL
} else if (length(population_sort_by) == 1 &&
           (is.na(population_sort_by) || population_sort_by == "NULL" || population_sort_by == "")) {
  population_sort_by_val <- NULL
} else {
  population_sort_by_val <- as.character(unlist(population_sort_by))[1]
}

if (is.null(site_order_val) && !is.null(population_sort_by_val)) {
  if (!(population_sort_by_val %in% colnames(indpopdata))) {
    stop(sprintf(
      "ERROR: population_sort_by column not found in indpopdata: %s",
      population_sort_by_val
    ))
  }

  sort_df <- tibble(
    Site = population_label_values,
    sort_value = indpopdata[[population_sort_by_val]]
  ) %>%
    distinct(Site, .keep_all = TRUE) %>%
    arrange(sort_value, Site)

  site_order_val <- sort_df$Site
  cat(sprintf(
    "Sorting populations by column '%s': %s\n",
    population_sort_by_val,
    paste(site_order_val, collapse = ", ")
  ))
}

# Create traditional structure barplot
structure_barplot <- structure_plot(
  admixture_df = qmatrix_with_data,
  type = "structure",
  cluster_cols = cluster_cols_val,
  site_dividers = site_dividers,
  divider_width = divider_width,
  site_order = site_order_val,
  labels = "site",
  display_site_labels = display_population_labels,
  flip_axis = flip_axis,
  site_ticks_size = site_ticks_size,
  site_labels_y = site_labels_y,
  site_labels_size = site_labels_size,
  site_labels_angle = site_labels_angle
) +
  coord_cartesian(clip = "off") +  # Allow labels to extend beyond plot area
  theme(
    axis.title.y = element_text(size = 8, hjust = 1),
    axis.text.y = element_text(size = 5),
    # Remove ALL facet strip "boxes" around site names (handles different ggplot2 versions)
    # Based on mapmixture source: site labels use facets, boxes come from strip backgrounds
    strip.background = element_blank(),
    strip.background.x = element_blank(),
    strip.background.y = element_blank(),
    strip.text = element_text(),
    strip.text.x = element_text(),  # Explicitly set strip text (no background)
    strip.text.y = element_text(),  # Explicitly set strip text (no background)
    panel.background = element_blank(),
    plot.background = element_blank(),
    # Remove any annotation backgrounds
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.key = element_blank(),
    # Add extra margin at bottom for rotated labels to prevent clipping
    plot.margin = margin(
      t = 5, 
      r = 5, 
      b = ifelse(site_labels_angle >= 45, max(30, site_labels_size * 5), 5), 
      l = 5, 
      unit = "pt"
    )
  )


# Calculate plot dimensions
n_individuals <- nrow(qmatrix_with_data)
cat(sprintf("Number of individuals: %d\n", n_individuals))

# Fixed height, width scales with number of individuals
height <- 1.2  # fixed at 1.2 inches
width <- max(8, 1 + 0.01 * n_individuals)  # at least 8 inches, scaling with individuals
cat(sprintf("Plot dimensions: %.2f x %.2f inches\n", width, height))

# Save barplot
ggsave_pdf(
  filename = output_barplot,
  plot = structure_barplot,
  width = width,
  height = height,
  dpi = dpi
)

# Save barplot ggplot object as RDS (facet script reads panel layout attributes).
bottom_margin_pt <- ifelse(
  site_labels_angle >= 45,
  max(30, site_labels_size * 5),
  5
)
top_margin_pt <- 5
site_label_height_in <- bottom_margin_pt / 72
bar_height_in <- height - (bottom_margin_pt + top_margin_pt) / 72
attr(structure_barplot, "panel_width") <- width
attr(structure_barplot, "panel_height") <- height
attr(structure_barplot, "panel_bottom_margin_pt") <- bottom_margin_pt
attr(structure_barplot, "panel_bar_height_in") <- bar_height_in
attr(structure_barplot, "panel_site_label_height_in") <- site_label_height_in
saveRDS(structure_barplot, file = output_barplot_rds)

cat("Successfully created structure barplot\n")
