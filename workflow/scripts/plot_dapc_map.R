#!/usr/bin/env Rscript
# Script to plot DAPC results on maps using mapmixture

library(tidyverse)
library(mapmixture)

# Prevent creation of Rplots.pdf
pdf(NULL)

# Redirect all output to log file
log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

# Snakemake inputs/outputs
dapc_results_rds <- snakemake@input[["dapc_results"]]
indpopdata_file <- snakemake@input[["indpopdata"]]
output_plot <- snakemake@output[["plot"]]
output_plot_rds <- snakemake@output[["plot_rds"]]

# Check if indpopdata file exists and is not NULL
if (is.null(indpopdata_file) || indpopdata_file == "NULL" || indpopdata_file == "" || !file.exists(indpopdata_file)) {
  stop("ERROR: indpopdata file not found or not specified. Map plotting requires a popdata file with geographic coordinates. Please specify 'popdata' in your config.yaml under parameters.")
}

# Plot dimension parameters
width <- as.numeric(snakemake@params[["width"]])
height <- as.numeric(snakemake@params[["height"]])
dpi <- as.numeric(snakemake@params[["dpi"]])

# Map and coordinate parameters
boundary <- snakemake@params[["boundary"]]
crs <- as.numeric(snakemake@params[["crs"]])
basemap <- snakemake@params[["basemap"]]

# Pie chart parameters
pie_size <- as.numeric(snakemake@params[["pie_size"]])
pie_border <- as.numeric(snakemake@params[["pie_border"]])
pie_border_col <- snakemake@params[["pie_border_col"]]
pie_opacity <- as.numeric(snakemake@params[["pie_opacity"]])

# Color parameters
land_colour <- snakemake@params[["land_colour"]]
sea_colour <- snakemake@params[["sea_colour"]]

# Map display parameters
expand <- as.logical(snakemake@params[["expand"]])
arrow <- as.logical(snakemake@params[["arrow"]])
arrow_size <- as.numeric(snakemake@params[["arrow_size"]])
arrow_position <- snakemake@params[["arrow_position"]]
scalebar <- as.logical(snakemake@params[["scalebar"]])
scalebar_size <- as.numeric(snakemake@params[["scalebar_size"]])
scalebar_position <- snakemake@params[["scalebar_position"]]

# Text parameters
plot_title <- snakemake@params[["plot_title"]]
axis_title_size <- as.numeric(snakemake@params[["axis_title_size"]])
axis_text_size <- as.numeric(snakemake@params[["axis_text_size"]])

# Border parameters
basemap_border <- as.logical(snakemake@params[["basemap_border"]])
basemap_border_col <- snakemake@params[["basemap_border_col"]]
basemap_border_lwd <- as.numeric(snakemake@params[["basemap_border_lwd"]])

# Legend parameter
legend <- as.logical(snakemake@params[["legend"]])

# Read DAPC results
message("\n=== READING DAPC RESULTS ===\n")
dapc_results <- readRDS(dapc_results_rds)
membership_probs <- dapc_results$membership_probs
k <- dapc_results$k  # Get K from results

message(sprintf("DAPC results file: %s\n", dapc_results_rds))
message(sprintf("K (number of clusters): %d\n", k))
message(sprintf("Membership probabilities dimensions: %d individuals x %d clusters\n", 
            nrow(membership_probs), ncol(membership_probs) - 1))  # -1 for Ind column

# Extract Q matrix (membership probabilities) - exclude the Ind column
qmatrix <- membership_probs[, !colnames(membership_probs) %in% "Ind"]
n_individuals <- nrow(qmatrix)
n_clusters <- ncol(qmatrix)

message(sprintf("Q matrix dimensions: %d individuals x %d clusters\n", n_individuals, n_clusters))

# Verify k matches number of clusters
if (k != n_clusters) {
  warning(sprintf("K (%d) does not match number of clusters in membership probabilities (%d). Using %d clusters.\n",
                  k, n_clusters, n_clusters))
}

# Get color palette from config - subset first n_clusters colors
message(sprintf("\n=== SETTING UP COLOR PALETTE ===\n"))
structure_colors_full <- unlist(snakemake@params[["structure_colors"]])
strcolors <- structure_colors_full[1:n_clusters]
message(sprintf("Using first %d colors from palette: %s\n", n_clusters, paste(strcolors, collapse = ", ")))

# Read indpopdata (generated file has columns: Ind, Site, Lat, Lon, ...)
message("\n=== READING INDPOPDATA ===\n")
indpopdata <- read.table(indpopdata_file, header = TRUE, sep = "\t")
message(sprintf("indpopdata file: %s\n", indpopdata_file))
message(sprintf("indpopdata dimensions: %d rows x %d columns\n", nrow(indpopdata), ncol(indpopdata)))
message(sprintf("Columns: %s\n", paste(colnames(indpopdata), collapse = ", ")))

# Match membership probabilities with Site assignments from indpopdata
message("\n=== MATCHING INDIVIDUALS ===\n")
site_assignments <- indpopdata %>%
  select(Ind, Site) %>%
  distinct(Ind, .keep_all = TRUE)

qmatrix_with_data <- merge(membership_probs, site_assignments, by = "Ind", all.x = TRUE)

if (any(is.na(qmatrix_with_data$Site))) {
  n_missing_site <- sum(is.na(qmatrix_with_data$Site))
  warning(sprintf("%d individuals do not have Site assignments and will be excluded\n", n_missing_site))
  qmatrix_with_data <- qmatrix_with_data[!is.na(qmatrix_with_data$Site), ]
}

cluster_cols <- colnames(qmatrix)
qmatrix_with_data <- qmatrix_with_data[, c("Site", "Ind", cluster_cols)]
colnames(qmatrix_with_data)[3:ncol(qmatrix_with_data)] <- paste0("Cluster", 1:n_clusters)

message(sprintf("qmatrix_with_data dimensions: %d rows x %d columns\n",
            nrow(qmatrix_with_data), ncol(qmatrix_with_data)))
message(sprintf("Columns: %s\n", paste(colnames(qmatrix_with_data), collapse = ", ")))

# Filter indpopdata to keep only individuals that are in qmatrix_with_data
message("\n=== FILTERING INDPOPDATA ===\n")
inds_in_qmatrix <- unique(qmatrix_with_data$Ind)
message(sprintf("Found %d unique individuals in qmatrix_with_data\n", length(inds_in_qmatrix)))

n_before_filter <- nrow(indpopdata)
indpopdata <- indpopdata %>%
  filter(Ind %in% inds_in_qmatrix)
n_after_filter <- nrow(indpopdata)

message(sprintf("After filtering, indpopdata contains %d individuals (removed %d)\n",
            n_after_filter, n_before_filter - n_after_filter))

# Create coords_df with unique sites and their coordinates from indpopdata
message("\n=== CREATING COORDINATES DATAFRAME ===\n")
coords_df <- indpopdata %>%
  select(Site, Lat, Lon) %>%
  distinct(Site, .keep_all = TRUE)

message(sprintf("Coords dataframe has %d unique sites\n", nrow(coords_df)))
message(sprintf("Coordinate range: Lat [%.2f, %.2f], Lon [%.2f, %.2f]\n",
            min(coords_df$Lat), max(coords_df$Lat),
            min(coords_df$Lon), max(coords_df$Lon)))

# Validate that all sites from qmatrix_with_data have coordinates
message("\n=== VALIDATING SITE COVERAGE ===\n")
sites_in_qmatrix <- unique(qmatrix_with_data$Site)
sites_in_coords <- unique(coords_df$Site)
missing_sites <- setdiff(sites_in_qmatrix, sites_in_coords)

if (length(missing_sites) > 0) {
  message(sprintf("ERROR: %d sites from qmatrix_with_data are missing in coords_df:\n",
              length(missing_sites)))
  cat(paste(missing_sites, collapse = ", "))
  cat("\n")
  stop("Sites in Q matrix are missing from coordinates file")
}

message(sprintf("Validation passed: All %d sites have coordinates\n", length(sites_in_qmatrix)))

# Run mapmixture with all parameters
message("\n=== RUNNING MAPMIXTURE ===\n")
message(sprintf("Input to mapmixture:\n"))
message(sprintf("  admixture_df: %d rows x %d columns\n", nrow(qmatrix_with_data), ncol(qmatrix_with_data)))
message(sprintf("  coords_df: %d sites\n", nrow(coords_df)))
message(sprintf("  Number of clusters: %d (optimal K from DAPC)\n", n_clusters))

# Parse boundary with error handling
boundary_parsed <- NULL
if (!(length(boundary) == 0 || is.null(boundary) || boundary == "NULL")) {
  tryCatch({
    boundary_parsed <- eval(parse(text = boundary))
  }, error = function(e) {
    stop(sprintf("Failed to parse boundary parameter '%s': %s", boundary, as.character(e)))
  })
}

# Set default title if not provided
# Always use empty title for DAPC maps
plot_title <- ""

p <- mapmixture(
  admixture_df = qmatrix_with_data,
  coords_df = coords_df,
  cluster_cols = strcolors,
  boundary = boundary_parsed,
  crs = crs,
  basemap = if (length(basemap) == 0 || is.null(basemap) || basemap == "NULL") NULL else basemap,
  pie_size = pie_size,
  pie_border = pie_border,
  pie_border_col = pie_border_col,
  pie_opacity = pie_opacity,
  land_colour = land_colour,
  sea_colour = sea_colour,
  expand = expand,
  arrow = arrow,
  arrow_size = arrow_size,
  arrow_position = arrow_position,
  scalebar = scalebar,
  scalebar_size = scalebar_size,
  scalebar_position = scalebar_position,
  plot_title = plot_title,
  axis_title_size = axis_title_size,
  axis_text_size = axis_text_size,
  basemap_border = basemap_border,
  basemap_border_col = basemap_border_col,
  basemap_border_lwd = basemap_border_lwd
)

# Remove legend if legend parameter is FALSE
if (!legend) {
  p <- p + theme(legend.position = "none")
}

message("\n=== SAVING OUTPUTS ===\n")
message(sprintf("Output PDF: %s\n", output_plot))
message(sprintf("Output RDS: %s\n", output_plot_rds))
message(sprintf("Plot dimensions: %.1f x %.1f inches @ %d dpi\n", width, height, dpi))

# Save the plot
dir.create(dirname(output_plot), recursive = TRUE, showWarnings = FALSE)
ggsave(
  filename = output_plot,
  plot = p,
  width = width,
  height = height,
  dpi = dpi,
  device = "pdf"
)

# Save the ggplot object as RDS
saveRDS(p, file = output_plot_rds)

message("\n=== COMPLETED SUCCESSFULLY ===\n")
message(sprintf("Final map includes:\n"))
message(sprintf("  - %d sites\n", nrow(coords_df)))
message(sprintf("  - %d individuals\n", nrow(qmatrix_with_data)))
message(sprintf("  - %d genetic clusters (K=%d)\n", n_clusters, k))

# Close log file sinks
sink(type = "message")
sink(type = "output")
close(log_file)

