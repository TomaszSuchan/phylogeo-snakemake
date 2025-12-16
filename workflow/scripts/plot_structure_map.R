#!/usr/bin/env Rscript
# Script to run mapmixture for plotting admixture results on maps

library(tidyverse)
library(mapmixture)

# Prevent creation of Rplots.pdf
pdf(NULL)

# Redirect all output to log file
log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

# Snakemake inputs/outputs
qmatrix_file <- snakemake@input[["qmatrix"]]
popmap_file <- snakemake@input[["popmap"]]
indpopdata_file <- snakemake@input[["indpopdata"]]  # Read from input (generated indpopdata)
output_plot <- snakemake@output[["plot"]]
output_plot_rds <- snakemake@output[["plot_rds"]]
output_prefix <- snakemake@params[["output_prefix"]]

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

# Read Q matrix (space or tab delimited, no headers)
message("\n=== READING Q MATRIX ===\n")
qmatrix <- read.table(qmatrix_file, header = FALSE, sep = "")
n_individuals <- nrow(qmatrix)
n_clusters <- ncol(qmatrix)
message(sprintf("Q matrix file: %s\n", qmatrix_file))
message(sprintf("Q matrix dimensions: %d individuals x %d clusters\n", n_individuals, n_clusters))

# Read popmap file, first column is individual, second population
message("\n=== READING POPMAP ===\n")
popmap <- read.table(popmap_file, header = FALSE, sep = "\t")
colnames(popmap) <- c("Ind", "Site")
message(sprintf("Popmap file: %s\n", popmap_file))
message(sprintf("Original popmap contains %d individuals with %d columns\n", nrow(popmap), ncol(popmap)))
message(sprintf("All individuals in popmap:\n"))
print(popmap)
message(sprintf("Unique sites in popmap: %d\n", length(unique(popmap$Site))))
message(sprintf("Sites: %s\n", paste(unique(popmap$Site), collapse = ", ")))

# For STRUCTURE, the Q matrix should be in the same order as the input file
# Filter popmap to match the number of individuals in Q matrix
if (nrow(popmap) > nrow(qmatrix)) {
  message(sprintf("WARNING: Filtering popmap from %d to %d individuals to match Q matrix\n",
              nrow(popmap), nrow(qmatrix)))
  popmap <- popmap[1:nrow(qmatrix), ]
}

# Check dimensions match
if (nrow(qmatrix) != nrow(popmap)) {
  stop(sprintf("ERROR: Q matrix has %d individuals but popmap has %d individuals",
               nrow(qmatrix), nrow(popmap)))
}
message(sprintf("Dimensions match: %d individuals\n", nrow(popmap)))

# Combine Q matrix with individual-level data, rename colnames
message("\n=== CREATING Q MATRIX WITH METADATA ===\n")
qmatrix_with_data <- cbind(popmap$Site, popmap$Ind, qmatrix)
colnames(qmatrix_with_data)[1] <- "Site"
colnames(qmatrix_with_data)[2] <- "Ind"
colnames(qmatrix_with_data) <- gsub("V", "Cluster", colnames(qmatrix_with_data))
message(sprintf("qmatrix_with_data dimensions: %d rows x %d columns\n",
            nrow(qmatrix_with_data), ncol(qmatrix_with_data)))
message(sprintf("Columns: %s\n", paste(colnames(qmatrix_with_data), collapse = ", ")))
message(sprintf("All qmatrix_with_data rows:\n"))
print(qmatrix_with_data)

# Read indpopdata (generated file has columns: Ind, Site, Lat, Lon, ...)
# The Site column from indpopdata already represents the correct individual-location mapping
message("\n=== READING INDPOPDATA ===\n")
indpopdata <- read.table(indpopdata_file, header = TRUE, sep = "\t")
message(sprintf("indpopdata file: %s\n", indpopdata_file))
message(sprintf("indpopdata dimensions: %d rows x %d columns\n", nrow(indpopdata), ncol(indpopdata)))
message(sprintf("Columns: %s\n", paste(colnames(indpopdata), collapse = ", ")))
message(sprintf("All indpopdata rows:\n"))
print(indpopdata)
message(sprintf("Unique sites in indpopdata: %d\n", length(unique(indpopdata$Site))))

# Get unique individuals from qmatrix_with_data
message("\n=== FILTERING INDPOPDATA ===\n")
inds_in_qmatrix <- unique(qmatrix_with_data$Ind)
message(sprintf("Found %d unique individuals in qmatrix_with_data\n", length(inds_in_qmatrix)))
message(sprintf("First 10 individuals: %s\n", paste(head(inds_in_qmatrix, 10), collapse = ", ")))

# Filter indpopdata to keep only individuals that are in qmatrix_with_data
n_before_filter <- nrow(indpopdata)
indpopdata <- indpopdata %>%
  filter(Ind %in% inds_in_qmatrix)
n_after_filter <- nrow(indpopdata)

message(sprintf("After filtering, indpopdata contains %d individuals (removed %d)\n",
            n_after_filter, n_before_filter - n_after_filter))

# Check for individuals in Q matrix but not in indpopdata
inds_in_indpopdata <- unique(indpopdata$Ind)
missing_in_indpopdata <- setdiff(inds_in_qmatrix, inds_in_indpopdata)
if (length(missing_in_indpopdata) > 0) {
  message(sprintf("WARNING: %d individuals in Q matrix are missing from indpopdata:\n",
              length(missing_in_indpopdata)))
  cat(paste(head(missing_in_indpopdata, 20), collapse = ", "))
  cat("\n")
}

# The qmatrix_with_data already has the correct Site from popmap (via indpopdata)
# Just verify the Site names match between qmatrix_with_data and indpopdata
message("\n=== SITE VERIFICATION ===\n")
sites_in_qmatrix <- unique(qmatrix_with_data$Site)
sites_in_indpopdata <- unique(indpopdata$Site)
message(sprintf("qmatrix_with_data has %d unique sites\n", length(sites_in_qmatrix)))
message(sprintf("indpopdata has %d unique sites\n", length(sites_in_indpopdata)))
message(sprintf("Sites in qmatrix: %s\n", paste(sites_in_qmatrix, collapse = ", ")))
message(sprintf("Sites in indpopdata: %s\n", paste(sites_in_indpopdata, collapse = ", ")))

# Create coords_df with unique sites and their coordinates from indpopdata
message("\n=== CREATING COORDINATES DATAFRAME ===\n")
coords_df <- indpopdata %>%
  select(Site, Lat, Lon) %>%
  distinct(Site, .keep_all = TRUE)

message(sprintf("Coords dataframe has %d unique sites\n", nrow(coords_df)))
message(sprintf("Coordinate range: Lat [%.2f, %.2f], Lon [%.2f, %.2f]\n",
            min(coords_df$Lat), max(coords_df$Lat),
            min(coords_df$Lon), max(coords_df$Lon)))
message("All sites with coordinates:\n")
print(coords_df)

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
message(sprintf("Total: %d unique sites with coordinates for %d individuals\n",
            nrow(coords_df), nrow(qmatrix_with_data)))

# Summary of individuals per site
message("\n=== INDIVIDUALS PER SITE ===\n")
inds_per_site <- qmatrix_with_data %>%
  group_by(Site) %>%
  summarise(n_individuals = n(), .groups = 'drop') %>%
  arrange(desc(n_individuals))
message(sprintf("All sites with individual counts:\n"))
print(inds_per_site)
message(sprintf("\nTotal sites: %d, Total individuals: %d\n",
            nrow(inds_per_site), sum(inds_per_site$n_individuals)))

# Run mapmixture with all parameters
message("\n=== RUNNING MAPMIXTURE ===\n")
message(sprintf("Input to mapmixture:\n"))
message(sprintf("  admixture_df: %d rows x %d columns\n", nrow(qmatrix_with_data), ncol(qmatrix_with_data)))
message(sprintf("  coords_df: %d sites\n", nrow(coords_df)))
message(sprintf("  Number of clusters: %d\n", n_clusters))
message(sprintf("All rows of admixture_df being passed to mapmixture:\n"))
print(qmatrix_with_data)

p <- mapmixture(
  admixture_df = qmatrix_with_data,
  coords_df = coords_df,
  boundary = if (length(boundary) == 0 || is.null(boundary) || boundary == "NULL") NULL else eval(parse(text = boundary)),
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
message(sprintf("  - %d genetic clusters (K)\n", n_clusters))

# Close log file sinks
sink(type = "message")
sink(type = "output")
close(log_file)

