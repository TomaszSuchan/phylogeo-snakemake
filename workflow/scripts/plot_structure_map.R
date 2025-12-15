#!/usr/bin/env Rscript
# Script to run mapmixture for plotting admixture results on maps

library(tidyverse)
library(mapmixture)

# Prevent creation of Rplots.pdf
pdf(NULL)

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
qmatrix <- read.table(qmatrix_file, header = FALSE, sep = "")
n_individuals <- nrow(qmatrix)
n_clusters <- ncol(qmatrix)
cat(sprintf("Q matrix dimensions: %d individuals x %d clusters\n", n_individuals, n_clusters))

# Read popmap file, first column is individual, second population
popmap <- read.table(popmap_file, header = FALSE, sep = "\t")
colnames(popmap) <- c("Ind", "Site")
cat(sprintf("Original popmap contains %d individuals with %d columns\n", nrow(popmap), ncol(popmap)))

# Get individual IDs from Q matrix file (first column if it exists in the file structure)
# Since Q matrix file might not have individual names, we need to read from a .indv or similar file
# For STRUCTURE output, individual order should match the input order
# We'll use the first n_individuals from popmap that match the analysis
# Actually, we need to ensure only individuals in the Q matrix are used

# For STRUCTURE, the Q matrix should be in the same order as the input file
# Filter popmap to match the number of individuals in Q matrix
if (nrow(popmap) > nrow(qmatrix)) {
  cat(sprintf("Filtering popmap from %d to %d individuals to match Q matrix\n", nrow(popmap), nrow(qmatrix)))
  popmap <- popmap[1:nrow(qmatrix), ]
}

# Check dimensions match
if (nrow(qmatrix) != nrow(popmap)) {
  stop(sprintf("Q matrix has %d individuals but popmap has %d individuals",
               nrow(qmatrix), nrow(popmap)))
}

# Combine Q matrix with individual-level data, rename colnames
qmatrix_with_data <- cbind(popmap$Site, popmap$Ind, qmatrix)
colnames(qmatrix_with_data)[1] <- "Site"
colnames(qmatrix_with_data)[2] <- "Ind"
colnames(qmatrix_with_data) <- gsub("V", "Cluster", colnames(qmatrix_with_data))

# Read indpopdata (generated file has columns: Ind, Site, Lat, Lon, ...)
# Select only the columns we need: Ind, Lat, Lon
indpopdata <- read.table(indpopdata_file, header = TRUE, sep = "\t")
popdata <- indpopdata[, c("Ind", "Lat", "Lon")]

# Get unique individuals from qmatrix_with_data
inds_in_qmatrix <- unique(qmatrix_with_data$Ind)
cat(sprintf("Found %d unique individuals in qmatrix_with_data\n", length(inds_in_qmatrix)))

# Filter popdata to keep only individuals that are in qmatrix_with_data
popdata <- popdata %>%
  filter(Ind %in% inds_in_qmatrix)

cat(sprintf("After filtering, popdata contains %d individuals\n", nrow(popdata)))

# Create Site names based on coordinates to group individuals at same location
# This allows mapmixture to aggregate individuals at the same coordinates into one pie chart
popdata <- popdata %>%
  mutate(Site_coord = paste0("Loc_", round(Lat, 6), "_", round(Lon, 6)))

# Update qmatrix_with_data to use the same coordinate-based Site names
# The left_join will create Site.x (original) and Site_coord, then we clean up
qmatrix_with_data <- qmatrix_with_data %>%
  left_join(popdata %>% select(Ind, Site_coord), by = "Ind") %>%
  select(-Site) %>%  # Remove original Site column
  rename(Site = Site_coord) %>%
  select(Site, Ind, everything())

# Also rename Site_coord to Site in popdata for consistency
popdata <- popdata %>%
  rename(Site = Site_coord)

cat(sprintf("Created %d unique location-based sites from %d individuals\n",
            length(unique(popdata$Site)), nrow(popdata)))

# Now create coords_df with unique sites and their coordinates
coords_df <- popdata %>%
  select(Site, Lat, Lon) %>%
  distinct(Site, .keep_all = TRUE)

cat(sprintf("Coords dataframe has %d unique sites\n", nrow(coords_df)))

# Validate that all sites from qmatrix_with_data have coordinates
sites_in_qmatrix <- unique(qmatrix_with_data$Site)
sites_in_coords <- unique(coords_df$Site)
missing_sites <- setdiff(sites_in_qmatrix, sites_in_coords)

if (length(missing_sites) > 0) {
  cat(sprintf("ERROR: %d sites from qmatrix_with_data are missing in coords_df:\n", length(missing_sites)))
  cat(paste(missing_sites, collapse = "\n"))
  cat("\n")
  stop("Sites in Q matrix are missing from coordinates file")
}

cat(sprintf("\nValidation passed: %d unique sites with coordinates for %d individuals\n",
            nrow(coords_df), nrow(qmatrix_with_data)))

# Run mapmixture with all parameters
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

