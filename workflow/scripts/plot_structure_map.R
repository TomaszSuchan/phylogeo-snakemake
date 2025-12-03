#!/usr/bin/env Rscript
# Script to run mapmixture for plotting admixture results on maps

library(tidyverse)
library(mapmixture)

# Prevent creation of Rplots.pdf
pdf(NULL)

# Snakemake inputs/outputs
qmatrix_file <- snakemake@input[["qmatrix"]]
popmap_file <- snakemake@input[["popmap"]]
popdata_file <- snakemake@params[["popdata"]]
output_plot <- snakemake@output[["plot"]]
output_plot_rds <- snakemake@output[["plot_rds"]]
output_prefix <- snakemake@params[["output_prefix"]]

# Check if popdata file exists and is not NULL
if (is.null(popdata_file) || popdata_file == "NULL" || popdata_file == "" || !file.exists(popdata_file)) {
  stop("ERROR: popdata file not found or not specified. Map plotting requires a popdata file with geographic coordinates. Please specify 'popdata' in your config.yaml under parameters.")
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

# Read popdata
popdata <- read.table(popdata_file, header = TRUE, sep = "\t")[,c(1,2,3)]

# Get unique sites from qmatrix_with_data
sites_in_qmatrix <- unique(qmatrix_with_data$Site)
cat(sprintf("Found %d unique sites in qmatrix_with_data\n", length(sites_in_qmatrix)))

# Filter popdata to keep only sites that are in qmatrix_with_data
popdata <- popdata %>%
  filter(Site %in% sites_in_qmatrix)

cat(sprintf("After filtering, popdata contains %d sites\n", nrow(popdata)))

# Check whether all sites from qmatrix_with_data have data in popdata
sites_in_popdata <- unique(popdata$Site)
missing_sites <- setdiff(sites_in_qmatrix, sites_in_popdata)

if (length(missing_sites) > 0) {
  warning(sprintf("Warning: %d sites from qmatrix_with_data are missing in popdata: %s",
                  length(missing_sites), paste(missing_sites, collapse = ", ")))
} else {
  cat("All sites from qmatrix_with_data have corresponding data in popdata\n")
}

# Run mapmixture with all parameters
p <- mapmixture(
  admixture_df = qmatrix_with_data,
  coords_df = popdata,
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

