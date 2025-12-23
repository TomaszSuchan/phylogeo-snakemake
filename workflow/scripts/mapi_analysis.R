#!/usr/bin/env Rscript
# MAPI Analysis for Population Genetic Data
# This script performs MAPI analysis using Euclidean genetic distance
# and geographic coordinates to identify spatial patterns of genetic diversity

library(mapi)
library(sf)

# ==============================================================================
# Read inputs from Snakemake
# ==============================================================================

# Input files
euclidean_dist_file <- snakemake@input[['euclidean_dist']]
indpopdata_file <- snakemake@input[['indpopdata']]

# Output files
mapi_gpkg <- snakemake@output[['mapi_gpkg']]
upper_tails_gpkg <- snakemake@output[['upper_tails_gpkg']]
lower_tails_gpkg <- snakemake@output[['lower_tails_gpkg']]

# Parameters
n_permutations <- snakemake@params[['n_permutations']]
grid_halfwidth <- snakemake@params[['grid_halfwidth']]
crs_projected <- snakemake@params[['crs_projected']]
crs_geographic <- snakemake@params[['crs_geographic']]
alpha <- snakemake@params[['alpha']]

# Debug output
cat("=== MAPI Analysis Script ===\n")
cat("Euclidean distance file:", euclidean_dist_file, "\n")
cat("Individual population data file:", indpopdata_file, "\n")
cat("Projected CRS:", crs_projected, "\n")
cat("Geographic CRS:", crs_geographic, "\n")
cat("Grid halfwidth:", grid_halfwidth, "\n")
cat("Permutations:", n_permutations, "\n")
cat("Alpha:", alpha, "\n")
cat("==============================\n\n")

# ==============================================================================
# 1. Load and prepare individual-level population/location data
# ==============================================================================

cat("Loading individual population data...\n")
indpopdata <- read.table(indpopdata_file, header = TRUE, sep = "\t")

# Debug: print column names and first few rows
cat("Column names:", paste(colnames(indpopdata), collapse=", "), "\n")
cat("First 3 rows:\n")
print(head(indpopdata, 3))

# Check required columns
required_cols <- c("Ind", "Lat", "Lon")
if (!all(required_cols %in% colnames(indpopdata))) {
  stop("indpopdata must contain columns: Ind, Lat, Lon. Found: ", paste(colnames(indpopdata), collapse=", "))
}

# Check for missing coordinates
if (any(is.na(indpopdata$Lat) | is.na(indpopdata$Lon))) {
  n_missing <- sum(is.na(indpopdata$Lat) | is.na(indpopdata$Lon))
  missing_inds <- indpopdata$Ind[is.na(indpopdata$Lat) | is.na(indpopdata$Lon)]
  cat("WARNING:", n_missing, "individuals have missing coordinates:\n")
  cat(paste(head(missing_inds, 10), collapse=", "), "\n")
  if (n_missing > 10) cat("... and", n_missing - 10, "more\n")
  stop("Missing coordinates detected. Please check the indpopdata file.")
}

# Create sf object with geographic coordinates
indpopdata_sf <- st_as_sf(indpopdata, coords = c("Lon", "Lat"), crs = crs_geographic)

# Transform to projected CRS for distance calculations
indpopdata_projected <- st_transform(indpopdata_sf, crs = crs_projected)

# Extract projected coordinates
coords <- st_coordinates(indpopdata_projected)
indpopdata$x <- coords[, 1]
indpopdata$y <- coords[, 2]

cat("Individual population data loaded:", nrow(indpopdata), "individuals\n")
cat("Coordinate range - Lat:", min(indpopdata$Lat), "to", max(indpopdata$Lat), "\n")
cat("Coordinate range - Lon:", min(indpopdata$Lon), "to", max(indpopdata$Lon), "\n")
cat("Projected X range:", min(indpopdata$x), "to", max(indpopdata$x), "\n")
cat("Projected Y range:", min(indpopdata$y), "to", max(indpopdata$y), "\n")

# ==============================================================================
# 2. Load genetic distance matrix and match with individual location data
# ==============================================================================

cat("Loading Euclidean genetic distance matrix...\n")

# Load Euclidean distance (check.names=FALSE prevents R from modifying column names)
euclidean_matrix <- as.matrix(read.table(euclidean_dist_file, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE))
cat("Euclidean distance matrix loaded:", nrow(euclidean_matrix), "x", ncol(euclidean_matrix), "\n")

# Get individual IDs from distance matrix
ind_ids <- rownames(euclidean_matrix)
cat("Number of individuals in distance matrix:", length(ind_ids), "\n")
cat("First 3 individual names from distance matrix:", head(ind_ids, 3), "\n")

# Match individuals from distance matrix with location data
cat("\nMatching individuals with location data...\n")
cat("Number of individuals in indpopdata:", nrow(indpopdata), "\n")
cat("First 3 individual names from indpopdata:", head(indpopdata$Ind, 3), "\n")

# Create individual data for MAPI (ind, x, y)
inddata <- indpopdata[, c("Ind", "x", "y")]
colnames(inddata) <- c("ind", "x", "y")

# Ensure order matches the distance matrix
cat("\nReordering inddata to match distance matrix...\n")
inddata <- inddata[match(ind_ids, inddata$ind), ]

# Check for missing individuals
if (any(is.na(inddata$ind))) {
  n_missing <- sum(is.na(inddata$ind))
  missing_inds <- ind_ids[is.na(inddata$ind)]
  cat("ERROR:", n_missing, "individuals from distance matrix not found in indpopdata:\n")
  cat(paste(head(missing_inds, 10), collapse=", "), "\n")
  if (n_missing > 10) cat("... and", n_missing - 10, "more\n")
  stop("Some individuals from distance matrix not found in location data")
}

# Check for duplicates
if (any(duplicated(inddata$ind))) {
  stop("Duplicate individual IDs found in location data")
}

cat("Individual location data prepared:", nrow(inddata), "individuals with coordinates\n")
cat("Data check - any NA in ind:", any(is.na(inddata$ind)), "\n")
cat("Data check - any NA in x:", any(is.na(inddata$x)), "\n")
cat("Data check - any NA in y:", any(is.na(inddata$y)), "\n")

# ==============================================================================
# 3. Create hexagonal grid
# ==============================================================================

cat("Creating hexagonal grid...\n")

# MAPI requires a data frame with ind, x, y for grid creation
grid_input <- inddata

# Determine grid halfwidth
if (is.null(grid_halfwidth) || is.na(grid_halfwidth)) {
  cat("Using MAPI auto-estimated grid...\n")
  mapi_grid <- MAPI_GridAuto(grid_input, crs_projected)
  cat("Grid created with", nrow(mapi_grid), "cells\n")
} else {
  cat("Creating hexagonal grid with specified halfwidth:", grid_halfwidth, "meters\n")
  mapi_grid <- MAPI_GridHexagonal(grid_input, crs_projected, grid_halfwidth)
  cat("Grid created with", nrow(mapi_grid), "cells\n")
}

# ==============================================================================
# 4. Convert distance matrix to long format for MAPI
# ==============================================================================

cat("Converting distance matrix to long format...\n")

# Convert Euclidean distance to long format
euclidean_long <- data.frame(
  ind1 = rep(rownames(euclidean_matrix), each = ncol(euclidean_matrix)),
  ind2 = rep(colnames(euclidean_matrix), times = nrow(euclidean_matrix)),
  value = as.vector(euclidean_matrix),
  stringsAsFactors = FALSE
)
# Remove NA values if any
euclidean_long <- euclidean_long[!is.na(euclidean_long$value), ]
cat("Euclidean long format:", nrow(euclidean_long), "pairwise distances\n")

cat("Distance matrix converted to long format\n")

# ==============================================================================
# 5. Run MAPI with Euclidean distance
# ==============================================================================

cat("\n=== Running MAPI with Euclidean distance ===\n")
cat("Permutations:", n_permutations, "\n")

mapi_results_euclidean <- MAPI_RunOnGrid(
  inddata,
  euclidean_long,
  mapi_grid,
  nbPermuts = n_permutations,
  nbCores = snakemake@threads
)

# Calculate significance tails
mapi_tails_euclidean <- MAPI_Tails(mapi_results_euclidean, alpha = alpha)

# Extract upper and lower tails
upper_tails_euclidean <- st_intersection(
  mapi_results_euclidean,
  mapi_tails_euclidean[mapi_tails_euclidean$tail == "upper", ]
)
upper_tails_euclidean <- upper_tails_euclidean[, colnames(mapi_results_euclidean)]

lower_tails_euclidean <- st_intersection(
  mapi_results_euclidean,
  mapi_tails_euclidean[mapi_tails_euclidean$tail == "lower", ]
)
lower_tails_euclidean <- lower_tails_euclidean[, colnames(mapi_results_euclidean)]

cat("MAPI analysis with Euclidean distance complete\n")
cat("Significant upper tail cells:", nrow(upper_tails_euclidean), "\n")
cat("Significant lower tail cells:", nrow(lower_tails_euclidean), "\n")

# ==============================================================================
# 6. Save results as GeoPackage files
# ==============================================================================

cat("\nSaving results to GeoPackage files...\n")

# Main results file
st_write(mapi_results_euclidean, dsn = mapi_gpkg, layer = "euclidean_results",
         driver = "GPKG", append = FALSE, quiet = TRUE)
st_write(mapi_tails_euclidean, dsn = mapi_gpkg, layer = "euclidean_tails",
         driver = "GPKG", append = TRUE, quiet = TRUE)

cat("Main results saved to:", mapi_gpkg, "\n")

# Upper tails
st_write(upper_tails_euclidean, dsn = upper_tails_gpkg, layer = "euclidean_upper",
         driver = "GPKG", append = FALSE, quiet = TRUE)

cat("Upper tails saved to:", upper_tails_gpkg, "\n")

# Lower tails
st_write(lower_tails_euclidean, dsn = lower_tails_gpkg, layer = "euclidean_lower",
         driver = "GPKG", append = FALSE, quiet = TRUE)

cat("Lower tails saved to:", lower_tails_gpkg, "\n")

# ==============================================================================
# 7. Summary statistics
# ==============================================================================

cat("\n=== MAPI Analysis Summary ===\n")
cat("\nEuclidean Distance:\n")
cat("  Mean avg_value:", mean(mapi_results_euclidean$avg_value, na.rm = TRUE), "\n")
cat("  SD avg_value:", sd(mapi_results_euclidean$avg_value, na.rm = TRUE), "\n")
cat("  Significant upper tails:", nrow(upper_tails_euclidean), "cells\n")
cat("  Significant lower tails:", nrow(lower_tails_euclidean), "cells\n")

cat("\nAnalysis complete!\n")
