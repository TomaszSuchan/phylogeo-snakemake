#!/usr/bin/env Rscript
# conStruct Analysis for Spatial Population Structure
# This script performs conStruct analysis using VCF data and geographic coordinates
# to infer spatial population structure with geography

# ==============================================================================
# Read inputs from Snakemake and redirect all output to log file
# ==============================================================================

# Get log file path from Snakemake
log_file_path <- snakemake@log[[1]]

# Redirect all output (stdout and stderr) to the log file
# This ensures all output from R, packages, and external processes (like Stan) goes to the log
sink(log_file_path, type = "output")
sink(log_file_path, append = TRUE, type = "message")

# Install conStruct if not available (from CRAN)
if (!require("conStruct", quietly = TRUE)) {
    cat("conStruct package not found. Installing from CRAN...\n")
    install.packages("conStruct", repos = "https://cloud.r-project.org", quiet = TRUE)
    library(conStruct)
}

library(vcfR)
library(adegenet)
library(geosphere)

# Input files
vcf_file <- snakemake@input[['vcf']]
indpopdata_file <- snakemake@input[['indpopdata']]

# Output files
results_rds <- snakemake@output[['results_rds']]
layer_proportions <- snakemake@output[['layer_proportions']]
log_file <- snakemake@output[['log_file']]

# Parameters
k <- as.integer(snakemake@params[['k']])
n_chains <- as.integer(snakemake@params[['n_chains']])
n_iterations <- as.integer(snakemake@params[['n_iterations']])
make.freqs <- as.logical(snakemake@params[['make_freqs']])
geoDist_param <- snakemake@params[['geoDist']]
coords_param <- snakemake@params[['coords']]
save.files <- as.logical(snakemake@params[['save_files']])

# Handle NULL/None values from Python (convert to R NULL)
if (is.null(geoDist_param) || geoDist_param == "None" || geoDist_param == "NULL") {
  geoDist_param <- NULL
}
if (is.null(coords_param) || coords_param == "None" || coords_param == "NULL") {
  coords_param <- NULL
}

# Debug output
cat("=== conStruct Analysis Script ===\n")
cat("VCF file:", vcf_file, "\n")
cat("Individual population data file:", indpopdata_file, "\n")
cat("K (number of layers):", k, "\n")
cat("Number of chains:", n_chains, "\n")
cat("Number of iterations:", n_iterations, "\n")
cat("================================\n\n")

# ==============================================================================
# 1. Load and prepare individual-level population/location data
# ==============================================================================

cat("Loading individual population data...\n")
indpopdata <- read.table(indpopdata_file, header = TRUE, sep = "\t")

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

# Extract coordinates
coords <- as.matrix(indpopdata[, c("Lon", "Lat")])
rownames(coords) <- indpopdata$Ind

cat("Individual population data loaded:", nrow(indpopdata), "individuals\n")
cat("Coordinate range - Lat:", min(indpopdata$Lat), "to", max(indpopdata$Lat), "\n")
cat("Coordinate range - Lon:", min(indpopdata$Lon), "to", max(indpopdata$Lon), "\n")

# ==============================================================================
# 2. Load VCF and convert to allele frequency format
# ==============================================================================

cat("Reading VCF file...\n")
vcf <- read.vcfR(vcf_file)

cat("Converting VCF to genind object...\n")
genind_obj <- vcfR2genind(vcf)

# Get individual names from genind
ind_names <- indNames(genind_obj)
cat("Number of individuals in VCF:", length(ind_names), "\n")
cat("First 3 individual names from VCF:", head(ind_names, 3), "\n")

# Match individuals from VCF with location data
cat("\nMatching individuals with location data...\n")
cat("Number of individuals in indpopdata:", nrow(indpopdata), "\n")
cat("First 3 individual names from indpopdata:", head(indpopdata$Ind, 3), "\n")

# Ensure order matches between VCF and coordinates
if (!all(ind_names %in% indpopdata$Ind)) {
  missing_inds <- ind_names[!ind_names %in% indpopdata$Ind]
  cat("ERROR: Individuals in VCF not found in indpopdata:\n")
  cat(paste(head(missing_inds, 10), collapse=", "), "\n")
  if (length(missing_inds) > 10) cat("... and", length(missing_inds) - 10, "more\n")
  stop("Some individuals from VCF not found in location data")
}

# Reorder coordinates to match VCF order
coords_ordered <- coords[ind_names, , drop = FALSE]

# Check for missing coordinates after reordering
if (any(is.na(coords_ordered))) {
  stop("Missing coordinates after matching with VCF individuals")
}

cat("Coordinates matched with VCF individuals\n")

# ==============================================================================
# 3. Calculate allele frequencies
# ==============================================================================

cat("Calculating allele frequencies...\n")

# Extract genotype matrix (allele counts: 0, 1, 2)
gen_matrix <- tab(genind_obj, NA.method = "mean")

# Convert to allele frequency format for conStruct
# conStruct expects allele frequencies (proportion of allele 1)
# For each locus, we need the frequency of the alternate allele

# Get number of loci and individuals
n_loci <- ncol(gen_matrix)
n_inds <- nrow(gen_matrix)

cat("Number of loci:", n_loci, "\n")
cat("Number of individuals:", n_inds, "\n")

# Convert genotype matrix (0, 1, 2) to allele frequencies
# Each value represents the count of alternate alleles (0=hom ref, 1=het, 2=hom alt)
# Divide by 2 to get frequency of alternate allele
allele_freqs <- gen_matrix / 2

# Ensure values are between 0 and 1
allele_freqs[allele_freqs < 0] <- 0
allele_freqs[allele_freqs > 1] <- 1

cat("Allele frequency matrix prepared:", nrow(allele_freqs), "individuals x", ncol(allele_freqs), "loci\n")

# ==============================================================================
# 4. Calculate geographic distance matrix
# ==============================================================================

cat("Calculating geographic distance matrix...\n")

# Calculate pairwise geographic distances in meters
geoDist <- distm(coords_ordered, fun = distHaversine)

cat("Geographic distance matrix calculated:", nrow(geoDist), "x", ncol(geoDist), "\n")
cat("Distance range:", min(geoDist[geoDist > 0]), "to", max(geoDist), "meters\n")

# ==============================================================================
# 5. Prepare data for conStruct
# ==============================================================================

cat("\n=== Preparing data for conStruct ===\n")

# conStruct requires:
# - allele.frequencies: matrix of allele frequencies (n_inds x n_loci)
# - coords: matrix of coordinates (n_inds x 2) with columns [longitude, latitude]
# - geoDist: matrix of pairwise geographic distances (optional but recommended)

# Ensure allele frequencies are in the correct format
# conStruct expects frequencies between 0 and 1
allele.frequencies <- as.matrix(allele_freqs)

# Ensure coordinates are in correct format [longitude, latitude]
coords_matrix <- coords_ordered

cat("Data prepared for conStruct:\n")
cat("  Allele frequencies:", nrow(allele.frequencies), "individuals x", ncol(allele.frequencies), "loci\n")
cat("  Coordinates:", nrow(coords_matrix), "individuals x", ncol(coords_matrix), "dimensions\n")
cat("  Geographic distances:", nrow(geoDist), "x", ncol(geoDist), "\n")

# ==============================================================================
# 6. Run conStruct analysis
# ==============================================================================

cat("\n=== Running conStruct analysis ===\n")
cat("K (number of layers):", k, "\n")
cat("Number of chains:", n_chains, "\n")
cat("Number of iterations:", n_iterations, "\n")
cat("Spatial model: TRUE\n")

# Create output directory if it doesn't exist
output_dir <- dirname(results_rds)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Run conStruct with spatial model
# conStruct can run in spatial or non-spatial mode
# We use spatial mode to incorporate geography
cat("Starting MCMC chains...\n")

# Run conStruct analysis
# conStruct function signature: conStruct(spatial, K, freqs, coords, geoDist, n.iter, ...)
# Note: conStruct writes files to disk by default
cat("Running conStruct with K =", k, "layers...\n")

# Use provided coords/geoDist if available, otherwise use computed ones
final_coords <- if (!is.null(coords_param)) coords_param else coords_matrix
final_geoDist <- if (!is.null(geoDist_param)) geoDist_param else geoDist

# Set up prefix for output files
output_prefix <- file.path(output_dir, paste0(basename(tools::file_path_sans_ext(results_rds)), ".K", k))

# Run conStruct analysis
# conStruct returns a list with results
# Note: make.freqs and save.files are not valid parameters for conStruct
# We already provide allele frequencies directly, so make.freqs is not needed
construct_results <- conStruct::conStruct(
  spatial = TRUE,
  K = k,
  freqs = allele.frequencies,
  coords = final_coords,
  geoDist = final_geoDist,
  n.iter = n_iterations,
  prefix = output_prefix
)

cat("conStruct analysis complete\n")

# ==============================================================================
# 7. Extract layer proportions
# ==============================================================================

cat("Extracting layer proportions...\n")

# Extract layer proportions (admixture proportions) from results
# conStruct results contain layer proportions in the "layer.proportions" component
if ("layer.proportions" %in% names(construct_results)) {
  layer_props <- construct_results$layer.proportions
} else if ("results" %in% names(construct_results) && "layer.proportions" %in% names(construct_results$results)) {
  layer_props <- construct_results$results$layer.proportions
} else {
  # Try to extract from the structure of results
  # conStruct may store results differently depending on version
  layer_props <- construct_results[[1]]$layer.proportions
}

# If layer_props is a list (multiple chains), take the first one or average
if (is.list(layer_props) && !is.data.frame(layer_props) && !is.matrix(layer_props)) {
  layer_props <- layer_props[[1]]
}

# Ensure it's a matrix or data frame
if (!is.matrix(layer_props) && !is.data.frame(layer_props)) {
  cat("WARNING: Could not extract layer proportions in expected format\n")
  cat("Available names in results:", paste(names(construct_results), collapse=", "), "\n")
  layer_props <- matrix(NA, nrow = n_inds, ncol = k)
  rownames(layer_props) <- ind_names
  colnames(layer_props) <- paste0("Layer", 1:k)
} else {
  # Add row names if missing
  if (is.null(rownames(layer_props))) {
    rownames(layer_props) <- ind_names[1:nrow(layer_props)]
  }
  # Add column names if missing
  if (is.null(colnames(layer_props))) {
    colnames(layer_props) <- paste0("Layer", 1:ncol(layer_props))
  }
}

# Write layer proportions to file
write.table(
  layer_props,
  file = layer_proportions,
  sep = "\t",
  col.names = NA,
  row.names = TRUE,
  quote = FALSE
)

cat("Layer proportions saved to:", layer_proportions, "\n")

# ==============================================================================
# 8. Save results as RDS
# ==============================================================================

cat("Saving results to RDS file...\n")
saveRDS(construct_results, file = results_rds)
cat("Results saved to:", results_rds, "\n")

# ==============================================================================
# 9. Write summary log
# ==============================================================================

cat("Writing summary log...\n")
log_conn <- file(log_file, "w")
cat("=== conStruct Analysis Summary ===\n", file = log_conn)
cat("K (number of layers):", k, "\n", file = log_conn)
cat("Number of chains:", n_chains, "\n", file = log_conn)
cat("Number of iterations:", n_iterations, "\n", file = log_conn)
cat("Number of individuals:", n_inds, "\n", file = log_conn)
cat("Number of loci:", n_loci, "\n", file = log_conn)
cat("Spatial model: TRUE\n", file = log_conn)
cat("\nLayer proportions extracted:", nrow(layer_props), "individuals x", ncol(layer_props), "layers\n", file = log_conn)
close(log_conn)

cat("\n=== conStruct Analysis Complete ===\n")
cat("Results saved to:", results_rds, "\n")
cat("Layer proportions saved to:", layer_proportions, "\n")
cat("Log file saved to:", log_file, "\n")

# Close sink to restore normal output
sink(type = "output")
sink(type = "message")

