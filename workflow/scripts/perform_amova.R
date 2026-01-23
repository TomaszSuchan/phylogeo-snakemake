#!/usr/bin/env Rscript
# Perform AMOVA (Analysis of Molecular Variance) analysis
# Works with VCF files and configurable hierarchical stratification
# Missing data is handled by mean imputation

library(vcfR)
library(adegenet)
library(poppr)
library(ggplot2)

# Get inputs from snakemake
vcf_file <- snakemake@input[["vcf"]]
popdata_file <- snakemake@input[["popdata"]]
output_amova <- snakemake@output[["amova"]]
output_plot <- snakemake@output[["plot"]]
output_plot_rds <- snakemake@output[["plot_rds"]]

# Get parameters
strata_list <- snakemake@params[["strata"]]
nperm <- snakemake@params[["nperm"]]

cat("Reading VCF file...\n")
vcf <- read.vcfR(vcf_file)

cat("Converting VCF to genind object...\n")
genind_obj <- vcfR2genind(vcf)

# Get individual names
ind_names <- indNames(genind_obj)
cat("Number of individuals:", length(ind_names), "\n")

# Track initial number of loci
n_loci_before <- nLoc(genind_obj)
cat("Initial number of loci:", n_loci_before, "\n")

# Read population data
cat("Reading population data...\n")
popdata <- read.table(popdata_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Match individuals with population data
if (!all(ind_names %in% popdata$Ind)) {
  missing_inds <- ind_names[!ind_names %in% popdata$Ind]
  cat("WARNING: Some individuals in VCF not found in popdata:\n")
  cat(paste(head(missing_inds, 10), collapse = ", "), "\n")
  if (length(missing_inds) > 10) cat("... and", length(missing_inds) - 10, "more\n")
  
  # Filter genind to only include individuals in popdata
  keep_inds <- ind_names[ind_names %in% popdata$Ind]
  genind_obj <- genind_obj[keep_inds, ]
  ind_names <- indNames(genind_obj)
  cat("Filtered to", length(ind_names), "individuals present in popdata\n")
}

# Read population data and match individuals
popdata_matched <- popdata[match(ind_names, popdata$Ind), ]

cat("Setting up hierarchical stratification...\n")

# Validate that strata_list is provided and not empty
if (is.null(strata_list) || length(strata_list) == 0) {
  stop("AMOVA requires at least one stratification level specified in config (e.g., ['Site'] or ['Region', 'Site'])")
}

# The last element in strata_list is the population level (most nested)
population_level <- strata_list[length(strata_list)]
higher_levels <- strata_list[-length(strata_list)]

cat("Population level:", population_level, "\n")
if (length(higher_levels) > 0) {
  cat("Higher levels:", paste(higher_levels, collapse = " > "), "\n")
}

# Validate that all strata columns exist in popdata
all_strata_cols <- strata_list
missing_cols <- all_strata_cols[!all_strata_cols %in% colnames(popdata_matched)]
if (length(missing_cols) > 0) {
  stop(paste("Stratification column(s) not found in popdata:", paste(missing_cols, collapse = ", "),
             "\nAvailable columns:", paste(colnames(popdata_matched), collapse = ", ")))
}

# Create strata data frame with all levels
# Initialize with the first column to set the number of rows
strata_df <- data.frame(
  temp = popdata_matched[[strata_list[1]]],
  stringsAsFactors = FALSE
)
colnames(strata_df)[1] <- strata_list[1]

# Add remaining stratification levels to strata_df
if (length(strata_list) > 1) {
  for (i in 2:length(strata_list)) {
    stratum_col <- strata_list[i]
    strata_df[[stratum_col]] <- popdata_matched[[stratum_col]]
  }
}

# Set the population level in genind object (required by adegenet)
pop(genind_obj) <- popdata_matched[[population_level]]
cat("Number of populations:", length(unique(pop(genind_obj))), "\n")

# Set strata in genind object
strata(genind_obj) <- strata_df

# Build hierarchical formula
# Formula uses all strata levels in order (highest to lowest/most nested)
# e.g., ["Region", "Site"] becomes ~Region/Site
hier_formula_str <- paste0("~", paste(strata_list, collapse = "/"))
hier_formula <- as.formula(hier_formula_str)
cat("Hierarchical formula:", hier_formula_str, "\n")

# Calculate genetic distance with mean imputation
cat("Calculating genetic distance matrix (with mean imputation of missing data)...\n")

# Extract genotype matrix with mean imputation
gen_matrix <- tab(genind_obj, NA.method = "mean")

# Check for missing data and calculate statistics
has_missing <- any(is.na(genind_obj@tab))
if (has_missing) {
  tab_raw <- genind_obj@tab
  ploidy <- genind_obj@ploidy[1]
  n_alleles_per_locus <- ncol(tab_raw) / n_loci_before
  
  # Count loci with missing data
  loci_with_missing <- logical(n_loci_before)
  for (i in 1:n_loci_before) {
    start_col <- (i - 1) * n_alleles_per_locus + 1
    end_col <- i * n_alleles_per_locus
    locus_cols <- start_col:end_col
    loci_with_missing[i] <- any(is.na(tab_raw[, locus_cols, drop = FALSE]))
  }
  n_loci_with_missing <- sum(loci_with_missing)
  cat("Found", n_loci_with_missing, "loci with missing data out of", n_loci_before, "total loci\n")
  
  # Calculate distribution of loci among individuals (before imputation)
  loci_per_ind <- numeric(nInd(genind_obj))
  for (j in 1:nInd(genind_obj)) {
    complete_loci <- 0
    for (i in 1:n_loci_before) {
      start_col <- (i - 1) * n_alleles_per_locus + 1
      end_col <- i * n_alleles_per_locus
      if (!any(is.na(tab_raw[j, start_col:end_col]))) {
        complete_loci <- complete_loci + 1
      }
    }
    loci_per_ind[j] <- complete_loci
  }
} else {
  n_loci_with_missing <- 0
  loci_per_ind <- rep(n_loci_before, nInd(genind_obj))
  cat("All loci complete (no missing data detected)\n")
}

n_loci_after <- n_loci_before
n_loci_removed <- 0

cat("Loci per individual - Mean:", round(mean(loci_per_ind), 2), 
    ", Median:", median(loci_per_ind),
    ", Min:", min(loci_per_ind),
    ", Max:", max(loci_per_ind), "\n")

genetic_dist <- dist.binary(gen_matrix, method = 1)

# Perform AMOVA
cat("Performing AMOVA with", nperm, "permutations...\n")

# Verify strata are properly set
if (is.null(strata(genind_obj))) {
  stop("Strata not properly set in genind object")
}

# Check that all individuals have strata
if (nrow(strata(genind_obj)) != nInd(genind_obj)) {
  stop("Mismatch between number of individuals (", nInd(genind_obj), 
       ") and strata rows (", nrow(strata(genind_obj)), ")")
}

# Verify we have individuals in all required strata levels
for (stratum_col in strata_list) {
  stratum_values <- unique(strata(genind_obj)[[stratum_col]])
  if (length(stratum_values) < 1) {
    stop("No individuals found in stratum level: ", stratum_col)
  }
  cat("Stratum level '", stratum_col, "' has ", length(stratum_values), " unique values\n", sep = "")
}

amova_result <- poppr.amova(
  genind_obj,
  hier = hier_formula,
  dist = genetic_dist,
  nperm = nperm
)

# Write summary to text file
cat("Writing AMOVA summary...\n")
sink(output_amova)
cat("AMOVA Results (with mean imputation of missing data)\n")
cat("====================================================\n\n")

# Data summary
cat("Data Summary\n")
cat("------------\n")
cat("Number of individuals:", length(ind_names), "\n")
cat("Number of populations:", length(unique(pop(genind_obj))), "\n")
cat("Initial number of loci:", n_loci_before, "\n")
cat("Loci after processing:", n_loci_after, "\n")
if (n_loci_removed > 0) {
  cat("Loci removed (contained missing data):", n_loci_removed, "\n")
  cat("Proportion of loci retained:", round(n_loci_after / n_loci_before * 100, 2), "%\n")
}
cat("\n")

# Loci distribution among individuals
cat("Loci Distribution Among Individuals\n")
cat("-----------------------------------\n")
cat("Mean loci per individual:", round(mean(loci_per_ind), 2), "\n")
cat("Median loci per individual:", median(loci_per_ind), "\n")
cat("Minimum loci per individual:", min(loci_per_ind), "\n")
cat("Maximum loci per individual:", max(loci_per_ind), "\n")
cat("Standard deviation:", round(sd(loci_per_ind), 2), "\n")
cat("\n")

# AMOVA parameters
cat("AMOVA Parameters\n")
cat("----------------\n")
cat("Hierarchical formula:", hier_formula_str, "\n")
cat("Number of permutations:", nperm, "\n")
cat("Missing data handling: Mean imputation\n")
cat("\n")

# AMOVA results
cat("AMOVA Results\n")
cat("-------------\n")
print(amova_result)
sink()

# Create variance components plot
cat("Creating variance components plot...\n")

# Extract variance components
var_data <- data.frame(
  Source = rownames(amova_result$componentsofcovariance),
  Percentage = amova_result$componentsofcovariance[, "%"],
  stringsAsFactors = FALSE
)

# Remove total row
var_data <- var_data[var_data$Source != "Total variations", ]

# Clean up source labels
var_data$Source <- gsub("Variations  ", "", var_data$Source)
var_data$Source <- gsub("Between samples Within ", "Within ", var_data$Source)
var_data$Source <- gsub("Between ", "Between ", var_data$Source)

# Set factor levels based on hierarchy (most nested to least nested)
if (nrow(var_data) > 0) {
  source_levels <- rev(unique(var_data$Source))
  var_data$Source <- factor(var_data$Source, levels = source_levels)
}

# Create plot
p <- ggplot(var_data, aes(x = Source, y = Percentage)) +
  geom_bar(stat = "identity", fill = "steelblue", color = "black") +
  labs(
    x = "Source of variation",
    y = "% of variation",
    title = ""
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )

# Save plot as PDF
ggsave(output_plot, plot = p, width = 5, height = 3, dpi = 300)

# Save plot as ggplot object (RDS)
saveRDS(p, file = output_plot_rds)

cat("AMOVA analysis complete!\n")
cat("Results saved to:", output_amova, "\n")
cat("Plot (PDF) saved to:", output_plot, "\n")
cat("Plot (RDS) saved to:", output_plot_rds, "\n")
