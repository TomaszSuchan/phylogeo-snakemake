#!/usr/bin/env Rscript
# Calculate Euclidean genetic distance from VCF file
# Euclidean distance is calculated from genotype matrix (allele counts: 0, 1, 2)

library(vcfR)
library(adegenet)

# Read VCF file
cat("Reading VCF file...\n")
vcf <- read.vcfR(snakemake@input[["vcf"]])

# Convert to genind object
cat("Converting to genind object...\n")
genind_obj <- vcfR2genind(vcf)

cat("Calculating Euclidean distance for", nInd(genind_obj), "individuals\n")

# Extract genotype matrix (coded as allele counts: 0, 1, 2)
cat("Extracting genotype matrix...\n")
gen_matrix <- tab(genind_obj, NA.method = "mean")  # Replace NAs with mean

# Calculate Euclidean distance
cat("Calculating Euclidean distances...\n")
euclidean_dist <- dist(gen_matrix, method = "euclidean")

# Convert to matrix
dist_matrix <- as.matrix(euclidean_dist)

# Use individual names from genind object (these will match VCF and Kosman)
rownames(dist_matrix) <- indNames(genind_obj)
colnames(dist_matrix) <- indNames(genind_obj)

cat("Sample names (first 3):", head(indNames(genind_obj), 3), "\n")

# Write to file
write.table(
  dist_matrix,
  file = snakemake@output[["dist"]],
  sep = "\t",
  col.names = NA,
  row.names = TRUE,
  quote = FALSE
)

cat("Euclidean distance calculation complete!\n")
cat("Output saved to:", snakemake@output[["dist"]], "\n")
