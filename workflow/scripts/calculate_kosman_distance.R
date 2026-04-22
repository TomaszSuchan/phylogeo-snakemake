library(PopGenReport)
library(vcfR)

# Read VCF and convert to genind
vcf <- read.vcfR(snakemake@input[['vcf']])
genind_obj <- vcfR2genind(vcf)

cat("Calculating Kosman distance for", nInd(genind_obj), "individuals\n")

# Calculate Kosman distance
kosman <- gd.kosman(genind_obj)

# Get distance matrix
dist_matrix <- as.matrix(kosman$geneticdist)

# Use individual names from genind object (these will have proper names from VCF)
rownames(dist_matrix) <- indNames(genind_obj)
colnames(dist_matrix) <- indNames(genind_obj)

cat("Sample names (first 3):", head(indNames(genind_obj), 3), "\n")

# Write with row and column names
write.table(dist_matrix,
            file = snakemake@output[['dist']],
            sep = '\t',
            quote = FALSE,
            col.names = NA)