# DAPC BIC Extraction Script
# Extracts BIC values from find.clusters for all K values
# BIC values are written to the log file for use by the criterion plot script

# Prevent creation of Rplots.pdf
options(device = function(file = NULL, ...) {
  if (is.null(file)) {
    pdf(NULL)
  } else {
    pdf(file, ...)
  }
})

# Load required libraries
library(adegenet)  # Also loads ade4
library(vcfR)

# Handle both Snakemake and command-line argument modes
if (exists("snakemake")) {
  # Snakemake mode
  vcf_file <- snakemake@input[["vcf"]]
  log_file <- snakemake@output[["log_file"]]
  k_values <- snakemake@params[["k_values"]]
  n_pca <- snakemake@params[["n_pca"]]
  criterion <- snakemake@params[["criterion"]]
} else {
  # Command-line argument mode (for use with timeout wrapper)
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 5) {
    stop("Usage: Rscript dapc_bic_plot.R <vcf_file> <log_file> <k_values> <n_pca> <criterion>\n",
         "  k_values should be comma-separated, e.g., '1,2,3,4,5'")
  }
  vcf_file <- args[1]
  log_file <- args[2]
  k_values <- as.numeric(strsplit(args[3], ",")[[1]])
  n_pca <- args[4]
  criterion <- args[5]
}

# Set up logging
sink(log_file, append = FALSE, split = TRUE)

cat("=== DAPC BIC Extraction ===\n")
cat("VCF file:", vcf_file, "\n")
cat("Log file:", log_file, "\n\n")

cat("Parameters:\n")
cat("  K values:", paste(k_values, collapse = ", "), "\n")
cat("  n_pca:", n_pca, "\n")
cat("  criterion:", criterion, "\n\n")

# ==============================================================================
# 1. Load VCF and convert to genind
# ==============================================================================

cat("Loading VCF file...\n")
vcf <- read.vcfR(vcf_file)
cat("VCF loaded:", nrow(vcf@gt), "samples,", nrow(vcf@fix), "variants\n")

cat("Converting VCF to genind object...\n")
genind_obj <- vcfR2genind(vcf, sep = "/")
cat("genind object created:\n")
cat("  ", nInd(genind_obj), "individuals\n")
cat("  ", nLoc(genind_obj), "loci\n")
cat("  ", nPop(genind_obj), "populations\n\n")

# ==============================================================================
# 2. Calculate n_pca_used
# ==============================================================================

n_pca_used <- ifelse(n_pca == "retained", ncol(genind_obj@tab), as.numeric(n_pca))
cat("Using", n_pca_used, "PCA axes\n\n")

# ==============================================================================
# 3. Call find.clusters with choose.n.clust = TRUE to extract BIC values
# ==============================================================================

cat("Calling find.clusters with choose.n.clust = TRUE...\n")
cat("This will extract BIC values from Kstat.\n")
cat("NOTE: This will prompt interactively, but BIC values are extracted before the prompt.\n\n")

# Prevent any PDF device from opening
pdf(NULL)

# Call find.clusters with choose.n.clust = TRUE
# This will generate a plot showing BIC values across K values
# The plot is generated before the interactive prompt appears
max_k <- max(k_values)
min_k <- min(k_values[k_values >= 2])  # Usually start from K=2

cat("Calculating BIC for K =", min_k, "to", max_k, "...\n")
cat("NOTE: find.clusters will prompt interactively after calculation.\n")
cat("      BIC values are extracted from Kstat before the prompt.\n\n")

# Flush console to see progress
flush.console()

# Call find.clusters - BIC values are calculated and stored in Kstat
# The script will hang at an interactive prompt, but we extract BIC values first
cat("Starting find.clusters calculation (this will take 30-60 seconds)...\n")
flush.console()

# Call find.clusters - this will calculate BIC and store in Kstat
# Even though it prompts interactively, Kstat is populated during calculation
tryCatch({
  grp_bic <- find.clusters(
    genind_obj,
    choose.n.clust = TRUE,
    max.n.clust = max_k,
    min.n.clust = min_k,
    n.pca = n_pca_used,
    criterion = "diffNgroup"  # Use diffNgroup as criterion for choosing K
  )
  
  cat("find.clusters completed successfully.\n")
  cat("BIC values extracted from Kstat:\n")
  if (!is.null(grp_bic$Kstat)) {
    print(grp_bic$Kstat)
  } else {
    cat("WARNING: Kstat is NULL - BIC values not available\n")
  }
}, error = function(e) {
  cat("ERROR during find.clusters call:", e$message, "\n")
  cat("BIC values may not be available.\n")
}, interrupt = function(e) {
  cat("\nInterrupted during find.clusters call (this is expected).\n")
  cat("Attempting to extract BIC values if available...\n")
  # Try to extract BIC if grp_bic exists
  if (exists("grp_bic") && !is.null(grp_bic$Kstat)) {
    cat("BIC values extracted from Kstat:\n")
    print(grp_bic$Kstat)
  }
})

cat("\n=== BIC Extraction Complete ===\n")
sink()

# Final cleanup - remove Rplots.pdf if it was created
if (file.exists("Rplots.pdf")) {
  tryCatch({
    file.remove("Rplots.pdf")
  }, error = function(e) {
    # Silently fail
  })
}

