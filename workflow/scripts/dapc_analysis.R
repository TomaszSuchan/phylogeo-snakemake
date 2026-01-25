# DAPC Analysis Script
# Performs Discriminant Analysis of Principal Components for population clustering
# Runs for a specific K value (provided as wildcard)

# Prevent creation of Rplots.pdf - set before loading libraries
# Set environment variable to redirect Rplots.pdf
Sys.setenv(R_INSTALL_PKG = "")
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

# Explicitly prevent Rplots.pdf creation
#pdf(NULL)
## Also try to delete it if it exists
#if (file.exists("Rplots.pdf")) {
#  tryCatch({
#    file.remove("Rplots.pdf")
#    cat("Removed existing Rplots.pdf\n")
#  }, error = function(e) {
#    cat("Could not remove Rplots.pdf:", e$message, "\n")
#  })
#}

# Set up logging
log_file <- snakemake@output[["log_file"]]
sink(log_file, append = FALSE, split = TRUE)

cat("=== DAPC Analysis (K =", snakemake@params[["k"]], ") ===\n")
cat("VCF file:", snakemake@input[["vcf"]], "\n")
cat("Output directory:", dirname(snakemake@output[["results_rds"]]), "\n\n")

# Read parameters
k <- as.numeric(snakemake@params[["k"]])
n_pca <- snakemake@params[["n_pca"]]
n_da <- snakemake@params[["n_da"]]
criterion <- snakemake@params[["criterion"]]

cat("Parameters:\n")
cat("  K:", k, "\n")
cat("  n_pca:", n_pca, "\n")
cat("  n_da:", n_da, "\n")
cat("  criterion:", criterion, "\n\n")

# ==============================================================================
# 1. Load VCF and convert to genind
# ==============================================================================

cat("Loading VCF file...\n")
vcf <- read.vcfR(snakemake@input[["vcf"]], verbose = FALSE)
cat("VCF loaded:", nrow(vcf@gt), "samples,", nrow(vcf@fix), "variants\n")

cat("Converting VCF to genind object...\n")
genind_obj <- vcfR2genind(vcf, sep = "/")
cat("genind object created:\n")
cat("  ", nInd(genind_obj), "individuals\n")
cat("  ", nLoc(genind_obj), "loci\n")
cat("  ", nPop(genind_obj), "populations\n\n")

# ==============================================================================
# 2. Find clusters for specified K
# ==============================================================================

cat("Finding clusters for K =", k, "...\n")
n_pca_used <- ifelse(n_pca == "retained", ncol(genind_obj@tab), as.numeric(n_pca))

# Initialize
criterion_value <- NA

# Get clustering for the target K
grp <- find.clusters(genind_obj, n.clust = k, choose.n.clust = FALSE,
                     n.pca = n_pca_used,
                     criterion = criterion)

# Check if find.clusters stores PCA coordinates internally
# The result object might have a 'stat' component or PCA information
# If available, use those coordinates instead of recalculating
if ("stat" %in% names(grp) && !is.null(grp$stat)) {
  cat("NOTE: find.clusters result has 'stat' component:", grp$stat, "\n")
  # The 'stat' component might contain useful information, but it's not the BIC we need
}

# Get the PCA coordinates that find.clusters used internally
# Use the genind tab and perform PCA with the same settings
genind_tab <- genind_obj@tab

# Clean data - replace infinite and NaN values with column means
genind_tab[is.infinite(genind_tab)] <- NA
genind_tab[is.nan(genind_tab)] <- NA

# Replace NA values with column means (or 0 if all NA)
for (j in 1:ncol(genind_tab)) {
  na_indices <- is.na(genind_tab[, j])
  if (any(na_indices)) {
    col_mean <- mean(genind_tab[, j], na.rm = TRUE)
    if (is.na(col_mean) || !is.finite(col_mean)) col_mean <- 0
    genind_tab[na_indices, j] <- col_mean
  }
}

# Perform PCA with the same number of axes as used in find.clusters
# Use dudi.pca from ade4 (which adegenet uses internally)
tryCatch({
  n_pca_final <- min(n_pca_used, nInd(genind_obj) - 1, ncol(genind_tab))
  
  if (n_pca_final < 1) {
    cat("WARNING: Cannot perform PCA (n_pca_final =", n_pca_final, ")\n")
    criterion_value <- NA
  } else {
    pca_dudi <- dudi.pca(genind_tab, center = TRUE, scale = FALSE, 
                         scannf = FALSE, nf = n_pca_final)
    pca_coords <- pca_dudi$li  # Individual coordinates
    
    n_individuals <- nrow(pca_coords)
    n_pca_axes <- ncol(pca_coords)
    
    cat("PCA coordinates: ", n_individuals, "individuals x", n_pca_axes, "axes\n")
    
    # Calculate within-cluster sum of squares (WSS) for BIC
    wss <- 0
    n_clusters_with_data <- 0
    for (cluster_id in unique(grp$grp)) {
      cluster_indices <- which(grp$grp == cluster_id)
      if (length(cluster_indices) > 1 && n_pca_axes > 0) {
        # Ensure indices are valid
        valid_indices <- cluster_indices[cluster_indices <= n_individuals]
        if (length(valid_indices) > 1) {
          cluster_coords <- pca_coords[valid_indices, , drop = FALSE]
          if (nrow(cluster_coords) > 1 && ncol(cluster_coords) > 0) {
            cluster_center <- colMeans(cluster_coords)
            if (all(is.finite(cluster_center))) {
              cluster_distances <- rowSums((cluster_coords - 
                                            matrix(cluster_center, nrow = nrow(cluster_coords), 
                                                   ncol = ncol(cluster_coords), byrow = TRUE))^2)
              wss <- wss + sum(cluster_distances)
              n_clusters_with_data <- n_clusters_with_data + 1
            }
          }
        }
      }
    }
    
    cat("WSS:", wss, ", Clusters with data:", n_clusters_with_data, "\n")
    
    # Calculate BIC using adegenet's formula: BIC = n * log(WSS/n) + k * log(n)
    if (wss > 0 && is.finite(wss) && n_individuals > 0) {
      criterion_value <- n_individuals * log(wss / n_individuals) + k * log(n_individuals)
      cat("BIC value:", criterion_value, "\n")
    } else {
      criterion_value <- NA
      cat("WARNING: Could not calculate BIC (WSS =", wss, ", n =", n_individuals, ")\n")
    }
  }
}, error = function(e) {
  cat("ERROR calculating BIC:", e$message, "\n")
  cat("Error traceback:\n")
  print(traceback())
  criterion_value <<- NA
})

# Get clustering for the target K (for DAPC analysis)
# Only do this if we haven't already obtained it from BIC extraction
if (!exists("grp") || is.null(grp)) {
  cat("Getting clustering for K =", k, "for DAPC analysis...\n")
  grp <- find.clusters(genind_obj, n.clust = k, choose.n.clust = FALSE,
                       n.pca = n_pca_used,
                       criterion = criterion)
}

# Criterion value calculation removed - using native BIC from find.clusters instead

# ==============================================================================
# 3. Perform DAPC
# ==============================================================================

cat("Performing DAPC...\n")
dapc_result <- dapc(genind_obj, grp$grp, 
                    n.pca = ifelse(n_pca == "retained", ncol(genind_obj@tab), as.numeric(n_pca)),
                    n.da = ifelse(n_da == "all", k - 1, as.numeric(n_da)))

cat("DAPC completed:\n")
cat("  ", ncol(dapc_result$ind.coord), "discriminant axes\n")
cat("  ", length(unique(grp$grp)), "clusters\n\n")

# ==============================================================================
# 4. Extract results
# ==============================================================================

# Cluster assignments
cluster_assignments <- data.frame(
  Ind = indNames(genind_obj),
  Cluster = as.character(grp$grp)
)

# Membership probabilities
membership_probs <- as.data.frame(dapc_result$posterior)
membership_probs$Ind <- indNames(genind_obj)
membership_probs <- membership_probs[, c("Ind", setdiff(names(membership_probs), "Ind"))]

# ==============================================================================
# 5. Calculate variance explained
# ==============================================================================

# Calculate variance explained by DA axes
da_eig <- dapc_result$eig
da_var_exp <- da_eig / sum(da_eig) * 100

cat("Variance explained by DA axes:\n")
for (i in 1:min(5, length(da_var_exp))) {
  cat("  DA", i, ":", round(da_var_exp[i], 2), "%\n")
}
cat("\n")

# ==============================================================================
# 6. Create scatter plot using native adegenet scatter() function
# ==============================================================================

cat("Creating scatter plot using native adegenet scatter()...\n")

# Save plot using PDF device
dir.create(dirname(snakemake@output[["scatter_plot"]]), recursive = TRUE, showWarnings = FALSE)
pdf(snakemake@output[["scatter_plot"]], width = 8, height = 6)
scatter(dapc_result)
dev.off()
cat("Scatter plot saved to:", snakemake@output[["scatter_plot"]], "\n")

# For RDS, we'll save the dapc_result object itself since we can't save base R plots easily
# Users can recreate the plot with scatter(dapc_result) if needed
saveRDS(dapc_result, snakemake@output[["scatter_plot_rds"]])
cat("DAPC result saved to RDS (use scatter() to recreate plot):", snakemake@output[["scatter_plot_rds"]], "\n")

# ==============================================================================
# 7. Save results
# ==============================================================================

cat("\nSaving results...\n")

# Create results list
results <- list(
  dapc_result = dapc_result,
  cluster_assignments = cluster_assignments,
  membership_probs = membership_probs,
  k = k,
  criterion_value = criterion_value,
  criterion_type = "BIC",
  da_var_exp = da_var_exp
)

# Save RDS
dir.create(dirname(snakemake@output[["results_rds"]]), recursive = TRUE, showWarnings = FALSE)
saveRDS(results, snakemake@output[["results_rds"]])
cat("Results RDS saved to:", snakemake@output[["results_rds"]], "\n")

# Save cluster assignments
write.table(cluster_assignments, snakemake@output[["cluster_assignments"]], 
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("Cluster assignments saved to:", snakemake@output[["cluster_assignments"]], "\n")

# Save membership probabilities
write.table(membership_probs, snakemake@output[["membership_probs"]], 
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("Membership probabilities saved to:", snakemake@output[["membership_probs"]], "\n")

cat("\n=== DAPC Analysis Complete (K =", k, ") ===\n")
sink()

# Final cleanup - remove Rplots.pdf if it was created
if (file.exists("Rplots.pdf")) {
  tryCatch({
    file.remove("Rplots.pdf")
  }, error = function(e) {
    # Silently fail
  })
}
