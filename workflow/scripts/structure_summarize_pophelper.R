library(pophelper)

# Snakemake inputs/outputs
runs <- unlist(snakemake@input)
output_qmatrix <- snakemake@output[["qmatrix"]]
output_aligned <- snakemake@output[["aligned_runs"]]

# Read STRUCTURE runs
slist <- readQ(runs)

# Align clusters across replicates
slist_aligned <- alignK(slist)

# Calculate mean Q matrix across aligned replicates
# Get dimensions
n_inds <- nrow(slist_aligned[[1]])
n_clusters <- ncol(slist_aligned[[1]])
n_runs <- length(slist_aligned)

# Initialize matrix for summing
qmatrix_sum <- matrix(0, nrow = n_inds, ncol = n_clusters)

# Sum across all aligned runs
for (i in 1:n_runs) {
  qmatrix_sum <- qmatrix_sum + as.matrix(slist_aligned[[i]])
}

# Calculate mean
qmatrix_mean <- qmatrix_sum / n_runs

# Write Q matrix (tab-delimited, no row/column names, for mapmixture compatibility)
write.table(qmatrix_mean,
            file = output_qmatrix,
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

# Save aligned runs for potential future use
saveRDS(slist_aligned, file = output_aligned)

# Log summary
cat(sprintf("Processed %d runs with %d individuals and %d clusters (K)\n",
            n_runs, n_inds, n_clusters))
cat(sprintf("Mean Q matrix written to: %s\n", output_qmatrix))
cat(sprintf("Aligned runs saved to: %s\n", output_aligned))
