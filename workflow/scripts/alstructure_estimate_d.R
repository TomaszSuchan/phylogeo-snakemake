#!/usr/bin/env Rscript
# ALStructure latent-dimension (K) estimation.
#
# estimate_d() (Cabreros & Storey 2019, after Leek 2011) estimates the rank d of
# the latent admixture subspace = the number of ancestral populations K. We also
# write the latent-subspace eigenvalue spectrum (top max(k_values) eigenvalues of
# the G matrix from lse()) so the shared choose-K plot can draw a scree; the
# elbow corroborates the estimate_d value.

log_file <- snakemake@log[[1]]
dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")
on.exit({
  while (sink.number(type = "message") > 0) sink(type = "message")
  while (sink.number(type = "output") > 0) sink(type = "output")
  close(log_con)
}, add = TRUE)

if (!requireNamespace("alstructure", quietly = TRUE)) {
  stop("alstructure is not installed in this rule environment.")
}
suppressPackageStartupMessages({
  library(alstructure)
})

geno_rds <- snakemake@input[["geno_rds"]]
choose_k_file <- snakemake@output[["choose_k_results"]]
cv_summary_file <- snakemake@output[["cv_summary"]]
k_values <- as.integer(unlist(snakemake@params[["k_values"]]))

geno_data <- readRDS(geno_rds)
# ALStructure expects X as m loci x n individuals; geno is individuals x loci.
X <- t(geno_data$geno)
m <- nrow(X)
n <- ncol(X)
cat("=== ALStructure estimate_d ===\n")
cat("SNP matrix X:", m, "loci x", n, "individuals\n")

# estimate_d()/lse() build the G matrix via eigen()/D_binomial() and do not
# tolerate NA, so impute missing dosages with the per-locus mean (the same
# strategy alstructure() applies internally via impute_mean).
na_cells <- sum(is.na(X))
if (na_cells > 0) {
  cat("Imputing", na_cells, "missing dosages with per-locus means...\n")
  row_means <- rowMeans(X, na.rm = TRUE)
  row_means[is.na(row_means)] <- 0
  na_idx <- which(is.na(X), arr.ind = TRUE)
  X[na_idx] <- row_means[na_idx[, 1]]
}

cat("Estimating latent dimension d (Leek 2011 plateau estimator)...\n")
d_hat <- as.integer(estimate_d(X))
cat("Estimated d (suggested K):", d_hat, "\n")

# Latent-subspace eigenvalue scree across the requested K range.
max_k <- max(k_values)
eig_d <- max(1L, min(as.integer(max_k), n - 1L))
eig_values <- lse(X, d = eig_d)$values
ks <- seq_len(eig_d)
cv_summary <- data.frame(
  K = ks,
  median = as.numeric(eig_values),
  min = as.numeric(eig_values),
  max = as.numeric(eig_values)
)
write.table(
  cv_summary,
  file = cv_summary_file,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

writeLines(
  c(
    "ALStructure estimates the number of ancestral populations as the dimension",
    "of the latent admixture subspace (estimate_d; Cabreros & Storey 2019, after",
    "Leek 2011). The eigenvalue scree (cv_summary) corroborates it: look for the",
    "elbow where latent-subspace eigenvalues drop toward the noise floor.",
    paste0("Suggested K (estimate_d): ", d_hat),
    "",
    "K\tlatent_subspace_eigenvalue",
    paste(cv_summary$K, signif(cv_summary$median, 6), sep = "\t")
  ),
  con = choose_k_file
)

cat("Eigenvalue scree written to:", cv_summary_file, "\n")
cat("choose-K results written to:", choose_k_file, "\n")
