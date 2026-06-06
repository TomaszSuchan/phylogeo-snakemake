#!/usr/bin/env Rscript
# ALStructure ancestry estimation for one K.
# Fits the PSD admixture model by latent-subspace estimation + alternating least
# squares (Cabreros & Storey 2019); deterministic up to the random ALS init,
# which is fixed with a seed for reproducibility.

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
  stop(
    "alstructure is not installed in this rule environment. The ",
    "install_alstructure Snakemake rule should install it before analysis."
  )
}
suppressPackageStartupMessages({
  library(alstructure)
})

geno_rds <- snakemake@input[["geno_rds"]]
qmatrix_file <- snakemake@output[["qmatrix"]]
results_rds <- snakemake@output[["results_rds"]]

k <- as.integer(snakemake@params[["k"]])
svd_method <- as.character(snakemake@params[["svd_method"]])
tolerance <- as.numeric(snakemake@params[["tolerance"]])
max_iters <- as.integer(snakemake@params[["max_iters"]])
order_method <- as.character(snakemake@params[["order_method"]])
seed <- as.integer(snakemake@params[["seed"]])

cat("=== ALStructure Analysis ===\n")
cat("K:", k, "\n")
cat("svd_method:", svd_method, "\n")
cat("tolerance:", tolerance, "\n")
cat("max_iters:", max_iters, "\n")
cat("order_method:", order_method, "\n")
cat("seed:", seed, "\n\n")

if (k < 2) {
  stop(
    "ALStructure requires at least 2 latent populations (K >= 2): estimate_d ",
    "enforces d >= 2 and the latent-subspace factorisation collapses at d = 1. ",
    "K = 1 is skipped automatically when targets are driven by k_values."
  )
}

dir.create(dirname(qmatrix_file), recursive = TRUE, showWarnings = FALSE)

geno_data <- readRDS(geno_rds)
sample_names <- geno_data$samples
# ALStructure expects X as m loci x n individuals; geno is individuals x loci.
X <- t(geno_data$geno)
cat("SNP matrix X:", nrow(X), "loci x", ncol(X), "individuals\n")
cat(
  "Missing genotypes:", sum(is.na(X)),
  " (", round(100 * mean(is.na(X)), 2), "%) - imputed internally by alstructure()\n",
  sep = ""
)

set.seed(seed)
fit <- alstructure(
  X = X,
  d_hat = k,
  svd_method = svd_method,
  tol = tolerance,
  max_iters = max_iters,
  order_method = order_method
)

# Q_hat is d x n (each column an individual); transpose to individuals x K.
qmat <- t(fit$Q_hat)
if (nrow(qmat) != length(sample_names)) {
  stop(
    "ALStructure Q matrix has ", nrow(qmat),
    " rows, but expected ", length(sample_names), " individuals."
  )
}
colnames(qmat) <- paste0("Cluster", seq_len(ncol(qmat)))
rownames(qmat) <- sample_names

write.table(
  qmat,
  file = qmatrix_file,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)
cat("Q matrix written to:", qmatrix_file, "\n")

saveRDS(
  list(
    Q_hat = qmat,
    P_hat = fit$P_hat,
    rowspace = fit$rowspace,
    iter = fit$iter,
    samples = sample_names,
    params = list(
      K = k,
      svd_method = svd_method,
      tolerance = tolerance,
      max_iters = max_iters,
      order_method = order_method,
      seed = seed
    )
  ),
  file = results_rds
)
cat("RDS written to:", results_rds, "\n")
cat("=== ALStructure Analysis Complete ===\n")
