#!/usr/bin/env Rscript
# Run sNMF (LEA) for a single K and export the ancestry (Q) matrix.
#
# For the requested K, several independent runs are fit (repetitions) and the
# run with the lowest cross-entropy is retained, as recommended for sNMF because
# the least-squares optimisation can converge to different local optima. The Q
# matrix is written without row/column names, in VCF/indpopdata sample order, so
# the shared mapmixture barplot and map scripts can align it to indpopdata rows.

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

if (!requireNamespace("LEA", quietly = TRUE)) {
  stop("The LEA package (which provides sNMF) is not installed in this rule environment.")
}
suppressPackageStartupMessages({
  library(LEA)
})

geno_file <- snakemake@input[["geno"]]
samples_file <- snakemake@input[["samples"]]
qmatrix_file <- snakemake@output[["qmatrix"]]
cross_entropy_file <- snakemake@output[["cross_entropy"]]
results_rds <- snakemake@output[["results_rds"]]

k <- as.integer(snakemake@params[["k"]])
repetitions <- as.integer(snakemake@params[["repetitions"]])
alpha <- as.numeric(snakemake@params[["alpha"]])
tolerance <- as.numeric(snakemake@params[["tolerance"]])
iterations <- as.integer(snakemake@params[["iterations"]])
ploidy <- as.integer(snakemake@params[["ploidy"]])
percentage <- as.numeric(snakemake@params[["percentage"]])
seed <- as.integer(snakemake@params[["seed"]])
threads <- as.integer(snakemake@threads)

cat("=== sNMF analysis ===\n")
cat("K:", k, "\n")
cat("repetitions:", repetitions, "\n")
cat("alpha:", alpha, "\n")
cat("tolerance:", tolerance, "\n")
cat("iterations:", iterations, "\n")
cat("ploidy:", ploidy, "\n")
cat("seed:", seed, "\n")
cat("threads:", threads, "\n")

dir.create(dirname(qmatrix_file), recursive = TRUE, showWarnings = FALSE)

sample_names <- readLines(samples_file)
sample_names <- sample_names[nzchar(sample_names)]

# Isolated working directory (see snmf_choose_k.R for the rationale).
work_dir <- file.path(dirname(qmatrix_file), paste0(".snmf_K", k, "_", Sys.getpid()))
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(work_dir, recursive = TRUE, force = TRUE), add = TRUE)
geno_copy <- file.path(work_dir, "snmf_input.geno")
file.copy(geno_file, geno_copy, overwrite = TRUE)

cat("Running sNMF...\n")
project <- snmf(
  geno_copy,
  K = k,
  project = "new",
  repetitions = repetitions,
  entropy = TRUE,
  percentage = percentage,
  alpha = alpha,
  tolerance = tolerance,
  iterations = iterations,
  ploidy = ploidy,
  CPU = threads,
  seed = seed
)

ce <- as.numeric(cross.entropy(project, K = k))
best_run <- which.min(ce)
cat("Cross-entropy per run:", paste(signif(ce, 6), collapse = ", "), "\n")
cat("Best run (lowest cross-entropy):", best_run, "\n")

qmat <- as.matrix(Q(project, K = k, run = best_run))
if (nrow(qmat) != length(sample_names)) {
  stop("sNMF Q matrix has ", nrow(qmat), " rows, but expected ", length(sample_names))
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

write.table(
  data.frame(K = k, run = best_run, cross_entropy = ce[best_run]),
  file = cross_entropy_file,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
cat("Cross-entropy summary written to:", cross_entropy_file, "\n")

saveRDS(
  list(
    qmatrix = qmat,
    samples = sample_names,
    cross_entropy = ce,
    best_run = best_run,
    params = list(
      K = k,
      repetitions = repetitions,
      alpha = alpha,
      tolerance = tolerance,
      iterations = iterations,
      ploidy = ploidy,
      percentage = percentage,
      seed = seed,
      threads = threads
    )
  ),
  file = results_rds
)
cat("RDS written to:", results_rds, "\n")
cat("=== sNMF analysis complete ===\n")
