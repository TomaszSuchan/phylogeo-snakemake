#!/usr/bin/env Rscript
# Run sNMF (LEA) across all configured K values and summarise the cross-entropy
# criterion used to choose the number of ancestral populations.
#
# sNMF (Frichot et al. 2014) estimates ancestry coefficients by sparse
# non-negative matrix factorisation and is robust to departures from
# Hardy-Weinberg/linkage equilibrium, including the allele-frequency gradients
# produced by isolation by distance. The number of clusters is chosen from the
# cross-entropy of masked (held-out) genotypes: a lower cross-entropy means the
# fitted ancestry/allele-frequency matrices predict masked genotypes better, and
# the best K is typically where cross-entropy reaches a minimum or plateau.

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
results_rds <- snakemake@output[["results_rds"]]
choose_k_file <- snakemake@output[["choose_k_results"]]
cv_summary_file <- snakemake@output[["cv_summary"]]

k_values <- as.integer(unlist(snakemake@params[["k_values"]]))
repetitions <- as.integer(snakemake@params[["repetitions"]])
alpha <- as.numeric(snakemake@params[["alpha"]])
tolerance <- as.numeric(snakemake@params[["tolerance"]])
iterations <- as.integer(snakemake@params[["iterations"]])
ploidy <- as.integer(snakemake@params[["ploidy"]])
percentage <- as.numeric(snakemake@params[["percentage"]])
seed <- as.integer(snakemake@params[["seed"]])
threads <- as.integer(snakemake@threads)

cat("=== sNMF choose-K ===\n")
cat("K values:", paste(k_values, collapse = ", "), "\n")
cat("repetitions:", repetitions, "\n")
cat("alpha:", alpha, "\n")
cat("tolerance:", tolerance, "\n")
cat("iterations:", iterations, "\n")
cat("ploidy:", ploidy, "\n")
cat("entropy masking percentage:", percentage, "\n")
cat("seed:", seed, "\n")
cat("threads:", threads, "\n")

# Run sNMF in an isolated working directory: snmf() writes a .snmf project tree
# next to its input .geno, so giving each rule its own copy avoids collisions
# when multiple sNMF rules run in parallel.
work_dir <- file.path(dirname(results_rds), paste0(".snmf_chooseK_", Sys.getpid()))
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(work_dir, recursive = TRUE, force = TRUE), add = TRUE)
geno_copy <- file.path(work_dir, "snmf_input.geno")
file.copy(geno_file, geno_copy, overwrite = TRUE)

cat("Running sNMF across K...\n")
project <- snmf(
  geno_copy,
  K = k_values,
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

# Summarise masked-genotype cross-entropy (median/min/max across repetitions)
# per K, mirroring the tess3r choose-K table layout.
med <- min_ce <- max_ce <- numeric(length(k_values))
for (i in seq_along(k_values)) {
  ce <- as.numeric(cross.entropy(project, K = k_values[i]))
  med[i] <- median(ce)
  min_ce[i] <- min(ce)
  max_ce[i] <- max(ce)
}
cv_summary <- data.frame(K = k_values, median = med, min = min_ce, max = max_ce)

write.table(
  cv_summary,
  file = cv_summary_file,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

best_k <- cv_summary$K[which.min(cv_summary$median)]

writeLines(
  c(
    "Choose K with the lowest masked-genotype cross-entropy across replicate runs.",
    "Smaller values indicate better prediction of held-out genotypes;",
    "look for a clear minimum or a plateau where cross-entropy stops decreasing.",
    paste0("Suggested K (minimum median cross-entropy): ", best_k),
    "",
    "K\tmedian_cross_entropy\tmin\tmax",
    paste(
      cv_summary$K,
      signif(cv_summary$median, 6),
      signif(cv_summary$min, 6),
      signif(cv_summary$max, 6),
      sep = "\t"
    )
  ),
  con = choose_k_file
)

saveRDS(
  list(
    cv_summary = cv_summary,
    params = list(
      K = k_values,
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

cat("Cross-entropy summary written to:", cv_summary_file, "\n")
cat("choose-K results written to:", choose_k_file, "\n")
cat("Suggested K:", best_k, "\n")
