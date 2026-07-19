#!/usr/bin/env Rscript
# Run the R related package (Coancestry) on a population-stratified VCF subset.
# Estimates Queller & Goodnight Rxy and optional additional moment estimators.

suppressPackageStartupMessages({
  if (!requireNamespace("related", quietly = TRUE)) {
    stop(
      "The related package is not installed. ",
      "The install_related Snakemake rule should install it before analysis."
    )
  }
  library(related)
  library(vcfR)
})

log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")
on.exit({
  sink(type = "message")
  sink(type = "output")
  close(log_file)
}, add = TRUE)
pdf(NULL)

vcf_file <- snakemake@input[["vcf"]]
output_file <- snakemake@output[["relatedness"]]
population_label <- snakemake@params[["population_label"]]
estimators <- snakemake@params[["estimators"]]

parse_gt_alleles <- function(gt_str) {
  if (is.na(gt_str) || gt_str %in% c("./.", ".|.", ".", "./.|.")) {
    return(c(0L, 0L))
  }
  parts <- strsplit(gt_str, "[/|]", fixed = FALSE)[[1]]
  alleles <- suppressWarnings(as.integer(parts))
  if (length(alleles) != 2L || any(is.na(alleles))) {
    return(c(0L, 0L))
  }
  alleles <- alleles + 1L
  sort(alleles)
}

cat("Reading VCF:", vcf_file, "\n")
vcf <- read.vcfR(vcf_file, verbose = FALSE)
gt <- extract.gt(vcf, element = "GT", as.numeric = FALSE)
inds <- colnames(gt)
nloci <- nrow(gt)
cat("Individuals:", length(inds), " Loci:", nloci, "\n")

if (length(inds) < 2L) {
  stop("related/coancestry requires at least 2 individuals in the population subset")
}
if (nloci < 1L) {
  stop("No variant loci found in the population subset VCF")
}

gmat <- matrix(0L, nrow = length(inds), ncol = 2L * nloci)
for (j in seq_along(inds)) {
  for (i in seq_len(nloci)) {
    alleles <- parse_gt_alleles(gt[i, j])
    gmat[j, 2L * i - 1L] <- alleles[1L]
    gmat[j, 2L * i] <- alleles[2L]
  }
}

count_observed_alleles <- function(locus_idx) {
  vals <- gmat[, c(2L * locus_idx - 1L, 2L * locus_idx)]
  vals <- vals[vals > 0L]
  length(unique(vals))
}

observed_alleles <- vapply(seq_len(nloci), count_observed_alleles, integer(1))
keep_loci <- observed_alleles > 0L
dropped_loci <- sum(!keep_loci)
if (dropped_loci > 0L) {
  cat(
    "Dropping", dropped_loci, "all-missing loci before coancestry",
    "(related requires 1-127 observed alleles per locus)\n"
  )
}
if (!any(keep_loci)) {
  stop(
    "No informative loci remain for coancestry in population subset ",
    snakemake@wildcards[["stratum"]]
  )
}
kept_cols <- as.vector(rbind(2L * which(keep_loci) - 1L, 2L * which(keep_loci)))
gmat <- gmat[, kept_cols, drop = FALSE]
nloci <- sum(keep_loci)
cat("Retained loci:", nloci, "\n")

gdata <- data.frame(ID = inds, gmat, check.names = FALSE, stringsAsFactors = FALSE)
colnames(gdata)[1] <- "ID"

estimator_args <- as.list(estimators)
if (length(estimator_args) == 0L) {
  estimator_args <- list(quellergt = 1L)
}
estimator_args <- lapply(estimator_args, as.integer)

cat("Running coancestry with estimators:", paste(names(estimator_args), collapse = ", "), "\n")
output <- do.call(
  coancestry,
  c(list(genotype.data = gdata), estimator_args)
)

rel <- as.data.frame(output$relatedness, stringsAsFactors = FALSE)
rel$population <- population_label
rel$stratum <- snakemake@wildcards[["stratum"]]

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
write.table(rel, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Wrote", nrow(rel), "pairwise estimates to", output_file, "\n")
