#!/usr/bin/env Rscript
# Convert the analysis VCF to an LEA/sNMF .geno genotype file.
#
# The .geno format (LEA package) is one row per SNP and one column per
# individual, with each cell a single character: 0, 1, or 2 alternate-allele
# copies, or 9 for a missing genotype, and no separators between columns
# (see ?LEA::read.geno). Individuals are written in VCF sample order, which is
# the same order used by indpopdata.txt, so the Q matrices written later line up
# with indpopdata rows in the shared barplot/map scripts.

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

suppressPackageStartupMessages({
  library(vcfR)
})

vcf_file <- snakemake@input[["vcf"]]
indpopdata_file <- snakemake@input[["indpopdata"]]
geno_file <- snakemake@output[["geno"]]
samples_file <- snakemake@output[["samples"]]
ploidy <- as.integer(snakemake@params[["ploidy"]])

dir.create(dirname(geno_file), recursive = TRUE, showWarnings = FALSE)

cat("=== sNMF/LEA genotype conversion ===\n")
cat("VCF:", vcf_file, "\n")
cat("ploidy:", ploidy, "\n")

vcf <- read.vcfR(vcf_file, verbose = FALSE)
sample_names <- colnames(vcf@gt)[-1]
cat("VCF samples:", length(sample_names), "\n")
cat("VCF variants:", nrow(vcf@fix), "\n")

indpopdata <- read.table(
  indpopdata_file,
  header = TRUE,
  sep = "\t",
  comment.char = "",
  fill = TRUE,
  blank.lines.skip = TRUE
)
if (!all(sample_names %in% indpopdata$Ind)) {
  missing_inds <- sample_names[!sample_names %in% indpopdata$Ind]
  stop(
    "Samples in VCF are missing from indpopdata: ",
    paste(head(missing_inds, 20), collapse = ", "),
    if (length(missing_inds) > 20) " ..."
  )
}

gt <- extract.gt(vcf, element = "GT", as.numeric = FALSE)
cat("Converting genotypes to alternate-allele dosages...\n")
# Individuals x loci dosage matrix (same convention as tess3_genotypes.R).
geno_dosage <- apply(gt, c(1, 2), function(x) {
  if (is.na(x) || x %in% c(".", "./.", ".|.")) {
    return(NA_real_)
  }
  alleles <- unlist(strsplit(x, "[/|]"))
  if (length(alleles) == 0 || any(alleles == ".")) {
    return(NA_real_)
  }
  sum(as.integer(alleles) > 0)
})
geno <- t(geno_dosage)
storage.mode(geno) <- "double"
rownames(geno) <- sample_names
colnames(geno) <- rownames(gt)

if (any(geno > ploidy, na.rm = TRUE)) {
  bad <- geno[geno > ploidy & !is.na(geno)]
  stop(
    "Genotype dosages exceed configured ploidy (", ploidy, "). ",
    "Maximum observed dosage: ", max(bad), ". ",
    "Increase snmf.ploidy if this is expected."
  )
}

cat("Genotype matrix:", nrow(geno), "individuals x", ncol(geno), "loci\n")
cat(
  "Missing genotypes:", sum(is.na(geno)),
  " (", round(100 * mean(is.na(geno)), 2), "%)\n",
  sep = ""
)

# LEA .geno: rows = SNPs, cols = individuals, missing coded as 9, no separators.
geno_int <- t(geno)                      # loci x individuals
geno_int[is.na(geno_int)] <- 9
storage.mode(geno_int) <- "integer"

# Warn about loci with no scored genotypes; sNMF tolerates them but they add no
# information and can slow convergence.
all_missing <- rowSums(geno_int != 9) == 0
if (any(all_missing)) {
  cat(
    "WARNING:", sum(all_missing),
    "loci have no scored genotypes (all individuals missing).\n"
  )
}

geno_lines <- apply(geno_int, 1, function(row) paste0(row, collapse = ""))
writeLines(geno_lines, con = geno_file)
cat("Wrote .geno file:", geno_file, "\n")
cat("  dimensions:", nrow(geno_int), "loci x", ncol(geno_int), "individuals\n")

writeLines(sample_names, con = samples_file)
cat("Wrote sample order file:", samples_file, "\n")
