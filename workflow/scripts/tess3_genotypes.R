#!/usr/bin/env Rscript
# Convert the analysis VCF to a tess3r-compatible genotype dosage matrix.

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
geno_rds <- snakemake@output[["geno_rds"]]
ploidy <- as.integer(snakemake@params[["ploidy"]])

dir.create(dirname(geno_rds), recursive = TRUE, showWarnings = FALSE)

cat("=== tess3r genotype conversion ===\n")
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
    "Increase tess3.ploidy if this is expected."
  )
}

cat("Genotype matrix:", nrow(geno), "individuals x", ncol(geno), "loci\n")
cat(
  "Missing genotypes:", sum(is.na(geno)),
  " (", round(100 * mean(is.na(geno)), 2), "%)\n",
  sep = ""
)
cat(
  "Dosage range:", min(geno, na.rm = TRUE), "to",
  max(geno, na.rm = TRUE), "\n"
)

saveRDS(
  list(
    geno = geno,
    samples = sample_names,
    ploidy = ploidy
  ),
  file = geno_rds
)
cat("Genotype RDS written to:", geno_rds, "\n")
