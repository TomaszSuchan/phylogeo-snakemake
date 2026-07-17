#!/usr/bin/env Rscript
# Run PC-Relate (GENESIS) on the full filtered VCF.
# PC-AiR identifies ancestry-representative PCs; PC-Relate conditions on them
# to estimate kinship and IBD probabilities robust to population structure.

suppressPackageStartupMessages({
  library(GENESIS)
  library(SeqArray)
  library(SNPRelate)
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
output_file <- snakemake@output[["kinship"]]
n_pcs <- as.integer(snakemake@params[["n_pcs"]])
ld_r2 <- as.numeric(snakemake@params[["ld_r2"]])
ld_window <- as.integer(snakemake@params[["ld_window"]])
maf <- as.numeric(snakemake@params[["maf"]])
return_ibd_probs <- isTRUE(snakemake@params[["return_ibd_probs"]])

gds_file <- tempfile(fileext = ".gds")
on.exit(unlink(gds_file), add = TRUE)

cat("Converting VCF to GDS:", vcf_file, "\n")
seqVCF2GDS(vcf_file, gds_file, verbose = FALSE)
gds <- seqOpen(gds_file)
on.exit(seqClose(gds), add = TRUE)

nsamp <- length(seqGetData(gds, "sample.id"))
cat("Individuals in GDS:", nsamp, "\n")
if (nsamp < 2L) {
  stop("PC-Relate requires at least 2 individuals (found ", nsamp, ")")
}

set.seed(100)
cat("LD pruning (r2 <=", ld_r2, ", window =", ld_window, ")\n")
snpset <- snpgdsLDpruning(
  gds,
  ld.threshold = ld_r2,
  window.size = ld_window,
  maf = maf,
  verbose = FALSE
)
snps.use <- unlist(snpset, use.names = FALSE)
cat("SNPs after LD pruning:", length(snps.use), "\n")
if (length(snps.use) < 50L) {
  stop("Too few SNPs (", length(snps.use), ") after LD pruning for PC-Relate")
}

cat("Computing KING kinship matrix for PC-AiR\n")
king <- snpgdsIBDKING(gds, snp.id = snps.use, verbose = FALSE)
kingMat <- snpKIN2kinship(king$kinship)

seqData <- SeqVarData(gds)
cat("Running PC-AiR\n")
pcair_res <- pcair(
  seqData,
  kinobj = kingMat,
  divobj = kingMat,
  snp.include = snps.use,
  verbose = FALSE
)

npcs <- min(n_pcs, ncol(pcair_res$vectors))
if (npcs < 1L) {
  stop("PC-AiR returned no usable principal components")
}
cat("Using", npcs, "ancestry PCs for PC-Relate\n")

training.set <- pcair_res$unrels
if (length(training.set) < 2L) {
  training.set <- rownames(pcair_res$vectors)
  cat(
    "WARNING: PC-AiR unrelated set has < 2 samples; using all",
    length(training.set), "individuals as training set\n"
  )
}

seqData <- SeqVarBlockIterator(seqData, variant.include = snps.use, verbose = FALSE)
cat("Running PC-Relate\n")
pcrel <- pcrelate(
  seqData,
  pcs = pcair_res$vectors[, seq_len(npcs), drop = FALSE],
  training.set = training.set,
  scaleKinship = TRUE,
  verbose = FALSE
)

kin <- as.data.frame(pcrel$kinBwtn, stringsAsFactors = FALSE)
if (!return_ibd_probs) {
  keep_cols <- intersect(
    colnames(kin),
    c("ID1", "ID2", "kin", "k0", "k1", "k2")
  )
  kin <- kin[, keep_cols, drop = FALSE]
}

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
write.table(kin, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Wrote", nrow(kin), "pairwise kinship estimates to", output_file, "\n")
