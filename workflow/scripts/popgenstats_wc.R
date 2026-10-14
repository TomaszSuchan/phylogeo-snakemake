#!/usr/bin/env Rscript
# Weir & Cockerham basic.stats, pairwise Fst (+ bootstrap CIs), rarefied allelic
# richness, and private alleles from the all-SNP VCF.

log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

suppressPackageStartupMessages({
  library(vcfR)
  library(adegenet)
  library(hierfstat)
  library(poppr)
})

vcf_file <- snakemake@input[["vcf"]]
popmap_file <- snakemake@input[["popmap"]]
out_overall <- snakemake@output[["overall"]]
out_popstats <- snakemake@output[["popstats"]]
out_pairwise <- snakemake@output[["pairwise_fst"]]
out_ci <- snakemake@output[["pairwise_fst_ci"]]

grouping_name <- as.character(snakemake@params[["grouping"]])
min_individuals <- suppressWarnings(as.integer(snakemake@params[["min_individuals"]]))
if (is.na(min_individuals)) min_individuals <- 2
nboot <- suppressWarnings(as.integer(snakemake@params[["nboot"]]))
if (is.na(nboot) || nboot < 1) nboot <- 1000

matrix_to_long <- function(mat, value_name) {
  if (is.null(dim(mat)) || nrow(mat) < 2) {
    return(data.frame(
      pop1 = character(0),
      pop2 = character(0),
      setNames(list(numeric(0)), value_name),
      stringsAsFactors = FALSE
    ))
  }
  rn <- rownames(mat)
  cn <- colnames(mat)
  if (is.null(rn)) rn <- as.character(seq_len(nrow(mat)))
  if (is.null(cn)) cn <- as.character(seq_len(ncol(mat)))
  ut <- which(upper.tri(mat), arr.ind = TRUE)
  data.frame(
    pop1 = rn[ut[, 1]],
    pop2 = cn[ut[, 2]],
    setNames(list(as.numeric(mat[ut])), value_name),
    stringsAsFactors = FALSE
  )
}

message("\n=== READING GENOTYPES (all-SNP) ===\n")
message(sprintf("VCF: %s\n", vcf_file))
vcf <- read.vcfR(vcf_file, verbose = FALSE)
genind_obj <- vcfR2genind(vcf)
message(sprintf("Read %d individuals x %d loci\n", nInd(genind_obj), nLoc(genind_obj)))

popmap <- read.table(
  popmap_file,
  header = FALSE,
  sep = "\t",
  col.names = c("individual", "population"),
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)

ind_names <- indNames(genind_obj)
idx <- match(ind_names, popmap$individual)
if (anyNA(idx)) {
  stop(
    "Samples present in the VCF but missing from the popmap: ",
    paste(ind_names[is.na(idx)], collapse = ", ")
  )
}
pop(genind_obj) <- factor(popmap$population[idx])

pop_sizes <- table(pop(genind_obj))
too_small <- names(pop_sizes)[pop_sizes < min_individuals]
if (length(too_small) > 0) {
  message(sprintf(
    "Dropping %d population(s) with < %d individuals: %s\n",
    length(too_small), min_individuals, paste(too_small, collapse = ", ")
  ))
  keep <- !(as.character(pop(genind_obj)) %in% too_small)
  genind_obj <- genind_obj[keep, ]
}
if (nPop(genind_obj) < 2) {
  stop("WC popgenstats require at least 2 populations with >= ",
       min_individuals, " individuals")
}
message(sprintf("Analysing %d individuals in %d populations (grouping: %s)\n",
                nInd(genind_obj), nPop(genind_obj), grouping_name))

message("\n=== CONVERTING TO hierfstat FORMAT ===\n")
data_hf <- genind2hierfstat(genind_obj)

message("\n=== basic.stats (Weir & Cockerham) ===\n")
bstats <- basic.stats(data_hf)
overall_vec <- bstats$overall
overall_df <- data.frame(
  statistic = names(overall_vec),
  value = as.numeric(overall_vec),
  stringsAsFactors = FALSE
)
message(sprintf("Overall Fis: %s\n", overall_vec["Fis"]))
message(sprintf("Overall Fst: %s\n", overall_vec["Fst"]))

pop_levels <- levels(pop(genind_obj))
pop_sizes <- table(pop(genind_obj))

ho_means <- colMeans(bstats$Ho, na.rm = TRUE)
he_means <- colMeans(bstats$Hs, na.rm = TRUE)
fis_means <- colMeans(bstats$Fis, na.rm = TRUE)

message("\n=== Allelic richness ===\n")
ar <- allelic.richness(data_hf)
ar_means <- colMeans(ar$Ar, na.rm = TRUE)
message(sprintf("Rarefaction min.n = %s\n", paste(ar$min.n, collapse = ", ")))

message("\n=== Private alleles ===\n")
pa <- private_alleles(genind_obj, count.alleles = TRUE)
private_counts <- rowSums(pa > 0)
# Align to population factor levels (private_alleles rows = pops)
private_aligned <- setNames(rep(0, length(pop_levels)), pop_levels)
private_aligned[names(private_counts)] <- as.numeric(private_counts)

popstats <- data.frame(
  population = pop_levels,
  n_individuals = as.integer(pop_sizes[pop_levels]),
  Ho = as.numeric(ho_means[pop_levels]),
  He = as.numeric(he_means[pop_levels]),
  Fis = as.numeric(fis_means[pop_levels]),
  Ar = as.numeric(ar_means[pop_levels]),
  private_alleles = as.numeric(private_aligned[pop_levels]),
  stringsAsFactors = FALSE
)
popstats <- popstats[order(popstats$population), , drop = FALSE]
message("Per-population summary:\n")
print(popstats, row.names = FALSE)

message("\n=== Pairwise WC Fst ===\n")
pw_fst <- pairwise.WCfst(data_hf)
pairwise_df <- matrix_to_long(pw_fst, "fst")
message(sprintf("Wrote %d pairwise comparisons\n", nrow(pairwise_df)))

message(sprintf("\n=== Bootstrap pairwise Fst CIs (nboot = %d) ===\n", nboot))
boot_fst <- boot.ppfst(data_hf, nboot = nboot)
ll_long <- matrix_to_long(boot_fst$ll, "ll")
ul_long <- matrix_to_long(boot_fst$ul, "ul")
ci_df <- merge(ll_long, ul_long, by = c("pop1", "pop2"), all = TRUE)
ci_df <- ci_df[order(ci_df$pop1, ci_df$pop2), , drop = FALSE]

message("\n=== WRITING OUTPUTS ===\n")
for (path in c(out_overall, out_popstats, out_pairwise, out_ci)) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}
write.table(overall_df, out_overall, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(popstats, out_popstats, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(pairwise_df, out_pairwise, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(ci_df, out_ci, sep = "\t", quote = FALSE, row.names = FALSE)

message("\n=== COMPLETED SUCCESSFULLY ===\n")

sink(type = "message")
sink(type = "output")
close(log_file)
