#!/usr/bin/env Rscript
# Population-specific Fst and Fis (Weir & Goudet 2017) plus Nei & Chesser (1983)
# Ho/Hs, computed from a biallelic-SNP dosage matrix.

log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

suppressPackageStartupMessages({
  library(gaston)
  library(hierfstat)
})

vcf_file <- snakemake@input[["vcf"]]
popmap_file <- snakemake@input[["popmap"]]
output_stats <- snakemake@output[["stats"]]

grouping_name <- as.character(snakemake@params[["grouping"]])
min_individuals <- suppressWarnings(as.integer(snakemake@params[["min_individuals"]]))
if (is.na(min_individuals)) min_individuals <- 2

message("\n=== READING GENOTYPES ===\n")
message(sprintf("VCF: %s\n", vcf_file))
bm <- hierfstat::read.VCF(vcf_file)
dos <- gaston::as.matrix(bm)
rownames(dos) <- bm@ped$id
message(sprintf("Read %d individuals x %d SNPs\n", nrow(dos), ncol(dos)))

popmap <- read.table(
  popmap_file,
  header = FALSE,
  sep = "\t",
  col.names = c("individual", "population"),
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)

idx <- match(rownames(dos), popmap$individual)
if (anyNA(idx)) {
  stop(
    "Samples present in the VCF but missing from the popmap: ",
    paste(rownames(dos)[is.na(idx)], collapse = ", ")
  )
}
pop <- popmap$population[idx]

pop_sizes <- table(pop)
too_small <- names(pop_sizes)[pop_sizes < min_individuals]
if (length(too_small) > 0) {
  message(sprintf(
    "Dropping %d population(s) with < %d individuals: %s\n",
    length(too_small), min_individuals, paste(too_small, collapse = ", ")
  ))
  keep <- !(pop %in% too_small)
  dos <- dos[keep, , drop = FALSE]
  pop <- pop[keep]
}
if (length(unique(pop)) < 2) {
  stop("Population-specific Fst requires at least 2 populations with >= ",
       min_individuals, " individuals")
}
message(sprintf("Analysing %d individuals in %d populations\n",
                nrow(dos), length(unique(pop))))

message("\n=== POPULATION-SPECIFIC Fst / Fis (Weir & Goudet 2017) ===\n")
fs <- hierfstat::fs.dosage(dos, pop = pop)
fs_mat <- fs$Fs
if (is.null(fs_mat) || !all(c("Fis", "Fst") %in% rownames(fs_mat))) {
  stop("Unexpected fs.dosage() output: missing Fis/Fst rows")
}
message("Overall (pooled) estimates:\n")
message(sprintf("  Fis = %.5f\n", fs_mat["Fis", "All"]))
message(sprintf("  Fst = %.5f\n", fs_mat["Fst", "All"]))

message("\n=== WITHIN-POPULATION HETEROZYGOSITY ===\n")
# Per locus Ho and Nei & Chesser (1983) sample-size corrected Hs, averaged over
# loci genotyped in at least two individuals of the population.
pop_heterozygosity <- function(sub) {
  n <- colSums(!is.na(sub))
  ho <- colSums(sub == 1, na.rm = TRUE) / n
  p <- colSums(sub, na.rm = TRUE) / (2 * n)
  hs <- (n / (n - 1)) * (1 - (p^2 + (1 - p)^2) - ho / (2 * n))
  usable <- n >= 2
  c(
    n_loci = sum(usable),
    Ho = mean(ho[usable]),
    Hs = mean(hs[usable])
  )
}

pops <- setdiff(colnames(fs_mat), "All")
het <- vapply(pops, function(p) pop_heterozygosity(dos[pop == p, , drop = FALSE]),
              numeric(3))

stats <- data.frame(
  population = pops,
  n_individuals = as.integer(pop_sizes[pops]),
  n_loci = as.integer(het["n_loci", ]),
  Ho = het["Ho", ],
  Hs = het["Hs", ],
  Fis = as.numeric(fs_mat["Fis", pops]),
  Fst = as.numeric(fs_mat["Fst", pops]),
  stringsAsFactors = FALSE
)
stats <- stats[order(stats$population), , drop = FALSE]

message(sprintf("\n=== WRITING %s (grouping: %s) ===\n", output_stats, grouping_name))
print(stats, row.names = FALSE)
dir.create(dirname(output_stats), recursive = TRUE, showWarnings = FALSE)
write.table(stats, output_stats, sep = "\t", quote = FALSE, row.names = FALSE)

message("\n=== COMPLETED SUCCESSFULLY ===\n")

sink(type = "message")
sink(type = "output")
close(log_file)
