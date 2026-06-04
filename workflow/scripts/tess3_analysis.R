#!/usr/bin/env Rscript
# tess3r spatial ancestry analysis for one K.

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

if (!requireNamespace("tess3r", quietly = TRUE)) {
  stop(
    "tess3r is not installed in this rule environment. The install_tess3 ",
    "Snakemake rule should install it before analysis."
  )
}
suppressPackageStartupMessages({
  library(tess3r)
})

vcf_file <- snakemake@input[["vcf"]]
indpopdata_file <- snakemake@input[["indpopdata"]]
qmatrix_file <- snakemake@output[["qmatrix"]]
results_rds <- snakemake@output[["results_rds"]]
cross_entropy_file <- snakemake@output[["cross_entropy"]]
max_cluster_png <- snakemake@output[["max_cluster_png"]]

k <- as.integer(snakemake@params[["k"]])
method <- as.character(snakemake@params[["method"]])
replicates <- as.integer(snakemake@params[["replicates"]])
max_iteration <- as.integer(snakemake@params[["max_iteration"]])
tolerance <- as.numeric(snakemake@params[["tolerance"]])
n_colors <- as.integer(snakemake@params[["n_colors"]])
ploidy <- as.integer(snakemake@params[["ploidy"]])
threads <- as.integer(snakemake@threads)

cat("=== tess3r Analysis ===\n")
cat("VCF:", vcf_file, "\n")
cat("indpopdata:", indpopdata_file, "\n")
cat("K:", k, "\n")
cat("method:", method, "\n")
cat("replicates:", replicates, "\n")
cat("max_iteration:", max_iteration, "\n")
cat("tolerance:", tolerance, "\n")
cat("threads:", threads, "\n\n")

dir.create(dirname(qmatrix_file), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(max_cluster_png), recursive = TRUE, showWarnings = FALSE)

indpopdata <- read.table(
  indpopdata_file,
  header = TRUE,
  sep = "\t",
  comment.char = "",
  fill = TRUE,
  blank.lines.skip = TRUE
)
required_cols <- c("Ind", "Lat", "Lon")
missing_cols <- setdiff(required_cols, colnames(indpopdata))
if (length(missing_cols) > 0) {
  stop("indpopdata is missing required columns: ", paste(missing_cols, collapse = ", "))
}
if (any(is.na(indpopdata$Lat) | is.na(indpopdata$Lon))) {
  stop("Missing Lat/Lon values in indpopdata; tess3r requires coordinates for every sample.")
}

cat("Reading VCF...\n")
vcf <- read.vcfR(vcf_file, verbose = FALSE)
sample_names <- colnames(vcf@gt)[-1]
cat("VCF samples:", length(sample_names), "\n")
cat("VCF variants:", nrow(vcf@fix), "\n")

if (!all(sample_names %in% indpopdata$Ind)) {
  missing_inds <- sample_names[!sample_names %in% indpopdata$Ind]
  stop(
    "Samples in VCF are missing from indpopdata: ",
    paste(head(missing_inds, 20), collapse = ", "),
    if (length(missing_inds) > 20) " ..."
  )
}

indpopdata_ordered <- indpopdata[match(sample_names, indpopdata$Ind), , drop = FALSE]
coords <- as.matrix(indpopdata_ordered[, c("Lon", "Lat")])
storage.mode(coords) <- "numeric"
rownames(coords) <- sample_names

gt <- extract.gt(vcf, element = "GT", as.numeric = FALSE)
cat("Converting genotypes to alternate-allele dosages...\n")
geno_dosage <- apply(gt, c(1, 2), function(x) {
  if (is.na(x) || x %in% c(".", "./.", ".|.")) {
    return(9)
  }
  alleles <- unlist(strsplit(x, "[/|]"))
  if (length(alleles) == 0 || any(alleles == ".")) {
    return(9)
  }
  sum(as.integer(alleles) > 0)
})
geno <- t(geno_dosage)
storage.mode(geno) <- "numeric"
rownames(geno) <- sample_names
colnames(geno) <- rownames(gt)

cat("Genotype matrix:", nrow(geno), "individuals x", ncol(geno), "loci\n")
cat("Coordinate matrix:", nrow(coords), "individuals x", ncol(coords), "columns\n")

cat("Running tess3r...\n")
tess3_obj <- tess3(
  X = geno,
  coord = coords,
  K = k,
  ploidy = ploidy,
  method = method,
  rep = replicates,
  max.iteration = max_iteration,
  tolerance = tolerance,
  openMP.core.num = threads
)

qmat <- qmatrix(tess3_obj, K = k)
qmat <- as.matrix(qmat)
if (nrow(qmat) != length(sample_names) && ncol(qmat) == length(sample_names)) {
  qmat <- t(qmat)
}
if (nrow(qmat) != length(sample_names)) {
  stop("tess3r Q matrix has ", nrow(qmat), " rows, but expected ", length(sample_names))
}
colnames(qmat) <- paste0("Cluster", seq_len(ncol(qmat)))
rownames(qmat) <- sample_names

# Shared map/barplot scripts read headerless Q matrices in filtered VCF sample order.
write.table(
  qmat,
  file = qmatrix_file,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)
cat("Q matrix written to:", qmatrix_file, "\n")

cross_entropy <- tryCatch(
  cross.entropy(tess3_obj, K = k),
  error = function(e) {
    cat("Could not extract cross-entropy:", conditionMessage(e), "\n")
    NA_real_
  }
)
write.table(
  data.frame(K = k, cross_entropy = cross_entropy),
  file = cross_entropy_file,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
cat("Cross-entropy summary written to:", cross_entropy_file, "\n")

saveRDS(
  list(
    tess3 = tess3_obj,
    qmatrix = qmat,
    samples = sample_names,
    coordinates = coords,
    params = list(
      K = k,
      method = method,
      replicates = replicates,
      max_iteration = max_iteration,
      tolerance = tolerance,
      ploidy = ploidy,
      threads = threads
    )
  ),
  file = results_rds
)
cat("RDS written to:", results_rds, "\n")

cat("Drawing tess3r native max-cluster map...\n")
tryCatch({
  png(max_cluster_png, width = 1800, height = 1400, res = 200)
  plot(
    tess3_obj,
    K = k,
    method = "map.max",
    resolution = c(300, 300),
    window = c(min(coords[, 1]), max(coords[, 1]), min(coords[, 2]), max(coords[, 2])),
    col.palette = CreatePalette(n_colors)
  )
  dev.off()
}, error = function(e) {
  cat("WARNING: tess3r native plot failed:", conditionMessage(e), "\n")
  if (dev.cur() != 1) dev.off()
  file.create(max_cluster_png)
})

cat("=== tess3r Analysis Complete ===\n")
