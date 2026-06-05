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

if (!requireNamespace("tess3r", quietly = TRUE)) {
  stop(
    "tess3r is not installed in this rule environment. The install_tess3 ",
    "Snakemake rule should install it before analysis."
  )
}
suppressPackageStartupMessages({
  library(tess3r)
})

geno_rds <- snakemake@input[["geno_rds"]]
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
ploidy <- as.integer(snakemake@params[["ploidy"]])
map_method <- as.character(snakemake@params[["map_method"]])
map_resolution <- as.integer(unlist(snakemake@params[["map_resolution"]]))
interpolation_knots <- as.integer(snakemake@params[["interpolation_knots"]])
threads <- as.integer(snakemake@threads)

cat("=== tess3r Analysis ===\n")
cat("K:", k, "\n")
cat("method:", method, "\n")
cat("ploidy:", ploidy, "\n")
cat("replicates:", replicates, "\n")
cat("max_iteration:", max_iteration, "\n")
cat("tolerance:", tolerance, "\n")
cat("threads:", threads, "\n\n")

dir.create(dirname(qmatrix_file), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(max_cluster_png), recursive = TRUE, showWarnings = FALSE)

geno_data <- readRDS(geno_rds)
geno <- geno_data$geno
sample_names <- geno_data$samples

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

indpopdata_ordered <- indpopdata[match(sample_names, indpopdata$Ind), , drop = FALSE]
coords <- as.matrix(indpopdata_ordered[, c("Lon", "Lat")])
storage.mode(coords) <- "numeric"
rownames(coords) <- sample_names

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

qmat_obj <- qmatrix(tess3_obj, K = k)
qmat <- as.matrix(qmat_obj)
if (nrow(qmat) != length(sample_names) && ncol(qmat) == length(sample_names)) {
  qmat <- t(qmat)
  qmat_obj <- t(qmat_obj)
}
if (nrow(qmat) != length(sample_names)) {
  stop("tess3r Q matrix has ", nrow(qmat), " rows, but expected ", length(sample_names))
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
      threads = threads,
      map_method = map_method,
      map_resolution = map_resolution,
      interpolation_knots = interpolation_knots
    )
  ),
  file = results_rds
)
cat("RDS written to:", results_rds, "\n")

cat("Drawing tess3r native max-cluster map...\n")
window <- c(
  min(coords[, 1], na.rm = TRUE),
  max(coords[, 1], na.rm = TRUE),
  min(coords[, 2], na.rm = TRUE),
  max(coords[, 2], na.rm = TRUE)
)
tryCatch({
  png(max_cluster_png, width = 1800, height = 1400, res = 200)
  plot_tess3Q_fn <- getS3method("plot", "tess3Q")
  plot_args <- list(
    x = qmat_obj,
    coord = coords,
    method = map_method,
    resolution = map_resolution,
    window = window,
    col.palette = CreatePalette(),
    xlab = "Longitude",
    ylab = "Latitude",
    main = "",
    cex = 0.4
  )
  if ("interpol" %in% names(formals(plot_tess3Q_fn))) {
    plot_args$interpol <- FieldsKrigModel(interpolation_knots)
  } else {
    plot_args$interpolation.model <- FieldsKrigModel(interpolation_knots)
  }
  do.call(plot_tess3Q_fn, plot_args)
  dev.off()
}, error = function(e) {
  cat("WARNING: tess3r native plot failed:", conditionMessage(e), "\n")
  if (dev.cur() != 1) dev.off()
  png(max_cluster_png, width = 1800, height = 1400, res = 200)
  plot.new()
  text(0.5, 0.5, paste("Max-cluster plot failed:", conditionMessage(e)))
  dev.off()
})

cat("=== tess3r Analysis Complete ===\n")
