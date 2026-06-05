#!/usr/bin/env Rscript
# Run tess3r across all K values and summarise model scores (plot.tess3 logic).

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
  stop("tess3r is not installed in this rule environment.")
}
suppressPackageStartupMessages({
  library(tess3r)
})

geno_rds <- snakemake@input[["geno_rds"]]
indpopdata_file <- snakemake@input[["indpopdata"]]
results_rds <- snakemake@output[["results_rds"]]
choose_k_file <- snakemake@output[["choose_k_results"]]
cv_summary_file <- snakemake@output[["cv_summary"]]

k_values <- as.integer(unlist(snakemake@params[["k_values"]]))
method <- as.character(snakemake@params[["method"]])
replicates <- as.integer(snakemake@params[["replicates"]])
max_iteration <- as.integer(snakemake@params[["max_iteration"]])
tolerance <- as.numeric(snakemake@params[["tolerance"]])
ploidy <- as.integer(snakemake@params[["ploidy"]])
mask <- as.numeric(snakemake@params[["mask"]])
crossvalid <- isTRUE(snakemake@params[["crossvalid"]])
crossentropy <- isTRUE(snakemake@params[["crossentropy"]])
threads <- as.integer(snakemake@threads)

if (crossvalid && (!is.finite(mask) || mask <= 0)) {
  stop(
    "tess3.crossvalid=TRUE requires tess3.mask > 0 ",
    "(tess3r package default mask=0 disables masked cross-validation)."
  )
}

cat("=== tess3r choose-K ===\n")
cat("K values:", paste(k_values, collapse = ", "), "\n")
cat("method:", method, "\n")
cat("replicates:", replicates, "\n")
cat("max_iteration:", max_iteration, "\n")
cat("mask:", mask, "\n")
cat("crossvalid:", crossvalid, "\n")
cat("crossentropy:", crossentropy, "\n")
cat("ploidy:", ploidy, "\n")
cat("threads:", threads, "\n")

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
indpopdata_ordered <- indpopdata[match(sample_names, indpopdata$Ind), , drop = FALSE]
coords <- as.matrix(indpopdata_ordered[, c("Lon", "Lat")])
storage.mode(coords) <- "numeric"
rownames(coords) <- sample_names

cat("Running tess3r across K...\n")
tess3_args <- list(
  X = geno,
  coord = coords,
  K = k_values,
  ploidy = ploidy,
  method = method,
  rep = replicates,
  max.iteration = max_iteration,
  tolerance = tolerance,
  openMP.core.num = threads
)
if (is.finite(mask) && mask > 0) {
  tess3_args$mask <- mask
}
tess3_obj <- do.call(tess3, tess3_args)

extract_scores <- function(obj, crossvalid = FALSE, crossentropy = FALSE) {
  med <- min <- max <- K <- numeric(length(obj))
  for (i in seq_along(obj)) {
    K[i] <- obj[[i]]$K
    if (!crossentropy) {
      if (!crossvalid) {
        med[i] <- median(obj[[i]]$rmse)
        min[i] <- min(obj[[i]]$rmse)
        max[i] <- max(obj[[i]]$rmse)
      } else {
        med[i] <- median(obj[[i]]$crossvalid.rmse)
        min[i] <- min(obj[[i]]$crossvalid.rmse)
        max[i] <- max(obj[[i]]$crossvalid.rmse)
      }
    } else if (!crossvalid) {
      med[i] <- median(obj[[i]]$crossentropy)
      min[i] <- min(obj[[i]]$crossentropy)
      max[i] <- max(obj[[i]]$crossentropy)
    } else {
      med[i] <- median(obj[[i]]$crossvalid.crossentropy)
      min[i] <- min(obj[[i]]$crossvalid.crossentropy)
      max[i] <- max(obj[[i]]$crossvalid.crossentropy)
    }
  }
  data.frame(K = K, median = med, min = min, max = max)
}

cv_summary <- extract_scores(
  tess3_obj,
  crossvalid = crossvalid,
  crossentropy = crossentropy
)

write.table(
  cv_summary,
  file = cv_summary_file,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

best_k <- cv_summary$K[which.min(cv_summary$median)]
score_label <- if (crossentropy) {
  if (crossvalid) "cross-validation cross-entropy" else "cross-entropy"
} else if (crossvalid) {
  "cross-validation RMSE"
} else {
  "RMSE"
}

writeLines(
  c(
    paste0("Choose K with the lowest ", score_label, " across replicate runs."),
    "Smaller values indicate better fit; look for a plateau or minimum.",
    paste0("Suggested K (minimum median score): ", best_k),
    "",
    paste0("K\tmedian_", gsub(" ", "_", score_label), "\tmin\tmax"),
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
    tess3 = tess3_obj,
    cv_summary = cv_summary,
    samples = sample_names,
    coordinates = coords,
    params = list(
      K = k_values,
      method = method,
      replicates = replicates,
      max_iteration = max_iteration,
      tolerance = tolerance,
      ploidy = ploidy,
      mask = mask,
      crossvalid = crossvalid,
      crossentropy = crossentropy,
      threads = threads
    )
  ),
  file = results_rds
)

cat("Score summary written to:", cv_summary_file, "\n")
cat("choose-K results written to:", choose_k_file, "\n")
cat("Suggested K:", best_k, "\n")
