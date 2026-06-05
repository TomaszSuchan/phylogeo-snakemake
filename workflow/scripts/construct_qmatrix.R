#!/usr/bin/env Rscript
# Export headerless layer-proportion matrix for mapmixture plotting.

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

layer_props <- read.table(
  snakemake@input[["layer_proportions"]],
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE
)

qmat <- as.matrix(layer_props)
if (any(is.na(qmat))) {
  stop("Layer proportions contain NA values; cannot write Q matrix for mapmixture.")
}

write.table(
  qmat,
  file = snakemake@output[["qmatrix"]],
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

cat("Wrote mapmixture Q matrix:", nrow(qmat), "individuals x", ncol(qmat), "layers\n")
