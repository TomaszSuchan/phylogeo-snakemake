#!/usr/bin/env Rscript
# Plot histogram of sequencing depth from vcftools .idepth or .ldepth.mean output.

library(ggplot2)

histogram_utils <- tryCatch({
  file.path(dirname(normalizePath(snakemake@script)), "plot_histogram_utils.R")
}, error = function(e) "workflow/scripts/plot_histogram_utils.R")
if (file.exists(histogram_utils)) {
  source(histogram_utils)
} else {
  source("workflow/scripts/plot_histogram_utils.R")
}

pdf(NULL)

log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")
on.exit({
  while (sink.number(type = "message") > 0) sink(type = "message")
  while (sink.number(type = "output") > 0) sink(type = "output")
  close(log_file)
}, add = TRUE)

if ("idepth" %in% names(snakemake@input)) {
  input_file <- snakemake@input[["idepth"]]
  depth_type <- "individual"
  x_label <- "Mean sequencing depth per individual"
} else if ("ldepth" %in% names(snakemake@input)) {
  input_file <- snakemake@input[["ldepth"]]
  depth_type <- "locus"
  x_label <- "Mean sequencing depth per locus"
} else {
  input_file <- snakemake@input[[1]]
  depth_type <- ifelse(grepl("\\.idepth$", input_file), "individual", "locus")
  x_label <- ifelse(
    depth_type == "individual",
    "Mean sequencing depth per individual",
    "Mean sequencing depth per locus"
  )
}

output_pdf <- snakemake@output[["pdf"]]
output_rds <- snakemake@output[["rds"]]
summary_path <- snakemake@output[["summary"]]

message(sprintf("\n=== READING %s DEPTH DATA ===\n", toupper(depth_type)))
depth_df <- read.table(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
if (!"MEAN_DEPTH" %in% names(depth_df)) {
  stop("Depth input must contain a MEAN_DEPTH column. Found: ", paste(names(depth_df), collapse = ", "))
}

depth <- depth_df$MEAN_DEPTH
depth <- depth[is.finite(depth)]
message(sprintf("Loaded %d finite depth values\n", length(depth)))
message(sprintf("Depth range: %.4f to %.4f\n", min(depth), max(depth)))

message("\n=== CREATING HISTOGRAM ===\n")
p <- ggplot(data.frame(mean_depth = depth), aes(x = mean_depth)) +
  geom_histogram_styled() +
  labs(
    x = x_label,
    y = "Frequency"
  ) +
  histogram_plot_theme()

message("\n=== CALCULATING SUMMARY STATISTICS ===\n")
n <- length(depth)
mean_x <- mean(depth)
median_x <- median(depth)
sd_x <- sd(depth)
min_x <- min(depth)
max_x <- max(depth)
quantiles <- quantile(depth, probs = c(0.05, 0.25, 0.75, 0.95), names = FALSE)

dir.create(dirname(summary_path), recursive = TRUE, showWarnings = FALSE)
con <- file(summary_path, open = "wt")
writeLines(sprintf("Depth summary (%s level)", depth_type), con)
writeLines(sprintf("N: %d", n), con)
writeLines(sprintf("Mean depth: %.6f", mean_x), con)
writeLines(sprintf("Median depth: %.6f", median_x), con)
writeLines(sprintf("SD depth: %.6f", sd_x), con)
writeLines(sprintf("Range: %.6f-%.6f", min_x, max_x), con)
writeLines(sprintf("P05 depth: %.6f", quantiles[1]), con)
writeLines(sprintf("P25 depth: %.6f", quantiles[2]), con)
writeLines(sprintf("P75 depth: %.6f", quantiles[3]), con)
writeLines(sprintf("P95 depth: %.6f", quantiles[4]), con)
close(con)

message("\n=== SAVING OUTPUT ===\n")
message(sprintf("Output PDF: %s\n", output_pdf))
message(sprintf("Output RDS: %s\n", output_rds))

dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)
ggsave_pdf(
  filename = output_pdf,
  plot = p,
  width = 10,
  height = 6,
  dpi = 300
)

dir.create(dirname(output_rds), recursive = TRUE, showWarnings = FALSE)
saveRDS(p, output_rds)

message("\n=== COMPLETED SUCCESSFULLY ===\n")
