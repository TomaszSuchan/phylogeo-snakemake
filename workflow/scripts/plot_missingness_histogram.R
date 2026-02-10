#!/usr/bin/env Rscript
# Plot histogram of missingness (individual or locus level)

library(ggplot2)

# Prevent creation of Rplots.pdf
pdf(NULL)

# Redirect all output to log file
log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

# Determine if this is imiss or lmiss based on input file
# Check which input was provided
if ("imiss" %in% names(snakemake@input)) {
  input_file <- snakemake@input[["imiss"]]
  is_imiss <- TRUE
} else if ("lmiss" %in% names(snakemake@input)) {
  input_file <- snakemake@input[["lmiss"]]
  is_imiss <- FALSE
} else {
  # Fallback to first input
  input_file <- snakemake@input[[1]]
  is_imiss <- grepl("\\.imiss$", input_file)
}

output_pdf <- snakemake@output[["pdf"]]
output_rds <- snakemake@output[["rds"]]

data_type <- ifelse(is_imiss, "Individual", "Locus")

message(sprintf("\n=== READING %s MISSINGNESS DATA ===\n", toupper(data_type)))
missing_df <- read.table(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
message(sprintf("Loaded %d %s\n", nrow(missing_df), tolower(data_type)))

# Determine missingness column name
if (is_imiss) {
  missing_col <- "F_MISS"
  x_label <- "Proportion of missing data per individual"
} else {
  missing_col <- "F_MISS"
  x_label <- "Proportion of missing data per locus"
}

message(sprintf("Missingness column: %s\n", missing_col))
message(sprintf("Missingness range: %.4f to %.4f\n", 
                min(missing_df[[missing_col]], na.rm = TRUE),
                max(missing_df[[missing_col]], na.rm = TRUE)))

# Create histogram (no title)
message("\n=== CREATING HISTOGRAM ===\n")
p <- ggplot(missing_df, aes(x = .data[[missing_col]])) +
  geom_histogram(bins = 50, fill = "steelblue", color = "black", alpha = 0.7) +
  labs(
    x = x_label,
    y = "Frequency"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

# Create summary statistics and save to text file (if output defined)
summary_path <- NULL
if (!is.null(snakemake@output[["summary"]])) {
  summary_path <- snakemake@output[["summary"]]
}

if (!is.null(summary_path)) {
  message("\n=== CALCULATING SUMMARY STATISTICS ===\n")
  x <- missing_df[[missing_col]]
  x <- x[is.finite(x)]

  n <- length(x)
  mean_x <- mean(x)
  median_x <- median(x)
  sd_x <- sd(x)
  min_x <- min(x)
  max_x <- max(x)

  # Bin into 10% intervals over [0, 1]: [0.0–0.1], (0.1–0.2], ..., (0.9–1.0]
  breaks <- seq(0, 1, by = 0.1)
  bin_counts <- integer(length(breaks) - 1)
  bin_perc <- numeric(length(breaks) - 1)
  bin_cum_perc <- numeric(length(breaks) - 1)
  bin_labels <- character(length(breaks) - 1)

  for (i in seq_len(length(breaks) - 1)) {
    lower <- breaks[i]
    upper <- breaks[i + 1]
    # First bin includes lower bound, others are (lower, upper]
    if (i == 1) {
      in_bin <- x >= lower & x <= upper
    } else {
      in_bin <- x > lower & x <= upper
    }
    bin_counts[i] <- sum(in_bin)
    bin_perc[i] <- if (n > 0) 100 * bin_counts[i] / n else 0
    bin_labels[i] <- sprintf("(%.1f-%.1f]", lower, upper)
  }
  # Cumulative percentage across bins
  if (n > 0) {
    bin_cum_perc <- cumsum(bin_perc)
  } else {
    bin_cum_perc[] <- 0
  }

  dir.create(dirname(summary_path), recursive = TRUE, showWarnings = FALSE)
  con <- file(summary_path, open = "wt")

  writeLines(sprintf("Missingness summary (%s level)", tolower(data_type)), con)
  writeLines(sprintf("N: %d", n), con)
  writeLines(sprintf("Mean missingness: %.6f (%.2f%%)", mean_x, 100 * mean_x), con)
  writeLines(sprintf("Median missingness: %.6f (%.2f%%)", median_x, 100 * median_x), con)
  writeLines(sprintf("SD missingness: %.6f (%.2f%%)", sd_x, 100 * sd_x), con)
  writeLines(sprintf("Range: %.6f-%.6f (%.2f%%-%.2f%%)", min_x, max_x, 100 * min_x, 100 * max_x), con)
  writeLines("", con)
  writeLines("Bin\tCount\tPercentage_of_samples\tCumulative_percentage", con)

  for (i in seq_along(bin_labels)) {
    writeLines(sprintf("%s\t%d\t%.2f\t%.2f",
                       bin_labels[i],
                       bin_counts[i],
                       bin_perc[i],
                       bin_cum_perc[i]),
               con)
  }

  close(con)
}

# Save plots
message("\n=== SAVING OUTPUT ===\n")
message(sprintf("Output PDF: %s\n", output_pdf))
message(sprintf("Output RDS: %s\n", output_rds))

dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)
ggsave(
  filename = output_pdf,
  plot = p,
  width = 10,
  height = 6,
  dpi = 300,
  device = "pdf"
)

dir.create(dirname(output_rds), recursive = TRUE, showWarnings = FALSE)
saveRDS(p, output_rds)

message("\n=== COMPLETED SUCCESSFULLY ===\n")

# Close log file sinks
sink(type = "message")
sink(type = "output")
close(log_file)

