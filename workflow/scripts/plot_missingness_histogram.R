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
  title <- "Individual Missingness Distribution"
} else {
  missing_col <- "F_MISS"
  x_label <- "Proportion of missing data per locus"
  title <- "Locus Missingness Distribution"
}

message(sprintf("Missingness column: %s\n", missing_col))
message(sprintf("Missingness range: %.4f to %.4f\n", 
                min(missing_df[[missing_col]], na.rm = TRUE),
                max(missing_df[[missing_col]], na.rm = TRUE)))

# Create histogram
message("\n=== CREATING HISTOGRAM ===\n")
p <- ggplot(missing_df, aes(x = .data[[missing_col]])) +
  geom_histogram(bins = 50, fill = "steelblue", color = "black", alpha = 0.7) +
  labs(
    x = x_label,
    y = "Frequency",
    title = title
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

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

