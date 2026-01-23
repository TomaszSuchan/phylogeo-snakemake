#!/usr/bin/env Rscript
# Create lmiss barplot (plain, no grouping)

library(ggplot2)

# Prevent creation of Rplots.pdf
pdf(NULL)

# Redirect all output to log file
log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

# Snakemake inputs/outputs
lmiss_file <- snakemake@input[["lmiss"]]
output_pdf <- snakemake@output[["pdf"]]
output_rds <- snakemake@output[["rds"]]

message("\n=== READING LMISS DATA ===\n")
lmiss_df <- read.table(lmiss_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
message(sprintf("Loaded %d loci\n", nrow(lmiss_df)))

# Sort by missingness
lmiss_df <- lmiss_df[order(lmiss_df$F_MISS), ]
# Create a simple index for x-axis (too many loci to show individual labels)
lmiss_df$locus_index <- 1:nrow(lmiss_df)

message("\n=== CREATING BARPLOT ===\n")
message(sprintf("Plotting %d loci\n", nrow(lmiss_df)))

# Create barplot
p <- ggplot(lmiss_df, aes(x = locus_index, y = F_MISS)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7, color = "black", linewidth = 0.1) +
  labs(
    x = "Locus (sorted by missingness)",
    y = "Proportion of missing data",
    title = "Locus Missingness Distribution"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),  # Too many loci to show labels
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  )

# Save plots
message("\n=== SAVING OUTPUT ===\n")
message(sprintf("Output PDF: %s\n", output_pdf))
message(sprintf("Output RDS: %s\n", output_rds))

dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)
ggsave(
  filename = output_pdf,
  plot = p,
  width = 12,
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

