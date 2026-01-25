#!/usr/bin/env Rscript
# Create boxplots of pixy values with bootstrapping
# Uses bootstrap confidence intervals from summary files

library(tidyverse)

# Prevent creation of Rplots.pdf
pdf(NULL)

# Redirect all output to log file
log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

# Snakemake inputs/outputs
stat_type <- snakemake@params[["stat_type"]]  # "pi", "dxy", or "fst"
summary_file <- snakemake@input[["summary"]]
output_pdf <- snakemake@output[["pdf"]]
output_rds <- snakemake@output[["rds"]]

message("\n=== READING SUMMARY FILE ===\n")
summary_df <- read.table(summary_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
message(sprintf("Loaded %d rows from summary file\n", nrow(summary_df)))

# Prepare data based on stat_type
if (stat_type == "pi") {
  # For pi: one value per population
  plot_df <- summary_df %>%
    select(population, mean_pi, ci_low, ci_high, bootstrap_se) %>%
    rename(
      Population = population,
      Value = mean_pi,
      CI_low = ci_low,
      CI_high = ci_high,
      SE = bootstrap_se
    )
  plot_title <- "Nucleotide Diversity (π) by Population"
  y_label <- "π"
  
} else if (stat_type == "dxy") {
  # For dxy: pairwise comparisons
  # Create a long format with both pop1 and pop2
  plot_df <- bind_rows(
    summary_df %>%
      select(pop1, mean_dxy, ci_low, ci_high, bootstrap_se) %>%
      rename(Population = pop1, Value = mean_dxy, CI_low = ci_low, CI_high = ci_high, SE = bootstrap_se),
    summary_df %>%
      select(pop2, mean_dxy, ci_low, ci_high, bootstrap_se) %>%
      rename(Population = pop2, Value = mean_dxy, CI_low = ci_low, CI_high = ci_high, SE = bootstrap_se)
  ) %>%
    group_by(Population) %>%
    summarise(
      Value = mean(Value, na.rm = TRUE),
      CI_low = mean(CI_low, na.rm = TRUE),
      CI_high = mean(CI_high, na.rm = TRUE),
      SE = mean(SE, na.rm = TRUE),
      .groups = "drop"
    )
  plot_title <- "Nucleotide Divergence (Dxy) by Population"
  y_label <- "Dxy"
  
} else if (stat_type == "fst") {
  # For fst: pairwise comparisons
  # Create a long format with both pop1 and pop2
  plot_df <- bind_rows(
    summary_df %>%
      select(pop1, mean_fst, ci_low, ci_high, bootstrap_se) %>%
      rename(Population = pop1, Value = mean_fst, CI_low = ci_low, CI_high = ci_high, SE = bootstrap_se),
    summary_df %>%
      select(pop2, mean_fst, ci_low, ci_high, bootstrap_se) %>%
      rename(Population = pop2, Value = mean_fst, CI_low = ci_low, CI_high = ci_high, SE = bootstrap_se)
  ) %>%
    group_by(Population) %>%
    summarise(
      Value = mean(Value, na.rm = TRUE),
      CI_low = mean(CI_low, na.rm = TRUE),
      CI_high = mean(CI_high, na.rm = TRUE),
      SE = mean(SE, na.rm = TRUE),
      .groups = "drop"
    )
  plot_title <- "FST by Population"
  y_label <- "FST"
  
} else {
  stop(sprintf("Unknown stat_type: %s. Must be 'pi', 'dxy', or 'fst'", stat_type))
}

message(sprintf("Prepared data: %d populations\n", nrow(plot_df)))
message(sprintf("Value range: [%.6f, %.6f]\n", 
                min(plot_df$Value, na.rm = TRUE),
                max(plot_df$Value, na.rm = TRUE)))

# Order populations by value
plot_df <- plot_df %>%
  arrange(Value) %>%
  mutate(Population = factor(Population, levels = Population))

# Create boxplot with error bars (using CI)
message("\n=== CREATING BOXPLOT ===\n")

p <- ggplot(plot_df, aes(x = Population, y = Value)) +
  geom_point(size = 3, color = "steelblue") +
  geom_errorbar(
    aes(ymin = CI_low, ymax = CI_high),
    width = 0.2,
    color = "steelblue",
    alpha = 0.7
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.5) +
  labs(
    title = plot_title,
    x = "Population",
    y = y_label
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

# Save plots
message("\n=== SAVING OUTPUT ===\n")
message(sprintf("Output PDF: %s\n", output_pdf))
message(sprintf("Output RDS: %s\n", output_rds))

dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)
ggsave(
  filename = output_pdf,
  plot = p,
  width = max(8, nrow(plot_df) * 0.5),
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




