#!/usr/bin/env Rscript
# Pi barplot with confidence intervals, one plot per pixy group_by column.

library(ggplot2)
library(dplyr)

ggsave_utils <- tryCatch(
  file.path(dirname(normalizePath(snakemake@script)), "plot_ggsave_utils.R"),
  error = function(e) "workflow/scripts/plot_ggsave_utils.R"
)
if (file.exists(ggsave_utils)) {
  source(ggsave_utils)
} else {
  source("workflow/scripts/plot_ggsave_utils.R")
}

plot_group_utils <- tryCatch(
  file.path(dirname(normalizePath(snakemake@script)), "plot_group_utils.R"),
  error = function(e) "workflow/scripts/plot_group_utils.R"
)
if (file.exists(plot_group_utils)) {
  source(plot_group_utils)
} else {
  source("workflow/scripts/plot_group_utils.R")
}

DEFAULT_FILL <- "steelblue"

pdf(NULL)

log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

pi_summary_file <- snakemake@input[["pi_summary"]]
popdata_file <- snakemake@input[["popdata"]]
output_pdf <- snakemake@output[["pdf"]]
output_rds <- snakemake@output[["rds"]]

grouping_name <- as.character(snakemake@params[["grouping"]])
group_colors <- group_fill_values(snakemake@params[["group_colors"]])
population_sort_by <- snakemake@params[["population_sort_by"]]

message("\n=== READING PI SUMMARY ===\n")
pi_df <- read.table(pi_summary_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
message(sprintf("Loaded %d populations for grouping '%s'\n", nrow(pi_df), grouping_name))

popdata <- read.table(popdata_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
site_col <- if ("Site" %in% colnames(popdata)) "Site" else colnames(popdata)[1]

if (is.null(group_sort_by(population_sort_by))) {
  message(sprintf("Using alphabetical order for %s\n", grouping_name))
} else if (length(group_sort_by(population_sort_by)) == 1 &&
           group_sort_by(population_sort_by) %in% colnames(popdata)) {
  message(sprintf(
    "Sorting %s by column '%s'\n",
    grouping_name,
    group_sort_by(population_sort_by)[1]
  ))
} else {
  message(sprintf("Using configured level order for %s\n", grouping_name))
}

pop_order <- population_levels(pi_df$population, popdata, population_sort_by, site_col)
pi_df$population <- factor(pi_df$population, levels = pop_order)
pi_df <- pi_df[order(match(as.character(pi_df$population), pop_order)), , drop = FALSE]

if (!is.null(group_colors) && length(group_colors) > 0) {
  message(sprintf(
    "Using configured fill colors for %s (%d colors)\n",
    grouping_name,
    length(group_colors)
  ))
}

message("\n=== CREATING BARPLOT ===\n")

if (!is.null(group_colors) && length(group_colors) > 0) {
  p <- ggplot(pi_df, aes(x = population, y = mean_pi, fill = population)) +
    geom_bar(stat = "identity", alpha = 0.7, color = "black", linewidth = 0.3) +
    scale_fill_manual(values = group_colors)
} else {
  p <- ggplot(pi_df, aes(x = population, y = mean_pi)) +
    geom_bar(stat = "identity", fill = DEFAULT_FILL, alpha = 0.7, color = "black", linewidth = 0.3)
}

p <- p +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.2, color = "black", linewidth = 0.5) +
  labs(
    x = grouping_name,
    y = expression(paste("Nucleotide Diversity (", pi, ")"))
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    panel.grid.major.x = element_blank(),
    legend.position = "none"
  )

message("\n=== SAVING OUTPUT ===\n")
message(sprintf("Output PDF: %s\n", output_pdf))
message(sprintf("Output RDS: %s\n", output_rds))

plot_width <- min(50, max(8, nrow(pi_df) * 0.3))
dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)
ggsave_pdf(
  filename = output_pdf,
  plot = p,
  width = plot_width,
  height = 6,
  dpi = 300
)

dir.create(dirname(output_rds), recursive = TRUE, showWarnings = FALSE)
saveRDS(p, output_rds)

message("\n=== COMPLETED SUCCESSFULLY ===\n")

sink(type = "message")
sink(type = "output")
close(log_file)
