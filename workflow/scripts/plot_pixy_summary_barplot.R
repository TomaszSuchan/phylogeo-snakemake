#!/usr/bin/env Rscript
# Generic pixy summary barplot with bootstrap CIs (pi / watterson_theta / tajima_d).

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

ylab_expression <- function(stat) {
  switch(
    stat,
    pi = expression(paste("Nucleotide Diversity (", pi, ")")),
    watterson_theta = expression(paste("Watterson's ", theta)),
    tajima_d = expression("Tajima's " * italic(D)),
    stop(sprintf("Unsupported pixy barplot statistic: %s", stat))
  )
}

value_column <- function(stat) {
  switch(
    stat,
    pi = "mean_pi",
    watterson_theta = "mean_watterson_theta",
    tajima_d = "mean_tajima_d",
    stop(sprintf("Unsupported pixy barplot statistic: %s", stat))
  )
}

pdf(NULL)

log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

summary_file <- snakemake@input[["summary"]]
popdata_file <- snakemake@input[["popdata"]]
output_pdf <- snakemake@output[["pdf"]]
output_rds <- snakemake@output[["rds"]]

stat <- as.character(snakemake@params[["stat"]])
grouping_name <- as.character(snakemake@params[["grouping"]])
group_colors <- group_fill_values(snakemake@params[["group_colors"]])
population_sort_by <- snakemake@params[["population_sort_by"]]
plot_width <- as.numeric(snakemake@params[["width"]])
plot_height <- as.numeric(snakemake@params[["height"]])
axis_title_size <- as.numeric(snakemake@params[["axis_title_size"]])
axis_text_size <- as.numeric(snakemake@params[["axis_text_size"]])
if (is.na(plot_width)) plot_width <- 8
if (is.na(plot_height)) plot_height <- 6
if (is.na(axis_title_size)) axis_title_size <- 10
if (is.na(axis_text_size)) axis_text_size <- 8

value_col <- value_column(stat)
y_lab <- ylab_expression(stat)

message(sprintf("\n=== READING %s SUMMARY ===\n", toupper(stat)))
sum_df <- read.table(summary_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
message(sprintf("Loaded %d populations for grouping '%s'\n", nrow(sum_df), grouping_name))

if (!(value_col %in% colnames(sum_df))) {
  stop(sprintf("Summary file missing value column '%s'", value_col))
}
for (req in c("ci_low", "ci_high", "population")) {
  if (!(req %in% colnames(sum_df))) {
    stop(sprintf("Summary file missing column '%s'", req))
  }
}

popdata <- read.table(popdata_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
site_col <- if ("Site" %in% colnames(popdata)) "Site" else colnames(popdata)[1]

if (is.null(group_sort_by(population_sort_by))) {
  message(sprintf("Using alphabetical order for %s\n", grouping_name))
} else if (length(group_sort_by(population_sort_by)) == 1 &&
           group_sort_by(population_sort_by)[1] %in% colnames(popdata)) {
  message(sprintf(
    "Sorting %s by column '%s'\n",
    grouping_name,
    group_sort_by(population_sort_by)[1]
  ))
} else {
  message(sprintf("Using configured level order for %s\n", grouping_name))
}

pop_order <- population_levels(sum_df$population, popdata, population_sort_by, site_col)
sum_df$population <- factor(sum_df$population, levels = pop_order)
sum_df <- sum_df[order(match(as.character(sum_df$population), pop_order)), , drop = FALSE]

if (!is.null(group_colors) && length(group_colors) > 0) {
  message(sprintf(
    "Using configured fill colors for %s (%d colors)\n",
    grouping_name,
    length(group_colors)
  ))
}

message("\n=== CREATING BARPLOT ===\n")

if (!is.null(group_colors) && length(group_colors) > 0) {
  p <- ggplot(sum_df, aes(x = population, y = .data[[value_col]], fill = population)) +
    geom_bar(stat = "identity", alpha = 0.7, color = "black", linewidth = 0.3) +
    scale_fill_manual(values = group_colors)
} else {
  p <- ggplot(sum_df, aes(x = population, y = .data[[value_col]])) +
    geom_bar(stat = "identity", fill = DEFAULT_FILL, alpha = 0.7, color = "black", linewidth = 0.3)
}

p <- p +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.2, color = "black", linewidth = 0.5) +
  labs(x = grouping_name, y = y_lab) +
  theme_bw() +
  theme(
    axis.title = element_text(size = axis_title_size),
    axis.text = element_text(size = axis_text_size),
    axis.text.x = element_text(angle = 45, hjust = 1, size = axis_text_size),
    panel.grid.major.x = element_blank(),
    legend.position = "none"
  )

message("\n=== SAVING OUTPUT ===\n")
message(sprintf("Output PDF: %s\n", output_pdf))
message(sprintf("Output RDS: %s\n", output_rds))

out_width <- min(50, max(plot_width, nrow(sum_df) * 0.3))
dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)
ggsave_pdf(
  filename = output_pdf,
  plot = p,
  width = out_width,
  height = plot_height,
  dpi = 300
)

dir.create(dirname(output_rds), recursive = TRUE, showWarnings = FALSE)
saveRDS(p, output_rds)

message("\n=== COMPLETED SUCCESSFULLY ===\n")

sink(type = "message")
sink(type = "output")
close(log_file)
