#!/usr/bin/env Rscript
# Combined pi / Watterson's theta / Tajima's D barplots (one panel each),
# faceted like ROH length-class plots (ncol = 1, free_y).

library(ggplot2)
library(dplyr)
library(tidyr)

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

STAT_LEVELS <- c("pi", "watterson_theta", "tajima_d")
STAT_LABELS <- c(
  pi = "Nucleotide diversity (π)",
  watterson_theta = "Watterson's θ",
  tajima_d = "Tajima's D"
)

pdf(NULL)

log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

pi_file <- snakemake@input[["pi_summary"]]
theta_file <- snakemake@input[["watterson_theta_summary"]]
tajima_file <- snakemake@input[["tajima_d_summary"]]
popdata_file <- snakemake@input[["popdata"]]
output_pdf <- snakemake@output[["pdf"]]
output_rds <- snakemake@output[["rds"]]

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

n_panels <- length(STAT_LEVELS)
# Match ROH classed-plot relative height: ~4/6 of single-panel height per panel.
panel_height <- plot_height * n_panels * (4 / 6)

read_stat_summary <- function(path, value_col, stat_id) {
  df <- read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  needed <- c("population", value_col, "ci_low", "ci_high")
  missing <- setdiff(needed, colnames(df))
  if (length(missing) > 0) {
    stop(sprintf("%s missing columns: %s", path, paste(missing, collapse = ", ")))
  }
  tibble(
    population = as.character(df$population),
    value = as.numeric(df[[value_col]]),
    ci_low = as.numeric(df$ci_low),
    ci_high = as.numeric(df$ci_high),
    statistic = stat_id
  )
}

message("\n=== READING DIVERSITY SUMMARIES ===\n")
long_df <- bind_rows(
  read_stat_summary(pi_file, "mean_pi", "pi"),
  read_stat_summary(theta_file, "mean_watterson_theta", "watterson_theta"),
  read_stat_summary(tajima_file, "mean_tajima_d", "tajima_d")
)
message(sprintf("Combined %d rows across %d statistics\n", nrow(long_df), n_panels))

popdata <- read.table(popdata_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
site_col <- if ("Site" %in% colnames(popdata)) "Site" else colnames(popdata)[1]

pop_order <- population_levels(unique(long_df$population), popdata, population_sort_by, site_col)
long_df$population <- factor(long_df$population, levels = pop_order)
long_df$statistic <- factor(long_df$statistic, levels = STAT_LEVELS, labels = STAT_LABELS[STAT_LEVELS])
long_df <- long_df[order(match(as.character(long_df$population), pop_order)), , drop = FALSE]

if (!is.null(group_colors) && length(group_colors) > 0) {
  message(sprintf(
    "Using configured fill colors for %s (%d colors)\n",
    grouping_name,
    length(group_colors)
  ))
}

message("\n=== CREATING COMBINED PANEL PLOT ===\n")

if (!is.null(group_colors) && length(group_colors) > 0) {
  p <- ggplot(long_df, aes(x = population, y = value, fill = population)) +
    geom_bar(stat = "identity", alpha = 0.7, color = "black", linewidth = 0.3) +
    scale_fill_manual(values = group_colors)
} else {
  p <- ggplot(long_df, aes(x = population, y = value)) +
    geom_bar(stat = "identity", fill = DEFAULT_FILL, alpha = 0.7, color = "black", linewidth = 0.3)
}

p <- p +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.2, color = "black", linewidth = 0.5) +
  facet_wrap(~ statistic, ncol = 1, scales = "free_y") +
  labs(x = grouping_name, y = NULL) +
  theme_bw() +
  theme(
    axis.title = element_text(size = axis_title_size),
    axis.text = element_text(size = axis_text_size),
    axis.text.x = element_text(angle = 45, hjust = 1, size = axis_text_size),
    strip.text = element_text(size = axis_title_size),
    panel.grid.major.x = element_blank(),
    legend.position = "none"
  )

message("\n=== SAVING OUTPUT ===\n")
message(sprintf("Output PDF: %s\n", output_pdf))
message(sprintf("Output RDS: %s\n", output_rds))

out_width <- min(49, max(plot_width, length(pop_order) * 0.3))
dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)
ggsave_pdf(
  filename = output_pdf,
  plot = p,
  width = out_width,
  height = panel_height,
  dpi = 300
)

dir.create(dirname(output_rds), recursive = TRUE, showWarnings = FALSE)
saveRDS(p, output_rds)

message("\n=== COMPLETED SUCCESSFULLY ===\n")

sink(type = "message")
sink(type = "output")
close(log_file)
