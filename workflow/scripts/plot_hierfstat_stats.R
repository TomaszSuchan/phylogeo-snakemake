#!/usr/bin/env Rscript
# Combined Ho / Hs / Fis / population-specific Fst barplots (one panel each),
# faceted like the pixy diversity panels (ncol = 1, free_y).

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

STAT_LEVELS <- c("Ho", "Hs", "Fis", "Fst")
STAT_LABELS <- c(
  Ho = "Observed heterozygosity (Ho)",
  Hs = "Expected heterozygosity (Hs)",
  Fis = "Inbreeding coefficient (Fis)",
  Fst = "Population-specific Fst"
)

pdf(NULL)

log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

stats_file <- snakemake@input[["stats"]]
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
# Match the pixy combined-panel relative height: ~4/6 of single-panel height per panel.
panel_height <- plot_height * n_panels * (4 / 6)

message("\n=== READING HIERFSTAT SUMMARY ===\n")
stats_df <- read.table(stats_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
missing <- setdiff(c("population", STAT_LEVELS), colnames(stats_df))
if (length(missing) > 0) {
  stop(sprintf("%s missing columns: %s", stats_file, paste(missing, collapse = ", ")))
}

long_df <- stats_df %>%
  select(population, all_of(STAT_LEVELS)) %>%
  pivot_longer(all_of(STAT_LEVELS), names_to = "statistic", values_to = "value")
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
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey40") +
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
