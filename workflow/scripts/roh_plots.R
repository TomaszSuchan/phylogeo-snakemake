#!/usr/bin/env Rscript
# ROH summary and group comparison plots

library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)

roh_utils <- tryCatch({
  file.path(dirname(normalizePath(snakemake@script)), "roh_utils.R")
}, error = function(e) "workflow/scripts/roh_utils.R")
if (file.exists(roh_utils)) {
  source(roh_utils)
} else {
  source("workflow/scripts/roh_utils.R")
}

histogram_utils <- tryCatch({
  file.path(dirname(normalizePath(snakemake@script)), "plot_histogram_utils.R")
}, error = function(e) "workflow/scripts/plot_histogram_utils.R")
if (file.exists(histogram_utils)) {
  source(histogram_utils)
} else {
  source("workflow/scripts/plot_histogram_utils.R")
}

group_col <- snakemake@params[["group_col"]]
plot_summary <- is.null(group_col)
per_ind_file <- snakemake@input[["per_ind"]]
log_file <- snakemake@log[[1]]

log_con <- file(log_file, open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")
on.exit({
  while (sink.number(type = "message") > 0) sink(type = "message")
  while (sink.number(type = "output") > 0) sink(type = "output")
  close(log_con)
}, add = TRUE)

out_pdf <- function(stem) snakemake@output[[stem]]
out_rds <- function(stem) snakemake@output[[paste0(stem, "_rds")]]

cat("=== ROH plots", if (plot_summary) "(summary)" else paste0("(group:", group_col, ")"), "===\n")

per_ind <- fread(per_ind_file, header = TRUE)

froh_class_cols <- paste0("F_ROH_", ROH_CLASS_LEVELS)
nroh_class_cols <- paste0("N_segments_", ROH_CLASS_LEVELS)
n_class_panels <- length(ROH_CLASS_LEVELS)
classes_panel_height <- 4 * n_class_panels

froh_by_class_long <- per_ind %>%
  dplyr::select(Sample, all_of(froh_class_cols)) %>%
  pivot_longer(
    cols = all_of(froh_class_cols),
    names_to = "ROH_class",
    names_prefix = "F_ROH_",
    values_to = "F_ROH"
  ) %>%
  mutate(ROH_class = factor(ROH_class, levels = ROH_CLASS_LEVELS))

nroh_by_class_long <- per_ind %>%
  dplyr::select(Sample, all_of(nroh_class_cols)) %>%
  pivot_longer(
    cols = all_of(nroh_class_cols),
    names_to = "ROH_class",
    names_prefix = "N_segments_",
    values_to = "N_ROH_segments"
  ) %>%
  mutate(ROH_class = factor(ROH_class, levels = ROH_CLASS_LEVELS))

group_boxplot <- function(data, col, y, fill_values = NULL) {
  if (!is.null(fill_values) && length(fill_values) > 0) {
    p <- ggplot(data, aes(x = .data[[col]], y = .data[[y]], fill = .data[[col]])) +
      geom_boxplot_styled() +
      scale_fill_manual(values = fill_values)
  } else {
    p <- ggplot(data, aes(x = .data[[col]], y = .data[[y]])) +
      geom_boxplot_styled(fill = HISTOGRAM_FILL)
  }
  p +
    geom_jitter(width = 0.2, alpha = 0.5) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
}

apply_group_levels <- function(data, col, sort_by = NULL) {
  level_order <- roh_group_levels(data, col, sort_by)
  data[[col]] <- factor(data[[col]], levels = level_order)
  data
}

plot_group_boxplots <- function(col, fill_values = NULL, sort_by = NULL) {
  if (!(col %in% colnames(per_ind))) {
    stop("ROH group_by column '", col, "' not found in per-individual summary.")
  }
  if (sum(!is.na(per_ind[[col]])) == 0) {
    stop("ROH group_by column '", col, "' has no non-NA values.")
  }

  level_order <- roh_group_levels(per_ind, col, sort_by)
  if (is.null(roh_group_sort_by(sort_by))) {
    cat("Using alphabetical order for", col, "\n")
  } else if (length(roh_group_sort_by(sort_by)) == 1 &&
             roh_group_sort_by(sort_by) %in% colnames(per_ind)) {
    cat("Sorting", col, "by column", roh_group_sort_by(sort_by)[1], "\n")
  } else {
    cat("Using configured level order for", col, "\n")
  }

  plot_data <- apply_group_levels(per_ind, col, sort_by)
  n_groups <- length(level_order)
  # ggplot2 rejects width/height >= 50 inches; cap below that limit.
  plot_width <- min(49, max(8, n_groups * 0.8))

  if (!is.null(fill_values) && length(fill_values) > 0) {
    cat("Using configured fill colors for", col, "\n")
  }

  p_froh <- group_boxplot(plot_data, col, "F_ROH", fill_values) +
    labs(x = col, y = "F_ROH")
  ggsave_pdf(out_pdf("froh"), p_froh, width = plot_width, height = 6)
  saveRDS(p_froh, out_rds("froh"))

  froh_group_long <- plot_data %>%
    dplyr::select(Sample, all_of(col), all_of(froh_class_cols)) %>%
    filter(!is.na(.data[[col]])) %>%
    pivot_longer(
      cols = all_of(froh_class_cols),
      names_to = "ROH_class",
      names_prefix = "F_ROH_",
      values_to = "F_ROH"
    ) %>%
    mutate(ROH_class = factor(ROH_class, levels = ROH_CLASS_LEVELS))

  p_froh_classes <- group_boxplot(froh_group_long, col, "F_ROH", fill_values) +
    facet_wrap(~ ROH_class, ncol = 1, scales = "free_y") +
    labs(x = col, y = "F_ROH")
  ggsave_pdf(out_pdf("froh_classes"), p_froh_classes, width = plot_width, height = classes_panel_height)
  saveRDS(p_froh_classes, out_rds("froh_classes"))

  p_nroh <- group_boxplot(plot_data, col, "N_ROH_segments", fill_values) +
    labs(x = col, y = "Number of ROH Segments")
  ggsave_pdf(out_pdf("nroh"), p_nroh, width = plot_width, height = 6)
  saveRDS(p_nroh, out_rds("nroh"))

  nroh_group_long <- plot_data %>%
    dplyr::select(Sample, all_of(col), all_of(nroh_class_cols)) %>%
    filter(!is.na(.data[[col]])) %>%
    pivot_longer(
      cols = all_of(nroh_class_cols),
      names_to = "ROH_class",
      names_prefix = "N_segments_",
      values_to = "N_ROH_segments"
    ) %>%
    mutate(ROH_class = factor(ROH_class, levels = ROH_CLASS_LEVELS))

  p_nroh_classes <- group_boxplot(nroh_group_long, col, "N_ROH_segments", fill_values) +
    facet_wrap(~ ROH_class, ncol = 1, scales = "free_y") +
    labs(x = col, y = "Number of ROH Segments")
  ggsave_pdf(out_pdf("nroh_classes"), p_nroh_classes, width = plot_width, height = classes_panel_height)
  saveRDS(p_nroh_classes, out_rds("nroh_classes"))
}

if (!plot_summary) {
  dir.create(dirname(out_pdf("froh")), recursive = TRUE, showWarnings = FALSE)
  group_fill_values <- roh_group_fill_values(snakemake@params[["group_colors"]])
  plot_group_boxplots(
    group_col,
    fill_values = group_fill_values,
    sort_by = snakemake@params[["group_sort_by"]]
  )
  cat("=== ROH plots complete ===\n")
  quit(save = "no", status = 0)
}

roh_file <- snakemake@input[["roh"]]
roh_data <- fread(roh_file, skip = 4, header = FALSE)
colnames(roh_data) <- c("Type", "Sample", "Chromosome", "Start", "End", "Length_bp", "N_markers", "Quality")
roh_data <- roh_data[Type == "RG", ]
roh_data$Length_Mb <- as.numeric(roh_data$Length_bp) / 1e6
roh_data$ROH_class <- factor(assign_roh_class(roh_data$Length_Mb), levels = ROH_CLASS_LEVELS)

dir.create(dirname(out_pdf("froh_histogram")), recursive = TRUE, showWarnings = FALSE)

cat("Creating summary plots...\n")

p1 <- ggplot(per_ind, aes(x = F_ROH)) +
  geom_histogram_styled(bins = 30) +
  labs(x = "F_ROH", y = "Number of Individuals") +
  histogram_plot_theme()
ggsave_pdf(out_pdf("froh_histogram"), p1, width = 8, height = 6)
saveRDS(p1, out_rds("froh_histogram"))

p2 <- ggplot(per_ind, aes(x = N_ROH_segments)) +
  geom_histogram_styled(bins = 30) +
  labs(x = "Number of ROH Segments", y = "Number of Individuals") +
  histogram_plot_theme()
ggsave_pdf(out_pdf("nroh_histogram"), p2, width = 8, height = 6)
saveRDS(p2, out_rds("nroh_histogram"))

p3 <- ggplot(per_ind, aes(x = Total_ROH_length_Mb)) +
  geom_histogram_styled(bins = 30) +
  labs(x = "Total ROH Length (Mb)", y = "Number of Individuals") +
  histogram_plot_theme()
ggsave_pdf(out_pdf("length_histogram"), p3, width = 8, height = 6)
saveRDS(p3, out_rds("length_histogram"))

p4 <- ggplot(roh_data, aes(x = Length_Mb)) +
  geom_histogram_styled() +
  labs(x = "ROH Length (Mb)", y = "Number of Segments") +
  histogram_plot_theme()
ggsave_pdf(out_pdf("length_segments_histogram"), p4, width = 8, height = 6)
saveRDS(p4, out_rds("length_segments_histogram"))

froh_class_summary <- froh_by_class_long %>%
  group_by(ROH_class) %>%
  summarise(Mean_F_ROH = mean(F_ROH, na.rm = TRUE), .groups = "drop")

p5 <- ggplot(froh_class_summary, aes(x = ROH_class, y = Mean_F_ROH)) +
  geom_bar(stat = "identity") +
  scale_x_discrete(limits = ROH_CLASS_LEVELS) +
  labs(x = "ROH Length Class", y = "Mean F_ROH") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave_pdf(out_pdf("froh_classes"), p5, width = 10, height = 6.5)
saveRDS(p5, out_rds("froh_classes"))

p6 <- ggplot(froh_by_class_long, aes(x = F_ROH)) +
  geom_histogram_styled(bins = 30) +
  facet_wrap(~ ROH_class, ncol = 1, scales = "free_y") +
  labs(x = "F_ROH", y = "Number of Individuals") +
  histogram_plot_theme()
ggsave_pdf(out_pdf("froh_classes_histogram"), p6, width = 8, height = classes_panel_height)
saveRDS(p6, out_rds("froh_classes_histogram"))

nroh_class_counts <- roh_data %>%
  group_by(ROH_class) %>%
  summarise(N_segments = n(), .groups = "drop")

p7 <- ggplot(nroh_class_counts, aes(x = ROH_class, y = N_segments)) +
  geom_bar(stat = "identity") +
  scale_x_discrete(limits = ROH_CLASS_LEVELS) +
  labs(x = "ROH Length Class", y = "Number of Segments") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave_pdf(out_pdf("nroh_classes"), p7, width = 10, height = 6.5)
saveRDS(p7, out_rds("nroh_classes"))

p8 <- ggplot(nroh_by_class_long, aes(x = N_ROH_segments)) +
  geom_histogram_styled(bins = 30) +
  facet_wrap(~ ROH_class, ncol = 1, scales = "free_y") +
  labs(x = "Number of ROH Segments", y = "Number of Individuals") +
  histogram_plot_theme()
ggsave_pdf(out_pdf("nroh_classes_histogram"), p8, width = 8, height = classes_panel_height)
saveRDS(p8, out_rds("nroh_classes_histogram"))

p9 <- ggplot(per_ind, aes(x = N_ROH_segments, y = F_ROH)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(x = "Number of ROH Segments", y = "F_ROH") +
  theme_bw()
ggsave_pdf(out_pdf("froh_vs_nroh"), p9, width = 8, height = 6)
saveRDS(p9, out_rds("froh_vs_nroh"))

cat("=== ROH plots complete ===\n")
