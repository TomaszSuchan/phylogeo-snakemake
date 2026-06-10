#!/usr/bin/env Rscript
# ROH group comparison plots for one configured group_by column

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

roh_file <- snakemake@input[["roh"]]
indpopdata_file <- snakemake@input[["indpopdata"]]
per_ind_file <- snakemake@input[["per_ind"]]
output_froh <- snakemake@output[["froh"]]
output_froh_class <- snakemake@output[["froh_class_panel"]]
output_nseg <- snakemake@output[["nseg"]]
output_class <- snakemake@output[["roh_class"]]
group_col <- snakemake@params[["group_col"]]
log_file <- snakemake@log[[1]]

log_con <- file(log_file, open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")
on.exit({
  while (sink.number(type = "message") > 0) sink(type = "message")
  while (sink.number(type = "output") > 0) sink(type = "output")
  close(log_con)
}, add = TRUE)

cat("=== ROH group plots:", group_col, "===\n")

per_ind <- fread(per_ind_file, header = TRUE)
if (!(group_col %in% colnames(per_ind))) {
  stop(
    "ROH group_by column '", group_col,
    "' not found in per-individual summary. Available columns: ",
    paste(colnames(per_ind), collapse = ", ")
  )
}
if (sum(!is.na(per_ind[[group_col]])) == 0) {
  stop("ROH group_by column '", group_col, "' has no non-NA values.")
}

n_groups <- length(unique(per_ind[[group_col]][!is.na(per_ind[[group_col]])]))
plot_width <- max(8, n_groups * 0.8)

dir.create(dirname(output_froh), recursive = TRUE, showWarnings = FALSE)

p_group_froh <- ggplot(per_ind, aes(x = .data[[group_col]], y = F_ROH, fill = .data[[group_col]])) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(
    title = paste("F_ROH by", group_col),
    x = group_col,
    y = "F_ROH"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )
ggsave(output_froh, p_group_froh, width = plot_width, height = 6)

froh_class_cols <- paste0("F_ROH_", ROH_CLASS_LEVELS)
froh_class_cols <- froh_class_cols[froh_class_cols %in% colnames(per_ind)]
if (length(froh_class_cols) == 0) {
  stop(
    "Per-individual ROH summary is missing per-class F_ROH columns (",
    paste(paste0("F_ROH_", ROH_CLASS_LEVELS), collapse = ", "),
    "). Re-run roh_summary after updating the workflow."
  )
}

froh_by_class_long <- per_ind %>%
  filter(!is.na(.data[[group_col]])) %>%
  select(Sample, all_of(group_col), all_of(froh_class_cols)) %>%
  pivot_longer(
    cols = all_of(froh_class_cols),
    names_to = "ROH_class",
    names_prefix = "F_ROH_",
    values_to = "F_ROH"
  ) %>%
  mutate(ROH_class = factor(ROH_class, levels = ROH_CLASS_LEVELS))

p_group_froh_class <- ggplot(
  froh_by_class_long,
  aes(x = .data[[group_col]], y = F_ROH, fill = .data[[group_col]])
) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~ ROH_class, nrow = 1, scales = "free_x") +
  labs(x = group_col, y = "F_ROH", caption = ROH_LENGTH_CLASS_CAPTION) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    strip.text = element_text(size = 10, face = "bold"),
    plot.caption = element_text(size = 8, hjust = 0)
  )
ggsave(output_froh_class, p_group_froh_class, width = max(12, n_groups * 2.4), height = 5.5)

p_group_nseg <- ggplot(per_ind, aes(x = .data[[group_col]], y = N_ROH_segments, fill = .data[[group_col]])) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(
    title = paste("Number of ROH Segments by", group_col),
    x = group_col,
    y = "Number of ROH Segments"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )
ggsave(output_nseg, p_group_nseg, width = plot_width, height = 6)

roh_data <- fread(roh_file, skip = 4, header = FALSE)
colnames(roh_data) <- c("Type", "Sample", "Chromosome", "Start", "End", "Length_bp", "N_markers", "Quality")
roh_data <- roh_data[Type == "RG", ]
roh_data$Length_Mb <- as.numeric(roh_data$Length_bp) / 1e6
roh_data$ROH_class <- factor(assign_roh_class(roh_data$Length_Mb), levels = ROH_CLASS_LEVELS)

indpopdata <- fread(indpopdata_file, header = TRUE)
ind_col <- if ("Ind" %in% colnames(indpopdata)) "Ind" else "Sample"
roh_data <- merge(roh_data, indpopdata, by.x = "Sample", by.y = ind_col, all.x = TRUE)

if (!(group_col %in% colnames(roh_data))) {
  stop("ROH group_by column '", group_col, "' not found after merging ROH segments with indpopdata.")
}

roh_class_summary <- roh_data %>%
  filter(!is.na(.data[[group_col]])) %>%
  group_by(.data[[group_col]], ROH_class) %>%
  summarise(N_segments = n(), .groups = "drop")

p_group_class <- ggplot(roh_class_summary, aes(x = .data[[group_col]], y = N_segments, fill = ROH_class)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  scale_fill_manual(values = ROH_CLASS_FILL, name = "ROH Class") +
  labs(
    title = paste("ROH Segments by Length Class and", group_col),
    x = group_col,
    y = "Number of Segments",
    caption = ROH_LENGTH_CLASS_CAPTION
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.caption = element_text(size = 8, hjust = 0),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
ggsave(
  output_class,
  p_group_class,
  width = max(10, length(unique(roh_class_summary[[group_col]])) * 0.8),
  height = 6
)

cat("=== Completed ROH group plots:", group_col, "===\n")
