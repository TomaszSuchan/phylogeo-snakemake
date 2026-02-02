#!/usr/bin/env Rscript
# ROH Summary and Visualization Script
# Processes bcftools roh output and generates summary statistics and visualizations

library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)

# Snakemake inputs/outputs
roh_file <- snakemake@input[["roh"]]
indpopdata_file <- snakemake@input[["indpopdata"]]
output_summary <- snakemake@output[["summary"]]
output_per_ind <- snakemake@output[["per_ind"]]
output_plots <- snakemake@output[["plots"]]
log_file <- snakemake@log[[1]]
group_by <- snakemake@params[["group_by"]]

# Redirect output to log
sink(log_file, type = "output", split = TRUE)

cat("=== ROH Summary Analysis ===\n")
cat("ROH file:", roh_file, "\n")
cat("Indpopdata file:", indpopdata_file, "\n")
cat("Group by columns:", paste(group_by, collapse = ", "), "\n")
cat("\n")

# Read ROH data
cat("Reading ROH data...\n")
roh_data <- fread(roh_file, skip = 4, header = FALSE)
colnames(roh_data) <- c("Type", "Sample", "Chromosome", "Start", "End", "Length_bp", "N_markers", "Quality")

# Filter to only RG (Run of Homozygosity) entries
roh_data <- roh_data[Type == "RG", ]

cat("  Total ROH segments:", nrow(roh_data), "\n")
cat("  Unique samples:", length(unique(roh_data$Sample)), "\n")
cat("  Unique chromosomes:", length(unique(roh_data$Chromosome)), "\n")
cat("\n")

# Convert Length_bp to numeric (in case it's character)
roh_data$Length_bp <- as.numeric(roh_data$Length_bp)
roh_data$Length_Mb <- roh_data$Length_bp / 1e6

# Read indpopdata if available
indpopdata <- NULL
if (file.exists(indpopdata_file)) {
  cat("Reading indpopdata...\n")
  indpopdata <- fread(indpopdata_file, header = TRUE)
  cat("  Columns in indpopdata:", paste(colnames(indpopdata), collapse = ", "), "\n")
  cat("  Samples in indpopdata:", nrow(indpopdata), "\n")
  
  # Merge with ROH data (use Ind column from indpopdata, Sample from roh_data)
  # First check if column is called "Ind" or "Sample"
  ind_col <- if ("Ind" %in% colnames(indpopdata)) "Ind" else "Sample"
  roh_data <- merge(roh_data, indpopdata, by.x = "Sample", by.y = ind_col, all.x = TRUE)
} else {
  cat("No indpopdata file found. Proceeding without population information.\n")
}

cat("\n")

# ==============================================================================
# 1. Per-individual ROH statistics
# ==============================================================================

cat("Calculating per-individual statistics...\n")

per_ind <- roh_data %>%
  group_by(Sample) %>%
  summarise(
    N_ROH_segments = n(),
    Total_ROH_length_bp = sum(Length_bp, na.rm = TRUE),
    Total_ROH_length_Mb = sum(Length_Mb, na.rm = TRUE),
    Mean_ROH_length_Mb = mean(Length_Mb, na.rm = TRUE),
    Median_ROH_length_Mb = median(Length_Mb, na.rm = TRUE),
    Max_ROH_length_Mb = max(Length_Mb, na.rm = TRUE),
    Min_ROH_length_Mb = min(Length_Mb, na.rm = TRUE),
    .groups = "drop"
  )

# Add population information if available
if (!is.null(indpopdata)) {
  ind_col <- if ("Ind" %in% colnames(indpopdata)) "Ind" else "Sample"
  per_ind <- merge(per_ind, indpopdata, by.x = "Sample", by.y = ind_col, all.x = TRUE)
}

# Calculate F_ROH (fraction of genome in ROH)
# This requires genome size - we'll estimate from the ROH data
# Sum of all unique chromosome lengths covered by ROH
genome_size <- roh_data %>%
  group_by(Chromosome) %>%
  summarise(
    chr_length = max(End, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  summarise(total_genome_size = sum(chr_length, na.rm = TRUE)) %>%
  pull(total_genome_size)

cat("  Estimated genome size (from ROH data):", genome_size / 1e6, "Mb\n")

per_ind$F_ROH <- per_ind$Total_ROH_length_bp / genome_size

# Categorize ROH by length classes
# Short: < 1 Mb (ancient inbreeding)
# Medium: 1-5 Mb (intermediate)
# Long: > 5 Mb (recent inbreeding)
roh_data$ROH_class <- case_when(
  roh_data$Length_Mb < 1 ~ "Short (<1 Mb)",
  roh_data$Length_Mb >= 1 & roh_data$Length_Mb < 5 ~ "Medium (1-5 Mb)",
  roh_data$Length_Mb >= 5 ~ "Long (>5 Mb)"
)

# Count ROH segments by class per individual
roh_by_class <- roh_data %>%
  group_by(Sample, ROH_class) %>%
  summarise(
    N_segments = n(),
    Total_length_Mb = sum(Length_Mb, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    id_cols = Sample,
    names_from = ROH_class,
    values_from = c(N_segments, Total_length_Mb),
    values_fill = 0
  )

# Merge with per_ind
per_ind <- merge(per_ind, roh_by_class, by = "Sample", all.x = TRUE)

# Fill NA with 0 for missing classes
class_cols <- grep("N_segments_|Total_length_Mb_", colnames(per_ind), value = TRUE)
if (length(class_cols) > 0) {
  for (col in class_cols) {
    per_ind[[col]][is.na(per_ind[[col]])] <- 0
  }
}

cat("  Processed", nrow(per_ind), "individuals\n")
cat("\n")

# ==============================================================================
# 2. Population-level summary (if populations available)
# ==============================================================================

# Use first group_by column for summary if available
pop_summary <- NULL
if (!is.null(indpopdata) && !is.null(group_by) && length(group_by) > 0 && !("none" %in% group_by)) {
  group_by <- group_by[group_by != "none"]
  if (length(group_by) > 0) {
    available_cols <- colnames(per_ind)
    valid_group_by <- group_by[group_by %in% available_cols]
    
    if (length(valid_group_by) > 0) {
      # Use first valid grouping column for summary
      summary_col <- valid_group_by[1]
      
      if (sum(!is.na(per_ind[[summary_col]])) > 0) {
        cat("Calculating population-level statistics by", summary_col, "...\n")
        
        pop_summary <- per_ind %>%
          group_by(.data[[summary_col]]) %>%
          summarise(
            N_individuals = n(),
            Mean_F_ROH = mean(F_ROH, na.rm = TRUE),
            Median_F_ROH = median(F_ROH, na.rm = TRUE),
            Mean_N_ROH = mean(N_ROH_segments, na.rm = TRUE),
            Mean_total_ROH_Mb = mean(Total_ROH_length_Mb, na.rm = TRUE),
            Mean_max_ROH_Mb = mean(Max_ROH_length_Mb, na.rm = TRUE),
            .groups = "drop"
          )
        
        # Rename the grouping column to "Population" for consistency
        colnames(pop_summary)[1] <- "Population"
        
        cat("  Processed", nrow(pop_summary), "groups\n")
        cat("\n")
      }
    }
  }
}

# ==============================================================================
# 3. Overall summary statistics
# ==============================================================================

cat("Overall summary statistics:\n")
cat("  Mean F_ROH:", mean(per_ind$F_ROH, na.rm = TRUE), "\n")
cat("  Median F_ROH:", median(per_ind$F_ROH, na.rm = TRUE), "\n")
cat("  Mean number of ROH segments:", mean(per_ind$N_ROH_segments, na.rm = TRUE), "\n")
cat("  Mean total ROH length:", mean(per_ind$Total_ROH_length_Mb, na.rm = TRUE), "Mb\n")
cat("  Mean max ROH length:", mean(per_ind$Max_ROH_length_Mb, na.rm = TRUE), "Mb\n")
cat("\n")

# ==============================================================================
# 4. Create visualizations
# ==============================================================================

cat("Creating visualizations...\n")

# Create output directory for plots (output_plots is already a directory path)
plot_dir <- output_plots
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}

# Plot 1: F_ROH distribution
p1 <- ggplot(per_ind, aes(x = F_ROH)) +
  geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7, color = "black") +
  labs(
    title = "Distribution of F_ROH (Fraction of Genome in ROH)",
    x = "F_ROH",
    y = "Number of Individuals"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(file.path(plot_dir, "froh_histogram.pdf"), p1, width = 8, height = 6)

# Plot 2: Number of ROH segments distribution
p2 <- ggplot(per_ind, aes(x = N_ROH_segments)) +
  geom_histogram(bins = 30, fill = "darkgreen", alpha = 0.7, color = "black") +
  labs(
    title = "Distribution of Number of ROH Segments",
    x = "Number of ROH Segments",
    y = "Number of Individuals"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(file.path(plot_dir, "n_roh_segments_histogram.pdf"), p2, width = 8, height = 6)

# Plot 3: Total ROH length distribution
p3 <- ggplot(per_ind, aes(x = Total_ROH_length_Mb)) +
  geom_histogram(bins = 30, fill = "purple", alpha = 0.7, color = "black") +
  labs(
    title = "Distribution of Total ROH Length",
    x = "Total ROH Length (Mb)",
    y = "Number of Individuals"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(file.path(plot_dir, "total_roh_length_histogram.pdf"), p3, width = 8, height = 6)

# Plot 4: ROH length distribution (all segments)
p4 <- ggplot(roh_data, aes(x = Length_Mb)) +
  geom_histogram(bins = 50, fill = "orange", alpha = 0.7, color = "black") +
  scale_x_log10() +
  labs(
    title = "Distribution of ROH Segment Lengths",
    x = "ROH Length (Mb, log10 scale)",
    y = "Number of Segments"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(file.path(plot_dir, "roh_segment_length_distribution.pdf"), p4, width = 8, height = 6)

# Plot 5: ROH by length class
p5 <- ggplot(roh_data, aes(x = ROH_class, fill = ROH_class)) +
  geom_bar(color = "black") +
  scale_fill_manual(values = c("Short (<1 Mb)" = "lightblue", 
                                "Medium (1-5 Mb)" = "lightgreen", 
                                "Long (>5 Mb)" = "salmon")) +
  labs(
    title = "ROH Segments by Length Class",
    x = "ROH Length Class",
    y = "Number of Segments"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")

ggsave(file.path(plot_dir, "roh_by_class.pdf"), p5, width = 8, height = 6)

# Plot 6: F_ROH vs Number of ROH segments
p6 <- ggplot(per_ind, aes(x = N_ROH_segments, y = F_ROH)) +
  geom_point(alpha = 0.6, color = "steelblue") +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  labs(
    title = "F_ROH vs Number of ROH Segments",
    x = "Number of ROH Segments",
    y = "F_ROH"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(file.path(plot_dir, "froh_vs_n_segments.pdf"), p6, width = 8, height = 6)

# Population-level plots if available and group_by is specified
if (!is.null(indpopdata) && !is.null(group_by) && length(group_by) > 0 && !("none" %in% group_by)) {
  # Filter out "none" if present
  group_by <- group_by[group_by != "none"]
  
  if (length(group_by) > 0) {
    # Check which grouping columns exist in the data
    available_cols <- colnames(per_ind)
    valid_group_by <- group_by[group_by %in% available_cols]
    
    if (length(valid_group_by) > 0) {
      cat("Generating population-level plots for grouping columns:", paste(valid_group_by, collapse = ", "), "\n")
      
      # Generate plots for each grouping column
      for (group_col in valid_group_by) {
        # Check if column has valid data (not all NA)
        if (sum(!is.na(per_ind[[group_col]])) > 0 && length(unique(per_ind[[group_col]][!is.na(per_ind[[group_col]])])) > 1) {
          
          # Plot: F_ROH by grouping column
          p_group_froh <- ggplot(per_ind, aes(x = .data[[group_col]], y = F_ROH, fill = .data[[group_col]])) +
            geom_boxplot(alpha = 0.7) +
            geom_jitter(width = 0.2, alpha = 0.5) +
            labs(
              title = paste("F_ROH by", group_col),
              x = group_col,
              y = "F_ROH"
            ) +
            theme_minimal() +
            theme(plot.title = element_text(hjust = 0.5),
                  axis.text.x = element_text(angle = 45, hjust = 1),
                  legend.position = "none")
          
          ggsave(file.path(plot_dir, paste0("froh_by_", group_col, ".pdf")), 
                 p_group_froh, width = max(8, length(unique(per_ind[[group_col]][!is.na(per_ind[[group_col]])])) * 0.8), 
                 height = 6)
          
          # Plot: Number of ROH segments by grouping column
          p_group_nseg <- ggplot(per_ind, aes(x = .data[[group_col]], y = N_ROH_segments, fill = .data[[group_col]])) +
            geom_boxplot(alpha = 0.7) +
            geom_jitter(width = 0.2, alpha = 0.5) +
            labs(
              title = paste("Number of ROH Segments by", group_col),
              x = group_col,
              y = "Number of ROH Segments"
            ) +
            theme_minimal() +
            theme(plot.title = element_text(hjust = 0.5),
                  axis.text.x = element_text(angle = 45, hjust = 1),
                  legend.position = "none")
          
          ggsave(file.path(plot_dir, paste0("n_roh_segments_by_", group_col, ".pdf")), 
                 p_group_nseg, width = max(8, length(unique(per_ind[[group_col]][!is.na(per_ind[[group_col]])])) * 0.8), 
                 height = 6)
          
          # Plot: ROH length classes by grouping column
          roh_class_summary <- roh_data %>%
            filter(!is.na(.data[[group_col]])) %>%
            group_by(.data[[group_col]], ROH_class) %>%
            summarise(N_segments = n(), .groups = "drop")
          
          p_group_class <- ggplot(roh_class_summary, aes(x = .data[[group_col]], y = N_segments, fill = ROH_class)) +
            geom_bar(stat = "identity", position = "stack", color = "black") +
            scale_fill_manual(values = c("Short (<1 Mb)" = "lightblue", 
                                          "Medium (1-5 Mb)" = "lightgreen", 
                                          "Long (>5 Mb)" = "salmon"),
                              name = "ROH Class") +
            labs(
              title = paste("ROH Segments by Length Class and", group_col),
              x = group_col,
              y = "Number of Segments"
            ) +
            theme_minimal() +
            theme(plot.title = element_text(hjust = 0.5),
                  axis.text.x = element_text(angle = 45, hjust = 1))
          
          ggsave(file.path(plot_dir, paste0("roh_class_by_", group_col, ".pdf")), 
                 p_group_class, width = max(10, length(unique(roh_class_summary[[group_col]])) * 0.8), 
                 height = 6)
          
          cat("  Generated plots for", group_col, "\n")
        } else {
          cat("  Skipping", group_col, "- insufficient data or only one group\n")
        }
      }
    } else {
      cat("  Warning: None of the specified group_by columns found in indpopdata\n")
      cat("  Available columns:", paste(available_cols, collapse = ", "), "\n")
    }
  }
}

cat("  Created", length(list.files(plot_dir, pattern = "\\.pdf$")), "plot files\n")
cat("\n")

# ==============================================================================
# 5. Write output files
# ==============================================================================

cat("Writing output files...\n")

# Overall summary
summary_stats <- data.frame(
  Metric = c(
    "Total_ROH_segments",
    "Unique_samples",
    "Unique_chromosomes",
    "Mean_F_ROH",
    "Median_F_ROH",
    "Mean_N_ROH_segments",
    "Mean_total_ROH_length_Mb",
    "Mean_max_ROH_length_Mb",
    "Genome_size_Mb"
  ),
  Value = c(
    nrow(roh_data),
    length(unique(roh_data$Sample)),
    length(unique(roh_data$Chromosome)),
    mean(per_ind$F_ROH, na.rm = TRUE),
    median(per_ind$F_ROH, na.rm = TRUE),
    mean(per_ind$N_ROH_segments, na.rm = TRUE),
    mean(per_ind$Total_ROH_length_Mb, na.rm = TRUE),
    mean(per_ind$Max_ROH_length_Mb, na.rm = TRUE),
    genome_size / 1e6
  )
)

write.table(summary_stats, file = output_summary, sep = "\t", quote = FALSE, row.names = FALSE)

# Per-individual summary
write.table(per_ind, file = output_per_ind, sep = "\t", quote = FALSE, row.names = FALSE)

# Population summary if available
if (!is.null(pop_summary)) {
  pop_summary_file <- gsub("_per_ind\\.txt$", "_pop_summary.txt", output_per_ind)
  write.table(pop_summary, file = pop_summary_file, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("  Population summary written to:", pop_summary_file, "\n")
}

cat("  Summary statistics written to:", output_summary, "\n")
cat("  Per-individual statistics written to:", output_per_ind, "\n")
cat("\n")

cat("=== Analysis Complete ===\n")

sink()

