#!/usr/bin/env Rscript
# ROH summary statistics (tables and tests; plots are in roh_plots.R)

library(data.table)
library(dplyr)
library(tidyr)
library(multcomp)

roh_utils <- tryCatch({
  file.path(dirname(normalizePath(snakemake@script)), "roh_utils.R")
}, error = function(e) "workflow/scripts/roh_utils.R")
if (file.exists(roh_utils)) {
  source(roh_utils)
} else {
  source("workflow/scripts/roh_utils.R")
}

# Snakemake inputs/outputs
roh_file <- snakemake@input[["roh"]]
indpopdata_file <- snakemake@input[["indpopdata"]]
output_summary <- snakemake@output[["summary"]]
output_per_ind <- snakemake@output[["per_ind"]]
output_stats <- snakemake@output[["stats"]]
log_file <- snakemake@log[[1]]
group_by <- snakemake@params[["group_by"]]

# Redirect output to log
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")
on.exit({
  while (sink.number(type = "message") > 0) sink(type = "message")
  while (sink.number(type = "output") > 0) sink(type = "output")
  close(log_con)
}, add = TRUE)

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

roh_data$ROH_class <- factor(assign_roh_class(roh_data$Length_Mb), levels = ROH_CLASS_LEVELS)
cat("  ROH length classes use G = 100/(2 × cM) at ", ROH_RECOMBINATION_RATE_CM_PER_MB,
    " cM/Mb: ", paste(ROH_CLASS_LEVELS, collapse = "; "), "\n", sep = "")

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

for (cls in ROH_CLASS_LEVELS) {
  len_col <- paste0("Total_length_Mb_", cls)
  froh_col <- paste0("F_ROH_", cls)
  if (len_col %in% colnames(per_ind)) {
    per_ind[[froh_col]] <- per_ind[[len_col]] * 1e6 / genome_size
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
# 4. Calculate statistics and perform multiple comparison tests
# ==============================================================================

cat("Calculating statistics and performing multiple comparison tests...\n")

# Initialize output for statistics
stats_output <- character()

# Process each grouping column
if (!is.null(indpopdata) && !is.null(group_by) && length(group_by) > 0 && !("none" %in% group_by)) {
  group_by <- group_by[group_by != "none"]
  
  if (length(group_by) > 0) {
    available_cols <- colnames(per_ind)
    valid_group_by <- group_by[group_by %in% available_cols]
    
    if (length(valid_group_by) > 0) {
      for (group_col in valid_group_by) {
        # Check if column has valid data
        if (sum(!is.na(per_ind[[group_col]])) > 0 && length(unique(per_ind[[group_col]][!is.na(per_ind[[group_col]])])) > 1) {
          
          # Filter data for this grouping column
          data_subset <- per_ind[!is.na(per_ind[[group_col]]), ]
          groups <- unique(data_subset[[group_col]])
          
          if (length(groups) > 1) {
            stats_output <- c(stats_output, 
              paste0("\n", "=", paste(rep("=", 79), collapse = ""), "\n"),
              paste0("STATISTICS FOR GROUPING COLUMN: ", group_col, "\n"),
              paste0("=", paste(rep("=", 79), collapse = ""), "\n")
            )
            
            # ============================================================
            # F_ROH Statistics
            # ============================================================
            stats_output <- c(stats_output, "\n### F_ROH Statistics ###\n\n")
            
            # Calculate statistics by group
            froh_stats <- data_subset %>%
              group_by(.data[[group_col]]) %>%
              summarise(
                n = n(),
                mean = mean(F_ROH, na.rm = TRUE),
                median = median(F_ROH, na.rm = TRUE),
                sd = sd(F_ROH, na.rm = TRUE),
                min = min(F_ROH, na.rm = TRUE),
                max = max(F_ROH, na.rm = TRUE),
                .groups = "drop"
              )
            
            # Format and add to output
            stats_output <- c(stats_output, 
              paste(sprintf("%-20s", c("Group", "n", "mean", "median", "sd", "min", "max")), collapse = "\t"),
              "\n"
            )
            for (i in 1:nrow(froh_stats)) {
              stats_output <- c(stats_output,
                paste(sprintf("%-20s", c(
                  as.character(froh_stats[[group_col]][i]),
                  froh_stats$n[i],
                  sprintf("%.6f", froh_stats$mean[i]),
                  sprintf("%.6f", froh_stats$median[i]),
                  sprintf("%.6f", froh_stats$sd[i]),
                  sprintf("%.6f", froh_stats$min[i]),
                  sprintf("%.6f", froh_stats$max[i])
                )), collapse = "\t"),
                "\n"
              )
            }
            
            # Perform ANOVA
            stats_output <- c(stats_output, "\n### F_ROH ANOVA ###\n\n")
            tryCatch({
              # Create formula dynamically
              formula_str <- paste("F_ROH ~ as.factor(", group_col, ")", sep = "")
              aov_froh <- aov(as.formula(formula_str), data = data_subset)
              aov_summary <- summary(aov_froh)
              stats_output <- c(stats_output, capture.output(print(aov_summary)), "\n")
              
              # Multiple comparison test (Tukey HSD) - always perform regardless of significance
              p_value <- aov_summary[[1]][["Pr(>F)"]][1]
              if (!is.na(p_value)) {
                if (p_value < 0.05) {
                  stats_output <- c(stats_output, paste0("\nANOVA p-value: ", sprintf("%.6f", p_value), " (significant at α=0.05)\n"))
                } else {
                  stats_output <- c(stats_output, paste0("\nANOVA p-value: ", sprintf("%.6f", p_value), " (not significant at α=0.05)\n"))
                }
                stats_output <- c(stats_output, "\n### F_ROH Multiple Comparisons (Tukey HSD) ###\n\n")
                tukey_froh <- TukeyHSD(aov_froh)
                stats_output <- c(stats_output, capture.output(print(tukey_froh)), "\n")
                
                # Also get pairwise comparisons using multcomp
                tryCatch({
                  data_subset$group_factor <- as.factor(data_subset[[group_col]])
                  aov_froh2 <- aov(F_ROH ~ group_factor, data = data_subset)
                  glht_froh <- glht(aov_froh2, linfct = mcp(group_factor = "Tukey"))
                  summary_glht <- summary(glht_froh)
                  stats_output <- c(stats_output, "\n### F_ROH Pairwise Comparisons (multcomp) ###\n\n")
                  stats_output <- c(stats_output, capture.output(print(summary_glht)), "\n")
                }, error = function(e) {
                  stats_output <<- c(stats_output, "Note: multcomp pairwise comparisons skipped (using TukeyHSD results above)\n")
                })
              } else {
                stats_output <- c(stats_output, "\nWarning: Could not calculate ANOVA p-value.\n\n")
              }
            }, error = function(e) {
              stats_output <<- c(stats_output, "Error in ANOVA:", e$message, "\n\n")
            })
            
            # ============================================================
            # Number of ROH Segments Statistics
            # ============================================================
            stats_output <- c(stats_output, "\n### Number of ROH Segments Statistics ###\n\n")
            
            # Calculate statistics by group
            nseg_stats <- data_subset %>%
              group_by(.data[[group_col]]) %>%
              summarise(
                n = n(),
                mean = mean(N_ROH_segments, na.rm = TRUE),
                median = median(N_ROH_segments, na.rm = TRUE),
                sd = sd(N_ROH_segments, na.rm = TRUE),
                min = min(N_ROH_segments, na.rm = TRUE),
                max = max(N_ROH_segments, na.rm = TRUE),
                .groups = "drop"
              )
            
            # Format and add to output
            stats_output <- c(stats_output,
              paste(sprintf("%-20s", c("Group", "n", "mean", "median", "sd", "min", "max")), collapse = "\t"),
              "\n"
            )
            for (i in 1:nrow(nseg_stats)) {
              stats_output <- c(stats_output,
                paste(sprintf("%-20s", c(
                  as.character(nseg_stats[[group_col]][i]),
                  nseg_stats$n[i],
                  sprintf("%.2f", nseg_stats$mean[i]),
                  sprintf("%.2f", nseg_stats$median[i]),
                  sprintf("%.2f", nseg_stats$sd[i]),
                  sprintf("%.2f", nseg_stats$min[i]),
                  sprintf("%.2f", nseg_stats$max[i])
                )), collapse = "\t"),
                "\n"
              )
            }
            
            # Perform ANOVA
            stats_output <- c(stats_output, "\n### Number of ROH Segments ANOVA ###\n\n")
            tryCatch({
              # Create formula dynamically
              formula_str <- paste("N_ROH_segments ~ as.factor(", group_col, ")", sep = "")
              aov_nseg <- aov(as.formula(formula_str), data = data_subset)
              aov_summary <- summary(aov_nseg)
              stats_output <- c(stats_output, capture.output(print(aov_summary)), "\n")
              
              # Multiple comparison test (Tukey HSD) - always perform regardless of significance
              p_value <- aov_summary[[1]][["Pr(>F)"]][1]
              if (!is.na(p_value)) {
                if (p_value < 0.05) {
                  stats_output <- c(stats_output, paste0("\nANOVA p-value: ", sprintf("%.6f", p_value), " (significant at α=0.05)\n"))
                } else {
                  stats_output <- c(stats_output, paste0("\nANOVA p-value: ", sprintf("%.6f", p_value), " (not significant at α=0.05)\n"))
                }
                stats_output <- c(stats_output, "\n### Number of ROH Segments Multiple Comparisons (Tukey HSD) ###\n\n")
                tukey_nseg <- TukeyHSD(aov_nseg)
                stats_output <- c(stats_output, capture.output(print(tukey_nseg)), "\n")
                
                # Also get pairwise comparisons using multcomp
                tryCatch({
                  data_subset$group_factor <- as.factor(data_subset[[group_col]])
                  aov_nseg2 <- aov(N_ROH_segments ~ group_factor, data = data_subset)
                  glht_nseg <- glht(aov_nseg2, linfct = mcp(group_factor = "Tukey"))
                  summary_glht <- summary(glht_nseg)
                  stats_output <- c(stats_output, "\n### Number of ROH Segments Pairwise Comparisons (multcomp) ###\n\n")
                  stats_output <- c(stats_output, capture.output(print(summary_glht)), "\n")
                }, error = function(e) {
                  stats_output <<- c(stats_output, "Note: multcomp pairwise comparisons skipped (using TukeyHSD results above)\n")
                })
              } else {
                stats_output <- c(stats_output, "\nWarning: Could not calculate ANOVA p-value.\n\n")
              }
            }, error = function(e) {
              stats_output <<- c(stats_output, "Error in ANOVA:", e$message, "\n\n")
            })
          }
        }
      }
    }
  }
}

# If no grouping data available, write a message
if (length(stats_output) == 0) {
  stats_output <- c("No grouping data available for statistical comparisons.\n",
                    "Statistics require grouping columns specified in group_by parameter.\n")
}

# Write statistics to file
cat("Writing statistics and comparison tests to:", output_stats, "\n")
writeLines(stats_output, con = output_stats)

cat("  Statistics written to:", output_stats, "\n")
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

