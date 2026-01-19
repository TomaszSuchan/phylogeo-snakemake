#!/usr/bin/env Rscript
# Create Pi barplot with confidence intervals, optionally grouped by stratification

library(ggplot2)
library(dplyr)

# Prevent creation of Rplots.pdf
pdf(NULL)

# Redirect all output to log file
log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

# Snakemake inputs/outputs
pi_summary_file <- snakemake@input[["pi_summary"]]
popdata_file <- snakemake@input[["popdata"]]
output_pdf <- snakemake@output[["pdf"]]
output_rds <- snakemake@output[["rds"]]

# Parameters
color_by_name <- NULL
pca_colors <- NULL
plot_type <- "grouped"  # "grouped" or "sorted"
if ("color_by" %in% names(snakemake@params)) {
  color_by_name <- as.character(snakemake@params[["color_by"]])
  if (color_by_name == "" || color_by_name == "none") {
    color_by_name <- NULL
  }
}
if ("pca_colors" %in% names(snakemake@params)) {
  pca_colors <- snakemake@params[["pca_colors"]]
  if (!is.null(pca_colors)) {
    pca_colors <- unlist(pca_colors)
  }
}
if ("plot_type" %in% names(snakemake@params)) {
  plot_type <- as.character(snakemake@params[["plot_type"]])
}

message("\n=== READING PI SUMMARY ===\n")
pi_df <- read.table(pi_summary_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
message(sprintf("Loaded %d populations\n", nrow(pi_df)))

# Read popdata if stratification is requested
if (!is.null(color_by_name)) {
  message("\n=== READING POPDATA ===\n")
  popdata <- read.table(popdata_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  message(sprintf("Popdata columns: %s\n", paste(names(popdata), collapse = ", ")))
  
  # Get first column name (should be Site/population identifier)
  site_col <- colnames(popdata)[1]
  
  # Get unique population-level data from popdata
  # Aggregate by Site to get population-level metadata
  popdata_unique <- popdata[, c(site_col, color_by_name), drop = FALSE]
  popdata_unique <- popdata_unique[!duplicated(popdata_unique[[site_col]]), ]
  
  message(sprintf("Unique populations in popdata: %d\n", nrow(popdata_unique)))
  
  # Try matching population names
  # First try: direct match with full population name (e.g., "Csr_BBE")
  pi_df_merged <- merge(pi_df, popdata_unique, 
                        by.x = "population", 
                        by.y = site_col, 
                        all.x = TRUE)
  
  # If that fails for some populations, try removing project prefix (e.g., "BBE")
  if (any(is.na(pi_df_merged[[color_by_name]]))) {
    # Only try alternative matching for populations that didn't match
    unmatched <- is.na(pi_df_merged[[color_by_name]])
    pi_df_unmatched <- pi_df[unmatched, ]
    
    if (nrow(pi_df_unmatched) > 0) {
      pi_df_unmatched$population_clean <- gsub("^[^_]+_", "", pi_df_unmatched$population)
      pi_df_unmatched_merged <- merge(pi_df_unmatched, popdata_unique, 
                                      by.x = "population_clean", 
                                      by.y = site_col, 
                                      all.x = TRUE)
      # Remove temporary column
      pi_df_unmatched_merged$population_clean <- NULL
      
      # Combine matched and unmatched results
      pi_df_merged[unmatched, ] <- pi_df_unmatched_merged
    }
  }
  
  # Check if color_by column exists
  if (!color_by_name %in% names(pi_df_merged)) {
    warning(sprintf("Column '%s' not found in popdata. Plotting without grouping.\n", color_by_name))
    color_by_name <- NULL
  } else {
    # Remove populations without stratification data
    na_mask <- is.na(pi_df_merged[[color_by_name]])
    if (any(na_mask)) {
      n_na <- sum(na_mask)
      warning(sprintf("%d populations have missing values in '%s' column and will be excluded\n", 
                     n_na, color_by_name))
      pi_df_merged <- pi_df_merged[!na_mask, ]
    }
    
    # Remove empty strings
    empty_mask <- pi_df_merged[[color_by_name]] == ""
    if (any(empty_mask)) {
      n_empty <- sum(empty_mask)
      warning(sprintf("%d populations have empty strings in '%s' column and will be excluded\n", 
                     n_empty, color_by_name))
      pi_df_merged <- pi_df_merged[!empty_mask, ]
    }
    
    pi_df <- pi_df_merged
    message(sprintf("After merging with popdata: %d populations\n", nrow(pi_df)))
  }
}

message("\n=== CREATING BARPLOT ===\n")

# Create base plot
if (is.null(color_by_name)) {
  # No grouping - simple barplot, ordered by pi value
  pi_df <- pi_df[order(pi_df$mean_pi), ]
  pi_df$population <- factor(pi_df$population, levels = pi_df$population)
  
  p <- ggplot(pi_df, aes(x = population, y = mean_pi)) +
    geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7, color = "black", linewidth = 0.3) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), 
                  width = 0.2, color = "black", linewidth = 0.5) +
    labs(
      x = "Population",
      y = expression(paste("Nucleotide Diversity (", pi, ")"))
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      panel.grid.major.x = element_blank()
    )
} else {
  # Colored by stratification
  # Determine colors to use - ensure consistent color mapping across plot types
  unique_categories <- sort(unique(pi_df[[color_by_name]]))  # Sort to ensure consistent order
  n_categories <- length(unique_categories)
  
  if (!is.null(pca_colors) && length(pca_colors) > 0) {
    # Use colors from config
    if (n_categories > length(pca_colors)) {
      # More categories than colors - interpolate using colorRampPalette
      message(sprintf("Pi barplot: %d categories but only %d colors defined. Interpolating additional colors.",
                     n_categories, length(pca_colors)))
      colors_all <- colorRampPalette(pca_colors)(n_categories)
    } else {
      # Use subset of defined colors
      colors_all <- pca_colors[1:n_categories]
    }
    # Create named vector for consistent color mapping
    names(colors_all) <- unique_categories
    colors <- colors_all
  } else {
    # Fallback to default ggplot2 colors if no colors provided
    message("Pi barplot: No colors defined in config. Using default ggplot2 colors.")
    colors <- NULL
  }
  
  if (plot_type == "sorted") {
    # Sorted by pi value (low to high), colored by stratification
    pi_df <- pi_df[order(pi_df$mean_pi), ]
    pi_df$population <- factor(pi_df$population, levels = pi_df$population)
    
    p <- ggplot(pi_df, aes(x = population, y = mean_pi, fill = .data[[color_by_name]])) +
      geom_bar(stat = "identity", alpha = 0.7, color = "black", linewidth = 0.3) +
      geom_errorbar(aes(ymin = ci_low, ymax = ci_high), 
                    width = 0.2, color = "black", linewidth = 0.5) +
      labs(
        x = "Population",
        y = expression(paste("Nucleotide Diversity (", pi, ")")),
        fill = color_by_name
      ) +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        panel.grid.major.x = element_blank(),
        legend.position = "right"
      )
    
    # Add color scale if colors were defined
    if (!is.null(colors)) {
      p <- p + scale_fill_manual(values = colors)
    }
    
  } else {
    # Grouped by stratification - order by group, then by pi within group
    pi_df <- pi_df[order(pi_df[[color_by_name]], pi_df$mean_pi), ]
    
    # Create factor levels to maintain order
    pi_df$population <- factor(pi_df$population, levels = pi_df$population)
    
    # Count populations per group for separators
    group_counts <- table(pi_df[[color_by_name]])
    group_boundaries <- cumsum(group_counts)
    
    p <- ggplot(pi_df, aes(x = population, y = mean_pi, fill = .data[[color_by_name]])) +
      geom_bar(stat = "identity", alpha = 0.7, color = "black", linewidth = 0.3) +
      geom_errorbar(aes(ymin = ci_low, ymax = ci_high), 
                    width = 0.2, color = "black", linewidth = 0.5) +
      labs(
        x = "Population",
        y = expression(paste("Nucleotide Diversity (", pi, ")")),
        fill = color_by_name
      ) +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        panel.grid.major.x = element_blank(),
        legend.position = "right"
      )
    
    # Add color scale if colors were defined (use same colors as sorted version)
    if (!is.null(colors)) {
      p <- p + scale_fill_manual(values = colors)
    }
    
    # Add vertical lines between groups
    for (boundary in group_boundaries[-length(group_boundaries)]) {
      p <- p + geom_vline(xintercept = boundary + 0.5, 
                         color = "gray50", linetype = "dashed", linewidth = 0.5)
    }
  }
}

# Save plots
message("\n=== SAVING OUTPUT ===\n")
message(sprintf("Output PDF: %s\n", output_pdf))
message(sprintf("Output RDS: %s\n", output_rds))

dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)
ggsave(
  filename = output_pdf,
  plot = p,
  width = max(8, nrow(pi_df) * 0.3),  # Adjust width based on number of populations
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

