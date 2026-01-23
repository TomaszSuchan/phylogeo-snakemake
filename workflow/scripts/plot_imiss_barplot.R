#!/usr/bin/env Rscript
# Create imiss barplot, optionally grouped by stratification (similar to pixy pi barplots)

library(ggplot2)
library(dplyr)

# Prevent creation of Rplots.pdf
pdf(NULL)

# Redirect all output to log file
log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

# Snakemake inputs/outputs
imiss_file <- snakemake@input[["imiss"]]
popdata_file <- snakemake@input[["popdata"]]
output_pdf <- snakemake@output[["pdf"]]
output_rds <- snakemake@output[["rds"]]

# Parameters
color_by_name <- NULL
pca_colors <- NULL
plot_type <- "grouped"  # "grouped", "sorted", or "plain"
if ("color_by" %in% names(snakemake@params)) {
  color_by_param <- snakemake@params[["color_by"]]
  if (!is.null(color_by_param) && !is.na(color_by_param)) {
    color_by_name <- as.character(color_by_param)
    if (length(color_by_name) == 0 || color_by_name == "" || color_by_name == "none") {
      color_by_name <- NULL
    }
  } else {
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

message("\n=== READING IMISS DATA ===\n")
imiss_df <- read.table(imiss_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
message(sprintf("Loaded %d individuals\n", nrow(imiss_df)))

# Read popdata if stratification is requested
if (!is.null(color_by_name) && plot_type != "plain") {
  message("\n=== READING POPDATA ===\n")
  popdata <- read.table(popdata_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  message(sprintf("Popdata columns: %s\n", paste(names(popdata), collapse = ", ")))
  
  # Merge imiss with popdata by individual ID
  # Try matching INDV column with Ind column
  imiss_df_merged <- merge(imiss_df, popdata, 
                          by.x = "INDV", 
                          by.y = "Ind", 
                          all.x = TRUE)
  
  # Check if color_by column exists
  if (!color_by_name %in% names(imiss_df_merged)) {
    warning(sprintf("Column '%s' not found in popdata. Plotting without grouping.\n", color_by_name))
    color_by_name <- NULL
  } else {
    # Remove individuals without stratification data
    na_mask <- is.na(imiss_df_merged[[color_by_name]])
    if (any(na_mask)) {
      n_na <- sum(na_mask)
      warning(sprintf("%d individuals have missing values in '%s' column and will be excluded\n", 
                     n_na, color_by_name))
      imiss_df_merged <- imiss_df_merged[!na_mask, ]
    }
    
    imiss_df <- imiss_df_merged
    message(sprintf("After merging with popdata: %d individuals\n", nrow(imiss_df)))
  }
}

message("\n=== CREATING BARPLOT ===\n")

# Create base plot
if (is.null(color_by_name) || plot_type == "plain") {
  # No grouping - simple barplot, ordered by missingness
  imiss_df <- imiss_df[order(imiss_df$F_MISS), ]
  imiss_df$INDV <- factor(imiss_df$INDV, levels = imiss_df$INDV)
  
  p <- ggplot(imiss_df, aes(x = INDV, y = F_MISS)) +
    geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7, color = "black", linewidth = 0.3) +
    labs(
      x = "Individual",
      y = "Proportion of missing data"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      panel.grid.major.x = element_blank()
    )
} else {
  # Colored by stratification
  unique_categories <- sort(unique(imiss_df[[color_by_name]]))
  n_categories <- length(unique_categories)
  
  if (!is.null(pca_colors) && length(pca_colors) > 0) {
    if (n_categories > length(pca_colors)) {
      message(sprintf("Imiss barplot: %d categories but only %d colors defined. Interpolating additional colors.",
                     n_categories, length(pca_colors)))
      colors_all <- colorRampPalette(pca_colors)(n_categories)
    } else {
      colors_all <- pca_colors[1:n_categories]
    }
    names(colors_all) <- unique_categories
    colors <- colors_all
  } else {
    message("Imiss barplot: No colors defined in config. Using default ggplot2 colors.")
    colors <- NULL
  }
  
  if (plot_type == "sorted") {
    # Sorted by missingness (low to high), colored by stratification
    imiss_df <- imiss_df[order(imiss_df$F_MISS), ]
    imiss_df$INDV <- factor(imiss_df$INDV, levels = imiss_df$INDV)
    
    p <- ggplot(imiss_df, aes(x = INDV, y = F_MISS, fill = .data[[color_by_name]])) +
      geom_bar(stat = "identity", alpha = 0.7, color = "black", linewidth = 0.3) +
      labs(
        x = "Individual",
        y = "Proportion of missing data",
        fill = color_by_name
      ) +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        panel.grid.major.x = element_blank(),
        legend.position = "right"
      )
    
    if (!is.null(colors)) {
      p <- p + scale_fill_manual(values = colors)
    }
    
  } else {
    # Grouped by stratification
    imiss_df <- imiss_df %>%
      group_by(.data[[color_by_name]]) %>%
      arrange(F_MISS) %>%
      ungroup()
    
    imiss_df$INDV <- factor(imiss_df$INDV, levels = unique(imiss_df$INDV))
    
    p <- ggplot(imiss_df, aes(x = INDV, y = F_MISS, fill = .data[[color_by_name]])) +
      geom_bar(stat = "identity", alpha = 0.7, color = "black", linewidth = 0.3) +
      facet_wrap(as.formula(paste("~", color_by_name)), scales = "free_x", ncol = 4) +
      labs(
        x = "Individual",
        y = "Proportion of missing data",
        fill = color_by_name
      ) +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        panel.grid.major.x = element_blank(),
        legend.position = "bottom",
        strip.text = element_text(size = 10, face = "bold")
      )
    
    if (!is.null(colors)) {
      p <- p + scale_fill_manual(values = colors)
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
  width = max(8, nrow(imiss_df) * 0.2),
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

