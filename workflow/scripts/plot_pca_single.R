# Debug: check if Snakemake object exists
if (!exists("snakemake")) {
  stop("Snakemake object not found! This script should be run via Snakemake.")
}

# Redirect all script output/messages to the rule log.
log_file <- snakemake@log[[1]]
dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")
on.exit({
  while (sink.number(type = "message") > 0) sink(type = "message")
  while (sink.number(type = "output") > 0) sink(type = "output")
  close(log_con)
}, add = TRUE)

# Load libraries
suppressPackageStartupMessages({
  library(RColorBrewer)
  library(ggplot2)
  library(scales)
})

ggsave_utils <- tryCatch(
  file.path(dirname(normalizePath(snakemake@script)), "plot_ggsave_utils.R"),
  error = function(e) "workflow/scripts/plot_ggsave_utils.R"
)
if (file.exists(ggsave_utils)) {
  source(ggsave_utils)
} else {
  source("workflow/scripts/plot_ggsave_utils.R")
}

group_utils <- tryCatch(
  file.path(dirname(normalizePath(snakemake@script)), "plot_group_utils.R"),
  error = function(e) "workflow/scripts/plot_group_utils.R"
)
if (file.exists(group_utils)) {
  source(group_utils)
} else {
  source("workflow/scripts/plot_group_utils.R")
}

# Debug: working directory
message("Working directory: ", getwd())

# Debug: print Snakemake inputs, outputs, and params
message("Snakemake inputs:")
print(snakemake@input)
message("Snakemake outputs:")
print(snakemake@output)
message("Snakemake params:")
print(snakemake@params)

# Enhanced PCA plotting function
plot_pca <- function(individuals, eigenvecs, eigenvals, popdata,
                     indmiss = NULL, color_by_name = NULL, pc1 = 1, pc2 = 2,
                     group_colors = NULL, plot_type = "colored",
                     point_size = 3, axis_title_size = 10, axis_text_size = 8) {

  # Create data frame with PC scores
  pca_df <- data.frame(
    Ind = individuals,
    PC1 = eigenvecs[, pc1],
    PC2 = eigenvecs[, pc2]
  )

  # Calculate variance explained
  var_exp <- eigenvals / sum(eigenvals) * 100

  axis_theme <- theme(
    axis.title = element_text(size = axis_title_size),
    axis.text = element_text(size = axis_text_size)
  )

  # Generate plot based on plot_type
  if (plot_type == "labeled") {
    # Labels only: gray points + sample labels, no colors
    message("Generating labeled plot (gray points + sample names)")
    p <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
      geom_point(size = point_size, alpha = 0.7, pch = 21, fill = "gray80", color = "black") +
      geom_text(aes(label = Ind), size = 1.5, check_overlap = FALSE) +
      labs(
        x = sprintf("PC%d (%.1f%%)", pc1, var_exp[pc1]),
        y = sprintf("PC%d (%.1f%%)", pc2, var_exp[pc2])
      ) +
      theme_bw() +
      axis_theme

  } else if (plot_type == "missing") {
    # Color by missing data rate
    message("Generating missing data plot")
    if (is.null(indmiss)) {
      stop("plot_type='missing' but no missing data provided")
    }

    plot_df <- merge(pca_df, indmiss[, c("INDV", "F_MISS")],
                     by.x = "Ind", by.y = "INDV", all.x = TRUE)

    # Check for samples without missing data info
    if (any(is.na(plot_df$F_MISS))) {
      warning(sprintf("%d samples lack missing data information", sum(is.na(plot_df$F_MISS))))
    }

    p <- ggplot(plot_df, aes(x = PC1, y = PC2, fill = F_MISS)) +
      geom_point(size = point_size, alpha = 0.7, pch = 21, color = "black") +
      scale_fill_distiller(
        name = "Missing Data",
        labels = scales::percent,
        palette = "YlOrBr",
        direction = 1
      ) +
      labs(
        x = sprintf("PC%d (%.1f%%)", pc1, var_exp[pc1]),
        y = sprintf("PC%d (%.1f%%)", pc2, var_exp[pc2])
      ) +
      theme_bw() +
      theme(legend.position = "right") +
      axis_theme

  } else {
    # Colored by population/metadata
    message("Generating colored plot")

    # Check if color_by is specified and valid
    if (is.null(color_by_name) || color_by_name == "" || color_by_name == "none") {
      # No coloring - just gray points
      message("No color_by specified - plotting with gray points only")
      p <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
        geom_point(size = point_size, alpha = 0.7, pch = 21, fill = "gray80", color = "black") +
        labs(
          x = sprintf("PC%d (%.1f%%)", pc1, var_exp[pc1]),
          y = sprintf("PC%d (%.1f%%)", pc2, var_exp[pc2])
        ) +
        theme_bw() +
        axis_theme
    } else {
      # Original categorical coloring by population metadata
      ind_col <- colnames(popdata)[1]

      # Debug: show sample counts before merge
      message(sprintf("PCA samples: %d, Popdata samples: %d", nrow(pca_df), nrow(popdata)))

      plot_df <- merge(pca_df, popdata, by.x = "Ind", by.y = ind_col, all.x = TRUE)

      message(sprintf("After merge: %d samples", nrow(plot_df)))

      # Check for samples without metadata
      if (any(is.na(plot_df[[2]]))) {
        n_missing <- sum(is.na(plot_df[[2]]))
        missing_samples <- plot_df$Ind[is.na(plot_df[[2]])]
        warning(sprintf("%d samples lack population metadata and will be excluded from plot:", n_missing))
        warning(paste(head(missing_samples, 10), collapse = ", "))
        if (n_missing > 10) warning(sprintf("... and %d more", n_missing - 10))
        # Remove samples without metadata
        plot_df <- plot_df[!is.na(plot_df[[2]]), ]
        message(sprintf("After removing samples without metadata: %d samples", nrow(plot_df)))
      }

      # Get color column - handle both name and index
      if (!is.null(color_by_name) && color_by_name %in% names(popdata)) {
        color_col <- color_by_name
      } else {
        stop(paste("Column", color_by_name, "not found in popdata"))
      }

      # Check unique values in color column (including NA and empty strings)
      unique_vals <- unique(plot_df[[color_col]])
      message(sprintf("Unique values in '%s': %s", color_col, paste(unique_vals, collapse = ", ")))
      message(sprintf("Number of unique values: %d", length(unique_vals)))

      # Remove NA values from the color column and warn user
      na_mask <- is.na(plot_df[[color_col]])
      if (any(na_mask)) {
        n_na <- sum(na_mask)
        na_samples <- plot_df$Ind[na_mask]
        warning(sprintf("%d samples have missing values in '%s' column and will be excluded:", n_na, color_col))
        warning(paste(head(na_samples, 10), collapse = ", "))
        if (n_na > 10) warning(sprintf("... and %d more", n_na - 10))
        plot_df <- plot_df[!na_mask, ]
      }

      # Remove empty strings from the color column
      empty_mask <- plot_df[[color_col]] == ""
      if (any(empty_mask)) {
        n_empty <- sum(empty_mask)
        empty_samples <- plot_df$Ind[empty_mask]
        warning(sprintf("%d samples have empty strings in '%s' column and will be excluded:", n_empty, color_col))
        warning(paste(head(empty_samples, 10), collapse = ", "))
        if (n_empty > 10) warning(sprintf("... and %d more", n_empty - 10))
        plot_df <- plot_df[!empty_mask, ]
      }

      color_is_numeric <- is.numeric(plot_df[[color_col]])
      if (!color_is_numeric) {
        plot_df[[color_col]] <- as.character(plot_df[[color_col]])
      }

      p <- ggplot(plot_df, aes(x = PC1, y = PC2, fill = .data[[color_col]])) +
        geom_point(size = point_size, alpha = 0.7, pch = 21, color = "black") +
        labs(
          x = sprintf("PC%d (%.1f%%)", pc1, var_exp[pc1]),
          y = sprintf("PC%d (%.1f%%)", pc2, var_exp[pc2]),
          fill = color_col
        ) +
        theme_bw() +
        axis_theme

      if (color_is_numeric) {
        message(sprintf("PCA plot: using continuous viridis scale for numeric '%s'", color_col))
        p <- p + scale_fill_viridis_c(name = color_col, option = "D")
      } else {
        colors <- group_fill_values(group_colors)
        if (is.null(colors)) {
          message("PCA plot: No colors defined in config. Using default ggplot2 colors.")
        } else {
          p <- p + scale_fill_manual(values = colors)
        }
      }
    }
  }

  return(p)
}

message("plot_pca function loaded")

# Read Snakemake inputs
eigvecs_file <- snakemake@input[["eigvecs"]]
eigvals_file <- snakemake@input[["eigvals"]]
popdata_file <- snakemake@input[["indpopdata"]]
output_pdf   <- snakemake@output[["pdf"]]
output_rds   <- snakemake@output[["rds"]]

# Read parameters
pc1 <- as.numeric(snakemake@params[["pc1"]])
pc2 <- as.numeric(snakemake@params[["pc2"]])
plot_type <- as.character(snakemake@params[["plot_type"]])
point_size <- as.numeric(snakemake@params[["point_size"]])
axis_title_size <- as.numeric(snakemake@params[["axis_title_size"]])
axis_text_size <- as.numeric(snakemake@params[["axis_text_size"]])
if (is.na(point_size)) point_size <- 3
if (is.na(axis_title_size)) axis_title_size <- 10
if (is.na(axis_text_size)) axis_text_size <- 8

# Get color_by and group_colors only for colored plots
color_by_name <- NULL
group_colors <- NULL
if (plot_type == "colored") {
  if ("color_by" %in% names(snakemake@params)) {
    color_by_name <- as.character(snakemake@params[["color_by"]])
  }
  if ("group_colors" %in% names(snakemake@params)) {
    group_colors <- snakemake@params[["group_colors"]]
  }
}

# Debug: check param values
message("pc1 = ", pc1)
message("pc2 = ", pc2)
message("plot_type = ", plot_type)
message("point_size = ", point_size)
message("axis_title_size = ", axis_title_size)
message("axis_text_size = ", axis_text_size)
if (!is.null(color_by_name)) {
  message("color_by_name = ", color_by_name)
}
if (!is.null(group_colors)) {
  message("group_colors configured for PCA plot")
}

# Read input data
# Read eigenvectors file - skip the #FID header line and use column positions
eigenvecs2 <- read.table(eigvecs_file, sep = "\t", skip = 1)
individuals <- eigenvecs2$V1
eigenvecs <- eigenvecs2[, -c(1,2)]
eigenvals <- scan(eigvals_file, quiet = TRUE)
popdata <- read.table(popdata_file, header = TRUE, sep = "\t")

# Debug: show first few individuals
message("First 5 individuals from eigvecs: ", paste(head(individuals, 5), collapse = ", "))
message("First 5 individuals from popdata: ", paste(head(popdata[[1]], 5), collapse = ", "))

# Read missing data if available
indmiss <- NULL
indmiss_file <- NULL
if ("indmiss" %in% names(snakemake@input)) {
  indmiss_input <- snakemake@input[["indmiss"]]
  if (length(indmiss_input) > 0 && indmiss_input != "") {
    indmiss_file <- indmiss_input
    message("Reading missing data from: ", indmiss_file)
    if (file.exists(indmiss_file)) {
      indmiss <- read.table(indmiss_file, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
      message("Missing data loaded: ", nrow(indmiss), " individuals")
      message("Columns: ", paste(names(indmiss), collapse = ", "))
    }
  }
}

# Debug: show columns in popdata
message("popdata columns: ", paste(names(popdata), collapse = ", "))

# Generate PCA plot
plt_pca <- plot_pca(
  individuals = individuals,
  eigenvecs = eigenvecs,
  eigenvals = eigenvals,
  popdata = popdata,
  indmiss = indmiss,
  color_by_name = color_by_name,
  pc1 = pc1,
  pc2 = pc2,
  group_colors = group_colors,
  plot_type = plot_type,
  point_size = point_size,
  axis_title_size = axis_title_size,
  axis_text_size = axis_text_size
)

# Debug: check plot object
message("Generated PCA plot object: ", class(plt_pca))

# Save PDF
dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)
ggsave_pdf(output_pdf, plt_pca, width = 6, height = 5, dpi = 300)
message("PDF saved to ", output_pdf)

# Save RDS (ggplot2 object)
dir.create(dirname(output_rds), recursive = TRUE, showWarnings = FALSE)
saveRDS(plt_pca, output_rds)
message("RDS saved to ", output_rds)

message("Script completed successfully")
