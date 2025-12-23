# Load libraries
library(RColorBrewer)
library(ggplot2)
library(scales)

# Debug: working directory
message("Working directory: ", getwd())

# Debug: check if Snakemake object exists
if (!exists("snakemake")) {
  stop("Snakemake object not found! This script should be run via Snakemake.")
}

# Debug: print Snakemake inputs, outputs, and params
message("Snakemake inputs:")
print(snakemake@input)
message("Snakemake outputs:")
print(snakemake@output)
message("Snakemake params:")
print(snakemake@params)

# Function to generate all PC combinations
get_pca_combinations <- function(pc_max) {
  combinations <- list()
  for (i in 1:(pc_max)) {
    for (j in (i+1):(pc_max+1)) {
      if (j <= pc_max) {
        combinations[[length(combinations) + 1]] <- c(i, j)
      }
    }
  }
  return(combinations)
}

# PCA facet plotting function using facet_wrap
plot_pca_facet <- function(individuals, eigenvecs, eigenvals, popdata,
                           indmiss = NULL, color_by_name = NULL, pc_max = 2,
                           pca_colors = NULL, plot_type = "colored") {

  # Generate all PC combinations
  pc_combinations <- get_pca_combinations(pc_max)

  # Calculate variance explained
  var_exp <- eigenvals / sum(eigenvals) * 100

  # Build a long-format data frame with all PC combinations
  plot_data_list <- list()

  for (combo in pc_combinations) {
    pc1 <- combo[1]
    pc2 <- combo[2]

    # Create data frame for this PC combination
    combo_df <- data.frame(
      Ind = individuals,
      PC_X = eigenvecs[, pc1],
      PC_Y = eigenvecs[, pc2],
      PC_X_label = sprintf("PC%d", pc1),
      PC_Y_label = sprintf("PC%d", pc2),
      PC_X_var = sprintf("%.1f%%", var_exp[pc1]),
      PC_Y_var = sprintf("%.1f%%", var_exp[pc2]),
      stringsAsFactors = FALSE
    )

    plot_data_list[[length(plot_data_list) + 1]] <- combo_df
  }

  # Combine all combinations into one data frame
  plot_data <- do.call(rbind, plot_data_list)

  # Create combined labels for faceting
  plot_data$PC_X_full <- paste0(plot_data$PC_X_label, " (", plot_data$PC_X_var, ")")
  plot_data$PC_Y_full <- paste0(plot_data$PC_Y_label, " (", plot_data$PC_Y_var, ")")

  # Convert to factors to control order
  plot_data$PC_X_full <- factor(plot_data$PC_X_full, levels = unique(plot_data$PC_X_full))
  plot_data$PC_Y_full <- factor(plot_data$PC_Y_full, levels = unique(plot_data$PC_Y_full))

  # Generate plot based on plot_type
  if (plot_type == "labeled") {
    # Labels only: gray points + sample labels, no colors
    p <- ggplot(plot_data, aes(x = PC_X, y = PC_Y)) +
      geom_point(size = 2, alpha = 0.7, pch = 21, fill = "gray80", color = "black") +
      geom_text(aes(label = Ind), size = 1, check_overlap = FALSE) +
      facet_grid(PC_Y_full ~ PC_X_full, scales = "free") +
      labs(x = NULL, y = NULL) +
      theme_bw() +
      theme(
        strip.text = element_text(size = 9, face = "bold"),
        strip.background = element_rect(fill = "gray90")
      )

  } else if (plot_type == "missing") {
    # Color by missing data rate
    if (is.null(indmiss)) {
      stop("plot_type='missing' but no missing data provided")
    }

    plot_data <- merge(plot_data, indmiss[, c("INDV", "F_MISS")],
                       by.x = "Ind", by.y = "INDV", all.x = TRUE)

    p <- ggplot(plot_data, aes(x = PC_X, y = PC_Y, fill = F_MISS)) +
      geom_point(size = 2, alpha = 0.7, pch = 21, color = "black") +
      scale_fill_distiller(
        name = "Missing Data",
        labels = scales::percent,
        palette = "YlOrBr",
        direction = 1
      ) +
      facet_grid(PC_Y_full ~ PC_X_full, scales = "free") +
      labs(x = NULL, y = NULL) +
      theme_bw() +
      theme(
        strip.text = element_text(size = 9, face = "bold"),
        strip.background = element_rect(fill = "gray90"),
        legend.position = "right"
      )

  } else {
    # Colored by population/metadata
    if (is.null(color_by_name) || color_by_name == "" || color_by_name == "none") {
      # No coloring - just gray points
      p <- ggplot(plot_data, aes(x = PC_X, y = PC_Y)) +
        geom_point(size = 2, alpha = 0.7, pch = 21, fill = "gray80", color = "black") +
        facet_grid(PC_Y_full ~ PC_X_full, scales = "free") +
        labs(x = NULL, y = NULL) +
        theme_bw() +
        theme(
          strip.text = element_text(size = 9, face = "bold"),
          strip.background = element_rect(fill = "gray90")
        )
    } else {
      # Categorical coloring by population metadata
      ind_col <- colnames(popdata)[1]
      plot_data <- merge(plot_data, popdata, by.x = "Ind", by.y = ind_col, all.x = TRUE)

      # Remove samples without metadata
      if (any(is.na(plot_data[[2]]))) {
        plot_data <- plot_data[!is.na(plot_data[[2]]), ]
      }

      # Get color column
      if (!is.null(color_by_name) && color_by_name %in% names(popdata)) {
        color_col <- color_by_name
      } else {
        stop(paste("Column", color_by_name, "not found in popdata"))
      }

      # Remove NA and empty values from color column
      na_mask <- is.na(plot_data[[color_col]])
      if (any(na_mask)) {
        plot_data <- plot_data[!na_mask, ]
      }
      empty_mask <- plot_data[[color_col]] == ""
      if (any(empty_mask)) {
        plot_data <- plot_data[!empty_mask, ]
      }

      # Determine colors to use
      n_categories <- length(unique(plot_data[[color_col]]))

      if (!is.null(pca_colors) && length(pca_colors) > 0) {
        if (n_categories > length(pca_colors)) {
          colors <- colorRampPalette(pca_colors)(n_categories)
        } else {
          colors <- pca_colors[1:n_categories]
        }
      } else {
        colors <- NULL
      }

      p <- ggplot(plot_data, aes(x = PC_X, y = PC_Y, fill = .data[[color_col]])) +
        geom_point(size = 2, alpha = 0.7, pch = 21, color = "black") +
        facet_grid(PC_Y_full ~ PC_X_full, scales = "free") +
        labs(x = NULL, y = NULL, fill = color_col) +
        theme_bw() +
        theme(
          strip.text = element_text(size = 9, face = "bold"),
          strip.background = element_rect(fill = "gray90"),
          legend.position = "right"
        )

      if (!is.null(colors)) {
        p <- p + scale_fill_manual(values = colors)
      }
    }
  }

  return(p)
}

message("plot_pca_facet function loaded")

# Read Snakemake inputs
eigvecs_file <- snakemake@input[["eigvecs"]]
eigvals_file <- snakemake@input[["eigvals"]]
popdata_file <- snakemake@input[["indpopdata"]]
output_pdf   <- snakemake@output[["pdf"]]
output_rds   <- snakemake@output[["rds"]]

# Read parameters
pc_max <- as.numeric(snakemake@params[["pc_max"]])
plot_type <- as.character(snakemake@params[["plot_type"]])

# Get color_by and pca_colors only for colored plots
color_by_name <- NULL
pca_colors <- NULL
if (plot_type == "colored") {
  if ("color_by" %in% names(snakemake@params)) {
    color_by_name <- as.character(snakemake@params[["color_by"]])
  }
  if ("pca_colors" %in% names(snakemake@params)) {
    pca_colors <- snakemake@params[["pca_colors"]]
    if (!is.null(pca_colors)) {
      pca_colors <- unlist(pca_colors)
    }
  }
}

# Debug: check param values
message("pc_max = ", pc_max)
message("plot_type = ", plot_type)
if (!is.null(color_by_name)) {
  message("color_by_name = ", color_by_name)
}
if (!is.null(pca_colors)) {
  message("pca_colors = ", paste(pca_colors, collapse = ", "))
}

# Read input data
eigenvecs2 <- read.table(eigvecs_file, sep = "\t", skip = 1)
individuals <- eigenvecs2$V1
eigenvecs <- eigenvecs2[, -c(1,2)]
eigenvals <- scan(eigvals_file, quiet = TRUE)
popdata <- read.table(popdata_file, header = TRUE, sep = "\t")

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
    }
  }
}

# Generate PCA facet plot
plt_pca_facet <- plot_pca_facet(
  individuals = individuals,
  eigenvecs = eigenvecs,
  eigenvals = eigenvals,
  popdata = popdata,
  indmiss = indmiss,
  color_by_name = color_by_name,
  pc_max = pc_max,
  pca_colors = pca_colors,
  plot_type = plot_type
)

# Debug: check plot object
message("Generated PCA facet plot object: ", class(plt_pca_facet))

# Calculate dimensions based on number of plots
n_plots <- length(get_pca_combinations(pc_max))
ncols <- ceiling(sqrt(n_plots))
nrows <- ceiling(n_plots / ncols)
plot_width <- ncols * 4
plot_height <- nrows * 3.5

# Save PDF
dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)
ggplot2::ggsave(output_pdf, plt_pca_facet, width = plot_width, height = plot_height, dpi = 300)
message("PDF saved to ", output_pdf)

# Save RDS (ggplot2 object)
dir.create(dirname(output_rds), recursive = TRUE, showWarnings = FALSE)
saveRDS(plt_pca_facet, output_rds)
message("RDS saved to ", output_rds)

message("Script completed successfully")
