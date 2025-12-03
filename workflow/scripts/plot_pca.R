# Load libraries
library(RColorBrewer)
library(ggplot2)
library(viridis)
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

# Enhanced PCA plotting function with missing data support
plot_pca <- function(individuals, eigenvecs, eigenvals, popdata,
                     indmiss = NULL, color_by_name = NULL, pc1 = 1, pc2 = 2) {

  # Create data frame with PC scores
  pca_df <- data.frame(
    Ind = individuals,
    PC1 = eigenvecs[, pc1],
    PC2 = eigenvecs[, pc2]
  )

  # Calculate variance explained
  var_exp <- eigenvals / sum(eigenvals) * 100

  # Check if we should color by missing data
  if (!is.null(color_by_name) && color_by_name == "missing") {
    # Color by missing data rate
    if (is.null(indmiss)) {
      stop("color_by='missing' but no missing data provided")
    }

    plot_df <- merge(pca_df, indmiss[, c("INDV", "F_MISS")],
                     by.x = "Ind", by.y = "INDV", all.x = TRUE)

    # Check for samples without missing data info
    if (any(is.na(plot_df$F_MISS))) {
      warning(sprintf("%d samples lack missing data information", sum(is.na(plot_df$F_MISS))))
    }

    p <- ggplot(plot_df, aes(x = PC1, y = PC2, color = F_MISS)) +
      geom_point(size = 3, alpha = 0.7) +
      scale_color_viridis_c(
        name = "Missing Data",
        labels = scales::percent,
        option = "plasma",
        direction = -1  # Higher missing = warmer colors
      ) +
      labs(
        x = sprintf("PC%d (%.1f%%)", pc1, var_exp[pc1]),
        y = sprintf("PC%d (%.1f%%)", pc2, var_exp[pc2])
      ) +
      theme_bw() +
      theme(legend.position = "right")

  } else {
    # Original categorical coloring by population metadata
    ind_col <- colnames(popdata)[1]
    plot_df <- merge(pca_df, popdata, by.x = "Ind", by.y = ind_col)

    # Get color column - handle both name and index
    if (!is.null(color_by_name) && color_by_name %in% names(popdata)) {
      color_col <- color_by_name
    } else {
      stop(paste("Column", color_by_name, "not found in popdata"))
    }

    p <- ggplot(plot_df, aes(x = PC1, y = PC2, color = .data[[color_col]])) +
      geom_point(size = 3, alpha = 0.7) +
      labs(
        x = sprintf("PC%d (%.1f%%)", pc1, var_exp[pc1]),
        y = sprintf("PC%d (%.1f%%)", pc2, var_exp[pc2]),
        color = color_col
      ) +
      theme_bw()
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

# Read optional missing data input
indmiss_file <- NULL

if ("indmiss" %in% names(snakemake@input)) {
  indmiss_input <- snakemake@input[["indmiss"]]
  # Check if input is non-empty (Snakemake passes empty list [] when not needed)
  if (length(indmiss_input) > 0 && indmiss_input != "") {
    indmiss_file <- indmiss_input
    message("Missing data file specified: ", indmiss_file)
  } else {
    message("No missing data file provided (empty input)")
  }
} else {
  message("No 'indmiss' input in snakemake object")
}

# Debug: check input files exist
message("eigvecs_file exists? ", file.exists(eigvecs_file))
message("eigvals_file exists? ", file.exists(eigvals_file))
message("popdata_file exists? ", file.exists(popdata_file))

# Read parameters and coerce types safely
pc1 <- as.numeric(snakemake@params[["pc1"]])
pc2 <- as.numeric(snakemake@params[["pc2"]])
color_by_name <- as.character(snakemake@params[["color_by"]])

# Debug: check param values
message("pc1 = ", pc1)
message("pc2 = ", pc2)
message("color_by_name = ", color_by_name)

# Read input data
eigenvecs2 <- read.table(eigvecs_file)
individuals <- eigenvecs2$V1
eigenvecs <- eigenvecs2[, -c(1,2)]
eigenvals <- scan(eigvals_file)
popdata <- read.table(popdata_file, header = TRUE)

# Read missing data if available
indmiss <- NULL
if (!is.null(indmiss_file) && file.exists(indmiss_file)) {
  message("Reading missing data from: ", indmiss_file)
  indmiss <- read.table(indmiss_file, header = TRUE, stringsAsFactors = FALSE)
  message("Missing data loaded: ", nrow(indmiss), " individuals")
  message("Columns: ", paste(names(indmiss), collapse = ", "))
}

# Debug: show columns in popdata
message("popdata columns: ", paste(names(popdata), collapse = ", "))
message("Requested color_by column: ", color_by_name)

# Generate PCA plot
plt_pca <- plot_pca(
  individuals = individuals,
  eigenvecs = eigenvecs,
  eigenvals = eigenvals,
  popdata = popdata,
  indmiss = indmiss,
  color_by_name = color_by_name,
  pc1 = pc1,
  pc2 = pc2
)

# Debug: check plot object
message("Generated PCA plot object: ", class(plt_pca))

# Save PDF
dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)
ggplot2::ggsave(output_pdf, plt_pca, width = 6, height = 5, dpi = 300)
message("PDF saved to ", output_pdf)

# Save RDS (ggplot2 object)
dir.create(dirname(output_rds), recursive = TRUE, showWarnings = FALSE)
saveRDS(plt_pca, output_rds)
message("RDS saved to ", output_rds)