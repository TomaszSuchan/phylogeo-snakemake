# Load libraries
library(RColorBrewer)
library(ggplot2)

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

# Simple PCA plotting function
plot_pca <- function(individuals, eigenvecs, eigenvals, popdata, color_by = 2, pc1 = 1, pc2 = 2) {

  # Create data frame with PC scores
  pca_df <- data.frame(
    Ind = individuals,
    PC1 = eigenvecs[, pc1],
    PC2 = eigenvecs[, pc2]
  )

  # Merge with population data (use first column of popdata for merging)
  ind_col <- colnames(popdata)[1]
  plot_df <- merge(pca_df, popdata, by.x = "Ind", by.y = ind_col)

  # Calculate variance explained
  var_exp <- eigenvals / sum(eigenvals) * 100

  # Create plot
  color_col <- colnames(popdata)[color_by]

  p <- ggplot(plot_df, aes(x = PC1, y = PC2, color = .data[[color_col]])) +
    geom_point(size = 3, alpha = 0.7) +
    labs(
      x = sprintf("PC%d (%.1f%%)", pc1, var_exp[pc1]),
      y = sprintf("PC%d (%.1f%%)", pc2, var_exp[pc2]),
      color = color_col
    ) +
    theme_bw()

  return(p)
}

message("plot_pca function loaded")

# Read Snakemake inputs
eigvecs_file <- snakemake@input[["eigvecs"]]
eigvals_file <- snakemake@input[["eigvals"]]
popdata_file <- snakemake@input[["indpopdata"]]
output_pdf   <- snakemake@output[["pdf"]]
output_rds   <- snakemake@output[["rds"]]

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

# Debug: show columns in popdata
message("popdata columns: ", paste(names(popdata), collapse = ", "))

# Convert color_by from column name to index
if (!color_by_name %in% names(popdata)) {
  stop(paste("Column", color_by_name, "not found in popdata"))
}
color_by_index <- which(names(popdata) == color_by_name)
message("color_by_index = ", color_by_index)

# Generate PCA plot
plt_pca <- plot_pca(
  individuals = individuals,
  eigenvecs = eigenvecs,
  eigenvals = eigenvals,
  popdata = popdata,
  color_by = color_by_index,
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