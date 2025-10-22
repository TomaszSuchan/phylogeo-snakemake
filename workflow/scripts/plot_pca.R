# Load R packages
library(RColorBrewer)
library(ggplot2)
source("external/phylogeographeR/R/plot_PCA.R")

# Snakemake inputs and params
eigvecs_file <- snakemake@input[["eigvecs"]]
eigvals_file <- snakemake@input[["eigvals"]]
popdata_file <- snakemake@input[["popdata"]]
output_file  <- snakemake@output[[1]]

pc1 <- as.numeric(snakemake@params[["pc1"]])
pc2 <- as.numeric(snakemake@params[["pc2"]])
color_by_name <- as.character(snakemake@params[["color_by"]])

# Read inputs
eigenvecs2 <- read.table(eigvecs_file)
individuals <- eigenvecs2$V1
eigenvecs <- eigenvecs2[, -c(1,2)]
eigenvals <- scan(eigvals_file)
popdata <- read.table(popdata_file, header = TRUE)

# Convert color_by from column name to index (as required by plot_pca)
if (!color_by_name %in% names(popdata)) {
  stop(paste("Column", color_by_name, "not found in popdata"))
}
color_by_index <- which(names(popdata) == color_by_name)

# Generate PCA plot using your package function
plt_pca <- plot_pca(
  individuals = individuals,
  eigenvecs = eigenvecs,
  eigenvals = eigenvals,
  popdata = popdata,
  color_by = color_by_index,
  pc1 = pc1,
  pc2 = pc2
)

# Save the plot
ggplot2::ggsave(output_file, plt_pca, width = 6, height = 5, dpi = 300)