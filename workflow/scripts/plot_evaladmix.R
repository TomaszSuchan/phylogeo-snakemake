#!/usr/bin/env Rscript
# Plot evalAdmix correlation matrix using official visFuns.R from evalAdmix

# Redirect all output to log file
log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

# Snakemake inputs/outputs
corres_file <- snakemake@input[["corres"]]
indpopdata_file <- snakemake@input[["indpopdata"]]
qfile <- snakemake@input[["qfile"]]
visfuns_file <- snakemake@input[["visfuns"]]
output_pdf <- snakemake@output[["pdf"]]
output_rds <- snakemake@output[["rds"]]
method_name <- snakemake@params[["method"]]
k_value <- snakemake@params[["k"]]

message("\n=== LOADING VISFUNS.R ===\n")
source(visfuns_file)

message("\n=== READING INPUT FILES ===\n")

# Read correlation matrix (r)
message("Reading correlation matrix from: ", corres_file)
r <- as.matrix(read.table(corres_file, header = FALSE, sep = " "))
message(sprintf("Correlation matrix dimensions: %d x %d", nrow(r), ncol(r)))

# Read admixture proportions (q) - optional but recommended for ordering
message("Reading admixture proportions from: ", qfile)
q <- as.matrix(read.table(qfile, header = FALSE, sep = " "))
message(sprintf("Q matrix dimensions: %d x %d", nrow(q), ncol(q)))

# Read indpopdata to get population assignments in filtered VCF order
message("Reading indpopdata from: ", indpopdata_file)
indpopdata <- read.table(indpopdata_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
indpopdata <- indpopdata[!duplicated(indpopdata$Ind), ]
if (nrow(indpopdata) > nrow(q)) {
  indpopdata <- indpopdata[1:nrow(q), ]
}
pop <- as.vector(indpopdata$Site)
message(sprintf("indpopdata contains %d individuals", length(pop)))

# Check dimensions match
if (nrow(r) != length(pop)) {
  stop(sprintf("Dimension mismatch: correlation matrix has %d rows but indpopdata has %d individuals",
               nrow(r), length(pop)))
}

if (nrow(q) != length(pop)) {
  stop(sprintf("Dimension mismatch: Q matrix has %d rows but indpopdata has %d individuals",
               nrow(q), length(pop)))
}

message("\n=== ORDERING INDIVIDUALS ===\n")
# Order individuals by population and admixture proportions
ord <- orderInds(pop = pop, q = q)
message("Individuals ordered for plotting")

message("\n=== CREATING PLOT ===\n")
# Create output directory
dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)

# Open PDF device
pdf(output_pdf, width = 10, height = 9)

# Plot correlation of residuals using visFuns.R function
plotCorRes(
  cor_mat = r,
  pop = pop,
  ord = ord,
  title = sprintf("Admixture evaluation as correlation of residuals\n%s, K=%s", method_name, k_value),
  max_z = 0.25,
  min_z = -0.25,
  cex.main = 1.5,
  cex.lab = 1.5,
  cex.legend = 1.5,
  color_palette = c("#001260", "#EAEDE9", "#601200"),
  pop_labels = c(TRUE, TRUE),
  plot_legend = TRUE,
  adjlab = 0.1,
  rotatelabpop = 0,
  rotatelabsuperpop = 0,
  lineswidth = 1,
  lineswidthsuperpop = 2,
  adjlabsuperpop = 0.16,
  cex.lab.2 = 1.5
)

dev.off()

# Save plot object to RDS (if plotCorRes returns a plot object, otherwise save the data)
plot_data <- list(
  cor_mat = r,
  pop = pop,
  ord = ord,
  q = q,
  method = method_name,
  k = k_value
)
saveRDS(plot_data, output_rds)

message("\n=== COMPLETED SUCCESSFULLY ===\n")
message(sprintf("Output PDF: %s", output_pdf))
message(sprintf("Output RDS: %s", output_rds))

# Close log file sinks
sink(type = "message")
sink(type = "output")
close(log_file)
