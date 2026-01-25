# Script to plot BIC values from dapc_bic_plot.log.txt
# Parses the native BIC values extracted from find.clusters and creates a plot
# matching the aesthetics of the original criterion plot

# Prevent creation of Rplots.pdf - set before loading libraries
Sys.setenv(R_INSTALL_PKG = "")
options(device = function(file = NULL, ...) {
  if (is.null(file)) {
    pdf(NULL)
  } else {
    pdf(file, ...)
  }
})

library(ggplot2)

# Explicitly prevent Rplots.pdf creation
pdf(NULL)
# Also try to delete it if it exists
if (file.exists("Rplots.pdf")) {
  tryCatch({
    file.remove("Rplots.pdf")
  }, error = function(e) {
    # Silently fail
  })
}

# Snakemake inputs/outputs
log_file <- snakemake@input[["log_file"]]
output_plot <- snakemake@output[["plot"]]
output_plot_rds <- snakemake@output[["plot_rds"]]

cat("=== DAPC BIC Plot from Log ===\n")
cat("Log file:", log_file, "\n")
cat("Output plot:", output_plot, "\n\n")

# Read the log file
log_content <- readLines(log_file)

# Find the line with "BIC values extracted from Kstat:"
bic_header_idx <- grep("BIC values extracted from Kstat:", log_content)
if (length(bic_header_idx) == 0) {
  stop("Could not find 'BIC values extracted from Kstat:' in log file")
}

# Get the next few lines (usually 2-3 lines contain all K and BIC values)
bic_lines <- log_content[(bic_header_idx + 1):min(bic_header_idx + 5, length(log_content))]
bic_lines <- bic_lines[bic_lines != ""]  # Remove empty lines

# Parse K values and BIC values
k_values <- integer()
bic_values <- numeric()

for (line in bic_lines) {
  line <- trimws(line)
  if (line == "") next
  
  # Check if line contains K= pattern
  if (grepl("K=\\d+", line)) {
    # Extract all K values from this line
    k_matches <- regmatches(line, gregexpr("K=(\\d+)", line))[[1]]
    k_vals <- as.integer(sub("K=(\\d+)", "\\1", k_matches))
    k_values <- c(k_values, k_vals)
  } else if (grepl("^\\d+\\.\\d+", line)) {
    # This line contains BIC values (numbers)
    # Extract all numbers
    bic_vals <- as.numeric(strsplit(line, "\\s+")[[1]])
    bic_values <- c(bic_values, bic_vals)
  }
}

# Create data frame - match K values with BIC values
# The format is: K values on one line, BIC values on the next line(s)
# They should be in the same order

if (length(k_values) == 0 || length(bic_values) == 0) {
  stop("Could not extract K values or BIC values from log file")
}

if (length(k_values) != length(bic_values)) {
  cat("WARNING: Mismatch between K values and BIC values\n")
  cat("K values (", length(k_values), "):", paste(k_values, collapse = ", "), "\n")
  cat("BIC values (", length(bic_values), "):", paste(bic_values, collapse = ", "), "\n")
  cat("BIC lines:\n")
  for (line in bic_lines) {
    cat("  ", line, "\n")
  }
  # Try to match what we can
  min_len <- min(length(k_values), length(bic_values))
  k_values <- k_values[1:min_len]
  bic_values <- bic_values[1:min_len]
  cat("Using first", min_len, "values\n")
}

bic_data <- data.frame(K = k_values, BIC = bic_values)

# Sort by K
bic_data <- bic_data[order(bic_data$K), ]

cat("Parsed BIC values:\n")
print(bic_data)

# Determine optimal K (minimum BIC)
optimal_k <- bic_data$K[which.min(bic_data$BIC)]
cat("\nOptimal K:", optimal_k, "(Minimum BIC)\n")

# Create plot with same aesthetics as criterion plot
p <- ggplot(bic_data, aes(x = K, y = BIC)) +
  geom_point() +
  geom_line() +
  labs(
    x = "K",
    y = "BIC",
    title = "",
    subtitle = ""
  ) +
  scale_x_continuous(breaks = bic_data$K) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(face = "italic")
  )

# Save plot
dir.create(dirname(output_plot), recursive = TRUE, showWarnings = FALSE)
ggsave(output_plot, p, width = 4, height = 3, dpi = 300)
cat("Plot saved to:", output_plot, "\n")

# Save RDS
saveRDS(p, output_plot_rds)
cat("Plot RDS saved to:", output_plot_rds, "\n")

cat("\n=== BIC Plot Complete ===\n")

# Final cleanup - remove Rplots.pdf if it was created
if (file.exists("Rplots.pdf")) {
  tryCatch({
    file.remove("Rplots.pdf")
  }, error = function(e) {
    # Silently fail
  })
}

