# Script to plot BIC values from per-K DAPC log files
# Parses "BIC value:" entries from each dapc.K*.log.txt file.
# Uses plot_choose_k_utils.R so DAPC criterion plots match Evanno choose-K styling.

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

script_dir <- tryCatch(
  dirname(normalizePath(snakemake@script)),
  error = function(e) "workflow/scripts"
)
source(file.path(script_dir, "plot_choose_k_utils.R"))
plot_dims <- read_choose_k_plot_dims(snakemake@params)

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
log_files <- as.character(snakemake@input)
output_plot <- snakemake@output[["plot"]]
output_plot_rds <- snakemake@output[["plot_rds"]]

cat("=== DAPC BIC Plot from Log ===\n")
cat("Input log files:\n")
for (lf in log_files) {
  cat("  ", lf, "\n")
}
cat("Output plot:", output_plot, "\n\n")

extract_k <- function(lines, path) {
  header_line <- grep("^=== DAPC Analysis \\(K\\s*=\\s*[0-9]+\\s*\\) ===$", lines, value = TRUE)
  if (length(header_line) > 0) {
    k_match <- regmatches(header_line[1], regexpr("[0-9]+", header_line[1]))
    if (length(k_match) > 0) return(as.integer(k_match))
  }
  params_line <- grep("^\\s*K:\\s*[0-9]+\\s*$", lines, value = TRUE)
  if (length(params_line) > 0) {
    k_match <- regmatches(params_line[1], regexpr("[0-9]+", params_line[1]))
    if (length(k_match) > 0) return(as.integer(k_match))
  }
  basename_k <- regmatches(basename(path), regexpr("K[0-9]+", basename(path)))
  if (length(basename_k) > 0) {
    return(as.integer(sub("K", "", basename_k)))
  }
  return(NA_integer_)
}

extract_bic <- function(lines) {
  bic_line <- grep("^\\s*BIC value:\\s*[-+]?[0-9]*\\.?[0-9]+\\s*$", lines, value = TRUE)
  if (length(bic_line) == 0) return(NA_real_)
  bic_match <- regmatches(bic_line[1], regexpr("[-+]?[0-9]*\\.?[0-9]+", bic_line[1]))
  as.numeric(bic_match)
}

parsed <- lapply(log_files, function(path) {
  if (!file.exists(path)) {
    return(data.frame(K = NA_integer_, BIC = NA_real_, source = path, stringsAsFactors = FALSE))
  }
  lines <- readLines(path, warn = FALSE)
  data.frame(
    K = extract_k(lines, path),
    BIC = extract_bic(lines),
    source = path,
    stringsAsFactors = FALSE
  )
})

bic_data <- do.call(rbind, parsed)
bic_data <- bic_data[is.finite(bic_data$BIC) & !is.na(bic_data$K), c("K", "BIC", "source")]
bic_data <- bic_data[bic_data$K >= 2, ]

if (nrow(bic_data) == 0) {
  stop("Could not parse any valid K/BIC pairs from DAPC log files")
}

if (any(duplicated(bic_data$K))) {
  cat("WARNING: Duplicate K values found in log files; keeping first occurrence per K.\n")
  bic_data <- bic_data[!duplicated(bic_data$K), ]
}

# Sort by K
bic_data <- bic_data[order(bic_data$K), ]

cat("Parsed BIC values:\n")
print(bic_data[, c("K", "BIC")])

# Determine optimal K (minimum BIC)
optimal_k <- bic_data$K[which.min(bic_data$BIC)]
cat("\nOptimal K:", optimal_k, "(Minimum BIC)\n")

p <- plot_choose_k_line(
  data = bic_data,
  x = "K",
  y = "BIC",
  ylab = "BIC"
)

# Save plot
dir.create(dirname(output_plot), recursive = TRUE, showWarnings = FALSE)
choose_k_ggsave(
  filename = output_plot,
  plot = p,
  width = plot_dims$width,
  height = plot_dims$height,
  dpi = plot_dims$dpi
)
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

