#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(OptM)
})

pdf(NULL)

log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")
on.exit({
  while (sink.number(type = "message") > 0) sink(type = "message")
  while (sink.number(type = "output") > 0) sink(type = "output")
  close(log_file)
}, add = TRUE)

extract_final_llik <- function(text) {
  matches <- regmatches(text, gregexpr("Exiting ln\\(likelihood\\) [-0-9.eE]+", text))[[1]]
  if (length(matches) == 0) {
    return(NA_real_)
  }
  as.numeric(sub(".* ", "", tail(matches, 1)))
}

extract_final_m <- function(text) {
  matches <- regmatches(text, gregexpr("Exiting ln\\(likelihood\\) [-0-9.eE]+ with [0-9]+ migration events", text))[[1]]
  if (length(matches) == 0) {
    return(NA_integer_)
  }
  as.integer(sub(".* with ([0-9]+) migration events", "\\1", tail(matches, 1)))
}

optm_run_summary <- function(folder) {
  files <- list.files(folder, pattern = "\\.llik$", full.names = TRUE)
  if (length(files) == 0) {
    return(data.frame())
  }
  rows <- lapply(files, function(path) {
    base <- basename(path)
    text <- paste(readLines(path, warn = FALSE), collapse = " ")
    data.frame(
      file = path,
      replicate = as.integer(sub(".*\\.r([0-9]+)\\.m[0-9]+\\.llik$", "\\1", base)),
      requested_m = as.integer(sub(".*\\.m([0-9]+)\\.llik$", "\\1", base)),
      achieved_m = extract_final_m(text),
      final_log_likelihood = extract_final_llik(text),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

format_range <- function(values) {
  values <- sort(unique(values[!is.na(values)]))
  if (length(values) == 0) {
    return("none")
  }
  paste(values, collapse = ",")
}

fallback_summary <- function(summary_path) {
  summary <- read.table(
    summary_path,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = ""
  )
  summary$migration_edges <- as.integer(summary$migration_edges)
  summary$final_log_likelihood <- vapply(summary$llik, extract_final_llik, numeric(1))
  summary[order(summary$migration_edges), ]
}

folder <- snakemake@params[["folder"]]
method <- as.character(snakemake@params[["method"]])
width <- as.numeric(snakemake@params[["width"]])
height <- as.numeric(snakemake@params[["height"]])
dpi <- as.numeric(snakemake@params[["dpi"]])

dir.create(dirname(snakemake@output[["pdf"]]), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(snakemake@output[["tsv"]]), recursive = TRUE, showWarnings = FALSE)

cat("Running OptM on folder:", folder, "\n")
cat("Method:", method, "\n")

optm_result <- tryCatch(
  OptM::optM(
    folder,
    orientagraph = TRUE,
    method = method,
    tsv = snakemake@output[["tsv"]]
  ),
  error = function(e) {
    cat("OptM failed; writing fallback likelihood summary.\n")
    cat("Reason:", conditionMessage(e), "\n")
    NULL
  }
)

if (!is.null(optm_result)) {
  grDevices::pdf(snakemake@output[["pdf"]], width = width, height = height)
  tryCatch(
    OptM::plot_optM(optm_result, method = method),
    error = function(e) {
      plot.new()
      title("OptM plot failed")
      text(0.5, 0.5, conditionMessage(e), cex = 0.8)
      cat("OptM plot failed:", conditionMessage(e), "\n")
    }
  )
  grDevices::dev.off()
  saveRDS(
    list(
      optm = optm_result,
      method = method,
      orientagraph = TRUE,
      fallback = FALSE
    ),
    snakemake@output[["rds"]]
  )
  cat("Saved OptM outputs.\n")
} else {
  summary <- fallback_summary(snakemake@input[["summary"]])
  optm_runs <- optm_run_summary(folder)
  n_replicates <- if (nrow(optm_runs) > 0) length(unique(optm_runs$replicate)) else NA_integer_
  requested_m <- if (nrow(optm_runs) > 0) optm_runs$requested_m else summary$migration_edges
  achieved_m <- if (nrow(optm_runs) > 0) optm_runs$achieved_m else summary$migration_edges
  requested_text <- format_range(requested_m)
  achieved_text <- format_range(achieved_m)
  replicate_text <- ifelse(is.na(n_replicates), "unknown", as.character(n_replicates))
  status_text <- paste0(
    "OptM could not estimate ", method,
    " statistics. Replicates present: ", replicate_text,
    "; requested m values: ", requested_text,
    "; achieved migration-event counts in .llik files: ", achieved_text,
    ". Evanno requires more than three usable migration-edge levels."
  )
  summary$optm_status <- paste(
    status_text,
    "If this is not a tiny test dataset, expand treemix.migration_edges and inspect",
    "whether higher-m runs actually add the requested number of migration events."
  )
  write.table(
    summary,
    file = snakemake@output[["tsv"]],
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  plot_obj <- ggplot2::ggplot(summary, ggplot2::aes(migration_edges, final_log_likelihood)) +
    ggplot2::geom_line(linewidth = 0.5) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_x_continuous(breaks = summary$migration_edges) +
    ggplot2::labs(
      title = "OptM fallback: final log likelihoods",
      subtitle = paste0("Delta-m unavailable; requested m=", requested_text, ", achieved m=", achieved_text),
      x = "Migration edges (m)",
      y = "Final log likelihood"
    ) +
    ggplot2::theme_minimal(base_size = 11)
  ggplot2::ggsave(snakemake@output[["pdf"]], plot_obj, width = width, height = height, dpi = dpi)
  saveRDS(
    list(
      plot = plot_obj,
      summary = summary,
      optm_runs = optm_runs,
      method = method,
      orientagraph = TRUE,
      fallback = TRUE
    ),
    snakemake@output[["rds"]]
  )
}

cat("OptM/fallback PDF:", snakemake@output[["pdf"]], "\n")
