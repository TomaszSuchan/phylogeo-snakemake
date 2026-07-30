#!/usr/bin/env Rscript

# OptM computes Evanno / linear / SiZer stats; the PDF uses the shared
# STRUCTURE Evanno-style ggplot helpers (plot_choose_k_utils.R).

suppressPackageStartupMessages({
  library(ggplot2)
  library(OptM)
})

pdf(NULL)

script_dir <- tryCatch(
  dirname(normalizePath(snakemake@script)),
  error = function(e) "workflow/scripts"
)
source(file.path(script_dir, "plot_choose_k_utils.R"))
plot_dims <- read_choose_k_plot_dims(snakemake@params)

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

# Build STRUCTURE Evanno-style facets from OptM Evanno table: L(m) ± range and Δm.
plot_optm_evanno_ggplot <- function(optm_df) {
  needed <- c("m", "mean(Lm)", "min(Lm)", "max(Lm)", "Deltam")
  missing <- setdiff(needed, names(optm_df))
  if (length(missing) > 0) {
    stop("OptM Evanno table missing columns: ", paste(missing, collapse = ", "))
  }

  lm_df <- data.frame(
    m = optm_df[["m"]],
    Parameter = "Lm",
    Value = optm_df[["mean(Lm)"]],
    Min = optm_df[["min(Lm)"]],
    Max = optm_df[["max(Lm)"]],
    stringsAsFactors = FALSE
  )
  delta_df <- data.frame(
    m = optm_df[["m"]],
    Parameter = "Deltam",
    Value = optm_df[["Deltam"]],
    Min = optm_df[["Deltam"]],
    Max = optm_df[["Deltam"]],
    stringsAsFactors = FALSE
  )
  plot_df <- rbind(lm_df, delta_df)
  plot_df$Parameter <- factor(plot_df$Parameter, levels = c("Lm", "Deltam"))

  facet_labels <- c(
    Lm = "italic(L)(italic(m))~\"\u00B1 SD\"",
    Deltam = "Delta*italic(m)"
  )

  ggplot2::ggplot(plot_df, ggplot2::aes(x = m, y = Value, ymin = Min, ymax = Max)) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::geom_errorbar(width = 0.3, linewidth = 0.3) +
    ggplot2::scale_x_continuous(breaks = sort(unique(plot_df$m))) +
    ggplot2::facet_wrap(
      ~Parameter,
      scales = "free",
      labeller = ggplot2::as_labeller(facet_labels, default = ggplot2::label_parsed)
    ) +
    ggplot2::xlab("m") +
    choose_k_plot_theme()
}

folder <- snakemake@params[["folder"]]
method <- as.character(snakemake@params[["method"]])

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
  if (identical(method, "Evanno") && is.data.frame(optm_result)) {
    plot_obj <- plot_optm_evanno_ggplot(optm_result)
    choose_k_ggsave(
      filename = snakemake@output[["pdf"]],
      plot = plot_obj,
      width = plot_dims$width,
      height = plot_dims$height,
      dpi = plot_dims$dpi
    )
  } else {
    # linear / SiZer: keep OptM's native plot
    OptM::plot_optM(
      optm_result,
      method = method,
      plot = FALSE,
      pdf = snakemake@output[["pdf"]]
    )
    plot_obj <- NULL
  }
  saveRDS(
    list(
      plot = plot_obj,
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
  plot_obj <- plot_choose_k_line(
    data = data.frame(
      m = summary$migration_edges,
      Value = summary$final_log_likelihood
    ),
    x = "m",
    y = "Value",
    ylab = "Final log likelihood",
    xlab = "m"
  )
  choose_k_ggsave(
    filename = snakemake@output[["pdf"]],
    plot = plot_obj,
    width = plot_dims$width,
    height = plot_dims$height,
    dpi = plot_dims$dpi
  )
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
