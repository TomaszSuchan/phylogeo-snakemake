#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
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

summary <- read.table(
  snakemake@input[["summary"]],
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE,
  quote = ""
)
summary$migration_edges <- as.integer(summary$migration_edges)
summary$final_log_likelihood <- vapply(summary$llik, extract_final_llik, numeric(1))
summary <- summary[order(summary$migration_edges), ]

plot_obj <- ggplot2::ggplot(summary, ggplot2::aes(migration_edges, final_log_likelihood)) +
  ggplot2::geom_line(linewidth = 0.5) +
  ggplot2::geom_point(size = 2) +
  ggplot2::scale_x_continuous(breaks = summary$migration_edges) +
  ggplot2::labs(
    title = "OrientAGraph / TreeMix-compatible model likelihoods",
    x = "Migration edges (m)",
    y = "Final log likelihood"
  ) +
  ggplot2::theme_minimal(base_size = 11)

dir.create(dirname(snakemake@output[["pdf"]]), recursive = TRUE, showWarnings = FALSE)
ggsave_pdf(
  snakemake@output[["pdf"]],
  plot_obj,
  width = as.numeric(snakemake@params[["width"]]),
  height = as.numeric(snakemake@params[["height"]]),
  dpi = as.numeric(snakemake@params[["dpi"]])
)
saveRDS(list(plot = plot_obj, summary = summary), snakemake@output[["rds"]])

cat("Saved likelihood plot:", snakemake@output[["pdf"]], "\n")
