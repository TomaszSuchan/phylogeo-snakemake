ggsave_utils <- tryCatch(
  file.path(dirname(normalizePath(snakemake@script)), "plot_ggsave_utils.R"),
  error = function(e) "workflow/scripts/plot_ggsave_utils.R"
)
if (file.exists(ggsave_utils)) {
  source(ggsave_utils)
} else {
  source("workflow/scripts/plot_ggsave_utils.R")
}

# Shared ggplot helpers for histogram plots (stats_vcf, ROH, genome_scan, etc.).

HISTOGRAM_FILL <- "steelblue"
HISTOGRAM_COLOR <- "black"
HISTOGRAM_ALPHA <- 0.7
HISTOGRAM_DEFAULT_BINS <- 50

histogram_plot_theme <- function() {
  ggplot2::theme_bw() +
    ggplot2::theme(
      axis.title = ggplot2::element_text(size = 12),
      axis.text = ggplot2::element_text(size = 10)
    )
}

geom_histogram_styled <- function(bins = HISTOGRAM_DEFAULT_BINS, ...) {
  ggplot2::geom_histogram(
    bins = bins,
    fill = HISTOGRAM_FILL,
    color = HISTOGRAM_COLOR,
    alpha = HISTOGRAM_ALPHA,
    ...
  )
}

geom_boxplot_styled <- function(fill = NULL, ...) {
  args <- list(
    alpha = HISTOGRAM_ALPHA,
    color = HISTOGRAM_COLOR,
    ...
  )
  if (!is.null(fill)) {
    args$fill <- fill
  }
  do.call(ggplot2::geom_boxplot, args)
}
