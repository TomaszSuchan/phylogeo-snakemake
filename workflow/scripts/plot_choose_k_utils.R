ggsave_utils <- tryCatch(
  file.path(dirname(normalizePath(snakemake@script)), "plot_ggsave_utils.R"),
  error = function(e) "workflow/scripts/plot_ggsave_utils.R"
)
if (file.exists(ggsave_utils)) {
  source(ggsave_utils)
} else {
  source("workflow/scripts/plot_ggsave_utils.R")
}

# Shared ggplot helpers for choose-K / model-comparison score plots.
#
# Why a separate file?
# Several analyses (STRUCTURE Evanno, ADMIXTURE CV, tess3r CV, conStruct MAP
# log-posterior, DAPC BIC) need the same publication-ready line-plot styling.
# plot_evanno.R is the reference implementation; this module keeps theme, geoms,
# and ggsave dimensions consistent without duplicating code in each script.
#
# Dimensions come from config parameters.choose_k_plot (passed via Snakemake
# rule params), not from map_background, because score plots are not maps.

DEFAULT_CHOOSE_K_PLOT_WIDTH <- 10
DEFAULT_CHOOSE_K_PLOT_HEIGHT <- 5
DEFAULT_CHOOSE_K_PLOT_DPI <- 300

read_choose_k_plot_dims <- function(params) {
  get_param <- function(name, default) {
    val <- params[[name]]
    if (is.null(val) || length(val) == 0 || any(is.na(val))) {
      return(default)
    }
    as.numeric(val)[1]
  }
  list(
    width = get_param("width", DEFAULT_CHOOSE_K_PLOT_WIDTH),
    height = get_param("height", DEFAULT_CHOOSE_K_PLOT_HEIGHT),
    dpi = get_param("dpi", DEFAULT_CHOOSE_K_PLOT_DPI)
  )
}

choose_k_plot_theme <- function() {
  ggplot2::theme_bw() +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(face = "italic"),
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank()
    )
}

choose_k_ggsave <- function(
  filename,
  plot,
  width = DEFAULT_CHOOSE_K_PLOT_WIDTH,
  height = DEFAULT_CHOOSE_K_PLOT_HEIGHT,
  dpi = DEFAULT_CHOOSE_K_PLOT_DPI,
  ...
) {
  ggsave_pdf(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    ...
  )
}

plot_choose_k_line <- function(data, x, y, ymin = NULL, ymax = NULL, ylab = NULL, xlab = "K") {
  mapping <- ggplot2::aes(x = .data[[x]], y = .data[[y]])
  if (!is.null(ymin) && !is.null(ymax) && ymin %in% names(data) && ymax %in% names(data)) {
    mapping <- ggplot2::aes(
      x = .data[[x]],
      y = .data[[y]],
      ymin = .data[[ymin]],
      ymax = .data[[ymax]]
    )
  }

  p <- ggplot2::ggplot(data, mapping) +
    ggplot2::geom_line() +
    ggplot2::geom_point()

  if (!is.null(ymin) && !is.null(ymax) && ymin %in% names(data) && ymax %in% names(data)) {
    p <- p + ggplot2::geom_errorbar(width = 0.3, linewidth = 0.3)
  }

  p +
    ggplot2::scale_x_continuous(breaks = sort(unique(data[[x]]))) +
    ggplot2::xlab(xlab) +
    {if (!is.null(ylab)) ggplot2::ylab(ylab) else ggplot2::waiver()} +
    choose_k_plot_theme()
}
