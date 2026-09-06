#!/usr/bin/env Rscript
# Arrange per-K ancestry barplots (K >= 2) side by side in one ROW, lowest K on
# the left and highest K on the right. Each panel is drawn as a *vertical*
# barplot (individuals run top-to-bottom, ancestry proportion 0..1 left-to-right)
# by flipping the shared horizontal per-K barplot RDS objects. Site labels are
# shown once in a shared strip on the far left; a single cluster colour key is
# placed on the far right.

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
  library(gtable)
})

pdf(NULL)

script_dir <- tryCatch(
  dirname(normalizePath(snakemake@script)),
  error = function(e) "workflow/scripts"
)
source(file.path(script_dir, "plot_ggsave_utils.R"))

log_file <- snakemake@log[[1]]
dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
# NOTE: only the *output* stream is redirected to the log connection. Sinking
# the *message* stream to a file connection deadlocks with grid::grid.grabExpr()
# under ggplot2 >= 4.0, so progress is logged via cat() (output stream) instead.
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output")
on.exit({
  while (sink.number(type = "output") > 0) sink(type = "output")
  close(log_con)
}, add = TRUE)

log_msg <- function(...) cat(sprintf(...), "\n", sep = "")

params <- snakemake@params
if (length(params) == 1L && is.list(params[[1L]]) && length(names(params)) == 0L) {
  params <- params[[1L]]
}

output_pdf <- snakemake@output[["pdf"]]
output_rds <- snakemake@output[["rds"]]
method_label <- params[["method_label"]]
label_width_in <- as.numeric(params[["label_width"]])
panel_gap_pt <- as.numeric(params[["panel_gap"]])
legend_pad_in <- as.numeric(params[["legend_pad"]])
if (!is.finite(legend_pad_in) || legend_pad_in < 0) {
  legend_pad_in <- 0.6
}
is_pos_num <- function(x) length(x) == 1L && is.finite(x) && x > 0
col_width_in <- suppressWarnings(as.numeric(params[["col_width"]]))
label_size <- suppressWarnings(as.numeric(params[["label_size"]]))
if (!is_pos_num(label_size)) {
  label_size <- 2.5
}

klabel_row_in <- 0.28

# ---------------------------------------------------------------------------
# Layer helpers: per-K barplot RDS objects share a common structure --
#   GeomBar      : the stacked ancestry bars
#   GeomSegment  : site dividers (span the full proportion axis) AND short site
#                  tick marks drawn just below the baseline (y <= 0)
#   GeomLabel    : site name labels drawn below the baseline (y < 0)
# In the horizontal source plot, individuals map to x and proportion to y.
# ---------------------------------------------------------------------------

is_tick_segment <- function(d) {
  is.data.frame(d) && nrow(d) > 0L &&
    all(c("y", "yend") %in% names(d)) &&
    is.finite(max(c(d$y, d$yend), na.rm = TRUE)) &&
    max(c(d$y, d$yend), na.rm = TRUE) <= 0
}

# Vertical bar panel: keep bars + dividers, drop site labels/ticks, flip to
# vertical and strip all axes/legend so panels tile cleanly.
build_bar_panel <- function(p, x_range) {
  p <- unserialize(serialize(p, NULL))
  for (i in seq_along(p@layers)) {
    d <- p@layers[[i]]$data
    if (!is.data.frame(d)) next
    geom <- class(p@layers[[i]]$geom)[1]
    if (geom == "GeomLabel" || (geom == "GeomSegment" && is_tick_segment(d))) {
      p@layers[[i]]$data <- d[0, , drop = FALSE]
    }
  }
  suppressMessages(
    p +
      coord_flip(xlim = x_range, ylim = c(0, 1), clip = "off") +
      labs(x = NULL, y = NULL, title = NULL) +
      theme(
        # Blank the specific .x/.y children: the source plot sets
        # axis.text.y explicitly, which would otherwise override a parent
        # axis.text = element_blank() and print per-individual labels after
        # the flip.
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.length = unit(0, "pt"),
        axis.line = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none",
        plot.margin = margin(0, panel_gap_pt / 2, 0, panel_gap_pt / 2, unit = "pt")
      )
  )
}

# Shared site-label strip: keep only site labels + tick marks, anchor labels
# just left of the baseline (right-justified) and flip to vertical.
build_label_strip <- function(p, x_range) {
  p <- unserialize(serialize(p, NULL))
  for (i in seq_along(p@layers)) {
    d <- p@layers[[i]]$data
    if (!is.data.frame(d)) next
    geom <- class(p@layers[[i]]$geom)[1]
    keep <- geom == "GeomLabel" || (geom == "GeomSegment" && is_tick_segment(d))
    if (!keep) {
      p@layers[[i]]$data <- d[0, , drop = FALSE]
    } else if (geom == "GeomLabel") {
      d$y <- -0.03
      p@layers[[i]]$data <- d
      p@layers[[i]]$aes_params$angle <- 0
      p@layers[[i]]$aes_params$hjust <- 1
      p@layers[[i]]$aes_params$size <- label_size
      p@layers[[i]]$aes_params$label.size <- 0
      p@layers[[i]]$aes_params$fill <- NA
    }
  }
  label_min <- -1.0
  suppressMessages(
    p +
      coord_flip(xlim = x_range, ylim = c(label_min, 0.02), clip = "off") +
      labs(x = NULL, y = NULL, title = NULL) +
      theme_void() +
      theme(
        legend.position = "none",
        plot.margin = margin(0, 2, 0, 2, unit = "pt")
      )
  )
}

# Pull the cluster colour key out of the richest (highest-K) panel.
extract_legend <- function(p) {
  g <- ggplotGrob(
    p + theme(legend.position = "right", legend.direction = "vertical")
  )
  idx <- which(g$layout$name == "guide-box-right")
  if (length(idx) == 0L) {
    idx <- grep("guide-box", g$layout$name)
  }
  if (length(idx) == 0L) {
    return(nullGrob())
  }
  g$grobs[[idx[[1]]]]
}

# ---------------------------------------------------------------------------
# Load per-K panels (sorted ascending: lowest K first -> leftmost).
# ---------------------------------------------------------------------------

barplot_input_names <- grep("^barplot_k", names(snakemake@input), value = TRUE)
if (length(barplot_input_names) == 0L) {
  stop("No per-K barplot RDS inputs found (expected names like barplot_k2, barplot_k3, ...).")
}
extract_k <- function(name) as.integer(sub("^barplot_k", "", name))
ord <- order(vapply(barplot_input_names, extract_k, integer(1)))
barplot_input_names <- barplot_input_names[ord]
barplot_ks <- vapply(barplot_input_names, extract_k, integer(1))

log_msg(
  "Building vertical %s barplot facet from K = %s (left -> right)",
  method_label,
  paste(barplot_ks, collapse = ", ")
)

barplot_paths <- snakemake@input[barplot_input_names]
panels_raw <- lapply(barplot_paths, function(path) {
  p <- readRDS(path)
  if (!inherits(p, "ggplot")) {
    stop("Expected ggplot in ", path)
  }
  p
})
n_panels <- length(panels_raw)

# Panel geometry (identical across K for a given project).
panel_width <- as.numeric(attr(panels_raw[[1L]], "panel_width", exact = TRUE))
bar_height_in <- suppressWarnings(
  as.numeric(attr(panels_raw[[1L]], "panel_bar_height_in", exact = TRUE))
)
if (!is.finite(panel_width) || panel_width <= 0) {
  stop("Barplot RDS missing panel_width; re-run barplot_* rules.")
}
# The individual axis (previously the wide dimension) becomes the panel height.
indiv_extent_in <- panel_width
# Each K column width defaults to the source bar thickness (kept tunable).
if (!is_pos_num(col_width_in)) {
  col_width_in <- if (is_pos_num(bar_height_in)) bar_height_in else 1.0
}

# Shared individual-axis range so every panel + label strip line up.
built1 <- ggplot_build(panels_raw[[1L]])
x_range <- built1$layout$panel_params[[1L]]$x.range

bar_panels <- lapply(panels_raw, build_bar_panel, x_range = x_range)
label_plot <- build_label_strip(panels_raw[[1L]], x_range = x_range)
legend_grob <- extract_legend(panels_raw[[n_panels]])

bar_grobs <- lapply(bar_panels, function(p) {
  grid::grid.grabExpr(
    suppressMessages(print(p)),
    width = col_width_in,
    height = indiv_extent_in,
    wrap = FALSE,
    wrap.grobs = FALSE
  )
})
label_grob <- grid::grid.grabExpr(
  suppressMessages(print(label_plot)),
  width = label_width_in,
  height = indiv_extent_in,
  wrap = FALSE,
  wrap.grobs = FALSE
)

# ---------------------------------------------------------------------------
# Assemble: 2 rows (K labels, panels) x (1 label col + n K cols + 1 legend col).
# ---------------------------------------------------------------------------

col_widths <- c(label_width_in, rep(col_width_in, n_panels), legend_pad_in)
row_heights <- c(klabel_row_in, indiv_extent_in)

gt <- gtable(
  widths = unit(col_widths, "in"),
  heights = unit(row_heights, "in")
)

gt <- gtable_add_grob(
  gt, label_grob, t = 2, l = 1, name = "site_labels", clip = "off"
)

for (i in seq_len(n_panels)) {
  col <- i + 1L
  gt <- gtable_add_grob(
    gt,
    textGrob(paste0("K = ", barplot_ks[[i]]), gp = gpar(fontsize = 9)),
    t = 1, l = col, name = paste0("k_", i)
  )
  gt <- gtable_add_grob(
    gt, bar_grobs[[i]], t = 2, l = col, name = paste0("bar_", i), clip = "off"
  )
}

gt <- gtable_add_grob(
  gt, legend_grob, t = 1, b = 2, l = n_panels + 2L, name = "legend", clip = "off"
)

page_width <- sum(col_widths)
page_height <- sum(row_heights)

log_msg(
  "Saving vertical barplot facet PDF: %s (%.2f x %.2f in)",
  output_pdf, page_width, page_height
)
ggsave_pdf(output_pdf, gt, width = page_width, height = page_height)

log_msg("Saving vertical barplot facet RDS: %s", output_rds)
saveRDS(gt, file = output_rds)

log_msg("Vertical barplot facet plot completed successfully.")
