#!/usr/bin/env Rscript
# Shared builder for multi-K ancestry barplot facet panels.
#
# build_barplot_facet_gtable() stacks per-K barplot ggplot objects (K >= 2) into
# a single-column gtable with side "K = n" labels and a shared site-label strip
# at the bottom. Both the column layout (plot_structure_barplot_facet.R) and the
# vertical layout (plot_structure_barplot_facet_vertical.R, which rotates the
# result 90 degrees) use this so the two facets stay visually consistent.
#
# `panels_raw` and `barplot_ks` must be supplied in the desired top-to-bottom
# order (ascending K for the column layout; descending K for the vertical layout
# so that, after a 90-degree clockwise rotation, the lowest K ends up on the
# left and the highest K on the right).

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
  library(gtable)
})

build_barplot_facet_gtable <- function(panels_raw, barplot_ks,
                                        flip_axis = FALSE,
                                        label_width_in,
                                        panel_gap_pt,
                                        legend_pad_in = 0.4) {
  flip_axis <- isTRUE(flip_axis)
  if (!is.finite(legend_pad_in) || legend_pad_in < 0) {
    legend_pad_in <- 0.4
  }
  ylab_width <- 0.35

  is_site_tick_segment <- function(layer_data) {
    nrow(layer_data) > 0L &&
      all(c("y", "yend") %in% names(layer_data)) &&
      is.finite(max(c(layer_data$y, layer_data$yend), na.rm = TRUE)) &&
      max(c(layer_data$y, layer_data$yend), na.rm = TRUE) <= 0
  }

  is_site_tick_segment_flipped <- function(layer_data) {
    nrow(layer_data) > 0L &&
      all(c("x", "xend") %in% names(layer_data)) &&
      is.finite(max(c(layer_data$x, layer_data$xend), na.rm = TRUE)) &&
      max(c(layer_data$x, layer_data$xend), na.rm = TRUE) <= 0
  }

  strip_for_facet <- function(p, compact_bottom = FALSE) {
    p <- unserialize(serialize(p, NULL))
    for (i in seq_along(p$layers)) {
      d <- p$layers[[i]]$data
      if (!is.data.frame(d)) {
        next
      }
      if (inherits(p$layers[[i]]$geom, "GeomLabel")) {
        p$layers[[i]]$data <- d[0, , drop = FALSE]
      } else if (inherits(p$layers[[i]]$geom, "GeomSegment")) {
        tick_segment <- if (flip_axis) {
          is_site_tick_segment_flipped(d)
        } else {
          is_site_tick_segment(d)
        }
        if (tick_segment) {
          p$layers[[i]]$data <- d[0, , drop = FALSE]
        }
      }
    }

    half_gap <- panel_gap_pt / 2
    bottom_margin <- if (compact_bottom) 0 else half_gap
    if (flip_axis) {
      p <- p + theme(
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.length.y = unit(0, "pt"),
        axis.title.x = element_blank()
      )
    } else {
      p <- p + theme(
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.length.x = unit(0, "pt"),
        axis.line.x = element_blank(),
        axis.title.y = element_blank()
      )
    }

    p + theme(
      plot.margin = margin(half_gap, 5, bottom_margin, 5, unit = "pt")
    )
  }

  build_site_label_strip <- function(p, bottom_margin_pt = 30) {
    p <- unserialize(serialize(p, NULL))
    built <- ggplot_build(p)
    if (flip_axis) {
      axis_min <- built$layout$panel_params[[1L]]$x.range[1L]
      label_lim <- c(axis_min - 0.18, -0.02)
    } else {
      axis_min <- built$layout$panel_params[[1L]]$y.range[1L]
      label_lim <- c(axis_min - 0.18, -0.02)
    }

    for (i in seq_along(p$layers)) {
      d <- p$layers[[i]]$data
      if (!is.data.frame(d)) {
        next
      }
      if (inherits(p$layers[[i]]$geom, "GeomBar")) {
        p$layers[[i]]$data <- d[0, , drop = FALSE]
      } else if (inherits(p$layers[[i]]$geom, "GeomSegment")) {
        tick_segment <- if (flip_axis) {
          is_site_tick_segment_flipped(d)
        } else {
          is_site_tick_segment(d)
        }
        if (!tick_segment) {
          p$layers[[i]]$data <- d[0, , drop = FALSE]
        }
      }
    }

    if (flip_axis) {
      p <- p + coord_cartesian(xlim = label_lim, clip = "off")
      p <- p + theme(
        axis.text.x = element_text(size = 2.2, colour = "transparent"),
        axis.ticks.x = element_line(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.length.y = unit(0, "pt"),
        axis.line.y = element_blank()
      )
    } else {
      p <- p + coord_cartesian(ylim = label_lim, clip = "off")
      p <- p + theme(
        axis.text.y = element_text(size = 5, colour = "transparent"),
        axis.ticks.y = element_line(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.length.x = unit(0, "pt"),
        axis.line.x = element_blank()
      )
    }

    p + theme(
      legend.position = "none",
      panel.grid = element_blank(),
      plot.margin = margin(22, 5, bottom_margin_pt, 5, unit = "pt")
    )
  }

  align_plot_grob <- function(grob, width_in, valign = c("center", "bottom"),
                              shift_down_in = 0) {
    valign <- match.arg(valign)
    just_y <- if (valign == "bottom") "bottom" else "center"
    y_pos <- if (valign == "bottom") unit(-shift_down_in, "in") else 0.5
    gTree(
      children = gList(grob),
      vp = viewport(
        x = unit(0, "in"),
        y = y_pos,
        width = unit(width_in, "in"),
        height = unit(1, "npc"),
        just = c("left", just_y),
        clip = "on"
      ),
      cl = "aligned_panel"
    )
  }

  has_site_labels <- function(p) {
    any(vapply(
      p$layers,
      function(layer) {
        is.data.frame(layer$data) &&
          inherits(layer$geom, "GeomLabel") &&
          nrow(layer$data) > 0L
      },
      logical(1)
    ))
  }

  panel_width <- attr(panels_raw[[1L]], "panel_width", exact = TRUE)
  panel_height <- attr(panels_raw[[1L]], "panel_height", exact = TRUE)
  if (is.null(panel_width) || is.null(panel_height)) {
    stop("Barplot RDS missing panel_width/panel_height; re-run barplot_* rules.")
  }
  panel_width <- as.numeric(panel_width)
  panel_height <- as.numeric(panel_height)

  bottom_margin_pt <- attr(panels_raw[[1L]], "panel_bottom_margin_pt", exact = TRUE)
  bar_height_in <- attr(panels_raw[[1L]], "panel_bar_height_in", exact = TRUE)
  site_label_height_in <- attr(panels_raw[[1L]], "panel_site_label_height_in", exact = TRUE)
  if (is.null(bottom_margin_pt) || !is.finite(as.numeric(bottom_margin_pt))) {
    bottom_margin_pt <- 30
  }
  if (is.null(bar_height_in) || !is.finite(as.numeric(bar_height_in))) {
    bar_height_in <- panel_height
  } else {
    bar_height_in <- as.numeric(bar_height_in)
  }
  if (is.null(site_label_height_in) || !is.finite(as.numeric(site_label_height_in))) {
    site_label_height_in <- bottom_margin_pt / 72
  } else {
    site_label_height_in <- as.numeric(site_label_height_in)
  }
  bottom_margin_pt <- as.numeric(bottom_margin_pt)
  label_bottom_pt <- bottom_margin_pt + 12

  n_panels <- length(panels_raw)
  show_site_label_row <- has_site_labels(panels_raw[[1L]])
  panels <- lapply(seq_len(n_panels), function(i) {
    strip_for_facet(
      panels_raw[[i]],
      compact_bottom = show_site_label_row && i == n_panels
    )
  })

  plot_width <- panel_width + legend_pad_in
  label_gap_in <- 1 / 72
  site_label_row_in <- label_bottom_pt / 72
  row_heights_in <- rep(bar_height_in, n_panels)
  if (show_site_label_row) {
    row_heights_in <- c(row_heights_in, label_gap_in, site_label_row_in)
  }
  label_row <- if (show_site_label_row) n_panels + 2L else NA_integer_

  axis_label <- if (flip_axis) {
    panels_raw[[1L]]$labels$x
  } else {
    panels_raw[[1L]]$labels$y
  }
  if (is.null(axis_label) || !nzchar(axis_label)) {
    axis_label <- NULL
  }

  bar_grobs <- lapply(panels, function(p) {
    grid::grid.grabExpr(
      print(p),
      width = panel_width,
      height = bar_height_in,
      wrap = FALSE,
      wrap.grobs = FALSE
    )
  })

  label_grob <- NULL
  if (show_site_label_row) {
    label_plot <- build_site_label_strip(panels_raw[[1L]], bottom_margin_pt = label_bottom_pt)
    label_grob <- grid::grid.grabExpr(
      print(label_plot),
      width = panel_width,
      height = site_label_row_in,
      wrap = FALSE,
      wrap.grobs = FALSE
    )
  }

  gt <- gtable(
    heights = unit(row_heights_in, "in"),
    widths = unit(c(ylab_width, label_width_in, plot_width), "in")
  )

  if (!is.null(axis_label)) {
    gt <- gtable_add_grob(
      gt,
      textGrob(axis_label, rot = 90, gp = gpar(fontsize = 8)),
      t = 1,
      b = n_panels,
      l = 1,
      name = "axis_label"
    )
  }

  for (i in seq_len(n_panels)) {
    gt <- gtable_add_grob(
      gt,
      textGrob(paste0("K = ", barplot_ks[[i]]), rot = 90, gp = gpar(fontsize = 9)),
      t = i,
      l = 2,
      name = paste0("k_", i)
    )
    gt <- gtable_add_grob(
      gt,
      align_plot_grob(bar_grobs[[i]], panel_width, shift_down_in = 0),
      t = i,
      l = 3,
      name = paste0("bar_", i),
      clip = if (show_site_label_row && i == n_panels) "on" else "off"
    )
  }

  if (show_site_label_row) {
    gt <- gtable_add_grob(
      gt,
      align_plot_grob(label_grob, panel_width, valign = "bottom", shift_down_in = 0.08),
      t = label_row,
      l = 3,
      name = "site_labels",
      clip = "on"
    )
  }

  page_width <- ylab_width + label_width_in + plot_width
  page_height <- sum(row_heights_in)

  list(gt = gt, width = page_width, height = page_height)
}

# Rotate a gtable/grob 90 degrees clockwise. Returns the rotated grob plus the
# swapped page dimensions to pass to ggsave().
rotate_grob_cw <- function(grob, width_in, height_in) {
  rotated <- grid::grobTree(
    grob,
    vp = grid::viewport(
      width = grid::unit(width_in, "in"),
      height = grid::unit(height_in, "in"),
      angle = -90
    )
  )
  list(grob = rotated, width = height_in, height = width_in)
}
