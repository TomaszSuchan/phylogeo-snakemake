#!/usr/bin/env Rscript

# NeighborNet plotting.
#
# Edge and tip coordinates are taken from phangorn::coords(net, dim = "2D") --
# the same planar split-graph coordinates the reference plot.networx renderer
# uses. NeighborNet always produces a circular (planar) split system, so the
# split "boxes" close by construction when drawn from these coordinates.
# tanggle::ggsplitnet was used previously but its own coordinate handling left
# the boxes open; we now draw the ggplot directly with geom_segment instead, and
# additionally emit the base phangorn plot.networx figure as a reference.

TIP_LABEL_OFFSET_FRAC <- 0.02 # fraction of x/y span; offset along edge, outward from parent
TIP_LABEL_SIZE_MM <- 2.3      # geom_text font size in mm (ggplot2)
EDGE_COLOUR <- "grey20"       # split-network edge colour (SplitsTree-like thin dark lines)
EDGE_LINEWIDTH <- 0.3         # split-network edge width

# Extra room so angled tip names are not clipped (data limits + panel + outer margins)
PLOT_EXPAND_MULT_WITH_LABELS <- 0.12   # scale_x/y expansion multipliers when labels are shown
PLOT_EXPAND_MULT_NO_LABELS <- 0.06     # smaller padding for the no-label PDF
PLOT_MARGIN_MM <- 14                   # theme plot.margin on all sides

suppressPackageStartupMessages({
  library(ggplot2)
  library(phangorn)
})

pdf(NULL)

log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

net_file <- snakemake@input[["net"]]
indpopdata_file <- snakemake@input[["indpopdata"]]
output_pdf_with <- snakemake@output[["pdf_with_tip_labels"]]
output_pdf_no <- snakemake@output[["pdf_no_tip_labels"]]
output_pdf_phangorn <- snakemake@output[["pdf_phangorn"]]
output_rds <- snakemake@output[["rds"]]

color_by <- as.character(snakemake@params[["color_by"]])
plot_width <- as.numeric(snakemake@params[["width"]])
plot_height <- as.numeric(snakemake@params[["height"]])
plot_dpi <- as.numeric(snakemake@params[["dpi"]])

user_colors <- NULL
if ("neighbornet_colors" %in% names(snakemake@params)) {
  user_colors <- snakemake@params[["neighbornet_colors"]]
  if (!is.null(user_colors)) {
    user_colors <- unlist(user_colors)
  }
}

cat("Reading NeighborNet object:", net_file, "\n")
net <- readRDS(net_file)

tip_labels <- net$tip.label
if (is.null(tip_labels) || length(tip_labels) == 0) {
  stop("Could not extract tip labels from NeighborNet object.")
}
ntip <- length(tip_labels)

cat("Reading indpopdata:", indpopdata_file, "\n")
indpopdata <- read.table(
  indpopdata_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!("Ind" %in% colnames(indpopdata))) {
  stop("indpopdata must contain an 'Ind' column.")
}

if (!(color_by %in% colnames(indpopdata))) {
  stop(
    "Configured neighbornet color_by column '", color_by,
    "' is not present in indpopdata. Available columns: ",
    paste(colnames(indpopdata), collapse = ", ")
  )
}

tip_meta <- data.frame(label = tip_labels, stringsAsFactors = FALSE)
tip_meta <- merge(
  tip_meta,
  indpopdata[, c("Ind", color_by), drop = FALSE],
  by.x = "label",
  by.y = "Ind",
  all.x = TRUE,
  sort = FALSE
)
colnames(tip_meta)[colnames(tip_meta) == color_by] <- "group"
tip_meta$group[is.na(tip_meta$group) | tip_meta$group == ""] <- "unassigned"
tip_meta$group <- as.factor(tip_meta$group)
# Keep tip_meta row order aligned with the network's tip ordering (1..ntip).
tip_meta <- tip_meta[match(tip_labels, tip_meta$label), , drop = FALSE]

n_groups <- length(levels(tip_meta$group))
cat("Number of tips:", ntip, "\n")
cat("Coloring tips by:", color_by, "\n")
cat("Groups:", paste(levels(tip_meta$group), collapse = ", "), "\n")

# Resolve a named colour palette (shared by the ggplot and the base-R figure).
palette_vals <- NULL
if (!is.null(user_colors) && length(user_colors) > 0) {
  if (n_groups > length(user_colors)) {
    palette_vals <- grDevices::colorRampPalette(user_colors)(n_groups)
  } else {
    palette_vals <- user_colors[seq_len(n_groups)]
  }
  names(palette_vals) <- levels(tip_meta$group)
}

# --- Planar split-graph coordinates from phangorn (boxes close by construction) ---
xy <- phangorn::coords(net, dim = "2D")
if (is.null(dim(xy)) || ncol(xy) < 2) {
  stop("phangorn::coords(net, dim = '2D') did not return a 2-column coordinate matrix.")
}
edge <- net$edge
if (is.null(edge) || ncol(edge) != 2) {
  stop("NeighborNet object has no usable $edge matrix.")
}
if (max(edge) > nrow(xy)) {
  stop("Edge indices exceed the number of coordinate vertices; cannot draw network.")
}

# Every split edge connects two graph vertices; drawing all of them reproduces
# the closed parallelogram "boxes" of the split network.
seg_df <- data.frame(
  x = xy[edge[, 1], 1],
  y = xy[edge[, 1], 2],
  xend = xy[edge[, 2], 1],
  yend = xy[edge[, 2], 2]
)

# Tip vertices are nodes 1..ntip (networx follows the ape tip-numbering convention).
tip_df <- data.frame(
  label = tip_labels,
  x = xy[seq_len(ntip), 1],
  y = xy[seq_len(ntip), 2],
  group = tip_meta$group,
  stringsAsFactors = FALSE
)

# Place each tip label just past the tip, along the direction of its terminal
# edge (neighbour vertex -> tip), and flip text in the left hemisphere so it
# stays readable -- same convention as ggtree::geom_tiplab2.
parent_of_tip <- vapply(seq_len(ntip), function(i) {
  er <- which(edge[, 1] == i | edge[, 2] == i)[1]
  if (is.na(er)) return(NA_integer_)
  if (edge[er, 1] == i) edge[er, 2] else edge[er, 1]
}, integer(1))

span <- max(diff(range(xy[, 1])), diff(range(xy[, 2])))
if (!is.finite(span) || span <= 0) span <- 1
eps <- span * TIP_LABEL_OFFSET_FRAC

dx <- tip_df$x - xy[parent_of_tip, 1]
dy <- tip_df$y - xy[parent_of_tip, 2]
ang <- atan2(dy, dx)
ang[is.na(ang)] <- 0
tip_df$x_lab <- tip_df$x + eps * cos(ang)
tip_df$y_lab <- tip_df$y + eps * sin(ang)

deg <- ang * 180 / pi
flip <- deg > 90 | deg < -90
tip_df$text_angle <- ifelse(flip, deg + 180, deg)
tip_df$hjust <- ifelse(flip, 1, 0)

# --- ggplot figures (edges from phangorn coordinates) ---
build_base <- function() {
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = seg_df,
      mapping = ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      colour = EDGE_COLOUR,
      linewidth = EDGE_LINEWIDTH
    ) +
    ggplot2::geom_point(
      data = tip_df,
      mapping = ggplot2::aes(x = x, y = y, color = group),
      size = 2
    ) +
    ggplot2::labs(color = color_by) +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "right")
  if (!is.null(palette_vals)) {
    p <- p + ggplot2::scale_color_manual(values = palette_vals)
  }
  p
}

p_base <- build_base()
p_with <- p_base +
  ggplot2::geom_text(
    data = tip_df,
    mapping = ggplot2::aes(
      x = x_lab,
      y = y_lab,
      label = label,
      color = group,
      angle = text_angle,
      hjust = hjust
    ),
    inherit.aes = FALSE,
    size = TIP_LABEL_SIZE_MM,
    lineheight = 0.9,
    show.legend = FALSE
  )

finalize_plot_margins <- function(p, expand_mult) {
  m <- ggplot2::margin(PLOT_MARGIN_MM, PLOT_MARGIN_MM, PLOT_MARGIN_MM, PLOT_MARGIN_MM, unit = "mm")
  p +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(expand_mult, expand_mult))) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(expand_mult, expand_mult))) +
    ggplot2::coord_fixed(clip = "off") +
    ggplot2::theme(plot.margin = m)
}

p_with <- finalize_plot_margins(p_with, PLOT_EXPAND_MULT_WITH_LABELS)
p_base <- finalize_plot_margins(p_base, PLOT_EXPAND_MULT_NO_LABELS)

dir.create(dirname(output_pdf_with), recursive = TRUE, showWarnings = FALSE)

ggsave(
  filename = output_pdf_with,
  plot = p_with,
  width = plot_width,
  height = plot_height,
  dpi = plot_dpi,
  units = "in",
  limitsize = FALSE
)

ggsave(
  filename = output_pdf_no,
  plot = p_base,
  width = plot_width,
  height = plot_height,
  dpi = plot_dpi,
  units = "in",
  limitsize = FALSE
)

# --- Reference base-R figure via phangorn::plot.networx ---
# This is phangorn's own renderer (the canonical NeighborNet drawing); it draws
# every split and closes the boxes, and serves as a cross-check on the ggplot.
tip_cols <- if (!is.null(palette_vals)) {
  unname(palette_vals[as.character(tip_meta$group)])
} else {
  rep("black", ntip)
}

pdf(output_pdf_phangorn, width = plot_width, height = plot_height)
plotted <- tryCatch({
  plot(net, type = "2D", show.tip.label = TRUE, tip.color = tip_cols, cex = 0.6)
  TRUE
}, error = function(e) {
  message("plot.networx with tip.color failed (", conditionMessage(e),
          "); retrying without tip colouring.")
  FALSE
})
if (!isTRUE(plotted)) {
  plot(net, type = "2D", show.tip.label = TRUE, cex = 0.6)
}
dev.off()

saveRDS(
  list(
    with_tip_labels = p_with,
    no_tip_labels = p_base,
    net = net,
    tip_metadata = tip_meta,
    tip_label_positions = tip_df,
    edges = seg_df,
    color_by = color_by
  ),
  output_rds
)

cat("Saved:", output_pdf_with, "\n")
cat("Saved:", output_pdf_no, "\n")
cat("Saved:", output_pdf_phangorn, "\n")
cat("Saved:", output_rds, "\n")
