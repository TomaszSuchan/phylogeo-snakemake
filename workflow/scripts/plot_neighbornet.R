#!/usr/bin/env Rscript

# Tip labels: place text along the edge direction (parent → tip), starting just past the tip.
# ggtree::geom_tiplab2 uses a horizontal x-nudge that does not follow slanted split-network edges.
TIP_LABEL_OFFSET_FRAC <- 0 # fraction of x/y span; offset along edge, outward from parent
TIP_LABEL_SIZE_MM <- 2.3      # geom_text font size in mm (ggplot2)

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggtree)
  library(tanggle)
})

pdf(NULL)

log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

net_file <- snakemake@input[["net"]]
indpopdata_file <- snakemake@input[["indpopdata"]]
output_pdf_with <- snakemake@output[["pdf_with_tip_labels"]]
output_pdf_no <- snakemake@output[["pdf_no_tip_labels"]]
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

n_groups <- length(levels(tip_meta$group))
cat("Number of tips:", length(tip_labels), "\n")
cat("Coloring tips by:", color_by, "\n")
cat("Groups:", paste(levels(tip_meta$group), collapse = ", "), "\n")

p_net <- tanggle::ggsplitnet(net) %<+% tip_meta

p_base <- p_net +
  ggtree::geom_tippoint(aes(color = group), size = 2) +
  ggplot2::labs(color = color_by) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "right")

if (!is.null(user_colors) && length(user_colors) > 0) {
  if (n_groups > length(user_colors)) {
    palette_vals <- grDevices::colorRampPalette(user_colors)(n_groups)
  } else {
    palette_vals <- user_colors[seq_len(n_groups)]
  }
  names(palette_vals) <- levels(tip_meta$group)
  p_base <- p_base + ggplot2::scale_color_manual(values = palette_vals)
}

# Tip label positions: follow edge direction (see tanggle fortify.networx: angle = atan2(y-yend, x-xend))
d <- p_net$data
need <- c("x", "y", "angle", "isTip", "label", "node")
if (!all(need %in% names(d))) {
  stop(
    "Unexpected ggsplitnet plot data (missing columns). Have: ",
    paste(names(d), collapse = ", ")
  )
}

tips <- d[as.logical(d$isTip) & !is.na(d$label) & nzchar(as.character(d$label)), , drop = FALSE]
tips <- tips[!duplicated(tips$node), , drop = FALSE]

span <- max(
  diff(range(d$x, na.rm = TRUE)),
  diff(range(d$y, na.rm = TRUE)),
  na.rm = TRUE
)
if (!is.finite(span) || span <= 0) span <- 1
eps <- span * TIP_LABEL_OFFSET_FRAC

theta_rad <- tips$angle * pi / 180
tips$x_lab <- tips$x + eps * cos(theta_rad)
tips$y_lab <- tips$y + eps * sin(theta_rad)

# Same hemisphere rule as ggtree::geom_tiplab2: readable text orientation
flip <- tips$angle >= 90 & tips$angle <= 270
tips$text_angle <- tips$angle
tips$text_angle[flip] <- tips$angle[flip] + 180
tips$hjust <- ifelse(flip, 1, 0)

p_with <- p_base +
  ggplot2::geom_text(
    data = tips,
    mapping = aes(
      x = x_lab,
      y = y_lab,
      label = label,
      color = group,
      angle = text_angle,
      hjust = hjust
    ),
    inherit.aes = FALSE,
    size = TIP_LABEL_SIZE_MM,
    lineheight = 0.9
  )

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

saveRDS(
  list(
    with_tip_labels = p_with,
    no_tip_labels = p_base,
    net = net,
    tip_metadata = tip_meta,
    tip_label_positions = tips,
    color_by = color_by
  ),
  output_rds
)

cat("Saved:", output_pdf_with, "\n")
cat("Saved:", output_pdf_no, "\n")
cat("Saved:", output_rds, "\n")
