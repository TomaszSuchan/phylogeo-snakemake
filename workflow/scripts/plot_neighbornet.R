#!/usr/bin/env Rscript

# NeighborNet plotting via tanggle::ggsplitnet (standard phangorn/ggtree workflow).
# See tanggle vignette: ggsplitnet(net) + geom_tiplab2().

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggtree)
  library(tanggle)
  library(phangorn)
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

group_utils <- tryCatch(
  file.path(dirname(normalizePath(snakemake@script)), "plot_group_utils.R"),
  error = function(e) "workflow/scripts/plot_group_utils.R"
)
if (file.exists(group_utils)) {
  source(group_utils)
} else {
  source("workflow/scripts/plot_group_utils.R")
}

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
group_colors_param <- snakemake@params[["group_colors"]]

cat("Reading NeighborNet object:", net_file, "\n")
net <- readRDS(net_file)

if (is.null(net$Ntip) || length(net$Ntip) != 1L || !is.finite(net$Ntip)) {
  net$Ntip <- length(net$tip.label)
}

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

tip_data <- data.frame(label = tip_labels, stringsAsFactors = FALSE)
tip_data <- merge(
  tip_data,
  indpopdata[, c("Ind", color_by), drop = FALSE],
  by.x = "label",
  by.y = "Ind",
  all.x = TRUE,
  sort = FALSE
)
colnames(tip_data)[colnames(tip_data) == color_by] <- "group"
tip_data$group[is.na(tip_data$group) | tip_data$group == ""] <- "unassigned"
tip_data <- tip_data[match(tip_labels, tip_data$label), , drop = FALSE]
tip_data$group <- factor(tip_data$group, levels = sort(unique(as.character(tip_data$group))))

n_groups <- length(levels(tip_data$group))
cat("Number of tips:", ntip, "\n")
cat("Coloring tips by:", color_by, "\n")
cat("Groups:", paste(levels(tip_data$group), collapse = ", "), "\n")

palette_vals <- group_fill_values(group_colors_param)
if (is.null(palette_vals) && !is.null(group_colors_param)) {
  unnamed <- unlist(group_colors_param, use.names = FALSE)
  unnamed <- unnamed[!is.na(unnamed) & unnamed != ""]
  if (length(unnamed) > 0) {
    if (n_groups > length(unnamed)) {
      palette_vals <- grDevices::colorRampPalette(unnamed)(n_groups)
    } else {
      palette_vals <- unnamed[seq_len(n_groups)]
    }
    names(palette_vals) <- levels(tip_data$group)
  }
}

tip_label_size <- if (ntip > 120) 1.4 else if (ntip > 60) 1.8 else 2.2
expand_frac <- if (ntip > 120) 0.18 else if (ntip > 60) 0.14 else 0.10

build_plot <- function(show_tip_labels) {
  p <- ggsplitnet(net) %<+% tip_data +
    theme_tree() +
    ggexpand(expand_frac) +
    ggexpand(expand_frac, direction = -1)

  if (show_tip_labels) {
    p <- p +
      geom_tiplab2(
        mapping = aes(color = group),
        size = tip_label_size
      ) +
      labs(color = color_by) +
      theme(legend.position = "right")
  } else {
    p <- p + theme(legend.position = "none")
  }

  if (!is.null(palette_vals)) {
    p <- p + scale_color_manual(values = palette_vals, drop = FALSE)
  }

  p
}

p_with <- build_plot(show_tip_labels = TRUE)
p_base <- build_plot(show_tip_labels = FALSE)

dir.create(dirname(output_pdf_with), recursive = TRUE, showWarnings = FALSE)

ggsave_pdf(
  filename = output_pdf_with,
  plot = p_with,
  width = plot_width,
  height = plot_height,
  dpi = plot_dpi,
  units = "in",
  limitsize = FALSE
)

ggsave_pdf(
  filename = output_pdf_no,
  plot = p_base,
  width = plot_width,
  height = plot_height,
  dpi = plot_dpi,
  units = "in",
  limitsize = FALSE
)

tip_cols <- if (!is.null(palette_vals)) {
  vals <- palette_vals[as.character(tip_data$group)]
  vals[is.na(vals)] <- "grey40"
  unname(vals)
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
    tip_metadata = tip_data,
    color_by = color_by
  ),
  output_rds
)

cat("Saved:", output_pdf_with, "\n")
cat("Saved:", output_pdf_no, "\n")
cat("Saved:", output_pdf_phangorn, "\n")
cat("Saved:", output_rds, "\n")
