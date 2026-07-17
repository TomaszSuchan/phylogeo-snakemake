#!/usr/bin/env Rscript
# Generic relatedness network plot for pairwise kinship/relatedness tables.

library(ggplot2)
library(dplyr)
library(igraph)
library(tidygraph)
library(ggraph)

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
  sink(type = "message")
  sink(type = "output")
  close(log_file)
}, add = TRUE)

pairwise_file <- snakemake@input[["pairwise"]]
indpopdata_file <- snakemake@input[["indpopdata"]]
output_pdf <- snakemake@output[["pdf"]]
output_rds <- snakemake@output[["rds"]]

threshold_profile <- snakemake@params[["threshold_profile"]]
weight_column <- snakemake@params[["weight_column"]]

color_by_name <- NULL
relatedness_colors <- NULL
plot_all <- FALSE
if ("color_by" %in% names(snakemake@params)) {
  color_by_name <- as.character(snakemake@params[["color_by"]])
}
if ("relatedness_colors" %in% names(snakemake@params)) {
  relatedness_colors <- snakemake@params[["relatedness_colors"]]
  if (!is.null(relatedness_colors)) {
    relatedness_colors <- unlist(relatedness_colors)
  }
}
if ("plot_all" %in% names(snakemake@params)) {
  plot_all <- isTRUE(as.logical(snakemake@params[["plot_all"]]))
}

THRESHOLD_PROFILES <- list(
  ajk = list(
    clone = 0.45,
    first = 0.225,
    second = 0.1125,
    min = 0.05,
    legend = "RELATEDNESS_AJK",
    labels = c(
      "clone" = "high (0.45, 1.00]",
      "1st-degree" = "elevated (0.225, 0.45]",
      "2nd-degree" = "moderate (0.1125, 0.225]",
      "other" = "low (0.05, 0.1125]"
    )
  ),
  manichaikul = list(
    clone = 0.45,
    first = 0.225,
    second = 0.1125,
    min = 0.05,
    legend = "RELATEDNESS_PHI",
    labels = c(
      "clone" = "duplicate / twin (0.45, 0.50]",
      "1st-degree" = "1st-degree (0.225, 0.45]",
      "2nd-degree" = "2nd-degree (0.1125, 0.225]",
      "other" = "<2nd-degree (0.05, 0.1125]"
    )
  ),
  king = list(
    clone = 0.354,
    first = 0.177,
    second = 0.0884,
    min = 0.05,
    legend = "kinship",
    labels = c(
      "clone" = "clone / duplicate (0.354, 0.50]",
      "1st-degree" = "1st-degree (0.177, 0.354]",
      "2nd-degree" = "2nd-degree (0.0884, 0.177]",
      "other" = "<2nd-degree (0.05, 0.0884]"
    )
  ),
  pihat = list(
    clone = 0.9,
    first = 0.45,
    second = 0.25,
    min = 0.05,
    legend = "PI_HAT",
    labels = c(
      "clone" = "clone / duplicate (0.90, 1.00]",
      "1st-degree" = "1st-degree (0.45, 0.90]",
      "2nd-degree" = "2nd-degree (0.25, 0.45]",
      "other" = "<2nd-degree (0.05, 0.25]"
    )
  )
)

if (!threshold_profile %in% names(THRESHOLD_PROFILES)) {
  stop("Unknown threshold_profile: ", threshold_profile)
}
thr <- THRESHOLD_PROFILES[[threshold_profile]]

pick_column <- function(df, candidates) {
  cols <- colnames(df)
  norm <- tolower(gsub("[^a-z0-9]", "", cols))
  for (candidate in candidates) {
    hit <- which(norm == tolower(gsub("[^a-z0-9]", "", candidate)))
    if (length(hit) > 0) {
      return(cols[hit[1]])
    }
  }
  for (candidate in candidates) {
    hit <- grep(candidate, cols, ignore.case = TRUE, value = TRUE)
    if (length(hit) > 0) {
      return(hit[1])
    }
  }
  NULL
}

read_pairwise_table <- function(path) {
  if (grepl("\\.tsv$", path, ignore.case = TRUE)) {
    return(read.delim(path, stringsAsFactors = FALSE, check.names = FALSE))
  }
  read.table(
    path,
    header = TRUE,
    sep = "",
    stringsAsFactors = FALSE,
    comment.char = "",
    check.names = FALSE,
    strip.white = TRUE
  )
}

message("\n=== READING PAIRWISE DATA ===\n")
pairwise_df <- read_pairwise_table(pairwise_file)
pairwise_df <- as.data.frame(pairwise_df, stringsAsFactors = FALSE)
message(sprintf("Loaded %d rows from %s\n", nrow(pairwise_df), pairwise_file))
message("Columns: ", paste(colnames(pairwise_df), collapse = ", "), "\n")

id1_col <- pick_column(pairwise_df, c("INDV1", "IID1", "ID1", "Ind1", "ind1"))
id2_col <- pick_column(pairwise_df, c("INDV2", "IID2", "ID2", "Ind2", "ind2"))
if (is.null(id1_col) || is.null(id2_col)) {
  stop("Could not identify individual ID columns in pairwise table")
}

if (is.null(weight_column) || !nzchar(weight_column)) {
  weight_column <- pick_column(
    pairwise_df,
    c(
      "RELATEDNESS_AJK", "RELATEDNESS_PHI", "PI_HAT", "KINSHIP", "kin",
      "Quellergt", "quellergt", "Rxy", "wang", "lynchrd", "lynchli",
      "ritland", "trioml", "dyadml", "RELATEDNESS"
    )
  )
}
if (is.null(weight_column) || !weight_column %in% colnames(pairwise_df)) {
  stop("Could not identify relatedness/kinship column in pairwise table")
}

pairwise_df[[weight_column]] <- suppressWarnings(as.numeric(pairwise_df[[weight_column]]))
pairwise_df <- pairwise_df %>%
  filter(
    !is.na(.data[[id1_col]]),
    !is.na(.data[[id2_col]]),
    .data[[id1_col]] != .data[[id2_col]],
    !is.na(.data[[weight_column]]),
    .data[[weight_column]] > thr$min
  )

message(sprintf(
  "After filtering (%s > %s): %d related pairs\n",
  weight_column, thr$min, nrow(pairwise_df)
))

indpopdata <- NULL
ind_col <- NULL
if (file.exists(indpopdata_file)) {
  indpopdata <- read.table(indpopdata_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  ind_col <- colnames(indpopdata)[1]
}

categorize_relationship <- function(value) {
  if (value > thr$clone) {
    "clone"
  } else if (value > thr$first) {
    "1st-degree"
  } else if (value > thr$second) {
    "2nd-degree"
  } else {
    "other"
  }
}

color_scheme <- c(
  "clone" = "black",
  "1st-degree" = "black",
  "2nd-degree" = "gray50",
  "other" = "gray85"
)
linetype_scheme <- c(
  "clone" = "solid",
  "1st-degree" = "dashed",
  "2nd-degree" = "dashed",
  "other" = "dotted"
)
edge_legend_title <- paste0("Relationship type (", thr$legend, ")")

if (nrow(pairwise_df) > 0) {
  pairwise_df$category <- vapply(pairwise_df[[weight_column]], categorize_relationship, character(1))
  pairwise_df$category <- factor(pairwise_df$category, levels = names(color_scheme))
} else {
  pairwise_df$category <- factor(character(0), levels = names(color_scheme))
}

related_individuals <- unique(c(pairwise_df[[id1_col]], pairwise_df[[id2_col]]))
if (plot_all && !is.null(indpopdata) && !is.null(ind_col)) {
  metadata_individuals <- unique(indpopdata[[ind_col]])
  metadata_individuals <- metadata_individuals[!is.na(metadata_individuals) & metadata_individuals != ""]
  all_individuals <- unique(c(metadata_individuals, related_individuals))
} else {
  all_individuals <- related_individuals
}

if (length(all_individuals) == 0) {
  p <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = "No related pairs found", size = 6) +
    theme_void()
  dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)
  ggsave_pdf(filename = output_pdf, plot = p, width = 10, height = 8, units = "in")
  saveRDS(list(plot = p, edges = pairwise_df), output_rds)
  quit(status = 0)
}

edges_df <- data.frame(
  from = pairwise_df[[id1_col]],
  to = pairwise_df[[id2_col]],
  weight = pairwise_df[[weight_column]],
  category = pairwise_df$category,
  stringsAsFactors = FALSE
)

nodes_df <- data.frame(name = all_individuals, stringsAsFactors = FALSE)
if (!is.null(indpopdata) && !is.null(ind_col)) {
  nodes_df <- merge(nodes_df, indpopdata, by.x = "name", by.y = ind_col, all.x = TRUE, sort = FALSE)
}

tbl_graph <- tbl_graph(nodes = nodes_df, edges = edges_df, directed = FALSE)
layout <- create_layout(tbl_graph, layout = "fr", niter = 2000)

node_color_col <- NULL
node_colors <- NULL
if (!is.null(color_by_name) && color_by_name != "" && color_by_name != "none" && !is.null(indpopdata)) {
  if (color_by_name %in% colnames(nodes_df)) {
    node_color_col <- color_by_name
    unique_vals <- unique(nodes_df[[node_color_col]][!is.na(nodes_df[[node_color_col]]) & nodes_df[[node_color_col]] != ""])
    n_categories <- length(unique_vals)
    if (!is.null(relatedness_colors) && length(relatedness_colors) > 0) {
      node_colors <- if (n_categories > length(relatedness_colors)) {
        colorRampPalette(relatedness_colors)(n_categories)
      } else {
        relatedness_colors[seq_len(n_categories)]
      }
    }
  }
}

build_plot <- function() {
  p <- ggraph(layout) +
    geom_edge_link(
      aes(color = category, linetype = category),
      width = 0.8,
      alpha = 0.75,
      lineend = "round"
    ) +
    scale_edge_color_manual(
      values = color_scheme,
      name = edge_legend_title,
      labels = thr$labels,
      drop = FALSE
    ) +
    scale_edge_linetype_manual(
      values = linetype_scheme,
      name = edge_legend_title,
      labels = thr$labels,
      drop = FALSE
    ) +
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9)
    )

  if (!is.null(node_color_col) && node_color_col %in% colnames(layout)) {
    p <- p +
      geom_node_point(
        aes(fill = .data[[node_color_col]]),
        color = "black",
        size = 2,
        shape = 21,
        stroke = 0.5
      ) +
      geom_node_text(aes(label = name), size = 2.5, vjust = 1.5, color = "black")
    if (!is.null(node_colors)) {
      p <- p + scale_fill_manual(
        values = node_colors,
        name = node_color_col,
        na.value = "gray90",
        guide = guide_legend(order = 2, override.aes = list(size = 3))
      )
    } else {
      p <- p + scale_fill_discrete(
        name = node_color_col,
        na.value = "gray90",
        guide = guide_legend(order = 2, override.aes = list(size = 3))
      )
    }
  } else {
    p <- p +
      geom_node_point(color = "lightblue", fill = "lightblue", size = 2, shape = 21, stroke = 0.5) +
      geom_node_text(aes(label = name), size = 2.5, vjust = 1.5, color = "black")
  }
  p
}

p <- build_plot()
dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)
ggsave_pdf(filename = output_pdf, plot = p, width = 12, height = 10, units = "in")
saveRDS(
  list(
    graph = tbl_graph,
    layout = layout,
    nodes = nodes_df,
    edges = edges_df,
    weight_column = weight_column,
    threshold_profile = threshold_profile,
    plot = p
  ),
  output_rds
)

message(sprintf("Plot saved to: %s\n", output_pdf))
message(sprintf("Graph data saved to: %s\n", output_rds))
