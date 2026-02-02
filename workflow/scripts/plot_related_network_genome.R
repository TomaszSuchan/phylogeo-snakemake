#!/usr/bin/env Rscript
# Plot relatedness network graph based on PLINK --genome results
# Uses ggraph + tidygraph for professional network visualization

library(ggplot2)
library(dplyr)
library(igraph)
library(tidygraph)
library(ggraph)

# Prevent creation of Rplots.pdf
pdf(NULL)

# Redirect all output to log file
log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

# Snakemake inputs/outputs
genome_file <- snakemake@input[["genome"]]
indpopdata_file <- snakemake@input[["indpopdata"]]
output_pdf <- snakemake@output[["pdf"]]
output_rds <- snakemake@output[["rds"]]

# Get color_by and relatedness_colors parameters
color_by_name <- NULL
relatedness_colors <- NULL
if ("color_by" %in% names(snakemake@params)) {
  color_by_name <- as.character(snakemake@params[["color_by"]])
}
if ("relatedness_colors" %in% names(snakemake@params)) {
  relatedness_colors <- snakemake@params[["relatedness_colors"]]
  if (!is.null(relatedness_colors)) {
    relatedness_colors <- unlist(relatedness_colors)
  }
}

# Relationship thresholds (PI_HAT values)
# Adjusted to better distinguish high PI_HAT values
clone_threshold <- 0.9      # PI_HAT > 0.9 for clones/twins (very high relatedness)
first_degree_threshold <- 0.45  # 0.45 < PI_HAT <= 0.9 for 1st-degree (high relatedness)
second_degree_threshold <- 0.25  # 0.25 < PI_HAT <= 0.45 for 2nd-degree
# 0.05 < PI_HAT <= 0.25 is "other" (but we'll filter out unrelated pairs with PI_HAT <= 0.05)

# Color scheme
color_scheme <- c(
  "clone" = "black",
  "1st-degree" = "black",
  "2nd-degree" = "gray50",
  "other" = "gray85"
)

# Linetype scheme
linetype_scheme <- c(
  "clone" = "solid",
  "1st-degree" = "dashed",
  "2nd-degree" = "dashed",
  "other" = "dotted"
)

message("\n=== READING RELATEDNESS DATA ===\n")
# PLINK --genome output uses space-separated columns with variable spacing
genome_df <- read.table(genome_file, header = TRUE, sep = "", stringsAsFactors = FALSE, 
                        comment.char = "", check.names = FALSE, strip.white = TRUE)
message(sprintf("Loaded %d pairwise relationships\n", nrow(genome_df)))

# Filter out self-comparisons and unrelated pairs (PI_HAT <= 0.05)
genome_df <- genome_df %>%
  filter(IID1 != IID2, PI_HAT > 0.05)

message(sprintf("After filtering: %d related pairs\n", nrow(genome_df)))

if (nrow(genome_df) == 0) {
  message("No related pairs found. Creating empty plot.\n")
  p <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = "No related pairs found", size = 6) +
    theme_void()
  
  dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)
  ggsave(
    filename = output_pdf,
    plot = p,
    width = 10,
    height = 8,
    units = "in"
  )
  saveRDS(p, output_rds)
  quit(status = 0)
}

# Categorize relationships based on PI_HAT
categorize_relationship <- function(pi_hat) {
  if (pi_hat > clone_threshold) {
    return("clone")
  } else if (pi_hat > first_degree_threshold) {
    return("1st-degree")
  } else if (pi_hat > second_degree_threshold) {
    return("2nd-degree")
  } else {
    return("other")
  }
}

genome_df$category <- sapply(genome_df$PI_HAT, categorize_relationship)
genome_df$category <- factor(genome_df$category, levels = names(color_scheme))

message("\n=== RELATIONSHIP CATEGORIES ===\n")
category_counts <- table(genome_df$category)
print(category_counts)

# Get all unique individuals
all_individuals <- unique(c(genome_df$IID1, genome_df$IID2))
n_nodes <- length(all_individuals)
message(sprintf("\nTotal individuals: %d\n", n_nodes))

# Create edge data frame
edges_df <- data.frame(
  from = genome_df$IID1,
  to = genome_df$IID2,
  weight = genome_df$PI_HAT,
  category = genome_df$category,
  stringsAsFactors = FALSE
)

# Create node data frame
nodes_df <- data.frame(
  name = all_individuals,
  stringsAsFactors = FALSE
)

# Read indpopdata and merge with nodes if available
indpopdata <- NULL
if (file.exists(indpopdata_file)) {
  message("\n=== READING INDPOPDATA ===\n")
  indpopdata <- read.table(indpopdata_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  message(sprintf("Loaded indpopdata with %d samples and %d columns\n", nrow(indpopdata), ncol(indpopdata)))
  message("Columns:", paste(colnames(indpopdata), collapse = ", "), "\n")
  
  # Get the first column name (should be "Ind" or "Sample")
  ind_col <- colnames(indpopdata)[1]
  
  # Merge with nodes
  nodes_df <- merge(nodes_df, indpopdata, by.x = "name", by.y = ind_col, all.x = TRUE)
  message(sprintf("Merged indpopdata with nodes: %d nodes\n", nrow(nodes_df)))
}

# Create tidygraph
message("Creating network graph...\n")
tbl_graph <- tbl_graph(nodes = nodes_df, edges = edges_df, directed = FALSE)

message(sprintf("Graph created: %d nodes, %d edges\n", 
                gorder(tbl_graph), gsize(tbl_graph)))

# Calculate layout using ggraph - automatically handles disconnected components
message("Calculating layout...\n")
layout <- create_layout(tbl_graph, layout = "fr", niter = 2000)

message("Creating plot...\n")

# Determine node coloring
node_color_col <- NULL
node_colors <- NULL
if (!is.null(color_by_name) && color_by_name != "" && color_by_name != "none" && !is.null(indpopdata)) {
  if (color_by_name %in% colnames(nodes_df)) {
    node_color_col <- color_by_name
    message(sprintf("Coloring nodes by: %s\n", node_color_col))
    
    # Remove NA and empty values from color column
    na_mask <- is.na(nodes_df[[node_color_col]])
    if (any(na_mask)) {
      message(sprintf("Warning: %d nodes have NA values in '%s', will be colored gray\n", sum(na_mask), node_color_col))
    }
    empty_mask <- nodes_df[[node_color_col]] == ""
    if (any(empty_mask)) {
      message(sprintf("Warning: %d nodes have empty values in '%s', will be colored gray\n", sum(empty_mask), node_color_col))
    }
    
    # Determine colors to use
    unique_vals <- unique(nodes_df[[node_color_col]][!is.na(nodes_df[[node_color_col]]) & nodes_df[[node_color_col]] != ""])
    n_categories <- length(unique_vals)
    message(sprintf("Found %d unique values in '%s': %s\n", n_categories, node_color_col, paste(unique_vals, collapse = ", ")))
    
    if (!is.null(relatedness_colors) && length(relatedness_colors) > 0) {
      if (n_categories > length(relatedness_colors)) {
        node_colors <- colorRampPalette(relatedness_colors)(n_categories)
      } else {
        node_colors <- relatedness_colors[1:n_categories]
      }
      message(sprintf("Using %d colors from config\n", length(node_colors)))
    } else {
      message("No colors defined in config. Using default ggplot2 colors.\n")
    }
  } else {
    message(sprintf("Warning: Column '%s' not found in indpopdata. Using default node colors.\n", color_by_name))
  }
}

# Create the plot using ggraph
# Check if color column exists in layout (layout should have all node attributes)
if (!is.null(node_color_col) && node_color_col %in% colnames(layout)) {
  # Color nodes by specified column
  p <- ggraph(layout) +
    # Draw edges
    geom_edge_link(
      aes(color = category, linetype = category),
      width = 0.8,
      alpha = 0.75,
      lineend = "round"
    ) +
    # Draw nodes with color mapping
    geom_node_point(
      aes(fill = .data[[node_color_col]]),
      color = "black",
      size = 2,
      shape = 21,
      stroke = 0.5
    ) +
    # Add node labels
    geom_node_text(
      aes(label = name),
      size = 2.5,
      hjust = 0.5,
      vjust = 1.5,
      color = "black",
      repel = FALSE
    ) +
    # Add color scale for nodes if colors were defined
    (if (!is.null(node_colors)) {
      scale_fill_manual(
        values = node_colors,
        name = node_color_col,
        na.value = "gray90",
        guide = guide_legend(override.aes = list(size = 3))
      )
    } else {
      scale_fill_discrete(
        name = node_color_col,
        na.value = "gray90",
        guide = guide_legend(override.aes = list(size = 3))
      )
    }) +
    # Color scale for edges
    scale_edge_color_manual(
      values = color_scheme,
      name = "Relationship Type",
      drop = FALSE
    ) +
    # Linetype scale for edges
    scale_edge_linetype_manual(
      values = linetype_scheme,
      name = "Relationship Type",
      drop = FALSE
    ) +
    labs(
      title = "Genome Network Graph (PI_HAT)",
      subtitle = paste("Nodes:", n_nodes, "| Edges:", nrow(edges_df))
    ) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      legend.position = "bottom",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9)
    )
} else {
  # Default: no coloring (gray nodes)
  p <- ggraph(layout) +
    # Draw edges
    geom_edge_link(
      aes(color = category, linetype = category),
      width = 0.8,
      alpha = 0.75,
      lineend = "round"
    ) +
    # Draw nodes
    geom_node_point(
      color = "lightblue",
      fill = "lightblue",
      size = 2,
      shape = 21,
      stroke = 0.5
    ) +
    # Add node labels
    geom_node_text(
      aes(label = name),
      size = 2.5,
      hjust = 0.5,
      vjust = 1.5,
      color = "black",
      repel = FALSE
    ) +
    # Color scale for edges
    scale_edge_color_manual(
      values = color_scheme,
      name = "Relationship Type",
      drop = FALSE
    ) +
    # Linetype scale for edges
    scale_edge_linetype_manual(
      values = linetype_scheme,
      name = "Relationship Type",
      drop = FALSE
    ) +
    labs(
      title = "Genome Network Graph (PI_HAT)",
      subtitle = paste("Nodes:", n_nodes, "| Edges:", nrow(edges_df))
    ) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      legend.position = "bottom",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9)
    )
}

# Create output directory
dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)

# Save plot
ggsave(
  filename = output_pdf,
  plot = p,
  width = 12,
  height = 10,
  units = "in"
)

message(sprintf("Plot saved to: %s\n", output_pdf))

# Save data for RDS
graph_data <- list(
  graph = tbl_graph,
  layout = layout,
  nodes = nodes_df,
  edges = edges_df,
  color_scheme = color_scheme,
  plot = p
)
saveRDS(graph_data, output_rds)

message(sprintf("Graph data saved to: %s\n", output_rds))
message("=== DONE ===\n")
