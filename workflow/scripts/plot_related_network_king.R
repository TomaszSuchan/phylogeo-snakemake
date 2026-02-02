#!/usr/bin/env Rscript
# Plot relatedness network graph based on PLINK2 KING kinship table
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
king_file <- snakemake@input[["king"]]
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

# KING kinship thresholds (standard values from your pipeline)
clone_threshold <- 0.354      # KING > 0.354 for duplicates/twins
first_degree_threshold <- 0.177  # 0.177 < KING <= 0.354 for 1st-degree
second_degree_threshold <- 0.0884  # 0.0884 < KING <= 0.177 for 2nd-degree
# KING <= 0.0884 is "other" (but we'll filter out unrelated pairs)

# Color scheme
color_scheme <- c(
  "clone" = "black",
  "1st-degree" = "gray20",
  "2nd-degree" = "gray50",
  "other" = "gray85"
)

# Linetype scheme
linetype_scheme <- c(
  "clone" = "solid",
  "1st-degree" = "longdash",
  "2nd-degree" = "dashed",
  "other" = "dotted"
)

message("\n=== READING KING DATA ===\n")
# PLINK2 KING table uses space-separated columns with variable spacing
king_df <- read.table(king_file, header = TRUE, sep = "", stringsAsFactors = FALSE, 
                      comment.char = "", check.names = FALSE, strip.white = TRUE)
message(sprintf("Loaded %d pairwise relationships\n", nrow(king_df)))
message("Columns:", paste(colnames(king_df), collapse = ", "), "\n")

# Filter out self-comparisons and unrelated pairs (KING <= 0.05)
king_df <- king_df %>%
  filter(IID1 != IID2, KINSHIP > 0.05)

message(sprintf("After filtering: %d related pairs\n", nrow(king_df)))

if (nrow(king_df) == 0) {
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

# Categorize relationships based on KING kinship
categorize_relationship <- function(kinship) {
  if (kinship > clone_threshold) {
    return("clone")
  } else if (kinship > first_degree_threshold) {
    return("1st-degree")
  } else if (kinship > second_degree_threshold) {
    return("2nd-degree")
  } else {
    return("other")
  }
}

king_df$category <- sapply(king_df$KINSHIP, categorize_relationship)
king_df$category <- factor(king_df$category, levels = names(color_scheme))

message("\n=== RELATIONSHIP CATEGORIES ===\n")
category_counts <- table(king_df$category)
print(category_counts)

# Get all unique individuals
all_individuals <- unique(c(king_df$IID1, king_df$IID2))
n_nodes <- length(all_individuals)
message(sprintf("\nTotal individuals: %d\n", n_nodes))

# Create edge data frame
edges_df <- data.frame(
  from = king_df$IID1,
  to = king_df$IID2,
  weight = king_df$KINSHIP,
  category = king_df$category,
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
      title = "KING Network Graph",
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
      title = "KING Network Graph",
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
