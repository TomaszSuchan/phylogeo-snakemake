#!/usr/bin/env Rscript
# Plot relatedness network graph based on PLINK --genome results
# Creates a network visualization where edges are colored by relationship category
# Uses ggplot2 only (no igraph dependency)

library(ggplot2)
library(dplyr)

# Prevent creation of Rplots.pdf
pdf(NULL)

# Redirect all output to log file
log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

# Snakemake inputs/outputs
genome_file <- snakemake@input[["genome"]]
output_pdf <- snakemake@output[["pdf"]]
output_rds <- snakemake@output[["rds"]]

# Relationship thresholds (PI_HAT values)
# Adjusted to better distinguish high PI_HAT values
clone_threshold <- 0.9      # PI_HAT > 0.9 for clones/twins (very high relatedness)
first_degree_threshold <- 0.45  # 0.45 < PI_HAT <= 0.9 for 1st-degree (high relatedness)
second_degree_threshold <- 0.25  # 0.25 < PI_HAT <= 0.45 for 2nd-degree
# 0.05 < PI_HAT <= 0.25 is "other" (but we'll filter out unrelated pairs with PI_HAT <= 0.05)

# Color scheme
color_scheme <- c(
  "clone" = "black",
  "1st-degree" = "gray20",
  "2nd-degree" = "gray50",
  "other" = "gray85"
)

# Linetype scheme (helps distinguish edges even when colors are similar)
linetype_scheme <- c(
  "clone" = "solid",
  "1st-degree" = "longdash",
  "2nd-degree" = "dashed",
  "other" = "dotted"
)

message("\n=== READING RELATEDNESS DATA ===\n")
# PLINK --genome output uses space-separated columns with variable spacing
# Read with skip=0 to handle header, strip.white=TRUE to handle leading spaces
genome_df <- read.table(genome_file, header = TRUE, sep = "", stringsAsFactors = FALSE, 
                        comment.char = "", check.names = FALSE, strip.white = TRUE)
message(sprintf("Loaded %d pairwise relationships\n", nrow(genome_df)))
message("Columns:", paste(colnames(genome_df), collapse = ", "), "\n")

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

# Create node data frame
nodes_df <- data.frame(
  id = all_individuals,
  stringsAsFactors = FALSE
)

# Create edge data frame
edges_df <- data.frame(
  from = genome_df$IID1,
  to = genome_df$IID2,
  weight = genome_df$PI_HAT,
  category = genome_df$category,
  stringsAsFactors = FALSE
)

# Simple force-directed layout algorithm (Fruchterman-Reingold inspired)
# This is a simplified version that works reasonably well for small to medium networks
calculate_layout <- function(nodes, edges, n_iter = 1000) {
  n <- nrow(nodes)
  
  # Initialize positions randomly
  pos <- matrix(runif(n * 2, -1, 1), ncol = 2)
  colnames(pos) <- c("x", "y")
  
  # Create adjacency matrix for faster lookups
  node_ids <- nodes$id
  names(node_ids) <- 1:n
  
  # Convert node IDs to indices
  from_idx <- match(edges$from, node_ids)
  to_idx <- match(edges$to, node_ids)
  
  # Temperature for cooling
  temp <- 1.0
  cooling_rate <- 0.99
  
  for (iter in 1:n_iter) {
    # Calculate repulsive forces (all nodes repel each other)
    forces <- matrix(0, nrow = n, ncol = 2)
    
    for (i in 1:(n-1)) {
      for (j in (i+1):n) {
        dx <- pos[i, "x"] - pos[j, "x"]
        dy <- pos[i, "y"] - pos[j, "y"]
        dist <- sqrt(dx^2 + dy^2)
        
        if (dist > 0.001) {  # Avoid division by zero
          # Repulsive force (inverse square)
          force <- temp / (dist^2)
          forces[i, ] <- forces[i, ] + c(dx, dy) * force / dist
          forces[j, ] <- forces[j, ] - c(dx, dy) * force / dist
        }
      }
    }
    
    # Calculate attractive forces (only for connected nodes)
    for (k in 1:nrow(edges)) {
      i <- from_idx[k]
      j <- to_idx[k]
      
      if (!is.na(i) && !is.na(j)) {
        dx <- pos[i, "x"] - pos[j, "x"]
        dy <- pos[i, "y"] - pos[j, "y"]
        dist <- sqrt(dx^2 + dy^2)
        
        if (dist > 0.001) {
          # Attractive force (proportional to distance, weighted by edge weight)
          weight <- edges$weight[k]
          force <- dist * weight * 0.1
          forces[i, ] <- forces[i, ] - c(dx, dy) * force / dist
          forces[j, ] <- forces[j, ] + c(dx, dy) * force / dist
        }
      }
    }
    
    # Update positions
    pos <- pos + forces * 0.1
    
    # Cool down
    temp <- temp * cooling_rate
  }
  
  # Normalize to fit in [0, 1] range
  pos[, "x"] <- (pos[, "x"] - min(pos[, "x"])) / (max(pos[, "x"]) - min(pos[, "x"]) + 0.001)
  pos[, "y"] <- (pos[, "y"] - min(pos[, "y"])) / (max(pos[, "y"]) - min(pos[, "y"]) + 0.001)
  
  nodes$x <- pos[, "x"]
  nodes$y <- pos[, "y"]
  
  return(nodes)
}

message("Calculating layout...\n")
nodes_df <- calculate_layout(nodes_df, edges_df, n_iter = 1000)

# Merge node positions into edges
edges_plot <- edges_df %>%
  left_join(nodes_df, by = c("from" = "id")) %>%
  rename(x1 = x, y1 = y) %>%
  left_join(nodes_df, by = c("to" = "id")) %>%
  rename(x2 = x, y2 = y)

message("Creating plot...\n")

# Create the plot
p <- ggplot() +
  # Draw edges first (so they appear behind nodes)
  geom_segment(data = edges_plot,
               aes(x = x1, y = y1, xend = x2, yend = y2,
                   color = category, linetype = category),
               size = 0.8,
               alpha = 0.75,
               lineend = "round") +
  # Draw nodes
  geom_point(data = nodes_df,
             aes(x = x, y = y),
             color = "lightblue",
             fill = "lightblue",
             size = 2,
             shape = 21,
             stroke = 0.5) +
  # Add node labels
  geom_text(data = nodes_df,
            aes(x = x, y = y, label = id),
            size = 2.5,
            hjust = 0.5,
            vjust = 1.5,
            color = "black",
            check_overlap = FALSE,
            angle = 0) +
  # Color scale for edges
  scale_color_manual(values = color_scheme,
                    name = "Relationship Type",
                    drop = FALSE) +
  # Linetype scale for edges
  scale_linetype_manual(values = linetype_scheme,
                        name = "Relationship Type",
                        drop = FALSE,
                        guide = guide_legend(override.aes = list(color = color_scheme))) +
  labs(title = "Relatedness Network Graph",
       subtitle = paste("Nodes:", n_nodes, "| Edges:", nrow(edges_df))) +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  ) +
  coord_fixed()

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
  nodes = nodes_df,
  edges = edges_df,
  color_scheme = color_scheme,
  plot = p
)
saveRDS(graph_data, output_rds)

message(sprintf("Graph data saved to: %s\n", output_rds))
message("=== DONE ===\n")
