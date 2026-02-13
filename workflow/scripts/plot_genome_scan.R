#!/usr/bin/env Rscript
# Plot FST, pi, and dXY along the genome (like Figure 3A from Nature Communications paper)
# Window size: 100 bp for FST (dot plot), 1 Mb for π and dXY (curve visualization)

library(tidyverse)
library(gridExtra)

# Prevent creation of Rplots.pdf
pdf(NULL)

# Redirect all output to log file
log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

# Snakemake inputs/outputs
fst_file <- snakemake@input[["fst"]]
pi_file <- snakemake@input[["pi"]]
dxy_file <- snakemake@input[["dxy"]]
output_pdf <- snakemake@output[["pdf"]]
output_rds <- snakemake@output[["rds"]]
pop1 <- snakemake@params[["pop1"]]
pop2 <- snakemake@params[["pop2"]]

message("\n=== READING PIXY OUTPUT FILES ===\n")

# Read FST data
message("Reading FST file...")
fst_df <- read.table(fst_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
message(sprintf("  Loaded %d rows", nrow(fst_df)))
message(sprintf("  Columns: %s", paste(colnames(fst_df), collapse = ", ")))

# Read pi data
message("\nReading pi file...")
pi_df <- read.table(pi_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
message(sprintf("  Loaded %d rows", nrow(pi_df)))
message(sprintf("  Columns: %s", paste(colnames(pi_df), collapse = ", ")))

# Read dXY data
message("\nReading dXY file...")
dxy_df <- read.table(dxy_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
message(sprintf("  Loaded %d rows", nrow(dxy_df)))
message(sprintf("  Columns: %s", paste(colnames(dxy_df), collapse = ", ")))

# Process FST data
# Pixy FST output has columns: pop1, pop2, chromosome, window_pos_1, window_pos_2, no_sites, avg_pi_pop1, avg_pi_pop2, avg_dxy, fst
message("\n=== PROCESSING FST DATA ===\n")
fst_plot <- fst_df %>%
  filter(!is.na(fst)) %>%
  mutate(
    midpoint = (window_pos_1 + window_pos_2) / 2,
    chromosome = as.character(chromosome)
  ) %>%
  select(chromosome, midpoint, fst) %>%
  arrange(chromosome, midpoint)

message(sprintf("FST data: %d windows", nrow(fst_plot)))
message(sprintf("Chromosomes: %s", paste(unique(fst_plot$chromosome), collapse = ", ")))

# Process pi data
# Pixy pi output has columns: pop, chromosome, window_pos_1, window_pos_2, no_sites, avg_pi
message("\n=== PROCESSING PI DATA ===\n")
pi_plot <- pi_df %>%
  filter(!is.na(avg_pi)) %>%
  mutate(
    midpoint = (window_pos_1 + window_pos_2) / 2,
    chromosome = as.character(chromosome),
    population = pop
  ) %>%
  select(chromosome, midpoint, population, avg_pi) %>%
  arrange(chromosome, midpoint, population)

message(sprintf("Pi data: %d windows", nrow(pi_plot)))
message(sprintf("Populations: %s", paste(unique(pi_plot$population), collapse = ", ")))

# Process dXY data
# Pixy dXY output has columns: pop1, pop2, chromosome, window_pos_1, window_pos_2, no_sites, avg_dxy
message("\n=== PROCESSING DXY DATA ===\n")
dxy_plot <- dxy_df %>%
  filter(!is.na(avg_dxy)) %>%
  mutate(
    midpoint = (window_pos_1 + window_pos_2) / 2,
    chromosome = as.character(chromosome)
  ) %>%
  select(chromosome, midpoint, avg_dxy) %>%
  arrange(chromosome, midpoint)

message(sprintf("dXY data: %d windows", nrow(dxy_plot)))

# Create cumulative positions for plotting across chromosomes
message("\n=== CREATING CUMULATIVE POSITIONS ===\n")
chromosomes <- unique(c(fst_plot$chromosome, pi_plot$chromosome, dxy_plot$chromosome))
chromosomes <- sort(chromosomes)

# Calculate chromosome lengths and cumulative offsets
chr_lengths <- list()
cumulative_offset <- 0
chr_offsets <- numeric(length(chromosomes))
names(chr_offsets) <- chromosomes

for (chr in chromosomes) {
  # Get max position for this chromosome from all datasets
  max_pos_fst <- if (any(fst_plot$chromosome == chr)) max(fst_plot$midpoint[fst_plot$chromosome == chr], na.rm = TRUE) else 0
  max_pos_pi <- if (any(pi_plot$chromosome == chr)) max(pi_plot$midpoint[pi_plot$chromosome == chr], na.rm = TRUE) else 0
  max_pos_dxy <- if (any(dxy_plot$chromosome == chr)) max(dxy_plot$midpoint[dxy_plot$chromosome == chr], na.rm = TRUE) else 0
  
  chr_length <- max(max_pos_fst, max_pos_pi, max_pos_dxy, na.rm = TRUE)
  chr_lengths[[chr]] <- chr_length
  chr_offsets[chr] <- cumulative_offset
  cumulative_offset <- cumulative_offset + chr_length
}

message(sprintf("Total genome length: %.0f bp", cumulative_offset))
message(sprintf("Number of chromosomes: %d", length(chromosomes)))

# Add cumulative positions
fst_plot <- fst_plot %>%
  mutate(cumulative_pos = midpoint + chr_offsets[chromosome])

pi_plot <- pi_plot %>%
  mutate(cumulative_pos = midpoint + chr_offsets[chromosome])

dxy_plot <- dxy_plot %>%
  mutate(cumulative_pos = midpoint + chr_offsets[chromosome])

# Create chromosome boundaries for plotting
chr_boundaries <- data.frame(
  chromosome = chromosomes,
  start = chr_offsets,
  end = chr_offsets + unlist(chr_lengths),
  midpoint = chr_offsets + unlist(chr_lengths) / 2
)

# Create the plot
message("\n=== CREATING PLOT ===\n")

# Create three panels: FST, pi, dXY
p1 <- ggplot(fst_plot, aes(x = cumulative_pos, y = fst)) +
  geom_point(size = 0.5, alpha = 0.6, color = "black") +
  geom_vline(xintercept = chr_boundaries$end[-nrow(chr_boundaries)], 
             linetype = "dashed", color = "gray60", linewidth = 0.3) +
  scale_x_continuous(
    breaks = chr_boundaries$midpoint,
    labels = chr_boundaries$chromosome,
    expand = c(0.01, 0)
  ) +
  labs(
    x = "Chromosome",
    y = expression(F[ST]),
    title = sprintf("FST between %s and %s", pop1, pop2)
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 10, hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 9),
    panel.grid.minor = element_blank()
  )

# Pi plot with both populations
p2 <- ggplot(pi_plot, aes(x = cumulative_pos, y = avg_pi, color = population)) +
  geom_line(linewidth = 0.5, alpha = 0.8) +
  geom_vline(xintercept = chr_boundaries$end[-nrow(chr_boundaries)], 
             linetype = "dashed", color = "gray60", linewidth = 0.3) +
  scale_x_continuous(
    breaks = chr_boundaries$midpoint,
    labels = chr_boundaries$chromosome,
    expand = c(0.01, 0)
  ) +
  scale_color_manual(
    values = c("#377EB8", "#E41A1C"),  # Blue and red
    name = "Population"
  ) +
  labs(
    x = "Chromosome",
    y = expression(pi),
    title = sprintf("Nucleotide diversity (π) in %s and %s", pop1, pop2)
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 10, hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 9),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

# dXY plot
p3 <- ggplot(dxy_plot, aes(x = cumulative_pos, y = avg_dxy)) +
  geom_line(linewidth = 0.5, alpha = 0.8, color = "black") +
  geom_vline(xintercept = chr_boundaries$end[-nrow(chr_boundaries)], 
             linetype = "dashed", color = "gray60", linewidth = 0.3) +
  scale_x_continuous(
    breaks = chr_boundaries$midpoint,
    labels = chr_boundaries$chromosome,
    expand = c(0.01, 0)
  ) +
  labs(
    x = "Chromosome",
    y = expression(d[XY]),
    title = sprintf("Absolute divergence (dXY) between %s and %s", pop1, pop2)
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 10, hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 9),
    panel.grid.minor = element_blank()
  )

# Combine plots using gridExtra
# Create combined plot
combined_plot <- grid.arrange(
  p1, p2, p3,
  ncol = 1,
  heights = c(1, 1, 1)
)

# Save plot
message(sprintf("\nSaving plot to %s...", output_pdf))
ggsave(output_pdf, combined_plot, width = 12, height = 10, dpi = 300)
message("Plot saved successfully")

# Save RDS
message(sprintf("Saving RDS to %s...", output_rds))
saveRDS(combined_plot, output_rds)
message("RDS saved successfully")

message("\n=== DONE ===\n")

