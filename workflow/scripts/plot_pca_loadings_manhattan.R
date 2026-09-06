#!/usr/bin/env Rscript
# Manhattan-style plot of PCAone per-SNP loadings along the genome.
# Shows whether a PC is driven by markers spread across the genome (a flat
# "lawn", i.e. neutral/polygenic structure) or concentrated in a region /
# chromosome (a peak). Aesthetics mirror plot_genome_scan.R: cumulative genomic
# position on x, dashed chromosome boundaries, chromosome labels at midpoints.
#
# Colouring rule (per user request): points are coloured (alternating by
# chromosome) only when real chromosome data is available; if the loadings table
# carries no chromosome information (e.g. a non-model assembly exported to PLINK
# where every contig collapses to "0"), points are drawn in a single colour and
# the x-axis becomes a plain genome-order index.

library(tidyverse)

ggsave_utils <- tryCatch(
  file.path(dirname(normalizePath(snakemake@script)), "plot_ggsave_utils.R"),
  error = function(e) "workflow/scripts/plot_ggsave_utils.R"
)
if (file.exists(ggsave_utils)) {
  source(ggsave_utils)
} else {
  source("workflow/scripts/plot_ggsave_utils.R")
}

# Prevent creation of Rplots.pdf
pdf(NULL)

# Redirect all output to log file
log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

# Snakemake inputs/outputs/params
table_file <- snakemake@input[["table"]]
output_pdf <- snakemake@output[["pdf"]]
output_rds <- snakemake@output[["rds"]]

n_pc <- suppressWarnings(as.integer(snakemake@params[["n_pc"]]))
plot_width <- suppressWarnings(as.numeric(snakemake@params[["width"]]))
plot_height <- suppressWarnings(as.numeric(snakemake@params[["height"]]))
axis_title_size <- suppressWarnings(as.numeric(snakemake@params[["axis_title_size"]]))
axis_text_size <- suppressWarnings(as.numeric(snakemake@params[["axis_text_size"]]))
point_size <- suppressWarnings(as.numeric(snakemake@params[["point_size"]]))
if (is.na(n_pc) || n_pc < 1) n_pc <- 2
if (is.na(plot_width)) plot_width <- 12
if (is.na(plot_height)) plot_height <- 5.5
if (is.na(axis_title_size)) axis_title_size <- 10
if (is.na(axis_text_size)) axis_text_size <- 8
if (is.na(point_size)) point_size <- 0.6

message("=== READING LOADINGS TABLE ===")
df <- read.table(table_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                 colClasses = c(chrom = "character"))
message(sprintf("  Loaded %d rows; columns: %s", nrow(df), paste(colnames(df), collapse = ", ")))

# Keep only the leading PCs (PC1 .. PCn_pc), preserving numeric order.
pc_levels <- paste0("PC", seq_len(n_pc))
df <- df %>%
  filter(.data$pc %in% pc_levels) %>%
  mutate(
    pc = factor(.data$pc, levels = pc_levels),
    pos = suppressWarnings(as.numeric(.data$pos)),
    loading_sq = suppressWarnings(as.numeric(.data$loading_sq))
  )
if (nrow(df) == 0) stop("No rows for the requested PCs; check n_pc and the input table.")

# Number of SNPs (per PC) and the neutral baseline for squared loadings (1/M):
# squared loadings sum to 1 across all M SNPs, so the per-SNP average is 1/M.
M <- df %>% filter(.data$pc == pc_levels[1]) %>% nrow()
neutral_baseline <- 1 / M
message(sprintf("  M = %d SNPs per PC; neutral baseline (1/M) = %.3g", M, neutral_baseline))

# Decide whether real chromosome information is available.
missing_tokens <- c("0", "", ".", "NA", "na", "chr0", "chrUn")
chrom_values <- unique(df$chrom[!is.na(df$chrom)])
has_chrom <- length(setdiff(chrom_values, missing_tokens)) > 1
message(sprintf("  Distinct chrom values: %d; using chromosome colouring: %s",
                length(chrom_values), has_chrom))

if (has_chrom) {
  # Natural-sorted chromosome order (LR...924.1 < LR...925.1 ...).
  chrom_order <- str_sort(chrom_values, numeric = TRUE)
  df <- df %>% mutate(chrom = factor(.data$chrom, levels = chrom_order))

  # Per-chromosome length (max position seen) and cumulative offsets, computed
  # once over all SNPs (positions are identical across PC panels).
  chr_len <- df %>%
    filter(.data$pc == pc_levels[1]) %>%
    group_by(.data$chrom, .drop = FALSE) %>%
    summarise(len = suppressWarnings(max(.data$pos, na.rm = TRUE)), .groups = "drop") %>%
    mutate(len = ifelse(is.finite(.data$len), .data$len, 0)) %>%
    arrange(.data$chrom) %>%
    mutate(offset = lag(cumsum(as.numeric(.data$len)), default = 0))

  df <- df %>%
    left_join(chr_len %>% select(.data$chrom, .data$offset), by = "chrom") %>%
    mutate(
      gpos = .data$pos + .data$offset,
      chrom_band = factor(as.integer(.data$chrom) %% 2)
    )

  # x-axis: chromosome labels at band midpoints, dashed lines at boundaries.
  chr_axis <- chr_len %>%
    mutate(mid = .data$offset + as.numeric(.data$len) / 2,
           end = .data$offset + as.numeric(.data$len))
  boundaries <- chr_axis$end[-nrow(chr_axis)]

  p <- ggplot(df, aes(x = .data$gpos, y = .data$loading_sq, color = .data$chrom_band)) +
    geom_vline(xintercept = boundaries, linetype = "dashed",
               color = "gray70", linewidth = 0.2) +
    geom_point(size = point_size, alpha = 0.6) +
    scale_color_manual(values = c("0" = "gray30", "1" = "#3B6FB6"), guide = "none") +
    scale_x_continuous(breaks = chr_axis$mid, labels = chr_axis$chrom,
                       expand = c(0.01, 0)) +
    labs(x = "Chromosome", y = expression(Loading^2))
} else {
  # No chromosome data: single colour, plain genome-order index on x.
  df <- df %>%
    group_by(.data$pc) %>%
    arrange(.data$pos, .by_group = TRUE) %>%
    mutate(gpos = row_number()) %>%
    ungroup()

  p <- ggplot(df, aes(x = .data$gpos, y = .data$loading_sq)) +
    geom_point(size = point_size, alpha = 0.6, color = "gray30") +
    scale_x_continuous(expand = c(0.01, 0)) +
    labs(x = "SNP (genome order)", y = expression(Loading^2))
}

# Shared elements: neutral baseline, one row per PC, genome-scan-like theme.
p <- p +
  geom_hline(yintercept = neutral_baseline, linetype = "dotted",
             color = "firebrick", linewidth = 0.3) +
  facet_grid(pc ~ ., scales = "free_y") +
  theme_bw() +
  theme(
    plot.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = axis_text_size),
    axis.text.y = element_text(size = axis_text_size),
    axis.title = element_text(size = axis_title_size),
    strip.text = element_text(size = axis_title_size, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

message(sprintf("Saving plot to %s ...", output_pdf))
ggsave_pdf(output_pdf, p, width = plot_width, height = plot_height, dpi = 300)
saveRDS(p, output_rds)
message("=== DONE ===")
