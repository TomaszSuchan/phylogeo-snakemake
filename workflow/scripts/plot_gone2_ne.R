#!/usr/bin/env Rscript

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")
on.exit({
  sink(type = "message")
  sink(type = "output")
  close(log_con)
}, add = TRUE)

ggsave_utils <- tryCatch(
  file.path(dirname(normalizePath(snakemake@script)), "plot_ggsave_utils.R"),
  error = function(e) "workflow/scripts/plot_ggsave_utils.R"
)
if (file.exists(ggsave_utils)) {
  source(ggsave_utils)
} else {
  source("workflow/scripts/plot_ggsave_utils.R")
}

suppressPackageStartupMessages({
  library(ggplot2)
})

summary_file <- snakemake@input[["summary"]]
seeds_file <- snakemake@input[["seeds"]]
width <- as.numeric(snakemake@params[["width"]])
height <- as.numeric(snakemake@params[["height"]])
show_seeds <- isTRUE(as.logical(snakemake@params[["show_seed_trajectories"]]))

summary_df <- read.delim(summary_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
summary_df$generation <- as.numeric(summary_df$generation)
summary_df$ne_mean <- as.numeric(summary_df$ne_mean)
summary_df$ne_sd <- as.numeric(summary_df$ne_sd)
summary_df <- summary_df[
  is.finite(summary_df$generation) & is.finite(summary_df$ne_mean) & summary_df$ne_mean > 0,
  ,
  drop = FALSE
]
summary_df$ne_sd[!is.finite(summary_df$ne_sd)] <- 0
summary_df$ymin <- pmax(summary_df$ne_mean - summary_df$ne_sd, .Machine$double.xmin)
summary_df$ymax <- summary_df$ne_mean + summary_df$ne_sd

p <- ggplot(summary_df, aes(x = .data$generation, y = .data$ne_mean))

if (show_seeds && file.exists(seeds_file)) {
  seeds_df <- read.delim(seeds_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  seeds_df$generation <- as.numeric(seeds_df$generation)
  seeds_df$ne <- as.numeric(seeds_df$ne)
  seeds_df$seed <- as.factor(seeds_df$seed)
  seeds_df <- seeds_df[
    is.finite(seeds_df$generation) & is.finite(seeds_df$ne) & seeds_df$ne > 0,
    ,
    drop = FALSE
  ]
  p <- p +
    geom_line(
      data = seeds_df,
      aes(x = .data$generation, y = .data$ne, group = .data$seed),
      colour = "#92C5DE",
      linewidth = 0.2,
      alpha = 0.25,
      inherit.aes = FALSE
    )
}

p <- p +
  geom_ribbon(
    aes(ymin = .data$ymin, ymax = .data$ymax),
    fill = "#2166AC",
    alpha = 0.25,
    colour = NA
  ) +
  geom_line(linewidth = 0.9, colour = "#2166AC") +
  scale_y_log10() +
  scale_x_log10() +
  labs(
    x = "Generations ago",
    y = expression(N[e]),
    subtitle = "Thick line = mean across seeds; ribbon = ±1 SD; faint lines = individual seeds"
  ) +
  theme_bw(base_size = 11)

if (!show_seeds) {
  p <- p + labs(subtitle = "Thick line = mean across seeds; ribbon = ±1 SD")
}

ggsave_pdf(snakemake@output[["pdf"]], plot = p, width = width, height = height)
saveRDS(p, snakemake@output[["rds"]])
