#!/usr/bin/env Rscript

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")
on.exit({
  sink(type = "message")
  sink(type = "output")
  close(log_con)
}, add = TRUE)

script_dir <- tryCatch(
  dirname(normalizePath(snakemake@script)),
  error = function(e) "workflow/scripts"
)
ggsave_utils <- file.path(script_dir, "plot_ggsave_utils.R")
group_utils <- file.path(script_dir, "plot_group_utils.R")
if (!file.exists(ggsave_utils)) ggsave_utils <- "workflow/scripts/plot_ggsave_utils.R"
if (!file.exists(group_utils)) group_utils <- "workflow/scripts/plot_group_utils.R"
source(ggsave_utils)
source(group_utils)

suppressPackageStartupMessages({
  library(ggplot2)
})

populations_file <- snakemake@input[["populations"]]
out_dir <- snakemake@params[["out_dir"]]
project <- snakemake@wildcards[["project"]]
width <- as.numeric(snakemake@params[["width"]])
height <- as.numeric(snakemake@params[["height"]])
group_colors_param <- snakemake@params[["group_colors"]]
group_sort_by_param <- snakemake@params[["group_sort_by"]]
show_seeds <- isTRUE(as.logical(snakemake@params[["show_seed_trajectories"]]))
legend_title <- as.character(snakemake@params[["legend_title"]])
if (is.na(legend_title) || !nzchar(legend_title)) {
  legend_title <- "Population"
}

pop_df <- read.delim(
  populations_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)
if (!all(c("population", "pop") %in% names(pop_df))) {
  stop("populations file must contain columns 'population' and 'pop'")
}

read_summary <- function(summary_file, population) {
  if (!file.exists(summary_file)) {
    stop("Missing GONE2 summary file: ", summary_file)
  }
  df <- read.delim(summary_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  out <- data.frame(
    population = population,
    generation = as.numeric(df$generation),
    ne_mean = as.numeric(df$ne_mean),
    ne_sd = as.numeric(df$ne_sd),
    stringsAsFactors = FALSE
  )
  out$ne_sd[!is.finite(out$ne_sd)] <- 0
  out <- out[is.finite(out$generation) & is.finite(out$ne_mean) & out$ne_mean > 0, , drop = FALSE]
  out$ymin <- pmax(out$ne_mean - out$ne_sd, .Machine$double.xmin)
  out$ymax <- out$ne_mean + out$ne_sd
  out
}

read_seeds <- function(seeds_file, population) {
  if (!file.exists(seeds_file)) {
    return(NULL)
  }
  df <- read.delim(seeds_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  out <- data.frame(
    population = population,
    seed = as.factor(df$seed),
    generation = as.numeric(df$generation),
    ne = as.numeric(df$ne),
    stringsAsFactors = FALSE
  )
  out[is.finite(out$generation) & is.finite(out$ne) & out$ne > 0, , drop = FALSE]
}

summary_list <- lapply(seq_len(nrow(pop_df)), function(i) {
  summary_file <- file.path(out_dir, paste0(project, ".", pop_df$pop[i], "_GONE2_Ne_summary.tsv"))
  message("Reading ", summary_file)
  read_summary(summary_file, pop_df$population[i])
})
summary_df <- do.call(rbind, summary_list)

levels_order <- group_levels(summary_df, "population", group_sort_by(group_sort_by_param))
summary_df$population <- factor(summary_df$population, levels = levels_order)

palette_vals <- group_fill_values(group_colors_param)

p <- ggplot(summary_df, aes(x = .data$generation, y = .data$ne_mean, colour = .data$population, fill = .data$population))

if (show_seeds) {
  seeds_list <- lapply(seq_len(nrow(pop_df)), function(i) {
    seeds_file <- file.path(out_dir, paste0(project, ".", pop_df$pop[i], "_GONE2_Ne_seeds.tsv"))
    read_seeds(seeds_file, pop_df$population[i])
  })
  seeds_df <- do.call(rbind, Filter(Negate(is.null), seeds_list))
  if (!is.null(seeds_df) && nrow(seeds_df) > 0) {
    seeds_df$population <- factor(seeds_df$population, levels = levels_order)
    p <- p +
      geom_line(
        data = seeds_df,
        aes(x = .data$generation, y = .data$ne, group = interaction(.data$population, .data$seed)),
        linewidth = 0.15,
        alpha = 0.15,
        show.legend = FALSE
      )
  }
}

p <- p +
  geom_ribbon(aes(ymin = .data$ymin, ymax = .data$ymax), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 0.8) +
  scale_y_log10() +
  scale_x_log10() +
  labs(
    x = "Generations ago",
    y = expression(N[e]),
    colour = legend_title,
    fill = legend_title,
    subtitle = if (show_seeds) {
      "Thick lines = mean across seeds; ribbons = ±1 SD; faint lines = individual seeds"
    } else {
      "Thick lines = mean across seeds; ribbons = ±1 SD"
    }
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "right")

if (!is.null(palette_vals)) {
  p <- p +
    scale_colour_manual(values = palette_vals, drop = FALSE) +
    scale_fill_manual(values = palette_vals, drop = FALSE)
}

ggsave_pdf(snakemake@output[["pdf"]], plot = p, width = width, height = height)
saveRDS(p, snakemake@output[["rds"]])
message("Wrote ", snakemake@output[["pdf"]])
