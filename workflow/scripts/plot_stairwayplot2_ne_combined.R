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

read_ne <- function(ne_file, population) {
  if (!file.exists(ne_file)) {
    stop("Missing Stairway Plot 2 Ne file: ", ne_file)
  }
  df <- read.delim(ne_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  out <- data.frame(
    population = population,
    generation = as.numeric(df$generation),
    ne_median = as.numeric(df$ne_median),
    ne_lower_95 = as.numeric(df$ne_lower_95),
    ne_upper_95 = as.numeric(df$ne_upper_95),
    stringsAsFactors = FALSE
  )
  out <- out[is.finite(out$generation) & is.finite(out$ne_median) & out$ne_median > 0, , drop = FALSE]
  out$ne_lower_95[!is.finite(out$ne_lower_95)] <- out$ne_median[!is.finite(out$ne_lower_95)]
  out$ne_upper_95[!is.finite(out$ne_upper_95)] <- out$ne_median[!is.finite(out$ne_upper_95)]
  out$ymin <- pmax(out$ne_lower_95, .Machine$double.xmin)
  out$ymax <- out$ne_upper_95
  out
}

summary_list <- lapply(seq_len(nrow(pop_df)), function(i) {
  ne_file <- file.path(out_dir, paste0(project, ".", pop_df$pop[i], ".ne.tsv"))
  message("Reading ", ne_file)
  read_ne(ne_file, pop_df$population[i])
})
summary_df <- do.call(rbind, summary_list)

levels_order <- group_levels(summary_df, "population", group_sort_by(group_sort_by_param))
summary_df$population <- factor(summary_df$population, levels = levels_order)
palette_vals <- group_fill_values(group_colors_param)

p <- ggplot(
  summary_df,
  aes(x = .data$generation, y = .data$ne_median, colour = .data$population, fill = .data$population)
) +
  geom_ribbon(aes(ymin = .data$ymin, ymax = .data$ymax), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.8) +
  labs(
    x = "Generations ago",
    y = expression(N[e]),
    colour = legend_title,
    fill = legend_title
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "right")

if (!is.null(palette_vals)) {
  p <- p +
    scale_colour_manual(values = palette_vals, drop = FALSE) +
    scale_fill_manual(values = palette_vals, drop = FALSE)
}

p_log <- p + scale_y_log10() + scale_x_log10()
p_linear <- p
p_xlinear_ylog <- p + scale_y_log10()

ggsave_pdf(snakemake@output[["pdf"]], plot = p_log, width = width, height = height)
saveRDS(p_log, snakemake@output[["rds"]])
ggsave_pdf(snakemake@output[["pdf_linear"]], plot = p_linear, width = width, height = height)
saveRDS(p_linear, snakemake@output[["rds_linear"]])
ggsave_pdf(snakemake@output[["pdf_xlinear_ylog"]], plot = p_xlinear_ylog, width = width, height = height)
saveRDS(p_xlinear_ylog, snakemake@output[["rds_xlinear_ylog"]])
message(
  "Wrote ", snakemake@output[["pdf"]], ", ",
  snakemake@output[["pdf_linear"]], ", and ",
  snakemake@output[["pdf_xlinear_ylog"]]
)
