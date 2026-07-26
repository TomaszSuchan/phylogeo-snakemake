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
currentne2_utils <- file.path(script_dir, "currentne2_parse_utils.R")
if (!file.exists(ggsave_utils)) ggsave_utils <- "workflow/scripts/plot_ggsave_utils.R"
if (!file.exists(group_utils)) group_utils <- "workflow/scripts/plot_group_utils.R"
if (!file.exists(currentne2_utils)) currentne2_utils <- "workflow/scripts/currentne2_parse_utils.R"
source(ggsave_utils)
source(group_utils)
source(currentne2_utils)

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
y_log10 <- isTRUE(as.logical(snakemake@params[["y_log10"]]))
if (length(y_log10) == 0 || is.na(y_log10)) y_log10 <- TRUE

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

ne_list <- lapply(seq_len(nrow(pop_df)), function(i) {
  out_file <- file.path(out_dir, paste0(project, ".", pop_df$pop[i], "_currentNe2_OUTPUT.txt"))
  if (!file.exists(out_file)) {
    mix <- file.path(out_dir, paste0(project, ".", pop_df$pop[i], "_currentNe2_mix_OUTPUT.txt"))
    if (file.exists(mix)) out_file <- mix
  }
  message("Reading ", out_file)
  parse_currentne2(out_file, pop_df$population[i])
})
ne_df <- do.call(rbind, ne_list)

levels_order <- group_levels(ne_df, "population", group_sort_by(group_sort_by_param))
ne_df$population <- factor(ne_df$population, levels = levels_order)
ne_df$ci50_low[!is.finite(ne_df$ci50_low)] <- ne_df$ne[!is.finite(ne_df$ci50_low)]
ne_df$ci50_high[!is.finite(ne_df$ci50_high)] <- ne_df$ne[!is.finite(ne_df$ci50_high)]
ne_df$ci90_low[!is.finite(ne_df$ci90_low)] <- ne_df$ne[!is.finite(ne_df$ci90_low)]
ne_df$ci90_high[!is.finite(ne_df$ci90_high)] <- ne_df$ne[!is.finite(ne_df$ci90_high)]

palette_vals <- group_fill_values(group_colors_param)

p <- ggplot(ne_df, aes(x = .data$population, y = .data$ne, fill = .data$population)) +
  # Outer whiskers: 90% CI (thinner)
  geom_errorbar(
    aes(ymin = .data$ci90_low, ymax = .data$ci90_high, colour = .data$population),
    width = 0.18,
    linewidth = 0.45,
    na.rm = TRUE
  ) +
  # Inner whiskers: 50% CI (normal)
  geom_errorbar(
    aes(ymin = .data$ci50_low, ymax = .data$ci50_high, colour = .data$population),
    width = 0.12,
    linewidth = 0.85,
    na.rm = TRUE
  ) +
  geom_point(
    shape = 21,
    size = 3.2,
    colour = "black",
    stroke = 0.35
  ) +
  labs(
    x = legend_title,
    y = expression(N[e])
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 35, hjust = 1),
    panel.grid.minor = element_blank()
  )

if (y_log10) {
  p <- p + scale_y_log10()
}

if (!is.null(palette_vals)) {
  p <- p +
    scale_fill_manual(values = palette_vals, drop = FALSE) +
    scale_colour_manual(values = palette_vals, drop = FALSE)
}

ggsave_pdf(snakemake@output[["pdf"]], plot = p, width = width, height = height)
saveRDS(p, snakemake@output[["rds"]])
message("Wrote ", snakemake@output[["pdf"]])
