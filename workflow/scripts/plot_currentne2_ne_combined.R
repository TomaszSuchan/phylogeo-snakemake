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

parse_currentne2 <- function(path, population) {
  if (!file.exists(path)) {
    stop("Missing currentNe2 output: ", path)
  }
  lines <- readLines(path, warn = FALSE)
  get_after <- function(pattern) {
    idx <- grep(pattern, lines, perl = TRUE)
    if (length(idx) == 0) return(NA_real_)
    for (i in idx) {
      if (i >= length(lines)) next
      val <- suppressWarnings(as.numeric(trimws(lines[[i + 1]])))
      if (is.finite(val)) return(val)
    }
    NA_real_
  }
  wg_start <- grep("integration over the whole genome", lines, fixed = TRUE)
  bc_start <- grep("LD between chromosomes", lines, fixed = TRUE)
  if (length(wg_start) > 0) {
    end <- if (length(bc_start) > 0 && bc_start[1] > wg_start[1]) bc_start[1] - 1 else length(lines)
    lines <- lines[wg_start[1]:end]
  } else if (length(bc_start) > 0) {
    lines <- lines[bc_start[1]:length(lines)]
  }
  out <- data.frame(
    population = population,
    ne = get_after("^# Ne point estimate:"),
    ci50_low = get_after("^# Lower limit of 50% CI:"),
    ci50_high = get_after("^# Upper (bound|limit) of 50% CI:"),
    ci90_low = get_after("^# Lower limit of 90% CI:"),
    ci90_high = get_after("^# Upper limit of 90% CI:"),
    stringsAsFactors = FALSE
  )
  if (!is.finite(out$ne[1])) {
    stop("Could not parse Ne point estimate from ", path)
  }
  out
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

palette_vals <- group_fill_values(group_colors_param)

p <- ggplot(ne_df, aes(x = .data$population, y = .data$ne, colour = .data$population)) +
  geom_errorbar(
    aes(ymin = .data$ci90_low, ymax = .data$ci90_high),
    width = 0.2,
    linewidth = 0.6,
    na.rm = TRUE
  ) +
  geom_point(size = 2.8) +
  scale_y_log10() +
  labs(
    x = legend_title,
    y = expression(N[e]),
    colour = legend_title,
    subtitle = "Points = Ne point estimate; error bars = 90% CI"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 35, hjust = 1)
  )

if (!is.null(palette_vals)) {
  p <- p + scale_colour_manual(values = palette_vals, drop = FALSE)
}

ggsave_pdf(snakemake@output[["pdf"]], plot = p, width = width, height = height)
saveRDS(p, snakemake@output[["rds"]])
message("Wrote ", snakemake@output[["pdf"]])
