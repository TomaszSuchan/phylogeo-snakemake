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

out_file <- snakemake@input[["output"]]
width <- as.numeric(snakemake@params[["width"]])
height <- as.numeric(snakemake@params[["height"]])
stratum <- snakemake@wildcards[["stratum"]]

parse_currentne2 <- function(path) {
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
  # Prefer whole-genome block when present.
  wg_start <- grep("integration over the whole genome", lines, fixed = TRUE)
  bc_start <- grep("LD between chromosomes", lines, fixed = TRUE)
  if (length(wg_start) > 0) {
    end <- if (length(bc_start) > 0 && bc_start[1] > wg_start[1]) bc_start[1] - 1 else length(lines)
    lines <- lines[wg_start[1]:end]
    estimate_type <- "whole_genome"
  } else if (length(bc_start) > 0) {
    lines <- lines[bc_start[1]:length(lines)]
    estimate_type <- "between_chromosomes"
  } else {
    estimate_type <- "primary"
  }
  data.frame(
    estimate_type = estimate_type,
    ne = get_after("^# Ne point estimate:"),
    ci50_low = get_after("^# Lower limit of 50% CI:"),
    ci50_high = get_after("^# Upper (bound|limit) of 50% CI:"),
    ci90_low = get_after("^# Lower limit of 90% CI:"),
    ci90_high = get_after("^# Upper limit of 90% CI:"),
    stringsAsFactors = FALSE
  )
}

df <- parse_currentne2(out_file)
if (!is.finite(df$ne[1])) {
  stop("Could not parse Ne point estimate from ", out_file)
}
df$population <- stratum

p <- ggplot(df, aes(x = .data$population, y = .data$ne)) +
  geom_errorbar(
    aes(ymin = .data$ci90_low, ymax = .data$ci90_high),
    width = 0.15,
    colour = "#92C5DE",
    linewidth = 0.7,
    na.rm = TRUE
  ) +
  geom_errorbar(
    aes(ymin = .data$ci50_low, ymax = .data$ci50_high),
    width = 0.08,
    colour = "#2166AC",
    linewidth = 1.0,
    na.rm = TRUE
  ) +
  geom_point(size = 2.5, colour = "#2166AC") +
  scale_y_log10() +
  labs(
    x = NULL,
    y = expression(N[e]),
    subtitle = paste0("Estimate: ", df$estimate_type[1], "; thick = 50% CI, thin = 90% CI")
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave_pdf(snakemake@output[["pdf"]], plot = p, width = width, height = height)
saveRDS(p, snakemake@output[["rds"]])
