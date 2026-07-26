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

ne_file <- snakemake@input[["ne"]]
width <- as.numeric(snakemake@params[["width"]])
height <- as.numeric(snakemake@params[["height"]])

df <- read.delim(ne_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
df$generation <- as.numeric(df$generation)
df$ne_median <- as.numeric(df$ne_median)
df$ne_lower_95 <- as.numeric(df$ne_lower_95)
df$ne_upper_95 <- as.numeric(df$ne_upper_95)
df$ne_lower_75 <- as.numeric(df$ne_lower_75)
df$ne_upper_75 <- as.numeric(df$ne_upper_75)
df <- df[
  is.finite(df$generation) & is.finite(df$ne_median) & df$ne_median > 0,
  ,
  drop = FALSE
]
df$ne_lower_95[!is.finite(df$ne_lower_95)] <- df$ne_median[!is.finite(df$ne_lower_95)]
df$ne_upper_95[!is.finite(df$ne_upper_95)] <- df$ne_median[!is.finite(df$ne_upper_95)]
df$ne_lower_75[!is.finite(df$ne_lower_75)] <- df$ne_median[!is.finite(df$ne_lower_75)]
df$ne_upper_75[!is.finite(df$ne_upper_75)] <- df$ne_median[!is.finite(df$ne_upper_75)]

p <- ggplot(df, aes(x = .data$generation, y = .data$ne_median)) +
  geom_ribbon(
    aes(ymin = pmax(.data$ne_lower_95, .Machine$double.xmin), ymax = .data$ne_upper_95),
    fill = "#BDBDBD",
    alpha = 0.45,
    colour = NA
  ) +
  geom_ribbon(
    aes(ymin = pmax(.data$ne_lower_75, .Machine$double.xmin), ymax = .data$ne_upper_75),
    fill = "#636363",
    alpha = 0.35,
    colour = NA
  ) +
  geom_line(linewidth = 0.9, colour = "#B2182B") +
  labs(
    x = "Generations ago",
    y = expression(N[e])
  ) +
  theme_bw(base_size = 11)

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
