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

lines <- readLines(ne_file, warn = FALSE)
data_lines <- lines[!grepl("^#", lines)]
df <- read.delim(
  text = paste(data_lines, collapse = "\n"),
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)
gen_col <- grep("^Generation$", names(df), ignore.case = TRUE, value = TRUE)[1]
ne_col <- grep("^Ne_", names(df), value = TRUE)[1]
ne_df <- data.frame(
  generation = as.numeric(df[[gen_col]]),
  ne = as.numeric(df[[ne_col]]),
  stringsAsFactors = FALSE
)
ne_df <- ne_df[is.finite(ne_df$generation) & is.finite(ne_df$ne) & ne_df$ne > 0, , drop = FALSE]

p <- ggplot(ne_df, aes(x = .data$generation, y = .data$ne)) +
  geom_line(linewidth = 0.6, colour = "#2166AC") +
  geom_point(size = 1.2, colour = "#2166AC") +
  scale_y_log10() +
  scale_x_log10() +
  labs(x = "Generations ago", y = expression(N[e])) +
  theme_bw(base_size = 11)

ggsave_pdf(snakemake@output[["pdf"]], plot = p, width = width, height = height)
saveRDS(p, snakemake@output[["rds"]])
