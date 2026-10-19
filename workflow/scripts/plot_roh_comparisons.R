#!/usr/bin/env Rscript
# Mean ± SE and Tukey HSD heatmaps from the ROH comparison text report.

suppressPackageStartupMessages(library(ggplot2))

ggsave_utils <- tryCatch(
  file.path(dirname(normalizePath(snakemake@script)), "plot_ggsave_utils.R"),
  error = function(e) "workflow/scripts/plot_ggsave_utils.R"
)
if (file.exists(ggsave_utils)) {
  source(ggsave_utils)
} else {
  source("workflow/scripts/plot_ggsave_utils.R")
}

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")
on.exit({
  while (sink.number(type = "message") > 0) sink(type = "message")
  while (sink.number(type = "output") > 0) sink(type = "output")
  close(log_con)
}, add = TRUE)

group_col <- snakemake@params[["group_col"]]
axis_title_size <- as.numeric(snakemake@params[["axis_title_size"]])
axis_text_size <- as.numeric(snakemake@params[["axis_text_size"]])
plot_width <- as.numeric(snakemake@params[["width"]])
plot_height <- as.numeric(snakemake@params[["height"]])
plot_dpi <- as.numeric(snakemake@params[["dpi"]])
if (is.na(axis_title_size)) axis_title_size <- 10
if (is.na(axis_text_size)) axis_text_size <- 8
if (is.na(plot_width)) plot_width <- 8
if (is.na(plot_height)) plot_height <- 6
if (is.na(plot_dpi)) plot_dpi <- 300

metrics <- list(
  F_ROH = list(
    stats = "### F_ROH Statistics ###",
    tukey = "### F_ROH Multiple Comparisons (Tukey HSD) ###",
    ylab = "Mean F_ROH"
  ),
  N_ROH = list(
    stats = "### Number of ROH Segments Statistics ###",
    tukey = "### Number of ROH Segments Multiple Comparisons (Tukey HSD) ###",
    ylab = "Mean number of ROH segments"
  )
)

section_lines <- function(lines, start_idx, end_idx) {
  if (start_idx >= end_idx) return(character())
  lines[(start_idx + 1):(end_idx - 1)]
}

find_heading <- function(lines, heading, from = 1) {
  hit <- which(trimws(lines) == heading & seq_along(lines) >= from)
  if (length(hit) == 0) NA_integer_ else hit[1]
}

parse_summary <- function(block) {
  rows <- grep("^[^\t#].*\t", block, value = TRUE)
  rows <- rows[!grepl("^Group[[:space:]]", rows)]
  if (length(rows) == 0) {
    return(data.frame(group = character(), n = integer(), mean = numeric(), sd = numeric()))
  }
  parsed <- lapply(strsplit(rows, "\t", fixed = TRUE), function(parts) {
    parts <- trimws(parts)
    parts <- parts[nzchar(parts)]
    if (length(parts) < 4) return(NULL)
    data.frame(
      group = parts[1],
      n = as.integer(parts[2]),
      mean = as.numeric(parts[3]),
      sd = as.numeric(parts[4]),
      stringsAsFactors = FALSE
    )
  })
  parsed <- do.call(rbind, parsed)
  parsed$se <- parsed$sd / sqrt(parsed$n)
  parsed
}

trailing_numbers <- function(text) {
  hits <- regmatches(text, gregexpr("-?[0-9]+(?:\\.[0-9]+)?(?:[eE][+-]?[0-9]+)?", text))[[1]]
  if (length(hits) == 0) return(numeric())
  as.numeric(hits)
}

split_contrast <- function(label, groups) {
  groups <- groups[order(nchar(groups), decreasing = TRUE)]
  for (group in groups) {
    suffix <- paste0("-", group)
    if (endsWith(label, suffix)) {
      left <- substring(label, 1, nchar(label) - nchar(suffix))
      if (left %in% groups) return(c(left, group))
    }
  }
  c(NA_character_, NA_character_)
}

parse_tukey <- function(block, groups) {
  body <- block[!grepl("^\\s*$|^#|^---|^Fit:|^\\s*Tukey |^\\s*95% |^\\$`", block)]
  body <- body[!grepl("p adj$", trimws(body)) & trimws(body) != "p adj"]
  body <- body[!grepl("^\\s*diff\\s", body)]
  contrasts <- list()
  pending_p <- character()
  for (line in body) {
    nums <- trailing_numbers(line)
    if (length(nums) < 1) next
    label <- sub("\\s+-?[0-9].*$", "", line)
    label <- trimws(label)
    if (!nzchar(label)) next
    if (length(nums) >= 4) {
      contrasts[[length(contrasts) + 1]] <- data.frame(
        label = label, diff = nums[1], lwr = nums[2], upr = nums[3], p = nums[4],
        stringsAsFactors = FALSE
      )
    } else if (length(nums) == 3) {
      contrasts[[length(contrasts) + 1]] <- data.frame(
        label = label, diff = nums[1], lwr = nums[2], upr = nums[3], p = NA_real_,
        stringsAsFactors = FALSE
      )
    } else if (length(nums) == 1) {
      pending_p <- c(pending_p, nums[1])
      # store label alignment later; p-only rows follow the same order
      attr(pending_p, "labels") <- c(attr(pending_p, "labels"), label)
    }
  }
  if (length(contrasts) == 0) {
    return(data.frame(group_a = character(), group_b = character(), diff = numeric(), p = numeric()))
  }
  tukey <- do.call(rbind, contrasts)
  p_labels <- attr(pending_p, "labels")
  if (anyNA(tukey$p) && length(pending_p) == nrow(tukey)) {
    tukey$p <- as.numeric(pending_p)
  }
  split_groups <- lapply(tukey$label, split_contrast, groups = groups)
  tukey$group_b <- vapply(split_groups, `[`, character(1), 1)
  tukey$group_a <- vapply(split_groups, `[`, character(1), 2)
  tukey <- tukey[!is.na(tukey$group_a) & !is.na(tukey$group_b), c("group_a", "group_b", "diff", "p")]
  tukey
}

mean_plot <- function(summary_df, ylab) {
  summary_df$group <- factor(summary_df$group, levels = summary_df$group[order(summary_df$mean)])
  ggplot(summary_df, aes(x = mean, y = group)) +
    geom_errorbarh(aes(xmin = mean - se, xmax = mean + se), height = 0.2, linewidth = 0.3) +
    geom_point(size = 1.8) +
    labs(x = ylab, y = NULL, caption = "Point is the group mean. Whiskers are ± 1 SE.") +
    theme_bw() +
    theme(
      axis.title = element_text(size = axis_title_size),
      axis.text = element_text(size = axis_text_size)
    )
}

tukey_plot <- function(summary_df, tukey_df) {
  groups <- summary_df$group[order(summary_df$mean)]
  grid <- expand.grid(group_a = groups, group_b = groups, stringsAsFactors = FALSE)
  direct <- match(paste(grid$group_a, grid$group_b), paste(tukey_df$group_a, tukey_df$group_b))
  grid$diff <- tukey_df$diff[direct]
  grid$p <- tukey_df$p[direct]
  reverse <- match(paste(grid$group_b, grid$group_a), paste(tukey_df$group_a, tukey_df$group_b))
  missing <- is.na(grid$diff) & !is.na(reverse)
  grid$diff[missing] <- -tukey_df$diff[reverse[missing]]
  grid$p[missing] <- tukey_df$p[reverse[missing]]
  same <- grid$group_a == grid$group_b
  grid$diff[same] <- 0
  grid$p[same] <- 1
  grid$group_a <- factor(grid$group_a, levels = groups)
  grid$group_b <- factor(grid$group_b, levels = groups)
  grid$star <- ifelse(!is.na(grid$p) & grid$p < 0.05, "*", "")
  ggplot(grid, aes(x = group_a, y = group_b, fill = diff)) +
    geom_tile(color = "white", linewidth = 0.2) +
    geom_text(aes(label = star), size = 2.5, na.rm = TRUE) +
    scale_fill_gradient2(na.value = "grey92", name = "Difference") +
    labs(
      x = NULL,
      y = NULL,
      caption = "Tukey HSD difference of means. A star marks an adjusted p below 0.05."
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = axis_text_size),
      axis.text.y = element_text(size = axis_text_size),
      axis.title = element_text(size = axis_title_size)
    )
}

lines <- readLines(snakemake@input[["stats"]], warn = FALSE)
group_starts <- grep("^STATISTICS FOR GROUPING COLUMN:", lines)
group_names <- sub("^STATISTICS FOR GROUPING COLUMN:\\s*", "", lines[group_starts])
this_group <- which(group_names == group_col)
if (length(this_group) != 1) {
  stop("Grouping column not found in the comparison report: ", group_col)
}
block_start <- group_starts[this_group]
block_end <- if (this_group < length(group_starts)) group_starts[this_group + 1] else length(lines) + 1
block <- section_lines(lines, block_start - 1, block_end)

plots <- list()
n_groups <- NA_integer_
for (metric_name in names(metrics)) {
  spec <- metrics[[metric_name]]
  stats_at <- find_heading(block, spec$stats)
  tukey_at <- find_heading(block, spec$tukey)
  next_heading <- grep("^### ", block)
  next_heading <- next_heading[next_heading > stats_at]
  stats_end <- if (length(next_heading)) next_heading[1] else length(block) + 1
  tukey_next <- next_heading[next_heading > tukey_at]
  tukey_end <- if (length(tukey_next)) tukey_next[1] else length(block) + 1
  if (is.na(stats_at) || is.na(tukey_at)) {
    stop("Missing ", metric_name, " section for ", group_col)
  }
  summary_df <- parse_summary(section_lines(block, stats_at, stats_end))
  tukey_df <- parse_tukey(section_lines(block, tukey_at, tukey_end), summary_df$group)
  if (is.na(n_groups)) n_groups <- nrow(summary_df)
  cat(metric_name, "groups:", nrow(summary_df), "pairwise tests:", nrow(tukey_df), "\n")
  plots[[paste0(metric_name, "_mean")]] <- mean_plot(summary_df, spec$ylab)
  plots[[paste0(metric_name, "_tukey")]] <- tukey_plot(summary_df, tukey_df)
}

panel <- max(plot_width, plot_height, 0.32 * n_groups + 2)

dir.create(dirname(snakemake@output[["pdf"]]), recursive = TRUE, showWarnings = FALSE)
open_pdf <- cairo_pdf_device()
if (is.null(open_pdf)) {
  grDevices::pdf(snakemake@output[["pdf"]], width = panel, height = panel, onefile = TRUE)
} else {
  open_pdf(snakemake@output[["pdf"]], width = panel, height = panel, onefile = TRUE)
}
for (plot_obj in plots) print(plot_obj)
grDevices::dev.off()
saveRDS(plots, snakemake@output[["rds"]])
cat("Wrote", snakemake@output[["pdf"]], "\n")
