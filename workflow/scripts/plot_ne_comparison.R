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

GONE2_LABEL <- "GONE2"
STAIRWAYPLOT2_LABEL <- "Stairway Plot 2"
CURRENTNE2_LABEL <- "CurrentNe2"
METHOD_COLOURS <- c("#D62728", "#1F77B4")
names(METHOD_COLOURS) <- c(GONE2_LABEL, STAIRWAYPLOT2_LABEL)
CURRENTNE2_COLOUR <- "#2CA02C"

populations_file <- snakemake@input[["populations"]]
gone2_files <- as.character(snakemake@input[["gone2"]])
stairwayplot2_files <- as.character(snakemake@input[["stairwayplot2"]])
currentne2_files <- if ("currentne2" %in% names(snakemake@input)) {
  as.character(snakemake@input[["currentne2"]])
} else {
  character(0)
}
pops <- as.character(unlist(snakemake@params[["pops"]], use.names = FALSE))
width <- as.numeric(snakemake@params[["width"]])
height <- as.numeric(snakemake@params[["height"]])
group_sort_by_param <- snakemake@params[["group_sort_by"]]

if (length(pops) == 0) {
  stop("No population has both a GONE2 and a Stairway Plot 2 Ne trajectory")
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
labels <- pop_df$population[match(pops, pop_df$pop)]
labels[is.na(labels)] <- pops[is.na(labels)]

# GONE2 reports one Ne per discrete generation; the seed spread is not drawn here
# because the Stairway Plot ribbon already carries the uncertainty comparison.
read_gone2 <- function(path, population) {
  message("Reading ", path)
  df <- read.delim(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  out <- data.frame(
    population = population,
    method = GONE2_LABEL,
    generation = as.numeric(df$generation),
    ne = as.numeric(df$ne_mean),
    ymin = NA_real_,
    ymax = NA_real_,
    stringsAsFactors = FALSE
  )
  out[is.finite(out$generation) & is.finite(out$ne) & out$ne > 0, , drop = FALSE]
}

read_stairwayplot2 <- function(path, population) {
  message("Reading ", path)
  df <- read.delim(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  out <- data.frame(
    population = population,
    method = STAIRWAYPLOT2_LABEL,
    generation = as.numeric(df$generation),
    ne = as.numeric(df$ne_median),
    ymin = as.numeric(df$ne_lower_95),
    ymax = as.numeric(df$ne_upper_95),
    stringsAsFactors = FALSE
  )
  out <- out[is.finite(out$generation) & is.finite(out$ne) & out$ne > 0, , drop = FALSE]
  out$ymin[!is.finite(out$ymin)] <- out$ne[!is.finite(out$ymin)]
  out$ymax[!is.finite(out$ymax)] <- out$ne[!is.finite(out$ymax)]
  out$ymin <- pmax(out$ymin, .Machine$double.xmin)
  out
}

gone2_df <- do.call(rbind, Map(read_gone2, gone2_files, labels))
stairwayplot2_df <- do.call(rbind, Map(read_stairwayplot2, stairwayplot2_files, labels))
traj_df <- rbind(gone2_df, stairwayplot2_df)

# The two methods resolve different epochs; shading the GONE2 span makes the
# extent of the comparable window explicit in every facet.
gone2_window <- do.call(rbind, lapply(split(gone2_df, gone2_df$population), function(d) {
  data.frame(
    population = d$population[1],
    xmin = min(d$generation),
    xmax = max(d$generation),
    stringsAsFactors = FALSE
  )
}))

levels_order <- group_levels(traj_df, "population", group_sort_by(group_sort_by_param))
traj_df$population <- factor(traj_df$population, levels = levels_order)
gone2_window$population <- factor(gone2_window$population, levels = levels_order)
traj_df$method <- factor(traj_df$method, levels = names(METHOD_COLOURS))

p <- ggplot(
  traj_df,
  aes(x = .data$generation, y = .data$ne, colour = .data$method, fill = .data$method)
) +
  geom_rect(
    data = gone2_window,
    mapping = aes(xmin = .data$xmin, xmax = .data$xmax, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    fill = "grey50",
    alpha = 0.12
  ) +
  geom_ribbon(aes(ymin = .data$ymin, ymax = .data$ymax), alpha = 0.15, colour = NA, na.rm = TRUE) +
  geom_line(linewidth = 0.7)

if (length(currentne2_files) > 0) {
  currentne2_df <- do.call(rbind, Map(function(path, population) {
    message("Reading ", path)
    parse_currentne2(path, population)[, c("population", "ne"), drop = FALSE]
  }, currentne2_files, labels))
  currentne2_df$population <- factor(currentne2_df$population, levels = levels_order)
  # Carried as a column rather than referenced inside aes() so the saved RDS stays
  # self-contained when reopened outside this script.
  currentne2_df$reference <- CURRENTNE2_LABEL
  p <- p +
    geom_hline(
      data = currentne2_df,
      mapping = aes(yintercept = .data$ne, linetype = .data$reference),
      inherit.aes = FALSE,
      colour = CURRENTNE2_COLOUR,
      linewidth = 0.5
    ) +
    scale_linetype_manual(name = NULL, values = setNames("dashed", CURRENTNE2_LABEL))
}

p <- p +
  facet_wrap(~population) +
  scale_x_log10() +
  scale_y_log10() +
  scale_colour_manual(values = METHOD_COLOURS, drop = FALSE) +
  scale_fill_manual(values = METHOD_COLOURS, drop = FALSE) +
  labs(
    x = "Generations ago",
    y = expression(N[e]),
    colour = "Method",
    fill = "Method"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

ggsave_pdf(snakemake@output[["pdf"]], plot = p, width = width, height = height)
saveRDS(p, snakemake@output[["rds"]])
message("Wrote ", snakemake@output[["pdf"]])
