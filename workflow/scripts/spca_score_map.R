#!/usr/bin/env Rscript
# Map retained sPCA scores on the shared mapmixture basemap.

log_file <- snakemake@log[[1]]
dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")
on.exit({
  while (sink.number(type = "message") > 0) sink(type = "message")
  while (sink.number(type = "output") > 0) sink(type = "output")
  close(log_con)
}, add = TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(mapmixture)
  library(terra)
})

common_functions <- tryCatch({
  script_dir <- dirname(normalizePath(snakemake@script))
  file.path(script_dir, "common_map_functions.R")
}, error = function(e) "workflow/scripts/common_map_functions.R")
if (file.exists(common_functions)) {
  source(common_functions)
} else {
  source("workflow/scripts/common_map_functions.R")
}

`%||%` <- function(x, y) if (is.null(x)) y else x

unwrap_spca <- function(path) {
  obj <- readRDS(path)
  if (inherits(obj, "spca")) {
    return(list(spca = obj, nfposi = obj$nfposi, nfnega = obj$nfnega))
  }
  if (is.list(obj) && !is.null(obj$spca)) {
    return(list(
      spca = obj$spca,
      nfposi = obj$nfposi %||% obj$spca$nfposi,
      nfnega = obj$nfnega %||% obj$spca$nfnega
    ))
  }
  stop("Unrecognized sPCA results RDS format: ", path)
}

score_type <- as.character(snakemake@params[["score_type"]])
axis <- as.integer(snakemake@params[["axis"]])

message("=== sPCA score map ===")
message("score_type: ", score_type)
message("axis: ", axis)

res <- unwrap_spca(snakemake@input[["results_rds"]])
spca_obj <- res$spca
nfposi <- as.integer(res$nfposi)
nfnega <- as.integer(res$nfnega)

if (score_type == "global") {
  if (axis > nfposi) {
    stop("Requested global axis ", axis, " but nfposi is ", nfposi)
  }
  spca_col <- axis
  score_label <- paste0("Global axis ", axis)
} else if (score_type == "local") {
  if (axis > nfnega) {
    stop("Requested local axis ", axis, " but nfnega is ", nfnega)
  }
  spca_col <- nfposi + axis
  score_label <- paste0("Local axis ", axis)
} else {
  stop("score_type must be 'global' or 'local'; got ", score_type)
}

scores <- as.data.frame(spca_obj$li)
if (spca_col > ncol(scores)) {
  stop("Requested sPCA column ", spca_col, " but only ", ncol(scores), " score columns exist")
}
score_df <- data.frame(
  Ind = rownames(scores),
  spca_score = as.numeric(scores[[spca_col]]),
  stringsAsFactors = FALSE
)

indpopdata <- read.table(
  snakemake@input[["indpopdata"]],
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)
indpopdata <- coerce_indpopdata_lat_lon(indpopdata)
required_cols <- c("Ind", "Site", "Lat", "Lon")
if (!all(required_cols %in% names(indpopdata))) {
  stop(
    "indpopdata must contain columns: ",
    paste(required_cols, collapse = ", "),
    ". Found: ",
    paste(names(indpopdata), collapse = ", ")
  )
}

plot_df <- inner_join(indpopdata, score_df, by = "Ind")
if (nrow(plot_df) == 0) {
  stop("No overlapping individuals between sPCA scores and indpopdata")
}
missing_scores <- setdiff(score_df$Ind, plot_df$Ind)
if (length(missing_scores) > 0) {
  warning(length(missing_scores), " sPCA individuals were missing from indpopdata and not mapped")
}

coords_df <- indpopdata %>%
  select(Site, Lat, Lon) %>%
  distinct(Site, .keep_all = TRUE)

width <- as.numeric(snakemake@params[["width"]])
height <- as.numeric(snakemake@params[["height"]])
dpi <- as.numeric(snakemake@params[["dpi"]])
crs <- as.numeric(snakemake@params[["crs"]])
basemap <- snakemake@params[["basemap"]]
use_elevation_bg <- isTRUE(snakemake@params[["use_elevation_bg"]])

plot_title <- snakemake@params[["plot_title"]]
if (is.null(plot_title) || plot_title == "") {
  plot_title <- paste("sPCA", score_label)
}

map_params <- list(
  boundary = snakemake@params[["boundary"]],
  crs = crs,
  basemap = basemap,
  land_colour = snakemake@params[["land_colour"]],
  sea_colour = snakemake@params[["sea_colour"]],
  expand = as.logical(snakemake@params[["expand"]]),
  arrow = as.logical(snakemake@params[["arrow"]]),
  arrow_size = as.numeric(snakemake@params[["arrow_size"]]),
  arrow_position = snakemake@params[["arrow_position"]],
  scalebar = as.logical(snakemake@params[["scalebar"]]),
  scalebar_size = as.numeric(snakemake@params[["scalebar_size"]]),
  scalebar_position = snakemake@params[["scalebar_position"]],
  plot_title = plot_title,
  axis_title_size = as.numeric(snakemake@params[["axis_title_size"]]),
  axis_text_size = as.numeric(snakemake@params[["axis_text_size"]]),
  basemap_border = as.logical(snakemake@params[["basemap_border"]]),
  basemap_border_col = snakemake@params[["basemap_border_col"]],
  basemap_border_lwd = as.numeric(snakemake@params[["basemap_border_lwd"]])
)
map_params$basemap <- resolve_map_basemap(use_elevation_bg, snakemake@input, basemap)
map_params$raster_is_elevation_dem <- isTRUE(use_elevation_bg)

coords_transformed <- transform_coordinates(coords_df, crs)
plot_df_transformed <- transform_coordinates(plot_df, crs)

p <- create_basemap(coords_df, map_params) +
  geom_point(
    data = plot_df_transformed,
    aes(x = Lon, y = Lat, color = spca_score),
    size = as.numeric(snakemake@params[["point_size"]]),
    alpha = as.numeric(snakemake@params[["point_alpha"]])
  ) +
  scale_color_gradient2(
    name = score_label,
    low = snakemake@params[["low_colour"]],
    mid = snakemake@params[["mid_colour"]],
    high = snakemake@params[["high_colour"]],
    midpoint = 0
  ) +
  theme(legend.position = "right")

dir.create(dirname(snakemake@output[["plot"]]), recursive = TRUE, showWarnings = FALSE)
ggsave(
  filename = snakemake@output[["plot"]],
  plot = p,
  width = width,
  height = height,
  dpi = dpi,
  device = "pdf"
)
saveRDS(p, snakemake@output[["plot_rds"]])

message("Wrote ", snakemake@output[["plot"]])
