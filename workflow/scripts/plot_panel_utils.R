# Helpers for saving and reusing pre-rendered ggplot panel grobs.
# Map plots with ggspatial/terra layers must be rendered once at plot time;
# the resulting gtable can then be arranged in facet panels without rebuilding.

ggplot_to_panel_grob <- function(plot) {
  ggplot2::ggplotGrob(ggplot2::ggplot_build(plot))
}

save_plot_panel_rds <- function(plot, path) {
  panel <- ggplot_to_panel_grob(plot)
  saveRDS(panel, file = path)
}

read_panel_grob <- function(path, panel_kind = c("auto", "map")) {
  panel_kind <- match.arg(panel_kind)
  if (!file.exists(path)) {
    stop("Panel RDS not found: ", path)
  }
  obj <- readRDS(path)

  if (inherits(obj, "gtable")) {
    return(obj)
  }
  if (inherits(obj, "ggplot")) {
    if (identical(panel_kind, "map")) {
      stop(
        "Map panel RDS at ", path, " stores an unrendered ggplot object. ",
        "Re-run the corresponding mapmixture_* rule to regenerate panel RDS."
      )
    }
    return(ggplot_to_panel_grob(obj))
  }
  if (inherits(obj, "grob")) {
    return(obj)
  }

  stop(
    "Unsupported panel RDS class in ", path, ": ",
    paste(class(obj), collapse = "/")
  )
}
