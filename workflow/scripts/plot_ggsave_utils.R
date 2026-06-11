# Shared ggsave helper for PDF output.
#
# The default R pdf() device uses Type 1 fonts (Helvetica, etc.) that lack many
# Unicode glyphs (e.g. Slovak ľ/č, Greek Δ). Missing glyphs render as dots.
# cairo_pdf embeds system fonts with broader Unicode coverage when Cairo is available.

cairo_pdf_device <- function() {
  if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else NULL
}

ggsave_pdf <- function(
  filename,
  plot = ggplot2::last_plot(),
  ...,
  device = cairo_pdf_device()
) {
  ggplot2::ggsave(
    filename = filename,
    plot = plot,
    device = device,
    ...
  )
}
