# Shared parser for currentNe v2 text output.

# Point estimate and 50%/90% CI bounds from the whole-genome block, falling back
# to the between-chromosome block when the former is absent.
parse_currentne2 <- function(path, population) {
  if (!file.exists(path)) {
    stop("Missing currentNe2 output: ", path)
  }
  lines <- readLines(path, warn = FALSE)
  get_after <- function(pattern, lines_use = lines) {
    idx <- grep(pattern, lines_use, perl = TRUE)
    if (length(idx) == 0) return(NA_real_)
    for (i in idx) {
      if (i >= length(lines_use)) next
      val <- suppressWarnings(as.numeric(trimws(lines_use[[i + 1]])))
      if (is.finite(val)) return(val)
    }
    NA_real_
  }
  wg_start <- grep("integration over the whole genome", lines, fixed = TRUE)
  bc_start <- grep("LD between chromosomes", lines, fixed = TRUE)
  block <- lines
  if (length(wg_start) > 0) {
    end <- if (length(bc_start) > 0 && bc_start[1] > wg_start[1]) bc_start[1] - 1 else length(lines)
    block <- lines[wg_start[1]:end]
  } else if (length(bc_start) > 0) {
    block <- lines[bc_start[1]:length(lines)]
  }
  out <- data.frame(
    population = population,
    ne = get_after("^# Ne point estimate:", block),
    ci50_low = get_after("^# Lower limit of 50% CI:", block),
    ci50_high = get_after("^# Upper (bound|limit) of 50% CI:", block),
    ci90_low = get_after("^# Lower limit of 90% CI:", block),
    ci90_high = get_after("^# Upper limit of 90% CI:", block),
    stringsAsFactors = FALSE
  )
  if (!is.finite(out$ne[1])) {
    stop("Could not parse Ne point estimate from ", path)
  }
  out
}
