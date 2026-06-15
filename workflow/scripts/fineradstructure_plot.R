#!/usr/bin/env Rscript
# Plots fineRADstructure output using the workflow from millanek/fineRADstructure
# fineRADstructurePlot.R (Daniel Lawson / Milan Malinsky; GPL-3).
# See FinestructureLibrary.R in this directory.

pdf(NULL)

log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")
on.exit(
  {
    sink(type = "message")
    sink(type = "output")
    close(log_file)
  },
  add = TRUE
)

snakemake@source("FinestructureLibrary.R")

chunkfile <- snakemake@input[["chunks"]]
mcmcfile <- snakemake@input[["mcmc"]]
treefile <- snakemake@input[["mcmcTree"]]
indpopdata_file <- snakemake@input[["indpopdata"]]

simple_pdf <- snakemake@output[["simple_pdf"]]
popavg_pdf <- snakemake@output[["popavg_pdf"]]
labeled_pdf <- snakemake@output[["labeled_pdf"]]
output_rds <- snakemake@output[["rds"]]

maxIndv <- as.numeric(snakemake@params[["max_indv"]])
maxPop <- as.numeric(snakemake@params[["max_pop"]])
population_columns <- snakemake@params[["population_columns"]]
population_colors  <- snakemake@params[["population_colors"]]

# Normalise: YAML list or single string → character vector; NULL stays NULL
if (!is.null(population_columns)) {
  population_columns <- as.character(unlist(population_columns))
  if (length(population_columns) == 0 ||
      all(population_columns %in% c("", "NULL", "null"))) {
    population_columns <- NULL
  }
}

cat("=== fineRADstructure plotting (upstream-style) ===\n")
cat("chunks:", chunkfile, "\n")
cat("mcmc:", mcmcfile, "\n")
cat("tree:", treefile, "\n")
cat("indpopdata:", indpopdata_file, "\n")
cat("population_columns:", if (is.null(population_columns)) "NULL (labeled plot will be skipped)"
                           else paste(population_columns, collapse = ", "), "\n")
cat("================================================\n\n")

some.colors_end <- MakeColorYRP(final = c(0.2, 0.2, 0.2))

dataraw <- as.matrix(read.table(
  chunkfile,
  row.names = 1,
  header = TRUE,
  skip = 1,
  check.names = FALSE,
  comment.char = "",
  quote = ""
))

mcmcxml <- xmlTreeParse(mcmcfile)
mcmcdata <- as.data.frame.myres(mcmcxml)

treexml <- xmlTreeParse(treefile)
ttree <- extractTree(treexml)

if (!is.null(ttree$node.label)) {
  idx <- nzchar(ttree$node.label)
  if (any(idx)) {
    suppressWarnings({
      num <- as.numeric(ttree$node.label[idx])
      ok <- !is.na(num)
      w <- which(idx)
      ttree$node.label[w[ok]] <- format(num[ok], digits = 2)
    })
  }
}

tdend <- myapetodend(ttree, factor = 1)

mapstate <- extractValue(treexml, "Pop")
if (!length(mapstate)) {
  stop("No Pop states in tree XML (extractValue).")
}
mapstate_use <- mapstate[length(mapstate)]
mapstatelist <- popAsList(mapstate_use)
# popAsList() uses sapply(); equal-sized populations become a matrix and break
# getPopMeanMatrix() / plotFinestructure() — force a list of character vectors.
if (is.matrix(mapstatelist)) {
  mapstatelist <- lapply(
    seq_len(ncol(mapstatelist)),
    function(j) as.character(mapstatelist[, j])
  )
} else {
  mapstatelist <- lapply(mapstatelist, function(z) as.character(unlist(z)))
}

popnames <- lapply(mapstatelist, NameSummary)
popnamesplot <- lapply(mapstatelist, NameMoreSummary)
names(popnames) <- popnamesplot
names(popnamesplot) <- popnamesplot

popdend <- makemydend(tdend, mapstatelist)
popdend <- fixMidpointMembers(popdend)
popdendclear <- makemydend(tdend, mapstatelist, "NameMoreSummary")
popdendclear <- fixMidpointMembers(popdendclear)

fullorder <- labels(tdend)
miss <- setdiff(fullorder, rownames(dataraw))
if (length(miss)) {
  stop(
    "Tree tips missing from chunk matrix: ",
    paste(head(miss, 30), collapse = ", ")
  )
}
datamatrix <- dataraw[fullorder, fullorder, drop = FALSE]

tmpmat <- datamatrix
tmpmat[tmpmat > maxIndv] <- maxIndv
pdf(file = simple_pdf, height = 25, width = 25)
plotFinestructure(
  tmpmat,
  dimnames(tmpmat)[[1]],
  dend = tdend,
  cols = some.colors_end,
  cex.axis = 1.1,
  edgePar = list(p.lwd = 0, t.srt = 90, t.off = -0.1, t.cex = 1.2)
)
dev.off()
cat("Wrote simple coancestry plot\n")

popmeanmatrix <- getPopMeanMatrix(datamatrix, mapstatelist)
tmpmat <- popmeanmatrix
tmpmat[tmpmat > maxPop] <- maxPop
pdf(file = popavg_pdf, height = 20, width = 20)
plotFinestructure(
  tmpmat,
  dimnames(tmpmat)[[1]],
  dend = tdend,
  cols = some.colors_end,
  cex.axis = 1.1,
  edgePar = list(p.lwd = 0, t.srt = 90, t.off = -0.1, t.cex = 1.2)
)
dev.off()
cat("Wrote population-averaged coancestry plot\n")

# ── Plot 3: labeled with optional population annotation bars ─────────────────
if (is.null(population_columns)) {
  cat("population_columns is NULL — writing placeholder for labeled PDF\n")
  pdf(file = labeled_pdf, height = 4, width = 6)
  par(mar = c(0, 0, 0, 0))
  plot.new()
  text(0.5, 0.5, "population_columns not configured\n— labeled plot skipped",
       adj = 0.5, cex = 1.2)
  dev.off()
} else {

  # --- read indpopdata and validate columns -----------------------------------
  indpopdata_df <- read.table(
    indpopdata_file, header = TRUE, sep = "\t",
    stringsAsFactors = FALSE, check.names = FALSE, quote = ""
  )
  if (!"Ind" %in% colnames(indpopdata_df))
    stop("indpopdata must have an 'Ind' column.")
  missing_cols <- setdiff(population_columns, colnames(indpopdata_df))
  if (length(missing_cols) > 0)
    stop(sprintf(
      "population_columns not found in indpopdata: %s. Available: %s",
      paste(missing_cols, collapse = ", "),
      paste(colnames(indpopdata_df), collapse = ", ")
    ))

  # per-column individual → value lookup
  ind_to_pop <- lapply(setNames(population_columns, population_columns), function(col)
    setNames(indpopdata_df[[col]], indpopdata_df[["Ind"]]))

  # --- build color maps -------------------------------------------------------
  # 20-color qualitative palette (no external package required)
  qual_pal <- c(
    "#E41A1C","#377EB8","#4DAF4A","#FF7F00","#984EA3",
    "#A65628","#F781BF","#999999","#66C2A5","#FC8D62",
    "#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494",
    "#B3B3B3","#1B9E77","#D95F02","#7570B3","#E7298A"
  )
  user_colors <- if (!is.null(population_colors))
    unlist(population_colors, use.names = TRUE) else NULL

  col_color_maps <- lapply(population_columns, function(col) {
    vals <- sort(unique(na.omit(indpopdata_df[[col]])))
    n <- length(vals)
    base_pal <- if (n <= length(qual_pal)) qual_pal[seq_len(n)]
                else rainbow(n, s = 0.8, v = 0.85)
    auto_map <- setNames(base_pal, vals)
    if (!is.null(user_colors))
      auto_map[names(user_colors)[names(user_colors) %in% vals]] <-
        user_colors[names(user_colors) %in% vals]
    auto_map
  })
  names(col_color_maps) <- population_columns

  # --- relabel popdendclear with the first population column -----------------
  first_col <- population_columns[[1]]
  first_lookup <- ind_to_pop[[first_col]]

  make_pop_label <- function(members) {
    pops <- first_lookup[members]
    pops <- pops[!is.na(pops)]
    if (length(pops) == 0) return(paste(members, collapse = ";"))
    tab <- sort(table(pops), decreasing = TRUE)
    paste(paste0(tab, names(tab)), collapse = ";")
  }
  nms_to_custom <- setNames(sapply(mapstatelist, make_pop_label),
                            sapply(mapstatelist, NameMoreSummary))
  popdendclear_labeled <- dendrapply(popdendclear, function(node) {
    lbl <- attr(node, "label")
    if (!is.null(lbl) && lbl %in% names(nms_to_custom))
      attr(node, "label") <- nms_to_custom[[lbl]]
    node
  })
  cat("Relabelled clusters using indpopdata column:", first_col, "\n")

  # --- label positions (individual scale → K centres) ------------------------
  mappopcorrectorder <- NameExpand(labels(popdend))
  mappopsizes <- sapply(mappopcorrectorder, length)
  labellocs <- PopCenters(mappopsizes)

  # --- annotation bar layout --------------------------------------------------
  # All coordinates are in heatmap user units (1 unit = 1 individual spacing).
  # Bars are drawn below the x-axis (y < 0) and left of the y-axis (x < 0)
  # by calling rect() with xpd=TRUE after plotFinestructure returns.
  n_bars   <- length(population_columns)
  BAR_H    <- 1.0   # bar height (user units)
  BAR_GAP  <- 0.35  # gap between consecutive bars
  FIRST_GAP <- 0.5  # gap from axis edge to first bar

  # y-top of bar i (1-indexed), measured downward from y=0
  bar_y_top <- function(i) -(FIRST_GAP + (i - 1) * (BAR_H + BAR_GAP))
  bar_y_bot <- function(i) bar_y_top(i) - BAR_H

  text_labelsoff <- -bar_y_bot(n_bars) + 0.8   # push labels below last bar

  # --- draw the labeled plot --------------------------------------------------
  pdf(file = labeled_pdf, height = 25, width = 25)

  plotFinestructure(
    tmpmat,
    dimnames(tmpmat)[[1]],
    labelsx  = labels(popdendclear_labeled),
    labelsatx = labellocs,
    labelsaty = labellocs,
    labelsoff = c(text_labelsoff, text_labelsoff),
    xcrt = 0,
    ycrt = 45,
    cols = some.colors_end,
    dend = tdend,
    cex.axis = 1.1,
    edgePar = list(p.lwd = 0, t.srt = 90, t.off = -0.1, t.cex = 1.2),
    hmmar = c(3, 0, 0, 1)
  )

  # Draw per-individual colored annotation bars along both axes
  for (bi in seq_along(population_columns)) {
    col_name <- population_columns[[bi]]
    col_map  <- col_color_maps[[col_name]]
    lookup   <- ind_to_pop[[col_name]]
    yt <- bar_y_top(bi)
    yb <- bar_y_bot(bi)

    for (i in seq_along(fullorder)) {
      v    <- lookup[fullorder[i]]
      fill <- if (!is.na(v) && v %in% names(col_map)) col_map[[v]] else "#CCCCCC"
      # bottom annotation (x-axis)
      rect(i - 0.5, yb, i + 0.5, yt, col = fill, border = NA, xpd = TRUE)
      # left annotation (y-axis), coordinates mirrored
      rect(yb, i - 0.5, yt, i + 0.5, col = fill, border = NA, xpd = TRUE)
    }

    # column name rotated vertically, centred alongside the y-axis strip
    text(yb - 0.15, (length(fullorder) + 1) / 2,
         col_name, srt = 90, adj = 0.5, cex = 0.75, xpd = TRUE)
  }

  # Legend in the top-right empty panel (panel 1 in the layout)
  par(mfg = c(1, 2))
  par(mar = c(0, 0, 0.5, 0.5))
  plot.new()
  y_leg <- 0.97
  leg_cex <- 0.78
  sq <- 0.035  # legend square half-height
  for (bi in seq_along(population_columns)) {
    col_name <- population_columns[[bi]]
    col_map  <- col_color_maps[[col_name]]
    text(0.5, y_leg, col_name, font = 2, adj = 0.5, cex = leg_cex)
    y_leg <- y_leg - 0.045
    for (val in names(col_map)) {
      if (y_leg < 0.02) break
      rect(0.04, y_leg - sq, 0.12, y_leg + sq,
           col = col_map[[val]], border = NA, xpd = FALSE)
      text(0.16, y_leg, val, adj = 0, cex = leg_cex * 0.88)
      y_leg <- y_leg - 0.042
    }
    y_leg <- y_leg - 0.02
  }

  dev.off()
  cat("Wrote labeled population-averaged coancestry plot\n")
}

saveRDS(
  list(
    tree = ttree,
    tdend = tdend,
    popdend = popdend,
    popdendclear = popdendclear,
    mcmcdata = mcmcdata,
    datamatrix = datamatrix,
    popmeanmatrix = popmeanmatrix,
    mapstatelist = mapstatelist
  ),
  output_rds
)

cat("Done.\n")
