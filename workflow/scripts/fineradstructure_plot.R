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
population_column <- snakemake@params[["population_column"]]

cat("=== fineRADstructure plotting (upstream-style) ===\n")
cat("chunks:", chunkfile, "\n")
cat("mcmc:", mcmcfile, "\n")
cat("tree:", treefile, "\n")
cat("indpopdata:", indpopdata_file, "\n")
cat("population_column:", if (is.null(population_column)) "NULL (using individual-name summaries)" else population_column, "\n")
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

# Optionally relabel popdendclear leaves using an indpopdata column instead of
# the individual-name NameMoreSummary format (e.g. "3popA;2popB").
if (!is.null(population_column) && nchar(population_column) > 0) {
  indpopdata_df <- read.table(
    indpopdata_file, header = TRUE, sep = "\t",
    stringsAsFactors = FALSE, check.names = FALSE, quote = ""
  )
  if (!"Ind" %in% colnames(indpopdata_df)) {
    stop("indpopdata file must contain an 'Ind' column.")
  }
  if (!population_column %in% colnames(indpopdata_df)) {
    stop(sprintf(
      "population_column '%s' not found in indpopdata. Available columns: %s",
      population_column, paste(colnames(indpopdata_df), collapse = ", ")
    ))
  }
  ind_to_pop <- setNames(indpopdata_df[[population_column]], indpopdata_df[["Ind"]])

  make_pop_label <- function(members) {
    pops <- ind_to_pop[members]
    pops <- pops[!is.na(pops)]
    if (length(pops) == 0) return(paste(members, collapse = ";"))
    tab <- sort(table(pops), decreasing = TRUE)
    paste(paste0(tab, names(tab)), collapse = ";")
  }

  # Map NameMoreSummary(cluster) -> population-column label for each cluster
  label_map <- setNames(
    sapply(mapstatelist, make_pop_label),
    sapply(mapstatelist, NameMoreSummary)
  )

  popdendclear <- dendrapply(popdendclear, function(node) {
    lbl <- attr(node, "label")
    if (!is.null(lbl) && lbl %in% names(label_map)) {
      attr(node, "label") <- label_map[[lbl]]
    }
    node
  })
  cat("Relabelled clusters using indpopdata column:", population_column, "\n")
}

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

mappopcorrectorder <- NameExpand(labels(popdend))
mappopsizes <- sapply(mappopcorrectorder, length)
labellocs <- PopCenters(mappopsizes)
xcrt <- 0
ycrt <- 45

pdf(file = labeled_pdf, height = 25, width = 25)
plotFinestructure(
  tmpmat,
  dimnames(tmpmat)[[1]],
  labelsx = labels(popdendclear),
  labelsatx = labellocs,
  labelsaty = labellocs,  # mirror x: K labels at K population centres, not recycled over N ticks
  xcrt = xcrt,
  cols = some.colors_end,
  ycrt = ycrt,
  dend = tdend,
  cex.axis = 1.1,
  edgePar = list(p.lwd = 0, t.srt = 90, t.off = -0.1, t.cex = 1.2),
  hmmar = c(3, 0, 0, 1)
)
dev.off()
cat("Wrote labeled population-averaged coancestry plot\n")

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
