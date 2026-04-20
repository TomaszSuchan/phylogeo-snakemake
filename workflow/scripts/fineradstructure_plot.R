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

simple_pdf <- snakemake@output[["simple_pdf"]]
popavg_pdf <- snakemake@output[["popavg_pdf"]]
labeled_pdf <- snakemake@output[["labeled_pdf"]]
output_rds <- snakemake@output[["rds"]]

maxIndv <- as.numeric(snakemake@params[["max_indv"]])
maxPop <- as.numeric(snakemake@params[["max_pop"]])

cat("=== fineRADstructure plotting (upstream-style) ===\n")
cat("chunks:", chunkfile, "\n")
cat("mcmc:", mcmcfile, "\n")
cat("tree:", treefile, "\n")
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
