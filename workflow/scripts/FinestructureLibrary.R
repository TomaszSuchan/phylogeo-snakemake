#!/usr/bin/env Rscript
# Adapted from the upstream fineRADstructure/FinestructureLibrary.R helpers.
# This trimmed version keeps the functions needed for Snakemake plotting while
# avoiding the original NAME<number>-specific population naming helpers.

simplelabels <- function(x, labmode = "mcmc", oneminus = TRUE) {
  if (x == "") {
    return("")
  }
  if (labmode == "mcmc") {
    x <- strsplit(x, "-")[[1]][1]
  } else if (labmode == "pop") {
    x <- strsplit(x, "-")[[1]]
    if (length(x) > 1) {
      x <- x[2]
    } else {
      return("")
    }
  } else {
    return(x)
  }
  if (oneminus) {
    return(1 - as.numeric(x))
  }
  as.numeric(x)
}

extractTree <- function(txml, labmode = "mcmc", oneminus = TRUE, hidecertain = FALSE) {
  output_children <- txml$doc$children$outputFile
  tree_nodes <- output_children[which(names(output_children) == "Tree")]
  if (length(tree_nodes) == 0) {
    stop("No <Tree> node found in fineSTRUCTURE XML.")
  }
  tree <- ape::read.tree(text = XML::xmlValue(tree_nodes[[length(tree_nodes)]]))
  if (!is.null(tree$node.label)) {
    tree$node.label <- as.vector(
      sapply(tree$node.label, simplelabels, labmode = labmode, oneminus = oneminus)
    )
    if (hidecertain) {
      tree$node.label[which(tree$node.label == "1")] <- ""
    }
  }
  tree
}

getV <- function(it, v) {
  XML::xmlValue(XML::xmlChildren(it)[[v]])
}

extractValue <- function(txml, v, getNames = FALSE) {
  tmpits <- txml$doc$children$outputFile[which(names(txml$doc$children$outputFile) == "Iteration")]
  if (getNames) {
    num <- sapply(tmpits, getV, v = "Number")
  }
  res <- sapply(tmpits, getV, v = v)
  if (!getNames) {
    names(res) <- NULL
  } else {
    names(res) <- num
  }
  res
}

popAsList <- function(s) {
  s <- gsub(")", "", s, fixed = TRUE)
  parts <- strsplit(s, c("(", ")"), fixed = TRUE)[[1]]
  parts <- parts[nzchar(parts)]
  lapply(parts, function(x) strsplit(x, ",", fixed = TRUE)[[1]])
}

matColSums <- function(mat, poplist) {
  res <- matrix(0, nrow = nrow(mat), ncol = length(poplist))
  colindex <- lapply(poplist, function(x) which(colnames(mat) %in% x))
  res <- t(apply(mat, 1, function(x) {
    sapply(colindex, function(y) sum(x[y]))
  }))
  if (nrow(res) == 1) {
    res <- t(res)
  }
  colnames(res) <- names(poplist)
  rownames(res) <- rownames(mat)
  res
}

getPopCountMatrix <- function(datamatrix, mapstatelist, remdiag = TRUE) {
  popnmatrix <- datamatrix
  popnmatrix[] <- 1
  if (remdiag) {
    diag(popnmatrix) <- 0
  }
  popcountmatrix <- matColSums(popnmatrix, mapstatelist)
  matColSums(t(popcountmatrix), mapstatelist)
}

getPopMatrix <- function(datamatrix, mapstatelist, correction = 1, remdiag = TRUE) {
  datamatrix[is.nan(datamatrix[])] <- 0
  popcountmatrix <- getPopCountMatrix(datamatrix, mapstatelist, remdiag = remdiag)
  popmatrix <- matColSums(datamatrix, mapstatelist)
  popmatrix <- t(matColSums(t(popmatrix), mapstatelist))
  popmatrix <- popmatrix / popcountmatrix
  popmatrix[is.nan(popmatrix[])] <- 0
  popmatrix / correction
}

getPopIndices <- function(indnames, mapstatelist) {
  sapply(indnames, function(x) {
    which(sapply(mapstatelist, function(y) x %in% y))
  })
}

getPopMeanMatrix <- function(datamatrix, mapstatelist) {
  popmeanmatrix <- datamatrix
  popindices <- getPopIndices(rownames(datamatrix), mapstatelist)
  popmatrix <- getPopMatrix(datamatrix, mapstatelist)
  for (i in seq_len(nrow(popmeanmatrix))) {
    for (j in seq_len(ncol(popmeanmatrix))) {
      popmeanmatrix[i, j] <- popmatrix[popindices[i], popindices[j]]
    }
  }
  popmeanmatrix
}

fixMidpointMembers <- function(x) {
  attr(x, "x.member") <- NULL
  attr(x, "value") <- NULL
  if (is.leaf(x)) {
    attr(x, "label") <- attr(x, "label")[[1]]
    attr(x, "midpoint") <- NULL
    attr(x, "members") <- 1
  }
  if (!is.leaf(x)) {
    attr(x, "members") <- length(labels(x))
  }
  x
}

setNodeLabels <- function(din, mylabs) {
  i <- 0
  edgeLab <- function(n, lablist) {
    if (!is.leaf(n)) {
      i <<- i + 1
      if (nchar(lablist[i]) > 0 && i > 1) {
        attr(n, "edgetext") <- lablist[i]
      }
    }
    n
  }
  dendrapply(din, edgeLab, mylabs)
}

my.as.hclust.phylo <- function(x, tol = 0.01, ...) {
  if (!ape::is.ultrametric(x, tol)) {
    stop("the tree is not ultrametric")
  }
  if (!ape::is.binary.tree(x)) {
    stop("the tree is not binary")
  }
  n <- length(x$tip.label)
  bt <- rev(ape::branching.times(x))
  N <- length(bt)
  nm <- as.numeric(names(bt))
  merge <- matrix(NA, N, 2)
  for (i in seq_len(N)) {
    ind <- which(x$edge[, 1] == nm[i])
    for (k in 1:2) {
      merge[i, k] <- if (x$edge[ind[k], 2] <= n) {
        -x$edge[ind[k], 2]
      } else {
        which(nm == x$edge[ind[k], 2])
      }
    }
  }
  names(bt) <- NULL
  obj <- list(
    merge = merge,
    height = bt,
    order = 1:(N + 1),
    labels = x$tip.label,
    call = match.call(),
    method = "unknown"
  )
  class(obj) <- "hclust"
  obj
}

positiveheights <- function(x) {
  if (attr(x, "height") < 0) {
    attr(x, "height") <- 0.0001
  }
  x
}

flattenheights <- function(x, factor = 0.25) {
  if (attr(x, "height") > 0) {
    attr(x, "height") <- attr(x, "height")^factor
  }
  x
}

myapetodend <- function(ttree, simplify = FALSE, tol = 0.1, factor = 0.25) {
  nodelab <- ttree$node.label
  ttree$node.label <- NULL
  htree <- my.as.hclust.phylo(ttree, tol)
  dend <- as.dendrogram(htree)
  dend <- dendrapply(dend, positiveheights)
  dend <- dendrapply(dend, flattenheights, factor = factor)
  if (!is.null(nodelab)) {
    dend <- setNodeLabels(dend, nodelab)
  }
  dend
}

makemydend <- function(tdend, lablist, summary = "NameLessSummary") {
  for (i in 1:2) {
    j <- i + 1
    if (j == 3) {
      j <- 1
    }
    test <- which(sapply(lablist, function(x) all(labels(tdend[[i]]) %in% x)))
    if (length(test) > 0) {
      if (all(lablist[[test]] %in% labels(tdend[[i]]))) {
        testj <- cut(tdend, attr(tdend, "height"))$upper
        attr(testj[[i]], "height") <- 0
        attr(testj[[i]], "label") <- get(summary)(lablist[[test]])
        testj[[j]] <- tdend[[j]]
        tdend <- testj
      }
    } else if (length(labels(tdend[[i]]) > 1)) {
      tdend[[i]] <- makemydend(tdend[[i]], lablist, summary)
    }
  }
  tdend <- dendrapply(tdend, fixMidpointMembers)
  suppressWarnings(midcache.dendrogram(tdend))
}

PopCenters <- function(popsizes) {
  csum <- ret <- rep(0, length(popsizes))
  csum[1] <- popsizes[1]
  ret[1] <- (popsizes[1] + 1) / 2
  for (i in 2:length(popsizes)) {
    csum[i] <- csum[i - 1] + popsizes[i]
    ret[i] <- csum[i - 1] + (popsizes[i] + 1) / 2
  }
  ret
}

MakeColorYRP <- function(colby = 0.05, final = NULL) {
  tmp <- c(
    rgb(1, seq(1, 0, -colby), 0),
    rgb(1, 0, seq(colby, 1, colby)),
    rgb(seq(1 - colby, 0, -colby), 0, 1.0)
  )
  if (is.null(final)) {
    return(tmp)
  }
  c(tmp, rgb(final[1], final[2], final[3]))
}

fs.plot.dendrogram <- function(
  x, type = c("rectangle", "triangle"), center = FALSE,
  edge.root = is.leaf(x) || !is.null(attr(x, "edgetext")),
  nodePar = NULL, edgePar = list(), leaflab = c("perpendicular", "textlike", "none"),
  dLeaf = NULL, xlab = "", ylab = "", xaxt = "n", yaxt = "s",
  horiz = FALSE, frame.plot = FALSE, xlim, ylim, height = 0, ...
) {
  if (height > 0) {
    x <- dendrapply(x, function(x, h) {
      attr(x, "height") <- attr(x, "height") + h
      x
    }, height)
  }
  type <- match.arg(type)
  leaflab <- match.arg(leaflab)
  hgt <- attr(x, "height")
  if (edge.root && is.logical(edge.root)) {
    edge.root <- 0.0625 * if (is.leaf(x)) 1 else hgt
  }
  mem.x <- stats:::.memberDend(x)
  yTop <- hgt + edge.root
  if (center) {
    x1 <- 0.5
    x2 <- mem.x + 0.5
  } else {
    x1 <- 1
    x2 <- mem.x
  }
  xl. <- c(x1 - 0.5, x2 + 0.5)
  yl. <- c(0, yTop)
  if (horiz) {
    tmp <- xl.
    xl. <- rev(yl.)
    yl. <- tmp
    tmp <- xaxt
    xaxt <- yaxt
    yaxt <- tmp
  }
  if (missing(xlim) || is.null(xlim)) {
    xlim <- xl.
  }
  if (missing(ylim) || is.null(ylim)) {
    ylim <- yl.
  }
  plot(0, xlim = xlim, ylim = ylim, type = "n", xlab = xlab, ylab = ylab,
       xaxt = xaxt, yaxt = yaxt, frame.plot = frame.plot, ...)
  if (is.null(dLeaf)) {
    dLeaf <- 0.75 * if (horiz) strwidth("w") else strheight("x")
  }
  if (edge.root) {
    x0 <- fs.plotNodeLimit(x1, x2, x, center)$x
    if (horiz) {
      segments(hgt, x0, yTop, x0)
    } else {
      segments(x0, hgt, x0, yTop)
    }
  }
  fs.plotNode(x1, x2, x, type = type, center = center, leaflab = leaflab,
              dLeaf = dLeaf, nodePar = nodePar, edgePar = edgePar, horiz = horiz)
}

fs.plotNode <- function(x1, x2, subtree, type, center, leaflab, dLeaf,
                        nodePar, edgePar, horiz = FALSE) {
  inner <- !is.leaf(subtree) && x1 != x2
  yTop <- attr(subtree, "height")
  bx <- fs.plotNodeLimit(x1, x2, subtree, center)
  xTop <- bx$x

  Xtract <- function(nam, L, default, indx) {
    rep(if (nam %in% names(L)) L[[nam]] else default, length.out = indx)[indx]
  }
  asTxt <- function(x) if (is.character(x) || is.expression(x) || is.null(x)) x else as.character(x)

  hasP <- !is.null(nPar <- attr(subtree, "nodePar"))
  if (!hasP) {
    nPar <- nodePar
  }
  i <- if (inner || hasP) 1 else 2
  lab.cex <- Xtract("lab.cex", nPar %||% list(), default = c(1, 1), i)
  lab.col <- Xtract("lab.col", nPar %||% list(), default = par("col"), i)
  lab.font <- Xtract("lab.font", nPar %||% list(), default = par("font"), i)

  if (is.leaf(subtree)) {
    if (leaflab == "perpendicular") {
      if (horiz) {
        text(yTop + dLeaf * lab.cex, xTop, asTxt(attr(subtree, "label")),
             xpd = TRUE, adj = c(0, 0.5), cex = lab.cex, col = lab.col, font = lab.font)
      } else {
        text(xTop, yTop - dLeaf * lab.cex, asTxt(attr(subtree, "label")),
             xpd = TRUE, srt = 90, adj = 1, cex = lab.cex, col = lab.col, font = lab.font)
      }
    }
  } else if (inner) {
    for (k in seq_along(subtree)) {
      child <- subtree[[k]]
      yBot <- attr(child, "height")
      if (is.null(yBot)) {
        yBot <- 0
      }
      xBot <- if (center) mean(bx$limit[k:(k + 1)]) else bx$limit[k] + stats:::.midDend(child)
      hasE <- !is.null(ePar <- attr(child, "edgePar"))
      if (!hasE) {
        ePar <- edgePar
      }
      i <- if (!is.leaf(child) || hasE) 1 else 2
      col <- Xtract("col", ePar, default = par("col"), i)
      lty <- Xtract("lty", ePar, default = par("lty"), i)
      lwd <- Xtract("lwd", ePar, default = par("lwd"), i)
      if (horiz) {
        if (type == "triangle") {
          segments(yTop, xTop, yBot, xBot, col = col, lty = lty, lwd = lwd)
        } else {
          segments(yTop, xTop, yTop, xBot, col = col, lty = lty, lwd = lwd)
          segments(yTop, xBot, yBot, xBot, col = col, lty = lty, lwd = lwd)
        }
      } else {
        if (type == "triangle") {
          segments(xTop, yTop, xBot, yBot, col = col, lty = lty, lwd = lwd)
        } else {
          segments(xTop, yTop, xBot, yTop, col = col, lty = lty, lwd = lwd)
          segments(xBot, yTop, xBot, yBot, col = col, lty = lty, lwd = lwd)
        }
      }
      fs.plotNode(bx$limit[k], bx$limit[k + 1], child, type, center, leaflab, dLeaf,
                  nodePar, edgePar, horiz)
    }
  }
  invisible()
}

fs.plotNodeLimit <- function(x1, x2, subtree, center) {
  inner <- !is.leaf(subtree) && x1 != x2
  if (inner) {
    K <- length(subtree)
    mTop <- stats:::.memberDend(subtree)
    limit <- integer(K)
    xx1 <- x1
    for (k in 1L:K) {
      m <- stats:::.memberDend(subtree[[k]])
      xx1 <- xx1 + (if (center) (x2 - x1) * m / mTop else m)
      limit[k] <- xx1
    }
    limit <- c(x1, limit)
  } else {
    limit <- c(x1, x2)
  }
  mid <- attr(subtree, "midpoint")
  center <- center || (inner && !is.numeric(mid))
  x <- if (center) mean(c(x1, x2)) else x1 + if (inner) mid else 0
  list(x = x, limit = limit)
}

plotFinestructure <- function(
  tmpmat, labelsx, labelsy = NULL, labelsatx = NULL, labelsaty = NULL,
  cols = NULL, dend = NULL, labmargin = 8, layoutd = 0.2, layoutf = 0.1,
  cex.axis = 0.5, xcrt = 0, ycrt = 0, colscale = NULL, text.col = NULL,
  ignorebelow = 0, ignoreabove = Inf, nodePar = list(cex = 0, lab.cex = 0.5, las = 2),
  edgePar = list(p.lwd = 0, t.srt = 0, t.off = 0.2), dendmar = c(0, 0, 2, 1),
  scalemar = c(2, 5, 2, 1), hmmar = c(0, 0, 0, 1), cex.scale = 1,
  scalelocs = NULL, scalenum = 10, scalesignif = 3, scalelabel = "",
  optcex = 1.0, optpch = 20, optcol = "black", labelsoff = c(1, 1),
  tickmarks = 1, startpt = "bottomleft", dolayout = TRUE, main = "",
  cex.main = 1, adj = 0, optpts = NULL
) {
  if (is.null(cols)) {
    cols <- MakeColorYRP()
  }
  if (is.null(labelsy)) {
    labelsy <- labelsx
  }
  if (is.null(labelsatx)) {
    labelsatx <- seq_len(ncol(tmpmat))
  }
  if (is.null(labelsaty)) {
    labelsaty <- seq_len(nrow(tmpmat))
  }

  valid_vals <- as.numeric(tmpmat)
  valid_vals <- valid_vals[is.finite(valid_vals)]
  valid_vals <- valid_vals[valid_vals > ignorebelow & valid_vals < ignoreabove]
  if (!length(valid_vals)) {
    valid_vals <- c(0, 1)
  }
  if (is.null(colscale)) {
    colscale <- range(valid_vals)
  }
  breaks <- seq(colscale[1], colscale[2], length.out = length(cols) + 1)

  if (dolayout) {
    layout(matrix(c(2, 1, 4, 3), 2, 2, byrow = TRUE),
           widths = c(1 - layoutf, layoutf), heights = c(layoutd, 1 - layoutd))
  }

  par(mar = c(0, labmargin, 0, 0) + dendmar)
  if (is.null(dend)) {
    plot.new()
    title(main = main, cex.main = cex.main, adj = adj)
  } else {
    fs.plot.dendrogram(
      dend, horiz = FALSE, axes = FALSE, xaxs = "i", leaflab = "none",
      nodePar = nodePar, edgePar = edgePar, main = main, cex.main = cex.main, adj = adj
    )
  }

  par(mar = c(0, 0, 0, 0))
  plot.new()

  par(mar = c(labmargin, labmargin, 0, 1) + hmmar)
  n <- nrow(tmpmat)
  image(
    x = seq_len(ncol(tmpmat)),
    y = seq_len(nrow(tmpmat)),
    z = t(tmpmat[n:1, , drop = FALSE]),
    col = cols,
    breaks = breaks,
    axes = FALSE,
    xlab = "",
    ylab = "",
    useRaster = TRUE
  )
  box()

  if (length(labelsatx) > 0) {
    axis(1, at = labelsatx, labels = FALSE, tck = -0.01 * tickmarks)
    text(labelsatx, par("usr")[3] - labelsoff[1], labels = labelsx, srt = xcrt,
         xpd = TRUE, adj = if (xcrt == 0) c(0.5, 1) else c(1, 0.5),
         cex = cex.axis, col = text.col %||% par("col"))
  }
  if (length(labelsaty) > 0) {
    y_at <- n + 1 - labelsaty
    axis(2, at = y_at, labels = FALSE, tck = -0.01 * tickmarks)
    text(par("usr")[1] - labelsoff[2], y_at, labels = labelsy, srt = ycrt,
         xpd = TRUE, adj = if (ycrt == 0) c(1, 0.5) else c(0.5, 1),
         cex = cex.axis, col = text.col %||% par("col"))
  }

  if (!is.null(optpts)) {
    tmp <- sapply(seq_len(nrow(optpts)), function(x) {
      plist <- which(optpts[, x] == 1)
      points(rep(x, length(plist)), plist, pch = optpch, cex = optcex, col = optcol)
      invisible(NULL)
    })
  }

  par(mar = c(labmargin, 1, 0, 4) + scalemar)
  scale_vals <- matrix(seq(colscale[1], colscale[2], length.out = 200), ncol = 1)
  image(
    x = 1,
    y = seq(colscale[1], colscale[2], length.out = 200),
    z = scale_vals,
    col = cols,
    breaks = breaks,
    axes = FALSE,
    xlab = "",
    ylab = "",
    useRaster = TRUE
  )
  box()
  if (is.null(scalelocs)) {
    scalelocs <- pretty(colscale, n = scalenum)
  }
  axis(4, at = scalelocs, labels = signif(scalelocs, scalesignif), las = 2, cex.axis = cex.scale)
  if (nzchar(scalelabel)) {
    mtext(scalelabel, side = 4, line = 2.5, cex = cex.scale)
  }
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
