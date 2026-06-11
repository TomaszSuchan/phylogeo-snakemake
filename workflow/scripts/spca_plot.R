#!/usr/bin/env Rscript
# sPCA plotting following adegenet tutorial-spca.pdf

Sys.setenv(R_INSTALL_PKG = "")
options(device = function(file = NULL, ...) {
  if (is.null(file)) {
    pdf(NULL)
  } else {
    pdf(file, ...)
  }
})

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
  library(adegenet)
  library(ggplot2)
})

ggsave_utils <- tryCatch(
  file.path(dirname(normalizePath(snakemake@script)), "plot_ggsave_utils.R"),
  error = function(e) "workflow/scripts/plot_ggsave_utils.R"
)
if (file.exists(ggsave_utils)) {
  source(ggsave_utils)
} else {
  source("workflow/scripts/plot_ggsave_utils.R")
}

save_base_plot <- function(path, expr, width = 8, height = 6) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  pdf(path, width = width, height = height)
  on.exit(dev.off(), add = TRUE)
  eval(expr, envir = parent.frame())
  invisible(NULL)
}

unwrap_results <- function(rds_path) {
  obj <- readRDS(rds_path)
  if (inherits(obj, "spca")) {
    return(list(spca = obj, genind = NULL, nfposi = obj$nfposi, nfnega = obj$nfnega))
  }
  if (is.list(obj) && !is.null(obj$spca)) {
    nfnega <- obj$nfnega %||% obj$spca$nfnega
    return(list(
      spca = obj$spca,
      genind = obj$genind,
      nfposi = obj$nfposi %||% obj$spca$nfposi,
      nfnega = nfnega
    ))
  }
  stop("Unrecognized sPCA results RDS format: ", rds_path)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

spca_eigenvalues <- function(spca_obj) {
  if (!is.null(spca_obj$eig)) return(as.numeric(spca_obj$eig))
  if (!is.null(spca_obj$val)) return(as.numeric(spca_obj$val))
  stop("sPCA object has no eigenvalues.")
}

plot_type <- as.integer(snakemake@params[["type"]])
color_by <- as.character(snakemake@params[["color_by"]])
pc_max <- as.integer(snakemake@params[["pc_max"]])
loading_threshold <- snakemake@params[["loading_threshold"]]
plot_axis <- as.integer(snakemake@params[["plot_axis"]])
plot_loading <- as.integer(snakemake@params[["plot_loading"]])

run_main <- plot_axis == 0L && plot_loading == 0L

cat("=== sPCA Plotting (adegenet tutorial) ===\n")
cat("plot_axis:", plot_axis, " plot_loading:", plot_loading, "\n")

res <- unwrap_results(snakemake@input[["results_rds"]])
spca_obj <- res$spca
genind_obj <- res$genind
nfposi <- as.integer(res$nfposi)
nfnega <- as.integer(res$nfnega)
n_axes <- nfposi + nfnega
eig_vals <- spca_eigenvalues(spca_obj)

if (run_main) {
  rtest_bundle <- readRDS(snakemake@input[["rtest_rds"]])
  global_test <- rtest_bundle$global
  local_test <- rtest_bundle$local

  indpopdata <- read.table(
    snakemake@input[["indpopdata"]],
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
  )

  dir.create(dirname(snakemake@output[["map_plot"]]), recursive = TRUE, showWarnings = FALSE)

  cat("Eigenvalue barplot (barplot(mySpca$eig))...\n")
  n_retained <- min(nfposi + nfnega, length(eig_vals))
  bar_cols <- rep("grey", length(eig_vals))
  if (n_retained > 0) bar_cols[seq_len(n_retained)] <- "red"
  save_base_plot(
    snakemake@output[["eig_barplot"]],
    quote(barplot(eig_vals, main = "Eigenvalues of sPCA", col = bar_cols, las = 2)),
    width = 10,
    height = 6
  )

  cat("screeplot(mySpca)...\n")
  save_base_plot(
    snakemake@output[["screeplot"]],
    quote(screeplot(spca_obj)),
    width = 10,
    height = 7
  )

  cat("plot(mySpca) composite...\n")
  save_base_plot(
    snakemake@output[["composite_plot"]],
    quote(plot(spca_obj)),
    width = 12,
    height = 10
  )
  saveRDS(spca_obj, snakemake@output[["composite_plot_rds"]])

  cat("plot.spca(mySpca)...\n")
  save_base_plot(
    snakemake@output[["map_plot"]],
    quote(plot.spca(spca_obj, type = plot_type)),
    width = 10,
    height = 8
  )
  saveRDS(spca_obj, snakemake@output[["map_plot_rds"]])

  cat("colorplot (global)...\n")
  save_base_plot(
    snakemake@output[["colorplot_global"]],
    quote(colorplot(spca_obj, cex = 2, main = "sPCA colorplot — global scores")),
    width = 8,
    height = 7
  )

  if (ncol(spca_obj$li) >= 2) {
    cat("colorplot (local)...\n")
    save_base_plot(
      snakemake@output[["colorplot_local"]],
      quote(colorplot(
        spca_obj,
        axes = 1:min(2, ncol(spca_obj$li)),
        useLag = FALSE,
        main = "sPCA colorplot — local scores"
      )),
      width = 8,
      height = 7
    )
  } else {
    save_base_plot(
      snakemake@output[["colorplot_local"]],
      quote({
        plot.new()
        text(0.5, 0.5, "Not enough local axes for colorplot", cex = 1.2)
      }),
      width = 8,
      height = 7
    )
  }

  cat("plot(global.rtest) / plot(local.rtest)...\n")
  save_base_plot(
    snakemake@output[["global_rtest_plot"]],
    quote(plot(global_test)),
    width = 7,
    height = 5
  )
  save_base_plot(
    snakemake@output[["local_rtest_plot"]],
    quote(plot(local_test)),
    width = 7,
    height = 5
  )

  make_scatter <- function(scores, title_prefix, pdf_out, rds_out) {
    if (ncol(scores) < 2) {
      save_base_plot(
        pdf_out,
        quote({
          plot.new()
          text(0.5, 0.5, paste(title_prefix, ": need >= 2 axes"), cex = 1.1)
        }),
        width = 8,
        height = 6
      )
      saveRDS(list(plot = NULL, spca_obj = spca_obj), rds_out)
      return(invisible(NULL))
    }
    pc1 <- 1
    pc2 <- min(2, min(ncol(scores), pc_max))
    plot_df <- data.frame(
      Ind = rownames(scores),
      PC1 = scores[, pc1],
      PC2 = scores[, pc2],
      stringsAsFactors = FALSE
    )
    plot_df <- merge(plot_df, indpopdata, by = "Ind", all.x = TRUE)
    use_color <- !is.null(color_by) && color_by != "" && toupper(color_by) != "NONE"
    if (use_color && !(color_by %in% colnames(plot_df))) {
      use_color <- FALSE
    }
    p <- ggplot(plot_df, aes(x = PC1, y = PC2))
    if (use_color) {
      p <- p + geom_point(aes(color = .data[[color_by]]), size = 2.5, alpha = 0.85) +
        labs(color = color_by)
    } else {
      p <- p + geom_point(size = 2.5, alpha = 0.85)
    }
    p <- p +
      labs(
        title = sprintf("%s axes %d vs %d", title_prefix, pc1, pc2),
        x = sprintf("Axis %d", pc1),
        y = sprintf("Axis %d", pc2)
      ) +
      theme_bw()
    ggsave_pdf(pdf_out, plot = p, width = 8, height = 6)
    saveRDS(list(plot = p, spca_obj = spca_obj, scores = scores), rds_out)
  }

  cat("ggplot scatters...\n")
  make_scatter(
    as.matrix(spca_obj$ls),
    "Global sPCA (lagged scores)",
    snakemake@output[["scatter_global"]],
    snakemake@output[["scatter_global_rds"]]
  )
  make_scatter(
    as.matrix(spca_obj$li),
    "Local sPCA scores",
    snakemake@output[["scatter_local"]],
    snakemake@output[["scatter_local_rds"]]
  )
}

if (plot_axis > 0L) {
  axis <- plot_axis
  if (axis > n_axes) {
    stop("Requested axis ", axis, " exceeds retained axes (", n_axes, ").")
  }
  cat("plot(mySpca, axis =", axis, ")...\n")
  save_base_plot(
    snakemake@output[["axis_plot"]],
    quote(plot(spca_obj, axis = axis)),
    width = 10,
    height = 8
  )
}

if (plot_loading > 0L) {
  axis <- plot_loading
  if (axis > ncol(spca_obj$c1)) {
    stop("Requested loading axis ", axis, " exceeds ncol(c1) = ", ncol(spca_obj$c1))
  }
  loc_fac <- if (!is.null(genind_obj)) genind_obj@loc.fac else NULL
  thresh <- loading_threshold
  if (is.null(thresh) || is.na(thresh) || as.character(thresh) %in% c("NULL", "null", "NA", "")) {
    thresh <- NULL
  } else {
    thresh <- as.numeric(thresh)
  }
  loadings <- spca_obj$c1[, axis, drop = TRUE]^2
  names(loadings) <- rownames(spca_obj$c1)
  cat("loadingplot axis", axis, "...\n")
  save_base_plot(
    snakemake@output[["loading_plot"]],
    quote({
      if (is.null(thresh)) {
        loadingplot(
          loadings,
          xlab = "Alleles",
          ylab = "Squared loading",
          main = paste("Allele contributions — axis", axis),
          fac = loc_fac
        )
      } else {
        loadingplot(
          loadings,
          threshold = thresh,
          xlab = "Alleles",
          ylab = "Squared loading",
          main = paste("Allele contributions — axis", axis),
          fac = loc_fac
        )
      }
    }),
    width = 12,
    height = 6
  )
}

cat("=== sPCA Plotting Complete ===\n")

if (file.exists("Rplots.pdf")) {
  tryCatch(file.remove("Rplots.pdf"), error = function(e) invisible(NULL))
}
