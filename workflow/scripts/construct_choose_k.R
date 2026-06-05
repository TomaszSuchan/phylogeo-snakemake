#!/usr/bin/env Rscript
# Summarise conStruct model-comparison scores across K values.

suppressPackageStartupMessages({
  library(conStruct)
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

results_files <- snakemake@input[["results_rds"]]
choose_k_file <- snakemake@output[["choose_k_results"]]
lpd_summary_file <- snakemake@output[["lpd_summary"]]
layer_summary_file <- snakemake@output[["layer_contribution_summary"]]
results_rds_out <- snakemake@output[["results_rds"]]

parse_k_from_path <- function(path) {
  k <- sub(".*\\.K([0-9]+)\\..*", "\\1", path)
  if (identical(k, path) || !nzchar(k)) {
    stop("Could not parse K from results path: ", path)
  }
  as.integer(k)
}

unwrap_results <- function(obj) {
  if (!is.null(obj$conStruct.results)) {
    return(list(results = obj$conStruct.results, data.block = obj$data.block))
  }
  list(results = obj, data.block = obj$data.block)
}

get_chain_names <- function(results) {
  chains <- grep("^chain_", names(results), value = TRUE)
  if (length(chains) > 0) {
    return(chains)
  }
  idx <- vapply(results, function(x) is.list(x) && !is.null(x$MAP), logical(1))
  names(results)[idx]
}

extract_lpd <- function(results) {
  chains <- get_chain_names(results)
  if (length(chains) == 0) {
    stop("No conStruct chain results found.")
  }
  vapply(chains, function(chain) as.numeric(results[[chain]]$MAP$lpd), numeric(1))
}

load_construct_robj <- function(path) {
  env <- new.env()
  loaded <- load(path, envir = env)
  if (length(loaded) != 1) {
    stop("Expected one object in ", path, ", found: ", paste(loaded, collapse = ", "))
  }
  env[[loaded[1]]]
}

load_data_block <- function(obj, results_path, k) {
  if (!is.null(obj$data.block)) {
    return(obj$data.block)
  }

  output_dir <- dirname(results_path)
  base_name <- sub("\\.results\\.rds$", "", basename(results_path))
  candidates <- c(
    file.path(output_dir, paste0(base_name, ".results.K", k, "_data.block.Robj")),
    file.path(output_dir, paste0(base_name, ".K", k, "_data.block.Robj")),
    file.path(output_dir, paste0(base_name, ".K", k, "_sp_data.block.Robj")),
    file.path(output_dir, paste0(base_name, "_data.block.Robj")),
    file.path(output_dir, paste0(base_name, "_sp_data.block.Robj"))
  )
  for (candidate in candidates) {
    if (file.exists(candidate)) {
      return(load_construct_robj(candidate))
    }
  }
  NULL
}

k_values <- vapply(results_files, parse_k_from_path, integer(1))
order_idx <- order(k_values)
results_files <- results_files[order_idx]
k_values <- k_values[order_idx]

lpd_rows <- list()
layer_rows <- list()
max_layers <- 0L
reference_props <- NULL

for (i in seq_along(results_files)) {
  k <- k_values[i]
  path <- results_files[i]
  cat("Processing conStruct results for K =", k, ":", path, "\n")

  obj <- readRDS(path)
  unwrapped <- unwrap_results(obj)
  results <- unwrapped$results
  data_block <- load_data_block(unwrapped, path, k)

  if (is.null(data_block)) {
    stop("data.block not found for K = ", k, ". Re-run construct_analysis with save.files = TRUE.")
  }

  lpd_values <- extract_lpd(results)
  lpd_rows[[length(lpd_rows) + 1L]] <- data.frame(
    K = k,
    median = median(lpd_values),
    min = min(lpd_values),
    max = max(lpd_values)
  )

  chains <- get_chain_names(results)
  layer_order <- NULL
  if (!is.null(reference_props)) {
    current_props <- results[[chains[1]]]$MAP$admix.proportions
    layer_order <- tryCatch(
      match.layers.x.runs(reference_props, current_props),
      error = function(e) NULL
    )
  }

  layer_contrib <- calculate.layer.contribution(
    conStruct.results = results[[chains[1]]],
    data.block = data_block,
    layer.order = layer_order
  )
  layer_contrib <- as.numeric(layer_contrib)
  max_layers <- max(max_layers, length(layer_contrib))

  for (layer_idx in seq_along(layer_contrib)) {
    layer_rows[[length(layer_rows) + 1L]] <- data.frame(
      K = k,
      Layer = layer_idx,
      Contribution = layer_contrib[layer_idx]
    )
  }

  if (is.null(reference_props)) {
    current_props <- results[[chains[1]]]$MAP$admix.proportions
    if (!is.null(layer_order)) {
      current_props <- current_props[, layer_order, drop = FALSE]
    }
    reference_props <- current_props
  }
}

lpd_summary <- do.call(rbind, lpd_rows)
layer_summary <- do.call(rbind, layer_rows)

best_k_lpd <- lpd_summary$K[which.max(lpd_summary$median)]
last_layer <- aggregate(Contribution ~ K, layer_summary, function(x) x[length(x)])
low_contrib_k <- last_layer$K[last_layer$Contribution < 0.05]
suggested_k_layer <- if (length(low_contrib_k) > 0) min(low_contrib_k) - 1L else max(k_values)
if (suggested_k_layer < min(k_values)) {
  suggested_k_layer <- min(k_values)
}

write.table(
  lpd_summary,
  file = lpd_summary_file,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
write.table(
  layer_summary,
  file = layer_summary_file,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

writeLines(
  c(
    "conStruct model comparison across K values.",
    "Higher MAP log-posterior (lpd) indicates better MCMC fit.",
    "Layer contributions near zero suggest excess model complexity; choose K below the first trivial layer.",
    paste0("Suggested K (maximum median MAP lpd): ", best_k_lpd),
    paste0("Suggested K (last layer contribution < 0.05 starts at K = ",
           if (length(low_contrib_k) > 0) min(low_contrib_k) else "NA",
           "; consider K = ", suggested_k_layer, ")"),
    "",
    "MAP log-posterior summary:",
    paste(c("K", "median", "min", "max"), collapse = "\t"),
    paste(
      lpd_summary$K,
      signif(lpd_summary$median, 6),
      signif(lpd_summary$min, 6),
      signif(lpd_summary$max, 6),
      sep = "\t"
    ),
    "",
    "Layer contributions (relative):",
    paste(c("K", "Layer", "Contribution"), collapse = "\t"),
    paste(
      layer_summary$K,
      layer_summary$Layer,
      signif(layer_summary$Contribution, 6),
      sep = "\t"
    )
  ),
  con = choose_k_file
)

saveRDS(
  list(
    lpd_summary = lpd_summary,
    layer_summary = layer_summary,
    suggested_k_lpd = best_k_lpd,
    suggested_k_layer = suggested_k_layer,
    k_values = k_values
  ),
  file = results_rds_out
)

cat("Wrote choose-K summary:", choose_k_file, "\n")
cat("Wrote MAP lpd summary:", lpd_summary_file, "\n")
cat("Wrote layer contribution summary:", layer_summary_file, "\n")
