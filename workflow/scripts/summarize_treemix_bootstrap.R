#!/usr/bin/env Rscript

suppressPackageStartupMessages({
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

pdf(NULL)

log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")
on.exit({
  while (sink.number(type = "message") > 0) sink(type = "message")
  while (sink.number(type = "output") > 0) sink(type = "output")
  close(log_file)
}, add = TRUE)

extract_final_llik <- function(text) {
  matches <- regmatches(text, gregexpr("Exiting ln\\(likelihood\\) [-0-9.eE]+", text))[[1]]
  if (length(matches) == 0) {
    return(NA_real_)
  }
  as.numeric(sub(".* ", "", tail(matches, 1)))
}

extract_final_m <- function(text) {
  matches <- regmatches(text, gregexpr("Exiting ln\\(likelihood\\) [-0-9.eE]+ with [0-9]+ migration events", text))[[1]]
  if (length(matches) == 0) {
    return(NA_integer_)
  }
  as.integer(sub(".* with ([0-9]+) migration events", "\\1", tail(matches, 1)))
}

parse_run_id <- function(path) {
  base <- basename(path)
  data.frame(
    replicate = as.integer(sub(".*\\.r([0-9]+)\\.m[0-9]+\\..*$", "\\1", base)),
    requested_m = as.integer(sub(".*\\.m([0-9]+)\\..*$", "\\1", base)),
    stringsAsFactors = FALSE
  )
}

read_vertices <- function(path) {
  lines <- readLines(gzfile(path), warn = FALSE)
  if (length(lines) == 0) {
    return(data.frame())
  }
  rows <- lapply(lines, function(line) {
    parts <- strsplit(line, "\\s+")[[1]]
    if (length(parts) < 10) {
      parts <- c(parts, rep(NA_character_, 10 - length(parts)))
    }
    data.frame(
      node = parts[1],
      label = parts[2],
      root = parts[3],
      mig = parts[4],
      tip = parts[5],
      parent = parts[6],
      child1 = parts[7],
      child2 = parts[8],
      child3 = parts[9],
      subtree = paste(parts[10:length(parts)], collapse = " "),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

tip_clade_label <- function(vertex) {
  if (nrow(vertex) == 0) {
    return(NA_character_)
  }
  if (!is.na(vertex$label) && vertex$label != "NA" && vertex$label != "") {
    return(vertex$label)
  }
  subtree <- vertex$subtree
  tips <- unique(regmatches(subtree, gregexpr("[A-Za-z0-9_.-]+(?=:[0-9.eE+-]+)", subtree, perl = TRUE))[[1]])
  tips <- tips[!is.na(tips) & tips != "" & tips != "NA"]
  if (length(tips) > 0) {
    return(paste(sort(tips), collapse = "|"))
  }
  paste0("node_", vertex$node)
}

read_migration_edges <- function(edges_path, vertices_path) {
  run_id <- parse_run_id(edges_path)
  vertices <- read_vertices(vertices_path)
  edges <- read.table(
    gzfile(edges_path),
    header = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )
  if (nrow(edges) == 0 || ncol(edges) < 5) {
    return(data.frame())
  }
  names(edges)[1:5] <- c("source_node", "target_node", "length", "weight", "edge_type")
  edges <- edges[edges$edge_type == "MIG", , drop = FALSE]
  if (nrow(edges) == 0) {
    return(data.frame())
  }
  edge_rows <- lapply(seq_len(nrow(edges)), function(i) {
    source_vertex <- vertices[vertices$node == as.character(edges$source_node[i]), , drop = FALSE]
    target_vertex <- vertices[vertices$node == as.character(edges$target_node[i]), , drop = FALSE]
    data.frame(
      requested_m = run_id$requested_m,
      replicate = run_id$replicate,
      source_node = as.character(edges$source_node[i]),
      target_node = as.character(edges$target_node[i]),
      source_clade = tip_clade_label(source_vertex),
      target_clade = tip_clade_label(target_vertex),
      weight = suppressWarnings(as.numeric(edges$weight[i])),
      edge_file = edges_path,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, edge_rows)
}

summarize_support <- function(edge_rows, replicate_summary) {
  if (nrow(edge_rows) == 0) {
    return(data.frame(
      requested_m = integer(),
      source_clade = character(),
      target_clade = character(),
      n_replicates_with_edge = integer(),
      n_replicates_total = integer(),
      support = numeric(),
      mean_weight = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  totals <- aggregate(replicate ~ requested_m, replicate_summary, function(x) length(unique(x)))
  names(totals)[2] <- "n_replicates_total"
  counts <- aggregate(
    replicate ~ requested_m + source_clade + target_clade,
    edge_rows,
    function(x) length(unique(x))
  )
  names(counts)[4] <- "n_replicates_with_edge"
  weights <- aggregate(
    weight ~ requested_m + source_clade + target_clade,
    edge_rows,
    function(x) mean(x, na.rm = TRUE)
  )
  names(weights)[4] <- "mean_weight"
  out <- merge(counts, totals, by = "requested_m", all.x = TRUE)
  out <- merge(out, weights, by = c("requested_m", "source_clade", "target_clade"), all.x = TRUE)
  out$support <- out$n_replicates_with_edge / out$n_replicates_total
  out[order(out$requested_m, -out$support, out$source_clade, out$target_clade), ]
}

llik_files <- unlist(snakemake@input[["llik"]])
edge_files <- unlist(snakemake@input[["edges"]])
vertex_files <- unlist(snakemake@input[["vertices"]])

cat("Summarizing TreeMix bootstrap replicates.\n")
cat("Likelihood files:", length(llik_files), "\n")

replicate_rows <- lapply(llik_files, function(path) {
  run_id <- parse_run_id(path)
  text <- paste(readLines(path, warn = FALSE), collapse = " ")
  data.frame(
    requested_m = run_id$requested_m,
    replicate = run_id$replicate,
    achieved_m = extract_final_m(text),
    final_log_likelihood = extract_final_llik(text),
    llik_file = path,
    stringsAsFactors = FALSE
  )
})
replicate_summary <- do.call(rbind, replicate_rows)
replicate_summary <- replicate_summary[order(replicate_summary$requested_m, replicate_summary$replicate), ]

edge_rows <- lapply(seq_along(edge_files), function(i) read_migration_edges(edge_files[i], vertex_files[i]))
edge_rows <- Filter(function(x) nrow(x) > 0, edge_rows)
if (length(edge_rows) > 0) {
  edge_rows <- do.call(rbind, edge_rows)
} else {
  edge_rows <- data.frame(
    requested_m = integer(),
    replicate = integer(),
    source_node = character(),
    target_node = character(),
    source_clade = character(),
    target_clade = character(),
    weight = numeric(),
    edge_file = character(),
    stringsAsFactors = FALSE
  )
}
support_summary <- summarize_support(edge_rows, replicate_summary)

dir.create(dirname(snakemake@output[["summary"]]), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(snakemake@output[["pdf"]]), recursive = TRUE, showWarnings = FALSE)
write.table(replicate_summary, snakemake@output[["summary"]], sep = "\t", row.names = FALSE, quote = FALSE)
write.table(support_summary, snakemake@output[["migration_edges"]], sep = "\t", row.names = FALSE, quote = FALSE)

plot_obj <- ggplot2::ggplot(
  replicate_summary,
  ggplot2::aes(x = factor(requested_m), y = final_log_likelihood)
) +
  ggplot2::geom_boxplot(outlier.shape = NA, width = 0.5) +
  ggplot2::geom_jitter(width = 0.12, height = 0, size = 1.5, alpha = 0.7) +
  ggplot2::labs(
    x = "Requested migration edges (m)",
    y = "Final log likelihood"
  ) +
  ggplot2::theme_minimal(base_size = 11)

ggsave_pdf(
  snakemake@output[["pdf"]],
  plot_obj,
  width = as.numeric(snakemake@params[["width"]]),
  height = as.numeric(snakemake@params[["height"]]),
  dpi = as.numeric(snakemake@params[["dpi"]])
)

saveRDS(
  list(
    replicate_summary = replicate_summary,
    migration_edges = edge_rows,
    support_summary = support_summary,
    plot = plot_obj
  ),
  snakemake@output[["rds"]]
)

cat("Bootstrap summary:", snakemake@output[["summary"]], "\n")
cat("Bootstrap migration-edge support:", snakemake@output[["migration_edges"]], "\n")
