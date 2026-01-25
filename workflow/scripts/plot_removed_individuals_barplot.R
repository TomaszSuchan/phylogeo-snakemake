#!/usr/bin/env Rscript
# Create barplot showing counts of removed individuals by relationship category,
# grouped by Site (population) with patterns instead of colors

library(ggplot2)
library(dplyr)

# Prevent creation of Rplots.pdf
pdf(NULL)

# Redirect all output to log file
log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

# Snakemake inputs/outputs
removed_file <- snakemake@input[["removed_individuals"]]
popdata_file <- snakemake@input[["popdata"]]
output_pdf <- snakemake@output[["pdf"]]
output_rds <- snakemake@output[["rds"]]

# Parameters
group_by_name <- NULL
if ("group_by" %in% names(snakemake@params)) {
  group_by_param <- snakemake@params[["group_by"]]
  if (!is.null(group_by_param) && !is.na(group_by_param)) {
    group_by_name <- as.character(group_by_param)
    if (length(group_by_name) == 0 || group_by_name == "" || group_by_name == "none") {
      group_by_name <- NULL
    }
  } else {
    group_by_name <- NULL
  }
}

message("\n=== READING REMOVED INDIVIDUALS ===\n")
removed_df <- read.table(removed_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
message(sprintf("Loaded %d removed individuals\n", nrow(removed_df)))

# If no removed individuals, create empty plot
if (nrow(removed_df) == 0 || (nrow(removed_df) == 1 && removed_df$individual[1] == "")) {
  message("No removed individuals found. Creating empty plot.\n")
  p <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = "No removed individuals", size = 6) +
    theme_void()
  
  dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)
  ggsave(
    filename = output_pdf,
    plot = p,
    width = 8,
    height = 6,
    dpi = 300,
    device = "pdf",
    limitsize = FALSE
  )
  
  dir.create(dirname(output_rds), recursive = TRUE, showWarnings = FALSE)
  saveRDS(p, output_rds)
  
  message("\n=== COMPLETED SUCCESSFULLY ===\n")
  quit(status = 0)
}

# Read popdata to get Site information
message("\n=== READING POPDATA ===\n")
popdata <- read.table(popdata_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
message(sprintf("Popdata columns: %s\n", paste(names(popdata), collapse = ", ")))

# Get Individual column name (indpopdata uses "Ind")
if ("Ind" %in% colnames(popdata)) {
  ind_col <- "Ind"
} else if ("Individual" %in% colnames(popdata)) {
  ind_col <- "Individual"
} else if ("Site" %in% colnames(popdata)) {
  ind_col <- "Site"
} else {
  ind_col <- colnames(popdata)[1]
}

# Merge removed individuals with popdata to get Site
removed_df_merged <- merge(removed_df, popdata[, c(ind_col, "Site"), drop = FALSE],
                           by.x = "individual",
                           by.y = ind_col,
                           all.x = TRUE)

# Remove individuals without Site data
na_mask <- is.na(removed_df_merged$Site)
if (any(na_mask)) {
  n_na <- sum(na_mask)
  warning(sprintf("%d individuals have missing Site values and will be excluded\n", n_na))
  removed_df_merged <- removed_df_merged[!na_mask, ]
}

# Remove empty strings
empty_mask <- removed_df_merged$Site == ""
if (any(empty_mask)) {
  n_empty <- sum(empty_mask)
  warning(sprintf("%d individuals have empty Site strings and will be excluded\n", n_empty))
  removed_df_merged <- removed_df_merged[!empty_mask, ]
}

message(sprintf("After merging with popdata: %d removed individuals with Site information\n", nrow(removed_df_merged)))

# Count removed individuals by Site and category
removed_count_df <- removed_df_merged %>%
  group_by(Site, category) %>%
  summarise(count = n(), .groups = "drop") %>%
  as.data.frame()

# Get total individuals per Site from popdata
total_count_df <- popdata %>%
  group_by(Site) %>%
  summarise(total = n(), .groups = "drop") %>%
  as.data.frame()

# Calculate kept individuals per Site (total - removed)
removed_total_df <- removed_df_merged %>%
  group_by(Site) %>%
  summarise(removed_total = n(), .groups = "drop") %>%
  as.data.frame()

kept_count_df <- merge(total_count_df, removed_total_df, by = "Site", all.x = TRUE)
kept_count_df$removed_total[is.na(kept_count_df$removed_total)] <- 0
kept_count_df$kept <- kept_count_df$total - kept_count_df$removed_total
kept_count_df <- kept_count_df[, c("Site", "kept")]
kept_count_df$category <- "kept"
kept_count_df <- kept_count_df[, c("Site", "category", "kept")]
colnames(kept_count_df)[3] <- "count"

# Combine removed and kept counts
count_df <- rbind(removed_count_df, kept_count_df)

# Order categories (kept should be last so it appears on top when stacked)
category_order <- c("clone", "1st-degree", "2nd-degree", "other", "kept")
count_df$category <- factor(count_df$category, levels = category_order)

# Order sites by total individuals (descending) - use total_count_df which has all individuals
site_totals <- total_count_df %>%
  arrange(desc(total))
count_df$Site <- factor(count_df$Site, levels = site_totals$Site)

message("\n=== COUNTS BY SITE AND CATEGORY ===\n")
print(count_df)

# If grouping by a column (e.g., Region), add that information
if (!is.null(group_by_name) && group_by_name %in% colnames(popdata)) {
  # Get unique Site-Region mapping
  site_region <- popdata[, c("Site", group_by_name)] %>%
    distinct()
  
  # Merge with count_df
  count_df <- merge(count_df, site_region, by = "Site", all.x = TRUE)
  
  # Order sites within each region by total individuals
  site_region_totals <- merge(total_count_df, site_region, by = "Site", all.x = TRUE)
  site_totals_by_region <- site_region_totals %>%
    arrange(.data[[group_by_name]], desc(total))
  count_df$Site <- factor(count_df$Site, levels = site_totals_by_region$Site)
  
  # Create facet plot grouped by region with stacked bars
  p <- ggplot(count_df, aes(x = Site, y = count, fill = category)) +
    geom_bar(
      stat = "identity", 
      position = "stack", 
      color = "black", 
      linewidth = 0.3
    ) +
    facet_wrap(as.formula(paste("~", group_by_name)), scales = "free_x", ncol = 4) +
    labs(
      x = "Population",
      y = "Number of Individuals",
      fill = "Category"
    ) +
    scale_fill_manual(
      values = c("clone" = "black",
                 "1st-degree" = "gray30",
                 "2nd-degree" = "gray60",
                 "other" = "gray85",
                 "kept" = "white"),
      breaks = category_order,
      labels = c("clone" = "Clone", "1st-degree" = "1st-degree", 
                 "2nd-degree" = "2nd-degree", "other" = "Other", "kept" = "Kept")
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      panel.grid.major.x = element_blank(),
      legend.position = "bottom",
      strip.text = element_text(size = 10, face = "bold")
    )
} else {
  # No grouping - show all sites with stacked bars
  p <- ggplot(count_df, aes(x = Site, y = count, fill = category)) +
    geom_bar(
      stat = "identity", 
      position = "stack", 
      color = "black", 
      linewidth = 0.3
    ) +
    labs(
      x = "Population",
      y = "Number of Individuals",
      fill = "Category"
    ) +
    scale_fill_manual(
      values = c("clone" = "black",
                 "1st-degree" = "gray30",
                 "2nd-degree" = "gray60",
                 "other" = "gray85",
                 "kept" = "white"),
      breaks = category_order,
      labels = c("clone" = "Clone", "1st-degree" = "1st-degree", 
                 "2nd-degree" = "2nd-degree", "other" = "Other", "kept" = "Kept")
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      panel.grid.major.x = element_blank(),
      legend.position = "bottom"
    )
}

# Save plots
message("\n=== SAVING OUTPUT ===\n")
message(sprintf("Output PDF: %s\n", output_pdf))
message(sprintf("Output RDS: %s\n", output_rds))

dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)
ggsave(
  filename = output_pdf,
  plot = p,
  width = max(12, length(unique(count_df$Site)) * 0.4),
  height = 6,
  dpi = 300,
  device = "pdf",
  limitsize = FALSE
)

dir.create(dirname(output_rds), recursive = TRUE, showWarnings = FALSE)
saveRDS(p, output_rds)

message("\n=== COMPLETED SUCCESSFULLY ===\n")
