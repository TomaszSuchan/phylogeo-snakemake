# Shared ROH length-class definitions.
#
# Physical length (Mb) is converted to generations since common ancestry using
# G = 100/(2 × cM), with cM estimated from megabases assuming a uniform
# recombination rate of 1 cM/Mb.

ROH_RECOMBINATION_RATE_CM_PER_MB <- 1.0

ROH_CLASS_LEVELS <- c(
  "Short (<1 Mb; >50 gen)",
  "Medium (1-5 Mb; 10-50 gen)",
  "Long (>5 Mb; <10 gen)"
)

assign_roh_class <- function(length_mb) {
  dplyr::case_when(
    length_mb < 1 ~ ROH_CLASS_LEVELS[1],
    length_mb >= 1 & length_mb < 5 ~ ROH_CLASS_LEVELS[2],
    length_mb >= 5 ~ ROH_CLASS_LEVELS[3],
    TRUE ~ NA_character_
  )
}

plot_group_utils <- tryCatch({
  file.path(dirname(normalizePath(snakemake@script)), "plot_group_utils.R")
}, error = function(e) "workflow/scripts/plot_group_utils.R")
if (file.exists(plot_group_utils)) {
  source(plot_group_utils)
} else {
  source("workflow/scripts/plot_group_utils.R")
}

roh_group_fill_values <- group_fill_values
roh_group_sort_by <- group_sort_by
roh_group_levels <- group_levels
