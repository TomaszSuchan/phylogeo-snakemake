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

ROH_LENGTH_CLASS_CAPTION <- paste0(
  "ROH length classes approximate generations since common ancestry as ",
  "G = 100/(2 × cM), with cM estimated from physical length assuming ",
  ROH_RECOMBINATION_RATE_CM_PER_MB, " cM/Mb."
)

assign_roh_class <- function(length_mb) {
  dplyr::case_when(
    length_mb < 1 ~ ROH_CLASS_LEVELS[1],
    length_mb >= 1 & length_mb < 5 ~ ROH_CLASS_LEVELS[2],
    length_mb >= 5 ~ ROH_CLASS_LEVELS[3],
    TRUE ~ NA_character_
  )
}

# Named level -> color map from parameters.roh.colors.{group_col}, or NULL.
roh_group_fill_values <- function(group_colors_param) {
  if (is.null(group_colors_param)) {
    return(NULL)
  }
  colors <- unlist(group_colors_param, use.names = TRUE)
  if (length(colors) == 0 || all(names(colors) == "")) {
    return(NULL)
  }
  colors
}
