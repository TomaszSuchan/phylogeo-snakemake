# Shared ROH length-class definitions and helpers.
#
# Physical length (Mb) is converted to generations since common ancestry using
# G = 100/(2 × cM), with cM estimated from megabases assuming a uniform
# recombination rate of 1 cM/Mb.

ROH_RECOMBINATION_RATE_CM_PER_MB <- 1.0

ROH_CLASS_LEVELS <- c(
  "Short (<1 Mb; >50 gen)",
  "Medium (1-5 Mb; 10-50 gen)",
  "Long (>5 Mb; ≤10 gen)"
)

ROH_CLASS_FILL <- c(
  "Short (<1 Mb; >50 gen)" = "lightblue",
  "Medium (1-5 Mb; 10-50 gen)" = "lightgreen",
  "Long (>5 Mb; ≤10 gen)" = "salmon"
)

ROH_LENGTH_CLASS_CAPTION <- paste0(
  "ROH length classes approximate generations since common ancestry as ",
  "G = 100/(2 × cM), with cM estimated from physical length assuming ",
  ROH_RECOMBINATION_RATE_CM_PER_MB, " cM/Mb ",
  "(>5 Mb: ≤10 gen; 1-5 Mb: 10-50 gen; <1 Mb: >50 gen)."
)

assign_roh_class <- function(length_mb) {
  dplyr::case_when(
    length_mb < 1 ~ ROH_CLASS_LEVELS[1],
    length_mb >= 1 & length_mb < 5 ~ ROH_CLASS_LEVELS[2],
    length_mb >= 5 ~ ROH_CLASS_LEVELS[3],
    TRUE ~ NA_character_
  )
}
