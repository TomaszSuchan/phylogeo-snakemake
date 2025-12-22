library(pophelper)
library(ggplot2)
library(ggtext)
library(tidyr)

# Snakemake inputs/outputs
runs <- unlist(snakemake@input)
output_pdf <- snakemake@output$pdf
output_rds <- snakemake@output$rds

# Read STRUCTURE runs
slist <- readQ(runs, filetype="structure", indlabfromfile = TRUE)

# Get runs summary
stats_str <- summariseQ(tabulateQ(slist))

# Perform Evanno method
evanno <- evannoMethodStructure(data=stats_str, exportplot=FALSE, returnplot=FALSE, returndata=TRUE)

# Create full plot: deltaK and mean likelihood L(K)
evanno_mean <- pivot_longer(evanno[,-c(1,2,4,6,7,8,10,11,13,14)], -c(k), values_to = "Value", names_to = "Parameter")
evanno_mean$Parameter <- gsub("mean", "", evanno_mean$Parameter)

evanno_min <- pivot_longer(evanno[,c(3,7,11,14)], -c(k), values_to = "Min", names_to = "Parameter")
evanno_min$Parameter <- gsub("min", "", evanno_min$Parameter)

evanno_max <- pivot_longer(evanno[,c(3,8,10,13)], -c(k), values_to = "Max", names_to = "Parameter")
evanno_max$Parameter <- gsub("max", "", evanno_max$Parameter)

evanno_long <- merge(merge(evanno_mean, evanno_min, by = c("k", "Parameter"), all.x = TRUE), evanno_max, by = c("k", "Parameter"), all.x = TRUE)
evanno_long$Parameter <- factor(evanno_long$Parameter, levels = c("elpd", "lnk1", "lnk2", "deltaK"))

# Only keep deltaK and elpd
evanno_long <- evanno_long[evanno_long$Parameter == "deltaK" | evanno_long$Parameter == "elpd", ]

facet_labels <- c(
                `elpd` = "*L*(*K*) \u00B1 SD",
                `deltaK` = "\u0394*K* "
                )

evanno_plot <- ggplot(data=evanno_long, aes(x=k, y=Value, ymin=Min, ymax=Max)) +
               geom_line() +
               geom_point() +
               geom_errorbar(width=.3, linewidth=.3) +
               scale_x_continuous(breaks=evanno$k) +
               facet_wrap(~Parameter , scales = "free", labeller = as_labeller(facet_labels)) +
               xlab("K") +
               theme_bw() +
               theme(strip.text = ggtext::element_markdown(),
                     axis.title.x = element_text(face = "italic"))

# Save as PDF
ggsave(
  filename = output_pdf,
  plot = evanno_plot,
  width = 10,
  height = 5,
  dpi = 300
)

# Save as RDS
saveRDS(evanno_plot, file = output_rds)
