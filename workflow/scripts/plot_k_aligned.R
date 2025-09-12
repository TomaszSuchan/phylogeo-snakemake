library(pophelper)
library(gridExtra)
library(ggplot2)

# Snakemake inputs/outputs
runs <- unlist(snakemake@input)
print(runs)
output_plot <- snakemake@output[[1]]

# Read runs
slist <- readQ(runs)
print(slist)
# Align clusters
slist_aligned <- alignK(slist)

# Plot
p <- plotQ(slist_aligned, imgoutput="join", returnplot=TRUE, exportplot=FALSE, basesize=11)

# Save the plot
n_reps <- length(slist_aligned)
print(n_reps)
n_inds <- nrow(slist_aligned[[1]])
print(n_inds)

height <- max(1.2, 1.2 * n_reps)  # 1.2 inch per replicate
width <- max(8, 1 + 0.01 * n_inds)  # at least 8 inches, but scaling with individuals

ggsave(
  filename = output_plot,
  plot = p$plot[[1]],
  width = width,
  height = height,
  dpi = 300
)