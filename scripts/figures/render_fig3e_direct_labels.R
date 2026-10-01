#!/usr/bin/env Rscript
# KEY MESSAGE: F0-reference DEG counts and continuous NMF axes describe stage-associated remodeling.
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(ggrepel); library(patchwork)
})
root <- Sys.getenv("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
source(file.path(root, "scripts/figures/publication_theme.R"))
input <- file.path(root, "figures/main/fig3_bulk_transcriptomics/source_tables")
output <- commandArgs(trailingOnly = TRUE)[1]
stopifnot(!is.na(output), !file.exists(output))
counts <- fread(file.path(input, "fig3e_deg_counts.tsv"))[axis == "fibrosis"]
activity <- fread(file.path(input, "fig3e_prespecified_nmf_activity.tsv"))[axis == "fibrosis"]
stages <- paste0("F", 0:4)
counts[, display := factor(display, levels = stages)]
activity[, group := factor(group, levels = stages)]
colors <- setNames(cat_palette[seq_len(6)], paste0("P", 1:6))
ends <- activity[group == "F4"]
stopifnot(nrow(ends) == 6L, uniqueN(ends$program) == 6L)
bars <- ggplot(counts, aes(display, signed_n, fill = direction)) +
  geom_col(width = .44) + geom_hline(yintercept = 0, linewidth = .25) +
  scale_fill_manual(values = c(Higher = masld_colors$mash, Lower = masld_colors$down), guide = "none") +
  scale_x_discrete(drop = FALSE, expand = expansion(add = c(.6, .9))) +
  scale_y_continuous(labels = abs, expand = expansion(mult = c(.13, .16))) +
  labs(x = "Fibrosis stage", y = "DEGs vs F0") + theme_masld_compact() +
  theme(plot.margin = margin(2, 2, 1, 2))
lines <- ggplot(activity, aes(group, centered_activity, group = program, color = program)) +
  geom_hline(yintercept = 0, color = "grey75", linewidth = .25) +
  geom_line(linewidth = .55) + geom_point(size = 1.05) +
  geom_text_repel(data = ends, aes(label = program), family = "Helvetica",
                  size = 6 / .pt, fontface = "plain", color = "black", nudge_x = .32,
                  direction = "y", hjust = 0, box.padding = .12,
                  point.padding = .1, min.segment.length = 0,
                  segment.size = .2, seed = 20260930, max.overlaps = Inf,
                  show.legend = FALSE) +
  scale_color_manual(values = colors, guide = "none") +
  scale_x_discrete(drop = FALSE, expand = expansion(add = c(.6, .9))) +
  labs(x = "Fibrosis stage", y = "NMF activity") + theme_masld_compact() +
  theme(legend.position = "none", plot.margin = margin(1, 2, 2, 2))
panel <- bars / lines + plot_layout(heights = c(1, 1.05))
ggsave(output, panel, width = 2.65, height = 2.27, device = cairo_pdf)
writeLines(capture.output(sessionInfo()), paste0(output, ".sessionInfo.txt"))
