#!/usr/bin/env Rscript
# figS_immune_exclusion.R — figure for I2

suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
source("scripts/figures/publication_theme.R")
source("scripts/figures/load_figure_data.R")

IN <- file.path(BASE, "Analysis/Spatial/results/immune_exclusion")
per_sample <- fread(file.path(IN, "per_sample_exclusion_summary.csv"))
cooc  <- fread(file.path(IN, "celltype_cooccurrence_spearman.csv"))
by_cls <- fread(file.path(IN, "celltype_mean_by_spot_class.csv"))

# Panel A: per-sample immune vs stromal scatter
pA <- ggplot(per_sample, aes(mean_immune, mean_stromal, label = sample_id)) +
  geom_point(size = 3, color = "#4472C4") +
  geom_text(hjust = -0.3, size = 2.5) +
  labs(x = "Mean immune score (per sample)",
       y = "Mean stromal score (per sample)",
       title = "Per-sample immune vs. stromal balance") +
  theme_masld()

# Panel B: co-occurrence heatmap
setnames(cooc, names(cooc), c("ct1","ct2","rho"))
cooc <- cooc[!is.na(rho) & ct1 != ct2]
# Keep top cell types
top_cts <- unique(c(cooc[, sum(abs(rho), na.rm = TRUE), by = ct1][order(-V1)][1:12, ct1],
                    cooc[, sum(abs(rho), na.rm = TRUE), by = ct2][order(-V1)][1:12, ct2]))
sub <- cooc[ct1 %in% top_cts & ct2 %in% top_cts]

pB <- ggplot(sub, aes(ct1, ct2, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 2.2) +
  scale_fill_gradient2(low = "#2980B9", mid = "grey93", high = "#C0392B",
                       midpoint = 0, limits = c(-1, 1), name = "Spearman rho") +
  labs(x = NULL, y = NULL,
       title = "Cell-type spatial co-occurrence (cell2location)") +
  theme_masld() +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 6.5),
        axis.text.y = element_text(size = 6.5))

# Panel C: mean abundance by spot class
by_cls_long <- melt(by_cls, id.vars = "spot_class",
                     variable.name = "celltype", value.name = "mean_abundance")
by_cls_long[, celltype := as.character(celltype)]
pC <- ggplot(by_cls_long, aes(celltype, mean_abundance, fill = spot_class)) +
  geom_col(position = "dodge", width = 0.7, color = "white") +
  scale_fill_brewer(palette = "Set1", name = NULL) +
  labs(x = NULL, y = "Mean abundance",
       title = "Cell-type abundance by spot/sample class") +
  theme_masld() +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 6.5),
        legend.position = "bottom",
        legend.key.size = unit(3, "mm"),
        legend.text = element_text(size = 6))

fig <- (pA | pC) / pB + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 8, face = "bold"))

ggsave(file.path(FIGS_CELLTYPE_DIR, "figS_I2_immune_exclusion.pdf"),
       fig, width = 14, height = 12)
message("Saved I2 figure")
