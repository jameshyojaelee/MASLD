#!/usr/bin/env Rscript
# figS_crossspecies_ccc.R
# Figure for Analysis L3 — Cross-species CCC conservation.
#
# Panels:
#   A — Stacked bar of axis conservation across all LR pairs + MASLD-enriched
#   B — Top 20 fully-conserved MASLD-enriched CCC axes
#   C — Species-divergent CCC axes (caution zone: mouse model may not recap)
#   D — Cell-type pair distribution of fully-conserved axes
#
# Output: figures/supplementary/figS_celltype_biology/figS_L3_crossspecies_ccc.pdf

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})

source("scripts/figures/publication_theme.R")
source("scripts/figures/load_figure_data.R")

IN <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/crossspecies_ccc")
pri <- fread(file.path(IN, "conserved_masld_ccc_priority.csv"))
div <- fread(file.path(IN, "divergent_masld_ccc.csv"))
ann <- fread(file.path(IN, "liana_crossspecies_annotated.csv"))

# Panel A — conservation distribution, all vs MASLD-enriched
all_cnt <- ann[, .(set = "All LR pairs", N = .N), by = axis_conservation]
masld_cnt <- ann[abs(score_diff) > 0.1, .(set = "MASLD-enriched (|score_diff|>0.1)", N = .N),
                 by = axis_conservation]
bar <- rbind(all_cnt, masld_cnt)
bar[, axis_conservation := factor(axis_conservation,
    levels = c("Fully_Conserved","Partially_Conserved","Other","Divergent","Unclassified"))]
bar[, set := factor(set, levels = unique(bar$set))]

pA <- ggplot(bar, aes(set, N, fill = axis_conservation)) +
  geom_col(position = "fill", color = "white") +
  scale_fill_manual(values = c("Fully_Conserved" = "#27AE60",
                               "Partially_Conserved" = "#2ECC71",
                               "Other" = "#95A5A6",
                               "Divergent" = "#E74C3C",
                               "Unclassified" = "grey80"),
                    name = NULL) +
  scale_y_continuous(labels = scales::percent_format(),
                     expand = c(0, 0)) +
  labs(x = NULL, y = "Proportion") +
  theme_masld() +
  theme(legend.key.size = unit(3, "mm"),
        legend.text = element_text(size = 6),
        axis.text.x = element_text(size = 6))

# Panel B — top 20 fully-conserved MASLD-enriched
topB <- pri[1:min(20, nrow(pri))]
topB[, pair_label := paste0(source, " -> ", target, ": ",
                             ligand_complex, " -> ", receptor_complex)]
topB[, pair_label := factor(pair_label, levels = rev(unique(pair_label)))]
pB <- ggplot(topB, aes(score_diff, pair_label)) +
  geom_col(fill = "#27AE60", width = 0.7, color = "white") +
  geom_text(aes(label = sprintf("%.3f", score_diff)),
            hjust = -0.1, size = GEOM_TEXT_6PT) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(x = "LIANA score_diff", y = NULL) +
  theme_masld() +
  theme(axis.text.y = element_text(size = 6))

# Panel C — divergent (top by |score_diff|)
topC <- div[1:min(15, nrow(div))]
topC[, pair_label := paste0(source, " -> ", target, ": ",
                             ligand_complex, " -> ", receptor_complex)]
topC[, pair_label := factor(pair_label, levels = rev(unique(pair_label)))]
pC <- ggplot(topC, aes(score_diff, pair_label)) +
  geom_col(fill = "#E74C3C", width = 0.7, color = "white") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey50") +
  scale_x_continuous(expand = expansion(mult = 0.2)) +
  labs(x = "LIANA score_diff", y = NULL) +
  theme_masld() +
  theme(axis.text.y = element_text(size = 6))

# Panel D — cell-type pair distribution
ct_pri <- pri[, .N, by = .(source, target)][order(-N)][1:14]
ct_pri[, pair := paste0(source, " -> ", target)]
ct_pri[, pair := factor(pair, levels = rev(pair))]
pD <- ggplot(ct_pri, aes(N, pair)) +
  geom_col(fill = "#27AE60", width = 0.65, color = "white") +
  geom_text(aes(label = N), hjust = -0.2, size = GEOM_TEXT_6PT) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Fully-conserved MASLD-enriched LR pairs", y = NULL) +
  theme_masld() +
  theme(axis.text.y = element_text(size = 6))

fig <- (pA | pD) / (pB | pC) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 6, face = "plain"))

out_path <- file.path(FIGS_CELLTYPE_DIR, "figS_L3_crossspecies_ccc.pdf")
ggsave(out_path, fig, width = fig_full_width, height = fig_full_width * (12 / 14))
message("[caption] A: Axis conservation (LIANA LR pairs vs cross-species atlas). B: Top 20 fully-conserved MASLD-enriched CCC axes (of ", nrow(pri), "). C: Species-divergent CCC axes (mouse model caution). D: Cell-type pairs with most conserved CCC.")
message("Saved: ", out_path)
