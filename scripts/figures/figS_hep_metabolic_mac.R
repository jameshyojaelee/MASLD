#!/usr/bin/env Rscript
# figS_hep_metabolic_mac.R
# Figure for Analysis E1 — Hep metabolic state × macrophage CCC.
#
# Panels:
#   A — Heatmap: meta-subtype × metabolic module (FAO/Lipogenic/etc.)
#   B — Top Hep->Mac metabolic-module ligands with LIANA score_diff
#   C — Per meta-subtype top Hep->Mac ligands (faceted bars)
#
# Output: figures/supplementary/figS_celltype_biology/figS_E1_hep_metabolic_mac.pdf

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})

source("scripts/figures/publication_theme.R")
source("scripts/figures/load_figure_data.R")

IN <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/metabolic_ccc")
mod_w <- fread(file.path(IN, "hep_metabolic_module_scores_wide.csv"))
lig_m <- fread(file.path(IN, "hep_to_mac_ligands_by_meta.csv"))
met_l <- fread(file.path(IN, "hep_to_mac_metabolic_ligand_axes.csv"))

# Panel A heatmap
long <- melt(mod_w, id.vars = "meta_subtype", variable.name = "module",
             value.name = "mean_lfc")
long[is.na(mean_lfc), mean_lfc := 0]

pA <- ggplot(long, aes(module, meta_subtype, fill = mean_lfc)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", mean_lfc)), size = 2.2) +
  scale_fill_gradient2(low = "#2980B9", mid = "grey93", high = "#C0392B",
                       midpoint = 0, name = "mean_lfc\n(wilcoxon)") +
  labs(x = NULL, y = NULL,
       title = "Hep meta-subtypes × metabolic modules") +
  theme_masld() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 7))

# Panel B: metabolic-module Hep->Mac ligands
metab_valid <- met_l[!is.na(score_diff) & !is.na(module)][order(-score_diff)][1:min(20, .N)]
if (nrow(metab_valid) > 0) {
  metab_valid[, pair_label := paste0(ligand, " -> ", receptor, " [", module, "]")]
  metab_valid[, pair_label := factor(pair_label, levels = rev(unique(pair_label)))]
  pB <- ggplot(metab_valid, aes(score_diff, pair_label, fill = module)) +
    geom_col(width = 0.7, color = "white") +
    geom_text(aes(label = sprintf("%.2f", score_diff)), hjust = -0.1, size = 2) +
    scale_fill_brewer(palette = "Set2", name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
    labs(x = "LIANA score_diff", y = NULL,
         title = "Metabolic-module Hep->Mac ligands (LIANA MASLD-enriched)") +
    theme_masld() +
    theme(axis.text.y = element_text(size = 5.8),
          legend.position = "bottom",
          legend.key.size = unit(3, "mm"),
          legend.text = element_text(size = 6))
} else {
  pB <- ggplot() + labs(title = "No metabolic-module ligand overlap") + theme_masld()
}

# Panel C: per meta-subtype top Hep->Mac ligands
top_per <- lig_m[, .SD[1:min(.N, 6)], by = meta_subtype]
top_per[, label := paste0(meta_subtype, ": ", gene)]
top_per[, label := factor(label, levels = rev(unique(label)))]
pC <- ggplot(top_per, aes(meta_lfc, label, fill = meta_subtype)) +
  geom_col(width = 0.7, color = "white") +
  facet_wrap(~ meta_subtype, scales = "free_y", ncol = 1) +
  scale_fill_brewer(palette = "Set1", guide = "none") +
  labs(x = "Mean logfoldchange (subtype markers vs rest)",
       y = NULL,
       title = "Top Hep->Mac ligands per Hep meta-subtype") +
  theme_masld() +
  theme(axis.text.y = element_text(size = 5.8),
        strip.text = element_text(size = 6.5))

fig <- (pA / pB) | pC
fig <- fig + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 8, face = "bold"))

out_path <- file.path(FIGS_CELLTYPE_DIR, "figS_E1_hep_metabolic_mac.pdf")
ggsave(out_path, fig, width = 14, height = 12)
message("Saved: ", out_path)
