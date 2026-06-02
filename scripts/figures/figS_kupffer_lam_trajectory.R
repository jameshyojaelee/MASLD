#!/usr/bin/env Rscript
# figS_kupffer_lam_trajectory.R — figure for D1

suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
source("scripts/figures/publication_theme.R")
source("scripts/figures/load_figure_data.R")

IN <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/macrophage_trajectory")
pan <- fread(file.path(IN, "macrophage_marker_panels_pseudotime.csv"))
ax  <- fread(file.path(IN, "hep_to_mac_ligands_x_pseudotime_receptors.csv"))

# Panel A: marker panel box of Spearman rho
pan[, panel := factor(panel, levels = c("Kupffer","M2_reparative","M1_inflammatory","LAM"))]
pA <- ggplot(pan, aes(panel, spearman_rho, fill = panel)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey50", linetype = "dashed") +
  geom_boxplot(outlier.alpha = 0, width = 0.55) +
  geom_jitter(width = 0.15, size = 0.9, alpha = 0.7) +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  labs(x = NULL, y = "Spearman rho vs macrophage pseudotime",
       title = "Macrophage marker panels along Kupffer->LAM pseudotime") +
  theme_masld()

# Panel B: top LAM-associated Hep->Mac
ax_valid <- ax[!is.na(receptor_rho) & receptor_padj < 0.05]
lam <- ax_valid[receptor_rho > 0][order(-score_diff)][1:15]
lam[, lab := paste0(ligand, " -> ", receptor, " (rec rho=", sprintf("%.2f", receptor_rho), ")")]
lam[, lab := factor(lab, levels = rev(unique(lab)))]
pB <- ggplot(lam, aes(score_diff, lab)) +
  geom_col(fill = "#C0392B", width = 0.7) +
  geom_text(aes(label = sprintf("%.2f", score_diff)), hjust = -0.1, size = 2) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(x = "LIANA score_diff (MASLD-enriched)", y = NULL,
       title = "LAM-associated Hep->Mac ligands (receptor_rho>0, padj<0.05)") +
  theme_masld() + theme(axis.text.y = element_text(size = 5.8))

# Panel C: top Kupffer-associated Hep->Mac
kup <- ax_valid[receptor_rho < 0][order(-score_diff)][1:15]
kup[, lab := paste0(ligand, " -> ", receptor, " (rec rho=", sprintf("%.2f", receptor_rho), ")")]
kup[, lab := factor(lab, levels = rev(unique(lab)))]
pC <- ggplot(kup, aes(score_diff, lab)) +
  geom_col(fill = "#2980B9", width = 0.7) +
  geom_text(aes(label = sprintf("%.2f", score_diff)), hjust = -0.1, size = 2) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(x = "LIANA score_diff", y = NULL,
       title = "Kupffer-associated Hep->Mac ligands (receptor_rho<0, padj<0.05)") +
  theme_masld() + theme(axis.text.y = element_text(size = 5.8))

fig <- pA / (pB + pC) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 8, face = "bold"))

ggsave(file.path(FIGS_CELLTYPE_DIR, "figS_D1_kupffer_lam_trajectory.pdf"),
       fig, width = 14, height = 12)
message("Saved D1 figure")
