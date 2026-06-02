#!/usr/bin/env Rscript
# figS_chromatin_ccc_coupling.R
# Figure for Analysis G2 — GWAS-ATAC chromatin × LIANA CCC coupling.
#
# Panels:
#   A — Enrichment forest: ATAC-regulated × LIANA ligands/receptors/DEGs
#   B — Motif-disrupted × LIANA ligands/receptors
#   C — Top 20 priority variant-disrupted CCC axes (MASLD-enriched)
#   D — Cell-type pair distribution among priority axes
#
# Output: figures/supplementary/figS_celltype_biology/figS_G2_chromatin_ccc_coupling.pdf

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

source("scripts/figures/publication_theme.R")
source("scripts/figures/load_figure_data.R")

IN <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/chromatin_ccc")
enr <- fread(file.path(IN, "gwas_atac_ccc_enrichment.csv"))
priority <- fread(file.path(IN, "priority_variant_disrupted_CCC_axes.csv"))

enr[, test_label := factor(test, levels = rev(test))]
enr[, log10OR := log2(pmax(odds_ratio, 1e-3))]

# ---- Panel A + B combined in forest ----------------------------------------
enr[, sig := fifelse(pvalue < 0.05, "sig", "ns")]
pA <- ggplot(enr, aes(log10OR, test_label, fill = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_col(width = 0.65, color = "white") +
  geom_text(aes(label = sprintf("OR=%.2f, p=%.3g\n(%d/%d)",
                                odds_ratio, pvalue, k, n_b)),
            hjust = ifelse(enr$log10OR > 0, -0.05, 1.05),
            size = 2.1, color = "black") +
  scale_fill_manual(values = c("sig" = "#C0392B","ns" = "grey65"),
                    guide = "none") +
  scale_x_continuous(expand = expansion(mult = 0.4)) +
  labs(x = "log2 OR (GWAS-ATAC vs CCC overlap)", y = NULL,
       title = "GWAS-ATAC × LIANA CCC enrichment") +
  theme_masld() +
  theme(axis.text.y = element_text(size = 6.5))

# ---- Panel C: top 20 priority axes -----------------------------------------
top <- priority[!is.na(source)][1:min(20, .N)]
top[, pair_label := paste0(source, " -> ", target, ": ",
                           ligand_complex, " -> ", receptor_complex)]
top[, pair_label := factor(pair_label, levels = rev(unique(pair_label)))]
top[, disruption_type := fifelse(receptor_motif_disrupted & ligand_motif_disrupted,
                                  "both disrupted",
                                  fifelse(receptor_motif_disrupted, "receptor disrupted",
                                          fifelse(ligand_motif_disrupted, "ligand disrupted",
                                                  "ATAC-linked")))]

pC <- ggplot(top, aes(score_diff, pair_label, fill = disruption_type)) +
  geom_col(width = 0.7, color = "white") +
  geom_text(aes(label = sprintf("%.3f", score_diff)),
            hjust = -0.1, size = 2, color = "black") +
  scale_fill_brewer(palette = "Set1", name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(x = "LIANA score_diff (MASLD - Control)", y = NULL,
       title = sprintf("Top 20 variant-disrupted CCC axes (of %d)", nrow(priority))) +
  theme_masld() +
  theme(axis.text.y = element_text(size = 5.5),
        legend.position = "bottom",
        legend.key.size = unit(3, "mm"),
        legend.text = element_text(size = 6))

# ---- Panel D: cell-type pair distribution ---------------------------------
ct_pairs <- priority[!is.na(source), .(n = .N),
                     by = .(source, target)][order(-n)][1:14]
ct_pairs[, pair := paste0(source, " -> ", target)]
ct_pairs[, pair := factor(pair, levels = rev(pair))]

pD <- ggplot(ct_pairs, aes(n, pair)) +
  geom_col(fill = "#4472C4", width = 0.65, color = "white") +
  geom_text(aes(label = n), hjust = -0.2, size = 2.4) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = "Variant-disrupted LR pairs", y = NULL,
       title = "Cell-type pairs with most variant-perturbed CCC") +
  theme_masld() +
  theme(axis.text.y = element_text(size = 6.5))

fig <- (pA / pC) | (pD / plot_spacer())
fig <- fig + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 8, face = "bold"))

out_path <- file.path(FIGS_CELLTYPE_DIR, "figS_G2_chromatin_ccc_coupling.pdf")
ggsave(out_path, fig, width = 14, height = 11)
message("Saved: ", out_path)
