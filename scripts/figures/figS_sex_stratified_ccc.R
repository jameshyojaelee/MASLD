#!/usr/bin/env Rscript
# figS_sex_stratified_ccc.R
# Figure for Analysis L1 — Sex-stratified CCC.
#
# Panels:
#   A — Enrichment forest: sex-biased genes × CCC ligands/receptors
#   B — Counts of MASLD-enriched CCC axes by sex-bias class
#   C — Top 15 female-ligand MASLD-enriched CCC axes
#   D — Top 15 female-receptor MASLD-enriched CCC axes
#
# Output: figures/supplementary/figS_celltype_biology/figS_L1_sex_stratified_ccc.pdf

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})

source("scripts/figures/publication_theme.R")
source("scripts/figures/load_figure_data.R")

IN <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/sex_ccc")
enr <- fread(file.path(IN, "sex_ccc_enrichment.csv"))
ax  <- fread(file.path(IN, "sex_specific_ccc_axes.csv"))

# Panel A
enr[, sig := fifelse(pvalue < 0.05, "sig", "ns")]
enr[, test := factor(test, levels = rev(test))]
pA <- ggplot(enr, aes(odds_ratio, test, fill = sig)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_col(width = 0.65, color = "white") +
  geom_text(aes(label = sprintf("OR=%.2f, p=%.2g\n(%d/%d)", odds_ratio, pvalue, k, n_b)),
            hjust = -0.05, size = 2.2) +
  scale_fill_manual(values = c("sig" = "#C0392B","ns" = "grey65"), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.4))) +
  labs(x = "Odds ratio", y = NULL,
       title = "Sex-biased DEGs enriched in CCC ligand/receptor sets") +
  theme_masld() +
  theme(axis.text.y = element_text(size = 6))

# Panel B
cnt <- ax[axis_sex_bias != "Not_sex_biased",
          .N, by = axis_sex_bias][order(-N)]
cnt[, axis_sex_bias := factor(axis_sex_bias, levels = rev(axis_sex_bias))]
pB <- ggplot(cnt, aes(N, axis_sex_bias, fill = axis_sex_bias)) +
  geom_col(width = 0.6, color = "white") +
  geom_text(aes(label = N), hjust = -0.2, size = 2.6) +
  scale_fill_manual(values = c("Female_ligand" = "#E91E63","Female_receptor" = "#F06292",
                               "Male_ligand" = "#1565C0","Male_receptor" = "#42A5F5"),
                    guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = "MASLD-enriched CCC axes", y = NULL,
       title = "Sex-bias of MASLD-enriched CCC axes") +
  theme_masld()

make_top <- function(dt, bias_label, title_text, n = 15) {
  sub <- dt[axis_sex_bias == bias_label][order(-score_diff)][1:min(n, .N)]
  sub[, label := paste0(source, " -> ", target, ": ",
                        ligand_complex, " -> ", receptor_complex)]
  sub[, label := factor(label, levels = rev(unique(label)))]
  ggplot(sub, aes(score_diff, label)) +
    geom_col(fill = "#E91E63", width = 0.7, color = "white") +
    geom_text(aes(label = sprintf("%.3f", score_diff)),
              hjust = -0.1, size = 2) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
    labs(x = "LIANA score_diff", y = NULL, title = title_text) +
    theme_masld() +
    theme(axis.text.y = element_text(size = 5.2))
}
pC <- make_top(ax, "Female_ligand",   "Top female-ligand MASLD-enriched axes")
pD <- make_top(ax, "Female_receptor", "Top female-receptor MASLD-enriched axes")

fig <- (pA + pB) / (pC + pD) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 8, face = "plain"))

out_path <- file.path(FIGS_CELLTYPE_DIR, "figS_L1_sex_stratified_ccc.pdf")
ggsave(out_path, fig, width = 14, height = 11)
message("Saved: ", out_path)
