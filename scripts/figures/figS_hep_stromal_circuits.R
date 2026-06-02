#!/usr/bin/env Rscript
# figS_hep_stromal_circuits.R
# Figure for Analyses D2 + D3 — Hep-stromal circuits + stromal-stromal CCC.
# Spec: docs/superpowers/specs/2026-04-17-cell-type-resolved-masld-biology-design.md
#
# Panels:
#   A — Hep<->Fib top LR pairs (bidirectional, bulk-concordant highlighted)
#   B — Hep<->Mac top LR pairs
#   C — Hep<->LSEC top LR pairs
#   D — Stromal-stromal top LR pairs (Fib/LSEC/Mac inter-crosstalk)
#
# Output: figures/supplementary/figS_celltype_biology/figS_D2_D3_hep_stromal_circuits.pdf

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

source("scripts/figures/publication_theme.R")
source("scripts/figures/load_figure_data.R")

CCC <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/ccc")
hep_fib  <- fread(file.path(CCC, "D2_hep_HSC_bidirectional.csv"))
hep_mac  <- fread(file.path(CCC, "D2_hep_Mac_bidirectional.csv"))
hep_lsec <- fread(file.path(CCC, "D2_hep_LSEC_bidirectional.csv"))
stromal  <- fread(file.path(CCC, "D3_stromal_stromal.csv"))

make_panel <- function(dt, title, top_n = 14) {
  if (nrow(dt) == 0) return(ggplot() + labs(title = paste("No data:", title)) + theme_masld())
  dt_top <- dt[order(-score_diff)][1:min(top_n, nrow(dt))]
  # Include direction in label to disambiguate shared L-R
  dt_top[, pair_label := paste0(direction, ": ", ligand, " -> ", receptor)]
  dt_top[, pair_label := factor(pair_label, levels = rev(unique(pair_label)))]
  dt_top[, direction := as.character(direction)]
  dt_top[, bulk_ok := fifelse(!is.na(both_concordant) & both_concordant,
                               "bulk-concordant", "not concordant")]

  ggplot(dt_top, aes(score_diff, pair_label, fill = direction, color = bulk_ok)) +
    geom_col(width = 0.7, linewidth = 0.4) +
    geom_text(aes(label = sprintf("%.2f", score_diff)),
              hjust = -0.1, size = 2.1, color = "black") +
    scale_fill_brewer(palette = "Set2", name = NULL) +
    scale_color_manual(values = c("bulk-concordant" = "black",
                                  "not concordant"   = "grey75"),
                       name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(x = "LIANA score_diff (MASLD - Control)", y = NULL, title = title) +
    theme_masld() +
    theme(axis.text.y = element_text(size = 6),
          legend.position = "bottom",
          legend.key.size = unit(3, "mm"),
          legend.text = element_text(size = 6))
}

pA <- make_panel(hep_fib,  "Hepatocyte <-> Fibroblast (HSC proxy): top LR pairs")
pB <- make_panel(hep_mac,  "Hepatocyte <-> Macrophage: top LR pairs")
pC <- make_panel(hep_lsec, "Hepatocyte <-> LSEC: top LR pairs")

# Panel D: stromal-stromal top pairs
stromal_top <- stromal[order(-score_diff)][1:14]
stromal_top[, axis := paste0(source, " -> ", target)]
stromal_top[, pair_label := paste0(axis, ": ", ligand, " -> ", receptor)]
stromal_top[, pair_label := factor(pair_label, levels = rev(unique(pair_label)))]
stromal_top[, bulk_ok := fifelse(!is.na(both_concordant) & both_concordant,
                                 "bulk-concordant", "not concordant")]

pD <- ggplot(stromal_top, aes(score_diff, pair_label, fill = axis, color = bulk_ok)) +
  geom_col(width = 0.7, linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", score_diff)),
            hjust = -0.1, size = 2.1, color = "black") +
  scale_fill_brewer(palette = "Dark2", name = NULL) +
  scale_color_manual(values = c("bulk-concordant" = "black",
                                "not concordant"   = "grey75"),
                     name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "LIANA score_diff", y = NULL,
       title = "Stromal-stromal CCC (Fib/LSEC/Mac)") +
  theme_masld() +
  theme(axis.text.y = element_text(size = 6),
        legend.position = "bottom",
        legend.key.size = unit(3, "mm"),
        legend.text = element_text(size = 5.5))

fig <- (pA + pB) / (pC + pD) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 8, face = "bold"))

out_path <- file.path(FIGS_CELLTYPE_DIR, "figS_D2_D3_hep_stromal_circuits.pdf")
ggsave(out_path, fig, width = 14, height = 11)
message("Saved: ", out_path)
