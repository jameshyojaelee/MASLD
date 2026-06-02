#!/usr/bin/env Rscript
# figS_cellstate_signatures.R
# Figure for Analyses H2 + J1 + J2 + K1 — Cell-state signature scoring in bulk.
#
# Panels:
#   A — fgsea NES forest (all signatures, ordered by |NES|)
#   B — Top per-gene bulk LFC for top senescence/SASP signature
#   C — Top per-gene bulk LFC for ductular reaction (J2)
#   D — Top per-gene bulk LFC for T-cell exhaustion (J1)
#
# Output: figures/supplementary/figS_celltype_biology/figS_H2_J1_J2_K1_cellstate_signatures.pdf

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})

source("scripts/figures/publication_theme.R")
source("scripts/figures/load_figure_data.R")

IN <- file.path(BASE, "RNA-seq/results/celltype_attribution")
fr <- fread(file.path(IN, "cellstate_signatures_fgsea.csv"))
gene_lvl <- fread(file.path(IN, "cellstate_signatures_gene_level.csv"))

# Clean names
fr[, nice_name := gsub("_", " ", pathway)]
fr[, nice_name := factor(nice_name, levels = rev(nice_name[order(NES)]))]
fr[, sig := fifelse(padj < 0.05, "padj<0.05", "ns")]

pA <- ggplot(fr, aes(NES, nice_name, fill = sig)) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey50") +
  geom_col(width = 0.7, color = "white") +
  geom_text(aes(label = sprintf("NES=%.2f, p=%.2g\n(size=%d)", NES, pval, size),
                x = NES + sign(NES) * 0.03),
            hjust = ifelse(fr$NES > 0, 0, 1), size = 2.1) +
  scale_fill_manual(values = c("padj<0.05" = "#C0392B", "ns" = "grey65"), name = NULL) +
  scale_x_continuous(expand = expansion(mult = 0.4)) +
  labs(x = "NES (bulk dream t-stat ranking)", y = NULL,
       title = "Cell-state signatures enriched in MASLD bulk DE") +
  theme_masld() +
  theme(axis.text.y = element_text(size = 6))

make_gene_panel <- function(sig_name, title_text) {
  g <- gene_lvl[signature == sig_name][order(-bulk_t)]
  if (nrow(g) == 0) return(ggplot() + labs(title = paste("No data:", sig_name)) + theme_masld())
  top_up <- g[1:min(12, .N)]
  top_dn <- g[(max(1, .N - 11)):.N]
  combo <- unique(rbind(top_up, top_dn))
  combo <- combo[order(bulk_t)]
  combo[, symbol := factor(symbol, levels = unique(symbol))]
  combo[, direction := fifelse(bulk_t > 0, "up", "down")]

  ggplot(combo, aes(bulk_t, symbol, fill = direction)) +
    geom_col(width = 0.7, color = "white") +
    geom_text(aes(label = sprintf("%.1f", bulk_t)),
              hjust = ifelse(combo$bulk_t > 0, -0.1, 1.1), size = 2) +
    scale_fill_manual(values = c("up" = "#C0392B","down" = "#2980B9"), guide = "none") +
    scale_x_continuous(expand = expansion(mult = 0.2)) +
    labs(x = "Bulk dream t-stat", y = NULL, title = title_text) +
    theme_masld() +
    theme(axis.text.y = element_text(size = 6))
}

pB <- make_gene_panel("H2_SASP_profibrotic", "H2 SASP pro-fibrotic (bulk DE)")
pC <- make_gene_panel("J2_ductular_reaction", "J2 Ductular reaction (bulk DE)")
pD <- make_gene_panel("J1_Tcell_exhaustion", "J1 T-cell exhaustion (bulk DE)")

fig <- (pA / (pB + pC + pD)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 8, face = "bold"))

out_path <- file.path(FIGS_CELLTYPE_DIR, "figS_H2_J1_J2_K1_cellstate_signatures.pdf")
ggsave(out_path, fig, width = 14, height = 13)
message("Saved: ", out_path)
