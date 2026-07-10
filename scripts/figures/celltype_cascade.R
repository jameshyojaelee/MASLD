#!/usr/bin/env Rscript
# ============================================================================
# celltype_cascade.R
# Fig 2 panel c — Cell-type cascade per CRN transition
#
# Heatmap: 16 cell types (rows) x 4 CRN transitions (cols). Fill = log2FC of
# proportion at each transition; asterisks for padj<0.05; padj<0.01 = **. Hero
# annotations:
#   - T cells +1.34 log2FC at F0->F1 (adaptive-immune onset)
#   - Plasma cells +0.93 at F1->F2 (amplification)
#   - Cholangiocytes +3.31 at F3->F4 (ductular climax)
#
# Source: RNA-seq/results/granular_staging/two_transition_celltype_by_stage.csv
# (BayesPrism, all 4 adjacent CRN transitions).
# MuSiC-BayesPrism concord (Jaccard 0.774) cited as caption footer.
#
# Output: figures/main/fig3_RNAseq/panels/celltype_cascade.pdf
#   sized 180 x 50 mm (full width — 16 cell-types across 4 transitions).
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

PANEL_DIR <- file.path(FIG2_DIR, "panels")
DATA_DIR  <- file.path(PANEL_DIR, "data")
OUT_PDF   <- file.path(PANEL_DIR, "celltype_cascade.pdf")

ct <- fread(file.path(BASE,
  "RNA-seq/results/granular_staging/two_transition_celltype_by_stage.csv"))

# Normalize labels
ct[, celltype := gsub("\\.", " ", celltype)]
ct[, celltype := sub("Mono mono derived cells", "Mono", celltype)]
ct[, celltype := sub("Endothelial cells", "Endothelial", celltype)]
ct[, celltype := sub("Circulating NK NKT", "NK/NKT", celltype)]
ct[, celltype := sub("Resident NK", "NK resident", celltype)]
ct[, celltype := sub("Plasma cells", "Plasma", celltype)]
ct[, celltype := sub("B cells", "B", celltype)]
ct[, celltype := sub("T cells", "T", celltype)]

trans_levels <- c("F0→F1", "F1→F2", "F2→F3", "F3→F4")
ct[, transition_id := factor(sub("v", "→", transition), levels = trans_levels)]

# Order cell types by max |log2FC| across transitions (most dynamic at top)
ct_order <- ct[, .(max_abs = max(abs(log2FC_proportion), na.rm = TRUE)),
                by = celltype][order(-max_abs), celltype]
ct[, celltype := factor(celltype, levels = rev(ct_order))]

# Significance marker
ct[, sig_mark := fifelse(padj < 0.01, "**",
                  fifelse(padj < 0.05, "*", ""))]

# Hero numbers (preserve original values for annotation)
hero <- ct[(celltype == "T" & transition == "F0vF1") |
           (celltype == "Plasma" & transition == "F1vF2") |
           (celltype == "Cholangiocytes" & transition == "F3vF4")]
cat(sprintf("[hero] T@F0→F1 log2FC=%.2f; Plasma@F1→F2 log2FC=%.2f; Chol@F3→F4 log2FC=%.2f\n",
            ct[celltype == "T" & transition == "F0vF1", log2FC_proportion],
            ct[celltype == "Plasma" & transition == "F1vF2", log2FC_proportion],
            ct[celltype == "Cholangiocytes" & transition == "F3vF4", log2FC_proportion]))

p <- ggplot(ct, aes(x = transition_id, y = celltype,
                     fill = log2FC_proportion)) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_text(aes(label = sig_mark),
            size = GEOM_TEXT_6PT, color = "grey15", vjust = 0.7) +
  scale_fill_gradient2(name = "log2FC\nproportion",
                       low = "#1565C0", mid = "white",
                       high = masld_colors$nash,
                       midpoint = 0,
                       limits = c(-2, 3.5),
                       oob = scales::squish,
                       breaks = c(-2, 0, 2, 3)) +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_masld(base_size = 7) +
  theme(
    axis.text.x   = element_text(size = 6, face = "plain"),
    axis.text.y   = element_text(size = 6),
    legend.position = "right",
    legend.title  = element_text(size = 6),
    legend.text   = element_text(size = 6),
    legend.key.height = unit(0.4, "cm"),
    legend.key.width  = unit(0.18, "cm"),
    panel.grid    = element_blank()
  )

message("[caption] Cell-type cascade per CRN transition: T (F0->F1, +1.34) -> Plasma+Fib (F1->F2) -> Stellate+Chol (F2->F3) -> ductular (F3->F4, +3.31). BP x MuSiC Jaccard 0.774. * padj<0.05; ** padj<0.01.")

# ggsave(OUT_PDF, p,
#        width  = 180 / 25.4,
#        height = 50 / 25.4,
#        units  = "in",
#        device = cairo_pdf)
fwrite(ct, file.path(DATA_DIR, "celltype_cascade.csv"))
# cat(sprintf("[saved] %s\n", OUT_PDF))
