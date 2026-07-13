#!/usr/bin/env Rscript
# figS_celltype_attribution.R
# Figure for Analysis A1 — Per-cell-type intrinsic attribution of bulk DEGs.
# Spec: docs/superpowers/specs/2026-04-17-cell-type-resolved-masld-biology-design.md
#
# Panels:
#   A — Bar: counts of DEGs by primary-cell-type attribution class
#   B — Upset/heatmap: multi-cell-type DEG attribution matrix (top 40 genes)
#   C — Agreement with MuSiC hepatocyte class (confusion matrix)
#   D — Cell-type-specific DEG waterfall: bulk logFC split by primary celltype
#
# Output: figures/supplementary/figS_celltype_biology/figS_A1_celltype_attribution.pdf
# Env: rnaseq

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

source("scripts/figures/publication_theme.R")
source("scripts/figures/load_figure_data.R")

IN_FILE <- file.path(BASE, "RNA-seq/results/celltype_attribution",
                     "celltype_attribution_matrix.csv")
stopifnot(file.exists(IN_FILE))

dt <- fread(IN_FILE)
dt_deg <- dt[bulk_padj < 0.05 & abs(bulk_lfc) > 0.5]

# ---- Panel A: attribution class counts --------------------------------------
classA <- dt_deg[, .N, by = attribution_class][order(-N)][1:12]
classA[, attribution_class := factor(attribution_class, levels = attribution_class)]

pA <- ggplot(classA, aes(attribution_class, N)) +
  geom_col(fill = "#4472C4", color = "white", width = 0.65) +
  geom_text(aes(label = N), vjust = -0.3, size = GEOM_TEXT_6PT) +
  labs(x = NULL, y = "DEG count") +
  theme_masld() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 6))

# ---- Panel B: top multi-celltype DEGs heatmap -------------------------------
celltypes <- grep("^score_", names(dt), value = TRUE)
celltypes_clean <- sub("^score_", "", celltypes)

topN <- 40
# v2 schema: n_sig_concordant_ct + attribution_class == "multi_celltype"
if ("n_sig_concordant_ct" %in% names(dt_deg)) {
  top_genes <- dt_deg[attribution_class == "multi_celltype" |
                      grepl("_intrinsic_strong$", attribution_class)][
    order(-n_sig_concordant_ct, -abs(bulk_lfc))][1:topN]
} else {
  top_genes <- dt_deg[order(-abs(bulk_lfc))][1:topN]
}

if (nrow(top_genes) > 0) {
  mat_long <- melt(top_genes[, c("symbol", celltypes), with = FALSE],
                   id.vars = "symbol", variable.name = "celltype",
                   value.name = "score")
  mat_long[, celltype := sub("^score_", "", celltype)]
  mat_long[, score_clipped := pmax(pmin(score, 6), -6)]

  pB <- ggplot(mat_long, aes(celltype, symbol, fill = score_clipped)) +
    geom_tile(color = "white", linewidth = 0.2) +
    scale_fill_gradient2(low = "#2980B9", mid = "grey95", high = "#C0392B",
                         midpoint = 0, na.value = "grey85",
                         limits = c(-6, 6), name = "attribution\nscore (|t|)") +
    labs(x = NULL, y = NULL) +
    theme_masld() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 6),
          axis.text.y = element_text(size = 6))
} else {
  pB <- ggplot() + theme_masld()
}

# ---- Panel C: MuSiC agreement ------------------------------------------------
music_agree <- dt_deg[!is.na(music_category) & !is.na(primary_celltype),
                      .N, by = .(music_category, primary_celltype)]
if (nrow(music_agree) > 0) {
  pC <- ggplot(music_agree, aes(primary_celltype, music_category, fill = log10(N + 1))) +
    geom_tile(color = "white", linewidth = 0.2) +
    geom_text(aes(label = N), size = GEOM_TEXT_6PT) +
    scale_fill_gradient(low = "grey95", high = "#27AE60", name = "log10(N+1)") +
    labs(x = "Primary celltype (this analysis)",
         y = "MuSiC class (Script 25)") +
    theme_masld() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 6))
} else {
  pC <- ggplot() + theme_masld()
}

# ---- Panel D: waterfall per primary celltype --------------------------------
# Collapse *_strong + *_likely → CT label
dt_deg[, primary_celltype_label := sub("_intrinsic_(strong|likely)$", "", attribution_class)]
dt_deg[attribution_class %in% c("bulk_only","multi_celltype","NS_bulk"),
       primary_celltype_label := attribution_class]
top_ct_labels <- dt_deg[grepl("_intrinsic", attribution_class),
                        .N, by = primary_celltype_label][order(-N)][1:6, primary_celltype_label]
waterfall_dt <- dt_deg[primary_celltype_label %in% top_ct_labels]
waterfall_dt[, rank := rank(-bulk_lfc), by = primary_celltype_label]

pD <- ggplot(waterfall_dt, aes(rank, bulk_lfc, color = primary_celltype_label)) +
  geom_point(size = 0.35, alpha = 0.6) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey50") +
  facet_wrap(~ primary_celltype_label, scales = "free_x", nrow = 2) +
  labs(x = "Rank within cell type", y = "Bulk dream logFC") +
  theme_masld() +
  theme(legend.position = "none",
        strip.text = element_text(size = 6))

# ---- Assemble ---------------------------------------------------------------
message("[caption] A: Bulk DEGs by primary-cell-type attribution. B: Top multi-celltype DEGs (concordant in >=3 types). C: A1 attribution vs. existing MuSiC attribution. D: Bulk DEG fold-change landscape by primary cell type.")
fig <- (pA + pB) / (pC + pD) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 6, face = "plain"))

out_path <- file.path(FIGS_CELLTYPE_DIR, "figS_A1_celltype_attribution.pdf")
ggsave(out_path, fig, width = fig_full_width, height = 10 * (fig_full_width / 13))
message("Saved: ", out_path)
