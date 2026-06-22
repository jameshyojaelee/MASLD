#!/usr/bin/env Rscript
# fig3f — Finemapping (SuSiE-X PIP) vs. RNA-seq DEG.
# Mirrors fig3e structure (side-by-side disease | severity) but the y-axis is
# the per-gene SuSiE-X max PIP across the 28-GWAS portfolio (1KG EUR LD panel).
#
# PIP source: GWAS/finemapping/results/susiex_1kg/susiex_gene_summary_1kg.csv
# (the atlas does not currently carry susiex_max_pip).
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(ggrepel); library(patchwork)
})

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
PANEL_DIR <- file.path(BASE, "figures/main/fig2_genetics/panels")
fig_half_width <- 4.5

source(file.path(BASE, "scripts/figures/publication_theme.R"))

# ---- Anchor sets (shared with fig3e) ----
drug_target_genes <- c("THRB","RORA","NR1H4","PPARA","PPARG","GLP1R")
govaere_25 <- c("AKR1B10","DUSP6","GDF15","THBS2","A2M","CDH2",
                "COL1A1","COL3A1","COL4A1","COL4A2","COL6A3","DCN",
                "FBN1","FSTL1","IGFBP7","LUM","MFAP4","MMP2",
                "POSTN","SPARC","SPP1","TAGLN","THY1","TIMP1","VCAN")
govaere_fibrosis <- c("ACTA2","COL1A1","COL1A2","COL3A1","COL4A1","COL4A2",
                      "CCN2","IL8","DCN","FN1","IL1B","IL6","LOX","LOXL2",
                      "MMP2","MMP9","PDGFRB","SERPINE1","TGFB1","TGFB2",
                      "TIMP1","TIMP2","TNF","VIM","CCN4")
steatosite_15 <- c("MT1F","KCNH7","COL25A1","RASD2","CCN2","STC1",
                   "GDNF","PRRX1","FGF7","LCNL1","DPEP1","CHRDL2",
                   "LHX6","POU4F1","CDH16")
piras_key <- c("HSD17B13","PNPLA3","TM6SF2","MARC1","MBOAT7",
               "CIDEB","GPAM","APOH","SERPINA6","GHR",
               "CYP7A1","ABCB4","SLC22A1","SLC27A5","ACSM3",
               "ADH4","ADH1B","CYP2E1","CYP4A11","CES1",
               "AKR1D1","SLC10A1","HGF","CRP","LCAT",
               "CETP","ALDOB","PCK1","GLS2","HAO1")
feng_top <- c("EFHD1","MLIP","TREM2","SPP1","GPNMB","CCL2",
              "CCL20","CXCL1","CXCL6","IL1B","IL32",
              "COL1A1","COL1A2","COL3A1","FN1","LOX","LOXL2",
              "ACTA2","PDGFRB","TGFB1","SERPINE1","MMP9",
              "AKR1B10","GDF15","THY1","THBS2","LUM")
landmark_masld_genes <- c("LPL","FABP4","MMP9","CFLAR")
masld_associated_genes <- unique(c(govaere_25, govaere_fibrosis,
                                   steatosite_15, piras_key, feng_top,
                                   landmark_masld_genes))

classify_gene <- function(g) {
  fcase(
    g %in% drug_target_genes,      "Drug target",
    g %in% masld_associated_genes, "MASLD-associated",
    default                        = "Novel")
}

tier_colors_3e <- c(
  "Drug target"        = "#880E4F",
  "Novel"              = "#00695C",
  "MASLD-associated"   = "#7B1FA2")

# ---- Atlas + PIP merge ----
atlas_lite <- fread(
  file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"),
  select = c("human_symbol", "bulk_logFC", "bulk_padj",
             "nafl_vs_nash_logFC", "nafl_vs_nash_padj",
             "coloc_susie_best_pp4"))
stopifnot(all(c("bulk_padj", "bulk_logFC") %in% names(atlas_lite)))
setnames(atlas_lite, "human_symbol", "gene")

susiex_pip <- fread(
  file.path(BASE, "GWAS/finemapping/results/susiex_1kg/susiex_gene_summary_1kg.csv"),
  select = c("GeneSymbol", "susiex_max_pip"))
setnames(susiex_pip, "GeneSymbol", "gene")
atlas_lite <- merge(atlas_lite, susiex_pip, by = "gene", all.x = TRUE)

# Genes supported by BOTH methods — shared label set, identical to fig3e
dual_anchor_genes <- atlas_lite[
  !is.na(coloc_susie_best_pp4) & coloc_susie_best_pp4 >= 0.5 &
  !is.na(susiex_max_pip)       & susiex_max_pip >= 0.5, gene]

# ---- Panel builder ----
build_panel <- function(x_col, padj_col, x_label, title, show_y_label = TRUE) {
  dat <- copy(atlas_lite)
  setnames(dat, c(x_col, padj_col), c("x_logFC", "x_padj"))

  # Label anchors: dual-evidence genes that are also DEGs in this contrast
  atlas_anchor_genes <- dat[
    gene %in% dual_anchor_genes &
    !is.na(x_padj) & x_padj < 0.05 &
    !is.na(x_logFC) & abs(x_logFC) > 0.3, gene]
  all_anchors <- unique(c(drug_target_genes, landmark_masld_genes,
                          atlas_anchor_genes))

  gene_tbl <- dat[gene %in% all_anchors]
  gene_tbl[, pip := susiex_max_pip]
  gene_tbl[is.na(pip), pip := 0]
  gene_tbl[is.na(x_logFC), x_logFC := 0]

  bg <- dat[!is.na(x_logFC) & !is.na(x_padj) & x_padj < 0.05]
  bg[, pip := susiex_max_pip]
  bg[is.na(pip), pip := 0]
  bg <- bg[!gene %in% gene_tbl$gene]

  gene_tbl[, fig5_group := classify_gene(gene)]

  gene_tbl[, fig5_group := factor(fig5_group, levels = c(
    "Drug target", "Novel", "MASLD-associated"))]
  gene_tbl[, do_label := TRUE]
  gene_tbl[, is_sig := factor(ifelse(!is.na(x_padj) & x_padj < 0.05,
                                      "sig", "ns"), levels = c("sig", "ns"))]

  gene_tbl[, label_text := gene]

  p <- ggplot(gene_tbl, aes(x = x_logFC, y = pip, color = fig5_group)) +
    # shaded band highlighting high-confidence finemapping zone
    annotate("rect",
             xmin = -Inf, xmax = Inf, ymin = 0.9, ymax = 1.02,
             fill = "#E8F5E9", alpha = 0.45, color = NA) +
    geom_point(data = bg, inherit.aes = FALSE,
               aes(x = x_logFC, y = pip),
               color = "gray80", size = 0.4, alpha = 0.4, shape = 16) +
    geom_hline(yintercept = c(0.5, 0.9), linetype = "dashed",
               linewidth = 0.25, color = "gray55") +
    geom_vline(xintercept = 0, linewidth = 0.2, color = "gray70") +
    geom_point(aes(shape = is_sig), size = 2.2, alpha = 0.9) +
    ggrepel::geom_text_repel(data = gene_tbl[do_label == TRUE],
                     aes(label = label_text),
                     size = 2.1, max.overlaps = Inf,
                     force = 3, force_pull = 1,
                     box.padding = 0.3, point.padding = 0.15,
                     segment.size = 0.2, segment.color = "gray55",
                     min.segment.length = 0, fontface = "italic",
                     bg.color = "white", bg.r = 0.12,
                     show.legend = FALSE,
                     max.iter = 8000, seed = 1) +
    scale_color_manual(values = tier_colors_3e, name = NULL, drop = FALSE) +
    scale_x_continuous(expand = expansion(mult = c(0.08, 0.08))) +
    scale_y_continuous(limits = c(-0.02, 1.05),
                       breaks = c(0, 0.5, 0.9, 1.0),
                       labels = c("0", "0.5", "0.9", "1")) +
    labs(x = x_label,
         y = if (show_y_label) "Best SuSiE-X finemapping PIP" else NULL,
         title = title) +
    theme_masld() +
    scale_shape_manual(values = c("sig" = 16, "ns" = 1),
                       name = "DEG Significance",
                       labels = c("sig" = "FDR < 0.05", "ns" = "n.s."))
  if (!show_y_label) {
    p <- p + theme(axis.title.y = element_blank())
  }
  p
}

p_left <- build_panel(
  x_col       = "bulk_logFC",
  padj_col    = "bulk_padj",
  x_label     = expression("Transcript log"[2]*"FC (MASLD vs control)"),
  title       = "Finemapping vs. RNA-seq DEG",
  show_y_label = TRUE)

p_right <- build_panel(
  x_col       = "nafl_vs_nash_logFC",
  padj_col    = "nafl_vs_nash_padj",
  x_label     = expression("Transcript log"[2]*"FC (MASH vs MASL)"),
  title       = "Finemapping vs. MASH-vs-MASL DEG",
  show_y_label = FALSE)

combined <- (p_left | p_right) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom",
        legend.box = "vertical",
        legend.spacing.y = unit(0.05, "cm"),
        legend.key.size = unit(0.25, "cm")) &
  guides(color = guide_legend(override.aes = list(size = 2.5), nrow = 1),
         shape = guide_legend(nrow = 1))

ggsave(file.path(PANEL_DIR, "susiex_cross_ancestry.pdf"), combined,
       width = (fig_half_width + 0.5) * 2, height = 5.4,
       device = cairo_pdf)
cat("[render_fig3f] wrote", file.path(PANEL_DIR, "susiex_cross_ancestry.pdf"), "\n")
cat("[render_fig3f] DONE\n")
