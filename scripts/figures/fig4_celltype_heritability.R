#!/usr/bin/env Rscript
# fig4_celltype_heritability.R — Cell-type GWAS heritability enrichment
#
# Panel for Fig 4: which cell types' specifically expressed genes are
# enriched for MASLD GWAS signal? Two sub-panels:
# (a) Bar chart of enrichment per cell type (Wilcoxon + Fisher)
# (b) DEG quintile × COLOC trend
#
# Output: figures/main/fig4_causal_architecture/panel_celltype_heritability.pdf

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

# RETIRED 2026-07-07: all three outputs of this script (panel_celltype_heritability.pdf,
# panel_celltype_enrichment_bar.pdf, panel_geneset_enrichment_forest.pdf) are cell-type
# heritability panels that are NOT part of the Fig 2 (genetics) or FigS2 set — they were
# stale leftovers in fig2_genetics/panels/ (the FIG3_DIR constant resolves there). Retired
# wholesale (early quit) so they are never regenerated. Plotting code kept for provenance.
message("[fig4_celltype_heritability] RETIRED 2026-07-07 — outputs not in Fig2/FigS2 set; no panels written.")
quit(save = "no", status = 0)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))

source(file.path(BASE, "scripts/figures/load_figure_data.R"))
dir.create(file.path(FIG3_DIR, "panels"), showWarnings = FALSE, recursive = TRUE)
outdir <- file.path(FIG3_DIR, "panels")

# ===========================================================================
# 1. Load cell-type heritability results
# ===========================================================================
cat("Loading cell-type heritability results...\n")

herit <- fread(file.path(BASE,
  "RNA-seq/results/gwas_rna_integration/celltype_heritability_results.csv"))
cat("  Cell types:", nrow(herit), "\n")

# Clean cell type names (replace underscore with space)
herit[, cell_type_label := gsub("_", " ", cell_type)]

# Sort by enrichment
setorder(herit, -fold_enrichment_coloc)

# Significance stars
herit[, sig_label := fcase(
  fdr_coloc < 0.001, "***",
  fdr_coloc < 0.01, "**",
  fdr_coloc < 0.05, "*",
  default = "ns"
)]

# ===========================================================================
# 2. Panel (a): Cell-type enrichment bar chart
# ===========================================================================
cat("Generating enrichment bar chart...\n")

# Color by significance
herit[, bar_color := fifelse(fdr_coloc < 0.05, "Significant", "Not significant")]

# Map cell types to palette
herit[, ct_color := ct_palette[gsub("_", " ", cell_type)]]
herit[is.na(ct_color), ct_color := masld_colors$ns]

p_a <- ggplot(herit, aes(x = reorder(cell_type_label, fold_enrichment_coloc),
                          y = fold_enrichment_coloc)) +
  geom_col(aes(fill = bar_color), width = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  # Significance stars
  geom_text(aes(label = sig_label, y = fold_enrichment_coloc + 0.02),
            hjust = 0, size = GEOM_TEXT_6PT) +
  # Count annotation
  geom_text(aes(label = paste0("n=", n_coloc_in_set), y = 0.02),
            hjust = 0, size = GEOM_TEXT_6PT, color = "white") +
  scale_fill_manual(values = c("Significant" = masld_colors$up,
                                "Not significant" = masld_colors$ns),
                    name = NULL) +
  coord_flip(ylim = c(0, max(herit$fold_enrichment_coloc) * 1.15)) +
  labs(
    x = NULL,
    y = "Fold enrichment\n(COLOC PP.H4 in cell-type genes vs background)"
  ) +
  theme_masld(base_size = 7) +
  theme(legend.position = c(0.95, 0.05),
        legend.justification = c(1, 0),
        legend.background = element_rect(fill = "white", color = NA))
message("[caption] a: Cell-type GWAS enrichment — fold enrichment of COLOC PP.H4 in cell-type-specific genes vs background")

# ===========================================================================
# 3. Panel (b): DEG quintile × COLOC enrichment
# ===========================================================================
cat("Generating DEG quintile trend...\n")

quintile_file <- file.path(BASE,
  "RNA-seq/results/gwas_rna_integration/deg_quintile_coloc_enrichment.csv")
if (file.exists(quintile_file)) {
  quintile <- fread(quintile_file)

  p_b <- ggplot(quintile, aes(x = t_quintile, y = prop_coloc_05)) +
    geom_col(fill = masld_colors$gwas, width = 0.6) +
    geom_text(aes(label = sprintf("%.1f%%", prop_coloc_05 * 100)),
              vjust = -0.5, size = GEOM_TEXT_6PT) +
    scale_y_continuous(labels = scales::percent,
                       expand = expansion(mult = c(0, 0.15))) +
    labs(
      x = "DEG significance quintile (Q1=weakest, Q5=strongest)",
      y = "% genes with\nCOLOC PP.H4 > 0.5"
    ) +
    theme_masld(base_size = 7)
  message("[caption] b: GWAS signal by DEG strength — % genes with COLOC PP.H4 > 0.5 across DEG significance quintiles")
} else {
  p_b <- ggplot() + theme_void()
  message("[caption] b: DEG quintile data not available")
}

# ===========================================================================
# 4. Panel (c): Gene-set enrichment forest plot (from 202)
# ===========================================================================
cat("Generating gene-set enrichment forest...\n")

geneset_file <- file.path(BASE,
  "RNA-seq/results/gwas_rna_integration/geneset_enrichment_results.csv")
if (file.exists(geneset_file)) {
  geneset <- fread(geneset_file)
  geneset <- geneset[!is.na(fold_enrichment) & gene_set != ""]

  # Clean names
  geneset[, gene_set_label := gsub("_", " ", gene_set)]
  geneset[, gene_set_label := gsub("^sc ", "scRNA: ", gene_set_label)]

  # Sort
  setorder(geneset, fold_enrichment)

  # Significance
  geneset[, sig := fdr_coloc < 0.05]

  p_c <- ggplot(geneset, aes(x = fold_enrichment,
                              y = reorder(gene_set_label, fold_enrichment))) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.3) +
    geom_point(aes(color = sig, size = n_in_set)) +
    scale_color_manual(values = c("TRUE" = masld_colors$up,
                                  "FALSE" = masld_colors$ns),
                       labels = c("TRUE" = "FDR < 0.05", "FALSE" = "ns"),
                       name = NULL) +
    scale_size_continuous(range = c(1, 3), name = "Gene set\nsize") +
    labs(
      x = "Fold enrichment (COLOC in gene set vs background)",
      y = NULL
    ) +
    theme_masld(base_size = 7) +
    theme(legend.position = "right")
  message("[caption] c: Gene-set GWAS enrichment — fold enrichment of COLOC PP.H4 in curated gene sets vs background")
} else {
  p_c <- ggplot() + theme_void()
}

# ===========================================================================
# 5. Combine and save
# ===========================================================================
cat("Assembling composite figure...\n")

p_combined <- (p_a | p_b) / p_c +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = list(c("a", "b", "c")))

ggsave(file.path(outdir, "panel_celltype_heritability.pdf"),
       p_combined, width = 170, height = 140, units = "mm", device = cairo_pdf)
cat("Saved: panel_celltype_heritability.pdf\n")

# Individual panels
ggsave(file.path(outdir, "panel_celltype_enrichment_bar.pdf"),
       p_a, width = 85, height = 65, units = "mm", device = cairo_pdf)
ggsave(file.path(outdir, "panel_geneset_enrichment_forest.pdf"),
       p_c, width = 85, height = 80, units = "mm", device = cairo_pdf)

cat("Done.\n")
