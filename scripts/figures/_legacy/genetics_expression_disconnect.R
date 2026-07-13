#!/usr/bin/env Rscript
# genetics_expression_disconnect.R
# Fig 3E — the genetics<->expression cell-type disconnect (Fig 2 -> Fig 3 bridge).
# One point per COARSE liver LINEAGE (7; immune/stromal fine types grouped):
#   x = per-lineage GWAS-COLOC marker enrichment (genetic-causal signal)
#   y = per-lineage transcriptional disease effect (mean |log2FC|, scRNA pseudobulk
#       MASLD vs Healthy, lineage-aggregated)
# Headline: genetic-causal signal concentrates in immune/non-parenchymal lineages
# while the transcriptional disease response is dominated by hepatocytes — the two
# evidence layers diverge by cell type.
#
# LINEAGES (2026-07-02): Hepatocytes / Cholangiocytes / Endothelial /
#   Mesenchymal(fibro/stellate) / Myeloid(mac/mono/DC) / T_NK / B_Plasma.
#
# WHY scRNA pseudobulk, not bulk deconvolution: bulk deconvolution (TOAST) only
# yields a trustworthy cell-type DE for hepatocytes (92% of reads); for the rare
# non-parenchymal lineages the deconvolved betas are noise-inflated (Myeloid
# beta 3.5, 43% of genes "sig") — bulk cannot isolate non-hepatocyte expression.
# scRNA pseudobulk is the direct per-lineage measurement. The TOAST diagnostic is
# kept as a supp limitation (RNA-seq/results/celltype_attribution/coarse_lineage/).
#
# Output (PDF only): FIG2_DIR/panels/fig3e_genetics_expression_disconnect.pdf
# Sidecar: FIG2_DIR/panels/data/genetics_expression_disconnect.csv
#
# RIGOR: x traces to the canonical gene-level coloc lineage (marker-gene PP4
# Wilcoxon enrichment on each lineage's top-specificity genes — label
# "GWAS-COLOC marker enrichment", NOT LDSC heritability). y is a disease-vs-control
# contrast (MASLD_vs_Healthy), NOT scVI-inferred F-stage. Hepatocyte fold is
# version-unstable near 1.0; shown as the lowest/borderline lineage, never "depleted".

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(ggrepel)
})
set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

PANEL_DIR <- file.path(FIG2_DIR, "panels")
DATA_DIR  <- file.path(PANEL_DIR, "data")
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)
lab_size <- 6 / ggplot2::.pt

# ── X: per-lineage GWAS-COLOC marker enrichment (coarse, 201b) ────────────────
her <- fread(file.path(BASE,
  "RNA-seq/results/gwas_rna_integration/celltype_heritability_coarse_results.csv"),
  select = c("cell_type", "fold_enrichment_coloc", "fdr_coloc"))
setnames(her, c("fold_enrichment_coloc","fdr_coloc"), c("coloc_fold","coloc_fdr"))

# ── Y: per-lineage transcriptional disease effect (coarse scRNA pseudobulk) ───
sc_dir <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/pseudobulk_de_coarse")
yeff <- rbindlist(lapply(list.files(sc_dir, pattern = "_de\\.csv$", full.names = TRUE),
  function(f) {
    d <- fread(f); stopifnot(all(d$contrast == "MASLD_vs_Healthy"))
    data.table(cell_type = sub("_de\\.csv$", "", basename(f)),
               mean_abs_lfc = mean(abs(d$logFC), na.rm = TRUE), n_genes = nrow(d))
  }))

dt <- merge(her, yeff, by = "cell_type")
stopifnot(nrow(dt) == 7L)

# Display names + lineage colours (representative ct_palette hues)
disp <- c(Hepatocytes="Hepatocytes", Cholangiocytes="Cholangiocytes",
          Endothelial="Endothelial", Mesenchymal="Mesenchymal",
          Myeloid="Myeloid", T_NK="T / NK", B_Plasma="B / Plasma")
lin_col <- c("Hepatocytes"="#0D47A1", "Cholangiocytes"="#1565C0",
             "Endothelial"="#2E7D32", "Mesenchymal"="#F57F17",
             "Myeloid"="#C2185B", "T / NK"="#7B1FA2", "B / Plasma"="#5D4037")
dt[, lineage := disp[cell_type]]
dt[, coloc_sig := coloc_fdr < 0.05]

rho <- cor(dt$coloc_fold, dt$mean_abs_lfc, method = "spearman")
rho_lab <- sprintf("Spearman ρ = %.2f (n = 7 lineages)", rho)
cat(sprintf("[fig3E] coarse genetics<->expression Spearman rho = %.3f\n", rho))

# ── Plot ─────────────────────────────────────────────────────────────────────
xr <- range(dt$coloc_fold); yr <- range(dt$mean_abs_lfc)
p <- ggplot(dt, aes(x = coloc_fold, y = mean_abs_lfc)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.3, colour = "gray55") +
  geom_point(aes(colour = lineage, shape = coloc_sig), size = 2.6, stroke = 0.5, fill = "white") +
  ggrepel::geom_text_repel(aes(label = lineage, colour = lineage),
    size = lab_size, fontface = "plain", box.padding = 0.6, point.padding = 0.35,
    segment.size = 0.25, segment.color = "gray60", min.segment.length = 0,
    max.overlaps = Inf, seed = 42, bg.color = "white", bg.r = 0.12, show.legend = FALSE) +
  annotate("text", x = 1.0, y = mean(yr), hjust = -0.06, vjust = 1,
           label = "no genetic\nenrichment", size = lab_size, colour = "gray55", lineheight = 0.9) +
  annotate("text", x = xr[2], y = yr[2] * 1.03, hjust = 1, vjust = 1,
           label = rho_lab, size = lab_size, colour = "gray40") +
  scale_colour_manual(values = lin_col, guide = "none") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21), guide = "none") +
  scale_x_continuous(expand = expansion(mult = 0.12)) +
  scale_y_continuous(expand = expansion(mult = c(0.10, 0.12))) +
  labs(x = "GWAS-COLOC marker enrichment (fold)",
       y = "Transcriptional disease effect (mean |log₂FC|)") +
  theme_masld_compact()

message("CAPTION Fig 3E: Genetics<->expression disconnect across 7 liver lineages. ",
        "x = per-lineage GWAS-COLOC marker-gene enrichment (mean colocalization PP4 of a ",
        "lineage's top disease-specific genes vs all others; canonical gene-level coloc; ",
        "dashed line = no enrichment; open point = not FDR-significant [hepatocytes only]). ",
        "y = per-lineage transcriptional disease effect (mean |log2FC|, scRNA pseudobulk ",
        "MASLD vs Healthy, lineage-aggregated). Genetic-causal signal is enriched in ",
        "immune/non-parenchymal lineages (cholangiocyte, myeloid, T/NK) while the ",
        "transcriptional response is dominated by hepatocytes (highest |log2FC|, lowest/",
        "borderline genetic enrichment). Bulk deconvolution recovers cell-type DE only for ",
        "hepatocytes (see supp/limitation); scRNA pseudobulk is the per-lineage measure. ", rho_lab, ".")

out_pdf <- file.path(PANEL_DIR, "fig3e_genetics_expression_disconnect.pdf")
save_fig(p, out_pdf, width = fig_half_width + 0.3, height = 3.0)
cat("[fig3E] Saved:", out_pdf, "\n")

out_csv <- file.path(DATA_DIR, "genetics_expression_disconnect.csv")
fwrite(dt[order(-coloc_fold), .(lineage, cell_type,
        coloc_fold = round(coloc_fold,4), coloc_fdr = signif(coloc_fdr,3),
        coloc_sig, mean_abs_lfc = round(mean_abs_lfc,4), n_genes)], out_csv)
cat("[fig3E] Sidecar:", out_csv, "\n")
print(dt[order(-coloc_fold), .(lineage, coloc_fold = round(coloc_fold,3),
        coloc_sig, mean_abs_lfc = round(mean_abs_lfc,3))])
