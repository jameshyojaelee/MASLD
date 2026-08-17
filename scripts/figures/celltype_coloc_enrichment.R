#!/usr/bin/env Rscript
# =============================================================================
# Fig 4 (validation) panel — Cell-type COLOC-marker enrichment
# Output: figures/main/fig5_molecular_context/fig4g_celltype_coloc.pdf
# =============================================================================
# KEY MESSAGE:
#   GWAS-colocalized signal is IMMUNE / NON-PARENCHYMAL-leaning across liver
#   cell types, while transcriptomic dysregulation concentrates in hepatocytes
#   (Fig 4 elsewhere) — a genetics ⊥ expression DISCONNECT that makes
#   multi-evidence integration a necessity, not redundancy.
#
# METRIC (be precise): cell-type COLOC-MARKER enrichment. For each cell type we
#   take its marker-gene set and test whether their gene-level COLOC PP4 is
#   elevated vs all other genes (one-sided Wilcoxon). `fold_enrichment_coloc` =
#   mean COLOC of marker genes / mean COLOC of the rest. This is NOT partitioned
#   heritability (no LDSC), despite the source filename.
#
# Source (READ, never hardcoded):
#   RNA-seq/results/gwas_rna_integration/celltype_heritability_results.csv
#   columns: cell_type, fold_enrichment_coloc, wilcox_p_coloc, fdr_coloc, ...
#
# HONESTY (mirrors docs/technical/NUMBERS_HISTORY.md, corrected 2026-06-22):
#   - The ROBUST claim is immune / non-parenchymal enrichment (fold > 1).
#   - The HEPATOCYTE fold is VERSION-UNSTABLE (1.16× ↔ 0.93×, near 1.0) and is
#     NOT load-bearing. We show it as the lowest/borderline lineage but do NOT
#     label it "depleted." No on-panel "all FDR-significant" claim — the
#     FDR-corrected significance is a gradient (see fdr_coloc), so significance
#     is marked per cell type off fdr_coloc, not asserted in bulk.
#
# Env: rnaseq
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

# ── Data ─────────────────────────────────────────────────────────────────────
src <- file.path(BASE,
  "RNA-seq/results/gwas_rna_integration/celltype_heritability_results.csv")
ct <- fread(src)

# The file has a duplicated header row written as a data row; drop any row whose
# cell_type is literally "cell_type" and coerce the numeric columns.
ct <- ct[cell_type != "cell_type"]
num_cols <- c("fold_enrichment_coloc", "wilcox_p_coloc", "fdr_coloc")
ct[, (num_cols) := lapply(.SD, as.numeric), .SDcols = num_cols]
ct <- ct[!is.na(fold_enrichment_coloc)]
stopifnot(nrow(ct) == 11L)  # 11 true liver cell types

# ── Display names + compartment grouping ─────────────────────────────────────
disp_map <- c(
  Circulating_NK_NKT        = "Circulating NK / NKT",
  Cholangiocytes            = "Cholangiocytes",
  Macrophages               = "Macrophages",
  B_cells                   = "B cells",
  Resident_NK               = "Resident NK",
  T_cells                   = "T cells",
  `Mono+mono_derived_cells` = "Monocyte / mono-derived",
  Plasma_cells              = "Plasma cells",
  Endothelial_cells         = "Endothelial",
  Fibroblasts               = "Fibroblasts",
  Hepatocytes               = "Hepatocytes"
)
ct[, display := disp_map[cell_type]]
stopifnot(!any(is.na(ct$display)))

# Compartment: immune leukocytes vs non-parenchymal stromal/epithelial vs
# parenchymal hepatocyte. Hepatocyte gets its own neutral colour so it reads as
# the borderline/lowest lineage without being labelled "depleted".
immune <- c("Circulating_NK_NKT", "Resident_NK", "T_cells", "B_cells",
            "Plasma_cells", "Macrophages", "Mono+mono_derived_cells")
nonparenchymal <- c("Cholangiocytes", "Endothelial_cells", "Fibroblasts")
ct[, compartment := fifelse(cell_type %in% immune, "Immune",
                     fifelse(cell_type %in% nonparenchymal, "Non-parenchymal\n(stromal / epithelial)",
                             "Hepatocyte"))]

# ── FDR significance (honest, off the fdr_coloc column — NOT asserted in bulk) ─
# A cell type is marked significant only if fdr_coloc < 0.05. We add an asterisk
# to those bars; everything else (incl. hepatocytes) is unmarked.
ct[, fdr_sig := fdr_coloc < 0.05]

# ── Order: highest fold at top → hepatocytes at bottom ───────────────────────
setorder(ct, fold_enrichment_coloc)            # ascending -> bottom of plot = lowest after factor
ct[, display := factor(display, levels = display)]

# Value labels (e.g. "1.41×") just outside the bar end; asterisk if FDR-sig.
ct[, vlab := sprintf("%.2f×%s", fold_enrichment_coloc,
                     fifelse(fdr_sig, "*", ""))]

# ── Colours (compartment fills; ALL TEXT stays black) ────────────────────────
comp_cols <- c(
  "Immune"                                    = "#C9265E",  # Liang deep magenta
  "Non-parenchymal\n(stromal / epithelial)"   = "#F4A674",  # Liang warm peach
  "Hepatocyte"                                = "#9E9E9E"   # neutral gray (borderline lineage)
)
ct[, compartment := factor(compartment, levels = names(comp_cols))]

xr <- max(ct$fold_enrichment_coloc) * 1.28

# (Metric / honesty prose lives in the message() caption below, not on the panel.
#  The only on-panel annotation is the "* FDR < 0.05" note placed up by the legend.)

# ── Plot ─────────────────────────────────────────────────────────────────────
p <- ggplot(ct, aes(x = fold_enrichment_coloc, y = display, fill = compartment)) +
  geom_vline(xintercept = 1, linewidth = 0.3, color = "gray55") +
  geom_col(width = 0.68) +
  geom_text(aes(label = vlab), hjust = -0.12, size = GEOM_TEXT_6PT,
            color = "black") +
  scale_fill_manual(values = comp_cols, name = NULL) +
  scale_x_continuous(
    limits = c(0, xr),
    breaks = c(0, 0.5, 1.0, 1.5),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(x = "GWAS-COLOC marker enrichment (fold vs background)",
       y = NULL) +
  theme_masld_compact() +
  theme(
    axis.text.y    = element_text(color = "black"),
    legend.position = "top",
    legend.key.size = PUB_LEGEND_KEY,
    legend.text    = element_text(color = "black"),
    plot.margin    = margin(3, 8, 3, 4)
  )

# ── Save — RETIRED 2026-07-07 (do NOT re-create) ─────────────────────────────
# fig4g_celltype_coloc.pdf (cell-type COLOC-marker enrichment) was removed from
# Fig 4 as not relevant. The PDF is no longer written; the diagnostic messages
# below are kept for reference.
# out <- file.path(FIG4_DIR, "panels", "fig4g_celltype_coloc.pdf")
# dir.create(FIG4_DIR, recursive = TRUE, showWarnings = FALSE)
# pdf(out, width = 3.0, height = 2.3, useDingbats = FALSE)
# print(p)
# invisible(dev.off())

message("[caption] * = FDR < 0.05 (fdr_coloc)")
message("fig4g_celltype_coloc panel RETIRED — no PDF written.")
message("Source: ", src)
message(sprintf("  %d cell types | %d FDR-sig (fdr_coloc<0.05) | top = %s %.2fx | hep = %.2fx (fdr=%.2f)",
                nrow(ct), sum(ct$fdr_sig),
                ct[which.max(fold_enrichment_coloc), display],
                max(ct$fold_enrichment_coloc),
                ct[cell_type == "Hepatocytes", fold_enrichment_coloc],
                ct[cell_type == "Hepatocytes", fdr_coloc]))
# Per-cell-type folds plotted (for verification against NUMBERS.md §368):
for (i in order(-ct$fold_enrichment_coloc)) {
  message(sprintf("    %-26s %.3fx  wilcox_p=%.2e  fdr=%.3f%s",
                  ct$cell_type[i], ct$fold_enrichment_coloc[i],
                  ct$wilcox_p_coloc[i], ct$fdr_coloc[i],
                  ifelse(ct$fdr_sig[i], "  *", "")))
}
message("Metric: cell-type COLOC-MARKER enrichment (marker-gene COLOC PP4, 1-sided Wilcoxon vs other genes);")
message("  fold = mean marker PP4 / mean other PP4; vertical line = no enrichment (fold 1). NOT partitioned heritability.")
message("  Hepatocyte fold near 1.0 and version-unstable (1.16x vs 0.93x) — NOT load-bearing; robust claim = immune/")
message("  non-parenchymal enrichment, hepatocytes the lowest/borderline lineage.")
