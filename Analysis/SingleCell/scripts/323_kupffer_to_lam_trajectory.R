#!/usr/bin/env Rscript
# 323_kupffer_to_lam_trajectory.R
#
# Analysis D1 (v1) — Kupffer->LAM trajectory with hepatocyte-derived drivers.
#
# Spec: docs/superpowers/specs/2026-04-17-cell-type-resolved-masld-biology-design.md
#
# Strategy (v1, uses existing pseudotime data):
#   1. Load macrophage pseudotime correlation table (pseudotime_corr_Macrophages.csv)
#      — genes with Spearman rho vs macrophage pseudotime.
#   2. Score Kupffer markers (CLEC4F, TIMD4, MARCO, VSIG4, CD5L, CD163L1) vs
#      LAM markers (TREM2, CD9, GPNMB, SPP1, LIPA, LGALS3, FABP5, CTSB/D/L).
#   3. Expect Kupffer markers negatively correlated, LAM positively correlated
#      with pseudotime (Kupffer -> LAM direction).
#   4. Overlay LIANA Hep->Mac ligands + their bulk DE, ranking candidate
#      hepatocyte-derived drivers of LAM transition by their pseudotime correlation.
#
# v2 (needs NicheNet installed): NicheNet upstream analysis on LAM signature.
#
# Env: rnaseq
# Outputs: Analysis/SingleCell/results_gpu_v2/macrophage_trajectory/

suppressPackageStartupMessages({
  library(data.table)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
PT   <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/pseudotime")
CCC  <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/ccc")
OUTDIR <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/macrophage_trajectory")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

message("[1] Loading macrophage pseudotime correlations...")
pt_corr <- fread(file.path(PT, "pseudotime_corr_Macrophages.csv"))
message(sprintf("  Genes: %d", nrow(pt_corr)))

# ---------------------------------------------------------------------------
# DIRECTION SANITY CHECK
# ---------------------------------------------------------------------------
# Canonical biology (Remmerie 2020, Tran 2020, Jaitin 2019): resident Kupffer
# cells DEPLETE in MASLD and monocyte-derived LAMs EXPAND — so along the
# Kupffer->LAM pseudotime, canonical Kupffer markers (TIMD4, MARCO, VCAM1)
# must NEGATIVELY correlate with pseudotime. If we see a positive sum it means
# 301 upstream anchored the root wrong and every downstream label is inverted.
# Stop hard rather than write misleading CSVs.
kupffer_anchor <- c("TIMD4", "MARCO", "VCAM1")
anchor_rows    <- pt_corr[gene %in% kupffer_anchor]
missing_anchor <- setdiff(kupffer_anchor, anchor_rows$gene)
if (length(missing_anchor) > 0) {
  stop(sprintf(
    "Kupffer anchor markers missing from pseudotime_corr_Macrophages.csv: %s",
    paste(missing_anchor, collapse = ", ")
  ))
}
anchor_sum_rho <- sum(anchor_rows$spearman_rho, na.rm = TRUE)
message(sprintf(
  "  Kupffer anchor check: sum(rho) for TIMD4+MARCO+VCAM1 = %.3f (expect < 0)",
  anchor_sum_rho
))
if (anchor_sum_rho > 0) {
  stop(sprintf(paste0(
    "Kupffer->LAM pseudotime direction check FAILED: ",
    "sum(rho) for TIMD4+MARCO+VCAM1 = %.3f > 0. ",
    "Kupffer markers should anti-correlate with pseudotime (resident KC = ",
    "pseudotime 0). Re-run 301_palantir_dpt.py after confirming Kupffer ",
    "root anchor (KUPFFER_ROOT_MARKERS argmax) is active."
  ), anchor_sum_rho))
}

# Marker panels (curated)
kupffer_markers <- c("CLEC4F","TIMD4","MARCO","VSIG4","CD5L","CD163L1","FCN1",
                     "IL18","CETP","VCAM1","LYVE1","CD163","ITGAX")
lam_markers <- c("TREM2","CD9","GPNMB","SPP1","LIPA","LGALS3","FABP5","CTSB","CTSD",
                 "CTSL","APOE","APOC1","PLIN2","CD36")
m1_markers <- c("CCL2","TNF","IL6","IL1B","CXCL9","CXCL10","NOS2","IFIT1","MX1","IRF1")
m2_markers <- c("CD163","MRC1","STAB1","MARCO","CD206","IL10","TGFB1","IGF1")

score_panel <- function(markers, name) {
  pt <- pt_corr[gene %in% markers]
  if (nrow(pt) == 0) return(NULL)
  pt[, panel := name]
  pt
}
panels <- rbind(
  score_panel(kupffer_markers, "Kupffer"),
  score_panel(lam_markers, "LAM"),
  score_panel(m1_markers, "M1_inflammatory"),
  score_panel(m2_markers, "M2_reparative"),
  fill = TRUE
)
fwrite(panels, file.path(OUTDIR, "macrophage_marker_panels_pseudotime.csv"))

# Panel summary stats
panel_stats <- panels[, .(n_genes = .N,
                          mean_rho = mean(spearman_rho, na.rm = TRUE),
                          median_rho = median(spearman_rho, na.rm = TRUE),
                          n_positive = sum(spearman_rho > 0 & padj < 0.05),
                          n_negative = sum(spearman_rho < 0 & padj < 0.05)),
                      by = panel]

message("[2] Overlay LIANA Hep->Mac ligands...")
hep_mac <- fread(file.path(CCC, "D2_hep_Mac_bidirectional.csv"))
hep_to_mac <- hep_mac[direction == "Hep_to_Mac"]
# Receptor on macrophage: check their pseudotime correlation
receptor_pt <- merge(unique(hep_to_mac[, .(receptor)]),
                     pt_corr, by.x = "receptor", by.y = "gene", all.x = TRUE)
receptor_pt <- receptor_pt[!is.na(spearman_rho)]
# Ligand on hepatocyte (won't be in macrophage pseudotime), but we can
# correlate hepatocyte ligand expression with macrophage trajectory via bulk.

# For each LR pair: macrophage receptor pseudotime correlation tells us if
# receptor-expressing macrophages are "early" or "late" in the Kupffer->LAM axis
hep_to_mac_enriched <- merge(hep_to_mac[score_diff > 0.1],
                             pt_corr[, .(gene, receptor_rho = spearman_rho,
                                          receptor_padj = padj)],
                             by.x = "receptor", by.y = "gene", all.x = TRUE)
fwrite(hep_to_mac_enriched, file.path(OUTDIR, "hep_to_mac_ligands_x_pseudotime_receptors.csv"))

# LAM-associated receptors (receptor_rho > 0, sig)
lam_associated <- hep_to_mac_enriched[receptor_rho > 0 & receptor_padj < 0.05]
kupffer_associated <- hep_to_mac_enriched[receptor_rho < 0 & receptor_padj < 0.05]

message("[3] Summary + writing outputs...")
summary_lines <- c(
  sprintf("Macrophage pseudotime correlated genes: %d", nrow(pt_corr)),
  sprintf("Genes in panels: Kupffer=%d, LAM=%d, M1=%d, M2=%d",
          length(kupffer_markers), length(lam_markers),
          length(m1_markers), length(m2_markers)),
  "",
  "=== Panel pseudotime summary (expect Kupffer negative, LAM positive) ===",
  capture.output(print(panel_stats)),
  "",
  sprintf("Hep->Mac MASLD-enriched LR pairs with receptor on macrophage: %d",
          nrow(hep_to_mac_enriched)),
  sprintf("  LAM-associated (receptor_rho>0, padj<0.05): %d", nrow(lam_associated)),
  sprintf("  Kupffer-associated (receptor_rho<0, padj<0.05): %d", nrow(kupffer_associated)),
  "",
  "=== Top 20 LAM-associated Hep->Mac axes ===",
  capture.output(print(lam_associated[order(-score_diff)][1:min(20,.N),
                       .(source, target, ligand, receptor, score_diff,
                         receptor_rho, receptor_padj)], nrows = 20)),
  "",
  "=== Top 10 Kupffer-associated Hep->Mac axes ===",
  capture.output(print(kupffer_associated[order(-score_diff)][1:min(10,.N),
                       .(source, target, ligand, receptor, score_diff,
                         receptor_rho, receptor_padj)], nrows = 10))
)
writeLines(summary_lines, file.path(OUTDIR, "kupffer_to_lam_trajectory_summary.txt"))
writeLines(summary_lines)

message("Done. Outputs in: ", OUTDIR)
