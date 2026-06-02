#!/usr/bin/env Rscript
# ============================================================================
# 07f_canonical_lr_audit.R
#
# Master review M-P0-15: audit the canonical MASLD ligand-receptor signaling
# pairs that a Cell Metabolism reviewer will expect to see explicitly tested.
# Reports their LMM q-values across all four stage axes (coarse, documented
# F-stage, augmented F-stage, continuous macrophage pseudotime) so the manuscript
# can state whether they were tested and what their stage-LMM q-values were.
#
# Pairs audited (curated from Marra & Svegliati-Baroni 2018, Schwabe & Brenner
# 2014, Krenkel & Tacke 2017):
#   - TGFB1 -> TGFBR1   (canonical fibrogenic ligand-receptor)
#   - TGFB1 -> TGFBR2
#   - SPP1  -> CD44     (osteopontin, macrophage-derived fibrosis driver)
#   - SPP1  -> ITGB1    (alternative SPP1 receptor)
#   - CCL2  -> CCR2     (monocyte recruitment, MASLD-cardinal)
#   - PDGFA -> PDGFRA   (HSC activation)
#   - PDGFB -> PDGFRB
#   - IL6   -> IL6R     (inflammatory, JAK/STAT)
#   - TNF   -> TNFRSF1A (NASH cardinal cytokine)
#
# Output: stage_trajectory_v3/canonical_lr_audit.tsv
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
V3_DIR <- file.path(BASE,
  "Analysis/SingleCell/results_gpu_v2_phase05/ccc/stage_trajectory_v3")

CANONICAL_PAIRS <- data.table(
  ligand   = c("TGFB1", "TGFB1", "SPP1", "SPP1",
               "CCL2", "PDGFA", "PDGFB", "IL6", "TNF"),
  receptor = c("TGFBR1", "TGFBR2", "CD44", "ITGB1",
               "CCR2", "PDGFRA", "PDGFRB", "IL6R", "TNFRSF1A"),
  reference = c("Schwabe 2014", "Schwabe 2014",
                "Krenkel 2017", "Krenkel 2017",
                "Marra 2018", "Tsuchida 2017", "Tsuchida 2017",
                "Marra 2018", "Marra 2018")
)
CANONICAL_PAIRS[, lr_pair := paste(ligand, receptor, sep = "__")]

cat(sprintf("[audit] %d canonical MASLD LR pairs to look up\n",
            nrow(CANONICAL_PAIRS)))

axis_files <- list(
  coarse        = "stage_lr_lmm_coarse_v3.tsv",
  documented    = "stage_lr_lmm_fstage_documented_v3.tsv",
  augmented     = "stage_lr_lmm_fstage_augmented_v3.tsv",
  continuous    = "stage_lr_lmm_continuous_v3.tsv"
)

results <- list()
for (axis_name in names(axis_files)) {
  fpath <- file.path(V3_DIR, axis_files[[axis_name]])
  if (!file.exists(fpath)) {
    cat(sprintf("[warn] %s missing -> skipping\n", fpath))
    next
  }
  dt <- fread(fpath)
  if (!all(c("ligand_complex", "receptor_complex") %in% names(dt))) {
    cat(sprintf("[warn] %s missing LR columns\n", fpath))
    next
  }
  # Match on un-complexed ligand/receptor symbols
  dt[, ligand := sub("_.*", "", ligand_complex)]
  dt[, receptor := sub("_.*", "", receptor_complex)]
  dt[, lr_pair := paste(ligand, receptor, sep = "__")]
  hits <- dt[lr_pair %in% CANONICAL_PAIRS$lr_pair]
  hits[, axis := axis_name]
  # For each canonical pair, take its best-supporting ct_pair (smallest q)
  best <- hits[order(q_bonferroni_family),
               .SD[1],
               by = .(lr_pair, axis)]
  results[[axis_name]] <- best
}

if (length(results) == 0) {
  stop("[audit] no axis files found")
}

audit <- rbindlist(results, fill = TRUE)
# Reshape: one row per canonical pair, columns per axis
out <- merge(CANONICAL_PAIRS,
             audit[, .(lr_pair, axis, ct_pair, Estimate, pval, q_within_ct,
                       q_bonferroni_family, n_donors)],
             by = "lr_pair", all.x = TRUE)
out[, tested := !is.na(pval)]
setorder(out, ligand, receptor, axis)

OUT_TSV <- file.path(V3_DIR, "canonical_lr_audit.tsv")
fwrite(out, OUT_TSV, sep = "\t")
cat(sprintf("\n[output] -> %s\n", OUT_TSV))

cat("\n[summary] canonical MASLD LR pair audit:\n")
summary_tab <- out[, .(
  n_axes_tested = sum(tested),
  n_axes_significant = sum(!is.na(q_bonferroni_family) &
                           q_bonferroni_family < 0.05),
  best_q = min(q_bonferroni_family, na.rm = TRUE),
  best_axis = axis[which.min(q_bonferroni_family)]
), by = .(ligand, receptor)]
print(summary_tab)

cat("\n[interpretation guide]\n")
cat("  - tested but not significant (q > 0.05):  expected null in MASLD scRNA\n")
cat("  - NOT tested at all (n_axes_tested == 0): pair absent from LIANA universe\n")
cat("  - significant in 1+ axes: should be reported in methods/results\n")
cat("[done]\n")
