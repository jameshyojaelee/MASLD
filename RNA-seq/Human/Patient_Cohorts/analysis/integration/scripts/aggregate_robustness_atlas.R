#!/usr/bin/env Rscript
# aggregate_robustness_atlas.R
# ---------------------------------------------------------------------------
# Master per-gene robustness atlas. Joins per-pillar outputs into a single
# CSV that the eventual figures + main-text Methods table read from. This is
# tier-agnostic — the post-hoc "high-confidence" decision is made in figures.
#
# Inputs (all under audit_sensitivity/):
#   pillar_A_stability.csv          (CPSS + bootstrap + kfold)
#   pillar_C_empirical_null.csv     (permutation z + empirical p)
#   pillar_D_residual_varpart.csv   (per-gene residual variance fractions)
#   pillar_D_sva_concordance.csv    (per-gene SVA concordance flag)
# Plus:
#   integration/loo_cv/loo_cv_per_gene.csv  (5-fold LOCO recurrence)
#   integration/variance_partition.csv      (raw cohort variance fraction)
#   integration/dream_results.csv           (canonical DEG flag)
#
# Output:
#   audit_sensitivity/dream_robustness_atlas.csv
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({ library(data.table) })

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
RDIR <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
OUT_DIR <- file.path(BASE, "RNA-seq/results/audit_sensitivity")

PADJ_THR <- 0.05; LFC_THR <- 0.5

# --- Canonical full dream ---
full <- fread(file.path(RDIR, "dream_results.csv"))
atlas <- full[, .(gene,
                  dream_logFC = round(logFC, 4),  # C2-OK-sensitivity (dream robustness battery)
                  dream_padj  = signif(padj, 4),  # C2-OK-sensitivity (dream robustness battery)
                  is_canonical_DEG = (padj < PADJ_THR & abs(logFC) > LFC_THR))]
cat("Atlas seed rows:", nrow(atlas), " canonical DEGs:", sum(atlas$is_canonical_DEG), "\n")

# --- Pillar A ---
pA_path <- file.path(OUT_DIR, "pillar_A_stability.csv")
if (file.exists(pA_path)) {
  pA <- fread(pA_path)
  cols <- intersect(c("gene", "cpss_pi_hat", "cpss_pfer_at_pi_hat",
                      "bootstrap_freq", "bootstrap_logFC_mean", "bootstrap_logFC_sd",
                      "bootstrap_logFC_q025", "bootstrap_logFC_q975",
                      "kfold10_recur", "kfold50_recur"), names(pA))
  atlas <- merge(atlas, pA[, ..cols], by = "gene", all.x = TRUE)
} else {
  cat("WARN: pillar_A_stability.csv not found\n")
}

# --- Pillar C ---
pC_path <- file.path(OUT_DIR, "pillar_C_empirical_null.csv")
if (file.exists(pC_path)) {
  pC <- fread(pC_path, select = c("gene", "perm_z", "perm_emp_p"))
  atlas <- merge(atlas, pC, by = "gene", all.x = TRUE)
} else {
  cat("WARN: pillar_C_empirical_null.csv not found\n")
}

# --- Pillar D residual PVCA ---
pD_path <- file.path(OUT_DIR, "pillar_D_residual_varpart.csv")
if (file.exists(pD_path)) {
  pD <- fread(pD_path)
  if ("dataset" %in% names(pD)) {
    pD_thin <- pD[, .(gene,
                      residual_cohort_var = round(dataset, 4),
                      residual_disease_var = if ("group_binary" %in% names(pD))
                        round(group_binary, 4) else NA_real_)]
    atlas <- merge(atlas, pD_thin, by = "gene", all.x = TRUE)
  }
} else {
  cat("WARN: pillar_D_residual_varpart.csv not found\n")
}

# --- Pillar D SVA concordance ---
pD_sva_path <- file.path(OUT_DIR, "pillar_D_sva_concordance.csv")
if (file.exists(pD_sva_path)) {
  pDs <- fread(pD_sva_path, select = c("gene", "sva_concordant", "logFC_sva", "padj_sva"))
  setnames(pDs, c("logFC_sva", "padj_sva"), c("dream_sva_logFC", "dream_sva_padj"))
  pDs[, dream_sva_logFC := round(dream_sva_logFC, 4)]
  pDs[, dream_sva_padj  := signif(dream_sva_padj, 4)]
  atlas <- merge(atlas, pDs, by = "gene", all.x = TRUE)
} else {
  cat("WARN: pillar_D_sva_concordance.csv not found\n")
}

# --- Existing 5-fold LOCO recurrence ---
loco_path <- file.path(RDIR, "loo_cv/loo_cv_per_gene.csv")
if (file.exists(loco_path)) {
  loco <- fread(loco_path, select = c("gene", "n_loo_sig", "n_loo_tested", "robustness"))
  setnames(loco,
           c("n_loo_sig", "n_loo_tested", "robustness"),
           c("loco_n_sig", "loco_n_tested", "loco_robustness_class"))
  atlas <- merge(atlas, loco, by = "gene", all.x = TRUE)
} else {
  cat("WARN: loo_cv_per_gene.csv not found\n")
}

# --- Raw variance partition (cohort fraction) ---
vp_path <- file.path(RDIR, "variance_partition.csv")
if (file.exists(vp_path)) {
  vp <- fread(vp_path)
  if ("dataset" %in% names(vp)) {
    vp_thin <- vp[, .(gene, varpart_cohort_raw = round(dataset, 4),
                      varpart_disease_raw = if ("condition" %in% names(vp)) round(condition, 4) else NA_real_)]
    atlas <- merge(atlas, vp_thin, by = "gene", all.x = TRUE)
  }
} else {
  cat("WARN: variance_partition.csv not found\n")
}

# --- Derived n_pillars_passed ---
# Per-gene pillar pass: A passes if cpss_pi_hat>=0.7 AND bootstrap_freq>=0.9
#                       C passes if abs(perm_z)>=3
#                       D passes if residual_cohort_var<0.05 AND sva_concordant
# (Pillar B is signature-level, not per-gene; not counted here.)
# Ensure expected columns exist even if upstream pillar files are missing.
for (col in c("cpss_pi_hat", "bootstrap_freq", "perm_z",
              "residual_cohort_var", "sva_concordant")) {
  if (!col %in% names(atlas)) atlas[, (col) := NA]
}
atlas[, A_pass := !is.na(cpss_pi_hat) & cpss_pi_hat >= 0.7 &
                  !is.na(bootstrap_freq) & bootstrap_freq >= 0.9]
atlas[, C_pass := !is.na(perm_z) & abs(perm_z) >= 3]
atlas[, D_pass := !is.na(residual_cohort_var) & residual_cohort_var < 0.05 &
                  !is.na(sva_concordant) & sva_concordant == TRUE]
atlas[, n_pillars_passed := as.integer(A_pass) + as.integer(C_pass) + as.integer(D_pass)]

# Tier definition deferred — this is a placeholder column the figure script
# fills based on the post-hoc threshold decision.
atlas[, is_high_confidence := NA]

fwrite(atlas, file.path(OUT_DIR, "dream_robustness_atlas.csv"))
cat(sprintf("\nSaved master atlas: %d rows × %d columns\n", nrow(atlas), ncol(atlas)))

# Headline reporting
cat("\n=== Master atlas summary ===\n")
n_canon <- sum(atlas$is_canonical_DEG)
cat(sprintf("Canonical DEGs (4,370 expected, kallisto canonical): %d\n", n_canon))
deg <- atlas[is_canonical_DEG == TRUE]
if (n_canon > 0) {
  cat(sprintf("  Pillar A pass: %d (%.1f%%)\n", sum(deg$A_pass, na.rm = TRUE),
              100 * mean(deg$A_pass, na.rm = TRUE)))
  cat(sprintf("  Pillar C pass: %d (%.1f%%)\n", sum(deg$C_pass, na.rm = TRUE),
              100 * mean(deg$C_pass, na.rm = TRUE)))
  cat(sprintf("  Pillar D pass: %d (%.1f%%)\n", sum(deg$D_pass, na.rm = TRUE),
              100 * mean(deg$D_pass, na.rm = TRUE)))
  cat(sprintf("  All 3 pass:    %d (%.1f%%)\n", sum(deg$n_pillars_passed == 3, na.rm = TRUE),
              100 * mean(deg$n_pillars_passed == 3, na.rm = TRUE)))
}
cat("\nDone master atlas.\n")
