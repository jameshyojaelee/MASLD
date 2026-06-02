#!/usr/bin/env Rscript
# verify_robustness.R
# ---------------------------------------------------------------------------
# Sanity-check verifier for the 4-pillar robustness pipeline. Run AFTER each
# stage to fail-fast if expectations are violated. Mirrors the verification
# section of the design spec.
#
# Usage:
#   Rscript verify_robustness.R [stage]
# where stage is one of: pillar_A | pillar_B | pillar_C | pillar_D | atlas | all
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({ library(data.table) })

args <- commandArgs(trailingOnly = TRUE)
stage <- if (length(args) > 0) args[1] else "all"

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
RDIR <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
AUDIT <- file.path(BASE, "RNA-seq/results/audit_sensitivity")

# Known canonical MASLD genes — should be selected almost always under any
# reasonable subsampling. If freq << 0.5 for these, something is broken.
# (Top fibrosis / inflammation markers from the existing dream output.)
SENTINEL_GENES <- c("COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL6A3",
                    "PDGFRA", "PDGFRB", "THBS2", "IGFBP3", "LUM",
                    "SPP1", "MMP2", "VIM", "TIMP1", "ACTA2")

check <- function(condition, ok_msg, fail_msg) {
  if (isTRUE(condition)) cat("  ✓ ", ok_msg, "\n")
  else                   cat("  ✗ ", fail_msg, "\n")
  invisible(condition)
}

# ============================================================================
# Pillar A
# ============================================================================
verify_A <- function() {
  cat("\n────────────────  Pillar A — sample stability  ────────────────\n")
  f <- file.path(AUDIT, "pillar_A_stability.csv")
  if (!file.exists(f)) { cat("  pillar_A_stability.csv not found\n"); return(invisible()) }
  pA <- fread(f)
  cat(sprintf("  Genes: %d  |  canonical DEGs: %d\n",
              nrow(pA), sum(pA$is_canonical_DEG, na.rm = TRUE)))
  deg <- pA[is_canonical_DEG == TRUE]

  # Bimodality of CPSS pi-hat: most genes near 0, tail near 1.
  # The mean of canonical DEG pi-hats should be substantially > non-DEG.
  mean_deg_pi <- mean(deg$cpss_pi_hat, na.rm = TRUE)
  mean_non_pi <- mean(pA[is_canonical_DEG == FALSE]$cpss_pi_hat, na.rm = TRUE)
  check(!is.na(mean_deg_pi) && !is.na(mean_non_pi) && (mean_deg_pi - mean_non_pi) > 0.3,
        sprintf("CPSS π̂ separates DEG vs non-DEG (DEG mean=%.3f vs non-DEG=%.3f)",
                mean_deg_pi, mean_non_pi),
        sprintf("CPSS π̂ does NOT separate DEG vs non-DEG (DEG=%.3f, non-DEG=%.3f)",
                mean_deg_pi, mean_non_pi))

  # Bootstrap freq for sentinel genes — should be almost always selected
  sn <- pA[gene %in% SENTINEL_GENES, .(gene, bootstrap_freq, cpss_pi_hat)]
  if (nrow(sn) > 0) {
    n_high <- sum(sn$bootstrap_freq >= 0.9, na.rm = TRUE)
    cat(sprintf("  Sentinel genes (n=%d found): %d / %d with bootstrap freq >= 0.9\n",
                nrow(sn), n_high, nrow(sn)))
    check(n_high >= ceiling(0.6 * nrow(sn)),
          "Sentinel genes pass bootstrap threshold (>= 60%)",
          "Sentinel genes FAIL bootstrap threshold — investigate subsampling")
    print(sn[order(-bootstrap_freq)])
  }

  # Pi-hat distribution shape — bimodality test (rough)
  pi_vec <- pA$cpss_pi_hat[!is.na(pA$cpss_pi_hat)]
  if (length(pi_vec) > 0) {
    n_low  <- mean(pi_vec < 0.2)
    n_high <- mean(pi_vec > 0.8)
    n_mid  <- mean(pi_vec >= 0.4 & pi_vec <= 0.6)
    cat(sprintf("  π̂ distribution: %.1f%% < 0.2, %.1f%% > 0.8, %.1f%% in [0.4, 0.6]\n",
                100 * n_low, 100 * n_high, 100 * n_mid))
    check(n_low + n_high > 0.5,
          "π̂ distribution is bimodal (mass at extremes)",
          "π̂ distribution is unimodal — possible subsampling bug")
  }

  # Pillar A iter summary
  iter_f <- file.path(AUDIT, "pillar_A_iter_summary.csv")
  if (file.exists(iter_f)) {
    iter <- fread(iter_f)
    cat(sprintf("  Iter summary rows: %d\n", nrow(iter)))
    by_src <- iter[, .(mean_rho = round(mean(rho, na.rm = TRUE), 3),
                        median_rho = round(median(rho, na.rm = TRUE), 3),
                        mean_jacc = round(mean(jaccard, na.rm = TRUE), 3),
                        mean_recovery = round(mean(recovery_pct, na.rm = TRUE), 1),
                        n = .N), by = source]
    print(by_src)
    check(all(by_src$mean_rho > 0.7, na.rm = TRUE),
          "All sources have mean ρ > 0.7",
          "Some sources have mean ρ <= 0.7 — fold sample sizes too small?")
  }
}

# ============================================================================
# Pillar B
# ============================================================================
verify_B <- function() {
  cat("\n────────────────  Pillar B — transferability  ────────────────\n")
  f <- file.path(AUDIT, "pillar_B_loco_prediction.csv")
  if (!file.exists(f)) { cat("  pillar_B_loco_prediction.csv not found\n"); return(invisible()) }
  pB <- fread(f)
  cat(sprintf("  Rows: %d (expected = 2 × n_held_out)\n", nrow(pB)))

  m <- pB[, .(mean_auroc = round(mean(auroc, na.rm = TRUE), 4),
              min_auroc  = round(min(auroc, na.rm = TRUE), 4),
              mean_null  = round(mean(null_mean, na.rm = TRUE), 4),
              mean_z     = round(mean(z_vs_null, na.rm = TRUE), 2)),
          by = method]
  print(m)
  check(all(m$mean_auroc >= 0.75),
        sprintf("Mean AUROC >= 0.75 for all methods (mean across methods = %.3f)",
                mean(m$mean_auroc, na.rm = TRUE)),
        "Mean AUROC < 0.75 — transferability is weak; check signature size or scoring method")
  check(all(m$mean_null < 0.6),
        sprintf("Random-label null mean < 0.6 (mean across methods = %.3f)",
                mean(m$mean_null, na.rm = TRUE)),
        "Null AUROC > 0.6 — class-imbalance leakage suspected; double-check shuffle code")
  check(all(m$mean_z >= 3),
        sprintf("Mean z vs null >= 3 (mean across methods = %.2f)", mean(m$mean_z, na.rm = TRUE)),
        "Mean z vs null < 3 — observed AUROC indistinguishable from random")
  check(all(pB$auroc >= 0.65),
        "All per-cohort AUROCs >= 0.65",
        sprintf("Some per-cohort AUROCs below 0.65 (min = %.3f) — check the weakest fold",
                min(pB$auroc, na.rm = TRUE)))
}

# ============================================================================
# Pillar C
# ============================================================================
verify_C <- function() {
  cat("\n────────────────  Pillar C — empirical null  ────────────────\n")
  f_sum <- file.path(AUDIT, "pillar_C_summary.csv")
  if (!file.exists(f_sum)) { cat("  pillar_C_summary.csv not found\n"); return(invisible()) }
  s <- fread(f_sum)
  print(s)
  check(s$empirical_fdr[1] < 0.05,
        sprintf("Empirical FDR < 0.05 (= %.4f)", s$empirical_fdr[1]),
        sprintf("Empirical FDR >= 0.05 (= %.4f) — parametric FDR liberal; flag for rebuttal", s$empirical_fdr[1]))
  check(s$pct_canonical_perm_z_gt_3[1] >= 90,
        sprintf("%.1f%% canonical DEGs have |perm z| >= 3", s$pct_canonical_perm_z_gt_3[1]),
        sprintf("Only %.1f%% canonical DEGs have |perm z| >= 3 (target >= 90%%)",
                s$pct_canonical_perm_z_gt_3[1]))
  cat(sprintf("  Permutation count distribution: median=%.0f mean=%.0f q95=%.0f (vs %d observed)\n",
              s$perm_count_median[1], s$perm_count_mean[1],
              s$perm_count_q95[1], s$n_deg_observed[1]))
}

# ============================================================================
# Pillar D
# ============================================================================
verify_D <- function() {
  cat("\n────────────────  Pillar D — batch-RE diagnostics  ────────────────\n")

  # PVCA on residuals
  f_pvca <- file.path(AUDIT, "pillar_D_residual_varpart_summary.csv")
  if (file.exists(f_pvca)) {
    pvca <- fread(f_pvca)
    cat("  Residual variance partition:\n")
    print(pvca)
    cohort_med <- pvca[component == "dataset", median_var]
    if (length(cohort_med) > 0 && !is.na(cohort_med)) {
      check(cohort_med < 0.05,
            sprintf("Cohort residual variance < 5%% (median = %.3f)", cohort_med),
            sprintf("Cohort residual variance = %.3f >= 0.05 — random effect under-fit", cohort_med))
    }
  } else {
    cat("  pillar_D_residual_varpart_summary.csv not found\n")
  }

  # SVA concordance
  f_sva <- file.path(AUDIT, "pillar_D_sva_summary.csv")
  if (file.exists(f_sva)) {
    sva <- fread(f_sva)
    cat("  SVA summary:\n"); print(sva)
    check(!is.na(sva$jaccard[1]) && sva$jaccard[1] >= 0.85,
          sprintf("dream ∩ dream+SVA Jaccard = %.3f (>= 0.85)", sva$jaccard[1]),
          sprintf("Jaccard = %.3f < 0.85 — hidden batch present; reframe primary list as conservative subset",
                  sva$jaccard[1]))
    check(!is.na(sva$rho_logFC[1]) && sva$rho_logFC[1] >= 0.95,
          sprintf("logFC ρ between dream and dream+SVA = %.3f (>= 0.95)", sva$rho_logFC[1]),
          sprintf("logFC ρ = %.3f below 0.95 — sizable SV-driven shifts in effect estimates",
                  sva$rho_logFC[1]))
  } else {
    cat("  pillar_D_sva_summary.csv not found\n")
  }
}

# ============================================================================
# Master atlas
# ============================================================================
verify_atlas <- function() {
  cat("\n────────────────  Master atlas  ────────────────\n")
  f <- file.path(AUDIT, "dream_robustness_atlas.csv")
  if (!file.exists(f)) { cat("  dream_robustness_atlas.csv not found\n"); return(invisible()) }
  atlas <- fread(f)
  cat(sprintf("  Atlas dimensions: %d genes × %d cols\n", nrow(atlas), ncol(atlas)))
  cat(sprintf("  Canonical DEGs: %d\n", sum(atlas$is_canonical_DEG, na.rm = TRUE)))
  deg <- atlas[is_canonical_DEG == TRUE]
  if (nrow(deg) > 0) {
    cat(sprintf("  Pillar pass-summary on canonical DEGs:\n"))
    cat(sprintf("    A pass: %d (%.1f%%)\n",
                sum(deg$A_pass, na.rm = TRUE), 100 * mean(deg$A_pass, na.rm = TRUE)))
    cat(sprintf("    C pass: %d (%.1f%%)\n",
                sum(deg$C_pass, na.rm = TRUE), 100 * mean(deg$C_pass, na.rm = TRUE)))
    cat(sprintf("    D pass: %d (%.1f%%)\n",
                sum(deg$D_pass, na.rm = TRUE), 100 * mean(deg$D_pass, na.rm = TRUE)))
    cat(sprintf("    All 3 pass: %d (%.1f%%)\n",
                sum(deg$n_pillars_passed == 3, na.rm = TRUE),
                100 * mean(deg$n_pillars_passed == 3, na.rm = TRUE)))
  }
}

if (stage %in% c("pillar_A", "all")) verify_A()
if (stage %in% c("pillar_B", "all")) verify_B()
if (stage %in% c("pillar_C", "all")) verify_C()
if (stage %in% c("pillar_D", "all")) verify_D()
if (stage %in% c("atlas", "all"))    verify_atlas()
cat("\nDone.\n")
