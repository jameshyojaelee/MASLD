#!/usr/bin/env Rscript
# 05b_ashr_shrinkage.R
# ---------------------------------------------------------------------------
# Apply adaptive shrinkage (ashr) to dream mega-analysis effect sizes.
# Produces shrunk logFC estimates and local false sign rates (lfsr).
# This replaces arbitrary logFC cutoffs with data-driven effect size estimation.
#
# Input:  results/integration/dream_results.csv
# Output: results/integration/dream_results_ashr.csv
#         results/integration/ashr_diagnostics.csv
#
# Runs AFTER 05_dream_mega_analysis.R, BEFORE 07_consensus_degs.R
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ashr)
})

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
INT  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR <- file.path(INT, "results/integration")

cat("============================================================\n")
cat("05b: Adaptive shrinkage of dream effect sizes\n")
cat("============================================================\n\n")

# ── 1. Load dream results ────────────────────────────────────────────────────
dream <- fread(file.path(RDIR, "dream_results.csv"))
cat("Loaded:", nrow(dream), "genes from dream_results.csv\n")
cat("Columns:", paste(names(dream), collapse = ", "), "\n")

# ── 2. Use pre-computed SE from dream output ─────────────────────────────────
# Phase 1.4: dream now emits SE directly (moderated SE from Satterthwaite)
if ("SE" %in% names(dream)) {
  dream[, se := SE]
  cat("Using pre-computed SE column from dream output\n")
} else {
  dream[, se := abs(logFC / t)]
  cat("SE column not found; falling back to |logFC / t|\n")
}

# Handle edge cases
n_bad <- sum(is.na(dream$se) | is.infinite(dream$se) | dream$se <= 0)
if (n_bad > 0) {
  cat("WARNING:", n_bad, "genes with invalid SE (t=0 or NA) — excluded from shrinkage\n")
}
valid <- !is.na(dream$se) & is.finite(dream$se) & dream$se > 0

# ── 3. Run ashr ──────────────────────────────────────────────────────────────
cat("\nRunning ashr adaptive shrinkage on", sum(valid), "genes...\n")

ash_fit <- ash(
  betahat    = dream$logFC[valid],
  sebetahat  = dream$se[valid],
  mixcompdist = "halfuniform",
  method      = "shrink"
)

# Extract results
dream[, shrunk_logFC := NA_real_]
dream[, shrunk_se    := NA_real_]
dream[, lfsr         := NA_real_]
dream[, svalue       := NA_real_]

dream[valid, shrunk_logFC := ash_fit$result$PosteriorMean]
dream[valid, shrunk_se    := ash_fit$result$PosteriorSD]
dream[valid, lfsr         := ash_fit$result$lfsr]
dream[valid, svalue       := ash_fit$result$svalue]

cat("ashr complete.\n")

# ── 4. Summary statistics ────────────────────────────────────────────────────
cat("\n--- Shrinkage summary ---\n")
cat("  Genes with lfsr < 0.05:", sum(dream$lfsr < 0.05, na.rm = TRUE), "\n")
cat("  Genes with lfsr < 0.01:", sum(dream$lfsr < 0.01, na.rm = TRUE), "\n")

# Quantiles of |shrunk_logFC| among significant genes
sig <- dream[lfsr < 0.05]
cat("\n  |shrunk_logFC| quantiles (lfsr < 0.05, n=", nrow(sig), "):\n")
q <- quantile(abs(sig$shrunk_logFC), probs = c(0, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 1.0))
print(round(q, 4))

# Correlation between raw and shrunk
rho <- cor(dream$logFC[valid], dream$shrunk_logFC[valid], method = "spearman", use = "complete.obs")
cat(sprintf("\n  Spearman rho (raw vs shrunk logFC): %.4f\n", rho))

# ── 5. Positive control calibration ──────────────────────────────────────────
cat("\n--- Positive control calibration ---\n")

ctrl_file <- file.path(BASE, "results/library/positive_control.csv")
if (file.exists(ctrl_file)) {
  controls <- fread(ctrl_file)
  ctrl_symbols <- controls$`Gene symbol`

  # Map gene symbols via atlas
  atlas_map <- fread(file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"),
                     select = c("ensembl_id", "human_symbol"))
  atlas_map[, ensembl_base := sub("\\.[0-9]+$", "", ensembl_id)]
  atlas_map <- unique(atlas_map[human_symbol != "", .(ensembl_base, symbol = human_symbol)])

  dream[, ensembl_base := sub("\\.[0-9]+$", "", gene)]
  dream <- merge(dream, atlas_map, by = "ensembl_base", all.x = TRUE)
  dream[, is_control := symbol %in% ctrl_symbols]

  n_ctrl <- sum(dream$is_control, na.rm = TRUE)
  cat("  Controls found:", n_ctrl, "/", length(ctrl_symbols), "\n")

  # Test threshold calibration
  thresholds <- seq(0.05, 0.5, by = 0.025)
  calibration <- data.table(
    threshold = thresholds,
    n_degs = sapply(thresholds, function(th) sum(dream$lfsr < 0.05 & abs(dream$shrunk_logFC) > th, na.rm = TRUE)),
    ctrl_recovery = sapply(thresholds, function(th) {
      ctrl_data <- dream[is_control == TRUE]
      sum(abs(ctrl_data$shrunk_logFC) > th, na.rm = TRUE) / nrow(ctrl_data)
    })
  )
  calibration[, stringency := 1 - n_degs / max(n_degs)]
  calibration[, f1 := 2 * ctrl_recovery * stringency / (ctrl_recovery + stringency + 1e-10)]

  best_th <- calibration[which.max(f1), threshold]
  cat(sprintf("  Best F1-calibrated threshold: %.3f\n", best_th))
  cat(sprintf("  DEGs at best threshold: %d\n", calibration[threshold == best_th, n_degs]))
  cat(sprintf("  Control recovery: %.1f%%\n", calibration[threshold == best_th, ctrl_recovery * 100]))

  fwrite(calibration, file.path(RDIR, "ashr_threshold_calibration.csv"))
} else {
  cat("  Positive control file not found — skipping calibration\n")
  best_th <- 0.2  # default
}

# ── 6. Recommended threshold ─────────────────────────────────────────────────
# Use the 5th percentile of significant shrunk |logFC| as floor, calibrated against controls
empirical_floor <- round(q["5%"], 3)
recommended <- max(empirical_floor, best_th, na.rm = TRUE)

n_deg_new <- sum(dream$lfsr < 0.05 & abs(dream$shrunk_logFC) > recommended, na.rm = TRUE)
n_deg_old <- sum(dream$padj < 0.1, na.rm = TRUE)

cat("\n============================================================\n")
cat("RECOMMENDED logFC THRESHOLD:", recommended, "\n")
cat("  New DEG count (lfsr<0.05 + |shrunk_logFC|>", recommended, "):", n_deg_new, "\n")
cat("  Old DEG count (padj<0.1):", n_deg_old, "\n")
cat("  Reduction:", round(100 * (1 - n_deg_new / n_deg_old), 1), "%\n")
cat("============================================================\n")

# ── 7. Save outputs ──────────────────────────────────────────────────────────
# Main output: dream results with ashr columns
out_cols <- c("gene", "logFC", "AveExpr", "t", "P.Value", "padj", "z.std",
              "se", "shrunk_logFC", "shrunk_se", "lfsr", "svalue")
if ("symbol" %in% names(dream)) out_cols <- c(out_cols, "symbol")
out_cols <- intersect(out_cols, names(dream))

fwrite(dream[, ..out_cols], file.path(RDIR, "dream_results_ashr.csv"))
cat("\nSaved: dream_results_ashr.csv (", nrow(dream), "genes x", length(out_cols), "cols)\n")

# Diagnostics summary
diag <- data.table(
  metric = c("total_genes", "valid_for_ashr", "lfsr_lt_005", "lfsr_lt_001",
             "recommended_threshold", "n_deg_new", "n_deg_old",
             "spearman_raw_vs_shrunk", "empirical_5pct_floor"),
  value = c(nrow(dream), sum(valid), sum(dream$lfsr < 0.05, na.rm = TRUE),
            sum(dream$lfsr < 0.01, na.rm = TRUE),
            recommended, n_deg_new, n_deg_old, round(rho, 4), empirical_floor)
)
fwrite(diag, file.path(RDIR, "ashr_diagnostics.csv"))
cat("Saved: ashr_diagnostics.csv\n")

cat("\nScript 05b complete.\n")
