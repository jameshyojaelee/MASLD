#!/usr/bin/env Rscript
# 218_ashr_lfc_threshold.R
# Principled logFC threshold for DEGs using adaptive shrinkage (ashr)
# + validation against 65 positive control genes
#
# ashr shrinks noisy logFC toward zero. Genes with |shrunk_logFC| > 0
# at lfsr < 0.05 are "truly" differentially expressed with meaningful effect.
# The empirical floor from this analysis replaces arbitrary |logFC| > 0.25 or 0.5.

library(data.table)
library(ashr)
library(ggplot2)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))

OUT_DIR <- file.path(BASE, "RNA-seq/results/stratified_causal")
FIG_DIR <- file.path(BASE, "figures/supplementary/figS_methods_validation/sensitivity")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

cat("============================================================\n")
cat("218: Adaptive shrinkage logFC threshold\n")
cat("============================================================\n\n")

# ── 1. Load dream results ────────────────────────────────────────────────────
dream <- fread(file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/dream_results.csv"))
cat("Dream results:", nrow(dream), "genes\n")

# Derive SE from logFC and t-statistic
dream[, se := abs(logFC / t)]
# Handle edge cases
dream[is.na(se) | is.infinite(se), se := NA]
dream <- dream[!is.na(se) & se > 0 & !is.na(logFC)]
cat("After SE filtering:", nrow(dream), "genes\n")

# ── 2. Run ashr ──────────────────────────────────────────────────────────────
cat("\n--- Running ashr adaptive shrinkage ---\n")
ash_fit <- ash(betahat = dream$logFC, sebetahat = dream$se,
               mixcompdist = "halfuniform", method = "shrink")

dream[, shrunk_logFC := ash_fit$result$PosteriorMean]
dream[, shrunk_se    := ash_fit$result$PosteriorSD]
dream[, lfsr         := ash_fit$result$lfsr]
dream[, svalue       := ash_fit$result$svalue]

cat("ashr complete\n")
cat("  Genes with lfsr < 0.05:", sum(dream$lfsr < 0.05), "\n")
cat("  Genes with lfsr < 0.01:", sum(dream$lfsr < 0.01), "\n")

# ── 3. Determine empirical threshold ─────────────────────────────────────────
# Among genes with lfsr < 0.05, what's the distribution of |shrunk_logFC|?
sig_shrunk <- dream[lfsr < 0.05]
cat("\n--- Shrunk logFC distribution (lfsr < 0.05) ---\n")
cat("  n:", nrow(sig_shrunk), "\n")
cat("  |shrunk_logFC| quantiles:\n")
q <- quantile(abs(sig_shrunk$shrunk_logFC), probs = c(0, 0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 1.0))
print(round(q, 4))

# The 5th percentile of |shrunk_logFC| among significant genes = empirical floor
# Genes below this are borderline even after shrinkage
empirical_lfc <- round(q["5%"], 3)
cat("\n  >>> Empirical logFC threshold (5th percentile of shrunk significant):", empirical_lfc, "\n")

# Also check: what's the min |shrunk_logFC| where lfsr transitions from <0.05 to >0.05?
# This is the "natural" ashr boundary
dream_sorted <- dream[order(abs(shrunk_logFC))]
boundary_idx <- which(dream_sorted$lfsr < 0.05)[1]
if (!is.na(boundary_idx)) {
  natural_boundary <- abs(dream_sorted$shrunk_logFC[boundary_idx])
  cat("  Natural ashr boundary (min |shrunk_logFC| at lfsr<0.05):", round(natural_boundary, 4), "\n")
}

# ── 4. Positive control validation ───────────────────────────────────────────
cat("\n--- Positive control validation ---\n")

controls <- fread(file.path(BASE, "results/library/positive_control.csv"))
cat("Positive controls:", nrow(controls), "genes\n")

# Map control gene symbols to dream Ensembl IDs via atlas
atlas <- fread(file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"),
               select = c("ensembl_id", "human_symbol"))
atlas[, ensembl_base := sub("\\.[0-9]+$", "", ensembl_id)]
dream[, ensembl_base := sub("\\.[0-9]+$", "", gene)]

# Join symbols
dream <- merge(dream, unique(atlas[human_symbol != "", .(ensembl_base, symbol = human_symbol)]),
               by = "ensembl_base", all.x = TRUE)

# Find controls in dream
control_symbols <- controls$`Gene symbol`
dream[, is_control := symbol %in% control_symbols]
n_controls_found <- sum(dream$is_control, na.rm = TRUE)
cat("Controls found in dream:", n_controls_found, "/", length(control_symbols), "\n")

if (n_controls_found > 0) {
  ctrl <- dream[is_control == TRUE]
  cat("\nControl gene |shrunk_logFC| distribution:\n")
  ctrl_q <- quantile(abs(ctrl$shrunk_logFC), probs = c(0, 0.05, 0.10, 0.25, 0.50, 0.75, 1.0), na.rm = TRUE)
  print(round(ctrl_q, 4))

  # What fraction of controls pass the empirical threshold?
  n_pass <- sum(abs(ctrl$shrunk_logFC) > empirical_lfc, na.rm = TRUE)
  cat(sprintf("\nControls passing |shrunk_logFC| > %.3f: %d / %d (%.1f%%)\n",
              empirical_lfc, n_pass, nrow(ctrl), 100 * n_pass / nrow(ctrl)))

  # Test a range of thresholds against control recovery
  thresholds <- seq(0.05, 0.5, by = 0.025)
  recovery <- sapply(thresholds, function(th) {
    sum(abs(ctrl$shrunk_logFC) > th, na.rm = TRUE) / nrow(ctrl)
  })

  # Also compute total DEGs at each threshold
  total_degs <- sapply(thresholds, function(th) {
    sum(dream$lfsr < 0.05 & abs(dream$shrunk_logFC) > th, na.rm = TRUE)
  })

  threshold_dt <- data.table(
    threshold = thresholds,
    control_recovery = recovery,
    n_degs = total_degs
  )

  # Find optimal: max F1-like balance between recovery and stringency
  threshold_dt[, stringency := 1 - n_degs / max(n_degs)]
  threshold_dt[, f1 := 2 * control_recovery * stringency / (control_recovery + stringency + 1e-10)]
  best_threshold <- threshold_dt[which.max(f1), threshold]

  cat("\n--- Threshold calibration ---\n")
  print(threshold_dt[, .(threshold, control_recovery = round(control_recovery, 3),
                          n_degs, f1 = round(f1, 3))])
  cat("\n  >>> Best F1-calibrated threshold:", best_threshold, "\n")
  cat("  >>> Controls recovered at best threshold:",
      round(threshold_dt[threshold == best_threshold, control_recovery] * 100, 1), "%\n")
  cat("  >>> DEGs at best threshold:",
      threshold_dt[threshold == best_threshold, n_degs], "\n")
}

# ── 5. Final recommendation ──────────────────────────────────────────────────
# Use the more conservative of: empirical_lfc vs best_threshold
recommended <- max(empirical_lfc, best_threshold, na.rm = TRUE)
cat("\n============================================================\n")
cat("RECOMMENDED logFC THRESHOLD:", recommended, "\n")
cat("  Based on: ashr shrinkage 5th percentile =", empirical_lfc, "\n")
cat("  Calibrated against:", n_controls_found, "positive controls, best F1 =", best_threshold, "\n")
cat("  DEGs at this threshold (lfsr<0.05 + |shrunk_logFC|>", recommended, "):",
    sum(dream$lfsr < 0.05 & abs(dream$shrunk_logFC) > recommended, na.rm = TRUE), "\n")
cat("  Compare: original padj<0.1 DEGs:", sum(dream$padj < 0.1, na.rm = TRUE), "\n")
cat("============================================================\n")

# ── 6. Save outputs ──────────────────────────────────────────────────────────
# Full dream with shrunk logFC
fwrite(dream[, .(gene, ensembl_base, symbol, logFC, se, t, padj,
                  shrunk_logFC, shrunk_se, lfsr, svalue, is_control)],
       file.path(OUT_DIR, "dream_ashr_shrunk.csv"))

# Threshold summary
fwrite(data.table(
  method = c("ashr_5pct", "control_f1_calibrated", "recommended"),
  threshold = c(empirical_lfc, best_threshold, recommended),
  n_degs = c(
    sum(dream$lfsr < 0.05 & abs(dream$shrunk_logFC) > empirical_lfc, na.rm = TRUE),
    sum(dream$lfsr < 0.05 & abs(dream$shrunk_logFC) > best_threshold, na.rm = TRUE),
    sum(dream$lfsr < 0.05 & abs(dream$shrunk_logFC) > recommended, na.rm = TRUE)
  )
), file.path(OUT_DIR, "ashr_lfc_threshold.csv"))

if (exists("threshold_dt")) {
  fwrite(threshold_dt, file.path(OUT_DIR, "ashr_threshold_calibration.csv"))
}

# ── 7. Figure: threshold calibration curve ───────────────────────────────────
if (exists("threshold_dt")) {
  p1 <- ggplot(threshold_dt, aes(x = threshold)) +
    geom_line(aes(y = control_recovery, color = "Control recovery"), linewidth = 0.6) +
    geom_line(aes(y = n_degs / max(n_degs), color = "DEG fraction"), linewidth = 0.6) +
    geom_line(aes(y = f1, color = "F1 score"), linewidth = 0.8, linetype = "dashed") +
    geom_vline(xintercept = recommended, linetype = "dotted", color = "grey30") +
    annotate("text", x = recommended + 0.01, y = 0.95,
             label = paste0("threshold = ", recommended), hjust = 0, size = 2.5) +
    scale_color_manual(values = c("Control recovery" = "#C2185B",
                                   "DEG fraction" = "#1565C0",
                                   "F1 score" = "#00695C"),
                       name = NULL) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = "|shrunk logFC| threshold",
         y = "Fraction",
         title = "Data-driven logFC threshold calibration",
         subtitle = paste0("ashr shrinkage + ", n_controls_found, " positive control genes")) +
    theme_masld() +
    theme(legend.position = c(0.7, 0.8),
          plot.subtitle = element_text(size = 5.5, color = "grey40"))

  save_fig(p1, file.path(FIG_DIR, "figS_ashr_threshold_calibration.pdf"),
           width = fig_half_width, height = 3)
  cat("Calibration figure saved\n")
}

# ── 8. Figure: shrunk vs original logFC ──────────────────────────────────────
p2 <- ggplot(dream[!is.na(shrunk_logFC)], aes(x = logFC, y = shrunk_logFC)) +
  geom_hex(bins = 100) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.3, color = "grey50") +
  geom_hline(yintercept = c(-recommended, recommended), linetype = "dotted",
             linewidth = 0.3, color = "#C2185B") +
  scale_fill_gradient(low = "#E3F2FD", high = "#0D47A1", trans = "log10", name = "Count") +
  labs(x = "Original dream logFC", y = "ashr shrunk logFC",
       title = "Adaptive shrinkage of effect sizes",
       subtitle = paste0("Dotted lines: recommended threshold ±", recommended)) +
  theme_masld() +
  theme(plot.subtitle = element_text(size = 5.5, color = "grey40"))

save_fig(p2, file.path(FIG_DIR, "figS_ashr_shrinkage.pdf"),
         width = fig_half_width, height = fig_half_width)
cat("Shrinkage figure saved\n")

cat("\nDone.\n")
