#!/usr/bin/env Rscript
# 05b2_ashr_normal_sensitivity.R
# ---------------------------------------------------------------------------
# Sensitivity analysis: compare ashr mixcompdist = "halfuniform" (canonical)
# vs mixcompdist = "normal" to assess whether the prior choice materially
# affects shrinkage estimates. [R2 #8]
#
# Design:
#   1. Load dream_results.csv (canonical kallisto dream output)
#   2. Re-run ashr with mixcompdist = "normal"
#   3. Load existing ashr results (mixcompdist = "halfuniform") from dream_results_ashr.csv
#   4. Compare: Spearman rho of shrunk_logFC between halfuniform and normal
#   5. Compare: overlap of lfsr < 0.05 genes (Jaccard index)
#   6. If Spearman > 0.99 -> choice is immaterial
#
# Input:  results/integration/dream_results.csv
#         results/integration/dream_results_ashr.csv (halfuniform baseline)
# Output: RNA-seq/results/audit_sensitivity/ashr_mixcomp_sensitivity/
#           ashr_normal_results.csv
#           ashr_mixcomp_comparison.csv
#           REPORT.md
#
# SLURM: cpu, 4 CPU, 32G, 2h, --job-name=ashr-sensitivity
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ashr)
})

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
INT  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR <- file.path(INT, "results/integration")
ODIR <- file.path(BASE, "RNA-seq/results/audit_sensitivity/ashr_mixcomp_sensitivity")
dir.create(ODIR, showWarnings = FALSE, recursive = TRUE)

cat("============================================================\n")
cat("05b2: ashr mixcompdist sensitivity (halfuniform vs normal)\n")
cat("============================================================\n\n")

# ── 1. Load dream results ────────────────────────────────────────────────────
dream <- fread(file.path(RDIR, "dream_results.csv"))
cat("Loaded:", nrow(dream), "genes from dream_results.csv\n")

# ── 2. Compute SE ────────────────────────────────────────────────────────────
if ("SE" %in% names(dream)) {
  dream[, se := SE]
  cat("Using pre-computed SE column from dream output\n")
} else {
  dream[, se := abs(logFC / t)]
  cat("SE column not found; falling back to |logFC / t|\n")
}

# Handle edge cases
valid <- !is.na(dream$se) & is.finite(dream$se) & dream$se > 0
n_valid <- sum(valid)
n_bad <- sum(!valid)
cat(sprintf("Valid genes for ashr: %d (excluded: %d)\n", n_valid, n_bad))

# ── 3. Run ashr with mixcompdist = "normal" ──────────────────────────────────
cat("\nRunning ashr with mixcompdist = 'normal' on", n_valid, "genes...\n")
t0 <- Sys.time()

ash_normal <- ash(
  betahat     = dream$logFC[valid],
  sebetahat   = dream$se[valid],
  mixcompdist = "normal",
  method      = "shrink"
)

elapsed <- round(difftime(Sys.time(), t0, units = "secs"), 1)
cat(sprintf("ashr (normal) completed in %.1f seconds\n", elapsed))

# Store normal results
dream[, normal_shrunk_logFC := NA_real_]
dream[, normal_lfsr         := NA_real_]
dream[, normal_shrunk_se    := NA_real_]

dream[valid, normal_shrunk_logFC := ash_normal$result$PosteriorMean]
dream[valid, normal_lfsr         := ash_normal$result$lfsr]
dream[valid, normal_shrunk_se    := ash_normal$result$PosteriorSD]

# ── 4. Load canonical halfuniform results ────────────────────────────────────
ashr_file <- file.path(RDIR, "dream_results_ashr.csv")
if (!file.exists(ashr_file)) {
  stop("dream_results_ashr.csv not found — run 05b_ashr_shrinkage.R first")
}

ashr_hu <- fread(ashr_file)
cat("Loaded:", nrow(ashr_hu), "genes from dream_results_ashr.csv (halfuniform)\n")

# Match genes
if ("gene" %in% names(dream)) {
  key_col <- "gene"
} else {
  stop("Cannot find gene identifier column in dream_results.csv")
}

# Merge halfuniform results onto dream
ashr_hu_sub <- ashr_hu[, .(gene, hu_shrunk_logFC = shrunk_logFC, hu_lfsr = lfsr)]
merged <- merge(dream, ashr_hu_sub, by = "gene", all.x = TRUE)

# ── 5. Compute concordance metrics ──────────────────────────────────────────
cat("\n--- Concordance metrics ---\n")

# Restrict to genes with valid results in BOTH
both_valid <- !is.na(merged$hu_shrunk_logFC) & !is.na(merged$normal_shrunk_logFC)
n_both <- sum(both_valid)
cat(sprintf("Genes with valid results in both: %d\n", n_both))

# 5a. Spearman correlation of shrunk_logFC
rho_shrunk <- cor(merged$hu_shrunk_logFC[both_valid],
                  merged$normal_shrunk_logFC[both_valid],
                  method = "spearman", use = "complete.obs")
pearson_shrunk <- cor(merged$hu_shrunk_logFC[both_valid],
                      merged$normal_shrunk_logFC[both_valid],
                      method = "pearson", use = "complete.obs")
cat(sprintf("  Spearman rho (shrunk_logFC): %.6f\n", rho_shrunk))
cat(sprintf("  Pearson r (shrunk_logFC):    %.6f\n", pearson_shrunk))

# 5b. Direction concordance
dir_conc <- mean(sign(merged$hu_shrunk_logFC[both_valid]) ==
                 sign(merged$normal_shrunk_logFC[both_valid]), na.rm = TRUE)
cat(sprintf("  Direction concordance:       %.4f\n", dir_conc))

# 5c. Mean absolute difference in shrunk_logFC
mad_shrunk <- mean(abs(merged$hu_shrunk_logFC[both_valid] -
                       merged$normal_shrunk_logFC[both_valid]), na.rm = TRUE)
max_diff <- max(abs(merged$hu_shrunk_logFC[both_valid] -
                    merged$normal_shrunk_logFC[both_valid]), na.rm = TRUE)
cat(sprintf("  Mean |delta shrunk_logFC|:   %.6f\n", mad_shrunk))
cat(sprintf("  Max  |delta shrunk_logFC|:   %.6f\n", max_diff))

# 5d. lfsr overlap
hu_sig  <- merged$hu_lfsr < 0.05 & !is.na(merged$hu_lfsr)
nor_sig <- merged$normal_lfsr < 0.05 & !is.na(merged$normal_lfsr)
n_hu_sig  <- sum(hu_sig)
n_nor_sig <- sum(nor_sig)
n_intersect <- sum(hu_sig & nor_sig)
n_union     <- sum(hu_sig | nor_sig)
jaccard_lfsr <- ifelse(n_union > 0, n_intersect / n_union, NA)
cat(sprintf("\n  lfsr < 0.05 (halfuniform):   %d\n", n_hu_sig))
cat(sprintf("  lfsr < 0.05 (normal):        %d\n", n_nor_sig))
cat(sprintf("  Intersection:                %d\n", n_intersect))
cat(sprintf("  Union:                       %d\n", n_union))
cat(sprintf("  Jaccard (lfsr < 0.05):       %.4f\n", jaccard_lfsr))

# 5e. lfsr correlation
rho_lfsr <- cor(merged$hu_lfsr[both_valid],
                merged$normal_lfsr[both_valid],
                method = "spearman", use = "complete.obs")
cat(sprintf("  Spearman rho (lfsr):         %.6f\n", rho_lfsr))

# 5f. Also compare at the stringent lfsr < 0.05 + |shrunk_logFC| > 0.2 threshold
hu_strict  <- hu_sig & abs(merged$hu_shrunk_logFC) > 0.2
nor_strict <- nor_sig & abs(merged$normal_shrunk_logFC) > 0.2
n_hu_strict  <- sum(hu_strict, na.rm = TRUE)
n_nor_strict <- sum(nor_strict, na.rm = TRUE)
n_int_strict <- sum(hu_strict & nor_strict, na.rm = TRUE)
n_uni_strict <- sum(hu_strict | nor_strict, na.rm = TRUE)
jaccard_strict <- ifelse(n_uni_strict > 0, n_int_strict / n_uni_strict, NA)
cat(sprintf("\n  lfsr<0.05 + |shrunk|>0.2 (halfuniform): %d\n", n_hu_strict))
cat(sprintf("  lfsr<0.05 + |shrunk|>0.2 (normal):      %d\n", n_nor_strict))
cat(sprintf("  Jaccard (strict):                        %.4f\n", jaccard_strict))

# ── 6. Verdict ───────────────────────────────────────────────────────────────
cat("\n============================================================\n")
if (rho_shrunk > 0.99) {
  verdict <- "PASS"
  verdict_text <- sprintf(
    "PASS: Spearman rho = %.4f (> 0.99). mixcompdist choice is immaterial.",
    rho_shrunk)
} else if (rho_shrunk > 0.95) {
  verdict <- "PARTIAL"
  verdict_text <- sprintf(
    "PARTIAL: Spearman rho = %.4f (0.95-0.99). Minor sensitivity; report in supplement.",
    rho_shrunk)
} else {
  verdict <- "FAIL"
  verdict_text <- sprintf(
    "FAIL: Spearman rho = %.4f (< 0.95). mixcompdist choice materially affects results.",
    rho_shrunk)
}
cat(verdict_text, "\n")
cat("============================================================\n")

# ── 7. Save outputs ──────────────────────────────────────────────────────────

# 7a. Per-gene comparison table
comp <- merged[both_valid, .(
  gene,
  raw_logFC = logFC,
  hu_shrunk_logFC,
  normal_shrunk_logFC,
  delta_shrunk = normal_shrunk_logFC - hu_shrunk_logFC,
  hu_lfsr,
  normal_lfsr,
  delta_lfsr = normal_lfsr - hu_lfsr
)]
fwrite(comp, file.path(ODIR, "ashr_mixcomp_comparison.csv"))
cat("\nSaved: ashr_mixcomp_comparison.csv (", nrow(comp), "genes)\n")

# 7b. Normal-prior results
out_norm <- dream[, .(gene, logFC, normal_shrunk_logFC, normal_lfsr, normal_shrunk_se)]
fwrite(out_norm, file.path(ODIR, "ashr_normal_results.csv"))
cat("Saved: ashr_normal_results.csv\n")

# 7c. Summary metrics
metrics <- data.table(
  metric = c("n_genes_tested",
             "spearman_rho_shrunk_logFC",
             "pearson_r_shrunk_logFC",
             "direction_concordance",
             "mean_abs_delta_shrunk_logFC",
             "max_abs_delta_shrunk_logFC",
             "n_lfsr005_halfuniform",
             "n_lfsr005_normal",
             "jaccard_lfsr005",
             "jaccard_strict_lfsr005_lfc02",
             "spearman_rho_lfsr",
             "verdict"),
  value = c(as.character(n_both),
            sprintf("%.6f", rho_shrunk),
            sprintf("%.6f", pearson_shrunk),
            sprintf("%.4f", dir_conc),
            sprintf("%.6f", mad_shrunk),
            sprintf("%.6f", max_diff),
            as.character(n_hu_sig),
            as.character(n_nor_sig),
            sprintf("%.4f", jaccard_lfsr),
            sprintf("%.4f", jaccard_strict),
            sprintf("%.6f", rho_lfsr),
            verdict)
)
fwrite(metrics, file.path(ODIR, "ashr_mixcomp_metrics.csv"))
cat("Saved: ashr_mixcomp_metrics.csv\n")

# 7d. REPORT.md
report_lines <- c(
  "# ashr mixcompdist Sensitivity [R2 #8]",
  "",
  sprintf("**Date:** %s", Sys.time()),
  "",
  "## Design",
  "",
  "- **Question**: Does the ashr prior distribution choice (`mixcompdist`) materially",
  "  affect shrinkage estimates?",
  "- **Canonical**: `mixcompdist = \"halfuniform\"` (Script 05b)",
  "- **Sensitivity**: `mixcompdist = \"normal\"`",
  sprintf("- **N genes tested**: %d", n_both),
  "",
  "## Results",
  "",
  "### Shrunk logFC concordance",
  "",
  "| Metric | Value |",
  "|--------|-------|",
  sprintf("| Spearman rho | %.6f |", rho_shrunk),
  sprintf("| Pearson r | %.6f |", pearson_shrunk),
  sprintf("| Direction concordance | %.4f |", dir_conc),
  sprintf("| Mean |delta shrunk_logFC| | %.6f |", mad_shrunk),
  sprintf("| Max |delta shrunk_logFC| | %.6f |", max_diff),
  "",
  "### lfsr < 0.05 gene overlap",
  "",
  "| Metric | Value |",
  "|--------|-------|",
  sprintf("| N sig (halfuniform) | %d |", n_hu_sig),
  sprintf("| N sig (normal) | %d |", n_nor_sig),
  sprintf("| Intersection | %d |", n_intersect),
  sprintf("| Union | %d |", n_union),
  sprintf("| Jaccard | %.4f |", jaccard_lfsr),
  sprintf("| Spearman rho (lfsr) | %.6f |", rho_lfsr),
  "",
  "### Strict threshold (lfsr < 0.05 + |shrunk_logFC| > 0.2)",
  "",
  "| Metric | Value |",
  "|--------|-------|",
  sprintf("| N sig (halfuniform) | %d |", n_hu_strict),
  sprintf("| N sig (normal) | %d |", n_nor_strict),
  sprintf("| Jaccard (strict) | %.4f |", jaccard_strict),
  "",
  "## Verdict",
  "",
  sprintf("**%s**", verdict_text),
  "",
  "## Files",
  "",
  "- `ashr_normal_results.csv` -- per-gene normal-prior ashr output",
  "- `ashr_mixcomp_comparison.csv` -- per-gene halfuniform vs normal comparison",
  "- `ashr_mixcomp_metrics.csv` -- summary metrics",
  ""
)
writeLines(report_lines, file.path(ODIR, "REPORT.md"))
cat("Saved: REPORT.md\n")

cat("\nScript 05b2 complete.\n")
