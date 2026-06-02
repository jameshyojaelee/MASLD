#!/usr/bin/env Rscript
# M03b_per_diet_classification.R
# ---------------------------------------------------------------------------
# Per-diet heterogeneity classification using tau-squared from metafor rma().
#
# Motivation: Reviewers R2 #11 and R3 P1-2 flagged that pooling 5 radically
# different diet models (MCD, HFD, CDAHFD, FPC, LIDPAD) as exchangeable is
# scientifically invalid. This script quantifies between-diet heterogeneity
# per gene and classifies genes as diet-conserved, intermediate, or
# diet-specific based on tau-squared quantiles.
#
# Input:  per_diet/*_de_results.csv (from M02)
# Output: results/mouse_per_diet_classification.csv
# ---------------------------------------------------------------------------

# ---- Seed pinning (T2.4, 2026-04-22) -----
set.seed(42)
Sys.setenv(R_PARALLEL_SEED = "42")

suppressPackageStartupMessages({
  library(data.table)
  library(metafor)
})

# ---- Paths ----
MOUSE  <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Mouse"
INT    <- file.path(MOUSE, "Unified_Integration")
RDIR   <- file.path(INT, "results")
DEDIR  <- file.path(RDIR, "per_diet")
OUTDIR <- RDIR

cat("=== M03b: Per-Diet Tau-Squared Classification ===\n")
cat("Start:", format(Sys.time()), "\n\n")

# ============================================================
# 1. Load per-diet DE results (from M02)
# ============================================================
# Auto-discover diet models from M02 output files in per_diet/
de_files <- list.files(DEDIR, pattern = "_de_results\\.csv$", full.names = TRUE)
if (length(de_files) == 0L) {
  stop("No *_de_results.csv files found in ", DEDIR,
       "\n  Run M02_mouse_per_diet_de.R first.")
}
diet_models <- sub("_de_results\\.csv$", "", basename(de_files))
cat("Auto-discovered", length(diet_models), "diet models:",
    paste(diet_models, collapse = ", "), "\n")

de_list <- list()
for (dm in diet_models) {
  f <- file.path(DEDIR, paste0(dm, "_de_results.csv"))
  de_list[[dm]] <- fread(f)
  cat("  Loaded:", dm, "-", nrow(de_list[[dm]]), "genes\n")
}

# Load the summary to annotate sample sizes
de_summary <- fread(file.path(DEDIR, "de_summary.csv"))
cat("\nPer-diet sample sizes:\n")
print(de_summary[, .(diet_model, n_samples, n_disease, n_control)])
cat("\n")

# ============================================================
# 2. Find genes present in >= 2 diet models
#    (minimum for rma; M03 uses >= 3, we keep >= 2 here to
#    also classify genes testable in exactly 2 diets)
# ============================================================
all_unique_genes <- unique(unlist(lapply(de_list, function(x) x$gene)))
gene_diet_counts <- sapply(all_unique_genes, function(g) {
  sum(vapply(de_list, function(x) g %in% x$gene, logical(1)))
})

min_diets_required <- 2L
testable_genes <- all_unique_genes[gene_diet_counts >= min_diets_required]

cat("Gene coverage across diet models:\n")
cat("  Total unique genes:", length(all_unique_genes), "\n")
cat("  Genes in >=", min_diets_required, "models:", length(testable_genes),
    sprintf(" (%.1f%%)\n", 100 * length(testable_genes) / length(all_unique_genes)))
cat("  Distribution:\n")
print(table(gene_diet_counts))
cat("\n")

# ============================================================
# 3. Run rma() per gene, extracting tau-squared + Q-test
# ============================================================
cat("Running rma() for", length(testable_genes), "genes to extract tau-squared...\n")

meta_results <- rbindlist(lapply(testable_genes, function(g) {
  effects <- vapply(de_list, function(x) {
    v <- x[gene == g, logFC]
    if (length(v) == 0L) NA_real_ else v[1L]
  }, numeric(1))
  ses <- vapply(de_list, function(x) {
    v <- x[gene == g, SE_unmoderated]
    if (length(v) == 0L) NA_real_ else v[1L]
  }, numeric(1))

  valid <- !is.na(effects) & !is.na(ses) & is.finite(ses) & ses > 0
  if (sum(valid) < 2L) return(NULL)

  tryCatch({
    fit <- rma(yi = effects[valid], sei = ses[valid], method = "REML")

    # Per-diet logFC for the driving-diet annotation
    lfc_vec <- effects[valid]
    names(lfc_vec) <- names(effects)[valid]

    # Identify direction-consistent diets (same sign as pooled estimate)
    pooled_sign <- sign(as.numeric(fit$b))
    same_dir <- names(lfc_vec)[sign(lfc_vec) == pooled_sign]
    opp_dir  <- names(lfc_vec)[sign(lfc_vec) != pooled_sign & lfc_vec != 0]

    # Per-diet significance (from original per-diet padj)
    per_diet_sig <- vapply(names(lfc_vec), function(d) {
      pv <- de_list[[d]][gene == g, adj.P.Val]
      if (length(pv) == 0L) NA_real_ else pv[1L]
    }, numeric(1))
    sig_diets <- names(per_diet_sig)[!is.na(per_diet_sig) & per_diet_sig < 0.05]

    # Strongest driving diet: highest |logFC| among direction-consistent
    if (length(same_dir) > 0) {
      driver <- same_dir[which.max(abs(lfc_vec[same_dir]))]
    } else {
      driver <- names(which.max(abs(lfc_vec)))
    }

    # Build per-diet logFC string: "MCD=1.23;HFD=-0.45;..."
    lfc_str <- paste(paste0(names(lfc_vec), "=", round(lfc_vec, 4)),
                     collapse = ";")
    padj_str <- paste(paste0(names(per_diet_sig), "=",
                             ifelse(is.na(per_diet_sig), "NA",
                                    formatC(per_diet_sig, format = "e", digits = 2))),
                      collapse = ";")

    data.table(
      gene           = g,
      meta_logFC     = as.numeric(fit$b),
      meta_se        = as.numeric(fit$se),
      meta_pval      = as.numeric(fit$pval),
      tau2           = as.numeric(fit$tau2),
      tau2_se        = as.numeric(fit$se.tau2),
      I2             = as.numeric(fit$I2),
      H2             = as.numeric(fit$H2),
      QE             = as.numeric(fit$QE),
      QEp            = as.numeric(fit$QEp),
      n_diets        = sum(valid),
      n_diets_sig    = length(sig_diets),
      n_diets_same_dir = length(same_dir),
      n_diets_opp_dir  = length(opp_dir),
      strongest_driver = driver,
      strongest_driver_lfc = as.numeric(lfc_vec[driver]),
      sig_diets      = paste(sig_diets, collapse = ";"),
      same_dir_diets = paste(same_dir, collapse = ";"),
      opp_dir_diets  = paste(opp_dir, collapse = ";"),
      per_diet_lfc   = lfc_str,
      per_diet_padj  = padj_str
    )
  }, error = function(e) NULL)
}))

cat("  rma() succeeded for", nrow(meta_results), "/", length(testable_genes), "genes\n")
if (nrow(meta_results) == 0L) stop("No genes passed rma() -- check per-diet inputs.")

# BH-adjust meta p-values
meta_results[, meta_padj := p.adjust(meta_pval, method = "BH")]

# ============================================================
# 4. Classify genes by tau-squared quantile
# ============================================================
tau2_vals <- meta_results$tau2

# Compute quantile thresholds on tau2 (excluding tau2 == 0 for a more
# informative distribution, but classification uses all genes)
cat("\nTau-squared summary (all genes):\n")
print(summary(tau2_vals))

q50 <- quantile(tau2_vals, 0.50, na.rm = TRUE)
q75 <- quantile(tau2_vals, 0.75, na.rm = TRUE)

cat("\nClassification thresholds:\n")
cat("  Median (50th pctl) tau2:", round(q50, 6), "\n")
cat("  75th percentile tau2:  ", round(q75, 6), "\n\n")

meta_results[, diet_class := fcase(
  tau2 < q50,                       "diet_conserved",
  tau2 >= q50 & tau2 <= q75,        "intermediate",
  tau2 > q75,                       "diet_specific"
)]

# Additionally flag by Cochran's Q test (formal heterogeneity test)
meta_results[, cochran_q_sig := QEp < 0.05]

# ============================================================
# 5. For diet-specific genes: identify driving diets
# ============================================================
# Already annotated above (strongest_driver, sig_diets, per_diet_lfc).
# Now compute a diet_driver_summary: which diet is named as strongest_driver
# most often among diet-specific genes.

diet_spec <- meta_results[diet_class == "diet_specific"]

cat("=== CLASSIFICATION DISTRIBUTION ===\n\n")
class_tbl <- meta_results[, .N, by = diet_class][order(-N)]
class_tbl[, pct := round(100 * N / sum(N), 1)]
print(class_tbl)
cat("\n")

cat("Cochran's Q significant (p<0.05):", sum(meta_results$cochran_q_sig, na.rm = TRUE),
    "/", nrow(meta_results), sprintf("(%.1f%%)\n\n",
    100 * sum(meta_results$cochran_q_sig, na.rm = TRUE) / nrow(meta_results)))

# Direction consistency
cat("Direction consistency:\n")
cat("  All diets same direction:", sum(meta_results$n_diets_opp_dir == 0),
    sprintf("(%.1f%%)\n", 100 * sum(meta_results$n_diets_opp_dir == 0) / nrow(meta_results)))
cat("  At least 1 diet opposite:", sum(meta_results$n_diets_opp_dir > 0),
    sprintf("(%.1f%%)\n\n", 100 * sum(meta_results$n_diets_opp_dir > 0) / nrow(meta_results)))

# Per-diet driver breakdown among diet-specific genes
cat("=== DIET-SPECIFIC GENES: DRIVING DIETS ===\n\n")
if (nrow(diet_spec) > 0) {
  driver_tbl <- diet_spec[, .N, by = strongest_driver][order(-N)]
  driver_tbl[, pct := round(100 * N / sum(N), 1)]
  print(driver_tbl)
  cat("\n")

  # Also count how many diets are significant per gene
  cat("Number of diets significant (padj<0.05) among diet-specific genes:\n")
  print(diet_spec[, .N, by = n_diets_sig][order(n_diets_sig)])
  cat("\n")
}

# ============================================================
# 6. Cross-tabulate: diet_class vs meta-analysis significance
# ============================================================
cat("=== CROSS-TABULATION: diet_class x meta-significance ===\n\n")
meta_results[, meta_sig := meta_padj < 0.05]
xtab <- meta_results[, .N, by = .(diet_class, meta_sig)]
xtab_wide <- dcast(xtab, diet_class ~ meta_sig, value.var = "N", fill = 0)
setnames(xtab_wide, c("diet_class", "not_sig", "sig"))
xtab_wide[, total := not_sig + sig]
xtab_wide[, pct_sig := round(100 * sig / total, 1)]
print(xtab_wide)
cat("\n")

# ============================================================
# 7. Per-diet contribution summary
# ============================================================
cat("=== PER-DIET CONTRIBUTION SUMMARY ===\n")
cat("(Number of genes where each diet is significant at padj<0.05)\n\n")

for (dm in diet_models) {
  n_tested <- sum(vapply(seq_len(nrow(meta_results)), function(i) {
    grepl(dm, meta_results$per_diet_lfc[i])
  }, logical(1)))
  n_sig <- sum(vapply(seq_len(nrow(meta_results)), function(i) {
    dm %in% strsplit(meta_results$sig_diets[i], ";")[[1]]
  }, logical(1)))
  cat(sprintf("  %-8s: %5d tested, %5d sig (%.1f%%)\n",
              dm, n_tested, n_sig, 100 * n_sig / max(n_tested, 1)))
}
cat("\n")

# ============================================================
# 8. Top diet-conserved and diet-specific genes
# ============================================================
cat("=== TOP 20 DIET-CONSERVED GENES (lowest tau2, meta-sig) ===\n")
top_conserved <- meta_results[diet_class == "diet_conserved" & meta_padj < 0.05
                             ][order(tau2)][1:min(20, .N)]
if (nrow(top_conserved) > 0) {
  print(top_conserved[, .(gene, meta_logFC = round(meta_logFC, 3),
                          meta_padj = formatC(meta_padj, format = "e", digits = 2),
                          tau2 = round(tau2, 6), I2 = round(I2, 1),
                          n_diets, n_diets_sig)])
}
cat("\n")

cat("=== TOP 20 DIET-SPECIFIC GENES (highest tau2, meta-sig) ===\n")
top_specific <- meta_results[diet_class == "diet_specific" & meta_padj < 0.05
                            ][order(-tau2)][1:min(20, .N)]
if (nrow(top_specific) > 0) {
  print(top_specific[, .(gene, meta_logFC = round(meta_logFC, 3),
                         meta_padj = formatC(meta_padj, format = "e", digits = 2),
                         tau2 = round(tau2, 4), I2 = round(I2, 1),
                         n_diets, strongest_driver)])
}
cat("\n")

# ============================================================
# 9. Save output
# ============================================================
meta_results <- meta_results[order(meta_padj)]

out_file <- file.path(OUTDIR, "mouse_per_diet_classification.csv")
fwrite(meta_results, out_file)
cat("Saved:", out_file, "\n")
cat("  Rows:", nrow(meta_results), " Cols:", ncol(meta_results), "\n")

cat("\n=== M03b complete:", format(Sys.time()), "===\n")
