#!/usr/bin/env Rscript
# ===========================================================================
# 05b_ordinal_lvqw.R
# ---------------------------------------------------------------------------
# HARMONIZED replacement for 05b_dream_ordinal.R: same ORDINAL/fibrosis
# analysis (numeric trend + per-transition successive-stage contrasts), but
# driven by the SINGLE vetted limma-voom-quality-weighted (LVQW) engine
# (de_engine_lvqw.R) instead of dream(). The random intercept (1|dataset) is
# replaced by a FIXED dataset effect; inferred_sex stays a nuisance covariate.
#
# Cohort-set principle (WIDER disease-internal set — consistent with 14)
# ---------------------------------------------------------------------
# These are disease-INTERNAL progression contrasts, so we do NOT impose
# include_in_mega. We mirror 14_fibrosis_progression_lvqw.R EXACTLY:
#   - QC gate: pass_technical == TRUE (qc/sample_qc_report.csv).
#   - Universe: EVERY cohort carrying a defined fibrosis_stage in {0..4}
#     (7 cohorts: GSE130970, GSE135251, GSE162694, GSE174478, GSE193066,
#      GSE213621, GSE240729; n = 1,077).
#   - ORDINAL trend: every cohort with >=2 distinct stages (all 7).
#   - EACH TRANSITION: every cohort containing BOTH adjacent stages of the pair
#     (build_design_guarded drops a `dataset` lacking a stage, and drops
#      `inferred_sex` if it collapses to a single level). This is why F3->F4
#      uses 6 cohorts (GSE213621 has no F4).
# This makes 05b's fibrosis ordinal consistent with 14's (same numeric trend on
# the same 7-cohort set), per PI 2026-06-08.
#
# Input loading mirrors 05h / 14_lvqw (merged_dge.rds + meta_matched.rds; subset
# only, NO re-TMM / NO re-filterByExpr — the merged DGE is already normalized).
#
# Two fits (same biology as the original 05b)
# -------------------------------------------
#   (1) Numeric trend:  ~ dataset + inferred_sex + fib_num   (coef = fib_num)
#       fib_num = as.numeric(fibrosis_stage); per-stage LFC slope. Pooled across
#       all 7 ordinal cohorts. build_design_guarded.
#   (2) Successive-stage transitions F0->F1, F1->F2, F2->F3, F3->F4: each a binary
#       low/high contrast pooled across the cohorts that contain BOTH stages.
#       Design ~ dataset + inferred_sex + fib_group ; coef = "fib_grouphigh"
#       (mirrors 14's parameterization). This REPLACES makeContrastsDream and
#       lets the per-transition cohort set differ (F3->F4 drops GSE213621).
#
# STAGING ONLY. Output:
#   results/integration/ordinal_lvqw.csv         (long, two-tier schema)
#   results/integration/ordinal_lvqw_sanity.csv  (per-term summary / dream gate)
# Backs up dream_results_ordinal.csv -> *_premigration_backup.csv if present.
# DOES NOT overwrite any canonical file.
# ===========================================================================

set.seed(42)
Sys.setenv(R_PARALLEL_SEED = "42")

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(ashr)
  library(data.table)
})

INPUT_ROOT  <- Sys.getenv("MASLD_INPUT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
OUTPUT_ROOT <- Sys.getenv("MASLD_OUTPUT_ROOT", INPUT_ROOT)
INT      <- file.path(INPUT_ROOT,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration")
IN_RDIR  <- file.path(INT, "results/integration")
OUT_RDIR <- file.path(OUTPUT_ROOT,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
SCRIPT_DIR <- file.path(INT, "scripts")
dir.create(OUT_RDIR, recursive = TRUE, showWarnings = FALSE)

source(file.path(SCRIPT_DIR, "de_engine_lvqw.R"))

cat("=== 05b_ordinal_lvqw.R (LVQW harmonization of 05b_dream_ordinal) ===\n")
cat("Started:", as.character(Sys.time()), "\n\n")

# ---------------------------------------------------------------------------
# Input loading (05h / 14_lvqw template) — subset only, NO re-norm/filter.
# Apply 14's QC gate (pass_technical) and the fibrosis-annotated universe.
# ---------------------------------------------------------------------------
dge_all <- readRDS(file.path(IN_RDIR, "merged_dge.rds"))
meta    <- as.data.table(readRDS(file.path(IN_RDIR, "meta_matched.rds")))
qc      <- fread(file.path(INT, "qc/sample_qc_report.csv"))

keep_ids <- qc[pass_technical == TRUE, sample_id]
meta <- meta[sample_id %in% keep_ids & sample_id %in% colnames(dge_all)]

# fibrosis-annotated, restricted to integer stages 0..4
meta[, fib_int := suppressWarnings(as.integer(as.character(fibrosis_stage)))]
meta_fib <- meta[!is.na(fib_int) & fib_int %in% 0:4]
cat(sprintf("[in] %d genes x %d samples in merged_dge; %d fibrosis-annotated {0..4} (pass_technical)\n",
            nrow(dge_all), ncol(dge_all), nrow(meta_fib)))
cat("Fibrosis stage distribution:\n"); print(table(meta_fib$fib_int, useNA = "always"))
cat("Cohort x stage:\n"); print(table(meta_fib$dataset, meta_fib$fib_int))

# helper: subset the (already normalized) merged DGE to a sample set and attach
# sample covariates from meta. NO calcNormFactors / NO filterByExpr re-run.
subset_dge <- function(sample_ids, cov_dt) {
  idx <- colnames(dge_all) %in% sample_ids
  d   <- dge_all[, idx]
  cov <- cov_dt[match(colnames(d), cov_dt$sample_id)]
  d$samples$dataset      <- factor(cov$dataset)
  d$samples$inferred_sex <- factor(cov$inferred_sex)
  d$samples$fib_int      <- cov$fib_int
  d
}

results_long <- list()

# ===========================================================================
# Fit 1: FIBROSIS ORDINAL SLOPE — numeric trend, dataset + sex fixed
#   ~ dataset + inferred_sex + fib_num ; coef = fib_num
#   pooled across every cohort with >=2 distinct stages (all 7 here).
# ===========================================================================
cat("\n--- Fit 1: numeric fibrosis trend (LVQW) ---\n")
ord_cohorts <- meta_fib[, .(n = uniqueN(fib_int)), by = dataset][n >= 2, dataset]
m_ord <- meta_fib[dataset %in% ord_cohorts]
cat(sprintf("Ordinal cohort set (>=2 stages): {%s}\n",
            paste(sort(unique(as.character(m_ord$dataset))), collapse = ", ")))
cat(sprintf("Ordinal samples: %d across %d cohorts\n",
            nrow(m_ord), uniqueN(m_ord$dataset)))

dge_ord  <- subset_dge(m_ord$sample_id, m_ord)
info_ord <- data.frame(
  dataset      = factor(dge_ord$samples$dataset),
  inferred_sex = factor(dge_ord$samples$inferred_sex),
  fib_num      = as.numeric(dge_ord$samples$fib_int)
)
rownames(info_ord) <- colnames(dge_ord)

bg_ord <- build_design_guarded(info_ord, c("dataset", "inferred_sex", "fib_num"))
cat("Trend formula used:", bg_ord$formula_used, "\n")
if (length(bg_ord$dropped))
  cat("  dropped:", paste(bg_ord$dropped, collapse = ", "), "\n")
stopifnot("fib_num" %in% colnames(bg_ord$design))

dt_trend <- fit_lvqw(dge_ord, bg_ord$design, coef = "fib_num", do_ashr = TRUE)
dt_trend[, transition := "linear_trend"]
dt_trend[, term       := "fib_num"]
dt_trend[, n_cohorts  := uniqueN(m_ord$dataset)]
dt_trend[, n_samples  := nrow(m_ord)]
dt_trend[, method     := "limma_voom_qw_C2"]
results_long[["linear_trend"]] <- dt_trend
cat("Linear trend DEGs (padj<0.1):", sum(dt_trend$padj < 0.1, na.rm = TRUE), "\n")

# ===========================================================================
# Fit 2: SUCCESSIVE-STAGE TRANSITIONS — binary low/high, dataset + sex fixed
#   per pair: ~ dataset + inferred_sex + fib_group ; coef = fib_grouphigh
#   pooled across cohorts containing BOTH adjacent stages (build_design_guarded
#   drops a dataset lacking a stage; F3->F4 -> 6 cohorts, GSE213621 has no F4).
# ===========================================================================
cat("\n--- Fit 2: successive-stage transitions (LVQW) ---\n")
stage_pairs <- list(
  F0vF1 = c(lo = 0L, hi = 1L),
  F1vF2 = c(lo = 1L, hi = 2L),
  F2vF3 = c(lo = 2L, hi = 3L),
  F3vF4 = c(lo = 3L, hi = 4L)
)

for (tname in names(stage_pairs)) {
  f_lo <- stage_pairs[[tname]]["lo"]; f_hi <- stage_pairs[[tname]]["hi"]

  both <- meta_fib[fib_int %in% c(f_lo, f_hi),
                   .(ns = uniqueN(fib_int)), by = dataset][ns == 2L, dataset]
  m_pair <- meta_fib[dataset %in% both & fib_int %in% c(f_lo, f_hi)]

  if (nrow(m_pair) < 6L || uniqueN(m_pair$fib_int) < 2L) {
    cat(sprintf("[skip] %s (n=%d, cohorts={%s})\n",
                tname, nrow(m_pair), paste(sort(both), collapse = ",")))
    next
  }
  cat(sprintf("%s: %d cohorts {%s}, %d samples (low=%d, high=%d)\n",
              tname, length(both), paste(sort(both), collapse = ","), nrow(m_pair),
              sum(m_pair$fib_int == f_lo), sum(m_pair$fib_int == f_hi)))

  dge_p <- subset_dge(m_pair$sample_id, m_pair)
  fib_group <- factor(ifelse(dge_p$samples$fib_int == f_hi, "high", "low"),
                      levels = c("low", "high"))
  info_p <- data.frame(
    dataset      = factor(dge_p$samples$dataset),
    inferred_sex = factor(dge_p$samples$inferred_sex),
    fib_group    = fib_group
  )
  rownames(info_p) <- colnames(dge_p)

  bg_p <- build_design_guarded(info_p, c("dataset", "inferred_sex", "fib_group"))
  if (length(bg_p$dropped))
    cat(sprintf("  design: %s ; dropped: %s\n",
                bg_p$formula_used, paste(bg_p$dropped, collapse = ", ")))
  else
    cat(sprintf("  design: %s\n", bg_p$formula_used))
  stopifnot("fib_grouphigh" %in% colnames(bg_p$design))

  dt_tr <- fit_lvqw(dge_p, bg_p$design, coef = "fib_grouphigh", do_ashr = TRUE)
  dt_tr[, transition := tname]
  dt_tr[, term       := "fib_grouphigh"]
  dt_tr[, n_cohorts  := length(both)]
  dt_tr[, n_samples  := nrow(m_pair)]
  dt_tr[, method     := "limma_voom_qw_C2"]
  results_long[[tname]] <- dt_tr
  cat("  DEGs (padj<0.1):", sum(dt_tr$padj < 0.1, na.rm = TRUE), "\n")
}

# ---------------------------------------------------------------------------
# Combine + save (STAGING ONLY)
# Schema: gene, logFC, SE, t, P.Value, padj, shrunk_logFC, lfsr, AveExpr,
#         transition, term, method, n_cohorts, n_samples
# ---------------------------------------------------------------------------
out <- rbindlist(results_long, use.names = TRUE, fill = TRUE)
setcolorder(out, c("gene", "logFC", "SE", "t", "P.Value", "padj",
                   "shrunk_logFC", "lfsr", "AveExpr",
                   "transition", "term", "method", "n_cohorts", "n_samples"))
out_csv <- file.path(OUT_RDIR, "ordinal_lvqw.csv")
fwrite(out, out_csv)
cat("\nSaved (STAGING):", out_csv, "\n")
cat("Rows:", nrow(out), "(", length(unique(out$gene)), "unique genes x",
    length(unique(out$transition)), "terms )\n")

cat("\n--- per-term DEG counts (Tier-2 padj<0.1) ---\n")
print(out[, .(n_cohorts = n_cohorts[1], n_samples = n_samples[1],
              n_tested = .N,
              n_DEG_padj0.1 = sum(padj < 0.1, na.rm = TRUE),
              n_up = sum(padj < 0.1 & logFC > 0, na.rm = TRUE),
              n_dn = sum(padj < 0.1 & logFC < 0, na.rm = TRUE)),
          by = transition])

# ---------------------------------------------------------------------------
# Backup canonical dream_results_ordinal.csv (copy, if present) — never overwrite
# ---------------------------------------------------------------------------
dream_csv  <- file.path(IN_RDIR, "dream_results_ordinal.csv")
backup_csv <- file.path(IN_RDIR, "dream_results_ordinal_premigration_backup.csv")
if (file.exists(dream_csv)) {
  if (!file.exists(backup_csv)) {
    file.copy(dream_csv, backup_csv, overwrite = FALSE)
    cat("\nBacked up canonical ->", backup_csv, "\n")
  } else {
    cat("\nBackup already exists, leaving as-is:", backup_csv, "\n")
  }
} else {
  cat("\n[note] dream_results_ordinal.csv not present on disk — nothing to back up.\n")
}

# ===========================================================================
# SANITY GATE vs dream_results_ordinal.csv (per term n DEG, direction
# concordance, logFC Spearman). If the dream reference is absent, emit a
# self-contained LVQW summary so the gate file is never silently empty.
# ===========================================================================
cat("\n=== SANITY GATE ===\n")
san_rows <- list()
if (file.exists(dream_csv)) {
  dr <- fread(dream_csv)
  for (tr in unique(out$transition)) {
    a <- out[transition == tr, .(gene, lvqw_lfc = logFC, lvqw_padj = padj)]
    b <- dr[transition == tr, .(gene, dr_lfc = logFC, dr_padj = padj)]
    m <- merge(a, b, by = "gene")
    if (!nrow(m)) next
    rho  <- suppressWarnings(cor(m$lvqw_lfc, m$dr_lfc, method = "spearman",
                                 use = "complete.obs"))
    pear <- suppressWarnings(cor(m$lvqw_lfc, m$dr_lfc, method = "pearson",
                                 use = "complete.obs"))
    dirA <- mean(sign(m$lvqw_lfc) == sign(m$dr_lfc), na.rm = TRUE)
    n_lvqw <- sum(a$lvqw_padj < 0.1, na.rm = TRUE)
    n_dr   <- sum(b$dr_padj   < 0.1, na.rm = TRUE)
    jac <- length(intersect(a[lvqw_padj < 0.1, gene], b[dr_padj < 0.1, gene])) /
           max(1L, length(union(a[lvqw_padj < 0.1, gene], b[dr_padj < 0.1, gene])))
    san_rows[[tr]] <- data.table(
      transition = tr, n_genes_matched = nrow(m),
      n_DEG_lvqw = n_lvqw, n_DEG_dream = n_dr,
      logFC_spearman = rho, logFC_pearson = pear,
      direction_agree = dirA, jaccard_DEG = jac,
      reference = "dream_results_ordinal.csv")
  }
} else {
  cat("[gate] No dream_results_ordinal.csv reference on disk — ",
      "emitting LVQW-only per-term summary.\n", sep = "")
  for (tr in unique(out$transition)) {
    a <- out[transition == tr]
    san_rows[[tr]] <- data.table(
      transition = tr, n_genes_matched = nrow(a),
      n_DEG_lvqw = sum(a$padj < 0.1, na.rm = TRUE), n_DEG_dream = NA_integer_,
      logFC_spearman = NA_real_, logFC_pearson = NA_real_,
      direction_agree = NA_real_, jaccard_DEG = NA_real_,
      reference = "NONE (dream ordinal not on disk)")
  }
}
san <- rbindlist(san_rows, use.names = TRUE, fill = TRUE)
san_csv <- file.path(OUT_RDIR, "ordinal_lvqw_sanity.csv")
fwrite(san, san_csv)
cat("Saved sanity:", san_csv, "\n")
print(san)

cat("\nFinished:", as.character(Sys.time()), "\n")
