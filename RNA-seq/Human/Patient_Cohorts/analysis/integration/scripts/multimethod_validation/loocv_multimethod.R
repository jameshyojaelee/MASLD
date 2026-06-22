#!/usr/bin/env Rscript
# loocv_multimethod.R
# ===========================================================================
# Leave-one-cohort-out (LOO-CV) fold of the multi-method DE validation harness.
#
# One fold = hold out ONE mega cohort (env HELD_OUT), refit every ACTIVE method
# on the remaining (k-1) training cohorts, then score each method's training
# result against the held-out cohort's per-study DE (= an independent cohort
# never seen during training).
#
# All methods run on the IDENTICAL training split, so per-method LOO metrics are
# directly comparable fold-by-fold (the aggregator's paired Wilcoxon relies on
# this).
#
# Methods (gated by env VALIDATION_METHODS, default "dream,deseq2"):
#   dream   — run_dream(counts_tr, meta_tr, bp)   [refit on training split]
#   deseq2  — run_deseq2(counts_tr, meta_tr, bp)  [refit on training split]
#   metafor — PRECOMPUTED per-study DE reuse (gated off by default).
#
# Run via: HELD_OUT=GSE126848 Rscript loocv_multimethod.R
# Driver:  run_loocv_multimethod.sh (array 1-5, one cohort per task)
# ===========================================================================

set.seed(42)
Sys.setenv(R_PARALLEL_SEED = "42")

t_start <- Sys.time()

suppressPackageStartupMessages({
  library(data.table)
})

# --- Shared helpers (load_mega_data, run_*, compute_loocv_metrics, ...) ------
SCRIPT_DIR <- file.path(
  Sys.getenv("MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"),
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/multimethod_validation")
source(file.path(SCRIPT_DIR, "de_validation_helpers.R"))

# ===========================================================================
# 1. Resolve the fold (held-out cohort + training set)
# ===========================================================================
HELD_OUT <- Sys.getenv("HELD_OUT", "")
COHORTS  <- mega_cohorts()
if (!nzchar(HELD_OUT))
  stop("HELD_OUT env var not set. Use: HELD_OUT=<cohort> Rscript loocv_multimethod.R")
if (!HELD_OUT %in% COHORTS)
  stop("HELD_OUT='", HELD_OUT, "' is not a mega cohort. Valid: ",
       paste(COHORTS, collapse = ", "))

training <- setdiff(COHORTS, HELD_OUT)

cat("===== LOO-CV multi-method fold =====\n")
cat("Started:        ", format(t_start), "\n")
cat("Held-out cohort:", HELD_OUT, "\n")
cat("Training set:   ", paste(training, collapse = ", "),
    sprintf("(k = %d)\n", length(training)))
cat("Active methods: ", paste(ACTIVE_METHODS(), collapse = ", "), "\n\n")

# ===========================================================================
# 2. Load mega data, subset to training cohorts
# ===========================================================================
d         <- load_mega_data()
train_idx <- d$meta$dataset %in% training
counts_tr <- d$counts[, train_idx, drop = FALSE]
meta_tr   <- d$meta[train_idx]
stopifnot(all(meta_tr$sample_id == colnames(counts_tr)))

n_samples_training <- ncol(counts_tr)
cat(sprintf("Training data: %d genes x %d samples\n",
            nrow(counts_tr), n_samples_training))
cat("  per-cohort n:\n")
print(table(meta_tr$dataset))
cat("\n")

# ===========================================================================
# 3. Run each active method on the training split
# ===========================================================================
bp <- bp_param()
method_results <- list()

for (m in ACTIVE_METHODS()) {
  cat(sprintf("── Running method '%s' on training split ──\n", m))
  res <- tryCatch({
    if (m == "dream") {
      run_dream(counts_tr, meta_tr, bp)
    } else if (m == "deseq2") {
      run_deseq2(counts_tr, meta_tr, bp)
    } else if (m == "metafor") {
      # RECOMPUTE per-study DE in BINARY mode on the training cohorts (NOT the
      # on-disk per_study/*.csv, which are the yaml averaged-substage estimand).
      # The metafor investigation (M3) requires the binary Disease-vs-Control
      # contrast so metafor targets the SAME estimand as dream/DESeq2. Holding a
      # cohort out doesn't change the others' samples, but the on-disk CSVs are
      # the WRONG (config) estimand — so we refit binary here.
      per_study_list <- lapply(training, function(ds) {
        sel <- meta_tr$dataset == ds
        run_per_study_voom(ds, counts_tr[, sel, drop = FALSE],
                           meta_tr[sel], contrast_mode = "binary")
      })
      run_metafor(per_study_list, K = length(training), bp)
    } else if (m == "limma_voom") {
      run_limma_voom(counts_tr, meta_tr, bp)
    } else if (m == "limma_voom_qw") {
      run_limma_voom_qw(counts_tr, meta_tr, bp)
    } else if (m == "limma_trend") {
      run_limma_trend(counts_tr, meta_tr, bp)
    } else if (m == "edger_qlf") {
      run_edger_qlf(counts_tr, meta_tr, bp)
    } else if (m == "edger_qlf_robust") {
      run_edger_qlf_robust(counts_tr, meta_tr, bp)
    } else if (m == "edger_lrt") {
      run_edger_lrt(counts_tr, meta_tr, bp)
    } else {
      stop("Unknown method '", m, "' in VALIDATION_METHODS")
    }
  }, error = function(e) {
    warning(sprintf("Method '%s' FAILED on fold %s: %s", m, HELD_OUT,
                    conditionMessage(e)))
    NULL
  })
  if (is.null(res)) {
    cat(sprintf("  '%s' -> NULL (failed); skipping in metrics.\n", m))
  } else {
    cat(sprintf("  '%s' -> %d genes (%d with padj<0.05)\n",
                m, nrow(res), sum(!is.na(res$padj) & res$padj < 0.05)))
  }
  method_results[[m]] <- res
}

# ===========================================================================
# 4. Held-out target: the held-out cohort's per-study DE
# ===========================================================================
ho_file <- file.path(PERSTUDY_DIR, paste0(HELD_OUT, "_de_results.csv"))
if (!file.exists(ho_file))
  stop("Held-out per-study DE not found: ", ho_file)
ho    <- fread(ho_file)
ho_dt <- ho[, .(gene, ho_logFC = logFC, ho_padj = adj.P.Val, ho_pval = P.Value)]
cat(sprintf("\nHeld-out target (%s): %d genes (%d with padj<0.05)\n",
            HELD_OUT, nrow(ho_dt),
            sum(!is.na(ho_dt$ho_padj) & ho_dt$ho_padj < 0.05)))

# ===========================================================================
# 5. Per-method LOO metrics (training vs held-out)
# ===========================================================================
OUT_DIR <- file.path(RDIR, "multimethod_validation", "loocv")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Method tag so per-(method x cohort) array tasks write distinct files and don't
# clobber one another. Default (legacy dream,deseq2[,metafor]) => no tag, so the
# original loocv_<cohort>.csv naming is preserved. A single-method subset (the
# new per-method array) => loocv_<cohort>__<method>.csv. The aggregator globs
# loocv_<cohort>*.csv and rbinds all method rows.
.active <- ACTIVE_METHODS()
METHOD_TAG <- if (identical(sort(.active), sort(c("dream", "deseq2"))) ||
                  identical(sort(.active), sort(c("dream", "deseq2", "metafor"))))
                "" else paste0("__", paste(.active, collapse = "_"))

metric_rows <- list()
for (m in ACTIVE_METHODS()) {
  if (is.null(method_results[[m]])) next
  row <- compute_loocv_metrics(method_results[[m]], ho_dt, m)
  row[, `:=`(held_out_cohort    = HELD_OUT,
             n_training_cohorts = length(training),
             n_samples_training = n_samples_training)]
  metric_rows[[m]] <- row
}
if (length(metric_rows) == 0)
  stop("No method produced metrics for fold ", HELD_OUT, " (all methods failed).")

metrics <- rbindlist(metric_rows, fill = TRUE)
metrics_file <- file.path(OUT_DIR, paste0("loocv_", HELD_OUT, METHOD_TAG, ".csv"))
fwrite(metrics, metrics_file)
cat("\nSaved per-fold metrics:", metrics_file, "\n")

# ===========================================================================
# 6. Wide per-gene file: every method's (logFC, padj) + held-out, full outer
# ===========================================================================
pergene <- NULL
for (m in ACTIVE_METHODS()) {
  res <- method_results[[m]]
  if (is.null(res)) next
  cols <- res[, .(gene, logFC, padj)]
  setnames(cols, c("logFC", "padj"), c(paste0(m, "_logFC"), paste0(m, "_padj")))
  pergene <- if (is.null(pergene)) cols else merge(pergene, cols, by = "gene", all = TRUE)
}
# Held-out columns (full outer so the held-out universe is fully represented)
if (is.null(pergene)) {
  pergene <- ho_dt[, .(gene, ho_logFC, ho_padj)]
} else {
  pergene <- merge(pergene, ho_dt[, .(gene, ho_logFC, ho_padj)],
                   by = "gene", all = TRUE)
}
pergene_file <- file.path(OUT_DIR, paste0("loocv_pergene_", HELD_OUT, METHOD_TAG, ".csv"))
fwrite(pergene, pergene_file)
cat("Saved per-gene wide table:", pergene_file,
    sprintf("(%d genes)\n", nrow(pergene)))

# ===========================================================================
# 7. Summary table + timing
# ===========================================================================
cat("\n===== Fold summary:", HELD_OUT, "=====\n")
print(metrics[, .(method, n_genes_common, n_degs_train_005, n_degs_ho_005,
                  jaccard_01, lfc_spearman, direction_concordance,
                  auc_replication, fisher_or, pi1_heldout)], digits = 4)

t_end <- Sys.time()
cat("\nFinished:", format(t_end), "\n")
cat("Elapsed: ", round(as.numeric(difftime(t_end, t_start, units = "mins")), 1),
    "min\n")
