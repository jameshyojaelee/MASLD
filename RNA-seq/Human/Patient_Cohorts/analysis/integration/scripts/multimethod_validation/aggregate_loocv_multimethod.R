#!/usr/bin/env Rscript
# aggregate_loocv_multimethod.R
# ===========================================================================
# Aggregate the per-fold multi-method LOO-CV results (loocv_<cohort>.csv, one
# per held-out mega cohort) into:
#
#   loocv/loocv_perfold.csv     — raw rbind of all folds (method x fold rows)
#   loocv/loocv_summary.csv     — per-method mean +/- sd of each metric across
#                                 folds (one row per method)
#   loocv/method_comparison.csv — head-to-head: per metric, each method's mean
#                                 plus a PAIRED Wilcoxon p (dream-vs-deseq2, and
#                                 dream-vs-metafor when metafor is present).
#                                 Paired because all methods ran on identical
#                                 folds (same held-out cohort -> same split).
#
# Run AFTER all 5 folds complete (chained via --dependency=afterok).
# Driver: run_loocv_multimethod.sh
# ===========================================================================

suppressPackageStartupMessages({
  library(data.table)
})

PROJECT_ROOT <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
RDIR    <- file.path(PROJECT_ROOT,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
LOO_DIR <- file.path(RDIR, "multimethod_validation", "loocv")

# ===========================================================================
# 1. Collect per-fold metrics files (loocv_<cohort>.csv, NOT pergene / summary)
# ===========================================================================
all_files <- list.files(LOO_DIR, pattern = "^loocv_.*\\.csv$", full.names = TRUE)
fold_files <- grep("loocv_pergene_|loocv_summary|loocv_perfold|method_comparison",
                   all_files, value = TRUE, invert = TRUE)
if (length(fold_files) == 0)
  stop("No per-fold LOO-CV files found in ", LOO_DIR)

cat("===== Aggregating multi-method LOO-CV =====\n")
cat(sprintf("Found %d per-fold file(s):\n", length(fold_files)))
for (f in fold_files) cat("  ", basename(f), "\n")

perfold <- rbindlist(lapply(fold_files, fread), fill = TRUE)
fwrite(perfold, file.path(LOO_DIR, "loocv_perfold.csv"))
cat("\nSaved raw per-fold rbind:", file.path(LOO_DIR, "loocv_perfold.csv"),
    sprintf("(%d rows)\n", nrow(perfold)))

methods <- sort(unique(perfold$method))
cat("Methods present:", paste(methods, collapse = ", "), "\n")
cat("Folds present:  ", paste(sort(unique(perfold$held_out_cohort)),
                              collapse = ", "), "\n")

# Numeric metric columns (exclude identifiers)
id_cols     <- c("method", "held_out_cohort", "n_training_cohorts")
metric_cols <- setdiff(names(perfold), id_cols)
metric_cols <- metric_cols[vapply(metric_cols,
                                  function(c) is.numeric(perfold[[c]]), logical(1))]

# Quality metrics that get a paired-Wilcoxon test (R1 P2-2: do NOT test count /
# universe-size columns like n_genes_common / n_degs_* / n_samples_training —
# a paired test on those is meaningless). Counts still get a mean in the table.
quality_metrics <- intersect(
  c("jaccard_01", "lfc_spearman", "lfc_spearman_sig", "direction_concordance",
    "auc_replication", "fisher_or", "pi1_heldout"),
  metric_cols)

# ===========================================================================
# 2. Per-method mean +/- sd across folds -> loocv_summary.csv
# ===========================================================================
summary_dt <- perfold[, {
  out <- list()
  for (mc in metric_cols) {
    v <- get(mc)
    out[[paste0(mc, "_mean")]] <- mean(v, na.rm = TRUE)
    out[[paste0(mc, "_sd")]]   <- if (sum(!is.na(v)) > 1) sd(v, na.rm = TRUE) else NA_real_
  }
  out
}, by = method]

fwrite(summary_dt, file.path(LOO_DIR, "loocv_summary.csv"))
cat("\nSaved per-method summary:", file.path(LOO_DIR, "loocv_summary.csv"), "\n")

cat("\n── Per-method mean across folds ──\n")
mean_only <- summary_dt[, c("method", paste0(metric_cols, "_mean")), with = FALSE]
print(mean_only, digits = 4)

# ===========================================================================
# 3. Head-to-head method comparison -> method_comparison.csv
#   For each metric: each method's mean + PAIRED Wilcoxon p (dream vs each other
#   method) across the shared folds. Paired because every method was run on the
#   identical held-out splits.
# ===========================================================================
# Paired comparisons: dream against every other present method.
pairs <- if ("dream" %in% methods)
  lapply(setdiff(methods, "dream"), function(o) c("dream", o)) else list()

paired_wilcox <- function(metric, a, b) {
  wide <- dcast(perfold[method %in% c(a, b)],
                held_out_cohort ~ method, value.var = metric)
  if (!all(c(a, b) %in% names(wide))) return(NA_real_)
  ok <- complete.cases(wide[[a]], wide[[b]])
  # Need >= 2 paired non-NA folds, with at least one non-zero difference.
  if (sum(ok) < 2) return(NA_real_)
  va <- wide[[a]][ok]; vb <- wide[[b]][ok]
  if (all(va - vb == 0)) return(NA_real_)
  tryCatch(
    suppressWarnings(wilcox.test(va, vb, paired = TRUE))$p.value,
    error = function(e) NA_real_)
}

comp <- rbindlist(lapply(metric_cols, function(mc) {
  row <- data.table(metric = mc)
  for (mth in methods)
    row[[paste0(mth, "_mean")]] <- mean(perfold[method == mth][[mc]], na.rm = TRUE)
  if (mc %in% quality_metrics)
    for (pr in pairs)
      row[[paste0("p_", pr[1], "_vs_", pr[2])]] <- paired_wilcox(mc, pr[1], pr[2])
  row
}), fill = TRUE)

fwrite(comp, file.path(LOO_DIR, "method_comparison.csv"))
cat("\nSaved head-to-head comparison:",
    file.path(LOO_DIR, "method_comparison.csv"), "\n")
cat("\n── Head-to-head method comparison (own-universe) ──\n")
print(comp, digits = 4)

# ===========================================================================
# 4. COMMON-UNIVERSE arm (R3 P1-1): recompute the metrics on the gene set
#    tested by ALL present methods, per fold, from the loocv_pergene_<cohort>.csv
#    files. Confirms the own-universe head-to-head isn't an artifact of the
#    methods' differing tested universes. (pi1 needs ho_pval, absent from the
#    pergene files, so pi1 is NA in the common-universe arm — the other 5
#    metrics are recomputed exactly.)
# ===========================================================================
HELPERS <- file.path(PROJECT_ROOT,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/multimethod_validation/de_validation_helpers.R")
pg_files <- list.files(LOO_DIR, pattern = "^loocv_pergene_.*\\.csv$", full.names = TRUE)
if (length(pg_files) > 0 && file.exists(HELPERS)) {
  suppressWarnings(suppressMessages(source(HELPERS)))   # for compute_loocv_metrics
  meths_cu <- intersect(methods, c("dream", "deseq2", "metafor"))
  common_perfold <- rbindlist(lapply(pg_files, function(f) {
    pg <- fread(f)
    cohort <- sub("^loocv_pergene_(.*)\\.csv$", "\\1", basename(f))
    padj_cols <- paste0(meths_cu, "_padj")
    if (!all(c(padj_cols, "ho_padj", "ho_logFC") %in% names(pg))) return(NULL)
    keep <- pg[, Reduce(`&`, lapply(padj_cols, function(pc) !is.na(get(pc)))) & !is.na(ho_padj)]
    cg <- pg[keep]
    if (nrow(cg) < 50) return(NULL)
    ho_dt <- cg[, .(gene, ho_logFC, ho_padj, ho_pval = NA_real_)]   # ho_pval absent -> pi1 NA
    rbindlist(lapply(meths_cu, function(m) {
      md <- cg[, .(gene, logFC = get(paste0(m, "_logFC")), padj = get(paste0(m, "_padj")))]
      r  <- compute_loocv_metrics(md, ho_dt, m)
      r[, `:=`(held_out_cohort = cohort, n_genes_common_universe = nrow(cg))]
      r
    }))
  }))
  if (!is.null(common_perfold) && nrow(common_perfold) > 0) {
    fwrite(common_perfold, file.path(LOO_DIR, "loocv_perfold_common.csv"))
    cu_metrics <- intersect(quality_metrics, names(common_perfold))
    comp_cu <- rbindlist(lapply(cu_metrics, function(mc) {
      row <- data.table(metric = mc)
      for (mth in meths_cu)
        row[[paste0(mth, "_mean")]] <- mean(common_perfold[method == mth][[mc]], na.rm = TRUE)
      row
    }), fill = TRUE)
    fwrite(comp_cu, file.path(LOO_DIR, "method_comparison_common.csv"))
    cat("\n── Head-to-head method comparison (COMMON universe, ",
        nrow(common_perfold) / length(meths_cu), " folds) ──\n", sep = "")
    print(comp_cu, digits = 4)
  }
}

cat("\nDone!\n")
