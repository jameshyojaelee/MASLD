#!/usr/bin/env Rscript
# =============================================================================
# 88_leave_two_out.R
#
# Leave-two-cohorts-out (LTO) stress test for the V3 fibrosis staging
# classifier (F>=3 binary). With 6 fibrosis-eligible cohorts, tests all
# (6 choose 2) = 15 held-out pairs to assess robustness when 25-40% of
# training data is removed simultaneously.
#
# Pipeline per LTO combination:
#   1. Hold out 2 cohorts (test set)
#   2. Train on remaining 4 using within-fold variance selection (top 3K) +
#      elastic net for F>=3 (same as Script 76)
#   3. Predict combined held-out set, compute AUROC
#
# Inputs:
#   - results/integration/merged_dge.rds (34,453 genes x 1,444 samples)
#   - results/staging_classifier/modeling_metadata.csv
#   - results/staging_classifier/v3_proper_cv_results.csv (reference AUROC)
#
# Outputs (to results/staging_classifier/):
#   - lto_results.csv: 15 rows with pair, auroc, n_train, n_test
#   - lto_summary.csv: Aggregate statistics (mean, min, max, worst pair)
#
# Usage: Rscript 88_leave_two_out.R
# SLURM: bigmem, 16 CPUs, 200GB RAM, 48h
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(glmnet)
  library(pROC)
})

set.seed(42)

# --- Paths -------------------------------------------------------------------
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT    <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR   <- file.path(INT, "results")
OUTDIR <- file.path(RDIR, "staging_classifier")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

cat("=== 88: Leave-Two-Cohorts-Out (LTO) Stress Test ===\n")
cat("Started:", as.character(Sys.time()), "\n\n")

K_FEATURES <- 3000
ALPHA <- 0.5

# =============================================================================
# LOAD DATA
# =============================================================================
cat("Loading merged DGE...\n")
dge <- readRDS(file.path(RDIR, "integration/merged_dge.rds"))
dge <- calcNormFactors(dge, method = "TMM")
logcpm <- cpm(dge, log = TRUE, prior.count = 1)
cat("  Expression:", nrow(logcpm), "genes x", ncol(logcpm), "samples\n")

cat("Loading modeling metadata...\n")
meta <- fread(file.path(OUTDIR, "modeling_metadata.csv"))

# Align
common <- intersect(colnames(logcpm), meta$sample_id)
meta <- meta[match(common, sample_id)]
logcpm <- logcpm[, common]
cat("  Aligned:", ncol(logcpm), "samples\n")

# Subset to fibrosis-eligible (fib_ge3 target, non-excluded, non-NA)
fib_valid <- meta$loco_fold_fibrosis != "excluded" &
             !is.na(meta$fib_ge3) &
             meta$fib_ge3 >= 0
meta_fib <- meta[fib_valid]
logcpm_fib <- logcpm[, meta_fib$sample_id]

cohorts <- sort(unique(meta_fib$loco_fold_fibrosis))
cohorts <- cohorts[cohorts != "excluded" & cohorts != "NA"]
cat("  Fibrosis-eligible:", nrow(meta_fib), "samples\n")
cat("  Cohorts:", paste(cohorts, collapse = ", "), "\n")

# Per-cohort sample counts
for (coh in cohorts) {
  n_coh <- sum(meta_fib$loco_fold_fibrosis == coh)
  n_pos <- sum(meta_fib$fib_ge3[meta_fib$loco_fold_fibrosis == coh] == 1)
  cat("    ", coh, ":", n_coh, "samples (", n_pos, "F>=3)\n")
}

# Reference AUROC from V3
v3_file <- file.path(OUTDIR, "v3_proper_cv_results.csv")
ref_auroc <- NA_real_
if (file.exists(v3_file)) {
  v3_preds <- fread(v3_file)
  v3_fib <- v3_preds[target == "fib_ge3"]
  ref_auroc <- tryCatch(
    as.numeric(auc(roc(v3_fib$true_label, v3_fib$probability, quiet = TRUE))),
    error = function(e) NA_real_
  )
  cat("\n  V3 reference AUROC (LOCO, 6-fold):", round(ref_auroc, 4), "\n")
}

# =============================================================================
# HELPER: Within-fold feature selection + rank transform
# (Identical to Script 76)
# =============================================================================
select_features_and_transform <- function(logcpm_train, logcpm_test, k = K_FEATURES) {
  gene_var <- apply(logcpm_train, 1, var)
  top_k_genes <- names(sort(gene_var, decreasing = TRUE))[1:min(k, length(gene_var))]

  rank_train <- apply(logcpm_train[top_k_genes, , drop = FALSE], 2, function(x) {
    r <- rank(x, ties.method = "average"); r / length(r)
  })
  rank_test <- apply(logcpm_test[top_k_genes, , drop = FALSE], 2, function(x) {
    r <- rank(x, ties.method = "average"); r / length(r)
  })

  list(
    X_train = t(rank_train),
    X_test  = t(rank_test),
    selected_genes = top_k_genes
  )
}

# =============================================================================
# GENERATE ALL (6 choose 2) = 15 LTO PAIRS
# =============================================================================
pairs <- combn(cohorts, 2, simplify = FALSE)
cat("\n  LTO pairs:", length(pairs), "\n\n")

# =============================================================================
# RUN LTO FOR EACH PAIR
# =============================================================================
lto_results <- list()

for (p_idx in seq_along(pairs)) {
  held_out <- pairs[[p_idx]]
  pair_label <- paste(held_out, collapse = " + ")
  cat("  Pair", p_idx, "/", length(pairs), ":", pair_label, "... ")

  # Split train / test
  test_mask  <- meta_fib$loco_fold_fibrosis %in% held_out
  train_mask <- !test_mask

  n_train <- sum(train_mask)
  n_test  <- sum(test_mask)

  y_train <- meta_fib$fib_ge3[train_mask]
  y_test  <- meta_fib$fib_ge3[test_mask]

  # Check viability
  if (n_test < 10 || n_train < 30) {
    cat("SKIP (too few samples: train=", n_train, ", test=", n_test, ")\n")
    lto_results[[p_idx]] <- data.table(
      pair = pair_label,
      cohort_1 = held_out[1],
      cohort_2 = held_out[2],
      auroc = NA_real_,
      n_train = n_train,
      n_test = n_test,
      n_train_pos = sum(y_train == 1),
      n_test_pos = sum(y_test == 1),
      pct_data_held_out = round(n_test / nrow(meta_fib) * 100, 1),
      status = "skipped_too_few"
    )
    next
  }

  if (length(unique(y_test)) < 2) {
    cat("SKIP (single class in test)\n")
    lto_results[[p_idx]] <- data.table(
      pair = pair_label,
      cohort_1 = held_out[1],
      cohort_2 = held_out[2],
      auroc = NA_real_,
      n_train = n_train,
      n_test = n_test,
      n_train_pos = sum(y_train == 1),
      n_test_pos = sum(y_test == 1),
      pct_data_held_out = round(n_test / nrow(meta_fib) * 100, 1),
      status = "skipped_single_class"
    )
    next
  }

  if (length(unique(y_train)) < 2) {
    cat("SKIP (single class in train)\n")
    lto_results[[p_idx]] <- data.table(
      pair = pair_label,
      cohort_1 = held_out[1],
      cohort_2 = held_out[2],
      auroc = NA_real_,
      n_train = n_train,
      n_test = n_test,
      n_train_pos = sum(y_train == 1),
      n_test_pos = sum(y_test == 1),
      pct_data_held_out = round(n_test / nrow(meta_fib) * 100, 1),
      status = "skipped_single_class_train"
    )
    next
  }

  # Within-fold feature selection on TRAINING data only
  feat <- select_features_and_transform(
    logcpm_fib[, meta_fib$sample_id[train_mask]],
    logcpm_fib[, meta_fib$sample_id[test_mask]],
    k = K_FEATURES
  )

  # Class-weighted elastic net
  class_counts <- table(y_train)
  class_wts <- 1 / class_counts
  class_wts <- class_wts / sum(class_wts) * length(class_wts)
  sample_wts <- class_wts[as.character(y_train)]

  fit <- tryCatch({
    cv.glmnet(
      x = feat$X_train, y = factor(y_train),
      family = "binomial", alpha = ALPHA, nfolds = 5,
      type.measure = "auc", weights = as.numeric(sample_wts)
    )
  }, error = function(e) { cat("ERROR:", e$message, "\n"); NULL })

  if (is.null(fit)) {
    lto_results[[p_idx]] <- data.table(
      pair = pair_label,
      cohort_1 = held_out[1],
      cohort_2 = held_out[2],
      auroc = NA_real_,
      n_train = n_train,
      n_test = n_test,
      n_train_pos = sum(y_train == 1),
      n_test_pos = sum(y_test == 1),
      pct_data_held_out = round(n_test / nrow(meta_fib) * 100, 1),
      status = "model_failed"
    )
    next
  }

  prob <- as.numeric(predict(fit, newx = feat$X_test, s = "lambda.min",
                              type = "response"))
  pred <- ifelse(prob > 0.5, 1, 0)

  auroc_lto <- tryCatch(
    as.numeric(auc(roc(y_test, prob, quiet = TRUE))),
    error = function(e) NA_real_
  )
  acc <- mean(pred == y_test)

  # Per-cohort AUROC within the held-out pair
  per_cohort <- list()
  for (coh in held_out) {
    coh_mask <- meta_fib$loco_fold_fibrosis[test_mask] == coh
    if (sum(coh_mask) >= 5 && length(unique(y_test[coh_mask])) >= 2) {
      per_cohort[[coh]] <- tryCatch(
        as.numeric(auc(roc(y_test[coh_mask], prob[coh_mask], quiet = TRUE))),
        error = function(e) NA_real_
      )
    } else {
      per_cohort[[coh]] <- NA_real_
    }
  }

  # Count non-zero coefficients
  coefs <- coef(fit, s = "lambda.min")
  n_nonzero <- sum(coefs[, 1] != 0) - 1  # subtract intercept

  lto_results[[p_idx]] <- data.table(
    pair = pair_label,
    cohort_1 = held_out[1],
    cohort_2 = held_out[2],
    auroc = auroc_lto,
    accuracy = acc,
    n_train = n_train,
    n_test = n_test,
    n_train_pos = sum(y_train == 1),
    n_test_pos = sum(y_test == 1),
    pct_data_held_out = round(n_test / nrow(meta_fib) * 100, 1),
    auroc_cohort_1 = per_cohort[[held_out[1]]],
    auroc_cohort_2 = per_cohort[[held_out[2]]],
    n_nonzero_coefs = n_nonzero,
    status = "success"
  )

  cat("AUROC=", round(auroc_lto, 3),
      " (", held_out[1], ":", round(per_cohort[[held_out[1]]], 3),
      ", ", held_out[2], ":", round(per_cohort[[held_out[2]]], 3),
      ") Acc=", round(acc, 3),
      " [train=", n_train, " test=", n_test, "]\n")
}

# =============================================================================
# AGGREGATE RESULTS
# =============================================================================
lto_dt <- rbindlist(lto_results, fill = TRUE)
fwrite(lto_dt, file.path(OUTDIR, "lto_results.csv"))

cat("\n=== LTO Results (15 pairs) ===\n")
print(lto_dt[, .(pair, auroc, n_train, n_test, pct_data_held_out, status)])

# Summary
valid_lto <- lto_dt[status == "success"]

if (nrow(valid_lto) > 0) {
  worst_pair <- valid_lto[which.min(auroc)]
  best_pair  <- valid_lto[which.max(auroc)]

  lto_summary <- data.table(
    n_pairs_total = nrow(lto_dt),
    n_pairs_valid = nrow(valid_lto),
    ref_auroc_v3_loco = ref_auroc,
    mean_lto_auroc = mean(valid_lto$auroc, na.rm = TRUE),
    sd_lto_auroc = sd(valid_lto$auroc, na.rm = TRUE),
    min_lto_auroc = min(valid_lto$auroc, na.rm = TRUE),
    max_lto_auroc = max(valid_lto$auroc, na.rm = TRUE),
    median_lto_auroc = median(valid_lto$auroc, na.rm = TRUE),
    worst_pair = worst_pair$pair,
    worst_auroc = worst_pair$auroc,
    worst_pct_held_out = worst_pair$pct_data_held_out,
    best_pair = best_pair$pair,
    best_auroc = best_pair$auroc,
    best_pct_held_out = best_pair$pct_data_held_out,
    mean_pct_held_out = mean(valid_lto$pct_data_held_out),
    mean_accuracy = mean(valid_lto$accuracy, na.rm = TRUE),
    auroc_drop_vs_v3 = ref_auroc - mean(valid_lto$auroc, na.rm = TRUE)
  )

  fwrite(lto_summary, file.path(OUTDIR, "lto_summary.csv"))

  cat("\n=== LTO Summary ===\n")
  cat("  Valid pairs:", nrow(valid_lto), "/", nrow(lto_dt), "\n")
  cat("  V3 LOCO reference AUROC:", round(ref_auroc, 4), "\n")
  cat("  LTO mean AUROC:", round(lto_summary$mean_lto_auroc, 4),
      " +/- ", round(lto_summary$sd_lto_auroc, 4), "\n")
  cat("  LTO range: [", round(lto_summary$min_lto_auroc, 4), ",",
      round(lto_summary$max_lto_auroc, 4), "]\n")
  cat("  AUROC drop vs V3:", round(lto_summary$auroc_drop_vs_v3, 4), "\n")
  cat("  Worst pair:", lto_summary$worst_pair,
      "(AUROC=", round(lto_summary$worst_auroc, 4),
      ", held out", lto_summary$worst_pct_held_out, "%)\n")
  cat("  Best pair:", lto_summary$best_pair,
      "(AUROC=", round(lto_summary$best_auroc, 4),
      ", held out", lto_summary$best_pct_held_out, "%)\n")
  cat("  Mean % data held out:", round(lto_summary$mean_pct_held_out, 1), "%\n")

  # Cohort-level degradation analysis: which cohort's removal hurts most?
  cat("\n=== Per-Cohort Impact ===\n")
  cohort_impact <- list()
  for (coh in cohorts) {
    pairs_with_coh <- valid_lto[cohort_1 == coh | cohort_2 == coh]
    pairs_without_coh <- valid_lto[cohort_1 != coh & cohort_2 != coh]
    if (nrow(pairs_with_coh) > 0 && nrow(pairs_without_coh) > 0) {
      cohort_impact[[coh]] <- data.table(
        cohort = coh,
        n_pairs_held_out = nrow(pairs_with_coh),
        mean_auroc_when_held_out = mean(pairs_with_coh$auroc, na.rm = TRUE),
        mean_auroc_when_retained = mean(pairs_without_coh$auroc, na.rm = TRUE),
        degradation = mean(pairs_without_coh$auroc, na.rm = TRUE) -
                      mean(pairs_with_coh$auroc, na.rm = TRUE)
      )
    }
  }

  if (length(cohort_impact) > 0) {
    cohort_dt <- rbindlist(cohort_impact)
    cohort_dt <- cohort_dt[order(-degradation)]
    print(cohort_dt)

    most_critical <- cohort_dt[1]
    cat("\n  Most critical cohort:", most_critical$cohort,
        "— removing it degrades AUROC by",
        round(most_critical$degradation, 4), "on average\n")
  }
} else {
  cat("\n  WARNING: No valid LTO pairs — all skipped or failed\n")
  lto_summary <- data.table(
    n_pairs_total = nrow(lto_dt),
    n_pairs_valid = 0
  )
  fwrite(lto_summary, file.path(OUTDIR, "lto_summary.csv"))
}

cat("\nOutputs written to:", OUTDIR, "\n")
cat("  - lto_results.csv (", nrow(lto_dt), "rows)\n")
cat("  - lto_summary.csv\n")

cat("\n=== 88_leave_two_out.R completed:", as.character(Sys.time()), "===\n")
