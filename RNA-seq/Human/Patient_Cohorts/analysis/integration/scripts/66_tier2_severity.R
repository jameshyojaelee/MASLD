#!/usr/bin/env Rscript
# 66_tier2_severity.R
# Tier 2 of the hierarchical staging cascade: 4-class severity classification.
#   Healthy(0) / MASL(1) / MASH-noFibrosis(2) / MASH-Fibrosis(3)
#
# Builds multinomial elastic net, XGBoost (via caret), and ensemble classifiers
# with leave-one-cohort-out (LOCO) cross-validation. Also runs a sensitivity
# analysis with an alternative NAS/fibrosis-based severity grouping.
#
# Inputs:
#   - rank_expression_matrix.rds: Top 3K genes (genes x samples)
#   - modeling_metadata.csv: Must contain severity4, dataset columns
#
# Outputs:
#   - tier2_loco_results.csv: Per-sample predictions across all folds
#   - tier2_model_summary.csv: Per-model per-fold metrics (macro AUROC, per-class F1)
#   - tier2_confusion.csv: Aggregated confusion matrix
#   - tier2_sensitivity_alt_grouping.csv: Comparison of primary vs alternative grouping
#
# Usage: Rscript 66_tier2_severity.R
# SLURM: cpu, 8 CPUs, 64GB RAM, 48h

suppressPackageStartupMessages({
  library(data.table)
  library(glmnet)
  library(caret)
  library(pROC)
  library(ranger)
  library(Matrix)
})

set.seed(42)

# --- Paths ---
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR <- file.path(INT, "results")
OUTDIR <- file.path(RDIR, "staging_classifier")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

cat("=== 66: Tier 2 Severity Classification ===\n")
cat("Started:", as.character(Sys.time()), "\n\n")

# --- Load data ---
cat("Loading rank expression matrix...\n")
rank_mat <- readRDS(file.path(OUTDIR, "rank_expression_matrix.rds"))
cat("  Expression:", nrow(rank_mat), "genes x", ncol(rank_mat), "samples\n")

cat("Loading modeling metadata...\n")
meta <- fread(file.path(OUTDIR, "modeling_metadata.csv"))
cat("  Metadata:", nrow(meta), "samples\n\n")

# --- Align samples ---
common_samples <- intersect(colnames(rank_mat), meta$sample_id)
meta <- meta[sample_id %in% common_samples]
meta <- meta[match(common_samples, sample_id)]
X_full <- t(rank_mat[, common_samples])
stopifnot(all(rownames(X_full) == meta$sample_id))

# --- Filter to samples with known severity4 ---
sev_mask <- meta$severity4 >= 0
meta_sev <- meta[sev_mask]
X_sev <- X_full[sev_mask, ]
y_sev <- meta_sev$severity4

cat("Severity distribution (severity4 >= 0):\n")
sev_tab <- table(y_sev)
print(sev_tab)
cat("  0 (Healthy):", sev_tab["0"], "\n")
cat("  1 (MASL):", sev_tab["1"], "\n")
cat("  2 (MASH-noFib):", sev_tab["2"], "\n")
cat("  3 (MASH-Fib):", sev_tab["3"], "\n")
cat("Total:", sum(sev_tab), "\n\n")

SEVERITY_LABELS <- c("Healthy", "MASL", "MASH_noFib", "MASH_Fib")

# --- Compute class weights (inverse frequency) ---
class_freq <- table(y_sev)
n_total <- length(y_sev)
n_classes <- length(class_freq)
class_weights <- n_total / (n_classes * class_freq)
sample_weights <- as.numeric(class_weights[as.character(y_sev)])

# ============================================================
# Helper: Compute multiclass metrics (macro AUROC, per-class F1)
# ============================================================
compute_multi_metrics <- function(y_true, prob_mat, class_labels = 0:3) {
  # prob_mat: n_samples x n_classes matrix of probabilities
  # Returns: macro AUROC, accuracy, per-class precision/recall/F1

  y_pred <- class_labels[apply(prob_mat, 1, which.max)]
  n <- length(y_true)

  # Accuracy
  accuracy <- mean(y_true == y_pred)

  # Per-class metrics
  per_class <- data.table()
  auc_vals <- numeric(length(class_labels))

  for (i in seq_along(class_labels)) {
    cl <- class_labels[i]
    cl_name <- SEVERITY_LABELS[i]

    # Binary: this class vs all others
    y_bin <- as.integer(y_true == cl)
    p_bin <- prob_mat[, i]

    # AUROC (one-vs-rest)
    roc_obj <- tryCatch(
      roc(y_bin, p_bin, levels = c(0, 1), direction = "<", quiet = TRUE),
      error = function(e) NULL
    )
    auc_vals[i] <- if (!is.null(roc_obj)) as.numeric(auc(roc_obj)) else NA_real_

    # Precision / Recall / F1
    tp <- sum(y_true == cl & y_pred == cl)
    fp <- sum(y_true != cl & y_pred == cl)
    fn <- sum(y_true == cl & y_pred != cl)
    prec <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
    rec  <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
    f1   <- if (!is.na(prec) && !is.na(rec) && (prec + rec) > 0) {
      2 * prec * rec / (prec + rec)
    } else { NA_real_ }

    per_class <- rbind(per_class, data.table(
      class = cl, class_name = cl_name,
      auroc = auc_vals[i], precision = prec, recall = rec, f1 = f1,
      n_true = sum(y_true == cl), n_pred = sum(y_pred == cl)
    ))
  }

  macro_auroc <- mean(auc_vals, na.rm = TRUE)
  macro_f1 <- mean(per_class$f1, na.rm = TRUE)

  list(
    macro_auroc = macro_auroc,
    macro_f1 = macro_f1,
    accuracy = accuracy,
    per_class = per_class
  )
}

# ============================================================
# LOCO Cross-Validation
# ============================================================
cat("=== LOCO Cross-Validation ===\n\n")

# Identify datasets with severity4 annotations
ds_counts <- meta_sev[, .N, by = dataset]
# Only use datasets that have at least 5 annotated samples
LOCO_DATASETS <- ds_counts[N >= 5, dataset]
cat("LOCO datasets:", paste(LOCO_DATASETS, collapse = ", "), "\n")

# Verify per-dataset class distribution
for (ds in LOCO_DATASETS) {
  ct <- meta_sev[dataset == ds, table(severity4)]
  cat("  ", ds, ":", paste(names(ct), ct, sep = "=", collapse = " "), "\n")
}
cat("\n")

# Filter to LOCO-eligible
loco_mask <- meta_sev$dataset %in% LOCO_DATASETS
meta_loco <- meta_sev[loco_mask]
X_loco <- X_sev[loco_mask, ]
y_loco <- meta_loco$severity4
w_loco <- sample_weights[loco_mask]

all_preds <- list()
all_metrics <- list()

for (fold_ds in LOCO_DATASETS) {
  cat("--- Fold:", fold_ds, "---\n")

  test_idx  <- which(meta_loco$dataset == fold_ds)
  train_idx <- which(meta_loco$dataset != fold_ds)

  if (length(test_idx) == 0) { cat("  SKIP: no test\n\n"); next }

  X_train <- X_loco[train_idx, ]
  y_train <- y_loco[train_idx]
  w_train <- w_loco[train_idx]
  X_test  <- X_loco[test_idx, ]
  y_test  <- y_loco[test_idx]

  cat("  Train:", length(y_train), " Test:", length(y_test), "\n")

  # Need at least 2 classes in training
  if (length(unique(y_train)) < 2) { cat("  SKIP: <2 classes\n\n"); next }

  # Factor labels for training
  y_train_fac <- factor(y_train, levels = 0:3)
  y_test_fac  <- factor(y_test, levels = 0:3)

  pred_dt <- data.table(
    sample_id = meta_loco$sample_id[test_idx],
    dataset = fold_ds,
    y_true = y_test
  )

  # ---- Model 1: Multinomial Elastic Net ----
  cat("  [1] Multinomial Elastic Net...")
  enet_prob <- matrix(NA_real_, nrow = length(y_test), ncol = 4)
  colnames(enet_prob) <- paste0("prob_enet_", 0:3)

  enet_fit <- tryCatch({
    cv.glmnet(X_train, y_train_fac, family = "multinomial", alpha = 0.5,
              weights = w_train, nfolds = 5, type.measure = "class",
              parallel = FALSE, type.multinomial = "grouped")
  }, error = function(e) { cat(" ERROR:", conditionMessage(e), "\n"); NULL })

  if (!is.null(enet_fit)) {
    # predict returns a list of class-specific arrays or a 3D array
    raw_pred <- predict(enet_fit, X_test, s = "lambda.min", type = "response")
    # Handle the 3D array: dim = (n_samples, n_classes, 1)
    if (is.array(raw_pred) && length(dim(raw_pred)) == 3) {
      enet_prob <- raw_pred[, , 1]
    } else {
      enet_prob <- as.matrix(raw_pred)
    }
    colnames(enet_prob) <- paste0("prob_enet_", 0:3)
    m1 <- compute_multi_metrics(y_test, enet_prob)
    cat(" macroAUROC=", round(m1$macro_auroc, 3), " acc=", round(m1$accuracy, 3), "\n")

    m1_row <- data.table(model = "multinomial_enet", fold = fold_ds,
                         macro_auroc = m1$macro_auroc, macro_f1 = m1$macro_f1,
                         accuracy = m1$accuracy)
    all_metrics <- c(all_metrics, list(m1_row))
  } else {
    cat(" FAILED\n")
  }
  pred_dt <- cbind(pred_dt, as.data.table(enet_prob))

  # ---- Model 2: XGBoost multiclass (via caret) ----
  cat("  [2] XGBoost (caret multi:softprob)...")
  xgb_prob <- matrix(NA_real_, nrow = length(y_test), ncol = 4)
  colnames(xgb_prob) <- paste0("prob_xgb_", 0:3)

  xgb_fit <- tryCatch({
    # caret needs factor levels as valid R names
    y_tr_fac <- factor(paste0("S", y_train), levels = paste0("S", 0:3))
    xgb_grid <- expand.grid(
      nrounds = c(100, 200), max_depth = c(3, 6), eta = c(0.05, 0.1),
      gamma = 0, colsample_bytree = 0.8, min_child_weight = 1, subsample = 0.8
    )
    ctrl <- trainControl(method = "cv", number = 3,
                         classProbs = TRUE, summaryFunction = multiClassSummary,
                         verboseIter = FALSE, allowParallel = FALSE)
    train(x = as.data.frame(X_train), y = y_tr_fac,
          method = "xgbTree", trControl = ctrl, tuneGrid = xgb_grid,
          metric = "AUC", weights = w_train, verbosity = 0)
  }, error = function(e) { cat(" ERROR:", conditionMessage(e), "\n"); NULL })

  if (!is.null(xgb_fit)) {
    xgb_pred_df <- predict(xgb_fit, newdata = as.data.frame(X_test), type = "prob")
    # Ensure correct column order: S0, S1, S2, S3
    xgb_prob <- as.matrix(xgb_pred_df[, paste0("S", 0:3)])
    colnames(xgb_prob) <- paste0("prob_xgb_", 0:3)
    m2 <- compute_multi_metrics(y_test, xgb_prob)
    cat(" macroAUROC=", round(m2$macro_auroc, 3), " acc=", round(m2$accuracy, 3), "\n")

    m2_row <- data.table(model = "xgboost", fold = fold_ds,
                         macro_auroc = m2$macro_auroc, macro_f1 = m2$macro_f1,
                         accuracy = m2$accuracy)
    all_metrics <- c(all_metrics, list(m2_row))
  } else {
    cat(" FAILED\n")
  }
  pred_dt <- cbind(pred_dt, as.data.table(xgb_prob))

  # ---- Ensemble (average probabilities) ----
  ens_prob <- matrix(NA_real_, nrow = length(y_test), ncol = 4)
  for (k in 1:4) {
    enet_col <- paste0("prob_enet_", k - 1)
    xgb_col  <- paste0("prob_xgb_", k - 1)
    vals <- cbind(
      if (enet_col %in% names(pred_dt)) pred_dt[[enet_col]] else NA,
      if (xgb_col %in% names(pred_dt)) pred_dt[[xgb_col]] else NA
    )
    ens_prob[, k] <- rowMeans(vals, na.rm = TRUE)
  }
  # Renormalize
  ens_prob <- ens_prob / rowSums(ens_prob)
  colnames(ens_prob) <- paste0("prob_ens_", 0:3)
  pred_dt <- cbind(pred_dt, as.data.table(ens_prob))
  pred_dt[, pred_class := apply(ens_prob, 1, which.max) - 1]

  me <- compute_multi_metrics(y_test, ens_prob)
  cat("  [E] Ensemble macroAUROC=", round(me$macro_auroc, 3),
      " acc=", round(me$accuracy, 3), "\n\n")
  me_row <- data.table(model = "ensemble", fold = fold_ds,
                       macro_auroc = me$macro_auroc, macro_f1 = me$macro_f1,
                       accuracy = me$accuracy)
  all_metrics <- c(all_metrics, list(me_row))

  all_preds <- c(all_preds, list(pred_dt))
}

# ============================================================
# Aggregate and Save
# ============================================================
cat("=== Aggregating Results ===\n\n")

pred_results <- rbindlist(all_preds, fill = TRUE)
fwrite(pred_results, file.path(OUTDIR, "tier2_loco_results.csv"))
cat("Per-sample predictions saved:", nrow(pred_results), "\n")

metric_results <- rbindlist(all_metrics, fill = TRUE)
fwrite(metric_results, file.path(OUTDIR, "tier2_model_summary.csv"))

# Model-level summary
cat("\nMean metrics by model:\n")
model_avg <- metric_results[, .(
  mean_macro_auroc = mean(macro_auroc, na.rm = TRUE),
  sd_macro_auroc   = sd(macro_auroc, na.rm = TRUE),
  mean_macro_f1    = mean(macro_f1, na.rm = TRUE),
  mean_accuracy    = mean(accuracy, na.rm = TRUE),
  n_folds = .N
), by = model]
print(model_avg)

# Confusion matrix (from ensemble predictions)
cat("\nConfusion matrix (ensemble, aggregated LOCO):\n")
cm_dt <- pred_results[, .N, by = .(y_true, pred_class)]
cm_dt[, y_true_label := SEVERITY_LABELS[y_true + 1]]
cm_dt[, pred_label   := SEVERITY_LABELS[pred_class + 1]]
cm_wide <- dcast(cm_dt, y_true_label ~ pred_label, value.var = "N", fill = 0)
print(cm_wide)
fwrite(cm_dt, file.path(OUTDIR, "tier2_confusion.csv"))

# ============================================================
# Sensitivity Analysis: Alternative Grouping
# ============================================================
cat("\n=== Sensitivity Analysis: Alternative Grouping ===\n")
cat("Alt: MASL = NAS<3, MASH-early = NAS 3-4 + F0-1, MASH-advanced = NAS 5+ OR F>=2\n\n")

# Recompute severity with alternative grouping
meta_alt <- copy(meta)
meta_alt[, severity4_alt := {
  sev <- rep(-1L, .N)
  # Healthy = controls
  sev[group_binary == "Control"] <- 0L
  # Among disease samples with NAS:
  dis_mask <- group_binary == "Disease"
  has_nas <- !is.na(nas_score) & nas_score >= 0
  has_fib <- fib_stage >= 0
  # MASL: NAS < 3
  sev[dis_mask & has_nas & nas_score < 3] <- 1L
  # MASH-early: NAS 3-4 AND F0-1
  sev[dis_mask & has_nas & has_fib & nas_score >= 3 & nas_score <= 4 & fib_stage <= 1] <- 2L
  # MASH-advanced: NAS >= 5 OR F >= 2
  sev[dis_mask & has_nas & nas_score >= 5] <- 3L
  sev[dis_mask & has_fib & fib_stage >= 2] <- 3L
  sev
}]

# Resolve conflicts: if both MASH-early and MASH-advanced were assigned, keep advanced
# (the last assignment wins above since F>=2 overrides NAS 3-4 F0-1)

# Run LOCO with alternative grouping
alt_mask <- meta_alt$severity4_alt >= 0 & meta_alt$sample_id %in% common_samples
meta_alt_sev <- meta_alt[alt_mask]
X_alt_sev <- X_full[match(meta_alt_sev$sample_id, common_samples), ]
y_alt_sev <- meta_alt_sev$severity4_alt

cat("Alternative grouping distribution:\n")
print(table(y_alt_sev))

# Compute class weights
alt_freq <- table(y_alt_sev)
alt_weights <- length(y_alt_sev) / (length(alt_freq) * alt_freq)
alt_sample_w <- as.numeric(alt_weights[as.character(y_alt_sev)])

ALT_LOCO_DS <- meta_alt_sev[, .N, by = dataset][N >= 5, dataset]
cat("Alt LOCO datasets:", paste(ALT_LOCO_DS, collapse = ", "), "\n\n")

alt_metrics <- list()

for (fold_ds in ALT_LOCO_DS) {
  loco_mask_alt <- meta_alt_sev$dataset %in% ALT_LOCO_DS
  test_idx  <- which(meta_alt_sev$dataset == fold_ds & loco_mask_alt)
  train_idx <- which(meta_alt_sev$dataset != fold_ds & loco_mask_alt)

  if (length(test_idx) == 0 || length(unique(y_alt_sev[train_idx])) < 2) next

  X_tr <- X_alt_sev[train_idx, ]
  y_tr <- factor(y_alt_sev[train_idx], levels = 0:3)
  w_tr <- alt_sample_w[train_idx]
  X_te <- X_alt_sev[test_idx, ]
  y_te <- y_alt_sev[test_idx]

  enet_fit <- tryCatch(
    cv.glmnet(X_tr, y_tr, family = "multinomial", alpha = 0.5,
              weights = w_tr, nfolds = 5, type.measure = "class",
              type.multinomial = "grouped"),
    error = function(e) NULL
  )

  if (!is.null(enet_fit)) {
    raw_p <- predict(enet_fit, X_te, s = "lambda.min", type = "response")
    if (is.array(raw_p) && length(dim(raw_p)) == 3) raw_p <- raw_p[, , 1]
    m <- compute_multi_metrics(y_te, as.matrix(raw_p))
    alt_metrics <- c(alt_metrics, list(data.table(
      grouping = "alternative", fold = fold_ds,
      macro_auroc = m$macro_auroc, macro_f1 = m$macro_f1, accuracy = m$accuracy
    )))
  }
}

# Combine with primary metrics
primary_enet <- metric_results[model == "multinomial_enet"]
primary_enet[, grouping := "primary"]

alt_agg <- rbindlist(alt_metrics, fill = TRUE)
if (nrow(alt_agg) > 0 && nrow(primary_enet) > 0) {
  sensitivity_cmp <- rbind(
    primary_enet[, .(grouping, fold, macro_auroc, macro_f1, accuracy)],
    alt_agg,
    fill = TRUE
  )

  cat("\nSensitivity comparison:\n")
  sens_summary <- sensitivity_cmp[, .(
    mean_macro_auroc = mean(macro_auroc, na.rm = TRUE),
    mean_macro_f1    = mean(macro_f1, na.rm = TRUE),
    mean_accuracy    = mean(accuracy, na.rm = TRUE),
    n_folds = .N
  ), by = grouping]
  print(sens_summary)

  fwrite(sensitivity_cmp, file.path(OUTDIR, "tier2_sensitivity_alt_grouping.csv"))
  cat("Sensitivity results saved\n")
} else {
  cat("Sensitivity analysis: insufficient data for comparison\n")
  fwrite(data.table(note = "insufficient_data"), file.path(OUTDIR, "tier2_sensitivity_alt_grouping.csv"))
}

# ============================================================
# Summary
# ============================================================
cat("\n=== FINAL SUMMARY ===\n")
cat("Total LOCO samples:", nrow(pred_results), "\n")
cat("Folds:", length(unique(pred_results$dataset)), "\n")
cat("Best model (macro AUROC):", model_avg[which.max(mean_macro_auroc), model], "\n")
cat("Best macro AUROC:", round(max(model_avg$mean_macro_auroc, na.rm = TRUE), 4), "\n")
cat("Overall accuracy (ensemble):", round(mean(pred_results$y_true == pred_results$pred_class, na.rm = TRUE), 4), "\n")

cat("\n=== 66_tier2_severity.R completed:", as.character(Sys.time()), "===\n")
