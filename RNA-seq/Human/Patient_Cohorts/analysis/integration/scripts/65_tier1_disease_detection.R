#!/usr/bin/env Rscript
# 65_tier1_disease_detection.R
# Tier 1 of the hierarchical staging cascade: Disease vs Control classifier.
#
# Builds logistic elastic net, random forest, and XGBoost (via caret) classifiers
# with leave-one-cohort-out (LOCO) cross-validation across 8 cohorts that have
# both disease and control samples.
#
# Inputs:
#   - rank_expression_matrix.rds: Within-sample rank-transformed expression (genes x samples)
#   - modeling_metadata.csv: Per-sample metadata with is_disease, dataset, etc.
#
# Outputs:
#   - tier1_loco_results.csv: Per-sample predictions across all folds
#   - tier1_model_summary.csv: Per-model per-fold metrics
#   - tier1_top_features.csv: Top 50 genes by mean |coefficient| from elastic net
#   - tier1_best_model.rds: Best model (highest mean LOCO AUROC) retrained on all data
#
# Usage: Rscript 65_tier1_disease_detection.R
# SLURM: cpu, 8 CPUs, 64GB RAM, 48h

suppressPackageStartupMessages({
  library(data.table)
  library(glmnet)
  library(ranger)
  library(caret)
  library(pROC)
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

cat("=== 65: Tier 1 Disease Detection Classifier ===\n")
cat("Started:", as.character(Sys.time()), "\n\n")

# --- Load data ---
cat("Loading rank expression matrix...\n")
rank_mat <- readRDS(file.path(OUTDIR, "rank_expression_matrix.rds"))
cat("  Expression matrix:", nrow(rank_mat), "genes x", ncol(rank_mat), "samples\n")

cat("Loading modeling metadata...\n")
meta <- fread(file.path(OUTDIR, "modeling_metadata.csv"))
cat("  Metadata:", nrow(meta), "samples\n\n")

# --- Align samples ---
common_samples <- intersect(colnames(rank_mat), meta$sample_id)
cat("Common samples:", length(common_samples), "\n")
meta <- meta[sample_id %in% common_samples]
meta <- meta[match(common_samples, sample_id)]
X_full <- t(rank_mat[, common_samples])  # samples x genes
stopifnot(all(rownames(X_full) == meta$sample_id))

# --- Define LOCO folds ---
# 8 cohorts with both disease and control samples
LOCO_DATASETS <- c("GSE130970", "GSE135251", "GSE162694", "GSE174478",
                    "GSE193066", "GSE240729", "GSE213621", "GSE126848")

# Verify each fold has both classes
cat("\nDataset composition:\n")
for (ds in LOCO_DATASETS) {
  n_dis <- sum(meta$dataset == ds & meta$is_disease == 1)
  n_ctl <- sum(meta$dataset == ds & meta$is_disease == 0)
  cat("  ", ds, ": Disease=", n_dis, " Control=", n_ctl, "\n")
}

# Filter to LOCO-eligible samples
loco_mask <- meta$dataset %in% LOCO_DATASETS
meta_loco <- meta[loco_mask]
X_loco <- X_full[loco_mask, ]
y_loco <- meta_loco$is_disease
cat("\nLOCO-eligible samples:", nrow(meta_loco),
    "(Disease:", sum(y_loco == 1), "Control:", sum(y_loco == 0), ")\n\n")

# ============================================================
# Helper: Compute classification metrics
# ============================================================
compute_metrics <- function(y_true, y_pred_prob, y_pred_class = NULL, pos_label = 1) {
  if (is.null(y_pred_class)) {
    y_pred_class <- ifelse(y_pred_prob >= 0.5, 1, 0)
  }

  # AUROC
  roc_obj <- tryCatch(
    roc(y_true, y_pred_prob, levels = c(0, 1), direction = "<", quiet = TRUE),
    error = function(e) NULL
  )
  auroc <- if (!is.null(roc_obj)) as.numeric(auc(roc_obj)) else NA_real_

  # Confusion matrix
  tp <- sum(y_true == 1 & y_pred_class == 1)
  tn <- sum(y_true == 0 & y_pred_class == 0)
  fp <- sum(y_true == 0 & y_pred_class == 1)
  fn <- sum(y_true == 1 & y_pred_class == 0)

  accuracy    <- (tp + tn) / length(y_true)
  sensitivity <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
  specificity <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
  precision   <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
  f1          <- if (!is.na(precision) && !is.na(sensitivity) && (precision + sensitivity) > 0) {
    2 * precision * sensitivity / (precision + sensitivity)
  } else { NA_real_ }

  data.table(
    auroc = auroc, accuracy = accuracy,
    sensitivity = sensitivity, specificity = specificity,
    precision = precision, f1 = f1,
    tp = tp, tn = tn, fp = fp, fn = fn
  )
}

# ============================================================
# LOCO Cross-Validation
# ============================================================
cat("=== LOCO Cross-Validation ===\n\n")

all_preds <- list()
all_metrics <- list()

for (fold_ds in LOCO_DATASETS) {
  cat("--- Fold:", fold_ds, "---\n")

  # Split
  test_idx  <- which(meta_loco$dataset == fold_ds)
  train_idx <- which(meta_loco$dataset != fold_ds)

  if (length(test_idx) == 0) {
    cat("  SKIP: no test samples\n\n")
    next
  }

  X_train <- X_loco[train_idx, ]
  y_train <- y_loco[train_idx]
  X_test  <- X_loco[test_idx, ]
  y_test  <- y_loco[test_idx]

  cat("  Train:", length(y_train), "(Dis:", sum(y_train == 1), "Ctl:", sum(y_train == 0), ")\n")
  cat("  Test: ", length(y_test),  "(Dis:", sum(y_test == 1),  "Ctl:", sum(y_test == 0),  ")\n")

  # Check that both classes exist in training
  if (length(unique(y_train)) < 2) {
    cat("  SKIP: training set has only one class\n\n")
    next
  }

  # --- Class weights for training ---
  n_pos <- sum(y_train == 1)
  n_neg <- sum(y_train == 0)
  w_train <- ifelse(y_train == 1, length(y_train) / (2 * n_pos),
                                   length(y_train) / (2 * n_neg))

  # ---- Model 1: Elastic Net (alpha = 0.5) ----
  cat("  [1] Elastic Net (glmnet)...")
  enet_fit <- tryCatch({
    cv.glmnet(X_train, y_train, family = "binomial", alpha = 0.5,
              weights = w_train, nfolds = 5, type.measure = "auc",
              parallel = FALSE)
  }, error = function(e) { cat(" ERROR:", conditionMessage(e), "\n"); NULL })

  if (!is.null(enet_fit)) {
    enet_prob <- as.numeric(predict(enet_fit, X_test, s = "lambda.min", type = "response"))
    enet_metrics <- compute_metrics(y_test, enet_prob)
    enet_metrics[, model := "elastic_net"]
    enet_metrics[, fold := fold_ds]
    all_metrics <- c(all_metrics, list(enet_metrics))

    pred_dt <- data.table(
      sample_id = meta_loco$sample_id[test_idx],
      dataset = fold_ds,
      y_true = y_test,
      prob_elastic_net = enet_prob
    )
    cat(" AUROC=", round(enet_metrics$auroc, 3), "\n")
  } else {
    pred_dt <- data.table(
      sample_id = meta_loco$sample_id[test_idx],
      dataset = fold_ds,
      y_true = y_test,
      prob_elastic_net = NA_real_
    )
  }

  # ---- Model 2: Random Forest (ranger) ----
  cat("  [2] Random Forest (ranger)...")
  rf_fit <- tryCatch({
    ranger(y = factor(y_train, levels = c(0, 1)),
           x = as.data.frame(X_train),
           num.trees = 500,
           mtry = floor(sqrt(ncol(X_train))),
           probability = TRUE,
           case.weights = w_train,
           num.threads = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4")),
           seed = 42)
  }, error = function(e) { cat(" ERROR:", conditionMessage(e), "\n"); NULL })

  if (!is.null(rf_fit)) {
    rf_pred <- predict(rf_fit, data = as.data.frame(X_test))
    rf_prob <- rf_pred$predictions[, "1"]
    rf_metrics <- compute_metrics(y_test, rf_prob)
    rf_metrics[, model := "random_forest"]
    rf_metrics[, fold := fold_ds]
    all_metrics <- c(all_metrics, list(rf_metrics))
    pred_dt[, prob_random_forest := rf_prob]
    cat(" AUROC=", round(rf_metrics$auroc, 3), "\n")
  } else {
    pred_dt[, prob_random_forest := NA_real_]
  }

  # ---- Model 3: XGBoost (via caret) ----
  cat("  [3] XGBoost (caret::xgbTree)...")
  xgb_fit <- tryCatch({
    # caret xgbTree with a grid search
    xgb_grid <- expand.grid(
      nrounds = c(100, 200),
      max_depth = c(3, 6),
      eta = c(0.05, 0.1),
      gamma = 0,
      colsample_bytree = 0.8,
      min_child_weight = 1,
      subsample = 0.8
    )
    ctrl <- trainControl(
      method = "cv", number = 3,
      classProbs = TRUE, summaryFunction = twoClassSummary,
      verboseIter = FALSE, allowParallel = FALSE
    )
    # caret requires factor levels as valid R names
    y_train_fac <- factor(ifelse(y_train == 1, "Disease", "Control"),
                          levels = c("Control", "Disease"))
    train(x = as.data.frame(X_train), y = y_train_fac,
          method = "xgbTree", trControl = ctrl, tuneGrid = xgb_grid,
          metric = "ROC", weights = w_train, verbosity = 0)
  }, error = function(e) { cat(" ERROR:", conditionMessage(e), "\n"); NULL })

  if (!is.null(xgb_fit)) {
    xgb_prob_df <- predict(xgb_fit, newdata = as.data.frame(X_test), type = "prob")
    xgb_prob <- xgb_prob_df[, "Disease"]
    xgb_metrics <- compute_metrics(y_test, xgb_prob)
    xgb_metrics[, model := "xgboost"]
    xgb_metrics[, fold := fold_ds]
    all_metrics <- c(all_metrics, list(xgb_metrics))
    pred_dt[, prob_xgboost := xgb_prob]
    cat(" AUROC=", round(xgb_metrics$auroc, 3), "\n")
  } else {
    pred_dt[, prob_xgboost := NA_real_]
  }

  # ---- Ensemble (average probabilities) ----
  prob_cols <- grep("^prob_", names(pred_dt), value = TRUE)
  prob_mat <- as.matrix(pred_dt[, ..prob_cols])
  pred_dt[, prob_ensemble := rowMeans(prob_mat, na.rm = TRUE)]

  ens_metrics <- compute_metrics(y_test, pred_dt$prob_ensemble)
  ens_metrics[, model := "ensemble"]
  ens_metrics[, fold := fold_ds]
  all_metrics <- c(all_metrics, list(ens_metrics))
  cat("  [E] Ensemble AUROC=", round(ens_metrics$auroc, 3), "\n\n")

  all_preds <- c(all_preds, list(pred_dt))
}

# ============================================================
# Aggregate Results
# ============================================================
cat("=== Aggregating Results ===\n\n")

# Per-sample predictions
pred_results <- rbindlist(all_preds, fill = TRUE)
pred_results[, pred_class := ifelse(prob_ensemble >= 0.5, 1, 0)]
fwrite(pred_results, file.path(OUTDIR, "tier1_loco_results.csv"))
cat("Per-sample predictions saved:", nrow(pred_results), "samples\n")

# Per-model per-fold metrics
metric_results <- rbindlist(all_metrics, fill = TRUE)
fwrite(metric_results, file.path(OUTDIR, "tier1_model_summary.csv"))

# Summary across folds
cat("\nMean LOCO AUROC by model:\n")
model_summary <- metric_results[, .(
  mean_auroc = mean(auroc, na.rm = TRUE),
  sd_auroc   = sd(auroc, na.rm = TRUE),
  mean_f1    = mean(f1, na.rm = TRUE),
  mean_acc   = mean(accuracy, na.rm = TRUE),
  mean_sens  = mean(sensitivity, na.rm = TRUE),
  mean_spec  = mean(specificity, na.rm = TRUE),
  n_folds    = .N
), by = model]
print(model_summary)

# ============================================================
# Feature Importance (Elastic Net Coefficients)
# ============================================================
cat("\n=== Feature Importance ===\n")

# Retrain elastic net on ALL LOCO data to extract coefficients
cat("Retraining elastic net on all LOCO data...\n")
w_all <- ifelse(y_loco == 1,
                length(y_loco) / (2 * sum(y_loco == 1)),
                length(y_loco) / (2 * sum(y_loco == 0)))

enet_all <- cv.glmnet(X_loco, y_loco, family = "binomial", alpha = 0.5,
                      weights = w_all, nfolds = 10, type.measure = "auc")

coefs <- as.matrix(coef(enet_all, s = "lambda.min"))
coef_dt <- data.table(
  gene = rownames(coefs),
  coefficient = as.numeric(coefs)
)
coef_dt <- coef_dt[gene != "(Intercept)"]
coef_dt[, abs_coef := abs(coefficient)]
coef_dt <- coef_dt[order(-abs_coef)]

# Top 50
top_features <- head(coef_dt, 50)
fwrite(top_features, file.path(OUTDIR, "tier1_top_features.csv"))
cat("Top 50 features saved. Top 10:\n")
print(head(top_features, 10))

# Non-zero coefficients
n_nonzero <- sum(coef_dt$abs_coef > 0)
cat("\nNon-zero coefficients:", n_nonzero, "of", nrow(coef_dt), "genes\n")

# ============================================================
# Save Best Model for Downstream Use
# ============================================================
cat("\n=== Saving Best Model ===\n")

# Identify best model by mean LOCO AUROC
best_model_name <- model_summary[which.max(mean_auroc), model]
cat("Best model:", best_model_name, "(mean AUROC =",
    round(model_summary[model == best_model_name, mean_auroc], 4), ")\n")

# Retrain best model on ALL data (including non-LOCO samples)
# Use all samples with known is_disease
all_mask <- meta$is_disease %in% c(0, 1)
X_all <- X_full[all_mask, ]
y_all <- meta$is_disease[all_mask]
w_all_full <- ifelse(y_all == 1,
                     length(y_all) / (2 * sum(y_all == 1)),
                     length(y_all) / (2 * sum(y_all == 0)))

best_model_obj <- list(
  model_type = best_model_name,
  loco_summary = model_summary,
  gene_names = colnames(X_all)
)

if (best_model_name == "elastic_net") {
  final_fit <- cv.glmnet(X_all, y_all, family = "binomial", alpha = 0.5,
                          weights = w_all_full, nfolds = 10, type.measure = "auc")
  best_model_obj$fit <- final_fit
  best_model_obj$lambda <- final_fit$lambda.min
} else if (best_model_name == "random_forest") {
  final_fit <- ranger(
    y = factor(y_all, levels = c(0, 1)),
    x = as.data.frame(X_all),
    num.trees = 500,
    mtry = floor(sqrt(ncol(X_all))),
    probability = TRUE,
    case.weights = w_all_full,
    num.threads = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4")),
    importance = "impurity",
    seed = 42
  )
  best_model_obj$fit <- final_fit
} else if (best_model_name == "xgboost") {
  y_fac <- factor(ifelse(y_all == 1, "Disease", "Control"),
                  levels = c("Control", "Disease"))
  ctrl <- trainControl(method = "cv", number = 5,
                       classProbs = TRUE, summaryFunction = twoClassSummary,
                       verboseIter = FALSE)
  xgb_grid <- expand.grid(
    nrounds = 200, max_depth = 6, eta = 0.1, gamma = 0,
    colsample_bytree = 0.8, min_child_weight = 1, subsample = 0.8
  )
  final_fit <- train(x = as.data.frame(X_all), y = y_fac,
                     method = "xgbTree", trControl = ctrl, tuneGrid = xgb_grid,
                     metric = "ROC", weights = w_all_full, verbosity = 0)
  best_model_obj$fit <- final_fit
} else if (best_model_name == "ensemble") {
  # For ensemble, save all three individual models retrained on all data
  enet_final <- cv.glmnet(X_all, y_all, family = "binomial", alpha = 0.5,
                           weights = w_all_full, nfolds = 10, type.measure = "auc")
  rf_final <- ranger(
    y = factor(y_all, levels = c(0, 1)),
    x = as.data.frame(X_all),
    num.trees = 500,
    mtry = floor(sqrt(ncol(X_all))),
    probability = TRUE,
    case.weights = w_all_full,
    num.threads = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4")),
    seed = 42
  )
  y_fac <- factor(ifelse(y_all == 1, "Disease", "Control"),
                  levels = c("Control", "Disease"))
  ctrl <- trainControl(method = "cv", number = 5,
                       classProbs = TRUE, summaryFunction = twoClassSummary,
                       verboseIter = FALSE)
  xgb_grid <- expand.grid(
    nrounds = 200, max_depth = 6, eta = 0.1, gamma = 0,
    colsample_bytree = 0.8, min_child_weight = 1, subsample = 0.8
  )
  xgb_final <- train(x = as.data.frame(X_all), y = y_fac,
                     method = "xgbTree", trControl = ctrl, tuneGrid = xgb_grid,
                     metric = "ROC", weights = w_all_full, verbosity = 0)
  best_model_obj$fit_enet <- enet_final
  best_model_obj$fit_rf <- rf_final
  best_model_obj$fit_xgb <- xgb_final
  best_model_obj$lambda <- enet_final$lambda.min
}

saveRDS(best_model_obj, file.path(OUTDIR, "tier1_best_model.rds"))
cat("Best model saved to tier1_best_model.rds\n")

# ============================================================
# Summary
# ============================================================
cat("\n=== FINAL SUMMARY ===\n")
cat("Total LOCO samples:", nrow(pred_results), "\n")
cat("Folds:", length(unique(pred_results$dataset)), "\n")
cat("Best model:", best_model_name, "\n")
cat("Mean LOCO AUROC:", round(model_summary[model == best_model_name, mean_auroc], 4), "\n")
cat("Mean LOCO F1:", round(model_summary[model == best_model_name, mean_f1], 4), "\n")
cat("Non-zero enet features:", n_nonzero, "\n")

# Overall confusion from ensemble predictions
cm <- pred_results[, .(
  accuracy = mean(y_true == pred_class, na.rm = TRUE),
  sensitivity = sum(y_true == 1 & pred_class == 1, na.rm = TRUE) /
                sum(y_true == 1, na.rm = TRUE),
  specificity = sum(y_true == 0 & pred_class == 0, na.rm = TRUE) /
                sum(y_true == 0, na.rm = TRUE)
)]
cat("\nOverall LOCO ensemble confusion:\n")
print(cm)

cat("\n=== 65_tier1_disease_detection.R completed:", as.character(Sys.time()), "===\n")
