#!/usr/bin/env Rscript
# 63_ordinal_models.R
# Ordinal elastic net + binary threshold models with LOCO cross-validation
#
# Plan 1 (Ordinal Evidence Convergence) — Phase 2, Script 1
#
# Models:
#   1. Ordinal elastic net for fibrosis (F0-F4), 6-fold LOCO
#   2. Ordinal elastic net for NAS (4 groups), 5-fold LOCO
#   3. Binary logistic elastic net for F>=3 vs F<3, 6-fold LOCO
#   4. Binary logistic elastic net for NAS>=5 vs NAS<5, 5-fold LOCO
#   5. Continuous regression (fibrosis as numeric), 6-fold LOCO
#
# Inputs:
#   - rank_expression_matrix.rds (genes x samples, top 3K)
#   - modeling_metadata.csv (per-sample metadata with folds)
#   - feature_candidates_3000.csv
#
# Outputs:
#   - ordinal_loco_fibrosis.csv (per-sample predictions)
#   - ordinal_loco_nas.csv
#   - binary_threshold_results.csv
#   - selected_features_per_fold.csv
#   - ordinal_model_summary.csv (aggregated metrics)
#   - continuous_regression_results.csv
#
# Usage: Rscript 63_ordinal_models.R
# SLURM: cpu, 8 CPUs, 64GB RAM, 48h

suppressPackageStartupMessages({
  library(data.table)
  library(glmnet)
  library(caret)
  library(pROC)
})

set.seed(42)

# =============================================================================
# PATHS
# =============================================================================
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR <- file.path(INT, "results")
OUTDIR <- file.path(RDIR, "staging_classifier")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

cat("=== 63: Ordinal Models with LOCO CV ===\n")
cat("Started:", as.character(Sys.time()), "\n\n")

# =============================================================================
# CHECK FOR ordinalNet
# =============================================================================
has_ordinalNet <- tryCatch({
  suppressPackageStartupMessages(library(ordinalNet))
  TRUE
}, error = function(e) {
  cat("ordinalNet not available; falling back to glmnet(family='multinomial')\n")
  FALSE
})

# =============================================================================
# LOAD DATA
# =============================================================================
cat("Loading rank expression matrix...\n")
rank_mat <- readRDS(file.path(OUTDIR, "rank_expression_matrix.rds"))
cat("  Dimensions:", nrow(rank_mat), "genes x", ncol(rank_mat), "samples\n")

cat("Loading modeling metadata...\n")
meta <- fread(file.path(OUTDIR, "modeling_metadata.csv"))
cat("  Metadata rows:", nrow(meta), "\n")

cat("Loading global feature candidates (fallback)...\n")
feat <- fread(file.path(OUTDIR, "feature_candidates_3000.csv"))
cat("  Feature candidates:", nrow(feat), "\n")
feature_genes_global <- feat$gene
GENE_LIST_DIR <- file.path(OUTDIR, "loco_gene_lists")

# Align samples: expression columns must match metadata rows
common_samples <- intersect(colnames(rank_mat), meta$sample_id)
cat("Common samples (expr & metadata):", length(common_samples), "\n")

# Reorder both
meta <- meta[match(common_samples, meta$sample_id), ]
X_full <- t(rank_mat[, common_samples])  # samples x genes
cat("Design matrix X:", nrow(X_full), "samples x", ncol(X_full), "genes\n\n")

stopifnot(nrow(X_full) == nrow(meta))
stopifnot(all(rownames(X_full) == meta$sample_id))

# =============================================================================
# HELPER: Compute class weights (inverse frequency)
# =============================================================================
compute_class_weights <- function(y) {
  tab <- table(y)
  n <- length(y)
  k <- length(tab)
  w <- n / (k * tab)
  sample_weights <- as.numeric(w[as.character(y)])
  return(sample_weights)
}

# =============================================================================
# HELPER: Compute ordinal metrics
# =============================================================================
compute_ordinal_metrics <- function(y_true, y_pred, y_prob = NULL, n_classes = NULL) {
  if (is.null(n_classes)) n_classes <- length(unique(y_true))

  # Accuracy
  acc <- mean(y_true == y_pred)

  # Adjacent accuracy (+/- 1)
  adj_acc <- mean(abs(as.numeric(y_true) - as.numeric(y_pred)) <= 1)

  # MAE
  mae <- mean(abs(as.numeric(y_true) - as.numeric(y_pred)))

  # Quadratic Weighted Kappa (using caret if available, else manual)
  qwk <- tryCatch({
    # Manual QWK implementation
    compute_qwk(as.numeric(y_true), as.numeric(y_pred), n_classes)
  }, error = function(e) NA_real_)

  return(data.table(
    accuracy = acc,
    adjacent_accuracy = adj_acc,
    mae = mae,
    qwk = qwk,
    n_test = length(y_true)
  ))
}

# =============================================================================
# HELPER: Quadratic Weighted Kappa (manual)
# =============================================================================
compute_qwk <- function(y_true, y_pred, n_classes) {
  min_val <- min(c(y_true, y_pred))
  max_val <- max(c(y_true, y_pred))
  labels <- min_val:max_val
  n <- length(labels)

  # Observed confusion matrix
  O <- matrix(0, nrow = n, ncol = n)
  for (i in seq_along(y_true)) {
    r <- which(labels == y_true[i])
    c <- which(labels == y_pred[i])
    if (length(r) == 1 && length(c) == 1) O[r, c] <- O[r, c] + 1
  }

  # Weight matrix (quadratic)
  W <- outer(seq_len(n), seq_len(n), function(i, j) (i - j)^2 / (n - 1)^2)

  # Expected confusion under independence
  row_sums <- rowSums(O)
  col_sums <- colSums(O)
  total <- sum(O)
  if (total == 0) return(NA_real_)

  E <- outer(row_sums, col_sums) / total

  # Kappa
  num <- sum(W * O)
  den <- sum(W * E)
  if (den == 0) return(1.0)

  return(1 - num / den)
}

# =============================================================================
# HELPER: Compute binary metrics
# =============================================================================
compute_binary_metrics <- function(y_true, y_pred, y_prob = NULL) {
  cm <- table(factor(y_true, levels = c(0, 1)), factor(y_pred, levels = c(0, 1)))
  tn <- cm[1, 1]; fp <- cm[1, 2]; fn <- cm[2, 1]; tp <- cm[2, 2]

  acc <- (tp + tn) / sum(cm)
  sens <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
  spec <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
  f1 <- if ((tp + fp + fn) > 0) 2 * tp / (2 * tp + fp + fn) else NA_real_

  auroc <- NA_real_
  if (!is.null(y_prob) && length(unique(y_true)) == 2) {
    auroc <- tryCatch({
      as.numeric(pROC::auc(pROC::roc(y_true, y_prob, quiet = TRUE)))
    }, error = function(e) NA_real_)
  }

  return(data.table(
    accuracy = acc,
    sensitivity = sens,
    specificity = spec,
    f1 = f1,
    auroc = auroc,
    n_test = length(y_true),
    n_pos = sum(y_true == 1),
    n_neg = sum(y_true == 0)
  ))
}

# =============================================================================
# HELPER: Run glmnet ordinal/multinomial LOCO for one fold
# =============================================================================
run_glmnet_ordinal_fold <- function(X_train, y_train, X_test, y_test,
                                     fold_name, target_type, alpha = 0.5) {
  cat("    Fold:", fold_name, " | Train:", nrow(X_train),
      "| Test:", nrow(X_test), "\n")

  # Class weights
  w_train <- compute_class_weights(y_train)

  levels_y <- sort(unique(c(y_train, y_test)))
  y_train_f <- factor(y_train, levels = levels_y)
  y_test_f  <- factor(y_test, levels = levels_y)

  # Attempt ordinalNet first, then fall back to glmnet multinomial
  model_type <- "glmnet_multinomial"
  fit <- NULL
  pred_class <- NULL
  pred_prob <- NULL
  coef_list <- NULL

  if (has_ordinalNet && target_type == "ordinal") {
    fit <- tryCatch({
      cat("      Trying ordinalNet...\n")
      ordinalNet(
        x = X_train,
        y = y_train_f,
        family = "cumulative",
        link = "logit",
        alpha = alpha,
        nFolds = 5,
        printProgress = FALSE,
        lambdaMinRatio = 0.01
      )
    }, error = function(e) {
      cat("      ordinalNet failed:", conditionMessage(e), "\n")
      NULL
    })

    if (!is.null(fit)) {
      model_type <- "ordinalNet"
      pred_prob_mat <- predict(fit, newx = X_test, type = "response")
      pred_class <- levels_y[apply(pred_prob_mat, 1, which.max)]
      pred_prob <- pred_prob_mat

      # Extract non-zero coefficients
      co <- coef(fit)
      if (is.matrix(co)) {
        nonzero <- which(co[-seq_len(length(levels_y) - 1), ] != 0, arr.ind = TRUE)
        if (nrow(nonzero) > 0) {
          gene_idx <- nonzero[, 1]
          coef_list <- data.table(
            gene = colnames(X_train)[gene_idx],
            coef_index = gene_idx
          )
        }
      }
    }
  }

  # Fallback: glmnet multinomial
  if (is.null(fit)) {
    cat("      Using glmnet multinomial (alpha=", alpha, ")...\n")
    fit <- tryCatch({
      cv.glmnet(
        x = X_train,
        y = y_train_f,
        family = "multinomial",
        type.multinomial = "grouped",
        alpha = alpha,
        weights = w_train,
        nfolds = 5,
        type.measure = "class",
        parallel = FALSE
      )
    }, error = function(e) {
      cat("      glmnet multinomial failed:", conditionMessage(e), "\n")
      NULL
    })

    if (!is.null(fit)) {
      model_type <- "glmnet_multinomial"
      pred_class_f <- predict(fit, newx = X_test, s = "lambda.min", type = "class")
      pred_class <- as.vector(pred_class_f)
      pred_prob_raw <- predict(fit, newx = X_test, s = "lambda.min", type = "response")
      # predict returns 3D array (samples x classes x 1), squeeze
      pred_prob <- pred_prob_raw[, , 1, drop = FALSE]
      dim(pred_prob) <- dim(pred_prob)[1:2]
      colnames(pred_prob) <- levels_y

      # Extract non-zero coefficients
      co <- coef(fit, s = "lambda.min")
      all_genes <- character(0)
      for (cl in names(co)) {
        nz <- which(co[[cl]][-1, 1] != 0)  # skip intercept
        all_genes <- union(all_genes, colnames(X_train)[nz])
      }
      if (length(all_genes) > 0) {
        coef_list <- data.table(gene = all_genes)
      }
    }
  }

  # Handle total failure
  if (is.null(fit)) {
    cat("      WARNING: All models failed for fold", fold_name, "\n")
    return(list(predictions = NULL, features = NULL, metrics = NULL))
  }

  # Convert pred_class to numeric for metric computation
  pred_numeric <- as.numeric(as.character(pred_class))
  y_numeric    <- as.numeric(as.character(y_test))

  # Metrics
  metrics <- compute_ordinal_metrics(y_numeric, pred_numeric,
                                      n_classes = length(levels_y))
  metrics[, fold := fold_name]
  metrics[, model_type := model_type]
  metrics[, n_features := if (!is.null(coef_list)) nrow(coef_list) else 0L]

  cat("      Acc:", round(metrics$accuracy, 3),
      "| QWK:", round(metrics$qwk, 3),
      "| MAE:", round(metrics$mae, 3),
      "| Adj:", round(metrics$adjacent_accuracy, 3),
      "| Feats:", metrics$n_features, "\n")

  # Per-sample predictions
  preds <- data.table(
    sample_id = rownames(X_test),
    fold = fold_name,
    y_true = y_numeric,
    y_pred = pred_numeric,
    model_type = model_type
  )

  # Attach class probabilities
  if (!is.null(pred_prob)) {
    prob_dt <- as.data.table(pred_prob)
    setnames(prob_dt, paste0("prob_", colnames(pred_prob)))
    preds <- cbind(preds, prob_dt)
  }

  # Features
  feats <- NULL
  if (!is.null(coef_list)) {
    feats <- copy(coef_list)
    feats[, fold := fold_name]
  }

  return(list(predictions = preds, features = feats, metrics = metrics))
}

# =============================================================================
# HELPER: Run binary logistic elastic net LOCO for one fold
# =============================================================================
run_glmnet_binary_fold <- function(X_train, y_train, X_test, y_test,
                                    fold_name, alpha = 0.5) {
  cat("    Fold:", fold_name, " | Train:", nrow(X_train),
      "| Test:", nrow(X_test),
      "| Pos:", sum(y_test == 1), "| Neg:", sum(y_test == 0), "\n")

  w_train <- compute_class_weights(y_train)

  fit <- tryCatch({
    cv.glmnet(
      x = X_train,
      y = y_train,
      family = "binomial",
      alpha = alpha,
      weights = w_train,
      nfolds = 5,
      type.measure = "auc",
      parallel = FALSE
    )
  }, error = function(e) {
    cat("      cv.glmnet binomial failed:", conditionMessage(e), "\n")
    NULL
  })

  if (is.null(fit)) {
    return(list(predictions = NULL, features = NULL, metrics = NULL))
  }

  pred_prob <- as.vector(predict(fit, newx = X_test, s = "lambda.min",
                                  type = "response"))
  pred_class <- as.integer(pred_prob >= 0.5)

  # Non-zero coefficients
  co <- coef(fit, s = "lambda.min")
  nz <- which(co[-1, 1] != 0)
  feats <- NULL
  if (length(nz) > 0) {
    feats <- data.table(
      gene = colnames(X_train)[nz],
      coef = co[-1, 1][nz],
      fold = fold_name
    )
  }

  metrics <- compute_binary_metrics(y_test, pred_class, pred_prob)
  metrics[, fold := fold_name]
  metrics[, n_features := length(nz)]

  cat("      Acc:", round(metrics$accuracy, 3),
      "| AUROC:", round(metrics$auroc, 3),
      "| Sens:", round(metrics$sensitivity, 3),
      "| Spec:", round(metrics$specificity, 3),
      "| Feats:", metrics$n_features, "\n")

  preds <- data.table(
    sample_id = rownames(X_test),
    fold = fold_name,
    y_true = y_test,
    y_pred = pred_class,
    y_prob = pred_prob
  )

  return(list(predictions = preds, features = feats, metrics = metrics))
}

# =============================================================================
# HELPER: Load per-fold gene list (falls back to global list)
# =============================================================================
load_fold_genes <- function(fold_id) {
  fold_gene_file <- file.path(GENE_LIST_DIR,
                               paste0("fold_", fold_id, "_feature_candidates_3000.csv"))
  if (file.exists(fold_gene_file)) {
    cat("    Using per-fold gene list for fold:", fold_id, "\n")
    fread(fold_gene_file)$gene
  } else {
    cat("    Fold gene list not found for", fold_id, "-- using global fallback\n")
    feature_genes_global
  }
}

# Subset X to fold-specific genes (only columns present in X)
subset_X_by_genes <- function(X, fold_genes) {
  keep <- intersect(fold_genes, colnames(X))
  X[, keep, drop = FALSE]
}

# =============================================================================
# MODEL 1: Ordinal elastic net for FIBROSIS (F0-F4), 6-fold LOCO
# =============================================================================
cat("\n========================================\n")
cat("MODEL 1: Ordinal Fibrosis (F0-F4)\n")
cat("========================================\n\n")

fib_mask <- meta$fib_stage >= 0 & meta$loco_fold_fibrosis != "excluded"
fib_idx  <- which(fib_mask)
fib_folds <- unique(meta$loco_fold_fibrosis[fib_idx])
fib_folds <- fib_folds[fib_folds != "excluded"]
cat("Fibrosis samples:", length(fib_idx), "| Folds:", length(fib_folds), "\n")
cat("Class distribution:\n")
print(table(meta$fib_stage[fib_idx]))

X_fib <- X_full[fib_idx, , drop = FALSE]
y_fib <- meta$fib_stage[fib_idx]
fold_fib <- meta$loco_fold_fibrosis[fib_idx]

fib_preds_list   <- list()
fib_feats_list   <- list()
fib_metrics_list <- list()

for (f in sort(fib_folds)) {
  test_idx  <- which(fold_fib == f)
  train_idx <- which(fold_fib != f)

  if (length(test_idx) < 5 || length(train_idx) < 20) {
    cat("  Skipping fold", f, "(too few samples)\n")
    next
  }

  # Check that training has at least 2 classes
  if (length(unique(y_fib[train_idx])) < 2) {
    cat("  Skipping fold", f, "(< 2 classes in training)\n")
    next
  }

  # Load per-fold gene list and subset features (removes leakage)
  fold_genes <- load_fold_genes(f)
  X_fib_fold <- subset_X_by_genes(X_fib, fold_genes)

  result <- run_glmnet_ordinal_fold(
    X_train = X_fib_fold[train_idx, , drop = FALSE],
    y_train = y_fib[train_idx],
    X_test  = X_fib_fold[test_idx, , drop = FALSE],
    y_test  = y_fib[test_idx],
    fold_name = f,
    target_type = "ordinal",
    alpha = 0.5
  )

  if (!is.null(result$predictions)) fib_preds_list[[f]]   <- result$predictions
  if (!is.null(result$features))    fib_feats_list[[f]]   <- result$features
  if (!is.null(result$metrics))     fib_metrics_list[[f]] <- result$metrics
}

fib_preds <- rbindlist(fib_preds_list, fill = TRUE)
fib_preds[, target := "fibrosis_ordinal"]
fib_feats <- rbindlist(fib_feats_list, fill = TRUE)
fib_feats[, target := "fibrosis_ordinal"]
fib_metrics <- rbindlist(fib_metrics_list, fill = TRUE)
fib_metrics[, target := "fibrosis_ordinal"]

cat("\nFibrosis ordinal aggregate:\n")
cat("  Mean accuracy:", round(mean(fib_metrics$accuracy, na.rm = TRUE), 3), "\n")
cat("  Mean QWK:", round(mean(fib_metrics$qwk, na.rm = TRUE), 3), "\n")
cat("  Mean MAE:", round(mean(fib_metrics$mae, na.rm = TRUE), 3), "\n")
cat("  Mean adj acc:", round(mean(fib_metrics$adjacent_accuracy, na.rm = TRUE), 3), "\n")

# =============================================================================
# MODEL 2: Ordinal elastic net for NAS (4 groups), 5-fold LOCO
# =============================================================================
cat("\n========================================\n")
cat("MODEL 2: Ordinal NAS (4 groups)\n")
cat("========================================\n\n")

nas_mask <- meta$nas_group4 >= 0 & meta$loco_fold_nas != "excluded"
nas_idx  <- which(nas_mask)
nas_folds <- unique(meta$loco_fold_nas[nas_idx])
nas_folds <- nas_folds[nas_folds != "excluded"]
cat("NAS samples:", length(nas_idx), "| Folds:", length(nas_folds), "\n")
cat("NAS group4 distribution:\n")
print(table(meta$nas_group4[nas_idx]))

X_nas <- X_full[nas_idx, , drop = FALSE]
y_nas <- meta$nas_group4[nas_idx]
fold_nas <- meta$loco_fold_nas[nas_idx]

nas_preds_list   <- list()
nas_feats_list   <- list()
nas_metrics_list <- list()

for (f in sort(nas_folds)) {
  test_idx  <- which(fold_nas == f)
  train_idx <- which(fold_nas != f)

  if (length(test_idx) < 5 || length(train_idx) < 20) {
    cat("  Skipping fold", f, "(too few samples)\n")
    next
  }

  if (length(unique(y_nas[train_idx])) < 2) {
    cat("  Skipping fold", f, "(< 2 classes in training)\n")
    next
  }

  # Load per-fold gene list and subset features (removes leakage)
  fold_genes <- load_fold_genes(f)
  X_nas_fold <- subset_X_by_genes(X_nas, fold_genes)

  result <- run_glmnet_ordinal_fold(
    X_train = X_nas_fold[train_idx, , drop = FALSE],
    y_train = y_nas[train_idx],
    X_test  = X_nas_fold[test_idx, , drop = FALSE],
    y_test  = y_nas[test_idx],
    fold_name = f,
    target_type = "ordinal",
    alpha = 0.5
  )

  if (!is.null(result$predictions)) nas_preds_list[[f]]   <- result$predictions
  if (!is.null(result$features))    nas_feats_list[[f]]   <- result$features
  if (!is.null(result$metrics))     nas_metrics_list[[f]] <- result$metrics
}

nas_preds <- rbindlist(nas_preds_list, fill = TRUE)
nas_preds[, target := "nas_ordinal"]
nas_feats <- rbindlist(nas_feats_list, fill = TRUE)
nas_feats[, target := "nas_ordinal"]
nas_metrics <- rbindlist(nas_metrics_list, fill = TRUE)
nas_metrics[, target := "nas_ordinal"]

cat("\nNAS ordinal aggregate:\n")
cat("  Mean accuracy:", round(mean(nas_metrics$accuracy, na.rm = TRUE), 3), "\n")
cat("  Mean QWK:", round(mean(nas_metrics$qwk, na.rm = TRUE), 3), "\n")
cat("  Mean MAE:", round(mean(nas_metrics$mae, na.rm = TRUE), 3), "\n")
cat("  Mean adj acc:", round(mean(nas_metrics$adjacent_accuracy, na.rm = TRUE), 3), "\n")

# =============================================================================
# MODEL 3: Binary F>=3 vs F<3, 6-fold LOCO
# =============================================================================
cat("\n========================================\n")
cat("MODEL 3: Binary F>=3 (Significant Fibrosis)\n")
cat("========================================\n\n")

fib_bin_mask <- meta$fib_ge3 >= 0 & meta$loco_fold_fibrosis != "excluded"
fib_bin_idx  <- which(fib_bin_mask)
cat("F>=3 binary samples:", length(fib_bin_idx), "\n")
cat("Class distribution:\n")
print(table(meta$fib_ge3[fib_bin_idx]))

X_fbin <- X_full[fib_bin_idx, , drop = FALSE]
y_fbin <- meta$fib_ge3[fib_bin_idx]
fold_fbin <- meta$loco_fold_fibrosis[fib_bin_idx]

fbin_preds_list   <- list()
fbin_feats_list   <- list()
fbin_metrics_list <- list()

for (f in sort(fib_folds)) {
  test_idx  <- which(fold_fbin == f)
  train_idx <- which(fold_fbin != f)

  if (length(test_idx) < 5 || length(train_idx) < 20) {
    cat("  Skipping fold", f, "\n")
    next
  }

  if (length(unique(y_fbin[train_idx])) < 2) {
    cat("  Skipping fold", f, "(single class in training)\n")
    next
  }

  # Load per-fold gene list and subset features (removes leakage)
  fold_genes <- load_fold_genes(f)
  X_fbin_fold <- subset_X_by_genes(X_fbin, fold_genes)

  result <- run_glmnet_binary_fold(
    X_train = X_fbin_fold[train_idx, , drop = FALSE],
    y_train = y_fbin[train_idx],
    X_test  = X_fbin_fold[test_idx, , drop = FALSE],
    y_test  = y_fbin[test_idx],
    fold_name = f,
    alpha = 0.5
  )

  if (!is.null(result$predictions)) fbin_preds_list[[f]]   <- result$predictions
  if (!is.null(result$features))    fbin_feats_list[[f]]   <- result$features
  if (!is.null(result$metrics))     fbin_metrics_list[[f]] <- result$metrics
}

fbin_preds <- rbindlist(fbin_preds_list, fill = TRUE)
fbin_preds[, target := "fib_ge3_binary"]
fbin_feats <- rbindlist(fbin_feats_list, fill = TRUE)
fbin_feats[, target := "fib_ge3_binary"]
fbin_metrics <- rbindlist(fbin_metrics_list, fill = TRUE)
fbin_metrics[, target := "fib_ge3_binary"]

cat("\nBinary F>=3 aggregate:\n")
cat("  Mean accuracy:", round(mean(fbin_metrics$accuracy, na.rm = TRUE), 3), "\n")
cat("  Mean AUROC:", round(mean(fbin_metrics$auroc, na.rm = TRUE), 3), "\n")
cat("  Mean sensitivity:", round(mean(fbin_metrics$sensitivity, na.rm = TRUE), 3), "\n")
cat("  Mean specificity:", round(mean(fbin_metrics$specificity, na.rm = TRUE), 3), "\n")

# =============================================================================
# MODEL 4: Binary NAS>=5 vs NAS<5, 5-fold LOCO
# =============================================================================
cat("\n========================================\n")
cat("MODEL 4: Binary NAS>=5 (Definite NASH)\n")
cat("========================================\n\n")

nas_bin_mask <- meta$nas_ge5 >= 0 & meta$loco_fold_nas != "excluded"
nas_bin_idx  <- which(nas_bin_mask)
cat("NAS>=5 binary samples:", length(nas_bin_idx), "\n")
cat("Class distribution:\n")
print(table(meta$nas_ge5[nas_bin_idx]))

X_nbin <- X_full[nas_bin_idx, , drop = FALSE]
y_nbin <- meta$nas_ge5[nas_bin_idx]
fold_nbin <- meta$loco_fold_nas[nas_bin_idx]

nbin_preds_list   <- list()
nbin_feats_list   <- list()
nbin_metrics_list <- list()

for (f in sort(nas_folds)) {
  test_idx  <- which(fold_nbin == f)
  train_idx <- which(fold_nbin != f)

  if (length(test_idx) < 5 || length(train_idx) < 20) {
    cat("  Skipping fold", f, "\n")
    next
  }

  if (length(unique(y_nbin[train_idx])) < 2) {
    cat("  Skipping fold", f, "(single class in training)\n")
    next
  }

  # Load per-fold gene list and subset features (removes leakage)
  fold_genes <- load_fold_genes(f)
  X_nbin_fold <- subset_X_by_genes(X_nbin, fold_genes)

  result <- run_glmnet_binary_fold(
    X_train = X_nbin_fold[train_idx, , drop = FALSE],
    y_train = y_nbin[train_idx],
    X_test  = X_nbin_fold[test_idx, , drop = FALSE],
    y_test  = y_nbin[test_idx],
    fold_name = f,
    alpha = 0.5
  )

  if (!is.null(result$predictions)) nbin_preds_list[[f]]   <- result$predictions
  if (!is.null(result$features))    nbin_feats_list[[f]]   <- result$features
  if (!is.null(result$metrics))     nbin_metrics_list[[f]] <- result$metrics
}

nbin_preds <- rbindlist(nbin_preds_list, fill = TRUE)
nbin_preds[, target := "nas_ge5_binary"]
nbin_feats <- rbindlist(nbin_feats_list, fill = TRUE)
nbin_feats[, target := "nas_ge5_binary"]
nbin_metrics <- rbindlist(nbin_metrics_list, fill = TRUE)
nbin_metrics[, target := "nas_ge5_binary"]

cat("\nBinary NAS>=5 aggregate:\n")
cat("  Mean accuracy:", round(mean(nbin_metrics$accuracy, na.rm = TRUE), 3), "\n")
cat("  Mean AUROC:", round(mean(nbin_metrics$auroc, na.rm = TRUE), 3), "\n")
cat("  Mean sensitivity:", round(mean(nbin_metrics$sensitivity, na.rm = TRUE), 3), "\n")
cat("  Mean specificity:", round(mean(nbin_metrics$specificity, na.rm = TRUE), 3), "\n")

# =============================================================================
# MODEL 5: Continuous regression (fibrosis stage as numeric), 6-fold LOCO
# =============================================================================
cat("\n========================================\n")
cat("MODEL 5: Continuous Regression (Fibrosis)\n")
cat("========================================\n\n")

cont_preds_list   <- list()
cont_feats_list   <- list()
cont_metrics_list <- list()

for (f in sort(fib_folds)) {
  test_idx  <- which(fold_fib == f)
  train_idx <- which(fold_fib != f)

  if (length(test_idx) < 5 || length(train_idx) < 20) {
    cat("  Skipping fold", f, "\n")
    next
  }

  cat("    Fold:", f, " | Train:", length(train_idx),
      "| Test:", length(test_idx), "\n")

  # Load per-fold gene list and subset features (removes leakage)
  fold_genes <- load_fold_genes(f)
  X_fib_fold <- subset_X_by_genes(X_fib, fold_genes)

  w_train <- compute_class_weights(y_fib[train_idx])

  fit <- tryCatch({
    cv.glmnet(
      x = X_fib_fold[train_idx, , drop = FALSE],
      y = y_fib[train_idx],
      family = "gaussian",
      alpha = 0.5,
      weights = w_train,
      nfolds = 5,
      type.measure = "mse",
      parallel = FALSE
    )
  }, error = function(e) {
    cat("      cv.glmnet gaussian failed:", conditionMessage(e), "\n")
    NULL
  })

  if (is.null(fit)) next

  pred_cont <- as.vector(predict(fit, newx = X_fib_fold[test_idx, , drop = FALSE],
                                  s = "lambda.min"))
  pred_rounded <- pmin(pmax(round(pred_cont), 0), 4)

  # Metrics
  mse <- mean((y_fib[test_idx] - pred_cont)^2)
  rmse <- sqrt(mse)
  mae_cont <- mean(abs(y_fib[test_idx] - pred_cont))
  r2 <- 1 - mse / var(y_fib[test_idx])

  # Ordinal metrics after rounding
  ord_met <- compute_ordinal_metrics(y_fib[test_idx], pred_rounded)

  metrics <- data.table(
    fold = f,
    rmse = rmse,
    mae_continuous = mae_cont,
    r_squared = r2,
    accuracy_rounded = ord_met$accuracy,
    qwk_rounded = ord_met$qwk,
    adj_acc_rounded = ord_met$adjacent_accuracy,
    n_test = length(test_idx)
  )

  # Non-zero features
  co <- coef(fit, s = "lambda.min")
  nz <- which(co[-1, 1] != 0)
  metrics[, n_features := length(nz)]

  feats <- NULL
  if (length(nz) > 0) {
    feats <- data.table(
      gene = colnames(X_fib_fold)[nz],
      coef = co[-1, 1][nz],
      fold = f,
      target = "fibrosis_continuous"
    )
  }

  cat("      RMSE:", round(rmse, 3),
      "| R2:", round(r2, 3),
      "| QWK (rounded):", round(ord_met$qwk, 3),
      "| Feats:", length(nz), "\n")

  cont_preds_list[[f]] <- data.table(
    sample_id = rownames(X_fib)[test_idx],
    fold = f,
    y_true = y_fib[test_idx],
    y_pred_continuous = pred_cont,
    y_pred_rounded = pred_rounded,
    target = "fibrosis_continuous"
  )
  cont_feats_list[[f]] <- feats
  cont_metrics_list[[f]] <- metrics
}

cont_preds <- rbindlist(cont_preds_list, fill = TRUE)
cont_feats <- rbindlist(cont_feats_list, fill = TRUE)
cont_metrics <- rbindlist(cont_metrics_list, fill = TRUE)
cont_metrics[, target := "fibrosis_continuous"]

cat("\nContinuous fibrosis aggregate:\n")
cat("  Mean RMSE:", round(mean(cont_metrics$rmse, na.rm = TRUE), 3), "\n")
cat("  Mean R2:", round(mean(cont_metrics$r_squared, na.rm = TRUE), 3), "\n")
cat("  Mean QWK (rounded):", round(mean(cont_metrics$qwk_rounded, na.rm = TRUE), 3), "\n")

# =============================================================================
# SAVE OUTPUTS
# =============================================================================
cat("\n========================================\n")
cat("SAVING OUTPUTS\n")
cat("========================================\n\n")

# 1. Per-sample predictions — fibrosis ordinal
fwrite(fib_preds, file.path(OUTDIR, "ordinal_loco_fibrosis.csv"))
cat("Saved ordinal_loco_fibrosis.csv:", nrow(fib_preds), "rows\n")

# 2. Per-sample predictions — NAS ordinal
fwrite(nas_preds, file.path(OUTDIR, "ordinal_loco_nas.csv"))
cat("Saved ordinal_loco_nas.csv:", nrow(nas_preds), "rows\n")

# 3. Binary threshold results (combine F>=3 and NAS>=5)
binary_preds <- rbindlist(list(fbin_preds, nbin_preds), fill = TRUE)
fwrite(binary_preds, file.path(OUTDIR, "binary_threshold_results.csv"))
cat("Saved binary_threshold_results.csv:", nrow(binary_preds), "rows\n")

# 4. Continuous regression predictions
fwrite(cont_preds, file.path(OUTDIR, "continuous_regression_results.csv"))
cat("Saved continuous_regression_results.csv:", nrow(cont_preds), "rows\n")

# 5. Selected features per fold (all models)
all_feats <- rbindlist(list(fib_feats, nas_feats, fbin_feats, nbin_feats, cont_feats),
                       fill = TRUE)
fwrite(all_feats, file.path(OUTDIR, "selected_features_per_fold.csv"))
cat("Saved selected_features_per_fold.csv:", nrow(all_feats), "rows\n")

# Feature frequency across folds
if (nrow(all_feats) > 0) {
  feat_freq <- all_feats[, .(
    n_folds = .N,
    targets = paste(unique(target), collapse = ";")
  ), by = gene]
  feat_freq <- feat_freq[order(-n_folds)]
  cat("\nTop 20 most frequently selected features:\n")
  print(head(feat_freq, 20))
}

# 6. Aggregated model summary
all_ordinal_metrics <- rbindlist(list(fib_metrics, nas_metrics), fill = TRUE)
all_binary_metrics <- rbindlist(list(fbin_metrics, nbin_metrics), fill = TRUE)

# Summarize ordinal models
ord_summary <- all_ordinal_metrics[, .(
  mean_accuracy = mean(accuracy, na.rm = TRUE),
  sd_accuracy = sd(accuracy, na.rm = TRUE),
  mean_qwk = mean(qwk, na.rm = TRUE),
  sd_qwk = sd(qwk, na.rm = TRUE),
  mean_mae = mean(mae, na.rm = TRUE),
  sd_mae = sd(mae, na.rm = TRUE),
  mean_adj_acc = mean(adjacent_accuracy, na.rm = TRUE),
  sd_adj_acc = sd(adjacent_accuracy, na.rm = TRUE),
  mean_n_features = mean(n_features, na.rm = TRUE),
  n_folds = .N,
  total_test_samples = sum(n_test, na.rm = TRUE)
), by = target]

# Summarize binary models
bin_summary <- all_binary_metrics[, .(
  mean_accuracy = mean(accuracy, na.rm = TRUE),
  sd_accuracy = sd(accuracy, na.rm = TRUE),
  mean_auroc = mean(auroc, na.rm = TRUE),
  sd_auroc = sd(auroc, na.rm = TRUE),
  mean_sensitivity = mean(sensitivity, na.rm = TRUE),
  sd_sensitivity = sd(sensitivity, na.rm = TRUE),
  mean_specificity = mean(specificity, na.rm = TRUE),
  sd_specificity = sd(specificity, na.rm = TRUE),
  mean_f1 = mean(f1, na.rm = TRUE),
  sd_f1 = sd(f1, na.rm = TRUE),
  mean_n_features = mean(n_features, na.rm = TRUE),
  n_folds = .N,
  total_test_samples = sum(n_test, na.rm = TRUE)
), by = target]

# Summarize continuous
cont_summary <- cont_metrics[, .(
  target = "fibrosis_continuous",
  mean_rmse = mean(rmse, na.rm = TRUE),
  sd_rmse = sd(rmse, na.rm = TRUE),
  mean_r_squared = mean(r_squared, na.rm = TRUE),
  sd_r_squared = sd(r_squared, na.rm = TRUE),
  mean_qwk_rounded = mean(qwk_rounded, na.rm = TRUE),
  sd_qwk_rounded = sd(qwk_rounded, na.rm = TRUE),
  mean_n_features = mean(n_features, na.rm = TRUE),
  n_folds = .N,
  total_test_samples = sum(n_test, na.rm = TRUE)
)]

# Combine all summaries into one table
model_summary <- rbindlist(list(
  ord_summary[, .(target, metric_type = "ordinal",
    mean_primary = mean_qwk, sd_primary = sd_qwk,
    mean_accuracy, sd_accuracy, mean_n_features, n_folds, total_test_samples)],
  bin_summary[, .(target, metric_type = "binary",
    mean_primary = mean_auroc, sd_primary = sd_auroc,
    mean_accuracy, sd_accuracy, mean_n_features, n_folds, total_test_samples)],
  cont_summary[, .(target, metric_type = "continuous",
    mean_primary = mean_qwk_rounded, sd_primary = sd_qwk_rounded,
    mean_accuracy = NA_real_, sd_accuracy = NA_real_,
    mean_n_features, n_folds, total_test_samples)]
), fill = TRUE)

fwrite(model_summary, file.path(OUTDIR, "ordinal_model_summary.csv"))
cat("\nSaved ordinal_model_summary.csv\n")

# Also save per-fold detail tables
fwrite(all_ordinal_metrics, file.path(OUTDIR, "ordinal_per_fold_metrics.csv"))
fwrite(all_binary_metrics, file.path(OUTDIR, "binary_per_fold_metrics.csv"))
fwrite(cont_metrics, file.path(OUTDIR, "continuous_per_fold_metrics.csv"))

# =============================================================================
# FINAL SUMMARY
# =============================================================================
cat("\n========================================\n")
cat("FINAL SUMMARY\n")
cat("========================================\n")
cat("\nModel summary:\n")
print(model_summary)

cat("\n=== 63_ordinal_models.R completed:", as.character(Sys.time()), "===\n")
