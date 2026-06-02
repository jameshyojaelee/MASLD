#!/usr/bin/env Rscript
# 67_tier3_staging.R
# Tier 3 of the hierarchical staging cascade: fine-grained staging within
# disease samples only (controls removed).
#
#   Tier 3a: NAS ordinal (4 groups: Low/Moderate/High/VeryHigh from nas_group4)
#   Tier 3b: Fibrosis ordinal (F0-F4 from fib_stage)
#   Tier 3c: Joint NAS + fibrosis (multi-output XGBoost via caret)
#
# Models: multinomial elastic net (glmnet) + ordinalNet (if available) per tier.
# LOCO across annotated cohorts.
#
# Inputs:
#   - rank_expression_matrix.rds: Top 3K genes (genes x samples)
#   - modeling_metadata.csv: Per-sample metadata with nas_group4, fib_stage, etc.
#
# Outputs:
#   - tier3a_nas_results.csv: Per-sample NAS ordinal predictions
#   - tier3b_fibrosis_results.csv: Per-sample fibrosis ordinal predictions
#   - tier3c_joint_results.csv: Joint NAS + fibrosis predictions (concordance)
#   - tier3_model_summary.csv: Per-tier per-model per-fold metrics
#
# Usage: Rscript 67_tier3_staging.R
# SLURM: cpu, 8 CPUs, 64GB RAM, 48h

suppressPackageStartupMessages({
  library(data.table)
  library(glmnet)
  library(caret)
  library(pROC)
  library(ranger)
  library(Matrix)
})

# Attempt to load ordinalNet (may not be installed)
has_ordinalNet <- tryCatch({
  suppressPackageStartupMessages(library(ordinalNet))
  TRUE
}, error = function(e) FALSE)
cat("ordinalNet available:", has_ordinalNet, "\n")

set.seed(42)

# --- Paths ---
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR <- file.path(INT, "results")
OUTDIR <- file.path(RDIR, "staging_classifier")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

cat("=== 67: Tier 3 Fine-Grained Staging ===\n")
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

# --- Only disease samples ---
dis_mask <- meta$is_disease == 1
meta_dis <- meta[dis_mask]
X_dis <- X_full[dis_mask, ]
cat("Disease samples:", nrow(meta_dis), "\n\n")

# NAS-relevant datasets
NAS_DATASETS <- c("GSE130970", "GSE135251", "GSE162694", "GSE174478", "GSE193066")
# Fibrosis-relevant datasets
FIB_DATASETS <- c("GSE130970", "GSE135251", "GSE162694", "GSE174478", "GSE193066", "GSE240729")

# ============================================================
# Helper: Ordinal classification metrics
# ============================================================
compute_ordinal_metrics <- function(y_true, y_pred, prob_mat = NULL,
                                     class_labels = NULL) {
  # Accuracy
  accuracy <- mean(y_true == y_pred)

  # Quadratic Weighted Kappa (QWK)
  # Compute manually to avoid extra dependencies
  n <- length(y_true)
  if (is.null(class_labels)) class_labels <- sort(unique(c(y_true, y_pred)))
  n_classes <- length(class_labels)

  # Confusion matrix
  cm <- table(factor(y_true, levels = class_labels),
              factor(y_pred, levels = class_labels))
  cm <- as.matrix(cm)

  # Weight matrix (quadratic)
  w_mat <- matrix(0, n_classes, n_classes)
  for (i in 1:n_classes) {
    for (j in 1:n_classes) {
      w_mat[i, j] <- (i - j)^2 / (n_classes - 1)^2
    }
  }

  # Expected confusion under independence
  row_sums <- rowSums(cm)
  col_sums <- colSums(cm)
  e_mat <- outer(row_sums, col_sums) / n

  # QWK
  num <- sum(w_mat * cm)
  den <- sum(w_mat * e_mat)
  qwk <- if (den > 0) 1 - num / den else 0

  # Mean Absolute Error (ordinal distance)
  mae <- mean(abs(y_true - y_pred))

  # Adjacent accuracy: correct or off-by-one
  adj_acc <- mean(abs(y_true - y_pred) <= 1)

  # Per-class F1
  per_class <- data.table()
  for (cl in class_labels) {
    tp <- sum(y_true == cl & y_pred == cl)
    fp <- sum(y_true != cl & y_pred == cl)
    fn <- sum(y_true == cl & y_pred != cl)
    prec <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
    rec  <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
    f1   <- if (!is.na(prec) && !is.na(rec) && (prec + rec) > 0) {
      2 * prec * rec / (prec + rec)
    } else { NA_real_ }
    per_class <- rbind(per_class, data.table(
      class = cl, precision = prec, recall = rec, f1 = f1,
      n_true = sum(y_true == cl)
    ))
  }

  # Macro AUROC if probability matrix provided
  macro_auroc <- NA_real_
  if (!is.null(prob_mat)) {
    auc_vals <- numeric(n_classes)
    for (i in seq_len(n_classes)) {
      y_bin <- as.integer(y_true == class_labels[i])
      if (sum(y_bin) > 0 && sum(y_bin) < length(y_bin)) {
        roc_obj <- tryCatch(
          roc(y_bin, prob_mat[, i], levels = c(0, 1), direction = "<", quiet = TRUE),
          error = function(e) NULL
        )
        auc_vals[i] <- if (!is.null(roc_obj)) as.numeric(auc(roc_obj)) else NA_real_
      } else {
        auc_vals[i] <- NA_real_
      }
    }
    macro_auroc <- mean(auc_vals, na.rm = TRUE)
  }

  list(
    accuracy = accuracy, qwk = qwk, mae = mae, adj_acc = adj_acc,
    macro_auroc = macro_auroc, macro_f1 = mean(per_class$f1, na.rm = TRUE),
    per_class = per_class
  )
}

# ============================================================
# TIER 3a: NAS Ordinal (4 groups)
# ============================================================
cat("=== TIER 3a: NAS Ordinal Classification ===\n\n")

nas_mask <- meta_dis$nas_group4 >= 0 & meta_dis$dataset %in% NAS_DATASETS
meta_nas <- meta_dis[nas_mask]
X_nas <- X_dis[nas_mask, ]
y_nas <- meta_nas$nas_group4

NAS_LABELS <- c("Low", "Moderate", "High", "VeryHigh")

cat("NAS group4 distribution:\n")
nas_tab <- table(y_nas)
for (i in seq_along(nas_tab)) cat("  ", i - 1, "(", NAS_LABELS[i], "):", nas_tab[i], "\n")

# Class weights
nas_freq <- table(y_nas)
nas_present <- sort(unique(y_nas))
nas_cw <- length(y_nas) / (length(nas_freq) * nas_freq)
nas_sw <- as.numeric(nas_cw[as.character(y_nas)])

# LOCO folds: NAS datasets with >= 5 disease samples with nas_group4
nas_ds_counts <- meta_nas[, .N, by = dataset][N >= 5]
NAS_LOCO <- nas_ds_counts$dataset
cat("NAS LOCO datasets:", paste(NAS_LOCO, collapse = ", "), "\n\n")

nas_preds <- list()
nas_metrics <- list()

for (fold_ds in NAS_LOCO) {
  cat("--- NAS fold:", fold_ds, "---\n")

  test_idx  <- which(meta_nas$dataset == fold_ds)
  train_idx <- which(meta_nas$dataset != fold_ds)
  if (length(test_idx) == 0) next

  X_tr <- X_nas[train_idx, ]
  y_tr <- y_nas[train_idx]
  w_tr <- nas_sw[train_idx]
  X_te <- X_nas[test_idx, ]
  y_te <- y_nas[test_idx]

  # Ensure at least 2 training classes
  if (length(unique(y_tr)) < 2) { cat("  SKIP: <2 classes\n\n"); next }

  y_tr_fac <- factor(y_tr, levels = nas_present)

  pred_dt <- data.table(
    sample_id = meta_nas$sample_id[test_idx],
    dataset = fold_ds, y_true = y_te
  )

  # ---- Multinomial Elastic Net ----
  cat("  [1] Multinomial Elastic Net...")
  enet_fit <- tryCatch(
    cv.glmnet(X_tr, y_tr_fac, family = "multinomial", alpha = 0.5,
              weights = w_tr, nfolds = 5, type.measure = "class",
              type.multinomial = "grouped"),
    error = function(e) { cat(" ERROR:", conditionMessage(e), "\n"); NULL }
  )

  if (!is.null(enet_fit)) {
    raw_p <- predict(enet_fit, X_te, s = "lambda.min", type = "response")
    if (is.array(raw_p) && length(dim(raw_p)) == 3) raw_p <- raw_p[, , 1]
    enet_prob <- as.matrix(raw_p)
    enet_pred <- nas_present[apply(enet_prob, 1, which.max)]
    m <- compute_ordinal_metrics(y_te, enet_pred, enet_prob, nas_present)
    cat(" acc=", round(m$accuracy, 3), " QWK=", round(m$qwk, 3),
        " MAE=", round(m$mae, 2), " adjAcc=", round(m$adj_acc, 3), "\n")

    nas_metrics <- c(nas_metrics, list(data.table(
      tier = "3a_NAS", model = "multinomial_enet", fold = fold_ds,
      accuracy = m$accuracy, qwk = m$qwk, mae = m$mae,
      adj_acc = m$adj_acc, macro_auroc = m$macro_auroc, macro_f1 = m$macro_f1
    )))

    for (k in seq_along(nas_present)) {
      pred_dt[, paste0("prob_enet_", nas_present[k]) := enet_prob[, k]]
    }
    pred_dt[, pred_enet := enet_pred]
  }

  # ---- ordinalNet (cumulative logit) ----
  if (has_ordinalNet) {
    cat("  [2] ordinalNet (cumulative logit)...")
    onet_fit <- tryCatch({
      ordinalNet(X_tr, factor(y_tr, levels = nas_present, ordered = TRUE),
                 family = "cumulative", link = "logit",
                 alpha = 0.5, lambdaVals = NULL,
                 nFolds = 5, printIter = FALSE)
    }, error = function(e) { cat(" ERROR:", conditionMessage(e), "\n"); NULL })

    if (!is.null(onet_fit)) {
      onet_prob <- predict(onet_fit, X_te, type = "response")
      onet_pred <- nas_present[apply(onet_prob, 1, which.max)]
      m2 <- compute_ordinal_metrics(y_te, onet_pred, onet_prob, nas_present)
      cat(" acc=", round(m2$accuracy, 3), " QWK=", round(m2$qwk, 3), "\n")

      nas_metrics <- c(nas_metrics, list(data.table(
        tier = "3a_NAS", model = "ordinalNet", fold = fold_ds,
        accuracy = m2$accuracy, qwk = m2$qwk, mae = m2$mae,
        adj_acc = m2$adj_acc, macro_auroc = m2$macro_auroc, macro_f1 = m2$macro_f1
      )))
      pred_dt[, pred_ordinal := onet_pred]
    }
  }

  # Final prediction: use enet as default
  if (!"pred_enet" %in% names(pred_dt)) pred_dt[, pred_enet := NA_integer_]
  pred_dt[, pred_final := pred_enet]

  nas_preds <- c(nas_preds, list(pred_dt))
  cat("\n")
}

nas_results <- rbindlist(nas_preds, fill = TRUE)
fwrite(nas_results, file.path(OUTDIR, "tier3a_nas_results.csv"))
cat("NAS predictions saved:", nrow(nas_results), "samples\n\n")

# ============================================================
# TIER 3b: Fibrosis Ordinal (F0-F4)
# ============================================================
cat("=== TIER 3b: Fibrosis Ordinal Classification ===\n\n")

fib_mask <- meta_dis$fib_stage >= 0 & meta_dis$dataset %in% FIB_DATASETS
meta_fib <- meta_dis[fib_mask]
X_fib <- X_dis[fib_mask, ]
y_fib <- meta_fib$fib_stage

FIB_LABELS <- paste0("F", 0:4)
fib_present <- sort(unique(y_fib))

cat("Fibrosis stage distribution:\n")
fib_tab <- table(y_fib)
for (s in names(fib_tab)) cat("  F", s, ":", fib_tab[s], "\n")

# Class weights
fib_freq <- table(y_fib)
fib_cw <- length(y_fib) / (length(fib_freq) * fib_freq)
fib_sw <- as.numeric(fib_cw[as.character(y_fib)])

# LOCO
fib_ds_counts <- meta_fib[, .N, by = dataset][N >= 5]
FIB_LOCO <- fib_ds_counts$dataset
cat("Fibrosis LOCO datasets:", paste(FIB_LOCO, collapse = ", "), "\n\n")

fib_preds <- list()
fib_metrics <- list()

for (fold_ds in FIB_LOCO) {
  cat("--- Fibrosis fold:", fold_ds, "---\n")

  test_idx  <- which(meta_fib$dataset == fold_ds)
  train_idx <- which(meta_fib$dataset != fold_ds)
  if (length(test_idx) == 0) next

  X_tr <- X_fib[train_idx, ]
  y_tr <- y_fib[train_idx]
  w_tr <- fib_sw[train_idx]
  X_te <- X_fib[test_idx, ]
  y_te <- y_fib[test_idx]

  if (length(unique(y_tr)) < 2) { cat("  SKIP: <2 classes\n\n"); next }

  # Use only classes present in training for factor levels
  train_present <- sort(unique(y_tr))
  y_tr_fac <- factor(y_tr, levels = train_present)

  pred_dt <- data.table(
    sample_id = meta_fib$sample_id[test_idx],
    dataset = fold_ds, y_true = y_te
  )

  # ---- Multinomial Elastic Net ----
  cat("  [1] Multinomial Elastic Net...")
  enet_fit <- tryCatch(
    cv.glmnet(X_tr, y_tr_fac, family = "multinomial", alpha = 0.5,
              weights = w_tr, nfolds = 5, type.measure = "class",
              type.multinomial = "grouped"),
    error = function(e) { cat(" ERROR:", conditionMessage(e), "\n"); NULL }
  )

  if (!is.null(enet_fit)) {
    raw_p <- predict(enet_fit, X_te, s = "lambda.min", type = "response")
    if (is.array(raw_p) && length(dim(raw_p)) == 3) raw_p <- raw_p[, , 1]
    enet_prob <- as.matrix(raw_p)

    # Map probabilities back to full F0-F4 range, filling missing classes with 0
    full_prob <- matrix(0, nrow = nrow(enet_prob), ncol = 5)
    colnames(full_prob) <- as.character(0:4)
    for (k in seq_along(train_present)) {
      full_prob[, as.character(train_present[k])] <- enet_prob[, k]
    }

    enet_pred <- as.integer(colnames(full_prob)[apply(full_prob, 1, which.max)])
    m <- compute_ordinal_metrics(y_te, enet_pred, full_prob, 0:4)
    cat(" acc=", round(m$accuracy, 3), " QWK=", round(m$qwk, 3),
        " MAE=", round(m$mae, 2), " adjAcc=", round(m$adj_acc, 3), "\n")

    fib_metrics <- c(fib_metrics, list(data.table(
      tier = "3b_fibrosis", model = "multinomial_enet", fold = fold_ds,
      accuracy = m$accuracy, qwk = m$qwk, mae = m$mae,
      adj_acc = m$adj_acc, macro_auroc = m$macro_auroc, macro_f1 = m$macro_f1
    )))

    for (k in 0:4) {
      pred_dt[, paste0("prob_F", k) := full_prob[, as.character(k)]]
    }
    pred_dt[, pred_enet := enet_pred]
  }

  # ---- ordinalNet ----
  if (has_ordinalNet) {
    cat("  [2] ordinalNet...")
    onet_fit <- tryCatch({
      ordinalNet(X_tr, factor(y_tr, levels = train_present, ordered = TRUE),
                 family = "cumulative", link = "logit",
                 alpha = 0.5, nFolds = 5, printIter = FALSE)
    }, error = function(e) { cat(" ERROR:", conditionMessage(e), "\n"); NULL })

    if (!is.null(onet_fit)) {
      onet_prob <- predict(onet_fit, X_te, type = "response")
      # Map to full F0-F4 range
      full_onet <- matrix(0, nrow = nrow(onet_prob), ncol = 5)
      colnames(full_onet) <- as.character(0:4)
      for (k in seq_along(train_present)) {
        full_onet[, as.character(train_present[k])] <- onet_prob[, k]
      }

      onet_pred <- as.integer(colnames(full_onet)[apply(full_onet, 1, which.max)])
      m2 <- compute_ordinal_metrics(y_te, onet_pred, full_onet, 0:4)
      cat(" acc=", round(m2$accuracy, 3), " QWK=", round(m2$qwk, 3), "\n")

      fib_metrics <- c(fib_metrics, list(data.table(
        tier = "3b_fibrosis", model = "ordinalNet", fold = fold_ds,
        accuracy = m2$accuracy, qwk = m2$qwk, mae = m2$mae,
        adj_acc = m2$adj_acc, macro_auroc = m2$macro_auroc, macro_f1 = m2$macro_f1
      )))
      pred_dt[, pred_ordinal := onet_pred]
    }
  }

  if (!"pred_enet" %in% names(pred_dt)) pred_dt[, pred_enet := NA_integer_]
  pred_dt[, pred_final := pred_enet]

  fib_preds <- c(fib_preds, list(pred_dt))
  cat("\n")
}

fib_results <- rbindlist(fib_preds, fill = TRUE)
fwrite(fib_results, file.path(OUTDIR, "tier3b_fibrosis_results.csv"))
cat("Fibrosis predictions saved:", nrow(fib_results), "samples\n\n")

# ============================================================
# TIER 3c: Joint NAS + Fibrosis (multi-output)
# ============================================================
cat("=== TIER 3c: Joint NAS + Fibrosis ===\n\n")

# Only samples with BOTH NAS group4 >= 0 AND fib_stage >= 0
joint_mask <- meta_dis$nas_group4 >= 0 & meta_dis$fib_stage >= 0
meta_joint <- meta_dis[joint_mask]
X_joint <- X_dis[joint_mask, ]
y_nas_joint <- meta_joint$nas_group4
y_fib_joint <- meta_joint$fib_stage

cat("Joint-annotated samples:", nrow(meta_joint), "\n")
cat("NAS distribution:\n"); print(table(y_nas_joint))
cat("Fibrosis distribution:\n"); print(table(y_fib_joint))

# LOCO across datasets with joint annotation
joint_ds_counts <- meta_joint[, .N, by = dataset][N >= 5]
JOINT_LOCO <- joint_ds_counts$dataset
cat("Joint LOCO datasets:", paste(JOINT_LOCO, collapse = ", "), "\n\n")

joint_preds <- list()
joint_metrics <- list()

for (fold_ds in JOINT_LOCO) {
  cat("--- Joint fold:", fold_ds, "---\n")

  test_idx  <- which(meta_joint$dataset == fold_ds)
  train_idx <- which(meta_joint$dataset != fold_ds)
  if (length(test_idx) == 0) next

  X_tr <- X_joint[train_idx, ]
  X_te <- X_joint[test_idx, ]

  y_nas_tr <- y_nas_joint[train_idx]
  y_nas_te <- y_nas_joint[test_idx]
  y_fib_tr <- y_fib_joint[train_idx]
  y_fib_te <- y_fib_joint[test_idx]

  if (length(unique(y_nas_tr)) < 2 || length(unique(y_fib_tr)) < 2) {
    cat("  SKIP: insufficient class variety\n\n")
    next
  }

  pred_dt <- data.table(
    sample_id = meta_joint$sample_id[test_idx],
    dataset = fold_ds,
    y_true_nas = y_nas_te,
    y_true_fib = y_fib_te
  )

  # ---- XGBoost for NAS (via caret) ----
  cat("  [NAS] XGBoost...")
  nas_train_present <- sort(unique(y_nas_tr))
  y_nas_fac <- factor(paste0("N", y_nas_tr), levels = paste0("N", nas_train_present))
  nas_cw_j <- length(y_nas_tr) / (length(nas_train_present) * table(y_nas_fac))
  nas_sw_j <- as.numeric(nas_cw_j[as.character(y_nas_fac)])

  xgb_nas <- tryCatch({
    ctrl <- trainControl(method = "cv", number = 3, classProbs = TRUE,
                         verboseIter = FALSE, allowParallel = FALSE)
    grid <- expand.grid(
      nrounds = 150, max_depth = 4, eta = 0.1, gamma = 0,
      colsample_bytree = 0.8, min_child_weight = 1, subsample = 0.8
    )
    train(x = as.data.frame(X_tr), y = y_nas_fac,
          method = "xgbTree", trControl = ctrl, tuneGrid = grid,
          weights = nas_sw_j, verbosity = 0)
  }, error = function(e) { cat(" ERROR:", conditionMessage(e), "\n"); NULL })

  if (!is.null(xgb_nas)) {
    nas_prob <- predict(xgb_nas, as.data.frame(X_te), type = "prob")
    # Map back to integer predictions
    nas_pred <- nas_train_present[apply(as.matrix(nas_prob), 1, which.max)]
    m_nas <- compute_ordinal_metrics(y_nas_te, nas_pred, class_labels = nas_train_present)
    cat(" acc=", round(m_nas$accuracy, 3), " QWK=", round(m_nas$qwk, 3))
    pred_dt[, pred_nas_joint := nas_pred]
  } else {
    pred_dt[, pred_nas_joint := NA_integer_]
    m_nas <- NULL
  }

  # ---- XGBoost for Fibrosis (via caret) ----
  cat("  [FIB] XGBoost...")
  fib_train_present <- sort(unique(y_fib_tr))
  y_fib_fac <- factor(paste0("F", y_fib_tr), levels = paste0("F", fib_train_present))
  fib_cw_j <- length(y_fib_tr) / (length(fib_train_present) * table(y_fib_fac))
  fib_sw_j <- as.numeric(fib_cw_j[as.character(y_fib_fac)])

  xgb_fib <- tryCatch({
    ctrl <- trainControl(method = "cv", number = 3, classProbs = TRUE,
                         verboseIter = FALSE, allowParallel = FALSE)
    grid <- expand.grid(
      nrounds = 150, max_depth = 4, eta = 0.1, gamma = 0,
      colsample_bytree = 0.8, min_child_weight = 1, subsample = 0.8
    )
    train(x = as.data.frame(X_tr), y = y_fib_fac,
          method = "xgbTree", trControl = ctrl, tuneGrid = grid,
          weights = fib_sw_j, verbosity = 0)
  }, error = function(e) { cat(" ERROR:", conditionMessage(e), "\n"); NULL })

  if (!is.null(xgb_fib)) {
    fib_prob <- predict(xgb_fib, as.data.frame(X_te), type = "prob")
    fib_pred <- fib_train_present[apply(as.matrix(fib_prob), 1, which.max)]
    m_fib <- compute_ordinal_metrics(y_fib_te, fib_pred, class_labels = fib_train_present)
    cat(" acc=", round(m_fib$accuracy, 3), " QWK=", round(m_fib$qwk, 3), "\n")
    pred_dt[, pred_fib_joint := fib_pred]
  } else {
    pred_dt[, pred_fib_joint := NA_integer_]
    m_fib <- NULL
  }

  # Record metrics
  if (!is.null(m_nas)) {
    joint_metrics <- c(joint_metrics, list(data.table(
      tier = "3c_joint_NAS", model = "xgboost_joint", fold = fold_ds,
      accuracy = m_nas$accuracy, qwk = m_nas$qwk, mae = m_nas$mae,
      adj_acc = m_nas$adj_acc, macro_auroc = m_nas$macro_auroc, macro_f1 = m_nas$macro_f1
    )))
  }
  if (!is.null(m_fib)) {
    joint_metrics <- c(joint_metrics, list(data.table(
      tier = "3c_joint_FIB", model = "xgboost_joint", fold = fold_ds,
      accuracy = m_fib$accuracy, qwk = m_fib$qwk, mae = m_fib$mae,
      adj_acc = m_fib$adj_acc, macro_auroc = m_fib$macro_auroc, macro_f1 = m_fib$macro_f1
    )))
  }

  joint_preds <- c(joint_preds, list(pred_dt))
  cat("\n")
}

joint_results <- rbindlist(joint_preds, fill = TRUE)
fwrite(joint_results, file.path(OUTDIR, "tier3c_joint_results.csv"))
cat("Joint predictions saved:", nrow(joint_results), "samples\n")

# --- Compare independent vs joint predictions ---
# Load Tier 3a/3b independent predictions to compare
if (nrow(joint_results) > 0 && nrow(nas_results) > 0 && nrow(fib_results) > 0) {
  cat("\n--- Concordance: independent vs joint ---\n")

  merged <- merge(
    joint_results[, .(sample_id, pred_nas_joint, pred_fib_joint, y_true_nas, y_true_fib)],
    nas_results[, .(sample_id, pred_nas_indep = pred_final)],
    by = "sample_id", all.x = TRUE
  )
  merged <- merge(
    merged,
    fib_results[, .(sample_id, pred_fib_indep = pred_final)],
    by = "sample_id", all.x = TRUE
  )

  # Concordance rate
  nas_concordance <- mean(merged$pred_nas_joint == merged$pred_nas_indep, na.rm = TRUE)
  fib_concordance <- mean(merged$pred_fib_joint == merged$pred_fib_indep, na.rm = TRUE)
  cat("  NAS concordance (joint vs independent):", round(nas_concordance, 3), "\n")
  cat("  Fibrosis concordance (joint vs independent):", round(fib_concordance, 3), "\n")

  joint_results[, nas_concordance := nas_concordance]
  joint_results[, fib_concordance := fib_concordance]
  fwrite(joint_results, file.path(OUTDIR, "tier3c_joint_results.csv"))
}

# ============================================================
# Aggregate all Tier 3 metrics
# ============================================================
cat("\n=== Aggregating Tier 3 Metrics ===\n\n")

all_t3_metrics <- rbindlist(c(nas_metrics, fib_metrics, joint_metrics), fill = TRUE)
fwrite(all_t3_metrics, file.path(OUTDIR, "tier3_model_summary.csv"))

# Summary table
t3_summary <- all_t3_metrics[, .(
  mean_accuracy = mean(accuracy, na.rm = TRUE),
  mean_qwk      = mean(qwk, na.rm = TRUE),
  mean_mae      = mean(mae, na.rm = TRUE),
  mean_adj_acc  = mean(adj_acc, na.rm = TRUE),
  mean_macro_f1 = mean(macro_f1, na.rm = TRUE),
  n_folds = .N
), by = .(tier, model)]

cat("Tier 3 Summary:\n")
print(t3_summary)

# ============================================================
# Final Summary
# ============================================================
cat("\n=== FINAL SUMMARY ===\n")
cat("Tier 3a (NAS) samples:", nrow(nas_results), "\n")
cat("Tier 3b (Fibrosis) samples:", nrow(fib_results), "\n")
cat("Tier 3c (Joint) samples:", nrow(joint_results), "\n")

if (nrow(all_t3_metrics) > 0) {
  best_nas <- all_t3_metrics[tier == "3a_NAS"][which.max(mean(qwk, na.rm = TRUE))]
  best_fib <- all_t3_metrics[tier == "3b_fibrosis"][which.max(mean(qwk, na.rm = TRUE))]
  if (nrow(best_nas) > 0)
    cat("Best NAS model:", best_nas$model[1], "mean QWK=", round(mean(all_t3_metrics[tier == "3a_NAS" & model == best_nas$model[1], qwk], na.rm = TRUE), 3), "\n")
  if (nrow(best_fib) > 0)
    cat("Best Fibrosis model:", best_fib$model[1], "mean QWK=", round(mean(all_t3_metrics[tier == "3b_fibrosis" & model == best_fib$model[1], qwk], na.rm = TRUE), 3), "\n")
}

cat("\n=== 67_tier3_staging.R completed:", as.character(Sys.time()), "===\n")
