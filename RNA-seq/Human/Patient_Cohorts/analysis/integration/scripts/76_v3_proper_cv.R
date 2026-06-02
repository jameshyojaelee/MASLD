#!/usr/bin/env Rscript
# 76_v3_proper_cv.R
# V3 Staging Classifier with PROPER within-fold feature selection
#
# Fixes the V2 information leakage: feature pre-filter (top 3K genes by dream
# t-stats) was computed on ALL samples including test folds. V3 performs
# feature selection INSIDE each LOCO fold using training-fold variance only.
#
# Models: Ordinal elastic net for fibrosis (F0-F4), binary elastic net (F>=3),
#         binary elastic net (NAS>=5), ordinal elastic net for NAS (4 groups)
#
# Validation: Leave-one-cohort-out (LOCO) with 6 fibrosis folds, 5 NAS folds
#
# Outputs:
#   - v3_proper_cv_results.csv: Per-sample predictions across all folds
#   - v3_feature_stability_per_fold.csv: Which genes selected per fold + Jaccard
#   - v3_vs_v2_comparison.csv: Head-to-head AUROC comparison V2 vs V3
#   - v3_model_summary.csv: Aggregated metrics per target
#
# Usage: Rscript 76_v3_proper_cv.R
# SLURM: bigmem, 16 CPUs, 200GB RAM, 48h

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(glmnet)
  library(pROC)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR <- file.path(INT, "results")
OUTDIR <- file.path(RDIR, "staging_classifier")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

cat("=== 76: V3 Proper Cross-Validation (No Feature Leakage) ===\n")
cat("Started:", as.character(Sys.time()), "\n\n")

K_FEATURES <- 3000  # Number of features to select per fold

# --- Load data ---
cat("Loading merged DGE (all 1,444 samples)...\n")
dge <- readRDS(file.path(RDIR, "integration/merged_dge.rds"))
dge <- calcNormFactors(dge, method = "TMM")
logcpm <- cpm(dge, log = TRUE, prior.count = 1)
cat("  Expression:", nrow(logcpm), "genes x", ncol(logcpm), "samples\n")

cat("Loading modeling metadata...\n")
meta <- fread(file.path(OUTDIR, "modeling_metadata.csv"))
cat("  Metadata:", nrow(meta), "samples\n")

# Align samples
common <- intersect(colnames(logcpm), meta$sample_id)
meta <- meta[match(common, sample_id)]
logcpm <- logcpm[, common]
cat("  Aligned:", ncol(logcpm), "samples\n\n")

# ============================================================
# Helper: Within-fold rank transform + feature selection
# ============================================================
select_features_and_transform <- function(logcpm_train, logcpm_test, k = K_FEATURES) {
  # 1. Compute per-gene variance on TRAINING data only (no label info → no leakage)
  gene_var <- apply(logcpm_train, 1, var)

  # 2. Select top K by variance
  top_k_genes <- names(sort(gene_var, decreasing = TRUE))[1:min(k, length(gene_var))]

  # 3. Rank transform both train and test (using only selected genes)
  rank_train <- apply(logcpm_train[top_k_genes, , drop = FALSE], 2, function(x) {
    r <- rank(x, ties.method = "average")
    r / length(r)
  })
  rank_test <- apply(logcpm_test[top_k_genes, , drop = FALSE], 2, function(x) {
    r <- rank(x, ties.method = "average")
    r / length(r)
  })

  list(
    X_train = t(rank_train),  # samples x genes
    X_test  = t(rank_test),
    selected_genes = top_k_genes
  )
}

# ============================================================
# Helper: Compute QWK
# ============================================================
compute_qwk <- function(pred, true, n_classes) {
  O <- table(factor(pred, levels = 0:(n_classes-1)),
             factor(true, levels = 0:(n_classes-1)))
  n <- sum(O)
  if (n == 0) return(NA_real_)
  W <- outer(0:(n_classes-1), 0:(n_classes-1), function(i,j) (i-j)^2 / (n_classes-1)^2)
  E <- outer(rowSums(O), colSums(O)) / n
  num <- sum(W * O); den <- sum(W * E)
  if (den == 0) return(1.0)
  1 - num / den
}

# ============================================================
# ANALYSIS: Proper LOCO-CV for each target
# ============================================================
targets <- list(
  fib_ge3 = list(
    col = "fib_ge3", fold_col = "loco_fold_fibrosis",
    type = "binary", desc = "F>=3 binary"
  ),
  nas_ge5 = list(
    col = "nas_ge5", fold_col = "loco_fold_nas",
    type = "binary", desc = "NAS>=5 binary"
  ),
  fib_5class = list(
    col = "fib_stage", fold_col = "loco_fold_fibrosis",
    type = "ordinal", n_classes = 5, desc = "Fibrosis F0-F4"
  ),
  nas_4group = list(
    col = "nas_group4", fold_col = "loco_fold_nas",
    type = "ordinal", n_classes = 4, desc = "NAS 4-group"
  )
)

all_results <- list()
all_features <- list()
all_summaries <- list()

for (tgt_name in names(targets)) {
  tgt <- targets[[tgt_name]]
  cat("\n=== Target:", tgt$desc, "===\n")

  # Get valid samples (non-excluded, non-NA labels)
  valid <- meta[[tgt$fold_col]] != "excluded" &
           !is.na(meta[[tgt$col]]) &
           meta[[tgt$col]] >= 0
  meta_valid <- meta[valid]
  logcpm_valid <- logcpm[, meta_valid$sample_id]

  folds <- unique(meta_valid[[tgt$fold_col]])
  folds <- folds[folds != "excluded" & folds != "NA"]
  cat("  Samples:", nrow(meta_valid), ", Folds:", length(folds), "\n")

  fold_results <- list()
  fold_features <- list()
  fold_aurocs <- numeric()
  fold_qwks <- numeric()

  for (fold in folds) {
    cat("  Fold:", fold, "... ")

    test_idx  <- which(meta_valid[[tgt$fold_col]] == fold)
    train_idx <- which(meta_valid[[tgt$fold_col]] != fold)

    if (length(test_idx) < 5 || length(train_idx) < 20) {
      cat("SKIP (too few samples)\n")
      next
    }

    y_train <- meta_valid[[tgt$col]][train_idx]
    y_test  <- meta_valid[[tgt$col]][test_idx]

    # Check if test fold has >=2 classes
    if (length(unique(y_test)) < 2 && tgt$type == "binary") {
      cat("SKIP (single class in test)\n")
      next
    }

    # WITHIN-FOLD feature selection + rank transform
    feat <- select_features_and_transform(
      logcpm_valid[, meta_valid$sample_id[train_idx]],
      logcpm_valid[, meta_valid$sample_id[test_idx]],
      k = K_FEATURES
    )

    # Compute class weights
    class_counts <- table(y_train)
    class_wts <- 1 / class_counts
    class_wts <- class_wts / sum(class_wts) * length(class_wts)
    sample_wts <- class_wts[as.character(y_train)]

    # Train model
    if (tgt$type == "binary") {
      fit <- tryCatch({
        cv.glmnet(
          x = feat$X_train, y = factor(y_train),
          family = "binomial", alpha = 0.5, nfolds = 5,
          type.measure = "auc", weights = as.numeric(sample_wts)
        )
      }, error = function(e) { cat("ERROR:", e$message, "\n"); NULL })

      if (is.null(fit)) next

      prob <- as.numeric(predict(fit, newx = feat$X_test, s = "lambda.min", type = "response"))
      pred <- ifelse(prob > 0.5, 1, 0)

      auroc <- tryCatch(as.numeric(auc(roc(y_test, prob, quiet = TRUE))), error = function(e) NA)
      qwk <- compute_qwk(pred, y_test, 2)
      acc <- mean(pred == y_test)

      fold_aurocs <- c(fold_aurocs, auroc)
      fold_qwks <- c(fold_qwks, qwk)
      cat("AUROC=", round(auroc, 3), " Acc=", round(acc, 3), "\n")

      # Save per-sample predictions
      fold_results[[fold]] <- data.table(
        sample_id = meta_valid$sample_id[test_idx],
        dataset = meta_valid$dataset[test_idx],
        true_label = y_test,
        predicted = pred,
        probability = prob,
        fold = fold,
        target = tgt_name
      )

    } else {
      # Ordinal: multinomial elastic net
      fit <- tryCatch({
        cv.glmnet(
          x = feat$X_train, y = factor(y_train),
          family = "multinomial", alpha = 0.5, nfolds = 5,
          type.measure = "class", weights = as.numeric(sample_wts),
          type.multinomial = "grouped"
        )
      }, error = function(e) { cat("ERROR:", e$message, "\n"); NULL })

      if (is.null(fit)) next

      prob_mat <- predict(fit, newx = feat$X_test, s = "lambda.min", type = "response")[,,1]
      pred <- as.integer(colnames(prob_mat)[apply(prob_mat, 1, which.max)])

      qwk <- compute_qwk(pred, y_test, tgt$n_classes)
      acc <- mean(pred == y_test)
      mae <- mean(abs(pred - y_test))
      adj_acc <- mean(abs(pred - y_test) <= 1)

      fold_qwks <- c(fold_qwks, qwk)
      cat("QWK=", round(qwk, 3), " Acc=", round(acc, 3), " MAE=", round(mae, 2), "\n")

      fold_results[[fold]] <- data.table(
        sample_id = meta_valid$sample_id[test_idx],
        dataset = meta_valid$dataset[test_idx],
        true_label = y_test,
        predicted = pred,
        fold = fold,
        target = tgt_name
      )
    }

    # Track selected features
    fold_features[[fold]] <- data.table(
      gene = feat$selected_genes,
      fold = fold,
      target = tgt_name
    )

    # Track non-zero coefficients
    coefs <- coef(fit, s = "lambda.min")
    if (is.list(coefs)) {
      # Multinomial: take union of non-zero coefs across classes
      nz_genes <- unique(unlist(lapply(coefs, function(c) {
        rownames(c)[which(c[,1] != 0)]
      })))
    } else {
      nz_genes <- rownames(coefs)[which(coefs[,1] != 0)]
    }
    nz_genes <- setdiff(nz_genes, "(Intercept)")

    fold_features[[fold]][, selected_by_model := gene %in% nz_genes]
  }

  # Aggregate
  if (length(fold_results) > 0) {
    tgt_results <- rbindlist(fold_results, fill = TRUE)
    all_results[[tgt_name]] <- tgt_results

    tgt_features <- rbindlist(fold_features, fill = TRUE)
    all_features[[tgt_name]] <- tgt_features

    # Summary
    all_summaries[[tgt_name]] <- data.table(
      target = tgt$desc,
      target_key = tgt_name,
      type = tgt$type,
      mean_auroc = if (length(fold_aurocs) > 0) mean(fold_aurocs, na.rm = TRUE) else NA,
      sd_auroc = if (length(fold_aurocs) > 1) sd(fold_aurocs, na.rm = TRUE) else NA,
      mean_qwk = if (length(fold_qwks) > 0) mean(fold_qwks, na.rm = TRUE) else NA,
      sd_qwk = if (length(fold_qwks) > 1) sd(fold_qwks, na.rm = TRUE) else NA,
      n_folds = length(fold_results),
      n_samples = nrow(tgt_results)
    )

    cat("  Summary: ", tgt$desc, " — ",
        if (!is.na(all_summaries[[tgt_name]]$mean_auroc))
          paste0("AUROC=", round(all_summaries[[tgt_name]]$mean_auroc, 3))
        else
          paste0("QWK=", round(all_summaries[[tgt_name]]$mean_qwk, 3)),
        " (", length(fold_results), " folds)\n")
  }
}

# ============================================================
# Save outputs
# ============================================================
cat("\n=== Saving V3 Results ===\n")

# Per-sample predictions
v3_results <- rbindlist(all_results, fill = TRUE)
fwrite(v3_results, file.path(OUTDIR, "v3_proper_cv_results.csv"))
cat("  v3_proper_cv_results.csv:", nrow(v3_results), "rows\n")

# Feature stability
v3_features <- rbindlist(all_features, fill = TRUE)
fwrite(v3_features, file.path(OUTDIR, "v3_feature_stability_per_fold.csv"))
cat("  v3_feature_stability_per_fold.csv:", nrow(v3_features), "rows\n")

# Compute Jaccard similarity between folds
cat("\n  Feature stability (Jaccard between folds):\n")
for (tgt_name in names(targets)) {
  feat_tgt <- v3_features[target == tgt_name]
  folds_present <- unique(feat_tgt$fold)
  if (length(folds_present) >= 2) {
    jaccards <- numeric()
    for (i in 1:(length(folds_present)-1)) {
      for (j in (i+1):length(folds_present)) {
        g1 <- feat_tgt[fold == folds_present[i], gene]
        g2 <- feat_tgt[fold == folds_present[j], gene]
        jac <- length(intersect(g1, g2)) / length(union(g1, g2))
        jaccards <- c(jaccards, jac)
      }
    }
    cat("    ", tgt_name, ": mean Jaccard =", round(mean(jaccards), 3),
        "(range", round(min(jaccards), 3), "-", round(max(jaccards), 3), ")\n")
  }
}

# Model summary
v3_summary <- rbindlist(all_summaries, fill = TRUE)
fwrite(v3_summary, file.path(OUTDIR, "v3_model_summary.csv"))
cat("\n  v3_model_summary.csv:\n")
print(v3_summary[, .(target, mean_auroc, mean_qwk, n_folds)])

# V2 vs V3 comparison
cat("\n=== V2 vs V3 Comparison ===\n")
v2_file <- file.path(OUTDIR, "ordinal_model_summary.csv")
if (file.exists(v2_file)) {
  v2_summary <- fread(v2_file)
  comparison <- data.table(
    target = v3_summary$target,
    v3_auroc = v3_summary$mean_auroc,
    v3_qwk = v3_summary$mean_qwk,
    v3_folds = v3_summary$n_folds
  )
  # Match V2 targets
  v2_map <- list(
    "F>=3 binary" = "fib_ge3_binary",
    "NAS>=5 binary" = "nas_ge5_binary",
    "Fibrosis F0-F4" = "fibrosis_ordinal",
    "NAS 4-group" = "nas_ordinal"
  )
  for (i in seq_len(nrow(comparison))) {
    v2_key <- v2_map[[comparison$target[i]]]
    if (!is.null(v2_key)) {
      v2_row <- v2_summary[target == v2_key]
      if (nrow(v2_row) > 0) {
        comparison[i, v2_auroc := v2_row$mean_primary[1]]
        comparison[i, v2_folds := v2_row$n_folds[1]]
      }
    }
  }
  comparison[, delta := fifelse(!is.na(v3_auroc) & !is.na(v2_auroc), v3_auroc - v2_auroc,
                         fifelse(!is.na(v3_qwk) & !is.na(v2_auroc), v3_qwk - v2_auroc, NA_real_))]

  fwrite(comparison, file.path(OUTDIR, "v3_vs_v2_comparison.csv"))
  cat("V2 vs V3 comparison saved\n")
  print(comparison)
}

cat("\n=== 76_v3_proper_cv.R completed:", as.character(Sys.time()), "===\n")
