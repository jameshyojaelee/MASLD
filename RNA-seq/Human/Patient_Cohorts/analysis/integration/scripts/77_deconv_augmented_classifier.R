#!/usr/bin/env Rscript
# =============================================================================
# 77_deconv_augmented_classifier.R
#
# Test whether adding cell-type deconvolution proportions improves fibrosis
# staging (F>=3 binary classification). Validates Kim 2025, Liver International
# at 10-cohort scale.
#
# Three models per LOCO fold:
#   A. Expression-only (top 3000 genes by training-fold variance, rank-transformed)
#   B. Expression + deconvolution proportions
#   C. Deconvolution-only
#
# Proper within-fold feature selection (no leakage).
#
# SLURM:
#   #SBATCH --partition=cpu
#   #SBATCH --cpus-per-task=8
#   #SBATCH --mem=64G
#   #SBATCH --time=48:00:00
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(glmnet)
  library(pROC)
})

set.seed(42)

# --- Paths -------------------------------------------------------------------
project_root <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

dge_path   <- file.path(project_root,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/merged_dge.rds")
meta_path  <- file.path(project_root,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/staging_classifier/modeling_metadata.csv")
deconv_dir <- file.path(project_root,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/deconvolution/bayesprism")
out_dir    <- file.path(project_root,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/staging_classifier")

# --- Locate deconvolution file -----------------------------------------------
# Primary: known path
deconv_file <- file.path(deconv_dir, "unified_bayesprism_proportions.csv")
if (!file.exists(deconv_file)) {
  # Fallback 1: recursive glob under results/deconvolution/
  deconv_base <- file.path(project_root,
    "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/deconvolution")
  hits <- list.files(deconv_base, pattern = "bayesprism.*proportions",
                     recursive = TRUE, full.names = TRUE)
  if (length(hits) == 0) {
    # Fallback 2: results/integration/deconvolution/
    deconv_base2 <- file.path(project_root,
      "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/deconvolution")
    hits <- list.files(deconv_base2, pattern = "bayesprism.*proportions",
                       recursive = TRUE, full.names = TRUE)
  }
  if (length(hits) == 0) {
    stop("ERROR: No file matching *bayesprism*proportions* found under ",
         deconv_base, " or ", deconv_base2,
         "\nCannot proceed without cell-type proportions.")
  }
  deconv_file <- hits[1]
  message("Using fallback deconvolution file: ", deconv_file)
}

# --- 1. Load data ------------------------------------------------------------
cat("=== Loading data ===\n")

# DGE -> TMM logCPM
cat("  Loading DGE ...\n")
dge <- readRDS(dge_path)
dge <- calcNormFactors(dge, method = "TMM")
logcpm <- cpm(dge, log = TRUE, prior.count = 1)  # 34453 x 1444
cat("  logCPM:", nrow(logcpm), "genes x", ncol(logcpm), "samples\n")

# Metadata
cat("  Loading metadata ...\n")
meta <- fread(meta_path)
cat("  Metadata:", nrow(meta), "samples\n")

# Deconvolution proportions
cat("  Loading deconvolution: ", basename(deconv_file), "\n")
deconv <- fread(deconv_file)
cat("  Deconvolution:", nrow(deconv), "samples x", ncol(deconv), "columns\n")

# Identify proportion columns (exclude sample_id, dataset)
prop_cols <- setdiff(names(deconv), c("sample_id", "dataset"))
cat("  Proportion columns:", paste(prop_cols, collapse = ", "), "\n")
n_deconv_features <- length(prop_cols)

# --- 2. Filter & align -------------------------------------------------------
cat("\n=== Filtering & aligning ===\n")

# Keep only samples with valid fibrosis labels and valid LOCO folds
meta_valid <- meta[fib_ge3 != -1 & loco_fold_fibrosis != "excluded"]
cat("  Samples with fibrosis labels:", nrow(meta_valid), "\n")

# Intersect all three data sources
common_ids <- Reduce(intersect, list(
  meta_valid$sample_id,
  colnames(logcpm),
  deconv$sample_id
))
cat("  Common samples across all sources:", length(common_ids), "\n")

# Align
meta_valid <- meta_valid[match(common_ids, sample_id)]
logcpm     <- logcpm[, common_ids, drop = FALSE]
deconv_dt  <- deconv[match(common_ids, sample_id)]

# Extract deconv matrix (samples x cell types)
deconv_mat <- as.matrix(deconv_dt[, ..prop_cols])
rownames(deconv_mat) <- common_ids

# Response
y <- meta_valid$fib_ge3
folds <- meta_valid$loco_fold_fibrosis
fold_names <- sort(unique(folds))
n_folds <- length(fold_names)

cat("  Outcome (fib_ge3): ", table(y)["0"], " negative, ", table(y)["1"], " positive\n")
cat("  LOCO folds (", n_folds, "): ", paste(fold_names, collapse = ", "), "\n")

# --- 3. LOCO-CV with proper within-fold feature selection --------------------
cat("\n=== Running LOCO-CV ===\n")

n_top_genes <- 3000

results_list <- list()
importance_list <- list()

for (fi in seq_along(fold_names)) {
  fold <- fold_names[fi]
  cat(sprintf("\n--- Fold %d/%d: %s (held out) ---\n", fi, n_folds, fold))

  test_idx  <- which(folds == fold)
  train_idx <- which(folds != fold)
  cat(sprintf("  Train: %d  Test: %d\n", length(train_idx), length(test_idx)))

  y_train <- y[train_idx]
  y_test  <- y[test_idx]

  # Check test fold has both classes
  if (length(unique(y_test)) < 2) {
    cat("  WARNING: Test fold has only one class, skipping AUROC.\n")
    next
  }

  # --- Inverse-frequency class weights ---
  tab <- table(y_train)
  w_train <- ifelse(y_train == 1,
                    length(y_train) / (2 * tab["1"]),
                    length(y_train) / (2 * tab["0"]))

  # --- Feature selection on TRAINING data only ---
  train_logcpm <- logcpm[, train_idx, drop = FALSE]
  gene_vars <- apply(train_logcpm, 1, var)
  top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:n_top_genes]

  # --- Rank-transform expression (within each set separately) ---
  # Training: rank across training samples
  expr_train_raw <- t(logcpm[top_genes, train_idx, drop = FALSE])  # n_train x 3000
  expr_train <- apply(expr_train_raw, 2, rank) / nrow(expr_train_raw)

  # Test: rank across test samples
  expr_test_raw <- t(logcpm[top_genes, test_idx, drop = FALSE])    # n_test x 3000
  expr_test <- apply(expr_test_raw, 2, rank) / nrow(expr_test_raw)

  # --- Deconv features ---
  deconv_train <- deconv_mat[train_idx, , drop = FALSE]
  deconv_test  <- deconv_mat[test_idx, , drop = FALSE]

  # --- Build model matrices ---
  # Model 1: Expression-only
  x_train_expr <- expr_train
  x_test_expr  <- expr_test

  # Model 2: Expression + deconv
  x_train_combo <- cbind(expr_train, deconv_train)
  x_test_combo  <- cbind(expr_test, deconv_test)

  # Model 3: Deconv-only
  x_train_deconv <- deconv_train
  x_test_deconv  <- deconv_test

  models <- list(
    list(name = "expression_only",   x_train = x_train_expr,   x_test = x_test_expr),
    list(name = "expression_deconv", x_train = x_train_combo,  x_test = x_test_combo),
    list(name = "deconv_only",       x_train = x_train_deconv, x_test = x_test_deconv)
  )

  for (m in models) {
    cat(sprintf("  Training %s (p=%d) ... ", m$name, ncol(m$x_train)))

    # cv.glmnet for lambda selection
    cvfit <- tryCatch(
      cv.glmnet(
        x = m$x_train,
        y = y_train,
        weights = w_train,
        family = "binomial",
        alpha = 0.5,
        nfolds = 5,
        type.measure = "auc"
      ),
      error = function(e) {
        cat("ERROR:", conditionMessage(e), "\n")
        return(NULL)
      }
    )

    if (is.null(cvfit)) {
      results_list[[length(results_list) + 1]] <- data.table(
        fold = fold, model_type = m$name, auroc = NA_real_, accuracy = NA_real_
      )
      next
    }

    # Predict on test
    pred_prob <- as.numeric(predict(cvfit, newx = m$x_test, s = "lambda.1se",
                                    type = "response"))
    pred_class <- ifelse(pred_prob >= 0.5, 1, 0)

    # AUROC
    roc_obj <- roc(y_test, pred_prob, quiet = TRUE)
    auc_val <- as.numeric(auc(roc_obj))
    acc_val <- mean(pred_class == y_test)

    cat(sprintf("AUROC=%.3f  Acc=%.3f\n", auc_val, acc_val))

    results_list[[length(results_list) + 1]] <- data.table(
      fold = fold, model_type = m$name, auroc = auc_val, accuracy = acc_val
    )

    # --- Extract deconv feature importance (non-zero coefficients) ---
    if (m$name %in% c("expression_deconv", "deconv_only")) {
      coefs <- as.matrix(coef(cvfit, s = "lambda.1se"))
      # Keep only deconv feature rows (exclude intercept)
      deconv_coef_names <- intersect(rownames(coefs), prop_cols)
      if (length(deconv_coef_names) > 0) {
        deconv_coefs <- coefs[deconv_coef_names, , drop = FALSE]
        nonzero <- deconv_coefs[deconv_coefs[, 1] != 0, , drop = FALSE]
        if (nrow(nonzero) > 0) {
          importance_list[[length(importance_list) + 1]] <- data.table(
            fold = fold,
            model_type = m$name,
            feature = rownames(nonzero),
            coefficient = nonzero[, 1]
          )
        }
      }
    }
  }
}

# --- 4. Assemble results ----------------------------------------------------
cat("\n=== Assembling results ===\n")

results <- rbindlist(results_list)
cat("Per-fold results:\n")
print(results)

# Summary: mean AUROC per model type
summary_dt <- results[, .(
  mean_auroc = mean(auroc, na.rm = TRUE),
  sd_auroc   = sd(auroc, na.rm = TRUE),
  mean_acc   = mean(accuracy, na.rm = TRUE),
  n_folds    = sum(!is.na(auroc))
), by = model_type]

# Compute delta from adding deconv
expr_auroc <- summary_dt[model_type == "expression_only", mean_auroc]
summary_dt[, delta_vs_expr := mean_auroc - expr_auroc]

cat("\nAblation summary:\n")
print(summary_dt)

# Feature importance
if (length(importance_list) > 0) {
  importance <- rbindlist(importance_list)
  cat("\nDeconv feature importance (non-zero coefficients):\n")
  print(importance)
} else {
  importance <- data.table(fold = character(), model_type = character(),
                           feature = character(), coefficient = numeric())
  cat("\nNo deconv features had non-zero coefficients.\n")
}

# --- 5. Write outputs --------------------------------------------------------
cat("\n=== Writing outputs ===\n")

fwrite(results,
       file.path(out_dir, "deconv_augmented_results.csv"))
cat("  Wrote deconv_augmented_results.csv\n")

fwrite(summary_dt,
       file.path(out_dir, "deconv_ablation.csv"))
cat("  Wrote deconv_ablation.csv\n")

fwrite(importance,
       file.path(out_dir, "deconv_feature_importance.csv"))
cat("  Wrote deconv_feature_importance.csv\n")

cat("\n=== Done ===\n")
cat(sprintf("Expression-only mean AUROC: %.3f\n", expr_auroc))
cat(sprintf("Expression+deconv mean AUROC: %.3f (delta: %+.3f)\n",
            summary_dt[model_type == "expression_deconv", mean_auroc],
            summary_dt[model_type == "expression_deconv", delta_vs_expr]))
cat(sprintf("Deconv-only mean AUROC: %.3f\n",
            summary_dt[model_type == "deconv_only", mean_auroc]))
