#!/usr/bin/env Rscript
# =============================================================================
# 78_sex_stratified_classifier.R
#
# Sex-stratified MASLD fibrosis classifiers (F>=3 binary).
# No published MASLD classifier is sex-stratified — fills urgent clinical gap.
#
# Four strategies, all with proper within-fold LOCO-CV:
#   A: Universal — all samples, top 3K genes by training-fold variance
#   B: Sex-specific — separate female-only and male-only LOCO classifiers
#   C: Sex-as-feature — all samples, top 3K genes + sex as binary covariate
#   D: Sex-specific gene sets — female classifier uses female-specific DEGs,
#      male classifier uses male-specific DEGs
#
# Elastic net (cv.glmnet, binomial, alpha=0.5), inverse-frequency class
# weights, AUROC via pROC.
#
# Inputs:
#   - merged_dge.rds: 34,453 genes x 1,444 samples (edgeR DGEList)
#   - modeling_metadata.csv: fib_ge3, loco_fold_fibrosis, sex columns
#   - unified_metadata.csv: fallback sex source
#   - sex_deg_classification.csv: sex_class column (Female_specific, etc.)
#
# Outputs (to results/staging_classifier/):
#   - sex_stratified_results.csv  (per-fold AUROC for each strategy x sex)
#   - sex_auroc_comparison.csv    (summary table)
#   - sex_specific_panels.csv     (top genes per sex)
#
# Usage: Rscript 78_sex_stratified_classifier.R
# SLURM: cpu, 8 CPUs, 64GB RAM, 48h
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

int_dir    <- file.path(project_root,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration")
res_dir    <- file.path(int_dir, "results/integration")
out_dir    <- file.path(int_dir, "results/staging_classifier")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

dge_path      <- file.path(res_dir, "merged_dge.rds")
meta_path     <- file.path(out_dir, "modeling_metadata.csv")
sex_deg_path  <- file.path(res_dir, "sex_deg_classification.csv")
unified_path  <- file.path(int_dir, "metadata/unified_metadata.csv")

N_TOP_GENES   <- 3000
MIN_TEST_N    <- 10
ALPHA         <- 0.5

cat("=== 78: Sex-Stratified Fibrosis Classifiers ===\n")
cat("Started:", as.character(Sys.time()), "\n\n")

# =============================================================================
# 1. Load data
# =============================================================================
cat("--- Loading data ---\n")

# DGE -> TMM logCPM
cat("  Loading DGE ...\n")
dge <- readRDS(dge_path)
dge <- calcNormFactors(dge, method = "TMM")
logcpm <- cpm(dge, log = TRUE, prior.count = 1)
cat("  logCPM:", nrow(logcpm), "genes x", ncol(logcpm), "samples\n")

# Metadata
cat("  Loading modeling metadata ...\n")
meta <- fread(meta_path)
cat("  Metadata:", nrow(meta), "rows\n")

# --- Resolve sex column ------------------------------------------------------
# modeling_metadata has sex as "M"/"F" or empty. Harmonize to "M"/"F".
harmonize_sex <- function(x) {
  x <- toupper(substr(trimws(as.character(x)), 1, 1))
  x[!x %in% c("F", "M")] <- NA_character_
  x
}

if ("sex" %in% names(meta)) {
  meta[, sex := harmonize_sex(sex)]
}

n_sex_available <- sum(!is.na(meta$sex) & meta$sex != "", na.rm = TRUE)
cat("  Samples with sex in modeling metadata:", n_sex_available, "\n")

# If sex is largely missing, merge from unified_metadata
if (n_sex_available < nrow(meta) * 0.3) {
  cat("  Sex largely missing — merging from unified_metadata.csv ...\n")
  if (file.exists(unified_path)) {
    unif <- fread(unified_path, select = c("sample_id", "sex"))
    unif[, sex := harmonize_sex(sex)]
    if ("sex" %in% names(meta)) meta[, sex := NULL]
    meta <- merge(meta, unif, by = "sample_id", all.x = TRUE)
    n_sex_available <- sum(!is.na(meta$sex) & meta$sex != "", na.rm = TRUE)
    cat("  After merge:", n_sex_available, "samples with sex\n")
  } else {
    cat("  WARNING: unified_metadata.csv not found at", unified_path, "\n")
  }
}

# Also try meta_matched.rds if available (inferred_sex covers imputed samples)
meta_rds_path <- file.path(res_dir, "meta_matched.rds")
if (file.exists(meta_rds_path)) {
  cat("  Loading inferred sex from meta_matched.rds ...\n")
  meta_matched <- readRDS(meta_rds_path)
  if ("inferred_sex" %in% names(meta_matched)) {
    sex_lookup <- data.table(
      sample_id = meta_matched$sample_id,
      inferred_sex = harmonize_sex(meta_matched$inferred_sex)
    )
    # Fill missing sex values
    still_na <- is.na(meta$sex) | meta$sex == ""
    if (any(still_na)) {
      idx <- match(meta$sample_id[still_na], sex_lookup$sample_id)
      meta$sex[still_na] <- sex_lookup$inferred_sex[idx]
      cat("  Filled", sum(!is.na(meta$sex[still_na])), "additional sex values from meta_matched\n")
    }
  }
}

cat("  Sex distribution (all samples):\n")
print(table(meta$sex, useNA = "always"))

# --- Load sex DEG classification (for Strategy D) ----------------------------
sex_degs <- NULL
if (file.exists(sex_deg_path)) {
  cat("  Loading sex DEG classification ...\n")
  sex_degs <- fread(sex_deg_path)
  cat("  Sex DEGs:", nrow(sex_degs), "genes\n")
  cat("  Sex classes:\n")
  print(table(sex_degs$sex_class))
} else {
  cat("  WARNING: sex_deg_classification.csv not found — Strategy D will be skipped\n")
}

# =============================================================================
# 2. Filter to fibrosis-eligible samples
# =============================================================================
cat("\n--- Filtering to fibrosis-eligible samples ---\n")

meta_valid <- meta[fib_ge3 != -1 & loco_fold_fibrosis != "excluded"]
cat("  Samples with fibrosis labels & LOCO folds:", nrow(meta_valid), "\n")

# Require sex for all strategies
meta_valid <- meta_valid[!is.na(sex) & sex %in% c("F", "M")]
cat("  After requiring sex:", nrow(meta_valid), "\n")

# Intersect with expression matrix
common_ids <- intersect(meta_valid$sample_id, colnames(logcpm))
meta_valid <- meta_valid[match(common_ids, sample_id)]
logcpm_sub <- logcpm[, common_ids, drop = FALSE]
cat("  Common samples:", length(common_ids), "\n")

y_all     <- meta_valid$fib_ge3
sex_all   <- meta_valid$sex
folds_all <- meta_valid$loco_fold_fibrosis
fold_names <- sort(unique(folds_all))
n_folds   <- length(fold_names)

cat("  Outcome (fib_ge3):", table(y_all)["0"], "neg,", table(y_all)["1"], "pos\n")
cat("  Sex:", table(sex_all)["F"], "F,", table(sex_all)["M"], "M\n")
cat("  LOCO folds (", n_folds, "):", paste(fold_names, collapse = ", "), "\n")

cat("\n  Sex x fib_ge3 x fold:\n")
for (fn in fold_names) {
  mask <- folds_all == fn
  cat("  ", fn, ": F(0=", sum(sex_all[mask] == "F" & y_all[mask] == 0),
      "/1=", sum(sex_all[mask] == "F" & y_all[mask] == 1),
      ") M(0=", sum(sex_all[mask] == "M" & y_all[mask] == 0),
      "/1=", sum(sex_all[mask] == "M" & y_all[mask] == 1), ")\n")
}

# =============================================================================
# 3. Helper functions
# =============================================================================

#' Run a single elastic net fold: feature-select, rank-transform, fit, predict
#'
#' @param logcpm_mat Full gene-by-sample logCPM matrix (all genes)
#' @param train_idx Integer indices into aligned sample dimension
#' @param test_idx Integer indices into aligned sample dimension
#' @param y_train Binary response for training samples
#' @param y_test Binary response for test samples
#' @param candidate_genes Character vector of gene IDs to restrict candidates to
#'   (NULL = use all rows of logcpm_mat). Feature selection picks top n_top from
#'   these candidates based on training variance.
#' @param extra_features_train Optional numeric matrix (n_train x p) to append
#' @param extra_features_test Optional numeric matrix (n_test x p) to append
#' @param n_top Number of top-variance genes to select from training fold
#' @return list(auroc, coefs, pred_prob, y_test, n_features, msg)
run_enet_fold <- function(logcpm_mat, train_idx, test_idx,
                          y_train, y_test,
                          candidate_genes = NULL,
                          extra_features_train = NULL,
                          extra_features_test = NULL,
                          n_top = N_TOP_GENES) {

  # Restrict to candidate genes if provided
  if (!is.null(candidate_genes)) {
    available <- intersect(candidate_genes, rownames(logcpm_mat))
    if (length(available) < 50) {
      return(list(auroc = NA_real_, coefs = data.table(), pred_prob = NULL,
                  y_test = y_test, n_features = 0,
                  msg = paste0("Only ", length(available), " candidate genes available")))
    }
    mat <- logcpm_mat[available, , drop = FALSE]
  } else {
    mat <- logcpm_mat
  }

  # Within-fold feature selection: top n_top by training variance
  train_mat <- mat[, train_idx, drop = FALSE]
  gene_vars <- apply(train_mat, 1, var)
  n_sel <- min(n_top, length(gene_vars))
  top_genes <- names(sort(gene_vars, decreasing = TRUE))[seq_len(n_sel)]

  # Rank-transform expression (within each set separately)
  expr_train_raw <- t(mat[top_genes, train_idx, drop = FALSE])
  expr_train <- apply(expr_train_raw, 2, rank) / nrow(expr_train_raw)

  expr_test_raw <- t(mat[top_genes, test_idx, drop = FALSE])
  expr_test <- apply(expr_test_raw, 2, rank) / nrow(expr_test_raw)

  # Append extra features if provided
  if (!is.null(extra_features_train) && !is.null(extra_features_test)) {
    expr_train <- cbind(expr_train, extra_features_train)
    expr_test  <- cbind(expr_test, extra_features_test)
  }

  # Inverse-frequency class weights
  tab <- table(y_train)
  w_train <- ifelse(y_train == 1,
                    length(y_train) / (2 * tab["1"]),
                    length(y_train) / (2 * tab["0"]))

  # Fit elastic net with internal 5-fold CV for lambda
  cvfit <- tryCatch(
    cv.glmnet(
      x = expr_train,
      y = y_train,
      weights = w_train,
      family = "binomial",
      alpha = ALPHA,
      nfolds = 5,
      type.measure = "auc"
    ),
    error = function(e) {
      return(NULL)
    }
  )

  if (is.null(cvfit)) {
    return(list(auroc = NA_real_, coefs = data.table(), pred_prob = NULL,
                y_test = y_test, n_features = n_sel, msg = "cv.glmnet failed"))
  }

  # Predict on test
  pred_prob <- as.numeric(predict(cvfit, newx = expr_test,
                                   s = "lambda.1se", type = "response"))

  # AUROC
  auroc_val <- tryCatch({
    roc_obj <- roc(y_test, pred_prob, quiet = TRUE)
    as.numeric(auc(roc_obj))
  }, error = function(e) NA_real_)

  # Extract non-zero coefficients
  coef_vec <- as.numeric(coef(cvfit, s = "lambda.1se"))
  coef_names <- c("(Intercept)", colnames(expr_train))
  nz_idx <- which(coef_vec != 0 & coef_names != "(Intercept)")

  coefs_dt <- data.table(
    gene = coef_names[nz_idx],
    coefficient = coef_vec[nz_idx]
  )

  list(auroc = auroc_val, coefs = coefs_dt, pred_prob = pred_prob,
       y_test = y_test, n_features = n_sel, msg = "OK")
}

#' Compute AUROC for a subset of test predictions, returning NA if insufficient
compute_subset_auroc <- function(y_test, pred_prob, subset_mask) {
  if (sum(subset_mask) < 5) return(NA_real_)
  y_sub <- y_test[subset_mask]
  p_sub <- pred_prob[subset_mask]
  if (length(unique(y_sub)) < 2) return(NA_real_)
  tryCatch(as.numeric(auc(roc(y_sub, p_sub, quiet = TRUE))),
           error = function(e) NA_real_)
}

# =============================================================================
# 4. Strategy A: Universal (all samples, top 3K genes)
# =============================================================================
cat("\n\n========================================\n")
cat("STRATEGY A: Universal (baseline)\n")
cat("========================================\n")

results_all <- list()
panel_list  <- list()

for (fi in seq_along(fold_names)) {
  fold <- fold_names[fi]
  test_idx  <- which(folds_all == fold)
  train_idx <- which(folds_all != fold)

  cat(sprintf("  Fold %d/%d: %s -- Train=%d, Test=%d\n",
              fi, n_folds, fold, length(train_idx), length(test_idx)))

  if (length(test_idx) < MIN_TEST_N) {
    cat("    SKIP: <", MIN_TEST_N, "test samples\n")
    next
  }
  if (length(unique(y_all[test_idx])) < 2) {
    cat("    SKIP: test fold has only one class\n")
    next
  }

  res <- run_enet_fold(logcpm_sub, train_idx, test_idx,
                       y_all[train_idx], y_all[test_idx])

  # Per-sex AUROC within this fold
  sex_test <- sex_all[test_idx]
  auroc_f <- compute_subset_auroc(y_all[test_idx], res$pred_prob,
                                  sex_test == "F")
  auroc_m <- compute_subset_auroc(y_all[test_idx], res$pred_prob,
                                  sex_test == "M")

  cat(sprintf("    AUROC(all)=%.3f  F=%.3f  M=%.3f  (%s)\n",
              res$auroc,
              ifelse(is.na(auroc_f), NaN, auroc_f),
              ifelse(is.na(auroc_m), NaN, auroc_m),
              res$msg))

  results_all[[length(results_all) + 1]] <- data.table(
    strategy = "A_universal",
    sex_subset = "all",
    fold = fold,
    auroc = res$auroc,
    auroc_female = auroc_f,
    auroc_male = auroc_m,
    n_train = length(train_idx),
    n_test = length(test_idx),
    n_test_female = sum(sex_test == "F"),
    n_test_male = sum(sex_test == "M"),
    n_features = res$n_features
  )
}

# =============================================================================
# 5. Strategy B: Sex-specific (separate classifiers)
# =============================================================================
cat("\n\n========================================\n")
cat("STRATEGY B: Sex-specific classifiers\n")
cat("========================================\n")

for (sx in c("F", "M")) {
  cat(sprintf("\n--- Sex: %s ---\n", sx))

  sex_mask <- sex_all == sx
  y_sex    <- y_all[sex_mask]
  folds_sex <- folds_all[sex_mask]
  full_col_idx <- which(sex_mask)

  for (fi in seq_along(fold_names)) {
    fold <- fold_names[fi]
    test_mask_local  <- folds_sex == fold
    train_mask_local <- folds_sex != fold

    test_full_idx  <- full_col_idx[test_mask_local]
    train_full_idx <- full_col_idx[train_mask_local]

    cat(sprintf("  Fold %d/%d: %s -- Train=%d, Test=%d\n",
                fi, n_folds, fold, length(train_full_idx), length(test_full_idx)))

    if (length(test_full_idx) < MIN_TEST_N) {
      cat("    SKIP: <", MIN_TEST_N, "test samples\n")
      next
    }
    if (length(unique(y_sex[test_mask_local])) < 2) {
      cat("    SKIP: test fold has only one class\n")
      next
    }
    if (length(unique(y_sex[train_mask_local])) < 2) {
      cat("    SKIP: train fold has only one class\n")
      next
    }

    res <- run_enet_fold(logcpm_sub, train_full_idx, test_full_idx,
                         y_sex[train_mask_local], y_sex[test_mask_local])

    cat(sprintf("    AUROC=%.3f (%s)\n", res$auroc, res$msg))

    results_all[[length(results_all) + 1]] <- data.table(
      strategy = "B_sex_specific",
      sex_subset = sx,
      fold = fold,
      auroc = res$auroc,
      auroc_female = if (sx == "F") res$auroc else NA_real_,
      auroc_male = if (sx == "M") res$auroc else NA_real_,
      n_train = length(train_full_idx),
      n_test = length(test_full_idx),
      n_test_female = if (sx == "F") length(test_full_idx) else 0L,
      n_test_male = if (sx == "M") length(test_full_idx) else 0L,
      n_features = res$n_features
    )

    # Record non-zero coefficients for panel output
    if (nrow(res$coefs) > 0) {
      panel_list[[length(panel_list) + 1]] <- cbind(
        data.table(strategy = "B_sex_specific", sex = sx, fold = fold),
        res$coefs
      )
    }
  }
}

# =============================================================================
# 6. Strategy C: Sex-as-feature (all samples, + sex binary)
# =============================================================================
cat("\n\n========================================\n")
cat("STRATEGY C: Sex-as-feature\n")
cat("========================================\n")

# Binary sex indicator: Female=0, Male=1
sex_numeric <- as.numeric(sex_all == "M")

for (fi in seq_along(fold_names)) {
  fold <- fold_names[fi]
  test_idx  <- which(folds_all == fold)
  train_idx <- which(folds_all != fold)

  cat(sprintf("  Fold %d/%d: %s -- Train=%d, Test=%d\n",
              fi, n_folds, fold, length(train_idx), length(test_idx)))

  if (length(test_idx) < MIN_TEST_N) {
    cat("    SKIP: <", MIN_TEST_N, "test samples\n")
    next
  }
  if (length(unique(y_all[test_idx])) < 2) {
    cat("    SKIP: test fold has only one class\n")
    next
  }

  # Sex feature matrices
  sex_train <- matrix(sex_numeric[train_idx], ncol = 1,
                      dimnames = list(NULL, "sex_male"))
  sex_test  <- matrix(sex_numeric[test_idx], ncol = 1,
                      dimnames = list(NULL, "sex_male"))

  res <- run_enet_fold(logcpm_sub, train_idx, test_idx,
                       y_all[train_idx], y_all[test_idx],
                       extra_features_train = sex_train,
                       extra_features_test  = sex_test)

  # Per-sex AUROC within this fold
  sex_test_vals <- sex_all[test_idx]
  auroc_f <- compute_subset_auroc(y_all[test_idx], res$pred_prob,
                                  sex_test_vals == "F")
  auroc_m <- compute_subset_auroc(y_all[test_idx], res$pred_prob,
                                  sex_test_vals == "M")

  cat(sprintf("    AUROC(all)=%.3f  F=%.3f  M=%.3f  (%s)\n",
              res$auroc,
              ifelse(is.na(auroc_f), NaN, auroc_f),
              ifelse(is.na(auroc_m), NaN, auroc_m),
              res$msg))

  results_all[[length(results_all) + 1]] <- data.table(
    strategy = "C_sex_as_feature",
    sex_subset = "all",
    fold = fold,
    auroc = res$auroc,
    auroc_female = auroc_f,
    auroc_male = auroc_m,
    n_train = length(train_idx),
    n_test = length(test_idx),
    n_test_female = sum(sex_test_vals == "F"),
    n_test_male = sum(sex_test_vals == "M"),
    n_features = res$n_features
  )
}

# =============================================================================
# 7. Strategy D: Sex-specific gene sets
# =============================================================================
cat("\n\n========================================\n")
cat("STRATEGY D: Sex-specific gene sets\n")
cat("========================================\n")

if (!is.null(sex_degs)) {
  # Identify female-specific and male-specific gene IDs
  female_genes <- sex_degs[sex_class %in% c("Female_biased", "Female_specific"), gene]
  male_genes   <- sex_degs[sex_class %in% c("Male_biased", "Male_specific"), gene]
  cat("  Female-biased DEGs:", length(female_genes), "\n")
  cat("  Male-biased DEGs:", length(male_genes), "\n")

  # Check overlap with expression matrix
  female_genes_avail <- intersect(female_genes, rownames(logcpm_sub))
  male_genes_avail   <- intersect(male_genes, rownames(logcpm_sub))
  cat("  Female DEGs in expression:", length(female_genes_avail), "\n")
  cat("  Male DEGs in expression:", length(male_genes_avail), "\n")

  for (sx in c("F", "M")) {
    cat(sprintf("\n--- Sex: %s (using %s-specific DEGs) ---\n", sx,
                ifelse(sx == "F", "female", "male")))

    candidate_genes <- if (sx == "F") female_genes_avail else male_genes_avail

    sex_mask <- sex_all == sx
    y_sex    <- y_all[sex_mask]
    folds_sex <- folds_all[sex_mask]
    full_col_idx <- which(sex_mask)

    for (fi in seq_along(fold_names)) {
      fold <- fold_names[fi]
      test_mask_local  <- folds_sex == fold
      train_mask_local <- folds_sex != fold

      test_full_idx  <- full_col_idx[test_mask_local]
      train_full_idx <- full_col_idx[train_mask_local]

      cat(sprintf("  Fold %d/%d: %s -- Train=%d, Test=%d\n",
                  fi, n_folds, fold, length(train_full_idx), length(test_full_idx)))

      if (length(test_full_idx) < MIN_TEST_N) {
        cat("    SKIP: <", MIN_TEST_N, "test samples\n")
        next
      }
      if (length(unique(y_sex[test_mask_local])) < 2) {
        cat("    SKIP: test fold has only one class\n")
        next
      }
      if (length(unique(y_sex[train_mask_local])) < 2) {
        cat("    SKIP: train fold has only one class\n")
        next
      }

      res <- run_enet_fold(logcpm_sub, train_full_idx, test_full_idx,
                           y_sex[train_mask_local], y_sex[test_mask_local],
                           candidate_genes = candidate_genes)

      cat(sprintf("    AUROC=%.3f (%s)\n", res$auroc, res$msg))

      results_all[[length(results_all) + 1]] <- data.table(
        strategy = "D_sex_gene_sets",
        sex_subset = sx,
        fold = fold,
        auroc = res$auroc,
        auroc_female = if (sx == "F") res$auroc else NA_real_,
        auroc_male = if (sx == "M") res$auroc else NA_real_,
        n_train = length(train_full_idx),
        n_test = length(test_full_idx),
        n_test_female = if (sx == "F") length(test_full_idx) else 0L,
        n_test_male = if (sx == "M") length(test_full_idx) else 0L,
        n_features = res$n_features
      )

      # Record non-zero coefficients for panel output
      if (nrow(res$coefs) > 0) {
        panel_list[[length(panel_list) + 1]] <- cbind(
          data.table(strategy = "D_sex_gene_sets", sex = sx, fold = fold),
          res$coefs
        )
      }
    }
  }
} else {
  cat("  SKIPPED: sex_deg_classification.csv not found\n")
}

# =============================================================================
# 8. Write outputs
# =============================================================================
cat("\n\n========================================\n")
cat("Writing outputs\n")
cat("========================================\n")

# --- sex_stratified_results.csv (per-fold) ---
results_dt <- rbindlist(results_all, fill = TRUE)
out_results <- file.path(out_dir, "sex_stratified_results.csv")
fwrite(results_dt, out_results)
cat("  Wrote:", out_results, "\n")
cat("  Rows:", nrow(results_dt), "\n")

# --- sex_specific_panels.csv (top genes per sex from B and D) ---
if (length(panel_list) > 0) {
  panels_dt <- rbindlist(panel_list, fill = TRUE)

  # Compute frequency across folds per strategy x sex x gene
  panel_summary <- panels_dt[, .(
    n_folds_selected = .N,
    mean_coefficient = mean(coefficient),
    mean_abs_coefficient = mean(abs(coefficient))
  ), by = .(strategy, sex, gene)]
  panel_summary <- panel_summary[order(strategy, sex, -n_folds_selected,
                                       -mean_abs_coefficient)]

  out_panels <- file.path(out_dir, "sex_specific_panels.csv")
  fwrite(panel_summary, out_panels)
  cat("  Wrote:", out_panels, "\n")
  cat("  Rows:", nrow(panel_summary), "\n")
} else {
  cat("  WARNING: No non-zero coefficients collected for panels.\n")
  panel_summary <- data.table(strategy = character(), sex = character(),
                              gene = character(), n_folds_selected = integer(),
                              mean_coefficient = numeric(),
                              mean_abs_coefficient = numeric())
  fwrite(panel_summary, file.path(out_dir, "sex_specific_panels.csv"))
}

# --- sex_auroc_comparison.csv (summary) ---
summary_dt <- results_dt[!is.na(auroc),
  .(mean_auroc      = mean(auroc, na.rm = TRUE),
    sd_auroc        = if (.N > 1) sd(auroc, na.rm = TRUE) else NA_real_,
    min_auroc       = min(auroc, na.rm = TRUE),
    max_auroc       = max(auroc, na.rm = TRUE),
    mean_auroc_female = mean(auroc_female, na.rm = TRUE),
    sd_auroc_female   = if (sum(!is.na(auroc_female)) > 1)
                          sd(auroc_female, na.rm = TRUE) else NA_real_,
    mean_auroc_male   = mean(auroc_male, na.rm = TRUE),
    sd_auroc_male     = if (sum(!is.na(auroc_male)) > 1)
                          sd(auroc_male, na.rm = TRUE) else NA_real_,
    n_folds         = .N,
    mean_n_train    = mean(n_train),
    mean_n_test     = mean(n_test)),
  by = .(strategy, sex_subset)]

out_summary <- file.path(out_dir, "sex_auroc_comparison.csv")
fwrite(summary_dt, out_summary)
cat("  Wrote:", out_summary, "\n")

# --- Print summary table ---
cat("\n=== AUROC Summary ===\n")
for (i in seq_len(nrow(summary_dt))) {
  r <- summary_dt[i]
  cat(sprintf("  %-25s [%s]  AUROC=%.3f (sd=%.3f)  F=%.3f  M=%.3f  [%d folds]\n",
    r$strategy, r$sex_subset,
    r$mean_auroc,
    ifelse(is.na(r$sd_auroc), 0, r$sd_auroc),
    ifelse(is.na(r$mean_auroc_female) | is.nan(r$mean_auroc_female), NA, r$mean_auroc_female),
    ifelse(is.na(r$mean_auroc_male) | is.nan(r$mean_auroc_male), NA, r$mean_auroc_male),
    r$n_folds))
}

# --- Panel overlap analysis ---
if (nrow(panel_summary) > 0) {
  cat("\n=== Gene Panel Overlap (B: Female vs Male) ===\n")
  b_female_genes <- unique(panel_summary[strategy == "B_sex_specific" & sex == "F", gene])
  b_male_genes   <- unique(panel_summary[strategy == "B_sex_specific" & sex == "M", gene])
  b_shared       <- intersect(b_female_genes, b_male_genes)
  cat("  Female-model genes:", length(b_female_genes), "\n")
  cat("  Male-model genes:", length(b_male_genes), "\n")
  cat("  Shared:", length(b_shared), "\n")
  cat("  Female-unique:", length(setdiff(b_female_genes, b_male_genes)), "\n")
  cat("  Male-unique:", length(setdiff(b_male_genes, b_female_genes)), "\n")

  # Print top 10 per sex
  for (sx in c("F", "M")) {
    top10 <- head(panel_summary[strategy == "B_sex_specific" & sex == sx], 10)
    if (nrow(top10) > 0) {
      cat(sprintf("\n  Top 10 genes for %s (Strategy B):\n",
                  ifelse(sx == "F", "female", "male")))
      for (j in seq_len(nrow(top10))) {
        cat(sprintf("    %2d. %s  (folds=%d, mean|coef|=%.4f)\n",
                    j, top10$gene[j], top10$n_folds_selected[j],
                    top10$mean_abs_coefficient[j]))
      }
    }
  }

  if (any(panel_summary$strategy == "D_sex_gene_sets")) {
    cat("\n=== Gene Panel Overlap (D: Female vs Male) ===\n")
    d_female_genes <- unique(panel_summary[strategy == "D_sex_gene_sets" & sex == "F", gene])
    d_male_genes   <- unique(panel_summary[strategy == "D_sex_gene_sets" & sex == "M", gene])
    d_shared       <- intersect(d_female_genes, d_male_genes)
    cat("  Female-model genes:", length(d_female_genes), "\n")
    cat("  Male-model genes:", length(d_male_genes), "\n")
    cat("  Shared:", length(d_shared), "\n")
  }
}

cat("\nCompleted:", as.character(Sys.time()), "\n")
cat("=== Done ===\n")
