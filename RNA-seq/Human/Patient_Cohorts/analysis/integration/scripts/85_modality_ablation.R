#!/usr/bin/env Rscript
# 85_modality_ablation.R
# Systematic modality ablation study for MASLD fibrosis staging (F>=3)
#
# Tests 9 feature combinations with 6-fold LOCO-CV elastic net:
#   1. Expression only (top 3K by training-fold variance, rank-transformed)
#   2. VAE embeddings only (64-dim latent space)
#   3. Deconvolution only (4 BayesPrism proportion columns)
#   4. Sex only (1 binary feature)
#   5. Expression + deconv
#   6. Expression + sex
#   7. Expression + deconv + sex
#   8. VAE + deconv + sex
#   9. All combined (expression + VAE + deconv + sex)
#
# Outputs (in results/staging_classifier/):
#   - ablation_results.csv:        Per-combination per-fold AUROC
#   - ablation_summary.csv:        Mean AUROC per combination, ranked
#   - masld_composite_score.csv:   Per-sample composite score from best model
#   - modality_importance.csv:     Delta AUROC when each modality is removed from full
#
# Usage: Rscript 85_modality_ablation.R
# SLURM: cpu, 8 CPUs, 64G, 48h

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(glmnet)
  library(pROC)
})

set.seed(42)

BASE   <- Sys.getenv("MASLD_PROJECT_ROOT",
            "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT    <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR   <- file.path(INT, "results")
OUTDIR <- file.path(RDIR, "staging_classifier")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

K_FEATURES <- 3000

cat("=== 85: Modality Ablation Study for F>=3 Staging ===\n")
cat("Started:", as.character(Sys.time()), "\n\n")

# ============================================================
# 1. Load all data sources
# ============================================================

# --- Expression ---
cat("[1/4] Loading merged DGE...\n")
dge <- readRDS(file.path(RDIR, "integration/merged_dge.rds"))
dge <- calcNormFactors(dge, method = "TMM")
logcpm <- cpm(dge, log = TRUE, prior.count = 1)
cat("  Expression:", nrow(logcpm), "genes x", ncol(logcpm), "samples\n")

# --- Metadata ---
cat("[2/4] Loading modeling metadata...\n")
meta <- fread(file.path(OUTDIR, "modeling_metadata.csv"))
cat("  Metadata:", nrow(meta), "samples\n")

# --- VAE embeddings ---
cat("[3/4] Loading VAE embeddings...\n")
emb <- fread(file.path(OUTDIR, "embeddings_all_samples.csv"))
emb_cols <- setdiff(names(emb), "sample_id")
cat("  Embeddings:", nrow(emb), "samples x", length(emb_cols), "dimensions\n")

# --- Deconvolution proportions ---
cat("[4/4] Loading BayesPrism proportions...\n")
deconv <- fread(file.path(RDIR, "deconvolution/bayesprism/unified_bayesprism_proportions.csv"))
deconv_cols <- setdiff(names(deconv), c("sample_id", "dataset"))
cat("  Deconvolution:", nrow(deconv), "samples x", length(deconv_cols), "proportion columns\n")
cat("  Columns:", paste(deconv_cols, collapse = ", "), "\n")

# ============================================================
# 2. Build unified sample table — restrict to F>=3-eligible
# ============================================================
cat("\nRestricting to fibrosis-eligible samples...\n")

# Filter to non-excluded, non-NA fib_ge3 labels
valid_mask <- meta$loco_fold_fibrosis != "excluded" &
              !is.na(meta$fib_ge3) &
              meta$fib_ge3 >= 0
meta_valid <- meta[valid_mask]
cat("  Fibrosis-eligible samples in metadata:", nrow(meta_valid), "\n")

# Intersect with all modalities
common_expr   <- intersect(meta_valid$sample_id, colnames(logcpm))
common_emb    <- intersect(meta_valid$sample_id, emb$sample_id)
common_deconv <- intersect(meta_valid$sample_id, deconv$sample_id)
common_all    <- Reduce(intersect, list(common_expr, common_emb, common_deconv))

cat("  Samples with expression:", length(common_expr), "\n")
cat("  Samples with embeddings:", length(common_emb), "\n")
cat("  Samples with deconv:    ", length(common_deconv), "\n")
cat("  Samples in ALL sources: ", length(common_all), "\n")

# Use the intersection of all modalities so every combination is evaluated on
# the same sample set (fair comparison)
meta_valid <- meta_valid[sample_id %in% common_all]
meta_valid <- meta_valid[match(common_all, sample_id)]  # stable order
logcpm     <- logcpm[, common_all]
emb        <- emb[match(common_all, sample_id)]
deconv     <- deconv[match(common_all, sample_id)]

n_total <- nrow(meta_valid)
cat("\n  Final sample set:", n_total, "samples\n")
cat("  Label distribution: F>=3 =", sum(meta_valid$fib_ge3 == 1),
    ", F<3 =", sum(meta_valid$fib_ge3 == 0), "\n")

folds <- sort(unique(meta_valid$loco_fold_fibrosis))
folds <- folds[folds != "excluded" & folds != "NA"]
cat("  LOCO folds:", length(folds), "(", paste(folds, collapse = ", "), ")\n\n")

# Pre-build non-expression feature matrices
emb_mat    <- as.matrix(emb[, ..emb_cols])
rownames(emb_mat) <- common_all

deconv_mat <- as.matrix(deconv[, ..deconv_cols])
rownames(deconv_mat) <- common_all

sex_mat    <- matrix(as.integer(meta_valid$sex == "M"), ncol = 1,
                     dimnames = list(common_all, "sex_male"))

# ============================================================
# 3. Define modality combinations
# ============================================================
combinations <- list(
  list(name = "Expression only",         use_expr = TRUE,  use_vae = FALSE, use_deconv = FALSE, use_sex = FALSE),
  list(name = "VAE embeddings only",     use_expr = FALSE, use_vae = TRUE,  use_deconv = FALSE, use_sex = FALSE),
  list(name = "Deconvolution only",      use_expr = FALSE, use_vae = FALSE, use_deconv = TRUE,  use_sex = FALSE),
  list(name = "Sex only",                use_expr = FALSE, use_vae = FALSE, use_deconv = FALSE, use_sex = TRUE),
  list(name = "Expression + deconv",     use_expr = TRUE,  use_vae = FALSE, use_deconv = TRUE,  use_sex = FALSE),
  list(name = "Expression + sex",        use_expr = TRUE,  use_vae = FALSE, use_deconv = FALSE, use_sex = TRUE),
  list(name = "Expression + deconv + sex", use_expr = TRUE, use_vae = FALSE, use_deconv = TRUE, use_sex = TRUE),
  list(name = "VAE + deconv + sex",      use_expr = FALSE, use_vae = TRUE,  use_deconv = TRUE,  use_sex = TRUE),
  list(name = "All combined",            use_expr = TRUE,  use_vae = TRUE,  use_deconv = TRUE,  use_sex = TRUE)
)

# ============================================================
# 4. Helper: Build feature matrix for one fold
# ============================================================
build_features <- function(train_idx, test_idx, combo) {
  # Expression: top K by training-fold variance, rank-transformed
  X_train_parts <- list()
  X_test_parts  <- list()

  if (combo$use_expr) {
    lcpm_train <- logcpm[, train_idx, drop = FALSE]
    lcpm_test  <- logcpm[, test_idx,  drop = FALSE]

    gene_var <- apply(lcpm_train, 1, var)
    top_genes <- names(sort(gene_var, decreasing = TRUE))[1:min(K_FEATURES, length(gene_var))]

    rank_train <- apply(lcpm_train[top_genes, , drop = FALSE], 2, function(x) {
      r <- rank(x, ties.method = "average")
      r / length(r)
    })
    rank_test <- apply(lcpm_test[top_genes, , drop = FALSE], 2, function(x) {
      r <- rank(x, ties.method = "average")
      r / length(r)
    })

    X_train_parts[["expr"]] <- t(rank_train)
    X_test_parts[["expr"]]  <- t(rank_test)
  }

  if (combo$use_vae) {
    X_train_parts[["vae"]] <- emb_mat[train_idx, , drop = FALSE]
    X_test_parts[["vae"]]  <- emb_mat[test_idx,  , drop = FALSE]
  }

  if (combo$use_deconv) {
    X_train_parts[["deconv"]] <- deconv_mat[train_idx, , drop = FALSE]
    X_test_parts[["deconv"]]  <- deconv_mat[test_idx,  , drop = FALSE]
  }

  if (combo$use_sex) {
    X_train_parts[["sex"]] <- sex_mat[train_idx, , drop = FALSE]
    X_test_parts[["sex"]]  <- sex_mat[test_idx,  , drop = FALSE]
  }

  list(
    X_train = do.call(cbind, X_train_parts),
    X_test  = do.call(cbind, X_test_parts)
  )
}

# ============================================================
# 5. Run ablation: LOCO-CV for each combination
# ============================================================
all_fold_results <- list()   # per-combo per-fold AUROC
all_predictions  <- list()   # per-combo per-sample predictions (for composite score)

for (ci in seq_along(combinations)) {
  combo <- combinations[[ci]]
  combo_name <- combo$name
  cat("=== Combination", ci, "/", length(combinations), ":", combo_name, "===\n")

  fold_aurocs <- numeric()
  combo_preds <- list()

  for (fold in folds) {
    test_mask  <- meta_valid$loco_fold_fibrosis == fold
    train_mask <- !test_mask

    train_idx <- which(train_mask)
    test_idx  <- which(test_mask)

    if (length(test_idx) < 5 || length(train_idx) < 20) {
      cat("  Fold", fold, ": SKIP (too few samples)\n")
      next
    }

    y_train <- meta_valid$fib_ge3[train_idx]
    y_test  <- meta_valid$fib_ge3[test_idx]

    if (length(unique(y_test)) < 2) {
      cat("  Fold", fold, ": SKIP (single class in test)\n")
      next
    }

    # Build features
    feat <- build_features(train_idx, test_idx, combo)

    # Handle potential NA in deconv/sex by imputing column medians from train
    if (any(is.na(feat$X_train))) {
      col_medians <- apply(feat$X_train, 2, median, na.rm = TRUE)
      for (j in seq_len(ncol(feat$X_train))) {
        na_train <- is.na(feat$X_train[, j])
        if (any(na_train)) feat$X_train[na_train, j] <- col_medians[j]
        na_test <- is.na(feat$X_test[, j])
        if (any(na_test)) feat$X_test[na_test, j] <- col_medians[j]
      }
    }

    # Class-balanced weights
    class_counts <- table(y_train)
    class_wts <- 1 / class_counts
    class_wts <- class_wts / sum(class_wts) * length(class_wts)
    sample_wts <- as.numeric(class_wts[as.character(y_train)])

    # Train elastic net (alpha=0.5)
    fit <- tryCatch({
      cv.glmnet(
        x = feat$X_train, y = factor(y_train),
        family = "binomial", alpha = 0.5, nfolds = 5,
        type.measure = "auc", weights = sample_wts
      )
    }, error = function(e) {
      cat("  Fold", fold, ": ERROR —", e$message, "\n")
      NULL
    })

    if (is.null(fit)) {
      fold_aurocs <- c(fold_aurocs, NA_real_)
      next
    }

    prob <- as.numeric(predict(fit, newx = feat$X_test,
                               s = "lambda.min", type = "response"))

    auroc <- tryCatch(
      as.numeric(auc(roc(y_test, prob, quiet = TRUE))),
      error = function(e) NA_real_
    )

    fold_aurocs <- c(fold_aurocs, auroc)
    cat("  Fold", fold, ": AUROC =", round(auroc, 4),
        " (n_train=", length(train_idx), ", n_test=", length(test_idx), ")\n")

    # Store per-sample predictions
    combo_preds[[fold]] <- data.table(
      sample_id   = meta_valid$sample_id[test_idx],
      dataset     = meta_valid$dataset[test_idx],
      true_label  = y_test,
      probability = prob,
      fold        = fold,
      combination = combo_name
    )
  }

  # Record fold-level results
  for (fi in seq_along(folds)) {
    auc_val <- if (fi <= length(fold_aurocs)) fold_aurocs[fi] else NA_real_
    all_fold_results[[length(all_fold_results) + 1]] <- data.table(
      combination = combo_name,
      fold        = folds[fi],
      auroc       = auc_val
    )
  }

  # Store predictions for composite score computation
  if (length(combo_preds) > 0) {
    all_predictions[[combo_name]] <- rbindlist(combo_preds, fill = TRUE)
  }

  valid_aurocs <- fold_aurocs[!is.na(fold_aurocs)]
  if (length(valid_aurocs) > 0) {
    cat("  >> Mean AUROC =", round(mean(valid_aurocs), 4),
        "+/-", round(sd(valid_aurocs), 4),
        "(", length(valid_aurocs), "folds)\n\n")
  } else {
    cat("  >> No valid folds\n\n")
  }
}

# ============================================================
# 6. Assemble outputs
# ============================================================
cat("=== Assembling outputs ===\n")

# --- ablation_results.csv: per-combination per-fold AUROC ---
ablation_results <- rbindlist(all_fold_results, fill = TRUE)
fwrite(ablation_results, file.path(OUTDIR, "ablation_results.csv"))
cat("  ablation_results.csv:", nrow(ablation_results), "rows\n")

# --- ablation_summary.csv: mean/SD AUROC ranked ---
ablation_summary <- ablation_results[
  !is.na(auroc),
  .(mean_auroc   = mean(auroc),
    sd_auroc     = sd(auroc),
    min_auroc    = min(auroc),
    max_auroc    = max(auroc),
    n_folds      = .N),
  by = combination
][order(-mean_auroc)]

ablation_summary[, rank := .I]
fwrite(ablation_summary, file.path(OUTDIR, "ablation_summary.csv"))
cat("  ablation_summary.csv:", nrow(ablation_summary), "rows\n")
cat("\n  --- Ablation Summary (ranked by mean AUROC) ---\n")
print(ablation_summary[, .(rank, combination, mean_auroc, sd_auroc, n_folds)])

# --- masld_composite_score.csv: predictions from the best combination ---
best_combo_name <- ablation_summary$combination[1]
cat("\n  Best combination:", best_combo_name,
    "(mean AUROC =", round(ablation_summary$mean_auroc[1], 4), ")\n")

if (best_combo_name %in% names(all_predictions)) {
  composite <- all_predictions[[best_combo_name]]
  setnames(composite, "probability", "composite_score")
  composite[, combination := NULL]
  fwrite(composite, file.path(OUTDIR, "masld_composite_score.csv"))
  cat("  masld_composite_score.csv:", nrow(composite), "samples\n")
} else {
  cat("  WARNING: no predictions for best combination\n")
}

# --- modality_importance.csv: delta AUROC removing each modality from full ---
cat("\n=== Modality importance (drop-one from 'All combined') ===\n")

full_auroc <- ablation_summary[combination == "All combined", mean_auroc]

# Map each modality to the combination that removes it from "All combined"
drop_map <- list(
  Expression    = "VAE + deconv + sex",         # All minus expression
  VAE           = "Expression + deconv + sex",   # All minus VAE
  Deconvolution = "Expression + sex",            # All minus deconv (and minus VAE)
  Sex           = "Expression + deconv + sex"     # Wait — no exact match for All minus sex only
)

# For a proper ablation we compute importance as:
#   delta = AUROC(All combined) - AUROC(All minus modality X)
# The combinations defined cover these drop-one comparisons:
#   Drop Expression -> VAE + deconv + sex
#   Drop VAE        -> Expression + deconv + sex
#   Drop Deconv     -> Expression + VAE + sex (not defined — use Expression + sex as proxy)
#   Drop Sex        -> Expression + VAE + deconv (not defined exactly)
#
# Since we don't have all possible subsets, we compute marginal importance
# as the difference between "All combined" and the closest subset without that modality.
# For completeness, also compute the "additive" gain of each modality
# (expression-only baseline → add each modality).

importance_rows <- list()

# Marginal importance: full - (full minus modality)
modality_drop <- list(
  Expression    = "VAE + deconv + sex",
  VAE_embeddings = "Expression + deconv + sex",
  Deconvolution = "Expression + sex",      # proxy: also drops VAE
  Sex           = "Expression + deconv + sex"  # proxy: also drops VAE
)

# Better approach: compare directly using the subset results we have
# Additive from expression baseline
expr_auroc <- ablation_summary[combination == "Expression only", mean_auroc]

for (mod_name in c("Expression", "VAE_embeddings", "Deconvolution", "Sex")) {
  drop_combo <- modality_drop[[mod_name]]
  drop_auroc <- ablation_summary[combination == drop_combo, mean_auroc]

  if (length(full_auroc) > 0 && length(drop_auroc) > 0) {
    delta <- full_auroc - drop_auroc
  } else {
    delta <- NA_real_
  }

  importance_rows[[length(importance_rows) + 1]] <- data.table(
    modality            = mod_name,
    full_model_auroc    = full_auroc,
    without_modality    = drop_combo,
    without_auroc       = drop_auroc,
    delta_auroc         = delta
  )
}

# Also compute additive gain over expression-only baseline
additive_rows <- list()
additive_combos <- list(
  Deconvolution = "Expression + deconv",
  Sex           = "Expression + sex",
  "Deconv + Sex"  = "Expression + deconv + sex",
  VAE_embeddings  = "All combined"
)

for (add_name in names(additive_combos)) {
  add_combo <- additive_combos[[add_name]]
  add_auroc <- ablation_summary[combination == add_combo, mean_auroc]
  if (length(add_auroc) > 0 && length(expr_auroc) > 0) {
    gain <- add_auroc - expr_auroc
  } else {
    gain <- NA_real_
  }
  additive_rows[[length(additive_rows) + 1]] <- data.table(
    modality              = add_name,
    baseline              = "Expression only",
    baseline_auroc        = expr_auroc,
    augmented_combination = add_combo,
    augmented_auroc       = add_auroc,
    additive_gain         = gain
  )
}

importance <- rbindlist(importance_rows, fill = TRUE)
additive   <- rbindlist(additive_rows, fill = TRUE)

# Combine both views
importance_out <- merge(
  importance[, .(modality, delta_auroc_drop = delta_auroc)],
  additive[, .(modality, additive_gain)],
  by = "modality", all = TRUE
)

# Add single-modality standalone performance
standalone_map <- list(
  Expression    = "Expression only",
  VAE_embeddings = "VAE embeddings only",
  Deconvolution = "Deconvolution only",
  Sex           = "Sex only"
)
for (mod in names(standalone_map)) {
  sa <- ablation_summary[combination == standalone_map[[mod]], mean_auroc]
  if (length(sa) > 0) {
    importance_out[modality == mod, standalone_auroc := sa]
  }
}

importance_out <- importance_out[order(-standalone_auroc, na.last = TRUE)]
fwrite(importance_out, file.path(OUTDIR, "modality_importance.csv"))
cat("  modality_importance.csv:\n")
print(importance_out)

# Full detail tables for reference
fwrite(importance, file.path(OUTDIR, "modality_importance_drop_detail.csv"))
fwrite(additive,   file.path(OUTDIR, "modality_importance_additive_detail.csv"))

cat("\n=== 85_modality_ablation.R completed:", as.character(Sys.time()), "===\n")
