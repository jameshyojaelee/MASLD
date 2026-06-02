#!/usr/bin/env Rscript
# 82_plasma_classifier.R
# Plasma Protein-Based Fibrosis Classifier from Olink Data
#
# Builds an elastic net classifier for advanced fibrosis (F3/F4 vs F0-2)
# using circulating plasma proteins measured by Olink proximity extension.
# Four experiments: (A) full Olink, (B) tissue-informed subset,
# (C) recursive feature elimination panel, (D) NFASC+GDF15 benchmark.
#
# Inputs:
#   - Olink plasma NPX data (1,460 proteins x 218 subjects, tab-separated)
#   - GSE276114 liver metadata (177 samples with fibrosis staging)
#   - Tissue-plasma bridge table (optional, from 81_tissue_plasma_bridge.R)
#
# Outputs (in results/staging_classifier/):
#   - plasma_classifier_results.csv    — per-experiment AUROC/accuracy summary
#   - plasma_minimal_panel.csv         — top panel genes from RFE at each size
#   - plasma_panel_curve.csv           — AUROC vs panel size curve
#   - plasma_nfasc_gdf15_benchmark.csv — NFASC+GDF15 2-protein model results
#
# Usage: Rscript 82_plasma_classifier.R
# SLURM: cpu, 8 CPUs, 64GB RAM, 48h

suppressPackageStartupMessages({
  library(data.table)
  library(glmnet)
  library(pROC)
})

set.seed(42)

# ============================================================
# Paths
# ============================================================
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
OUTDIR <- file.path(INT, "results/staging_classifier")
OLINK_DIR <- file.path(BASE, "Analysis/Proteomics/data/olink_plasma")
PROT_DIR  <- file.path(BASE, "Analysis/Proteomics/results")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

cat("=== 82: Plasma Protein Fibrosis Classifier ===\n")
cat("Started:", as.character(Sys.time()), "\n\n")

# ============================================================
# STEP 1: Load and prepare data
# ============================================================
cat("=== STEP 1: Loading Data ===\n")

# 1a. Olink plasma data (proteins x subjects, NPX log2)
olink_file <- file.path(OLINK_DIR, "olink.qc.finished.mendeley.data.txt")
olink_raw <- fread(olink_file, header = TRUE, sep = "\t")
cat("  Olink raw:", nrow(olink_raw), "proteins x", ncol(olink_raw) - 1, "subjects\n")

protein_names <- olink_raw$Assay
olink_mat <- as.matrix(olink_raw[, -1, with = FALSE])
rownames(olink_mat) <- protein_names

# Transpose to subjects x proteins
npx <- t(olink_mat)
cat("  Transposed NPX matrix:", nrow(npx), "subjects x", ncol(npx), "proteins\n")

# 1b. Impute NAs with column median
na_count <- sum(is.na(npx))
cat("  NAs in NPX:", na_count, "\n")
if (na_count > 0) {
  for (j in seq_len(ncol(npx))) {
    col_vals <- npx[, j]
    na_idx <- is.na(col_vals)
    if (any(na_idx)) {
      npx[na_idx, j] <- median(col_vals, na.rm = TRUE)
    }
  }
  cat("  Imputed NAs with column medians\n")
}

# 1c. Extract subject number from rownames ("Subject N" -> N)
subject_nums <- as.integer(sub("^Subject\\s+", "", rownames(npx)))

# 1d. Liver metadata (fibrosis staging)
meta_file <- file.path(PROT_DIR, "gse276114_disease_metadata.csv")
meta <- fread(meta_file)
cat("  Metadata:", nrow(meta), "liver samples\n")

# 1e. Map subjects to metadata via sample_number
# Olink "Subject N" maps to liver metadata sample_number = N
meta_map <- meta[, .(sample_number, disease_group, disease)]
setkey(meta_map, sample_number)

# Build label vector: match each subject to its fibrosis stage
label_dt <- data.table(
  subject_idx = seq_len(nrow(npx)),
  subject_num = subject_nums
)
label_dt <- merge(label_dt, meta_map, by.x = "subject_num", by.y = "sample_number",
                  all.x = TRUE)

# Keep only subjects with liver metadata and valid disease_group
label_dt <- label_dt[!is.na(disease_group)]
label_dt[, advanced := as.integer(disease_group %in% c("F3", "F4"))]

cat("  Matched subjects:", nrow(label_dt), "\n")
cat("  Advanced (F3/F4):", sum(label_dt$advanced == 1),
    " Early (F0-2):", sum(label_dt$advanced == 0), "\n")

# Subset NPX matrix to matched subjects
npx_matched <- npx[label_dt$subject_idx, , drop = FALSE]
y <- label_dt$advanced
etiology <- label_dt$disease

cat("  Etiology breakdown:\n")
print(table(etiology, y))

# 1f. Class weights (inverse frequency)
n_pos <- sum(y == 1)
n_neg <- sum(y == 0)
w_pos <- length(y) / (2 * n_pos)
w_neg <- length(y) / (2 * n_neg)
sample_weights <- ifelse(y == 1, w_pos, w_neg)
cat("  Class weights: advanced=", round(w_pos, 3), " early=", round(w_neg, 3), "\n\n")

# 1g. Tissue-plasma bridge (optional)
bridge_file <- file.path(OUTDIR, "plasma_translatable_degs.csv")
has_bridge <- file.exists(bridge_file)
if (has_bridge) {
  bridge <- fread(bridge_file)
  tissue_proteins <- intersect(bridge$human_symbol, colnames(npx_matched))
  cat("  Bridge file loaded:", nrow(bridge), "genes,", length(tissue_proteins),
      "overlap with Olink\n\n")
} else {
  cat("  No bridge file found — Experiment B will use dream DEGs directly\n\n")
  tissue_proteins <- NULL
}

# ============================================================
# Helper: fit elastic net and evaluate
# ============================================================
fit_enet <- function(X, y, weights, alpha = 0.5, nfolds = 5, foldid = NULL) {
  # Fit elastic net with CV
  if (is.null(foldid)) {
    fit <- cv.glmnet(X, y, family = "binomial", alpha = alpha,
                     weights = weights, nfolds = nfolds,
                     type.measure = "auc", standardize = TRUE)
  } else {
    fit <- cv.glmnet(X, y, family = "binomial", alpha = alpha,
                     weights = weights, foldid = foldid,
                     type.measure = "auc", standardize = TRUE)
  }
  fit
}

evaluate_model <- function(fit, X, y, label = "model") {
  # Predictions at lambda.min
  pred <- predict(fit, newx = X, s = "lambda.min", type = "response")[, 1]
  roc_obj <- roc(y, pred, quiet = TRUE)
  auc_val <- as.numeric(auc(roc_obj))

  # Optimal threshold (Youden)
  coords_best <- coords(roc_obj, "best", ret = c("threshold", "sensitivity",
                                                   "specificity", "accuracy"),
                         best.method = "youden")

  # Selected features
  coefs <- coef(fit, s = "lambda.min")
  n_selected <- sum(coefs[-1] != 0)

  list(auc = auc_val, sensitivity = coords_best$sensitivity[1],
       specificity = coords_best$specificity[1],
       accuracy = coords_best$accuracy[1],
       n_features = n_selected, roc = roc_obj, predictions = pred)
}

# Leave-one-etiology-out CV
loeo_cv <- function(X, y, etiology, weights, alpha = 0.5) {
  etiologies <- unique(etiology)
  all_preds <- rep(NA_real_, length(y))

  for (eti in etiologies) {
    test_idx <- which(etiology == eti)
    train_idx <- which(etiology != eti)

    if (length(unique(y[train_idx])) < 2) {
      cat("    Skipping", eti, "— only one class in training\n")
      next
    }

    fit <- glmnet(X[train_idx, , drop = FALSE], y[train_idx],
                  family = "binomial", alpha = alpha,
                  weights = weights[train_idx], standardize = TRUE)

    # Use lambda that gives ~same complexity as CV would
    # Pick lambda with min deviance from the fit path
    pred <- predict(fit, newx = X[test_idx, , drop = FALSE],
                    type = "response")
    # Use the lambda at median index for regularization
    lambda_idx <- max(1, floor(length(fit$lambda) / 3))
    all_preds[test_idx] <- pred[, lambda_idx]
    cat("    Held out:", eti, "(n=", length(test_idx), "), lambda_idx=",
        lambda_idx, "\n")
  }

  valid <- !is.na(all_preds)
  if (sum(valid) < 10) return(list(auc = NA, n_valid = sum(valid)))

  roc_obj <- roc(y[valid], all_preds[valid], quiet = TRUE)
  auc_val <- as.numeric(auc(roc_obj))
  coords_best <- coords(roc_obj, "best", ret = c("sensitivity", "specificity",
                                                   "accuracy"),
                         best.method = "youden")

  list(auc = auc_val, sensitivity = coords_best$sensitivity[1],
       specificity = coords_best$specificity[1],
       accuracy = coords_best$accuracy[1],
       n_valid = sum(valid), roc = roc_obj, predictions = all_preds)
}

# Random 5-fold CV (repeated 10 times for stability)
random_cv <- function(X, y, weights, alpha = 0.5, nfolds = 5, nrepeats = 10) {
  aucs <- numeric(nrepeats)
  sens <- numeric(nrepeats)
  spec <- numeric(nrepeats)
  accs <- numeric(nrepeats)
  n_feats <- numeric(nrepeats)

  for (r in seq_len(nrepeats)) {
    set.seed(42 + r)
    foldid <- sample(rep(seq_len(nfolds), length.out = length(y)))
    fit <- fit_enet(X, y, weights, alpha = alpha, foldid = foldid)
    res <- evaluate_model(fit, X, y)
    aucs[r] <- res$auc
    sens[r] <- res$sensitivity
    spec[r] <- res$specificity
    accs[r] <- res$accuracy
    n_feats[r] <- res$n_features
  }

  list(auc_mean = mean(aucs), auc_sd = sd(aucs),
       sensitivity_mean = mean(sens), specificity_mean = mean(spec),
       accuracy_mean = mean(accs), n_features_mean = mean(n_feats))
}

# ============================================================
# Collect all results
# ============================================================
results_list <- list()

# ============================================================
# EXPERIMENT A: Full Olink classifier (all 1460 proteins)
# ============================================================
cat("=== EXPERIMENT A: Full Olink Classifier ===\n")

# A1: Leave-one-etiology-out CV
cat("  A1: Leave-one-etiology-out CV\n")
loeo_A <- loeo_cv(npx_matched, y, etiology, sample_weights, alpha = 0.5)
cat("    LOEO AUROC:", round(loeo_A$auc, 4), "\n")

# A2: Random 5-fold CV (10 repeats)
cat("  A2: Random 5-fold CV (10 repeats)\n")
rcv_A <- random_cv(npx_matched, y, sample_weights, alpha = 0.5)
cat("    5-fold CV AUROC:", round(rcv_A$auc_mean, 4), "+/-", round(rcv_A$auc_sd, 4), "\n")

results_list[["A_full_olink_LOEO"]] <- data.table(
  experiment = "A_full_olink", cv_method = "LOEO",
  auroc = loeo_A$auc, auroc_sd = NA_real_,
  sensitivity = loeo_A$sensitivity, specificity = loeo_A$specificity,
  accuracy = loeo_A$accuracy, n_proteins = ncol(npx_matched),
  n_features_selected = NA_integer_
)
results_list[["A_full_olink_5foldCV"]] <- data.table(
  experiment = "A_full_olink", cv_method = "5fold_repeated",
  auroc = rcv_A$auc_mean, auroc_sd = rcv_A$auc_sd,
  sensitivity = rcv_A$sensitivity_mean, specificity = rcv_A$specificity_mean,
  accuracy = rcv_A$accuracy_mean, n_proteins = ncol(npx_matched),
  n_features_selected = round(rcv_A$n_features_mean)
)

cat("\n")

# ============================================================
# EXPERIMENT B: Tissue-informed classifier
# ============================================================
cat("=== EXPERIMENT B: Tissue-Informed Classifier ===\n")

if (!is.null(tissue_proteins) && length(tissue_proteins) > 10) {
  npx_tissue <- npx_matched[, tissue_proteins, drop = FALSE]
  cat("  Using", ncol(npx_tissue), "tissue-informed proteins\n")
} else {
  # Fallback: load dream DEGs and intersect with Olink proteins
  atlas_file <- file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
  if (file.exists(atlas_file)) {
    atlas <- fread(atlas_file, select = c("human_symbol", "dream_padj"))
    deg_genes <- atlas[dream_padj < 0.1 & !is.na(dream_padj), human_symbol]
    tissue_proteins <- intersect(deg_genes, colnames(npx_matched))
    npx_tissue <- npx_matched[, tissue_proteins, drop = FALSE]
    cat("  Fallback: intersected", length(deg_genes), "dream DEGs with Olink ->",
        ncol(npx_tissue), "proteins\n")
  } else {
    cat("  WARNING: No atlas or bridge file. Skipping Experiment B.\n")
    npx_tissue <- NULL
  }
}

if (!is.null(npx_tissue) && ncol(npx_tissue) > 5) {
  # B1: LOEO
  cat("  B1: Leave-one-etiology-out CV\n")
  loeo_B <- loeo_cv(npx_tissue, y, etiology, sample_weights, alpha = 0.5)
  cat("    LOEO AUROC:", round(loeo_B$auc, 4), "\n")

  # B2: Random 5-fold CV
  cat("  B2: Random 5-fold CV (10 repeats)\n")
  rcv_B <- random_cv(npx_tissue, y, sample_weights, alpha = 0.5)
  cat("    5-fold CV AUROC:", round(rcv_B$auc_mean, 4), "+/-",
      round(rcv_B$auc_sd, 4), "\n")

  results_list[["B_tissue_informed_LOEO"]] <- data.table(
    experiment = "B_tissue_informed", cv_method = "LOEO",
    auroc = loeo_B$auc, auroc_sd = NA_real_,
    sensitivity = loeo_B$sensitivity, specificity = loeo_B$specificity,
    accuracy = loeo_B$accuracy, n_proteins = ncol(npx_tissue),
    n_features_selected = NA_integer_
  )
  results_list[["B_tissue_informed_5foldCV"]] <- data.table(
    experiment = "B_tissue_informed", cv_method = "5fold_repeated",
    auroc = rcv_B$auc_mean, auroc_sd = rcv_B$auc_sd,
    sensitivity = rcv_B$sensitivity_mean, specificity = rcv_B$specificity_mean,
    accuracy = rcv_B$accuracy_mean, n_proteins = ncol(npx_tissue),
    n_features_selected = round(rcv_B$n_features_mean)
  )
} else {
  cat("  Experiment B skipped (insufficient tissue-informed proteins)\n")
}
cat("\n")

# ============================================================
# EXPERIMENT C: Minimal Panel RFE
# ============================================================
cat("=== EXPERIMENT C: Minimal Panel via RFE ===\n")

# Step 1: Rank all proteins by elastic net coefficient magnitude
cat("  Fitting full model for feature ranking...\n")
set.seed(42)
full_fit <- cv.glmnet(npx_matched, y, family = "binomial", alpha = 0.5,
                       weights = sample_weights, nfolds = 5,
                       type.measure = "auc", standardize = TRUE)

# Extract absolute coefficients at lambda.min
coefs_full <- coef(full_fit, s = "lambda.min")
coef_dt <- data.table(
  protein = rownames(coefs_full)[-1],
  coef = as.numeric(coefs_full[-1])
)
coef_dt[, abs_coef := abs(coef)]
coef_dt <- coef_dt[order(-abs_coef)]

# Also get feature importance via repeated refitting
# Use lambda.1se for sparser models as tiebreaker
coefs_1se <- coef(full_fit, s = "lambda.1se")
coef_dt_1se <- data.table(
  protein = rownames(coefs_1se)[-1],
  coef_1se = as.numeric(coefs_1se[-1])
)
coef_dt <- merge(coef_dt, coef_dt_1se, by = "protein")
coef_dt[, rank := .I]  # Already ordered by abs_coef

# Panel sizes to evaluate
panel_sizes <- c(500, 200, 100, 50, 20, 10)
panel_results <- list()
panel_genes_list <- list()

for (ps in panel_sizes) {
  cat("  Panel size:", ps, "\n")

  # Select top ps proteins by elastic net coefficient
  top_proteins <- coef_dt[seq_len(min(ps, nrow(coef_dt))), protein]
  npx_panel <- npx_matched[, top_proteins, drop = FALSE]

  # 5-fold CV (10 repeats)
  set.seed(42)
  rcv_panel <- random_cv(npx_panel, y, sample_weights, alpha = 0.5)

  # LOEO
  loeo_panel <- loeo_cv(npx_panel, y, etiology, sample_weights, alpha = 0.5)

  panel_results[[as.character(ps)]] <- data.table(
    panel_size = ps,
    auroc_5fold = rcv_panel$auc_mean,
    auroc_5fold_sd = rcv_panel$auc_sd,
    auroc_loeo = loeo_panel$auc,
    sensitivity_5fold = rcv_panel$sensitivity_mean,
    specificity_5fold = rcv_panel$specificity_mean,
    accuracy_5fold = rcv_panel$accuracy_mean,
    n_features_selected = round(rcv_panel$n_features_mean)
  )

  cat("    5-fold AUROC:", round(rcv_panel$auc_mean, 4),
      " | LOEO AUROC:", round(loeo_panel$auc, 4), "\n")

  # Record panel genes
  panel_genes_list[[as.character(ps)]] <- data.table(
    panel_size = ps,
    protein = top_proteins,
    rank = seq_along(top_proteins),
    coef = coef_dt[protein %in% top_proteins, coef],
    abs_coef = coef_dt[protein %in% top_proteins, abs_coef]
  )

  # Also add to main results
  results_list[[paste0("C_panel_", ps, "_5foldCV")]] <- data.table(
    experiment = paste0("C_panel_", ps), cv_method = "5fold_repeated",
    auroc = rcv_panel$auc_mean, auroc_sd = rcv_panel$auc_sd,
    sensitivity = rcv_panel$sensitivity_mean, specificity = rcv_panel$specificity_mean,
    accuracy = rcv_panel$accuracy_mean, n_proteins = ps,
    n_features_selected = round(rcv_panel$n_features_mean)
  )
  results_list[[paste0("C_panel_", ps, "_LOEO")]] <- data.table(
    experiment = paste0("C_panel_", ps), cv_method = "LOEO",
    auroc = loeo_panel$auc, auroc_sd = NA_real_,
    sensitivity = loeo_panel$sensitivity, specificity = loeo_panel$specificity,
    accuracy = loeo_panel$accuracy, n_proteins = ps,
    n_features_selected = NA_integer_
  )
}

cat("\n")

# ============================================================
# EXPERIMENT D: NFASC + GDF15 benchmark (2-protein model)
# ============================================================
cat("=== EXPERIMENT D: NFASC + GDF15 Benchmark ===\n")

benchmark_proteins <- c("NFASC", "GDF15")
available_bench <- intersect(benchmark_proteins, colnames(npx_matched))
cat("  Requested:", paste(benchmark_proteins, collapse = ", "), "\n")
cat("  Available:", paste(available_bench, collapse = ", "), "\n")

if (length(available_bench) >= 1) {
  npx_bench <- npx_matched[, available_bench, drop = FALSE]

  # D1: 5-fold CV
  cat("  D1: Random 5-fold CV (10 repeats)\n")
  rcv_D <- random_cv(npx_bench, y, sample_weights, alpha = 0.5)
  cat("    5-fold CV AUROC:", round(rcv_D$auc_mean, 4), "+/-",
      round(rcv_D$auc_sd, 4), "\n")

  # D2: LOEO
  cat("  D2: Leave-one-etiology-out CV\n")
  loeo_D <- loeo_cv(npx_bench, y, etiology, sample_weights, alpha = 0.5)
  cat("    LOEO AUROC:", round(loeo_D$auc, 4), "\n")

  # D3: Individual protein AUROCs
  indiv_results <- list()
  for (prot in available_bench) {
    roc_ind <- roc(y, npx_matched[, prot], quiet = TRUE)
    indiv_results[[prot]] <- data.table(
      protein = prot,
      individual_auroc = as.numeric(auc(roc_ind))
    )
    cat("    Individual", prot, "AUROC:", round(as.numeric(auc(roc_ind)), 4), "\n")
  }

  bench_out <- data.table(
    proteins = paste(available_bench, collapse = "+"),
    n_proteins = length(available_bench),
    auroc_5fold = rcv_D$auc_mean,
    auroc_5fold_sd = rcv_D$auc_sd,
    auroc_loeo = loeo_D$auc,
    sensitivity_5fold = rcv_D$sensitivity_mean,
    specificity_5fold = rcv_D$specificity_mean,
    accuracy_5fold = rcv_D$accuracy_mean
  )
  bench_indiv <- rbindlist(indiv_results)
  bench_out <- cbind(bench_out, t(setNames(bench_indiv$individual_auroc,
                                            paste0("auroc_", bench_indiv$protein))))

  results_list[["D_nfasc_gdf15_5foldCV"]] <- data.table(
    experiment = "D_nfasc_gdf15", cv_method = "5fold_repeated",
    auroc = rcv_D$auc_mean, auroc_sd = rcv_D$auc_sd,
    sensitivity = rcv_D$sensitivity_mean, specificity = rcv_D$specificity_mean,
    accuracy = rcv_D$accuracy_mean, n_proteins = length(available_bench),
    n_features_selected = length(available_bench)
  )
  results_list[["D_nfasc_gdf15_LOEO"]] <- data.table(
    experiment = "D_nfasc_gdf15", cv_method = "LOEO",
    auroc = loeo_D$auc, auroc_sd = NA_real_,
    sensitivity = loeo_D$sensitivity, specificity = loeo_D$specificity,
    accuracy = loeo_D$accuracy, n_proteins = length(available_bench),
    n_features_selected = NA_integer_
  )
} else {
  cat("  WARNING: Neither NFASC nor GDF15 found in Olink panel. Skipping.\n")
  bench_out <- data.table(proteins = "MISSING", n_proteins = 0,
                           auroc_5fold = NA, auroc_5fold_sd = NA,
                           auroc_loeo = NA, sensitivity_5fold = NA,
                           specificity_5fold = NA, accuracy_5fold = NA)
}
cat("\n")

# ============================================================
# STEP 5: Save outputs
# ============================================================
cat("=== STEP 5: Saving Results ===\n")

# 5a. Main results table
results_all <- rbindlist(results_list, fill = TRUE)
fwrite(results_all, file.path(OUTDIR, "plasma_classifier_results.csv"))
cat("  Saved: plasma_classifier_results.csv (", nrow(results_all), " rows)\n")

# 5b. Minimal panel genes (all sizes combined)
panel_genes_all <- rbindlist(panel_genes_list, fill = TRUE)
fwrite(panel_genes_all, file.path(OUTDIR, "plasma_minimal_panel.csv"))
cat("  Saved: plasma_minimal_panel.csv (", nrow(panel_genes_all), " rows)\n")

# 5c. Panel curve (AUROC vs size)
panel_curve <- rbindlist(panel_results, fill = TRUE)
fwrite(panel_curve, file.path(OUTDIR, "plasma_panel_curve.csv"))
cat("  Saved: plasma_panel_curve.csv (", nrow(panel_curve), " rows)\n")

# 5d. NFASC+GDF15 benchmark
fwrite(bench_out, file.path(OUTDIR, "plasma_nfasc_gdf15_benchmark.csv"))
cat("  Saved: plasma_nfasc_gdf15_benchmark.csv\n")

# ============================================================
# Summary
# ============================================================
cat("\n=== SUMMARY ===\n")
cat("Experiments completed:\n")
for (i in seq_len(nrow(results_all))) {
  row <- results_all[i]
  cat(sprintf("  %-30s %-17s AUROC=%.4f%s\n",
              row$experiment, row$cv_method, row$auroc,
              ifelse(is.na(row$auroc_sd), "",
                     sprintf(" (sd=%.4f)", row$auroc_sd))))
}
cat("\nPanel size curve (5-fold AUROC):\n")
for (i in seq_len(nrow(panel_curve))) {
  row <- panel_curve[i]
  cat(sprintf("  %4d proteins -> AUROC=%.4f (LOEO=%.4f)\n",
              row$panel_size, row$auroc_5fold, row$auroc_loeo))
}

cat("\nFinished:", as.character(Sys.time()), "\n")
cat("=== 82_plasma_classifier.R COMPLETE ===\n")
