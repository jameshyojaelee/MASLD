#!/usr/bin/env Rscript
# =============================================================================
# 148_deconv_gene_panel.R
#
# Cell-Type Composition Gene Panel: elastic net models predicting BayesPrism
# cell-type proportions from a minimal gene panel.
#
# Overview:
#   1. Load bulk DGE (34K genes x 1,444 samples) + BayesPrism proportions
#   2. Build candidate gene pool from scRNA markers + transition programs
#   3. LOCO-CV elastic net for 4 cell types (Hepatocyte, Macrophage,
#      Endothelial, Stellate) with logit-transformed proportions
#   4. Extract minimal panels (10, 25, 50 genes) ranked by cross-fold stability
#   5. Re-evaluate panels, stellate binary classification, benchmarks
#
# Follows Script 77 LOCO-CV + elastic net pattern.
#
# SLURM:
#   #SBATCH --partition=io
#   #SBATCH --cpus-per-task=8
#   #SBATCH --mem=32G
#   #SBATCH --time=48:00:00
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(glmnet)
  library(pROC)
  library(Matrix)
})

set.seed(42)

# --- Paths -------------------------------------------------------------------
project_root <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

int_root <- file.path(project_root,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration")

dge_path      <- file.path(int_root, "results/integration/merged_dge.rds")
meta_path     <- file.path(int_root, "results/staging_classifier/modeling_metadata.csv")
bp_path       <- file.path(int_root,
  "results/progression/cibersortx_celltype_expression/bayesprism_proportions.csv")
scrna_labels  <- file.path(int_root,
  "results/progression/scrna_reference/scrna_cell_labels.csv")
scrna_counts  <- file.path(int_root,
  "results/progression/scrna_reference/scrna_counts.csv.gz")
transition_path <- file.path(int_root,
  "results/progression/celltype_transition_programs.csv")
gene_annot    <- file.path(int_root,
  "results/gene_annotation/human_ensg_to_symbol.tsv")
panel50_path  <- file.path(int_root,
  "results/staging_classifier/minimal_panel_50.csv")
plasma_path   <- file.path(int_root,
  "results/staging_classifier/tissue_plasma_bridge.csv")

out_dir <- file.path(int_root, "results/deconvolution_panel")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Target cell types
target_types <- c("Hepatocyte", "Macrophage", "Endothelial", "Stellate")
stellate_threshold <- 0.04  # 4% threshold for binary stellate classification

# =============================================================================
# 1. LOAD DATA
# =============================================================================
cat("=== 1. Loading data ===\n")

# --- 1a. DGE -> TMM logCPM ---
cat("  Loading DGE ...\n")
dge <- readRDS(dge_path)
dge <- calcNormFactors(dge, method = "TMM")
logcpm <- cpm(dge, log = TRUE, prior.count = 1)
cat("  logCPM:", nrow(logcpm), "genes x", ncol(logcpm), "samples\n")

# --- 1b. Gene annotation (Ensembl -> symbol) ---
cat("  Loading gene annotation ...\n")
annot <- fread(gene_annot, sep = "\t")
# Create Ensembl-versioned -> symbol lookup
ensg_to_sym <- setNames(annot$symbol, annot$gene_id)
# Also create base Ensembl -> symbol (for matching scRNA symbols)
sym_to_ensg <- setNames(annot$gene_id, annot$symbol)
# Handle duplicates: keep first occurrence
sym_to_ensg <- sym_to_ensg[!duplicated(names(sym_to_ensg))]
cat("  Gene annotation:", length(ensg_to_sym), "mappings\n")

# --- 1c. Metadata with LOCO folds ---
cat("  Loading metadata ...\n")
meta <- fread(meta_path)
cat("  Metadata:", nrow(meta), "samples\n")

# --- 1d. BayesPrism proportions ---
cat("  Loading BayesPrism proportions ...\n")
bp <- fread(bp_path)
cat("  BayesPrism:", nrow(bp), "samples x", ncol(bp) - 1, "cell types\n")

# --- 1e. scRNA reference ---
cat("  Loading scRNA cell labels ...\n")
sc_labels <- fread(scrna_labels)
cat("  scRNA labels:", nrow(sc_labels), "cells,",
    length(unique(sc_labels$cell_type)), "cell types\n")

cat("  Loading scRNA counts (this may take a moment) ...\n")
sc_counts <- fread(scrna_counts)
# First column is cell_id (unnamed column V1)
sc_cell_ids <- sc_counts[[1]]
sc_gene_names <- colnames(sc_counts)[-1]
sc_mat <- as.matrix(sc_counts[, -1, with = FALSE])
rownames(sc_mat) <- sc_cell_ids
rm(sc_counts); gc(verbose = FALSE)
cat("  scRNA counts:", nrow(sc_mat), "cells x", ncol(sc_mat), "genes\n")

# --- 1f. Transition programs ---
cat("  Loading transition programs ...\n")
trans <- fread(transition_path)
cat("  Transition programs:", nrow(trans), "rows\n")

# =============================================================================
# 2. ALIGN DATA & BUILD GENE ID MAPPINGS
# =============================================================================
cat("\n=== 2. Aligning data ===\n")

# DGE rownames are versioned Ensembl IDs; scRNA has gene symbols
# Build bidirectional mapping for candidate genes
bulk_genes <- rownames(logcpm)  # ENSG00000310526.1 etc.
bulk_symbols <- ensg_to_sym[bulk_genes]

# Map bulk Ensembl -> symbol (keep valid only)
ensg_sym_map <- data.table(
  ensg = bulk_genes,
  symbol = as.character(bulk_symbols)
)
ensg_sym_map <- ensg_sym_map[!is.na(symbol) & symbol != ""]
# Handle symbol duplicates: keep first Ensembl ID per symbol
ensg_sym_map <- ensg_sym_map[!duplicated(symbol)]
cat("  Mapped", nrow(ensg_sym_map), "Ensembl IDs to unique symbols\n")

# Intersection of scRNA genes and bulk genes (by symbol)
shared_genes <- intersect(sc_gene_names, ensg_sym_map$symbol)
cat("  Shared genes (scRNA-bulk):", length(shared_genes), "\n")

# Align samples across bulk, BayesPrism, and metadata
common_ids <- Reduce(intersect, list(
  colnames(logcpm),
  bp$sample_id,
  meta$sample_id
))
cat("  Common samples:", length(common_ids), "\n")

# Subset and align
meta <- meta[match(common_ids, sample_id)]
bp   <- bp[match(common_ids, sample_id)]
logcpm <- logcpm[, common_ids, drop = FALSE]

# Extract proportion matrix for target types
prop_mat <- as.matrix(bp[, ..target_types])
rownames(prop_mat) <- common_ids
cat("  Proportion matrix:", nrow(prop_mat), "samples x", ncol(prop_mat),
    "cell types\n")
cat("  Proportion ranges:\n")
for (ct in target_types) {
  cat(sprintf("    %s: [%.4f, %.4f] mean=%.4f\n",
              ct, min(prop_mat[, ct]), max(prop_mat[, ct]), mean(prop_mat[, ct])))
}

# LOCO folds for fibrosis (6 cohort-based folds)
folds <- meta$loco_fold_fibrosis
valid_idx <- which(folds != "excluded")
cat("  Samples with valid LOCO folds:", length(valid_idx), "of", nrow(meta), "\n")

# =============================================================================
# 3. BUILD CANDIDATE GENE POOL
# =============================================================================
cat("\n=== 3. Building candidate gene pool ===\n")

# --- Pool A: scRNA cell-type markers (one-vs-rest Wilcoxon, top 200 per type)
cat("  Computing scRNA markers (Wilcoxon one-vs-rest) ...\n")

# Align scRNA labels to counts
sc_label_map <- setNames(sc_labels$cell_type, sc_labels$cell_id)
sc_labels_aligned <- sc_label_map[rownames(sc_mat)]

# Only compute for shared genes (present in both bulk and scRNA)
sc_shared_idx <- which(colnames(sc_mat) %in% shared_genes)
sc_mat_shared <- sc_mat[, sc_shared_idx, drop = FALSE]

marker_genes_list <- list()
n_markers_per_type <- 200

for (ct in target_types) {
  cat(sprintf("    %s: ", ct))
  is_type <- (sc_labels_aligned == ct)
  if (sum(is_type, na.rm = TRUE) < 10) {
    cat("SKIP (< 10 cells)\n")
    next
  }

  # Wilcoxon rank-sum test on each gene
  pvals <- numeric(ncol(sc_mat_shared))
  lfcs  <- numeric(ncol(sc_mat_shared))
  names(pvals) <- colnames(sc_mat_shared)
  names(lfcs)  <- colnames(sc_mat_shared)

  for (j in seq_len(ncol(sc_mat_shared))) {
    x_in  <- sc_mat_shared[which(is_type), j]
    x_out <- sc_mat_shared[which(!is_type), j]

    # Skip genes with no variation
    if (sd(c(x_in, x_out)) == 0) {
      pvals[j] <- 1
      lfcs[j]  <- 0
      next
    }

    wt <- suppressWarnings(wilcox.test(x_in, x_out, alternative = "greater"))
    pvals[j] <- wt$p.value
    lfcs[j]  <- log2(mean(x_in) + 1) - log2(mean(x_out) + 1)
  }

  # Take top 200 by p-value (upregulated in target type)
  top_idx <- order(pvals)[1:min(n_markers_per_type, length(pvals))]
  top_markers <- names(pvals)[top_idx]
  marker_genes_list[[ct]] <- top_markers
  cat(sprintf("%d markers (min p=%.2e)\n", length(top_markers), min(pvals)))
}

pool_a <- unique(unlist(marker_genes_list))
cat("  Pool A (scRNA markers):", length(pool_a), "unique genes\n")

# --- Pool B: Transition-specific genes (Stellate DEGs)
cat("  Extracting Stellate transition genes ...\n")
stellate_trans <- trans[cell_type == "Stellate" & padj < 0.01 & abs(lfc) > 0.3]
pool_b <- unique(stellate_trans$gene)
# pool_b genes are symbols; intersect with shared genes
pool_b <- intersect(pool_b, shared_genes)
cat("  Pool B (Stellate transition DEGs):", length(pool_b), "genes\n")

# --- Union ---
candidate_symbols <- unique(c(pool_a, pool_b))
# Map to Ensembl IDs (for subsetting bulk data)
cand_map <- ensg_sym_map[symbol %in% candidate_symbols]
candidate_ensg <- cand_map$ensg
cat("  Total candidate pool:", length(candidate_ensg), "genes mapped to Ensembl\n")

# Verify candidates exist in bulk
candidate_ensg <- intersect(candidate_ensg, rownames(logcpm))
cat("  Candidates in bulk data:", length(candidate_ensg), "\n")

# Create symbol lookup for candidates
cand_ensg_to_sym <- setNames(cand_map$symbol, cand_map$ensg)
cand_ensg_to_sym <- cand_ensg_to_sym[names(cand_ensg_to_sym) %in% candidate_ensg]

# Free scRNA memory
rm(sc_mat, sc_mat_shared, sc_labels_aligned, sc_label_map)
gc(verbose = FALSE)

# =============================================================================
# 4. LOCO-CV FOR COMPOSITION PREDICTION
# =============================================================================
cat("\n=== 4. LOCO-CV composition prediction ===\n")

fold_names <- sort(unique(folds[valid_idx]))
n_folds <- length(fold_names)
cat("  LOCO folds (", n_folds, "):", paste(fold_names, collapse = ", "), "\n")

# Use all samples (not just fibrosis-valid) since deconvolution targets
# are available for all 1,444 samples. We use the fibrosis LOCO folds
# for cross-cohort validation since they partition by cohort.
# However, we need valid fold assignments; use all non-excluded samples.

n_var_filter <- min(500, length(candidate_ensg))  # within-fold variance filter

# Storage
predictions_list <- list()
coef_list <- list()

logit_transform <- function(p, eps = 1e-6) {
  p <- pmin(pmax(p, eps), 1 - eps)
  log(p / (1 - p))
}

inv_logit <- function(x) {
  1 / (1 + exp(-x))
}

for (fi in seq_along(fold_names)) {
  fold <- fold_names[fi]
  cat(sprintf("\n--- Fold %d/%d: %s (held out) ---\n", fi, n_folds, fold))

  test_idx  <- which(folds == fold)
  train_idx <- which(folds != fold & folds != "excluded")
  cat(sprintf("  Train: %d  Test: %d\n", length(train_idx), length(test_idx)))

  if (length(test_idx) < 10) {
    cat("  SKIP: too few test samples\n")
    next
  }

  # --- Within-fold variance filter on candidate genes ---
  train_expr <- logcpm[candidate_ensg, train_idx, drop = FALSE]
  gene_vars <- apply(train_expr, 1, var)
  top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:n_var_filter]

  # --- Rank-transform expression ---
  expr_train <- t(logcpm[top_genes, train_idx, drop = FALSE])
  expr_train <- apply(expr_train, 2, rank) / nrow(expr_train)

  expr_test <- t(logcpm[top_genes, test_idx, drop = FALSE])
  expr_test <- apply(expr_test, 2, rank) / nrow(expr_test)

  # --- Fit elastic net for each cell type ---
  for (ct in target_types) {
    cat(sprintf("  %s: ", ct))

    y_train <- logit_transform(prop_mat[train_idx, ct])
    y_test_raw <- prop_mat[test_idx, ct]

    # Fit cv.glmnet (alpha=0.5 elastic net)
    cvfit <- tryCatch(
      cv.glmnet(
        x = expr_train,
        y = y_train,
        alpha = 0.5,
        nfolds = 5,
        type.measure = "mse"
      ),
      error = function(e) {
        cat("ERROR:", conditionMessage(e), "\n")
        return(NULL)
      }
    )

    if (is.null(cvfit)) next

    # Predict and inverse-logit
    pred_logit <- as.numeric(predict(cvfit, newx = expr_test, s = "lambda.1se"))
    pred_prop  <- inv_logit(pred_logit)

    # Spearman correlation
    sp <- cor(pred_prop, y_test_raw, method = "spearman")
    rmse <- sqrt(mean((pred_prop - y_test_raw)^2))
    cat(sprintf("rho=%.3f RMSE=%.4f ", sp, rmse))

    # Record predictions
    predictions_list[[length(predictions_list) + 1]] <- data.table(
      sample_id = common_ids[test_idx],
      cell_type = ct,
      fold = fold,
      predicted = pred_prop,
      actual = y_test_raw
    )

    # Record non-zero coefficients
    coefs <- as.matrix(coef(cvfit, s = "lambda.1se"))
    nonzero <- coefs[coefs[, 1] != 0, , drop = FALSE]
    nonzero <- nonzero[rownames(nonzero) != "(Intercept)", , drop = FALSE]

    if (nrow(nonzero) > 0) {
      # Map Ensembl back to symbols
      gene_syms <- cand_ensg_to_sym[rownames(nonzero)]
      gene_syms[is.na(gene_syms)] <- rownames(nonzero)[is.na(gene_syms)]

      coef_list[[length(coef_list) + 1]] <- data.table(
        gene_ensg = rownames(nonzero),
        gene_symbol = as.character(gene_syms),
        cell_type = ct,
        fold = fold,
        coefficient = nonzero[, 1]
      )
      cat(sprintf("(%d genes selected)\n", nrow(nonzero)))
    } else {
      cat("(0 genes selected)\n")
    }
  }
}

# Assemble predictions
predictions <- rbindlist(predictions_list)
cat("\n  Total predictions:", nrow(predictions), "\n")

# Assemble coefficients
all_coefs <- rbindlist(coef_list)
cat("  Total coefficient entries:", nrow(all_coefs), "\n")

# =============================================================================
# 5. PERFORMANCE SUMMARY (FULL CANDIDATE POOL)
# =============================================================================
cat("\n=== 5. Performance summary (full candidate pool) ===\n")

perf_full <- predictions[, .(
  spearman = cor(predicted, actual, method = "spearman"),
  rmse = sqrt(mean((predicted - actual)^2)),
  mae = mean(abs(predicted - actual)),
  n_samples = .N
), by = .(cell_type)]

cat("  Per-cell-type performance (full candidate pool):\n")
print(perf_full)

# =============================================================================
# 6. EXTRACT MINIMAL PANELS (10, 25, 50 genes)
# =============================================================================
cat("\n=== 6. Extracting minimal panels ===\n")

# Tally gene selection frequency across folds x cell types
gene_stability <- all_coefs[, .(
  n_selections = .N,
  n_folds = uniqueN(fold),
  n_celltypes = uniqueN(cell_type),
  mean_abs_coef = mean(abs(coefficient)),
  cell_types = paste(sort(unique(cell_type)), collapse = ";")
), by = .(gene_ensg, gene_symbol)]

# Rank by stability: primary = n_selections, secondary = mean_abs_coef
gene_stability <- gene_stability[order(-n_selections, -mean_abs_coef)]
gene_stability[, stability_rank := .I]

cat("  Total unique genes selected across folds:", nrow(gene_stability), "\n")
cat("  Top 10 most stable genes:\n")
print(head(gene_stability[, .(gene_symbol, n_selections, n_folds,
                               n_celltypes, mean_abs_coef)], 10))

# --- Re-evaluate panels at sizes 10, 25, 50 ---
panel_sizes <- c(10, 25, 50)
panel_results_list <- list()
panel_genes_list <- list()

for (ps in panel_sizes) {
  cat(sprintf("\n--- Panel size: %d ---\n", ps))

  if (ps > nrow(gene_stability)) {
    cat("  SKIP: only", nrow(gene_stability), "genes available\n")
    next
  }

  panel_ensg <- gene_stability$gene_ensg[1:ps]
  panel_sym  <- gene_stability$gene_symbol[1:ps]

  # Save panel
  panel_dt <- data.table(
    gene_ensg = panel_ensg,
    gene_symbol = panel_sym,
    stability_rank = 1:ps,
    n_selections = gene_stability$n_selections[1:ps],
    cell_types = gene_stability$cell_types[1:ps]
  )
  panel_genes_list[[as.character(ps)]] <- panel_dt

  # Re-fit LOCO-CV using only panel genes
  panel_pred_list <- list()

  for (fi in seq_along(fold_names)) {
    fold <- fold_names[fi]
    test_idx  <- which(folds == fold)
    train_idx <- which(folds != fold & folds != "excluded")

    if (length(test_idx) < 10) next

    # Rank-transform panel genes
    expr_train <- t(logcpm[panel_ensg, train_idx, drop = FALSE])
    expr_train <- apply(expr_train, 2, rank) / nrow(expr_train)

    expr_test <- t(logcpm[panel_ensg, test_idx, drop = FALSE])
    expr_test <- apply(expr_test, 2, rank) / nrow(expr_test)

    for (ct in target_types) {
      y_train <- logit_transform(prop_mat[train_idx, ct])
      y_test_raw <- prop_mat[test_idx, ct]

      cvfit <- tryCatch(
        cv.glmnet(
          x = expr_train,
          y = y_train,
          alpha = 0.5,
          nfolds = 5,
          type.measure = "mse"
        ),
        error = function(e) NULL
      )

      if (is.null(cvfit)) next

      pred_logit <- as.numeric(predict(cvfit, newx = expr_test, s = "lambda.1se"))
      pred_prop  <- inv_logit(pred_logit)

      panel_pred_list[[length(panel_pred_list) + 1]] <- data.table(
        sample_id = common_ids[test_idx],
        cell_type = ct,
        fold = fold,
        predicted = pred_prop,
        actual = y_test_raw,
        panel_size = ps
      )
    }
  }

  panel_preds <- rbindlist(panel_pred_list)

  # Performance per cell type
  panel_perf <- panel_preds[, .(
    spearman = cor(predicted, actual, method = "spearman"),
    rmse = sqrt(mean((predicted - actual)^2)),
    mae = mean(abs(predicted - actual)),
    n_samples = .N
  ), by = .(cell_type)]
  panel_perf[, panel_size := ps]

  panel_results_list[[as.character(ps)]] <- panel_perf
  cat("  Performance:\n")
  print(panel_perf[, .(cell_type, spearman, rmse)])
}

# Combine all panel performances
panel_perf_all <- rbindlist(panel_results_list)

# Add full candidate pool as "reference"
perf_full[, panel_size := length(candidate_ensg)]
panel_perf_all <- rbind(panel_perf_all, perf_full)

cat("\n  Combined panel performance:\n")
print(dcast(panel_perf_all, cell_type ~ panel_size, value.var = "spearman"))

# =============================================================================
# 7. STELLATE BINARY CLASSIFICATION
# =============================================================================
cat("\n=== 7. Stellate binary classification ===\n")

# Use full-pool predictions for stellate
stellate_preds <- predictions[cell_type == "Stellate"]

if (nrow(stellate_preds) > 0) {
  stellate_preds[, actual_high := as.integer(actual > stellate_threshold)]
  stellate_preds[, pred_high := as.integer(predicted > stellate_threshold)]

  cat("  Actual high stellate (>4%):", sum(stellate_preds$actual_high),
      "of", nrow(stellate_preds), "\n")
  cat("  Predicted high stellate:", sum(stellate_preds$pred_high), "\n")

  # AUROC
  if (length(unique(stellate_preds$actual_high)) == 2) {
    roc_stellate <- roc(stellate_preds$actual_high, stellate_preds$predicted,
                        quiet = TRUE)
    auc_stellate <- as.numeric(auc(roc_stellate))
    cat(sprintf("  Stellate >4%% AUROC: %.3f\n", auc_stellate))

    # Sensitivity / specificity at threshold
    sens <- sum(stellate_preds$pred_high == 1 & stellate_preds$actual_high == 1) /
            max(1, sum(stellate_preds$actual_high == 1))
    spec <- sum(stellate_preds$pred_high == 0 & stellate_preds$actual_high == 0) /
            max(1, sum(stellate_preds$actual_high == 0))
    cat(sprintf("  Sensitivity: %.3f  Specificity: %.3f\n", sens, spec))
  } else {
    auc_stellate <- NA
    sens <- NA
    spec <- NA
    cat("  Only one class present; cannot compute AUROC\n")
  }

  # --- Secondary: does predicted stellate >4% predict F>=3? ---
  cat("\n  Secondary: predicted Stellate >4%% as F>=3 predictor\n")
  fib_valid <- meta[fib_ge3 != -1]
  stellate_fib <- merge(stellate_preds, fib_valid[, .(sample_id, fib_ge3)],
                        by = "sample_id")

  if (nrow(stellate_fib) > 0 && length(unique(stellate_fib$fib_ge3)) == 2) {
    roc_fib <- roc(stellate_fib$fib_ge3, stellate_fib$predicted, quiet = TRUE)
    auc_fib <- as.numeric(auc(roc_fib))
    cat(sprintf("  Predicted stellate -> F>=3 AUROC: %.3f\n", auc_fib))
  } else {
    auc_fib <- NA
    cat("  Cannot compute F>=3 AUROC (insufficient data)\n")
  }

  # Save stellate binary results
  stellate_binary <- data.table(
    metric = c("auroc_stellate_binary", "sensitivity", "specificity",
               "n_actual_high", "n_pred_high", "n_total",
               "auroc_stellate_to_fib"),
    value = c(auc_stellate, sens, spec,
              sum(stellate_preds$actual_high),
              sum(stellate_preds$pred_high),
              nrow(stellate_preds),
              auc_fib)
  )
} else {
  stellate_binary <- data.table(metric = "error", value = NA)
  cat("  No stellate predictions available\n")
}

# =============================================================================
# 8. BENCHMARKS
# =============================================================================
cat("\n=== 8. Benchmarks ===\n")

benchmark_list <- list()

# --- 8a. Full BayesPrism (upper bound = perfect correlation with itself) ---
cat("  Benchmark A: Full BayesPrism (upper bound = 1.0)\n")
for (ct in target_types) {
  benchmark_list[[length(benchmark_list) + 1]] <- data.table(
    method = "BayesPrism_full", cell_type = ct,
    spearman = 1.0, rmse = 0.0, note = "upper_bound"
  )
}

# --- 8b. Staging panel (50 genes) ---
cat("  Benchmark B: Staging panel (50 genes) ...\n")
if (file.exists(panel50_path)) {
  staging_panel <- fread(panel50_path)
  staging_ensg <- staging_panel$gene
  # Ensure they exist in bulk
  staging_ensg <- intersect(staging_ensg, rownames(logcpm))
  cat("    Staging panel genes in bulk:", length(staging_ensg), "\n")

  if (length(staging_ensg) >= 5) {
    staging_pred_list <- list()

    for (fi in seq_along(fold_names)) {
      fold <- fold_names[fi]
      test_idx  <- which(folds == fold)
      train_idx <- which(folds != fold & folds != "excluded")
      if (length(test_idx) < 10) next

      expr_train <- t(logcpm[staging_ensg, train_idx, drop = FALSE])
      expr_train <- apply(expr_train, 2, rank) / nrow(expr_train)
      expr_test <- t(logcpm[staging_ensg, test_idx, drop = FALSE])
      expr_test <- apply(expr_test, 2, rank) / nrow(expr_test)

      for (ct in target_types) {
        y_train <- logit_transform(prop_mat[train_idx, ct])
        y_test_raw <- prop_mat[test_idx, ct]

        cvfit <- tryCatch(
          cv.glmnet(x = expr_train, y = y_train, alpha = 0.5,
                    nfolds = 5, type.measure = "mse"),
          error = function(e) NULL
        )
        if (is.null(cvfit)) next

        pred_logit <- as.numeric(predict(cvfit, newx = expr_test, s = "lambda.1se"))
        pred_prop <- inv_logit(pred_logit)

        staging_pred_list[[length(staging_pred_list) + 1]] <- data.table(
          sample_id = common_ids[test_idx], cell_type = ct, fold = fold,
          predicted = pred_prop, actual = y_test_raw
        )
      }
    }

    staging_preds <- rbindlist(staging_pred_list)
    staging_perf <- staging_preds[, .(
      spearman = cor(predicted, actual, method = "spearman"),
      rmse = sqrt(mean((predicted - actual)^2))
    ), by = cell_type]
    staging_perf[, `:=`(method = "staging_panel_50", note = "benchmark")]

    for (i in seq_len(nrow(staging_perf))) {
      benchmark_list[[length(benchmark_list) + 1]] <- staging_perf[i]
    }
    cat("    Staging panel performance:\n")
    print(staging_perf[, .(cell_type, spearman, rmse)])
  }
} else {
  cat("    Staging panel file not found, skipping\n")
}

# --- 8c. Unrestricted elastic net (top 3K genes by variance) ---
cat("  Benchmark C: Unrestricted elastic net (top 3K by variance) ...\n")

unrestricted_pred_list <- list()

for (fi in seq_along(fold_names)) {
  fold <- fold_names[fi]
  test_idx  <- which(folds == fold)
  train_idx <- which(folds != fold & folds != "excluded")
  if (length(test_idx) < 10) next

  # Within-fold top 3K by variance
  train_logcpm <- logcpm[, train_idx, drop = FALSE]
  gene_vars <- apply(train_logcpm, 1, var)
  top3k <- names(sort(gene_vars, decreasing = TRUE))[1:3000]

  expr_train <- t(logcpm[top3k, train_idx, drop = FALSE])
  expr_train <- apply(expr_train, 2, rank) / nrow(expr_train)
  expr_test <- t(logcpm[top3k, test_idx, drop = FALSE])
  expr_test <- apply(expr_test, 2, rank) / nrow(expr_test)

  for (ct in target_types) {
    y_train <- logit_transform(prop_mat[train_idx, ct])
    y_test_raw <- prop_mat[test_idx, ct]

    cvfit <- tryCatch(
      cv.glmnet(x = expr_train, y = y_train, alpha = 0.5,
                nfolds = 5, type.measure = "mse"),
      error = function(e) NULL
    )
    if (is.null(cvfit)) next

    pred_logit <- as.numeric(predict(cvfit, newx = expr_test, s = "lambda.1se"))
    pred_prop <- inv_logit(pred_logit)

    unrestricted_pred_list[[length(unrestricted_pred_list) + 1]] <- data.table(
      sample_id = common_ids[test_idx], cell_type = ct, fold = fold,
      predicted = pred_prop, actual = y_test_raw
    )
  }
}

unrestricted_preds <- rbindlist(unrestricted_pred_list)
unrestricted_perf <- unrestricted_preds[, .(
  spearman = cor(predicted, actual, method = "spearman"),
  rmse = sqrt(mean((predicted - actual)^2))
), by = cell_type]
unrestricted_perf[, `:=`(method = "unrestricted_3K", note = "benchmark")]

for (i in seq_len(nrow(unrestricted_perf))) {
  benchmark_list[[length(benchmark_list) + 1]] <- unrestricted_perf[i]
}
cat("    Unrestricted 3K performance:\n")
print(unrestricted_perf[, .(cell_type, spearman, rmse)])

# Add our panels to benchmark
for (ps in panel_sizes) {
  pp <- panel_perf_all[panel_size == ps]
  for (i in seq_len(nrow(pp))) {
    benchmark_list[[length(benchmark_list) + 1]] <- data.table(
      method = paste0("deconv_panel_", ps),
      cell_type = pp$cell_type[i],
      spearman = pp$spearman[i],
      rmse = pp$rmse[i],
      note = "ours"
    )
  }
}

# Full candidate pool
pp_full <- panel_perf_all[panel_size == length(candidate_ensg)]
for (i in seq_len(nrow(pp_full))) {
  benchmark_list[[length(benchmark_list) + 1]] <- data.table(
    method = "deconv_full_pool",
    cell_type = pp_full$cell_type[i],
    spearman = pp_full$spearman[i],
    rmse = pp_full$rmse[i],
    note = "ours"
  )
}

benchmarks <- rbindlist(benchmark_list, fill = TRUE)

cat("\n  Benchmark comparison:\n")
print(dcast(benchmarks, method ~ cell_type, value.var = "spearman"))

# =============================================================================
# 9. PLASMA TRANSLATABILITY
# =============================================================================
cat("\n=== 9. Plasma translatability ===\n")

if (file.exists(plasma_path)) {
  plasma <- fread(plasma_path)
  cat("  Loaded tissue_plasma_bridge:", nrow(plasma), "genes\n")

  # Check overlap with our panels
  for (ps in names(panel_genes_list)) {
    panel_dt <- panel_genes_list[[ps]]
    panel_sym <- panel_dt$gene_symbol

    # Match by symbol
    in_plasma <- panel_sym %in% plasma[in_plasma == TRUE, human_symbol]
    n_plasma <- sum(in_plasma)
    cat(sprintf("  Panel %s: %d/%d genes detectable in plasma (%.0f%%)\n",
                ps, n_plasma, length(panel_sym),
                100 * n_plasma / length(panel_sym)))
    if (n_plasma > 0) {
      cat("    Plasma-detectable:", paste(panel_sym[in_plasma], collapse = ", "), "\n")
    }
  }
} else {
  cat("  tissue_plasma_bridge.csv not found, skipping\n")
}

# =============================================================================
# 10. SAVE OUTPUTS
# =============================================================================
cat("\n=== 10. Saving outputs ===\n")

# Gene panels
for (ps in names(panel_genes_list)) {
  out_file <- file.path(out_dir, sprintf("deconv_panel_%s.csv", ps))
  fwrite(panel_genes_list[[ps]], out_file)
  cat("  Wrote", basename(out_file), "\n")
}

# Full gene stability ranking
fwrite(gene_stability, file.path(out_dir, "deconv_gene_stability.csv"))
cat("  Wrote deconv_gene_stability.csv\n")

# Per-sample predictions (full candidate pool)
fwrite(predictions, file.path(out_dir, "deconv_panel_predictions.csv"))
cat("  Wrote deconv_panel_predictions.csv\n")

# Performance table
fwrite(panel_perf_all, file.path(out_dir, "deconv_panel_performance.csv"))
cat("  Wrote deconv_panel_performance.csv\n")

# Stellate binary classification
fwrite(stellate_binary, file.path(out_dir, "deconv_panel_stellate_binary.csv"))
cat("  Wrote deconv_panel_stellate_binary.csv\n")

# Benchmark comparison
fwrite(benchmarks, file.path(out_dir, "deconv_panel_benchmark.csv"))
cat("  Wrote deconv_panel_benchmark.csv\n")

# =============================================================================
# SUMMARY
# =============================================================================
cat("\n")
cat("===========================================================================\n")
cat("  SCRIPT 148 COMPLETE: Cell-Type Composition Gene Panel\n")
cat("===========================================================================\n")
cat(sprintf("  Candidate gene pool: %d genes\n", length(candidate_ensg)))
cat(sprintf("  LOCO folds: %d\n", n_folds))
cat("\n  Full-pool Spearman correlations:\n")
for (i in seq_len(nrow(perf_full))) {
  cat(sprintf("    %s: rho=%.3f RMSE=%.4f\n",
              perf_full$cell_type[i], perf_full$spearman[i], perf_full$rmse[i]))
}
cat(sprintf("\n  Stellate >4%% AUROC: %.3f\n",
            stellate_binary[metric == "auroc_stellate_binary", value]))
cat(sprintf("  Outputs: %s\n", out_dir))
cat("===========================================================================\n")
