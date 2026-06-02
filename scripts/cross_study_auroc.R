#!/usr/bin/env Rscript
# cross_study_auroc.R
# ---------------------------------------------------------------------------
# Cross-study disease prediction using AUROC.
# A. Single-cohort train → held-out test (8x7 = 56 AUROCs)
# B. Integrated LOO train → held-out test (8 AUROCs)
#
# Algorithm: For each train signature, select top 500 genes by |t-stat|,
# compute dot-product activation score in test cohort, report AUROC.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(pROC)
})

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
RDIR <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
LOO_DIR <- file.path(RDIR, "loo_cv")
PS_DIR <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/per_study")

COHORTS <- c("GSE126848", "GSE130970", "GSE135251", "GSE162694",
             "GSE174478", "GSE193066", "GSE213621", "GSE240729")
# NOTE: GSE174478, GSE193066, GSE240729 have no healthy controls — their per-study
# DE contrast is advanced vs early fibrosis (not Disease vs Control). They are used
# as training signatures (fibrosis-associated genes overlap with disease genes) but
# CANNOT be used as test cohorts (AUROC requires both classes). The AUROC matrix is
# therefore 8 train x 5 test, not 8x8. Figure script annotates these accordingly.
FIBROSIS_ONLY <- c("GSE174478", "GSE193066", "GSE240729")
TOP_N <- 500  # Number of signature genes

cat("=== Cross-Study AUROC Prediction ===\n")
cat(sprintf("Signature size: top %d genes by |t-stat|\n\n", TOP_N))

# --- Load expression data ---
cat("Loading merged DGE...\n")
dge <- readRDS(file.path(RDIR, "merged_dge.rds"))

# Exclude permanently excluded datasets
excl_permanent <- c("GSE167523", "PRJNA512027")
keep <- !dge$samples$dataset %in% excl_permanent
dge <- dge[, keep]

# Compute log-CPM
cat("Computing log-CPM...\n")
logcpm <- cpm(dge, log = TRUE, prior.count = 1)
sample_meta <- data.table(
  sample_id = colnames(dge),
  dataset = dge$samples$dataset,
  group_binary = dge$samples$group_binary
)
cat(sprintf("Expression matrix: %d genes x %d samples\n", nrow(logcpm), ncol(logcpm)))

# --- Load per-study DE results ---
cat("\nLoading per-study DE results...\n")
ps_de <- list()
for (study in COHORTS) {
  f <- file.path(PS_DIR, paste0(study, "_de_results.csv"))
  if (file.exists(f)) {
    dt <- fread(f)
    if ("adj.P.Val" %in% names(dt)) setnames(dt, "adj.P.Val", "padj")
    ps_de[[study]] <- dt
    cat(sprintf("  %s: %d genes\n", study, nrow(dt)))
  } else {
    cat(sprintf("  WARNING: Missing %s\n", f))
  }
}

# --- Helper: compute activation score ---
compute_auroc <- function(signature_dt, test_samples, test_labels, expr_matrix, top_n = TOP_N) {
  # Select top genes by |t-stat|
  sig <- signature_dt[order(-abs(t))][1:min(top_n, nrow(signature_dt))]
  sig_genes <- sig$gene

  # Intersect with available genes in expression matrix
  avail <- intersect(sig_genes, rownames(expr_matrix))
  if (length(avail) < 10) return(NA_real_)

  sig <- sig[gene %in% avail]
  weights <- setNames(sig$logFC, sig$gene)  # Signed LFC weights

  # Extract test expression and z-score per gene within test cohort
  test_expr <- expr_matrix[avail, test_samples, drop = FALSE]
  test_z <- t(scale(t(test_expr)))  # Z-score across samples per gene
  test_z[is.na(test_z)] <- 0  # Handle zero-variance genes

  # Activation score: dot product of z-scored expression with signed weights
  w <- weights[avail]
  scores <- as.numeric(t(test_z) %*% w)

  # AUROC
  tryCatch({
    roc_obj <- roc(response = as.numeric(test_labels == "Disease"),
                   predictor = scores, quiet = TRUE)
    as.numeric(auc(roc_obj))
  }, error = function(e) NA_real_)
}

# ===========================================================================
# A. Single-cohort train → held-out test
# ===========================================================================
cat("\n── A. Single-cohort AUROC matrix ──\n")

auroc_matrix <- data.table(
  train = character(0), test = character(0), auroc = numeric(0),
  n_test = integer(0), n_disease = integer(0), n_control = integer(0),
  train_is_fibrosis_contrast = logical(0)
)

for (train_study in COHORTS) {
  if (is.null(ps_de[[train_study]])) next

  for (test_study in COHORTS) {
    # Get test samples
    test_idx <- sample_meta[dataset == test_study, sample_id]
    test_labels <- sample_meta[dataset == test_study, group_binary]

    if (length(unique(test_labels)) < 2) {
      cat(sprintf("  %s → %s: SKIPPED (only one class)\n", train_study, test_study))
      next
    }

    auc_val <- compute_auroc(ps_de[[train_study]], test_idx, test_labels, logcpm)

    auroc_matrix <- rbind(auroc_matrix, data.table(
      train = train_study, test = test_study, auroc = round(auc_val, 4),
      n_test = length(test_idx),
      n_disease = sum(test_labels == "Disease"),
      n_control = sum(test_labels == "Control"),
      train_is_fibrosis_contrast = train_study %in% FIBROSIS_ONLY
    ))

    cat(sprintf("  %s → %s: AUROC = %.3f (N=%d, D=%d, C=%d)\n",
                train_study, test_study, auc_val, length(test_idx),
                sum(test_labels == "Disease"), sum(test_labels == "Control")))
  }
}

# Save
fwrite(auroc_matrix, file.path(LOO_DIR, "auroc_matrix.csv"))
cat("\nSaved: auroc_matrix.csv\n")

# Summary stats
cross_study <- auroc_matrix[train != test]
self_study <- auroc_matrix[train == test]
cat(sprintf("\nSelf-prediction mean AUROC: %.3f (range: %.3f - %.3f)\n",
            mean(self_study$auroc, na.rm = TRUE),
            min(self_study$auroc, na.rm = TRUE),
            max(self_study$auroc, na.rm = TRUE)))
cat(sprintf("Cross-study mean AUROC: %.3f (range: %.3f - %.3f)\n",
            mean(cross_study$auroc, na.rm = TRUE),
            min(cross_study$auroc, na.rm = TRUE),
            max(cross_study$auroc, na.rm = TRUE)))

# Separate disease-contrast vs fibrosis-contrast training
disease_train <- cross_study[train_is_fibrosis_contrast == FALSE]
fibrosis_train <- cross_study[train_is_fibrosis_contrast == TRUE]
if (nrow(disease_train) > 0) {
  cat(sprintf("  Disease-contrast train mean: %.3f\n", mean(disease_train$auroc, na.rm = TRUE)))
}
if (nrow(fibrosis_train) > 0) {
  cat(sprintf("  Fibrosis-contrast train mean: %.3f\n", mean(fibrosis_train$auroc, na.rm = TRUE)))
}

# ===========================================================================
# B. Integrated LOO train → held-out test
# ===========================================================================
cat("\n── B. Integrated LOO AUROC ──\n")

auroc_integrated <- data.table(
  held_out = character(0), auroc = numeric(0),
  n_test = integer(0), n_disease = integer(0), n_control = integer(0)
)

for (held_out in COHORTS) {
  loo_file <- file.path(LOO_DIR, paste0("dream_loo_", held_out, ".csv"))
  if (!file.exists(loo_file)) {
    cat(sprintf("  Holding out %s: SKIPPED (LOO file missing)\n", held_out))
    next
  }

  loo_dream <- fread(loo_file)

  # Test on the held-out cohort
  test_idx <- sample_meta[dataset == held_out, sample_id]
  test_labels <- sample_meta[dataset == held_out, group_binary]

  if (length(unique(test_labels)) < 2) {
    cat(sprintf("  Holding out %s: SKIPPED (only one class in held-out)\n", held_out))
    next
  }

  auc_val <- compute_auroc(loo_dream, test_idx, test_labels, logcpm)

  auroc_integrated <- rbind(auroc_integrated, data.table(
    held_out = held_out, auroc = round(auc_val, 4),
    n_test = length(test_idx),
    n_disease = sum(test_labels == "Disease"),
    n_control = sum(test_labels == "Control")
  ))

  cat(sprintf("  Integrated (excl %s) → %s: AUROC = %.3f (N=%d)\n",
              held_out, held_out, auc_val, length(test_idx)))
}

# Save
fwrite(auroc_integrated, file.path(LOO_DIR, "auroc_integrated.csv"))
cat("\nSaved: auroc_integrated.csv\n")

# Summary
if (nrow(auroc_integrated) > 0) {
  cat(sprintf("\nIntegrated LOO mean AUROC: %.3f (range: %.3f - %.3f)\n",
              mean(auroc_integrated$auroc, na.rm = TRUE),
              min(auroc_integrated$auroc, na.rm = TRUE),
              max(auroc_integrated$auroc, na.rm = TRUE)))

  # Compare: per held-out cohort, integrated vs best single-study
  cat("\n── Per held-out cohort: Integrated vs best single-study ──\n")
  for (ho in auroc_integrated$held_out) {
    int_auc <- auroc_integrated[held_out == ho, auroc]
    single_aucs <- auroc_matrix[test == ho & train != ho, auroc]
    best_single <- if (length(single_aucs) > 0) max(single_aucs, na.rm = TRUE) else NA
    mean_single <- if (length(single_aucs) > 0) mean(single_aucs, na.rm = TRUE) else NA
    cat(sprintf("  %s: Integrated=%.3f, Best_single=%.3f, Mean_single=%.3f\n",
                ho, int_auc, best_single, mean_single))
  }
}

cat("\nDone!\n")
