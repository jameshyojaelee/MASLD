#!/usr/bin/env Rscript
# =============================================================================
# 80_classifier_sensitivity.R
#
# Sensitivity analyses for the V3 staging classifier (Script 76):
#   1. Bootstrap 95% CI on the observed F>=3 AUROC (1,000 resamples)
#   2. Permutation null (200 permutations of the full LOCO-CV pipeline)
#   3. Feature K sweep (K = 100, 200, 500, 1K, 2K, 3K, 5K)
#   4. DeLong test vs published panels (from Script 73 per-fold AUROCs)
#
# Standard reviewer-facing sensitivity battery for Cell Metabolism.
#
# Inputs:
#   - results/integration/merged_dge.rds (34,453 genes x 1,444 samples)
#   - results/staging_classifier/modeling_metadata.csv
#   - results/staging_classifier/v3_proper_cv_results.csv (V3 per-sample preds)
#   - results/staging_classifier/published_panel_benchmark.csv (Script 73)
#
# Outputs (to results/staging_classifier/):
#   - sensitivity_bootstrap_ci.csv
#   - sensitivity_permutation.csv
#   - sensitivity_k_sweep.csv
#   - delong_comparisons.csv
#
# Usage: Rscript 80_classifier_sensitivity.R
# SLURM: bigmem, 16 CPUs, 200GB RAM, 48h
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(glmnet)
  library(pROC)
  library(parallel)
})

set.seed(42)

# --- Paths -------------------------------------------------------------------
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT    <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR   <- file.path(INT, "results")
OUTDIR <- file.path(RDIR, "staging_classifier")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

NCORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "16"))

cat("=== 80: Classifier Sensitivity Analyses ===\n")
cat("Started:", as.character(Sys.time()), "\n")
cat("Cores:", NCORES, "\n\n")

# =============================================================================
# LOAD SHARED DATA
# =============================================================================
cat("Loading merged DGE (all 1,444 samples)...\n")
dge <- readRDS(file.path(RDIR, "integration/merged_dge.rds"))
dge <- calcNormFactors(dge, method = "TMM")
logcpm <- cpm(dge, log = TRUE, prior.count = 1)
cat("  Expression:", nrow(logcpm), "genes x", ncol(logcpm), "samples\n")

cat("Loading modeling metadata...\n")
meta <- fread(file.path(OUTDIR, "modeling_metadata.csv"))
cat("  Metadata:", nrow(meta), "samples\n")

# Align
common <- intersect(colnames(logcpm), meta$sample_id)
meta <- meta[match(common, sample_id)]
logcpm <- logcpm[, common]
cat("  Aligned:", ncol(logcpm), "samples\n")

# Subset to fibrosis-eligible samples (fib_ge3 target, same as Script 76)
fib_valid <- meta$loco_fold_fibrosis != "excluded" &
             !is.na(meta$fib_ge3) &
             meta$fib_ge3 >= 0
meta_fib <- meta[fib_valid]
logcpm_fib <- logcpm[, meta_fib$sample_id]
folds_fib <- unique(meta_fib$loco_fold_fibrosis)
folds_fib <- folds_fib[folds_fib != "excluded" & folds_fib != "NA"]
cat("  Fibrosis-eligible:", nrow(meta_fib), "samples,", length(folds_fib), "folds\n")
cat("  Folds:", paste(folds_fib, collapse = ", "), "\n\n")

# Load V3 per-sample predictions
v3_preds <- fread(file.path(OUTDIR, "v3_proper_cv_results.csv"))
v3_fib <- v3_preds[target == "fib_ge3"]
cat("  V3 fib_ge3 predictions:", nrow(v3_fib), "samples\n")

# Compute observed AUROC from V3 per-sample predictions
roc_obs <- roc(v3_fib$true_label, v3_fib$probability, quiet = TRUE)
AUROC_OBS <- as.numeric(auc(roc_obs))
cat("  Observed V3 F>=3 AUROC:", round(AUROC_OBS, 4), "\n\n")

# =============================================================================
# HELPER: Within-fold feature selection + rank transform
# (Identical to Script 76 to ensure reproducibility)
# =============================================================================
select_features_and_transform <- function(logcpm_train, logcpm_test, k = 3000) {
  gene_var <- apply(logcpm_train, 1, var)
  top_k_genes <- names(sort(gene_var, decreasing = TRUE))[1:min(k, length(gene_var))]

  rank_train <- apply(logcpm_train[top_k_genes, , drop = FALSE], 2, function(x) {
    r <- rank(x, ties.method = "average")
    r / length(r)
  })
  rank_test <- apply(logcpm_test[top_k_genes, , drop = FALSE], 2, function(x) {
    r <- rank(x, ties.method = "average")
    r / length(r)
  })

  list(
    X_train = t(rank_train),
    X_test  = t(rank_test),
    selected_genes = top_k_genes
  )
}

# =============================================================================
# HELPER: Run full LOCO-CV pipeline for fib_ge3
# Returns per-sample data.table with probability column, or NULL on failure
# =============================================================================
run_loco_fib_ge3 <- function(meta_v, logcpm_v, folds, k_features = 3000,
                              labels = NULL) {
  # If labels provided, use those (for permutation); otherwise use meta column
  fold_results <- list()

  for (fold in folds) {
    test_idx  <- which(meta_v$loco_fold_fibrosis == fold)
    train_idx <- which(meta_v$loco_fold_fibrosis != fold)

    if (length(test_idx) < 5 || length(train_idx) < 20) next

    y_train <- if (!is.null(labels)) labels[train_idx] else meta_v$fib_ge3[train_idx]
    y_test  <- if (!is.null(labels)) labels[test_idx]  else meta_v$fib_ge3[test_idx]

    if (length(unique(y_test)) < 2) next
    if (length(unique(y_train)) < 2) next

    feat <- select_features_and_transform(
      logcpm_v[, meta_v$sample_id[train_idx]],
      logcpm_v[, meta_v$sample_id[test_idx]],
      k = k_features
    )

    class_counts <- table(y_train)
    class_wts <- 1 / class_counts
    class_wts <- class_wts / sum(class_wts) * length(class_wts)
    sample_wts <- class_wts[as.character(y_train)]

    fit <- tryCatch({
      cv.glmnet(
        x = feat$X_train, y = factor(y_train),
        family = "binomial", alpha = 0.5, nfolds = 5,
        type.measure = "auc", weights = as.numeric(sample_wts)
      )
    }, error = function(e) NULL)

    if (is.null(fit)) next

    prob <- as.numeric(predict(fit, newx = feat$X_test, s = "lambda.min",
                                type = "response"))

    fold_results[[fold]] <- data.table(
      sample_id  = meta_v$sample_id[test_idx],
      true_label = y_test,
      probability = prob,
      fold = fold
    )
  }

  if (length(fold_results) == 0) return(NULL)
  rbindlist(fold_results, fill = TRUE)
}

# =============================================================================
# ANALYSIS 1: BOOTSTRAP 95% CI (1,000 resamples)
# =============================================================================
cat("=== Analysis 1: Bootstrap 95% CI ===\n")
N_BOOT <- 10000

boot_aurocs <- mclapply(seq_len(N_BOOT), function(b) {
  set.seed(42 + b)
  idx <- sample(nrow(v3_fib), replace = TRUE)
  boot_df <- v3_fib[idx]
  # Ensure both classes present
  if (length(unique(boot_df$true_label)) < 2) return(NA_real_)
  tryCatch(
    as.numeric(auc(roc(boot_df$true_label, boot_df$probability, quiet = TRUE))),
    error = function(e) NA_real_
  )
}, mc.cores = NCORES)

boot_aurocs <- unlist(boot_aurocs)
boot_aurocs <- boot_aurocs[!is.na(boot_aurocs)]

bootstrap_ci <- data.table(
  observed_auroc = AUROC_OBS,
  n_bootstrap = N_BOOT,
  n_valid = length(boot_aurocs),
  ci_lower = quantile(boot_aurocs, 0.025),
  ci_upper = quantile(boot_aurocs, 0.975),
  boot_mean = mean(boot_aurocs),
  boot_sd   = sd(boot_aurocs),
  boot_median = median(boot_aurocs)
)

fwrite(bootstrap_ci, file.path(OUTDIR, "sensitivity_bootstrap_ci.csv"))
cat("  Observed AUROC:", round(AUROC_OBS, 4), "\n")
cat("  Bootstrap 95% CI: [", round(bootstrap_ci$ci_lower, 4), ",",
    round(bootstrap_ci$ci_upper, 4), "]\n")
cat("  Bootstrap mean:", round(bootstrap_ci$boot_mean, 4),
    " SD:", round(bootstrap_ci$boot_sd, 4), "\n\n")

# =============================================================================
# ANALYSIS 2: PERMUTATION NULL (10,000 permutations)
# =============================================================================
cat("=== Analysis 2: Permutation Null (10000 permutations x 6 folds) ===\n")
N_PERM <- 10000

perm_aurocs <- mclapply(seq_len(N_PERM), function(p) {
  set.seed(42 + 10000 + p)
  # Shuffle fib_ge3 labels across all fibrosis-eligible samples
  shuffled_labels <- sample(meta_fib$fib_ge3)

  res <- run_loco_fib_ge3(meta_fib, logcpm_fib, folds_fib,
                            k_features = 3000, labels = shuffled_labels)
  if (is.null(res) || nrow(res) == 0) return(NA_real_)
  if (length(unique(res$true_label)) < 2) return(NA_real_)

  tryCatch(
    as.numeric(auc(roc(res$true_label, res$probability, quiet = TRUE))),
    error = function(e) NA_real_
  )
}, mc.cores = NCORES)

perm_aurocs <- unlist(perm_aurocs)
perm_valid <- perm_aurocs[!is.na(perm_aurocs)]

perm_pvalue <- mean(perm_valid >= AUROC_OBS)

perm_results <- data.table(
  observed_auroc = AUROC_OBS,
  n_permutations = N_PERM,
  n_valid = length(perm_valid),
  perm_mean_auroc = mean(perm_valid),
  perm_sd_auroc   = sd(perm_valid),
  perm_max_auroc  = max(perm_valid),
  perm_min_auroc  = min(perm_valid),
  perm_p_value = perm_pvalue
)

# Also save per-permutation AUROCs for histogram plotting
perm_detail <- data.table(
  permutation = seq_along(perm_aurocs),
  auroc = perm_aurocs,
  is_valid = !is.na(perm_aurocs)
)

perm_out <- rbind(
  perm_results,
  data.table(
    observed_auroc = NA, n_permutations = NA, n_valid = NA,
    perm_mean_auroc = NA, perm_sd_auroc = NA, perm_max_auroc = NA,
    perm_min_auroc = NA, perm_p_value = NA
  )[0]  # dummy to ensure schema
)

fwrite(perm_results, file.path(OUTDIR, "sensitivity_permutation.csv"))
fwrite(perm_detail, file.path(OUTDIR, "sensitivity_permutation_detail.csv"))

cat("  Observed AUROC:", round(AUROC_OBS, 4), "\n")
cat("  Permutation null mean:", round(mean(perm_valid), 4),
    " SD:", round(sd(perm_valid), 4), "\n")
cat("  Permutation null max:", round(max(perm_valid), 4), "\n")
cat("  Permutation p-value:", perm_pvalue, "\n")
cat("  (", length(perm_valid), "/", N_PERM, "valid permutations)\n\n")

# =============================================================================
# ANALYSIS 3: FEATURE K SWEEP
# =============================================================================
cat("=== Analysis 3: Feature K Sweep ===\n")
K_VALUES <- c(100, 200, 500, 1000, 2000, 3000, 5000)

k_sweep_results <- list()

for (k_val in K_VALUES) {
  cat("  K =", k_val, "... ")
  t_start <- Sys.time()

  res <- run_loco_fib_ge3(meta_fib, logcpm_fib, folds_fib,
                            k_features = k_val)

  if (is.null(res) || nrow(res) == 0) {
    cat("FAILED\n")
    k_sweep_results[[as.character(k_val)]] <- data.table(
      k = k_val, auroc = NA_real_, n_samples = 0, n_folds = 0,
      elapsed_sec = as.numeric(difftime(Sys.time(), t_start, units = "secs"))
    )
    next
  }

  auroc_k <- tryCatch(
    as.numeric(auc(roc(res$true_label, res$probability, quiet = TRUE))),
    error = function(e) NA_real_
  )

  # Per-fold AUROCs
  fold_aurocs_k <- res[, {
    if (length(unique(true_label)) >= 2) {
      .(auroc = tryCatch(as.numeric(auc(roc(true_label, probability, quiet = TRUE))),
                          error = function(e) NA_real_))
    } else {
      .(auroc = NA_real_)
    }
  }, by = fold]

  elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))

  k_sweep_results[[as.character(k_val)]] <- data.table(
    k = k_val,
    auroc = auroc_k,
    mean_fold_auroc = mean(fold_aurocs_k$auroc, na.rm = TRUE),
    sd_fold_auroc = sd(fold_aurocs_k$auroc, na.rm = TRUE),
    n_samples = nrow(res),
    n_folds = nrow(fold_aurocs_k),
    elapsed_sec = elapsed
  )

  cat("AUROC =", round(auroc_k, 4), " (", round(elapsed, 1), "s)\n")
}

k_sweep <- rbindlist(k_sweep_results, fill = TRUE)
fwrite(k_sweep, file.path(OUTDIR, "sensitivity_k_sweep.csv"))

cat("\n  K sweep summary:\n")
print(k_sweep[, .(k, auroc, mean_fold_auroc, sd_fold_auroc, n_folds)])
cat("\n")

# =============================================================================
# ANALYSIS 4: DELONG TEST VS PUBLISHED PANELS
# =============================================================================
cat("=== Analysis 4: DeLong Test vs Published Panels ===\n")

# The published_panel_benchmark.csv from Script 73 has mean_auroc and per_fold_auroc
# but NOT per-sample predictions. We re-run the published panels using the same
# LOCO fold structure and the same expression data so we can get per-sample
# probabilities for a valid DeLong test.

panel_file <- file.path(OUTDIR, "published_panel_benchmark.csv")
delong_results <- list()

if (file.exists(panel_file)) {
  panel_bench <- fread(panel_file)
  cat("  Published panels found:", nrow(panel_bench), "\n")

  # Load gene symbol -> Ensembl mapping from multi-evidence atlas
  atlas_file <- file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
  if (file.exists(atlas_file)) {
    atlas <- fread(atlas_file, select = c("gene", "gene_name"))
    # gene is Ensembl (versioned), gene_name is symbol
    # Also try unversioned
    atlas[, gene_base := sub("\\.[0-9]+$", "", gene)]
    symbol_map <- setNames(atlas$gene, atlas$gene_name)

    # Also map from expression matrix row names
    expr_genes <- rownames(logcpm_fib)
    expr_base <- sub("\\.[0-9]+$", "", expr_genes)
    names(expr_base) <- expr_genes
  } else {
    cat("  WARNING: Multi-evidence atlas not found, skipping DeLong\n")
    symbol_map <- character(0)
  }

  # Known HGNC alias -> GENCODE v49 name (from Script 73)
  aliases <- c(
    "CTGF" = "CCN2", "NOV" = "CCN3", "CYR61" = "CCN1",
    "WISP1" = "CCN4", "WISP2" = "CCN5"
  )

  # Helper: map gene symbols to expression matrix row IDs
  map_symbols_to_rows <- function(symbols) {
    mapped <- character(0)
    for (sym in symbols) {
      # Try alias first
      query_sym <- if (sym %in% names(aliases)) aliases[[sym]] else sym
      # Try direct symbol mapping
      if (query_sym %in% names(symbol_map)) {
        ensembl_id <- symbol_map[[query_sym]]
        if (ensembl_id %in% rownames(logcpm_fib)) {
          mapped <- c(mapped, ensembl_id)
          next
        }
      }
      # Try grep on rownames (partial Ensembl match)
      # Skip — too slow and unreliable
    }
    mapped
  }

  # For each published panel, run LOCO-CV to get per-sample predictions
  for (i in seq_len(nrow(panel_bench))) {
    panel_name <- panel_bench$panel_name[i]
    cat("  Panel:", panel_name, "... ")

    # Get panel gene symbols from the minimal_panel files
    panel_csv <- file.path(OUTDIR, paste0("minimal_panel_",
      gsub(".*_(\\d+)$", "\\1", panel_name), ".csv"))

    # For published panels (SteatoSITE, Govaere), we need the gene lists
    # These are defined in Script 73; we extract from the panel benchmark
    # For "Our_*" panels, they come from minimal_panel_*.csv
    if (grepl("^Our_", panel_name)) {
      n_genes <- as.integer(sub("Our_", "", panel_name))
      panel_csv <- file.path(OUTDIR, paste0("minimal_panel_", n_genes, ".csv"))
      if (file.exists(panel_csv)) {
        panel_genes_dt <- fread(panel_csv)
        # Try column names: gene, ensembl_id, gene_symbol
        if ("gene" %in% names(panel_genes_dt)) {
          panel_rows <- panel_genes_dt$gene[panel_genes_dt$gene %in% rownames(logcpm_fib)]
        } else {
          panel_rows <- character(0)
        }
      } else {
        cat("SKIP (panel CSV not found)\n")
        next
      }
    } else {
      # Published panels: we only have symbol info from benchmark CSV
      # Missing genes listed in the CSV; we need the full panel list
      # Reconstruct from panel_name
      if (panel_name == "SteatoSITE_15") {
        # 15-gene SteatoSITE panel (Govaere et al. 2020, Hepatology)
        panel_syms <- c("MT1F", "KCNH7", "RASD2", "GDNF", "PRRX1", "FGF7",
                         "LCNL1", "CHRDL2", "POU4F1", "CTGF", "CCL20", "PNPLA3",
                         "GDF15", "SEMA3E", "SPP1")
      } else if (panel_name == "Govaere_25") {
        # 25-gene Govaere NASH CRN staging (Govaere et al. 2020, J Hepatol)
        panel_syms <- c("DUSP6", "A2M", "CDH2", "POSTN", "SPARC", "TAGLN",
                         "THY1", "COL1A1", "COL1A2", "COL3A1", "COL6A3",
                         "ADAMTS2", "FAP", "TGFB1", "TGFBI", "TIMP1",
                         "MMP2", "MMP7", "CTGF", "GDF15", "SPP1",
                         "CCL20", "CXCL8", "IL32", "ICAM1")
      } else {
        cat("SKIP (unknown panel)\n")
        next
      }
      panel_rows <- map_symbols_to_rows(panel_syms)
    }

    n_found <- length(panel_rows)
    cat(n_found, "genes mapped ... ")

    if (n_found < 3) {
      cat("SKIP (too few genes)\n")
      next
    }

    # Run LOCO-CV with FIXED gene panel (no variance selection)
    panel_preds <- list()
    for (fold in folds_fib) {
      test_idx  <- which(meta_fib$loco_fold_fibrosis == fold)
      train_idx <- which(meta_fib$loco_fold_fibrosis != fold)
      if (length(test_idx) < 5 || length(train_idx) < 20) next

      y_train <- meta_fib$fib_ge3[train_idx]
      y_test  <- meta_fib$fib_ge3[test_idx]
      if (length(unique(y_test)) < 2 || length(unique(y_train)) < 2) next

      # Rank-transform the panel genes
      train_mat <- logcpm_fib[panel_rows, meta_fib$sample_id[train_idx], drop = FALSE]
      test_mat  <- logcpm_fib[panel_rows, meta_fib$sample_id[test_idx], drop = FALSE]

      rank_train <- apply(train_mat, 2, function(x) {
        r <- rank(x, ties.method = "average"); r / length(r)
      })
      rank_test <- apply(test_mat, 2, function(x) {
        r <- rank(x, ties.method = "average"); r / length(r)
      })

      X_train <- t(rank_train)
      X_test  <- t(rank_test)

      class_counts <- table(y_train)
      class_wts <- 1 / class_counts
      class_wts <- class_wts / sum(class_wts) * length(class_wts)
      sample_wts <- class_wts[as.character(y_train)]

      fit <- tryCatch({
        cv.glmnet(
          x = X_train, y = factor(y_train),
          family = "binomial", alpha = 0.5, nfolds = 5,
          type.measure = "auc", weights = as.numeric(sample_wts)
        )
      }, error = function(e) NULL)

      if (is.null(fit)) next

      prob <- as.numeric(predict(fit, newx = X_test, s = "lambda.min",
                                  type = "response"))
      panel_preds[[fold]] <- data.table(
        sample_id = meta_fib$sample_id[test_idx],
        true_label = y_test,
        probability = prob,
        fold = fold
      )
    }

    if (length(panel_preds) == 0) {
      cat("FAILED (no valid folds)\n")
      next
    }

    panel_df <- rbindlist(panel_preds)
    panel_auroc <- tryCatch(
      as.numeric(auc(roc(panel_df$true_label, panel_df$probability, quiet = TRUE))),
      error = function(e) NA_real_
    )

    # DeLong test: V3 vs this panel on overlapping samples
    overlap_ids <- intersect(v3_fib$sample_id, panel_df$sample_id)
    if (length(overlap_ids) >= 30) {
      v3_sub   <- v3_fib[match(overlap_ids, sample_id)]
      panel_sub <- panel_df[match(overlap_ids, sample_id)]

      roc_v3    <- roc(v3_sub$true_label, v3_sub$probability, quiet = TRUE)
      roc_panel <- roc(panel_sub$true_label, panel_sub$probability, quiet = TRUE)

      delong <- tryCatch({
        test_res <- roc.test(roc_v3, roc_panel, method = "delong")
        data.table(
          panel = panel_name,
          n_overlap = length(overlap_ids),
          v3_auroc = as.numeric(auc(roc_v3)),
          panel_auroc_refit = panel_auroc,
          panel_auroc_reported = panel_bench$mean_auroc[i],
          n_panel_genes = n_found,
          delong_z = test_res$statistic,
          delong_p = test_res$p.value,
          v3_better = as.numeric(auc(roc_v3)) > panel_auroc
        )
      }, error = function(e) {
        data.table(
          panel = panel_name,
          n_overlap = length(overlap_ids),
          v3_auroc = as.numeric(auc(roc_v3)),
          panel_auroc_refit = panel_auroc,
          panel_auroc_reported = panel_bench$mean_auroc[i],
          n_panel_genes = n_found,
          delong_z = NA_real_,
          delong_p = NA_real_,
          v3_better = as.numeric(auc(roc_v3)) > panel_auroc
        )
      })

      delong_results[[panel_name]] <- delong
      cat("AUROC=", round(panel_auroc, 3),
          " DeLong p=", signif(delong$delong_p, 3), "\n")
    } else {
      cat("AUROC=", round(panel_auroc, 3), " (too few overlap for DeLong)\n")
      delong_results[[panel_name]] <- data.table(
        panel = panel_name,
        n_overlap = length(overlap_ids),
        v3_auroc = AUROC_OBS,
        panel_auroc_refit = panel_auroc,
        panel_auroc_reported = panel_bench$mean_auroc[i],
        n_panel_genes = n_found,
        delong_z = NA_real_,
        delong_p = NA_real_,
        v3_better = AUROC_OBS > panel_auroc
      )
    }
  }
} else {
  cat("  published_panel_benchmark.csv not found — skipping DeLong\n")
}

if (length(delong_results) > 0) {
  delong_dt <- rbindlist(delong_results, fill = TRUE)
  fwrite(delong_dt, file.path(OUTDIR, "delong_comparisons.csv"))
  cat("\n  DeLong summary:\n")
  print(delong_dt[, .(panel, v3_auroc, panel_auroc_refit, delong_p, v3_better)])
} else {
  # Write empty file so downstream knows it ran
  fwrite(data.table(panel = character(0)), file.path(OUTDIR, "delong_comparisons.csv"))
  cat("  No DeLong comparisons computed\n")
}

# =============================================================================
# FINAL SUMMARY
# =============================================================================
cat("\n=== Summary ===\n")
cat("1. Bootstrap 95% CI: [", round(bootstrap_ci$ci_lower, 4), ",",
    round(bootstrap_ci$ci_upper, 4), "] (observed:", round(AUROC_OBS, 4), ")\n")
cat("2. Permutation p-value:", perm_pvalue, "(null mean:",
    round(mean(perm_valid), 4), ")\n")
cat("3. K sweep best:", k_sweep[which.max(auroc), .(k, auroc)], "\n")
cat("4. DeLong tests:", length(delong_results), "comparisons\n")

cat("\nOutputs written to:", OUTDIR, "\n")
cat("  - sensitivity_bootstrap_ci.csv\n")
cat("  - sensitivity_permutation.csv\n")
cat("  - sensitivity_permutation_detail.csv\n")
cat("  - sensitivity_k_sweep.csv\n")
cat("  - delong_comparisons.csv\n")

cat("\n=== 80_classifier_sensitivity.R completed:", as.character(Sys.time()), "===\n")
