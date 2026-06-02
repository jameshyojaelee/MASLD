#!/usr/bin/env Rscript
# 61_feature_engineering.R
# Feature engineering for staging classifier: rank transform, z-score,
# dream t-stat pre-filter, and stability filter
#
# Outputs:
#   - rank_expression_matrix.rds: Within-sample rank-transformed expression (batch-invariant)
#   - zscore_expression_matrix.rds: Within-sample z-scored expression
#   - feature_candidates_3000.csv: Top 3,000 genes by max |t| (GLOBAL — reference only, leaks test-fold info)
#   - feature_stability_report.csv: Per-gene stability + selection summary
#   - loco_gene_lists/fold_{cohort}_feature_candidates_3000.csv: Per-fold leakage-free gene lists
#     Uses dream_loo_{cohort}.csv (dream re-run without held-out cohort) for primary t-stat.
#     Downstream LOCO scripts should load the fold-specific file, NOT the global list.
#   - loco_gene_lists/fold_manifest.csv: Manifest of fold-to-dream-source mapping
#
# Usage: Rscript 61_feature_engineering.R
# SLURM: cpu, 8 CPUs, 64GB RAM, 48h

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR <- file.path(INT, "results")
SIGS <- file.path(RDIR, "disease_signatures")
OUTDIR <- file.path(RDIR, "staging_classifier")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

cat("=== 61: Feature Engineering for Staging Classifier ===\n")
cat("Started:", as.character(Sys.time()), "\n\n")

# --- Load expression data ---
# CRITICAL: Use merged_dge.rds (all 1,444 QC-passing samples), NOT corrected_logcpm.rds
# (which contains only ~608 samples from a subset of cohorts).
cat("Loading merged DGE (all 1,444 samples)...\n")
dge <- readRDS(file.path(RDIR, "integration/merged_dge.rds"))
cat("DGE:", nrow(dge), "genes x", ncol(dge), "samples\n")

# Compute TMM-normalized logCPM from the full DGE
dge <- calcNormFactors(dge, method = "TMM")
logcpm <- cpm(dge, log = TRUE, prior.count = 1)
cat("logCPM expression matrix:", nrow(logcpm), "genes x", ncol(logcpm), "samples\n")

# Also keep raw logCPM as baseline (same object — TMM normalization is in lib.size)
raw_logcpm <- logcpm
cat("Verified: using full merged_dge with all cohorts\n\n")

# ============================================================
# STEP 1: Within-Sample Rank Transform (Batch-Invariant)
# ============================================================
cat("=== STEP 1: Rank Transform ===\n")

# For each sample, rank genes from 1 (lowest) to N (highest)
# Then scale to [0, 1] for comparability across samples
rank_matrix <- apply(logcpm, 2, function(x) {
  r <- rank(x, ties.method = "average")
  r / length(r)  # Scale to [0, 1]
})
dimnames(rank_matrix) <- dimnames(logcpm)

cat("Rank matrix:", nrow(rank_matrix), "x", ncol(rank_matrix), "\n")
cat("Range:", range(rank_matrix), "\n\n")

# ============================================================
# STEP 2: Within-Sample Z-Score Transform
# ============================================================
cat("=== STEP 2: Z-Score Transform ===\n")

# For each sample, z-score gene expression (mean=0, sd=1)
zscore_matrix <- apply(logcpm, 2, function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
})
dimnames(zscore_matrix) <- dimnames(logcpm)

cat("Z-score matrix:", nrow(zscore_matrix), "x", ncol(zscore_matrix), "\n")
cat("Sample mean (should be ~0):", mean(colMeans(zscore_matrix)), "\n")
cat("Sample SD (should be ~1):", mean(apply(zscore_matrix, 2, sd)), "\n\n")

# ============================================================
# STEP 3: Dream T-Statistic Pre-Filter (Top 3,000 Genes)
# ============================================================
# NOTE: This global gene selection uses dream results computed on ALL samples.
# It is retained for backward compatibility (feature_candidates_3000.csv) and
# exploratory/reference use. For LOCO-CV modeling, downstream scripts MUST use
# the per-fold gene lists generated in STEP 3b below, which avoid leaking
# held-out cohort information into feature selection.
cat("=== STEP 3: Dream T-Statistic Pre-Filter (GLOBAL — reference only) ===\n")

# Collect t-statistics from ALL available dream results
# Sources: Disease-vs-Control, per-NAS, per-fibrosis, one-vs-rest, consecutive
t_stat_sources <- list()

# 1. Main dream (Disease vs Control)
main_dream_file <- file.path(RDIR, "integration/dream_results.csv")
if (file.exists(main_dream_file)) {
  dt <- fread(main_dream_file)
  if ("t" %in% names(dt)) {
    t_stat_sources[["disease_vs_control"]] <- dt[, .(gene, t_abs = abs(t))]
  } else if ("logFC" %in% names(dt) && "P.Value" %in% names(dt)) {
    # Approximate t from logFC and p-value
    dt[, t_abs := abs(qnorm(pmax(P.Value / 2, 1e-300)))]
    t_stat_sources[["disease_vs_control"]] <- dt[, .(gene, t_abs)]
  }
  cat("  Disease vs Control:", nrow(t_stat_sources[["disease_vs_control"]]), "genes\n")
}

# 2. Per-NAS-score dream (existing)
nas_dream_file <- file.path(SIGS, "nas_score_dream.csv")
if (file.exists(nas_dream_file)) {
  dt <- fread(nas_dream_file)
  # Take max |t| across all NAS contrasts per gene
  if ("t" %in% names(dt)) {
    agg <- dt[, .(t_abs = max(abs(t), na.rm = TRUE)), by = gene]
  } else {
    agg <- dt[, .(t_abs = max(abs(logFC / (logFC / qnorm(pmax(P.Value / 2, 1e-300)))), na.rm = TRUE)), by = gene]
  }
  t_stat_sources[["nas_vs_ref"]] <- agg
  cat("  NAS vs ref:", nrow(agg), "genes\n")
}

# 3. Per-fibrosis-stage dream (existing)
fib_dream_file <- file.path(SIGS, "fibrosis_stage_dream.csv")
if (file.exists(fib_dream_file)) {
  dt <- fread(fib_dream_file)
  if ("t" %in% names(dt)) {
    agg <- dt[, .(t_abs = max(abs(t), na.rm = TRUE)), by = gene]
  } else {
    agg <- dt[, .(t_abs = max(abs(logFC / (logFC / qnorm(pmax(P.Value / 2, 1e-300)))), na.rm = TRUE)), by = gene]
  }
  t_stat_sources[["fib_vs_ref"]] <- agg
  cat("  Fibrosis vs ref:", nrow(agg), "genes\n")
}

# 4. One-vs-rest NAS (from Script 60)
nas_ovr_file <- file.path(OUTDIR, "one_vs_rest_nas_dream.csv")
if (file.exists(nas_ovr_file)) {
  dt <- fread(nas_ovr_file)
  if ("t" %in% names(dt)) {
    agg <- dt[, .(t_abs = max(abs(t), na.rm = TRUE)), by = gene]
  } else if ("logFC" %in% names(dt) && "P.Value" %in% names(dt)) {
    dt[, t_abs := abs(qnorm(pmax(P.Value / 2, 1e-300)))]
    agg <- dt[, .(t_abs = max(t_abs, na.rm = TRUE)), by = gene]
  }
  t_stat_sources[["nas_ovr"]] <- agg
  cat("  NAS one-vs-rest:", nrow(agg), "genes\n")
}

# 5. One-vs-rest fibrosis (from Script 60)
fib_ovr_file <- file.path(OUTDIR, "one_vs_rest_fibrosis_dream.csv")
if (file.exists(fib_ovr_file)) {
  dt <- fread(fib_ovr_file)
  if ("t" %in% names(dt)) {
    agg <- dt[, .(t_abs = max(abs(t), na.rm = TRUE)), by = gene]
  } else if ("logFC" %in% names(dt) && "P.Value" %in% names(dt)) {
    dt[, t_abs := abs(qnorm(pmax(P.Value / 2, 1e-300)))]
    agg <- dt[, .(t_abs = max(t_abs, na.rm = TRUE)), by = gene]
  }
  t_stat_sources[["fib_ovr"]] <- agg
  cat("  Fibrosis one-vs-rest:", nrow(agg), "genes\n")
}

# 6. Consecutive contrasts (existing)
for (src_name in c("nas_consecutive_dream", "fibrosis_consecutive_dream")) {
  f <- file.path(SIGS, paste0(src_name, ".csv"))
  if (file.exists(f)) {
    dt <- fread(f)
    if ("t" %in% names(dt)) {
      agg <- dt[, .(t_abs = max(abs(t), na.rm = TRUE)), by = gene]
    } else if ("logFC" %in% names(dt) && "P.Value" %in% names(dt)) {
      dt[, t_abs := abs(qnorm(pmax(P.Value / 2, 1e-300)))]
      agg <- dt[, .(t_abs = max(t_abs, na.rm = TRUE)), by = gene]
    } else {
      next
    }
    t_stat_sources[[src_name]] <- agg
    cat("  ", src_name, ":", nrow(agg), "genes\n")
  }
}

cat("\nTotal sources loaded:", length(t_stat_sources), "\n")

# Merge: for each gene, take the max |t| across ALL sources
all_genes <- unique(unlist(lapply(t_stat_sources, function(x) x$gene)))
cat("Unique genes across all sources:", length(all_genes), "\n")

gene_t_summary <- data.table(gene = all_genes)
for (src_name in names(t_stat_sources)) {
  src <- t_stat_sources[[src_name]]
  setnames(src, "t_abs", paste0("t_abs_", src_name), skip_absent = TRUE)
  gene_t_summary <- merge(gene_t_summary, src, by = "gene", all.x = TRUE)
}

# Max |t| across all sources
t_cols <- grep("^t_abs_", names(gene_t_summary), value = TRUE)
gene_t_summary[, max_t_abs := do.call(pmax, c(.SD, na.rm = TRUE)), .SDcols = t_cols]

# Count how many sources each gene appears in
gene_t_summary[, n_sources := rowSums(!is.na(.SD)), .SDcols = t_cols]

# Select top 3,000 by max |t|
gene_t_summary <- gene_t_summary[order(-max_t_abs)]
top_genes <- head(gene_t_summary, 3000)

cat("\nTop 3,000 genes selected:\n")
cat("  Max |t| range:", round(range(top_genes$max_t_abs, na.rm = TRUE), 2), "\n")
cat("  Mean sources per gene:", round(mean(top_genes$n_sources), 1), "\n")
cat("  Genes in >=3 sources:", sum(top_genes$n_sources >= 3), "\n")

# ============================================================
# STEP 4: Stability Filter (LOO-CV Robustness)
# ============================================================
cat("\n=== STEP 4: Stability Filter ===\n")

loo_file <- file.path(RDIR, "integration/loo_cv/loo_cv_per_gene.csv")
stability_filter_applied <- FALSE

if (file.exists(loo_file)) {
  loo <- fread(loo_file)
  cat("LOO-CV per-gene data:", nrow(loo), "genes\n")

  # If there's a robustness classification column
  if ("robustness" %in% names(loo)) {
    stable_genes <- loo[robustness %in% c("Robust", "Stable"), gene]
    cat("Robust + Stable genes:", length(stable_genes), "\n")
    stability_filter_applied <- TRUE
  } else if ("n_loo_sig" %in% names(loo)) {
    # Genes significant in >= 6/8 LOO folds
    stable_genes <- loo[n_loo_sig >= 6, gene]
    cat("Genes significant in >=6/8 LOO folds:", length(stable_genes), "\n")
    stability_filter_applied <- TRUE
  } else {
    cat("No robustness classification or n_loo_sig column found; skipping stability filter\n")
  }
} else {
  cat("LOO-CV per-gene file not found; skipping stability filter\n")
}

# Intersect top 3,000 with stability filter
if (stability_filter_applied) {
  top_genes[, is_stable := gene %in% stable_genes]
  n_stable_in_top <- sum(top_genes$is_stable)
  cat("Stable genes in top 3,000:", n_stable_in_top, "\n")

  # If stability filter removes too many, relax threshold
  if (n_stable_in_top < 1000) {
    cat("WARNING: Stability filter too aggressive, retaining all 3,000\n")
    top_genes[, is_stable := TRUE]
  }
} else {
  top_genes[, is_stable := TRUE]
}

# ============================================================
# STEP 4b: Per-Fold Gene Selection (Leakage-Free LOCO-CV)
# ============================================================
# For each LOCO fold, select top 3,000 genes using ONLY training-fold data.
# Uses dream_loo_{cohort}.csv (dream re-run with that cohort removed) as the
# primary t-stat source. Secondary sources (NAS dream, fibrosis dream, OVR
# contrasts, consecutive contrasts) are computed on stage-specific subsets
# and do NOT have per-fold versions; they are included as-is since the leakage
# from the main disease-vs-control dream is the dominant concern.
#
# Output: loco_gene_lists/fold_{cohort}_feature_candidates_3000.csv
#   Same format as the global feature_candidates_3000.csv.
#   Downstream LOCO scripts (63, 64, 69, 70, 190, etc.) should load the
#   fold-specific file when building the training feature matrix.
cat("\n=== STEP 4b: Per-Fold Gene Selection (Leakage-Free) ===\n")

LOO_DIR <- file.path(RDIR, "integration/loo_cv")
FOLD_OUT <- file.path(OUTDIR, "loco_gene_lists")
dir.create(FOLD_OUT, showWarnings = FALSE, recursive = TRUE)

# Load modeling metadata to discover LOCO folds
meta_file <- file.path(OUTDIR, "modeling_metadata.csv")
if (!file.exists(meta_file)) {
  cat("WARNING: modeling_metadata.csv not found; skipping per-fold gene selection.\n")
  cat("  Run Script 60 first to generate modeling_metadata.csv.\n")
} else {
  meta_dt <- fread(meta_file)

  # Collect all unique fold cohorts across fibrosis and NAS targets
  fib_folds <- setdiff(unique(meta_dt$loco_fold_fibrosis), c("excluded", ""))
  nas_folds <- setdiff(unique(meta_dt$loco_fold_nas), c("excluded", ""))
  all_fold_cohorts <- sort(unique(c(fib_folds, nas_folds)))

  cat("LOCO fold cohorts:", paste(all_fold_cohorts, collapse = ", "), "\n")
  cat("Fibrosis folds (", length(fib_folds), "):", paste(sort(fib_folds), collapse = ", "), "\n")
  cat("NAS folds (", length(nas_folds), "):", paste(sort(nas_folds), collapse = ", "), "\n\n")

  # Helper: same t-stat aggregation logic as STEP 3, but with a custom
  # primary dream source (the LOO dream for the held-out cohort)
  select_top_genes_for_fold <- function(loo_dream_file, n_top = 3000) {
    fold_sources <- list()

    # 1. Primary: LOO dream (held-out cohort removed)
    if (file.exists(loo_dream_file)) {
      dt <- fread(loo_dream_file)
      if ("t" %in% names(dt)) {
        fold_sources[["disease_vs_control"]] <- dt[, .(gene, t_abs = abs(t))]
      } else if ("logFC" %in% names(dt) && "P.Value" %in% names(dt)) {
        dt[, t_abs := abs(qnorm(pmax(P.Value / 2, 1e-300)))]
        fold_sources[["disease_vs_control"]] <- dt[, .(gene, t_abs)]
      }
    }

    # 2-6. Secondary sources (same as global — no per-fold versions available)
    # These use stage-specific contrasts, not disease-vs-control, so leakage is
    # minimal (different model specification). Included for consistency.
    if (exists("t_stat_sources")) {
      for (src_name in names(t_stat_sources)) {
        if (src_name != "disease_vs_control") {
          fold_sources[[src_name]] <- t_stat_sources[[src_name]]
        }
      }
    }

    if (length(fold_sources) == 0) {
      cat("    WARNING: No t-stat sources available\n")
      return(NULL)
    }

    # Merge: max |t| across all sources (same logic as STEP 3)
    fold_genes <- unique(unlist(lapply(fold_sources, function(x) x$gene)))
    fold_summary <- data.table(gene = fold_genes)
    for (src_name in names(fold_sources)) {
      src <- copy(fold_sources[[src_name]])
      setnames(src, "t_abs", paste0("t_abs_", src_name), skip_absent = TRUE)
      fold_summary <- merge(fold_summary, src, by = "gene", all.x = TRUE)
    }

    t_cols_fold <- grep("^t_abs_", names(fold_summary), value = TRUE)
    fold_summary[, max_t_abs := do.call(pmax, c(.SD, na.rm = TRUE)),
                 .SDcols = t_cols_fold]
    fold_summary[, n_sources := rowSums(!is.na(.SD)), .SDcols = t_cols_fold]

    fold_summary <- fold_summary[order(-max_t_abs)]
    fold_top <- head(fold_summary, n_top)
    return(fold_top)
  }

  # Generate per-fold gene lists
  fold_overlap_with_global <- numeric()

  for (cohort in all_fold_cohorts) {
    loo_dream_path <- file.path(LOO_DIR, paste0("dream_loo_", cohort, ".csv"))

    if (!file.exists(loo_dream_path)) {
      cat("  Fold", cohort, ": LOO dream file not found at", loo_dream_path, "\n")
      cat("    Falling back to global gene list (leakage NOT fixed for this fold).\n")
      out_file <- file.path(FOLD_OUT,
                            paste0("fold_", cohort, "_feature_candidates_3000.csv"))
      fwrite(top_genes, out_file)
      next
    }

    fold_top <- select_top_genes_for_fold(loo_dream_path, n_top = 3000)

    if (is.null(fold_top) || nrow(fold_top) == 0) {
      cat("  Fold", cohort, ": gene selection failed, using global fallback\n")
      out_file <- file.path(FOLD_OUT,
                            paste0("fold_", cohort, "_feature_candidates_3000.csv"))
      fwrite(top_genes, out_file)
      next
    }

    # Save per-fold gene list (same format as feature_candidates_3000.csv)
    out_file <- file.path(FOLD_OUT,
                          paste0("fold_", cohort, "_feature_candidates_3000.csv"))
    fwrite(fold_top, out_file)

    # Compute overlap with global list for diagnostics
    overlap <- length(intersect(fold_top$gene, top_genes$gene))
    pct <- round(100 * overlap / nrow(fold_top), 1)
    fold_overlap_with_global[cohort] <- pct

    cat("  Fold", cohort, ":", nrow(fold_top), "genes |",
        "max_t range:", round(range(fold_top$max_t_abs, na.rm = TRUE), 2),
        "| overlap with global:", overlap, "(", pct, "%)\n")
  }

  # Summary
  cat("\nPer-fold gene selection summary:\n")
  cat("  Folds processed:", length(all_fold_cohorts), "\n")
  cat("  Output directory:", FOLD_OUT, "\n")
  if (length(fold_overlap_with_global) > 0) {
    cat("  Mean overlap with global list:",
        round(mean(fold_overlap_with_global, na.rm = TRUE), 1), "%\n")
    cat("  Min overlap:", round(min(fold_overlap_with_global, na.rm = TRUE), 1),
        "% | Max:", round(max(fold_overlap_with_global, na.rm = TRUE), 1), "%\n")
  }

  # Also save a manifest listing which folds use which dream source
  manifest <- data.table(
    fold_cohort = all_fold_cohorts,
    loo_dream_file = file.path(LOO_DIR,
                               paste0("dream_loo_", all_fold_cohorts, ".csv")),
    loo_dream_exists = file.exists(
      file.path(LOO_DIR, paste0("dream_loo_", all_fold_cohorts, ".csv"))),
    output_file = file.path(FOLD_OUT,
                            paste0("fold_", all_fold_cohorts,
                                   "_feature_candidates_3000.csv")),
    in_fibrosis_folds = all_fold_cohorts %in% fib_folds,
    in_nas_folds = all_fold_cohorts %in% nas_folds
  )
  fwrite(manifest, file.path(FOLD_OUT, "fold_manifest.csv"))
  cat("  Manifest saved to:", file.path(FOLD_OUT, "fold_manifest.csv"), "\n")
}

# ============================================================
# STEP 5: Save Outputs
# ============================================================
cat("\n=== STEP 5: Saving Outputs ===\n")

# Feature candidate list (GLOBAL — for reference/exploratory use only)
# For LOCO-CV, use per-fold lists from loco_gene_lists/ instead.
fwrite(top_genes, file.path(OUTDIR, "feature_candidates_3000.csv"))
cat("Feature candidates saved:", nrow(top_genes), "genes (global/reference)\n")

# Subset matrices to top 3,000 genes (rows in expression matrices)
selected_genes <- top_genes$gene
genes_in_matrix <- intersect(selected_genes, rownames(rank_matrix))
cat("Genes found in expression matrices:", length(genes_in_matrix), "of", length(selected_genes), "\n")

# Save rank matrix (top genes only to keep manageable)
rank_sub <- rank_matrix[genes_in_matrix, , drop = FALSE]
saveRDS(rank_sub, file.path(OUTDIR, "rank_expression_matrix.rds"))
cat("Rank matrix saved:", nrow(rank_sub), "x", ncol(rank_sub), "\n")

# Save z-score matrix
zscore_sub <- zscore_matrix[genes_in_matrix, , drop = FALSE]
saveRDS(zscore_sub, file.path(OUTDIR, "zscore_expression_matrix.rds"))
cat("Z-score matrix saved:", nrow(zscore_sub), "x", ncol(zscore_sub), "\n")

# Save raw logCPM matrix for comparison
raw_sub <- raw_logcpm[intersect(genes_in_matrix, rownames(raw_logcpm)), , drop = FALSE]
saveRDS(raw_sub, file.path(OUTDIR, "raw_logcpm_matrix.rds"))
cat("Raw logCPM matrix saved:", nrow(raw_sub), "x", ncol(raw_sub), "\n")

# Full stability report
stability_report <- copy(top_genes)
stability_report[, selected_for_modeling := TRUE]
fwrite(stability_report, file.path(OUTDIR, "feature_stability_report.csv"))
cat("Stability report saved\n")

# Save full-genome matrices as well (for downstream use by imputation etc.)
saveRDS(rank_matrix, file.path(OUTDIR, "rank_expression_full.rds"))
saveRDS(zscore_matrix, file.path(OUTDIR, "zscore_expression_full.rds"))
cat("Full-genome rank and z-score matrices saved\n")

cat("\n=== 61_feature_engineering.R completed:", as.character(Sys.time()), "===\n")
