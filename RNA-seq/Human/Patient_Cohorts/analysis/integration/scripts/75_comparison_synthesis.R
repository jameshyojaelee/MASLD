#!/usr/bin/env Rscript
# 75_comparison_synthesis.R
# Head-to-head comparison of all three staging classifier plans + visualization.
#
# Phase 6 (Comparison) -- aggregates and ranks all models, generates
# publication-quality multi-panel figures and final summary tables.
#
# Inputs (all from OUTDIR = results/staging_classifier/):
#   Plan 1: ordinal_loco_fibrosis.csv, ordinal_loco_nas.csv, binary_threshold_results.csv
#           plan1_sweep_fibrosis.csv, plan1_sweep_nas.csv, plan1_sweep_binary.csv
#   Plan 2: tier1_loco_results.csv, tier2_loco_results.csv, tier3a_nas_results.csv,
#           tier3b_fibrosis_results.csv, hierarchical_neural_results.csv,
#           svm_tier_results.csv, cascade_error_analysis.csv
#   Plan 3: embedding_model_results.csv, embedding_vs_raw_comparison.csv,
#           transformer_loco_results.csv, fusion_demo_results.csv
#   Phase 5: panel_performance_curve.csv, published_panel_benchmark.csv,
#            imputation_validation.csv, conformal_prediction_sets.csv,
#            stability_ranking_all_genes.csv
#
# Outputs:
#   - plan_comparison_summary.csv
#   - model_ranking.csv
#   - final_best_model.csv
#   - comparison_figures.pdf
#   - gene_panel_overlap.csv
#   - stage_classifier_final_report.csv
#
# Usage: Rscript 75_comparison_synthesis.R
# SLURM: cpu, 4 CPUs, 32GB RAM, 48h

suppressPackageStartupMessages({
  library(data.table)
  library(pROC)
  library(ggplot2)
  library(grid)
  library(gridExtra)
})

set.seed(42)

# =============================================================================
# PATHS
# =============================================================================
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR <- file.path(INT, "results")
OUTDIR <- file.path(RDIR, "staging_classifier")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

cat("=== 75: Comparison Synthesis Across All Plans ===\n")
cat("Started:", as.character(Sys.time()), "\n\n")

# =============================================================================
# PUBLICATION THEME (consistent with project convention)
# =============================================================================
theme_pub <- theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95"),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 11),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 9),
    legend.text = element_text(size = 8)
  )

# =============================================================================
# HELPER: Safe file reader (returns NULL if file missing)
# =============================================================================
safe_fread <- function(path, desc = "") {
  if (file.exists(path)) {
    dt <- fread(path)
    cat("  Loaded ", desc, ": ", nrow(dt), " rows x ", ncol(dt), " cols\n", sep = "")
    return(dt)
  } else {
    cat("  MISSING: ", desc, " (", basename(path), ")\n", sep = "")
    return(NULL)
  }
}

# =============================================================================
# HELPER: QWK (consistent with Scripts 63/67)
# =============================================================================
compute_qwk <- function(y_true, y_pred) {
  labels <- sort(unique(c(y_true, y_pred)))
  n <- length(labels)
  if (n < 2) return(NA_real_)

  O <- matrix(0, nrow = n, ncol = n)
  for (i in seq_along(y_true)) {
    r <- which(labels == y_true[i])
    c <- which(labels == y_pred[i])
    if (length(r) == 1 && length(c) == 1) O[r, c] <- O[r, c] + 1
  }
  W <- outer(seq_len(n), seq_len(n), function(i, j) (i - j)^2 / (n - 1)^2)
  row_sums <- rowSums(O); col_sums <- colSums(O); total <- sum(O)
  if (total == 0) return(NA_real_)
  E <- outer(row_sums, col_sums) / total
  num <- sum(W * O); den <- sum(W * E)
  if (den == 0) return(1.0)
  return(1 - num / den)
}

# =============================================================================
# SECTION 1: COLLECT ALL MODEL RESULTS INTO UNIFIED TABLE
# =============================================================================
cat("=== SECTION 1: Loading All Results ===\n\n")

unified <- list()

# --- Plan 1: Ordinal elastic net (Script 63) ---
p1_fib <- safe_fread(file.path(OUTDIR, "ordinal_loco_fibrosis.csv"), "P1 ordinal fibrosis")
p1_nas <- safe_fread(file.path(OUTDIR, "ordinal_loco_nas.csv"), "P1 ordinal NAS")
p1_bin <- safe_fread(file.path(OUTDIR, "binary_threshold_results.csv"), "P1 binary threshold")
p1_summary <- safe_fread(file.path(OUTDIR, "ordinal_model_summary.csv"), "P1 ordinal summary")

# --- Plan 1: Python model sweep (Script 64) ---
p1_sweep_fib <- safe_fread(file.path(OUTDIR, "plan1_sweep_fibrosis.csv"), "P1 sweep fibrosis")
p1_sweep_nas <- safe_fread(file.path(OUTDIR, "plan1_sweep_nas.csv"), "P1 sweep NAS")
p1_sweep_bin <- safe_fread(file.path(OUTDIR, "plan1_sweep_binary.csv"), "P1 sweep binary")

# --- Plan 2: Hierarchical cascade (Scripts 65-68) ---
p2_t1 <- safe_fread(file.path(OUTDIR, "tier1_loco_results.csv"), "P2 tier1")
p2_t2 <- safe_fread(file.path(OUTDIR, "tier2_loco_results.csv"), "P2 tier2")
p2_t3a <- safe_fread(file.path(OUTDIR, "tier3a_nas_results.csv"), "P2 tier3a NAS")
p2_t3b <- safe_fread(file.path(OUTDIR, "tier3b_fibrosis_results.csv"), "P2 tier3b fibrosis")
p2_neural <- safe_fread(file.path(OUTDIR, "hierarchical_neural_results.csv"), "P2 neural")
p2_svm <- safe_fread(file.path(OUTDIR, "svm_tier_results.csv"), "P2 SVM")
p2_cascade <- safe_fread(file.path(OUTDIR, "cascade_error_analysis.csv"), "P2 cascade error")

# --- Plan 3: DL embeddings (Scripts 69-72) ---
p3_emb <- safe_fread(file.path(OUTDIR, "embedding_model_results.csv"), "P3 embedding models")
p3_cmp <- safe_fread(file.path(OUTDIR, "embedding_vs_raw_comparison.csv"), "P3 emb vs raw")
p3_tfm <- safe_fread(file.path(OUTDIR, "transformer_loco_results.csv"), "P3 transformer")
p3_fus <- safe_fread(file.path(OUTDIR, "fusion_demo_results.csv"), "P3 fusion demo")

# --- Phase 5: Panels and imputation ---
panel_perf <- safe_fread(file.path(OUTDIR, "panel_performance_curve.csv"), "Panel performance")
pub_bench <- safe_fread(file.path(OUTDIR, "published_panel_benchmark.csv"), "Published benchmark")
imp_val <- safe_fread(file.path(OUTDIR, "imputation_validation.csv"), "Imputation validation")
conformal <- safe_fread(file.path(OUTDIR, "conformal_prediction_sets.csv"), "Conformal sets")
stability <- safe_fread(file.path(OUTDIR, "stability_ranking_all_genes.csv"), "Stability ranking")

cat("\n")

# =============================================================================
# HELPER: Extract metrics from per-sample prediction tables
# =============================================================================
extract_loco_metrics <- function(dt, model_name, plan, target_type,
                                  pred_col = "predicted", true_col = "true",
                                  prob_col = NULL, fold_col = "fold") {
  if (is.null(dt)) return(NULL)

  # Try to find the right column names
  possible_pred <- c(pred_col, "pred", "predicted_class", "y_pred", "pred_class")
  possible_true <- c(true_col, "actual", "true_class", "y_true", "true_label")
  possible_fold <- c(fold_col, "loco_fold", "dataset", "fold_name")
  possible_prob <- c(prob_col, "prob_1", "pred_prob", "probability")

  pred_col <- intersect(possible_pred, names(dt))[1]
  true_col <- intersect(possible_true, names(dt))[1]
  fold_col <- intersect(possible_fold, names(dt))[1]
  prob_col_found <- intersect(possible_prob, names(dt))
  prob_col <- if (length(prob_col_found) > 0) prob_col_found[1] else NULL

  if (is.na(pred_col) || is.na(true_col) || is.na(fold_col)) return(NULL)

  folds <- unique(dt[[fold_col]])
  folds <- folds[!folds %in% c("excluded", "NA", "")]

  results <- list()
  for (fld in folds) {
    sub <- dt[get(fold_col) == fld]
    y_true <- as.numeric(sub[[true_col]])
    y_pred <- as.numeric(sub[[pred_col]])

    valid <- !is.na(y_true) & !is.na(y_pred)
    y_true <- y_true[valid]
    y_pred <- y_pred[valid]
    if (length(y_true) < 5) next

    acc <- mean(y_true == y_pred)
    qwk <- tryCatch(compute_qwk(y_true, y_pred), error = function(e) NA_real_)
    mae <- mean(abs(y_true - y_pred))

    # AUROC for binary targets
    auroc <- NA_real_
    if (!is.null(prob_col) && prob_col %in% names(sub) &&
        length(unique(y_true)) == 2) {
      prob <- as.numeric(sub[[prob_col]][valid])
      auroc <- tryCatch({
        as.numeric(pROC::auc(pROC::roc(y_true, prob, quiet = TRUE)))
      }, error = function(e) NA_real_)
    }

    results[[length(results) + 1]] <- data.table(
      model = model_name, plan = plan, target = target_type,
      fold = fld, accuracy = acc, qwk = qwk, mae = mae, auroc = auroc,
      n_test = length(y_true)
    )
  }
  if (length(results) == 0) return(NULL)
  rbindlist(results)
}

# =============================================================================
# HELPER: Extract metrics from summary tables (already aggregated per fold)
# =============================================================================
extract_summary_metrics <- function(dt, plan_name) {
  if (is.null(dt)) return(NULL)

  # Try to standardize columns
  possible_model <- c("model", "model_name", "model_type", "algorithm")
  possible_target <- c("target", "target_type", "task", "endpoint")
  possible_fold <- c("fold", "loco_fold", "fold_name", "dataset")
  possible_acc <- c("accuracy", "acc", "test_accuracy")
  possible_qwk <- c("qwk", "kappa", "quadratic_weighted_kappa")
  possible_auroc <- c("auroc", "auc", "roc_auc", "mean_auroc", "test_auroc")

  model_col <- intersect(possible_model, names(dt))[1]
  target_col <- intersect(possible_target, names(dt))[1]
  fold_col <- intersect(possible_fold, names(dt))[1]
  acc_col <- intersect(possible_acc, names(dt))[1]
  qwk_col <- intersect(possible_qwk, names(dt))[1]
  auroc_col <- intersect(possible_auroc, names(dt))[1]

  if (is.na(model_col) && is.na(acc_col)) return(NULL)

  results <- copy(dt)
  results[, plan := plan_name]

  # Rename columns to standard names
  if (!is.na(model_col) && model_col != "model") setnames(results, model_col, "model")
  if (!is.na(target_col) && target_col != "target") setnames(results, target_col, "target")
  if (!is.na(fold_col) && fold_col != "fold") setnames(results, fold_col, "fold")
  if (!is.na(acc_col) && acc_col != "accuracy") setnames(results, acc_col, "accuracy")
  if (!is.na(qwk_col) && qwk_col != "qwk") setnames(results, qwk_col, "qwk")
  if (!is.na(auroc_col) && auroc_col != "auroc") setnames(results, auroc_col, "auroc")

  # Ensure required columns exist
  if (!"model" %in% names(results)) results[, model := "unknown"]
  if (!"target" %in% names(results)) results[, target := "unknown"]
  if (!"fold" %in% names(results)) results[, fold := "aggregated"]
  if (!"accuracy" %in% names(results)) results[, accuracy := NA_real_]
  if (!"qwk" %in% names(results)) results[, qwk := NA_real_]
  if (!"auroc" %in% names(results)) results[, auroc := NA_real_]
  if (!"mae" %in% names(results)) results[, mae := NA_real_]
  if (!"n_test" %in% names(results)) results[, n_test := NA_integer_]

  results[, .(model, plan, target, fold, accuracy, qwk, mae, auroc, n_test)]
}

# =============================================================================
# COLLECT METRICS
# =============================================================================
cat("=== SECTION 2: Extracting Metrics ===\n\n")

all_metrics <- list()

# Plan 1: Ordinal elastic net per-sample predictions
m <- extract_loco_metrics(p1_fib, "Elastic_Net_Ordinal", "Plan1", "fibrosis_5class")
if (!is.null(m)) all_metrics[[length(all_metrics) + 1]] <- m

m <- extract_loco_metrics(p1_nas, "Elastic_Net_Ordinal", "Plan1", "nas_4class")
if (!is.null(m)) all_metrics[[length(all_metrics) + 1]] <- m

m <- extract_loco_metrics(p1_bin, "Elastic_Net_Binary", "Plan1", "f_ge3_binary")
if (!is.null(m)) all_metrics[[length(all_metrics) + 1]] <- m

# Plan 1: Summary table if available
m <- extract_summary_metrics(p1_summary, "Plan1")
if (!is.null(m)) all_metrics[[length(all_metrics) + 1]] <- m

# Plan 1: Python sweep — parse per-fold results from sweep CSVs
for (info in list(
  list(dt = p1_sweep_fib, target = "fibrosis_5class"),
  list(dt = p1_sweep_nas, target = "nas_4class"),
  list(dt = p1_sweep_bin, target = "f_ge3_binary")
)) {
  if (!is.null(info$dt)) {
    m <- extract_summary_metrics(info$dt, "Plan1")
    if (!is.null(m)) {
      m[, target := info$target]
      all_metrics[[length(all_metrics) + 1]] <- m
    }
  }
}

# Plan 2: Per-sample predictions from tiers
m <- extract_loco_metrics(p2_t1, "Tier1_Disease", "Plan2", "disease_binary")
if (!is.null(m)) all_metrics[[length(all_metrics) + 1]] <- m

m <- extract_loco_metrics(p2_t2, "Tier2_Severity", "Plan2", "severity_4class")
if (!is.null(m)) all_metrics[[length(all_metrics) + 1]] <- m

m <- extract_loco_metrics(p2_t3a, "Tier3a_NAS", "Plan2", "nas_4class")
if (!is.null(m)) all_metrics[[length(all_metrics) + 1]] <- m

m <- extract_loco_metrics(p2_t3b, "Tier3b_Fibrosis", "Plan2", "fibrosis_5class")
if (!is.null(m)) all_metrics[[length(all_metrics) + 1]] <- m

# Plan 2: Neural and SVM summaries
m <- extract_summary_metrics(p2_neural, "Plan2")
if (!is.null(m)) all_metrics[[length(all_metrics) + 1]] <- m

m <- extract_summary_metrics(p2_svm, "Plan2")
if (!is.null(m)) all_metrics[[length(all_metrics) + 1]] <- m

# Plan 3: Embedding model results
m <- extract_summary_metrics(p3_emb, "Plan3")
if (!is.null(m)) all_metrics[[length(all_metrics) + 1]] <- m

# Plan 3: Transformer LOCO
m <- extract_loco_metrics(p3_tfm, "Transformer", "Plan3", "fibrosis_5class")
if (!is.null(m)) all_metrics[[length(all_metrics) + 1]] <- m

# Plan 3: Embedding vs raw comparison
m <- extract_summary_metrics(p3_cmp, "Plan3")
if (!is.null(m)) all_metrics[[length(all_metrics) + 1]] <- m

# Plan 3: Fusion
m <- extract_summary_metrics(p3_fus, "Plan3")
if (!is.null(m)) all_metrics[[length(all_metrics) + 1]] <- m

# Combine
if (length(all_metrics) == 0) {
  cat("\n  WARNING: No metrics loaded from any plan.\n")
  cat("  Creating empty output files.\n")

  fwrite(data.table(note = "No model results available"),
         file.path(OUTDIR, "plan_comparison_summary.csv"))
  fwrite(data.table(note = "No model results available"),
         file.path(OUTDIR, "model_ranking.csv"))
  fwrite(data.table(note = "No model results available"),
         file.path(OUTDIR, "final_best_model.csv"))
  fwrite(data.table(note = "No model results available"),
         file.path(OUTDIR, "stage_classifier_final_report.csv"))
  cat("\nFinished (no results):", as.character(Sys.time()), "\n")
  quit(save = "no", status = 0)
}

unified_dt <- rbindlist(all_metrics, fill = TRUE)
cat("\n  Unified metrics table:", nrow(unified_dt), "rows\n")
cat("  Plans:", paste(unique(unified_dt$plan), collapse = ", "), "\n")
cat("  Models:", length(unique(unified_dt$model)), "\n")
cat("  Targets:", paste(unique(unified_dt$target), collapse = ", "), "\n\n")

fwrite(unified_dt, file.path(OUTDIR, "plan_comparison_summary.csv"))
cat("  Saved plan_comparison_summary.csv\n")

# =============================================================================
# SECTION 3: MODEL RANKING
# =============================================================================
cat("\n=== SECTION 3: Model Ranking ===\n")

# Aggregate per model x plan x target
ranking <- unified_dt[, .(
  mean_accuracy = mean(accuracy, na.rm = TRUE),
  sd_accuracy   = sd(accuracy, na.rm = TRUE),
  mean_qwk      = mean(qwk, na.rm = TRUE),
  sd_qwk        = sd(qwk, na.rm = TRUE),
  mean_mae      = mean(mae, na.rm = TRUE),
  mean_auroc    = mean(auroc, na.rm = TRUE),
  sd_auroc      = sd(auroc, na.rm = TRUE),
  n_folds       = .N,
  n_test_total  = sum(n_test, na.rm = TRUE)
), by = .(model, plan, target)]

# Rank within each target
ranking_list <- list()
for (tgt in unique(ranking$target)) {
  sub <- ranking[target == tgt]

  # Primary ranking metric depends on target type
  if (grepl("binary", tgt, ignore.case = TRUE)) {
    setorder(sub, -mean_auroc)
    sub[, rank_by := "auroc"]
  } else {
    # Rank by QWK first, then accuracy as tiebreaker
    setorder(sub, -mean_qwk, -mean_accuracy)
    sub[, rank_by := "qwk"]
  }
  sub[, rank := .I]
  ranking_list[[length(ranking_list) + 1]] <- sub
}

ranking_dt <- rbindlist(ranking_list, fill = TRUE)
fwrite(ranking_dt, file.path(OUTDIR, "model_ranking.csv"))
cat("  Saved model_ranking.csv (", nrow(ranking_dt), " entries)\n")

# Print top models per target
for (tgt in unique(ranking_dt$target)) {
  cat("\n  Target:", tgt, "\n")
  sub <- ranking_dt[target == tgt & rank <= 3,
                    .(rank, model, plan, mean_accuracy, mean_qwk, mean_auroc)]
  print(sub)
}

# =============================================================================
# SECTION 4: BEST MODEL PER TARGET
# =============================================================================
cat("\n=== SECTION 4: Best Model Per Target ===\n")

best_models <- ranking_dt[rank == 1, .(
  target, model, plan, mean_accuracy, mean_qwk, mean_auroc, n_folds
)]

# Overall best (average rank across targets)
if (nrow(ranking_dt) > 0) {
  avg_rank <- ranking_dt[, .(
    mean_rank = mean(rank, na.rm = TRUE),
    targets_present = .N
  ), by = .(model, plan)]
  setorder(avg_rank, mean_rank)

  cat("  Overall ranking (avg rank across targets):\n")
  print(head(avg_rank, 5))
  best_models <- rbind(best_models,
    data.table(target = "OVERALL_BEST", model = avg_rank$model[1],
               plan = avg_rank$plan[1], mean_accuracy = NA_real_,
               mean_qwk = NA_real_, mean_auroc = NA_real_,
               n_folds = NA_integer_),
    fill = TRUE)
}

fwrite(best_models, file.path(OUTDIR, "final_best_model.csv"))
cat("  Saved final_best_model.csv\n")
print(best_models)

# =============================================================================
# SECTION 5: STATISTICAL COMPARISON (DeLong test for top 2 models)
# =============================================================================
cat("\n=== SECTION 5: Statistical Comparison (DeLong) ===\n")

# Only for binary targets with per-sample predictions
delong_results <- list()

# Try to find binary predictions from the top 2 models
binary_targets <- ranking_dt[grepl("binary", target, ignore.case = TRUE)]
if (nrow(binary_targets) >= 2) {
  top2 <- binary_targets[rank <= 2]
  cat("  Top 2 models for binary F>=3: ", top2$model[1], " vs ", top2$model[2], "\n")
  cat("  DeLong test requires paired per-sample probabilities from both models.\n")
  cat("  (Skipping if per-sample prob files not aligned)\n")

  # Attempt DeLong if we have per-sample predictions from both
  # This is a best-effort comparison
  if (!is.null(p1_bin)) {
    prob_cols <- grep("prob|probability", names(p1_bin), value = TRUE)
    true_cols <- intersect(c("true", "actual", "true_class", "y_true"), names(p1_bin))
    if (length(prob_cols) > 0 && length(true_cols) > 0) {
      y_true_bin <- as.numeric(p1_bin[[true_cols[1]]])
      y_prob_bin <- as.numeric(p1_bin[[prob_cols[1]]])
      valid <- !is.na(y_true_bin) & !is.na(y_prob_bin)
      if (sum(valid) > 20 && length(unique(y_true_bin[valid])) == 2) {
        roc1 <- tryCatch(
          pROC::roc(y_true_bin[valid], y_prob_bin[valid], quiet = TRUE),
          error = function(e) NULL)
        if (!is.null(roc1)) {
          delong_results[[1]] <- data.table(
            comparison = "Plan1_binary_self",
            auroc = as.numeric(pROC::auc(roc1)),
            ci_lower = as.numeric(pROC::ci.auc(roc1))[1],
            ci_upper = as.numeric(pROC::ci.auc(roc1))[3],
            n_samples = sum(valid)
          )
          cat("  Plan1 binary AUROC:", sprintf("%.3f [%.3f, %.3f]",
              delong_results[[1]]$auroc,
              delong_results[[1]]$ci_lower,
              delong_results[[1]]$ci_upper), "\n")
        }
      }
    }
  }
} else {
  cat("  Fewer than 2 binary models; skipping DeLong\n")
}

# =============================================================================
# SECTION 6: FEATURE SELECTION COMPARISON
# =============================================================================
cat("\n=== SECTION 6: Feature Selection Comparison ===\n")

gene_panel_overlap <- list()

# Load stability ranking (Plan 1 + 2)
stab_genes <- character(0)
if (!is.null(stability)) {
  stab_genes <- head(stability$gene, 50)
  cat("  Stability top 50:", length(stab_genes), "genes\n")
}

# Load attention-based panel (Plan 3)
attn_genes <- character(0)
attn_file <- file.path(OUTDIR, "attention_panel_50.csv")
if (file.exists(attn_file)) {
  attn_panel <- fread(attn_file)
  attn_genes <- attn_panel$gene
  cat("  Attention top 50:", length(attn_genes), "genes\n")
} else {
  # Try attention weights
  attn_w <- safe_fread(file.path(OUTDIR, "attention_weights_all.csv"), "attention weights")
  if (!is.null(attn_w) && "gene" %in% names(attn_w)) {
    num_cols <- names(attn_w)[sapply(attn_w, is.numeric)]
    if (length(num_cols) > 0) {
      setorderv(attn_w, num_cols[1], order = -1)
      attn_genes <- head(attn_w$gene, 50)
    }
  }
}

# Load tier-1 features (Plan 2)
tier1_genes <- character(0)
tier1_file <- file.path(OUTDIR, "tier1_top_features.csv")
if (file.exists(tier1_file)) {
  tier1_dt <- fread(tier1_file)
  if ("gene" %in% names(tier1_dt)) {
    tier1_genes <- tier1_dt$gene
    cat("  Tier-1 top features:", length(tier1_genes), "genes\n")
  }
}

# Build overlap table for UpSet plot
all_union <- unique(c(stab_genes, attn_genes, tier1_genes))
if (length(all_union) > 0) {
  overlap_dt <- data.table(
    gene = all_union,
    stability_top50 = all_union %in% stab_genes,
    attention_top50 = all_union %in% attn_genes,
    tier1_features  = all_union %in% tier1_genes
  )

  # Overlap stats
  if (length(stab_genes) > 0 && length(attn_genes) > 0) {
    overlap_sa <- length(intersect(stab_genes, attn_genes))
    cat("  Stability-Attention overlap:", overlap_sa, "/",
        length(union(stab_genes, attn_genes)), "\n")
  }
  if (length(stab_genes) > 0 && length(tier1_genes) > 0) {
    overlap_st <- length(intersect(stab_genes, tier1_genes))
    cat("  Stability-Tier1 overlap:", overlap_st, "/",
        length(union(stab_genes, tier1_genes)), "\n")
  }

  fwrite(overlap_dt, file.path(OUTDIR, "gene_panel_overlap.csv"))
  cat("  Saved gene_panel_overlap.csv (", nrow(overlap_dt), " genes)\n")
} else {
  cat("  No gene panels available for overlap analysis\n")
  overlap_dt <- data.table()
}

# =============================================================================
# SECTION 7: VISUALIZATION (PDF)
# =============================================================================
cat("\n=== SECTION 7: Generating Figures ===\n")

pdf_path <- file.path(OUTDIR, "comparison_figures.pdf")
pdf(pdf_path, width = 12, height = 10)

plots_created <- 0

# --- Panel A: LOCO performance heatmap (model x target) ---
if (nrow(ranking_dt) > 0) {
  cat("  Panel A: Performance heatmap\n")

  # Choose primary metric per target type
  hm_data <- ranking_dt[, .(model, plan, target, mean_accuracy, mean_qwk, mean_auroc)]
  hm_data[, display_metric := fifelse(
    grepl("binary", target, ignore.case = TRUE), mean_auroc, mean_qwk
  )]
  hm_data[, metric_label := fifelse(
    grepl("binary", target, ignore.case = TRUE), "AUROC", "QWK"
  )]
  hm_data[, model_plan := paste0(model, " (", plan, ")")]

  # Limit to top 8 models per target for readability
  # Merge rank from ranking_dt, then filter
  hm_data <- merge(hm_data, ranking_dt[, .(model, plan, target, model_rank = rank)],
                   by = c("model", "plan", "target"), all.x = TRUE)
  top_models <- hm_data[!is.na(model_rank) & model_rank <= 8]

  if (nrow(top_models) > 0) {
    p_hm <- ggplot(top_models,
                   aes(x = target, y = reorder(model_plan, display_metric),
                       fill = display_metric)) +
      geom_tile(color = "white", linewidth = 0.5) +
      geom_text(aes(label = sprintf("%.3f", display_metric)), size = 2.5) +
      scale_fill_gradient2(low = "#d73027", mid = "#fee08b", high = "#1a9850",
                           midpoint = 0.5, name = "Performance") +
      labs(title = "A) LOCO Performance (QWK for ordinal, AUROC for binary)",
           x = "Target", y = "Model (Plan)") +
      theme_pub +
      theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 7),
            axis.text.y = element_text(size = 7))

    print(p_hm)
    plots_created <- plots_created + 1
  }
}

# --- Panel B: Panel size vs AUROC performance curve ---
if (!is.null(panel_perf) && nrow(panel_perf) > 0) {
  cat("  Panel B: Panel size vs AUROC\n")

  p_panel <- ggplot(panel_perf, aes(x = panel_size, y = mean_auroc)) +
    geom_line(color = "#1f78b4", linewidth = 1) +
    geom_point(color = "#1f78b4", size = 3) +
    geom_errorbar(aes(ymin = mean_auroc - sd_auroc,
                      ymax = mean_auroc + sd_auroc),
                  width = 0.05 * max(panel_perf$panel_size, na.rm = TRUE),
                  color = "#1f78b4") +
    scale_x_continuous(breaks = panel_perf$panel_size, trans = "log10") +
    coord_cartesian(ylim = c(
      max(0.5, min(panel_perf$mean_auroc - panel_perf$sd_auroc, na.rm = TRUE) - 0.05),
      1.0
    )) +
    labs(title = "B) Gene Panel Size vs AUROC (F>=3 binary)",
         x = "Number of Genes (log scale)", y = "Mean LOCO AUROC") +
    theme_pub

  # Add reference line at 15 genes
  if (15 %in% panel_perf$panel_size) {
    p_panel <- p_panel +
      geom_vline(xintercept = 15, linetype = "dashed", color = "grey50")
  }

  print(p_panel)
  plots_created <- plots_created + 1
}

# --- Panel C: Published panel benchmark barplot ---
if (!is.null(pub_bench) && nrow(pub_bench) > 0 &&
    any(!is.na(pub_bench$mean_auroc))) {
  cat("  Panel C: Published panel benchmark\n")

  pub_plot <- pub_bench[!is.na(mean_auroc)]
  pub_plot[, panel_label := paste0(panel_name, "\n(", n_genes_found, " genes)")]

  p_bench <- ggplot(pub_plot,
                    aes(x = reorder(panel_label, mean_auroc), y = mean_auroc)) +
    geom_col(fill = "#a6cee3", color = "grey30", width = 0.6) +
    geom_errorbar(aes(ymin = mean_auroc - sd_auroc,
                      ymax = mean_auroc + sd_auroc),
                  width = 0.15) +
    geom_text(aes(label = sprintf("%.3f", mean_auroc)), vjust = -0.5, size = 3) +
    coord_cartesian(ylim = c(0.5, 1.0)) +
    labs(title = "C) Published vs Our Panel Benchmark (F>=3 binary)",
         x = "Gene Panel", y = "Mean LOCO AUROC") +
    theme_pub +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))

  print(p_bench)
  plots_created <- plots_created + 1
}

# --- Panel D: Gene panel overlap (UpSet-like horizontal barplot) ---
if (nrow(overlap_dt) > 0 && ncol(overlap_dt) > 1) {
  cat("  Panel D: Gene panel overlap\n")

  # Create intersection category
  set_cols <- setdiff(names(overlap_dt), "gene")
  overlap_dt[, category := apply(.SD, 1, function(x) {
    paste(set_cols[x], collapse = " + ")
  }), .SDcols = set_cols]
  overlap_dt[category == "", category := "None"]

  cat_counts <- overlap_dt[, .N, by = category]
  setorder(cat_counts, -N)

  p_upset <- ggplot(cat_counts, aes(x = reorder(category, N), y = N)) +
    geom_col(fill = "#b2df8a", color = "grey30", width = 0.6) +
    geom_text(aes(label = N), hjust = -0.2, size = 3) +
    coord_flip() +
    labs(title = "D) Gene Panel Overlap (Stability / Attention / Tier-1)",
         x = "Intersection", y = "Number of Genes") +
    theme_pub

  print(p_upset)
  plots_created <- plots_created + 1
}

# --- Panel E: Conformal prediction set width distribution ---
if (!is.null(conformal) && nrow(conformal) > 0 &&
    "mean_set_size" %in% names(conformal)) {
  cat("  Panel E: Conformal prediction set sizes\n")

  p_conf <- ggplot(conformal, aes(x = mean_set_size)) +
    geom_histogram(binwidth = 0.5, fill = "#cab2d6", color = "grey30",
                   boundary = 0.5) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
    labs(title = "E) Conformal Prediction Set Width (90% coverage)",
         x = "Mean Prediction Set Size", y = "Number of Samples") +
    theme_pub

  if ("confidence_level" %in% names(conformal)) {
    # Add confidence coloring
    conf_summary <- conformal[, .N, by = confidence_level]
    subtitle_text <- paste(
      conf_summary$confidence_level, ":", conf_summary$N,
      collapse = " | "
    )
    p_conf <- p_conf + labs(subtitle = subtitle_text)
  }

  print(p_conf)
  plots_created <- plots_created + 1
}

# --- Panel F: Plan comparison overview (bar chart of best per plan) ---
if (nrow(ranking_dt) > 0) {
  cat("  Panel F: Plan-level comparison\n")

  # Best model per plan per target
  plan_best <- ranking_dt[, .SD[which.min(rank)], by = .(plan, target)]
  plan_best[, display_metric := fifelse(
    grepl("binary", target, ignore.case = TRUE), mean_auroc, mean_qwk
  )]

  if (nrow(plan_best) > 0) {
    p_plan <- ggplot(plan_best,
                     aes(x = target, y = display_metric, fill = plan)) +
      geom_col(position = position_dodge(width = 0.7), width = 0.6,
               color = "grey30") +
      geom_text(aes(label = sprintf("%.3f", display_metric)),
                position = position_dodge(width = 0.7), vjust = -0.5, size = 2.5) +
      scale_fill_brewer(palette = "Set2", name = "Plan") +
      coord_cartesian(ylim = c(0, 1.05)) +
      labs(title = "F) Best Model per Plan (QWK for ordinal, AUROC for binary)",
           x = "Target", y = "Performance") +
      theme_pub +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))

    print(p_plan)
    plots_created <- plots_created + 1
  }
}

dev.off()
cat("  Saved comparison_figures.pdf (", plots_created, " panels)\n")

# =============================================================================
# SECTION 8: FINAL SUMMARY TABLE (for manuscript)
# =============================================================================
cat("\n=== SECTION 8: Final Report Table ===\n")

# Build one-row-per-model/plan summary
report_rows <- list()
for (tgt in unique(ranking_dt$target)) {
  sub <- ranking_dt[target == tgt]
  for (i in seq_len(nrow(sub))) {
    report_rows[[length(report_rows) + 1]] <- data.table(
      model = sub$model[i],
      plan = sub$plan[i],
      target = tgt,
      accuracy = sub$mean_accuracy[i],
      accuracy_sd = sub$sd_accuracy[i],
      qwk = sub$mean_qwk[i],
      qwk_sd = sub$sd_qwk[i],
      auroc = sub$mean_auroc[i],
      auroc_sd = sub$sd_auroc[i],
      mae = sub$mean_mae[i],
      n_folds = sub$n_folds[i],
      rank_in_target = sub$rank[i]
    )
  }
}

# Add panel information from Phase 5
if (!is.null(panel_perf) && nrow(panel_perf) > 0) {
  for (i in seq_len(nrow(panel_perf))) {
    report_rows[[length(report_rows) + 1]] <- data.table(
      model = paste0("Panel_", panel_perf$panel_size[i]),
      plan = "Phase5_RFE",
      target = "f_ge3_binary",
      accuracy = NA_real_, accuracy_sd = NA_real_,
      qwk = NA_real_, qwk_sd = NA_real_,
      auroc = panel_perf$mean_auroc[i],
      auroc_sd = panel_perf$sd_auroc[i],
      mae = NA_real_,
      n_folds = panel_perf$n_folds[i],
      rank_in_target = NA_integer_
    )
  }
}

# Add published benchmark
if (!is.null(pub_bench) && nrow(pub_bench) > 0) {
  for (i in seq_len(nrow(pub_bench))) {
    report_rows[[length(report_rows) + 1]] <- data.table(
      model = pub_bench$panel_name[i],
      plan = "Published",
      target = "f_ge3_binary",
      accuracy = NA_real_, accuracy_sd = NA_real_,
      qwk = NA_real_, qwk_sd = NA_real_,
      auroc = pub_bench$mean_auroc[i],
      auroc_sd = pub_bench$sd_auroc[i],
      mae = NA_real_,
      n_folds = NA_integer_,
      rank_in_target = NA_integer_
    )
  }
}

report_dt <- rbindlist(report_rows, fill = TRUE)

# Add n_genes_panel column (where applicable)
report_dt[, n_genes_panel := NA_integer_]
if (!is.null(panel_perf)) {
  for (i in seq_len(nrow(panel_perf))) {
    report_dt[model == paste0("Panel_", panel_perf$panel_size[i]),
              n_genes_panel := as.integer(panel_perf$panel_size[i])]
  }
}
if (!is.null(pub_bench)) {
  for (i in seq_len(nrow(pub_bench))) {
    report_dt[model == pub_bench$panel_name[i],
              n_genes_panel := as.integer(pub_bench$n_genes_found[i])]
  }
}

fwrite(report_dt, file.path(OUTDIR, "stage_classifier_final_report.csv"))
cat("  Saved stage_classifier_final_report.csv (", nrow(report_dt), " entries)\n")

# Print top results
cat("\n  === TOP RESULTS SUMMARY ===\n\n")
for (tgt in unique(ranking_dt$target)) {
  top <- ranking_dt[target == tgt & rank == 1]
  if (nrow(top) == 0) next
  metric_str <- if (grepl("binary", tgt)) {
    sprintf("AUROC=%.3f", top$mean_auroc)
  } else {
    sprintf("QWK=%.3f, Acc=%.3f", top$mean_qwk, top$mean_accuracy)
  }
  cat("  ", tgt, ": ", top$model, " (", top$plan, ") — ", metric_str, "\n", sep = "")
}

if (!is.null(pub_bench) && any(!is.na(pub_bench$mean_auroc))) {
  cat("\n  Published panel benchmark:\n")
  for (i in seq_len(nrow(pub_bench))) {
    cat("    ", pub_bench$panel_name[i], ": AUROC=",
        sprintf("%.3f", pub_bench$mean_auroc[i]), "\n", sep = "")
  }
}

cat("\n  Output files:\n")
cat("    plan_comparison_summary.csv\n")
cat("    model_ranking.csv\n")
cat("    final_best_model.csv\n")
cat("    comparison_figures.pdf\n")
cat("    gene_panel_overlap.csv\n")
cat("    stage_classifier_final_report.csv\n")

cat("\nFinished:", as.character(Sys.time()), "\n")
