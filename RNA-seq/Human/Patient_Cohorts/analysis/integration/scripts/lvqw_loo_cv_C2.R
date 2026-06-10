#!/usr/bin/env Rscript
# ===========================================================================
# lvqw_loo_cv_C2.R
# Item 2 of the C2 validation refresh: cross-cohort LOO-CV driven off the C2
# canonical (limma_voom_qw) DEG method instead of dream.
#
# Mirrors dream_loo_cv_v2.R EXACTLY (leakage-safe Fix A: per-fold filterByExpr +
# calcNormFactors), but swaps the dream() call for the vetted LVQW engine with
# the C2 production design:
#     voomWithQualityWeights -> lmFit -> eBayes,  ~ dataset + inferred_sex + group_binary
# (= 05h_limma_voom_qw_canonical.R, the producer of canonical_deg_results.csv).
#
# Per fold: hold out one mega cohort, re-fit LVQW on the remaining 4 cohorts,
# emit the training-fold DEG table. The aggregator (aggregate_lvqw_loo_cv_C2.R)
# then recovers each fold's training DEGs against the FULL C2 canonical.
#
# Parameterized by HELD_OUT env var.
# Output: results/integration/loo_cv_C2/lvqw_loo_C2_{HELD_OUT}.csv
#         results/integration/loo_cv_C2/loo_C2_metrics_{HELD_OUT}.csv
# ===========================================================================
set.seed(42)
Sys.setenv(R_PARALLEL_SEED = "42")

held_out <- Sys.getenv("HELD_OUT", "")
if (nchar(held_out) == 0) stop("HELD_OUT env var must be set (e.g., GSE126848)")
cat("=== LVQW LOO-CV (C2): holding out", held_out, "===\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(ashr)
  library(yaml)
})

BASE    <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
INT     <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR    <- file.path(INT, "results/integration")
SCRIPTS <- file.path(INT, "scripts")

# --- the single vetted LVQW engine (fit_lvqw) ---
source(file.path(SCRIPTS, "de_engine_lvqw.R"))

loo_dir <- file.path(RDIR, "loo_cv_C2")
dir.create(loo_dir, recursive = TRUE, showWarnings = FALSE)

# ===========================================================================
# 1. Cohort validation (yaml include_in_mega)
# ===========================================================================
ycfg <- yaml::read_yaml(file.path(BASE, "config/human_datasets.yaml"))$datasets
mega_cohorts <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))
if (!held_out %in% mega_cohorts)
  stop(held_out, " not in mega set. Choose: ", paste(mega_cohorts, collapse = ", "))
training_cohorts <- setdiff(mega_cohorts, held_out)
cat("Mega cohorts (k =", length(mega_cohorts), "):", paste(mega_cohorts, collapse = ", "), "\n")
cat("Training (k =", length(training_cohorts), "):", paste(training_cohorts, collapse = ", "), "\n")
cat("Held out:", held_out, "\n\n")

# ===========================================================================
# 2. Load raw counts + metadata, QC filter (pass_technical) — as 03_integrate
# ===========================================================================
cat("Loading merged counts (raw) + metadata...\n")
merged_raw <- readRDS(file.path(RDIR, "merged_counts_raw.rds"))
meta_all   <- readRDS(file.path(RDIR, "meta_matched.rds"))
setDT(meta_all)
qc <- fread(file.path(INT, "qc/sample_qc_report.csv"))
pass_ids <- qc[pass_technical == TRUE, sample_id]
merged_raw <- merged_raw[, colnames(merged_raw) %in% pass_ids]
meta_all   <- meta_all[sample_id %in% colnames(merged_raw)]
meta_all   <- meta_all[match(colnames(merged_raw), sample_id)]
stopifnot(all(meta_all$sample_id == colnames(merged_raw)))
cat("QC-passing samples:", ncol(merged_raw), "\n")

# ===========================================================================
# 3. FIX A: subset to training cohorts, recompute filterByExpr + normFactors
# ===========================================================================
cat("\n-- Fix A: independent gene filtering + normalization for this fold --\n")
train_idx    <- meta_all$dataset %in% training_cohorts
merged_train <- merged_raw[, train_idx]
meta_train   <- meta_all[train_idx]
cat("Training samples:", ncol(merged_train), "\n")

dge_train <- DGEList(counts = merged_train)
dge_train$samples$dataset      <- meta_train$dataset
dge_train$samples$group_binary <- factor(meta_train$group_binary, levels = c("Control", "Disease"))

# filterByExpr WITHIN this fold (same design as 05_dream / 05h producer filter)
design_filter <- model.matrix(~ 0 + group_binary, data = dge_train$samples)
keep <- filterByExpr(dge_train, design = design_filter)
n_genes_pre  <- nrow(dge_train)
n_genes_post <- sum(keep)
cat("Genes pre-filter:", n_genes_pre, "  post-filterByExpr:", n_genes_post, "\n")
dge_train <- dge_train[keep, , keep.lib.sizes = FALSE]

# calcNormFactors WITHIN this fold (RLE, matching 03_integrate_counts.R)
dge_train <- calcNormFactors(dge_train, method = "RLE")
cat("Normalization: RLE (recomputed for this fold)\n")

# ===========================================================================
# 4. Build C2 design info (cohort + sex FIXED effects)
# ===========================================================================
sex <- meta_train$inferred_sex[match(colnames(dge_train), meta_train$sample_id)]
info <- data.frame(
  group_binary = factor(dge_train$samples$group_binary, levels = c("Control", "Disease")),
  dataset      = droplevels(factor(dge_train$samples$dataset)),
  inferred_sex = factor(sex),
  stringsAsFactors = FALSE
)
rownames(info) <- colnames(dge_train)
cat("\nGroup x dataset (training):\n"); print(table(info$group_binary, info$dataset))

if (nlevels(info$dataset) < 2) stop("Need >= 2 datasets; have ", nlevels(info$dataset))

# C2 production design. Drop sex only if degenerate (LOO never triggers this).
design <- if (nlevels(droplevels(info$inferred_sex)) >= 2) {
            model.matrix(~ dataset + inferred_sex + group_binary, data = info)
          } else {
            model.matrix(~ dataset + group_binary, data = info)
          }
coef_name <- "group_binaryDisease"
stopifnot(coef_name %in% colnames(design))
cat("\nDesign:", paste(colnames(design), collapse = " + "), "\n")

# ===========================================================================
# 5. Fit LVQW engine on the training fold (voomWithQualityWeights -> eBayes + ashr)
# ===========================================================================
cat("\nRunning LVQW engine (voomWithQualityWeights -> lmFit -> eBayes -> ashr)...\n")
res_dt <- fit_lvqw(dge_train, design, coef = coef_name, do_ashr = TRUE,
                   mixcompdist = "normal", weights = TRUE)
# add symbol parity column (drop-in with canonical schema) — optional, skip if no metadata
# (kept lean; aggregator joins by gene only.)

train_degs_01    <- res_dt[padj < 0.1, gene]
train_degs_005   <- res_dt[padj < 0.05, gene]
train_degs_lfc   <- res_dt[padj < 0.05 & abs(logFC) >= 0.5, gene]
train_degs_ash   <- res_dt[lfsr < 0.05 & abs(shrunk_logFC) > 0.5, gene]

cat("\n===== TRAINING (LVQW, excluding", held_out, ") =====\n")
cat("Gene universe (fold):", nrow(res_dt), "\n")
cat("DEGs padj<0.1:", length(train_degs_01),
    " | padj<0.05:", length(train_degs_005),
    " | raw Tier1 (padj<.05,|lfc|>.5):", length(train_degs_lfc),
    " | ashr Tier1 (lfsr<.05,|slfc|>.5):", length(train_degs_ash), "\n")

out_file <- file.path(loo_dir, paste0("lvqw_loo_C2_", held_out, ".csv"))
fwrite(res_dt, out_file)
cat("Saved training results:", out_file, "\n")

# ===========================================================================
# 6. Held-out replication (Fix B): training DEGs vs held-out cohort per-study DE
# ===========================================================================
cat("\n-- Fix B: held-out replication --\n")
per_study_dir <- file.path(INT, "results/per_study")
heldout_file  <- file.path(per_study_dir, paste0(held_out, "_de_results.csv"))

metrics <- data.table(
  held_out_cohort = held_out, n_training_cohorts = length(training_cohorts),
  n_samples_training = ncol(dge_train), gene_universe_size = nrow(res_dt),
  n_genes_pre_filter = n_genes_pre,
  n_degs_training_01 = length(train_degs_01), n_degs_training_005 = length(train_degs_005),
  n_degs_training_lfc = length(train_degs_lfc), n_degs_training_ash = length(train_degs_ash),
  n_degs_heldout = NA_integer_, n_genes_common = NA_integer_,
  jaccard_01 = NA_real_, lfc_spearman = NA_real_, lfc_spearman_sig = NA_real_,
  auc_replication = NA_real_, fisher_or = NA_real_, fisher_pval = NA_real_,
  direction_concordance = NA_real_, pi1_heldout = NA_real_
)

if (file.exists(heldout_file)) {
  ho <- fread(heldout_file)
  setnames(ho, "adj.P.Val", "ho_padj", skip_absent = TRUE)
  if (!"ho_padj" %in% names(ho)) setnames(ho, "padj", "ho_padj", skip_absent = TRUE)
  setnames(ho, "logFC", "ho_logFC", skip_absent = TRUE)
  setnames(ho, "P.Value", "ho_pval", skip_absent = TRUE)

  common <- merge(res_dt[, .(gene, train_logFC = logFC, train_padj = padj)],
                  ho[, .(gene, ho_logFC, ho_padj, ho_pval)], by = "gene")
  n_common <- nrow(common)
  cat("Common genes (training & held-out):", n_common, "\n")

  train_sig_common <- common[train_padj < 0.1, gene]
  ho_sig_common    <- common[ho_padj < 0.1, gene]
  jaccard_01 <- if (length(union(train_sig_common, ho_sig_common)) > 0)
    length(intersect(train_sig_common, ho_sig_common)) /
      length(union(train_sig_common, ho_sig_common)) else 0

  lfc_spearman <- cor(common$train_logFC, common$ho_logFC, method = "spearman", use = "pairwise.complete.obs")
  train_sig_df <- common[train_padj < 0.1]
  lfc_spearman_sig <- if (nrow(train_sig_df) > 5)
    cor(train_sig_df$train_logFC, train_sig_df$ho_logFC, method = "spearman", use = "pairwise.complete.obs") else NA_real_
  dir_conc <- if (nrow(train_sig_df) > 0)
    mean(sign(train_sig_df$train_logFC) == sign(train_sig_df$ho_logFC)) * 100 else NA_real_

  ho_sig_binary <- as.integer(common$ho_padj < 0.05)
  train_abs_lfc <- abs(common$train_logFC)
  pos_scores <- train_abs_lfc[ho_sig_binary == 1]; neg_scores <- train_abs_lfc[ho_sig_binary == 0]
  auc_replication <- if (length(pos_scores) > 0 && length(neg_scores) > 0) {
    wt <- suppressWarnings(wilcox.test(pos_scores, neg_scores))
    as.numeric(wt$statistic) / (length(pos_scores) * length(neg_scores))
  } else NA_real_

  ts5 <- common$train_padj < 0.05; hs5 <- common$ho_padj < 0.05
  ct  <- matrix(c(sum(ts5 & hs5), sum(ts5 & !hs5), sum(!ts5 & hs5), sum(!ts5 & !hs5)), nrow = 2)
  ft  <- fisher.test(ct)

  ho_pv <- common[train_padj < 0.1, ho_pval]
  pi1 <- if (length(ho_pv) > 20) {
    lambda <- 0.5; pi0 <- min(1, sum(ho_pv > lambda, na.rm = TRUE) / (length(ho_pv) * (1 - lambda))); 1 - pi0
  } else NA_real_

  metrics[, `:=`(
    n_degs_heldout = length(common[ho_padj < 0.05, gene]), n_genes_common = n_common,
    jaccard_01 = round(jaccard_01, 4), lfc_spearman = round(lfc_spearman, 4),
    lfc_spearman_sig = round(lfc_spearman_sig, 4), auc_replication = round(auc_replication, 4),
    fisher_or = round(as.numeric(ft$estimate), 2), fisher_pval = ft$p.value,
    direction_concordance = round(dir_conc, 1), pi1_heldout = round(pi1, 4)
  )]
  cat(sprintf("Jaccard %.4f | rho %.4f | dir %.1f%% | AUC %.4f | OR %.2f\n",
              jaccard_01, lfc_spearman, dir_conc, auc_replication, ft$estimate))
} else {
  cat("WARNING: held-out per-study DE not found:", heldout_file, "\n")
}

metrics_file <- file.path(loo_dir, paste0("loo_C2_metrics_", held_out, ".csv"))
fwrite(metrics, metrics_file)
cat("\nSaved metrics:", metrics_file, "\n")
cat("\n===== LVQW LOO-CV (C2) COMPLETE (held out:", held_out, ") =====\n")
cat("Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
