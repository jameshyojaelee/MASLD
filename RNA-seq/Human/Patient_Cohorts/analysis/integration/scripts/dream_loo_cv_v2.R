#!/usr/bin/env Rscript
# dream_loo_cv_v2.R
# ---------------------------------------------------------------------------
# Leave-one-cohort-out cross-validation (v2) for dream mega-analysis.
#
# Reviewer-driven improvements over v1 (dream_loo_cv.R):
#   Fix A: Recompute filterByExpr() + calcNormFactors() per fold (no leakage).
#   Fix B: Held-out replication — test training DEGs against per-study DE from
#          the held-out cohort (AUC, Fisher enrichment, Spearman concordance).
#   Fix C: Structured output with per-fold gene-universe size, replication
#          metrics, and training stability metrics.
#
# Parameterized by HELD_OUT env var.
# Output: results/integration/loo_cv_v2/dream_loo_v2_{HELD_OUT}.csv
#                                        loo_v2_metrics_{HELD_OUT}.csv
# ---------------------------------------------------------------------------

# ---- Seed pinning (T2.4, 2026-04-22) -----
set.seed(42)
Sys.setenv(R_PARALLEL_SEED = "42")

held_out <- Sys.getenv("HELD_OUT", "")
if (nchar(held_out) == 0) stop("HELD_OUT env var must be set (e.g., GSE126848)")
cat("=== Dream LOO-CV v2: holding out", held_out, "===\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

suppressPackageStartupMessages({
  library(reformulas)
  library(lme4)
  library(data.table)
  library(edgeR)
  library(yaml)
})

# Force injection into lme4 namespace BEFORE loading variancePartition
ns_lme4 <- asNamespace("lme4")
for (fn in c("findbars", "nobars", "subbars", "rebuildFormula")) {
  if (exists(fn, envir = ns_lme4)) {
    try({
      unlockBinding(fn, ns_lme4)
      assign(fn, get(fn, asNamespace("reformulas")), envir = ns_lme4)
      lockBinding(fn, ns_lme4)
    }, silent = TRUE)
  }
}

suppressPackageStartupMessages({
  library(variancePartition)
  library(BiocParallel)
})

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts"
INT  <- file.path(BASE, "analysis/integration")
RDIR <- file.path(INT, "results/integration")

# --- Output directory (v2) ---
loo_dir <- file.path(RDIR, "loo_cv_v2")
dir.create(loo_dir, recursive = TRUE, showWarnings = FALSE)

# ===========================================================================
# 1. Cohort validation
# ===========================================================================
yaml_path <- file.path(
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
  "config/human_datasets.yaml")
ycfg <- yaml::read_yaml(yaml_path)$datasets
mega_cohorts <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))

if (!held_out %in% mega_cohorts) {
  stop(held_out, " is not in the mega-analysis cohort set (yaml include_in_mega). ",
       "Choose from: ", paste(mega_cohorts, collapse = ", "))
}

training_cohorts <- setdiff(mega_cohorts, held_out)
cat("Mega cohorts (k =", length(mega_cohorts), "):",
    paste(mega_cohorts, collapse = ", "), "\n")
cat("Training cohorts (k =", length(training_cohorts), "):",
    paste(training_cohorts, collapse = ", "), "\n")
cat("Held out:", held_out, "\n\n")

# ===========================================================================
# 2. Load raw counts + metadata (pre-filterByExpr)
# ===========================================================================
cat("Loading merged counts (raw) and metadata...\n")
merged_raw <- readRDS(file.path(RDIR, "merged_counts_raw.rds"))
meta_all   <- readRDS(file.path(RDIR, "meta_matched.rds"))

# Apply QC filter (pass_technical) — same as 03_integrate_counts.R
qc <- fread(file.path(INT, "qc/sample_qc_report.csv"))
pass_ids <- qc[pass_technical == TRUE, sample_id]
merged_raw <- merged_raw[, colnames(merged_raw) %in% pass_ids]
meta_all   <- meta_all[sample_id %in% colnames(merged_raw)]
meta_all   <- meta_all[match(colnames(merged_raw), meta_all$sample_id), ]
stopifnot(all(meta_all$sample_id == colnames(merged_raw)))
cat("QC-passing samples:", ncol(merged_raw), "\n")

# ===========================================================================
# 3. FIX A: Subset to training cohorts, recompute filterByExpr + normFactors
# ===========================================================================
cat("\n── Fix A: Independent gene filtering + normalization for this fold ──\n")

# Subset to training cohorts only
train_idx <- meta_all$dataset %in% training_cohorts
merged_train <- merged_raw[, train_idx]
meta_train   <- meta_all[train_idx, ]

cat("Training samples:", ncol(merged_train), "\n")
cat("Training datasets:", paste(unique(meta_train$dataset), collapse = ", "), "\n")

# Create DGEList from raw counts
dge_train <- DGEList(counts = merged_train)
dge_train$samples$dataset      <- meta_train$dataset
dge_train$samples$group_binary <- factor(meta_train$group_binary, levels = c("Control", "Disease"))

# filterByExpr WITHIN this fold
design_filter <- model.matrix(~ 0 + group_binary, data = dge_train$samples)
keep <- filterByExpr(dge_train, design = design_filter)
n_genes_pre  <- nrow(dge_train)
n_genes_post <- sum(keep)
cat("Genes pre-filter:", n_genes_pre, "\n")
cat("Genes post-filterByExpr:", n_genes_post, "\n")
dge_train <- dge_train[keep, , keep.lib.sizes = FALSE]

# calcNormFactors WITHIN this fold (RLE, matching 03_integrate_counts.R)
dge_train <- calcNormFactors(dge_train, method = "RLE")
cat("Normalization: RLE (recomputed for this fold)\n")

# ===========================================================================
# 4. Build metadata for dream
# ===========================================================================
matched_sex <- meta_train$inferred_sex
if (is.null(matched_sex)) {
  # Fallback: merge from meta_matched
  meta_new <- readRDS(file.path(RDIR, "meta_matched.rds"))
  matched_sex <- meta_new$inferred_sex[match(colnames(dge_train), meta_new$sample_id)]
}
n_na_sex <- sum(is.na(matched_sex))
if (n_na_sex > 0) {
  warning(n_na_sex, " samples have NA inferred_sex; these will be dropped by dream().")
}

info <- data.frame(
  group_binary = factor(dge_train$samples$group_binary, levels = c("Control", "Disease")),
  dataset = droplevels(factor(dge_train$samples$dataset)),
  inferred_sex = factor(matched_sex),
  stringsAsFactors = FALSE
)
rownames(info) <- colnames(dge_train)

cat("\nGroup distribution (training, excluding", held_out, "):\n")
print(table(info$group_binary, info$dataset))
cat("Sex distribution:\n")
print(table(info$inferred_sex, useNA = "always"))

# Check >= 2 datasets for random effect
n_datasets <- length(unique(info$dataset))
if (n_datasets < 2) stop("Need >= 2 datasets for random effect; only have ", n_datasets)

# ===========================================================================
# 5. Run dream on training fold
# ===========================================================================
form <- ~ group_binary + inferred_sex + (1|dataset)
cat("\nFormula:", deparse(form), "\n")

ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
cat("Using", ncpus, "CPU cores\n")
param <- if (ncpus > 1) MulticoreParam(ncpus) else SerialParam()

cat("Running voomWithDreamWeights...\n")
v <- suppressWarnings(voomWithDreamWeights(dge_train, form, info, BPPARAM = param))

cat("Running dream()...\n")
fit <- suppressWarnings(dream(v, form, info, BPPARAM = param))
# NOTE: do NOT call eBayes() after dream()

res <- topTable(fit, coef = "group_binaryDisease", number = Inf, sort.by = "none")
res$gene <- rownames(res)
res_dt <- as.data.table(res)
setnames(res_dt, "adj.P.Val", "padj")

# Training DEGs at two thresholds
train_degs_01    <- res_dt[padj < 0.1, gene]
train_degs_005   <- res_dt[padj < 0.05, gene]
train_degs_lfc   <- res_dt[padj < 0.05 & abs(logFC) >= 0.5, gene]

cat("\n===== TRAINING RESULTS (excluding", held_out, ") =====\n")
cat("Gene universe (this fold):", nrow(res_dt), "\n")
cat("DEGs (padj < 0.1):", length(train_degs_01), "\n")
cat("DEGs (padj < 0.05):", length(train_degs_005), "\n")
cat("DEGs (padj < 0.05, |logFC| >= 0.5):", length(train_degs_lfc), "\n")

# Save training-fold dream results
out_file <- file.path(loo_dir, paste0("dream_loo_v2_", held_out, ".csv"))
fwrite(res_dt, out_file)
cat("Saved training results:", out_file, "\n")

# ===========================================================================
# 6. FIX B: Held-out replication
# ===========================================================================
cat("\n── Fix B: Held-out replication ──\n")

# Load per-study DE for the held-out cohort (from Script 02)
per_study_dir <- file.path(INT, "results/per_study")
heldout_file <- file.path(per_study_dir, paste0(held_out, "_de_results.csv"))

if (!file.exists(heldout_file)) {
  cat("WARNING: Per-study DE results not found for", held_out, "\n")
  cat("  Expected:", heldout_file, "\n")
  cat("  Skipping held-out replication metrics.\n")
  # Write metrics with NA replication
  metrics <- data.table(
    held_out_cohort      = held_out,
    n_training_cohorts   = length(training_cohorts),
    n_samples_training   = ncol(dge_train),
    gene_universe_size   = nrow(res_dt),
    n_genes_pre_filter   = n_genes_pre,
    n_degs_training_01   = length(train_degs_01),
    n_degs_training_005  = length(train_degs_005),
    n_degs_training_lfc  = length(train_degs_lfc),
    n_degs_heldout       = NA_integer_,
    n_genes_common       = NA_integer_,
    jaccard_01           = NA_real_,
    lfc_spearman         = NA_real_,
    lfc_spearman_sig     = NA_real_,
    auc_replication      = NA_real_,
    fisher_or            = NA_real_,
    fisher_pval          = NA_real_,
    direction_concordance = NA_real_,
    pi1_heldout          = NA_real_
  )
} else {
  ho <- fread(heldout_file)
  cat("Held-out per-study DE:", held_out, "—", nrow(ho), "genes\n")

  # The per-study DE runs on the merged_dge gene universe (Script 02 uses
  # the same DGEList), so column names are the same.
  # Columns: logFC, AveExpr, t, P.Value, adj.P.Val, B, gene, dataset, SE, df.total, sex_source
  setnames(ho, "adj.P.Val", "ho_padj", skip_absent = TRUE)
  if (!"ho_padj" %in% names(ho)) setnames(ho, "padj", "ho_padj", skip_absent = TRUE)
  setnames(ho, "logFC", "ho_logFC", skip_absent = TRUE)
  setnames(ho, "P.Value", "ho_pval", skip_absent = TRUE)

  # Merge training and held-out on common genes
  common <- merge(
    res_dt[, .(gene, train_logFC = logFC, train_padj = padj, train_t = t)],
    ho[, .(gene, ho_logFC, ho_padj, ho_pval)],
    by = "gene"
  )
  n_common <- nrow(common)
  cat("Common genes (training & held-out):", n_common, "\n")

  ho_degs_01  <- common[ho_padj < 0.1, gene]
  ho_degs_005 <- common[ho_padj < 0.05, gene]
  cat("Held-out DEGs (padj < 0.1, among common):", length(ho_degs_01), "\n")
  cat("Held-out DEGs (padj < 0.05, among common):", length(ho_degs_005), "\n")

  # ------ Metric 1: Jaccard (training DEGs at padj<0.1 vs held-out padj<0.1) ------
  train_sig_common <- common[train_padj < 0.1, gene]
  ho_sig_common    <- common[ho_padj < 0.1, gene]
  jaccard_01 <- if (length(union(train_sig_common, ho_sig_common)) > 0) {
    length(intersect(train_sig_common, ho_sig_common)) /
      length(union(train_sig_common, ho_sig_common))
  } else 0
  cat(sprintf("Jaccard (train vs held-out, padj<0.1): %.4f\n", jaccard_01))

  # ------ Metric 2: Effect-size concordance (Spearman) ------
  # (a) All common genes
  lfc_spearman <- cor(common$train_logFC, common$ho_logFC, method = "spearman",
                      use = "pairwise.complete.obs")
  cat(sprintf("LFC Spearman (all common): %.4f\n", lfc_spearman))

  # (b) Restricted to training DEGs (padj<0.1)
  train_sig_df <- common[train_padj < 0.1]
  lfc_spearman_sig <- if (nrow(train_sig_df) > 5) {
    cor(train_sig_df$train_logFC, train_sig_df$ho_logFC, method = "spearman",
        use = "pairwise.complete.obs")
  } else NA_real_
  cat(sprintf("LFC Spearman (training DEGs only): %.4f\n", lfc_spearman_sig))

  # ------ Metric 3: Direction concordance ------
  # Among genes significant in training (padj<0.1), what fraction have same sign in held-out?
  dir_conc <- if (nrow(train_sig_df) > 0) {
    mean(sign(train_sig_df$train_logFC) == sign(train_sig_df$ho_logFC)) * 100
  } else NA_real_
  cat(sprintf("Direction concordance (training DEGs): %.1f%%\n", dir_conc))

  # ------ Metric 4: AUC — training |logFC| predicts held-out significance ------
  # Binary outcome: held-out padj < 0.05
  # Predictor: absolute training logFC (higher => more likely to replicate)
  ho_sig_binary <- as.integer(common$ho_padj < 0.05)
  train_abs_lfc <- abs(common$train_logFC)

  # Wilcoxon-Mann-Whitney AUC (concordance statistic) — no extra dependencies
  # AUC = P(score_positive > score_negative) estimated by the U-statistic
  pos_scores <- train_abs_lfc[ho_sig_binary == 1]
  neg_scores <- train_abs_lfc[ho_sig_binary == 0]

  auc_replication <- if (length(pos_scores) > 0 && length(neg_scores) > 0) {
    wt <- suppressWarnings(wilcox.test(pos_scores, neg_scores))
    wt$statistic / (length(pos_scores) * length(neg_scores))
  } else NA_real_
  cat(sprintf("AUC (train |logFC| -> held-out padj<0.05): %.4f\n", auc_replication))

  # ------ Metric 5: Fisher's exact test ------
  # Training DEGs (padj<0.05) enriched among held-out DEGs (padj<0.05)?
  train_sig_005 <- common$train_padj < 0.05
  ho_sig_005    <- common$ho_padj < 0.05
  ct <- matrix(c(
    sum( train_sig_005 &  ho_sig_005),  # both sig
    sum( train_sig_005 & !ho_sig_005),  # train only
    sum(!train_sig_005 &  ho_sig_005),  # heldout only
    sum(!train_sig_005 & !ho_sig_005)   # neither
  ), nrow = 2)
  ft <- fisher.test(ct)
  cat(sprintf("Fisher OR: %.2f (p = %.2e)\n", ft$estimate, ft$p.value))

  # ------ Metric 6: pi1 (Storey-style) ------
  # Estimate the fraction of true positives among held-out p-values for training DEGs
  # pi1 = 1 - pi0, where pi0 = #{p > lambda} / (n * (1 - lambda))
  ho_pvals_traindegs <- common[train_padj < 0.1, ho_pval]
  pi1 <- if (length(ho_pvals_traindegs) > 20) {
    lambda <- 0.5
    pi0 <- min(1, sum(ho_pvals_traindegs > lambda, na.rm = TRUE) /
                 (length(ho_pvals_traindegs) * (1 - lambda)))
    1 - pi0
  } else NA_real_
  cat(sprintf("pi1 (held-out, among training DEGs): %.4f\n", pi1))

  # ------ Assemble metrics ------
  metrics <- data.table(
    held_out_cohort      = held_out,
    n_training_cohorts   = length(training_cohorts),
    n_samples_training   = ncol(dge_train),
    gene_universe_size   = nrow(res_dt),
    n_genes_pre_filter   = n_genes_pre,
    n_degs_training_01   = length(train_degs_01),
    n_degs_training_005  = length(train_degs_005),
    n_degs_training_lfc  = length(train_degs_lfc),
    n_degs_heldout       = length(ho_degs_005),
    n_genes_common       = n_common,
    jaccard_01           = round(jaccard_01, 4),
    lfc_spearman         = round(lfc_spearman, 4),
    lfc_spearman_sig     = round(lfc_spearman_sig, 4),
    auc_replication      = round(auc_replication, 4),
    fisher_or            = round(ft$estimate, 2),
    fisher_pval          = ft$p.value,
    direction_concordance = round(dir_conc, 1),
    pi1_heldout          = round(pi1, 4)
  )
}

# Save per-fold metrics
metrics_file <- file.path(loo_dir, paste0("loo_v2_metrics_", held_out, ".csv"))
fwrite(metrics, metrics_file)
cat("\nSaved metrics:", metrics_file, "\n")

cat("\n===== LOO-CV v2 COMPLETE (held out:", held_out, ") =====\n")
cat("Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
