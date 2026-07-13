#!/usr/bin/env Rscript
# figS_steatosis_bmi_sensitivity.R
# BMI Sensitivity Analysis using Steatosis Grade as BMI Proxy
#
# Context: BMI data is unavailable across all 10 MASLD cohorts. Steatosis grade
# (0-3) is the histological measure most correlated with BMI (r~0.5-0.7 in
# published data). This script adds steatosis to the dream model as a covariate
# and compares results to the unadjusted model on the SAME samples.
#
# Data availability audit (Script 189 confirmed):
#   - Steatosis grade: 76 samples, GSE130970 only (single dataset, no random effect)
#   - NAS score: 660 samples, 5 datasets (supports random intercept for dataset)
#
# Strategy: Two complementary analyses
#   A) Steatosis-adjusted dream (N=76, GSE130970 only, fixed-effect model)
#   B) NAS-adjusted dream (N=660, 5 datasets, mixed-effect model)
#      NAS is used because steatosis is its largest component (0-3 of 0-8)
#      and NAS correlates with BMI (r~0.4-0.6 in published data)
#
# CRITICAL: Both formulas run on the SAME sample subset to isolate the covariate
# effect from sample size effects.
#
# Follows the pattern of Script 14.3 (age sensitivity: rho=0.998, Jaccard=0.938)
#
# Usage: sbatch run_steatosis_sensitivity.sbatch
# SLURM: cpu partition, 16 CPUs, 120GB RAM, 48h
# Env: micromamba activate rnaseq

suppressPackageStartupMessages({
  library(reformulas)
  library(lme4)
  library(data.table)
  library(edgeR)
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
  library(ggplot2)
  library(patchwork)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR <- file.path(INT, "results/integration")
STAGING_DIR <- file.path(INT, "results/staging_classifier")
OUTDIR <- file.path(BASE, "RNA-seq/results/audit_sensitivity/steatosis_sensitivity")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
# Canonical limma-voom-QW (C2) engine helper: build_design_guarded()
source(file.path(INT, "scripts/de_engine_lvqw.R"))

# Two DEG gates are used here, by design:
#  (1) PRIMARY canonical reference = TREAT (McCarthy & Smyth 2009): treat_fdr<0.05
#      at lfc=0.25 — the paper-wide canonical (canonical_deg_results.csv), cited as
#      the headline DEG count.
#  (2) SUBSET-INTERNAL significance gate = padj<0.05 (no |logFC| floor), used for
#      the unadjusted-vs-adjusted OVERLAP/Jaccard panels. TREAT@0.25 is underpowered
#      on these small covariate subsets (N=76 steatosis / N=407 NAS) and collapses
#      the overlap sets to 0; a plain significance gate keeps them interpretable.
#      The robustness conclusion rests on the threshold-free logFC concordance (rho).
TREAT_FDR_T <- 0.05   # (1) primary canonical TREAT gate
SIG_P       <- 0.05   # (2) subset-internal overlap significance gate (padj only)

FIGDIR <- FIGS_SENS_DIR
dir.create(file.path(FIGDIR, "panels"), showWarnings = FALSE, recursive = TRUE)

cat("=== Steatosis/BMI Sensitivity Analysis ===\n")
cat("Started:", as.character(Sys.time()), "\n\n")

# --- Setup parallel backend ---
ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
cat("Using", ncpus, "CPU cores\n")
param <- if (ncpus > 1) MulticoreParam(ncpus) else SerialParam()

# ============================================================
# DATA LOADING
# ============================================================
cat("Loading merged DGE...\n")
dge <- readRDS(file.path(RDIR, "merged_dge.rds"))
cat("  DGE dimensions:", nrow(dge), "genes x", ncol(dge), "samples\n")

cat("Loading matched metadata...\n")
meta_matched <- readRDS(file.path(RDIR, "meta_matched.rds"))

cat("Loading staging metadata (has steatosis_grade)...\n")
meta_staging <- fread(file.path(STAGING_DIR, "modeling_metadata.csv"))

cat("Loading QC report...\n")
qc <- fread(file.path(INT, "qc/sample_qc_report.csv"))
pass_samples <- qc[pass_technical == TRUE, sample_id]

cat("Loading primary dream results...\n")
dream_primary <- fread(file.path(RDIR, "canonical_deg_results.csv"))
if ("adj.P.Val" %in% names(dream_primary) && !"padj" %in% names(dream_primary))
  setnames(dream_primary, "adj.P.Val", "padj")
stopifnot("treat_fdr" %in% names(dream_primary))  # canonical TREAT gate
cat("  Primary canonical DEGs (TREAT FDR<0.05 @ lfc=0.25):",
    sum(dream_primary$treat_fdr < TREAT_FDR_T, na.rm = TRUE), "\n\n")

# ============================================================
# STEATOSIS DATA AUDIT
# ============================================================
cat("=== Steatosis Data Audit ===\n")

# Steatosis grade: -1 = missing in modeling_metadata
steat_samples <- meta_staging[steatosis_grade != -1, sample_id]
steat_datasets <- meta_staging[steatosis_grade != -1, unique(dataset)]
cat("Steatosis grade available:", length(steat_samples), "samples from",
    paste(steat_datasets, collapse = ", "), "\n")
cat("Steatosis distribution:\n")
print(meta_staging[steatosis_grade != -1, .N, by = steatosis_grade][order(steatosis_grade)])

# NAS score: available for 5 datasets
nas_samples <- meta_staging[!is.na(nas_score) & nas_score >= 0, sample_id]
nas_datasets <- meta_staging[!is.na(nas_score) & nas_score >= 0, unique(dataset)]
cat("\nNAS score available:", length(nas_samples), "samples from",
    paste(nas_datasets, collapse = ", "), "\n")
cat("NAS distribution:\n")
print(meta_staging[!is.na(nas_score) & nas_score >= 0, .N, by = nas_score][order(nas_score)])

# Also check which of those have controls
meta_staging_qc <- meta_staging[sample_id %in% pass_samples]
nas_with_groups <- meta_staging_qc[!is.na(nas_score) & nas_score >= 0]
cat("\nNAS samples by dataset and group_binary:\n")
print(nas_with_groups[, .N, by = .(dataset, group_binary)][order(dataset, group_binary)])


# ============================================================
# HELPER: run_dream_comparison
# ============================================================
run_dream_comparison <- function(dge, meta_matched, sample_ids, covar_name, covar_values,
                                  dataset_col, label, param) {
  cat(sprintf("\n=== %s: Dream Comparison (N=%d) ===\n", label, length(sample_ids)))

  # Intersect with DGE columns
  valid_samples <- intersect(sample_ids, colnames(dge))
  cat("Valid samples in DGE:", length(valid_samples), "\n")

  # Subset DGE
  dge_sub <- dge[, valid_samples]

  # Build info
  meta_idx <- match(colnames(dge_sub), meta_matched$sample_id)
  sex_sub <- meta_matched$inferred_sex[meta_idx]

  # Get dataset and covariate
  staging_idx <- match(colnames(dge_sub), meta_staging$sample_id)
  ds_sub <- meta_staging$dataset[staging_idx]
  covar_sub <- covar_values[match(colnames(dge_sub), sample_ids)]

  info_sub <- data.frame(
    group_binary = factor(
      meta_matched$group_binary[meta_idx],
      levels = c("Control", "Disease")
    ),
    dataset      = factor(ds_sub),
    inferred_sex = factor(sex_sub),
    covar        = as.numeric(covar_sub),
    row.names    = colnames(dge_sub),
    stringsAsFactors = FALSE
  )

  # Scale covariate
  info_sub$covar_scaled <- scale(info_sub$covar)[, 1]

  # Check dataset levels
  n_datasets <- nlevels(droplevels(info_sub$dataset))
  cat("Datasets in subset:", n_datasets, "\n")

  # Dataset group counts
  ds_groups <- data.table(
    dataset = ds_sub,
    group = meta_matched$group_binary[meta_idx]
  )[, .(n_ctrl = sum(group == "Control"), n_dis = sum(group == "Disease")), by = dataset]
  cat("Dataset group counts:\n")
  print(ds_groups)

  # Check for datasets with only one group (need controls for group_binary effect)
  ds_with_both <- ds_groups[n_ctrl > 0 & n_dis > 0, dataset]
  if (length(ds_with_both) < n_datasets) {
    cat("WARNING: Some datasets lack controls. Filtering to datasets with both groups.\n")
    keep_idx <- ds_sub %in% ds_with_both
    dge_sub <- dge_sub[, keep_idx]
    info_sub <- info_sub[keep_idx, ]
    info_sub$dataset <- droplevels(info_sub$dataset)
    n_datasets <- nlevels(info_sub$dataset)
    cat("  After filter:", ncol(dge_sub), "samples,", n_datasets, "datasets\n")
  }

  # Filter by expression
  keep_expr <- filterByExpr(dge_sub, group = info_sub$group_binary)
  dge_sub <- dge_sub[keep_expr, , keep.lib.sizes = FALSE]
  dge_sub <- calcNormFactors(dge_sub, method = "RLE")
  cat("After expression filter:", nrow(dge_sub), "genes\n")

  # --- limma-voom QW Run 1: WITHOUT covariate (dataset FIXED when >1) ---
  cat("\n--- limma-voom QW: WITHOUT", covar_name, "---\n")
  base_terms <- if (n_datasets > 1) {
    c("dataset", "group_binary", "inferred_sex")
  } else {
    c("group_binary", "inferred_sex")
  }
  des_base <- build_design_guarded(info_sub, base_terms)
  stopifnot("group_binaryDisease" %in% colnames(des_base$design))
  cat("Design terms:", paste(base_terms, collapse = " + "), "\n")

  v_base <- limma::voomWithQualityWeights(dge_sub, des_base$design)
  fit_base <- limma::eBayes(limma::lmFit(v_base, des_base$design))
  res_base <- topTable(fit_base, coef = "group_binaryDisease", number = Inf, sort.by = "none")
  res_base$gene <- rownames(res_base)
  res_base_dt <- as.data.table(res_base)
  if ("adj.P.Val" %in% names(res_base_dt)) setnames(res_base_dt, "adj.P.Val", "padj")

  n_deg_base <- sum(res_base_dt$padj < SIG_P, na.rm = TRUE)
  cat("DEGs (padj<0.05 sig. gate) without", covar_name, ":", n_deg_base, "\n")

  # --- limma-voom QW Run 2: WITH covariate (dataset FIXED when >1) ---
  cat("\n--- limma-voom QW: WITH", covar_name, "---\n")
  adj_terms <- c(base_terms, "covar_scaled")
  des_adj <- build_design_guarded(info_sub, adj_terms)
  stopifnot("group_binaryDisease" %in% colnames(des_adj$design))
  stopifnot("covar_scaled" %in% colnames(des_adj$design))
  cat("Design terms:", paste(adj_terms, collapse = " + "), "\n")

  v_adj <- limma::voomWithQualityWeights(dge_sub, des_adj$design)
  fit_adj <- limma::eBayes(limma::lmFit(v_adj, des_adj$design))
  res_adj <- topTable(fit_adj, coef = "group_binaryDisease", number = Inf, sort.by = "none")
  res_adj$gene <- rownames(res_adj)
  res_adj_dt <- as.data.table(res_adj)
  if ("adj.P.Val" %in% names(res_adj_dt)) setnames(res_adj_dt, "adj.P.Val", "padj")

  n_deg_adj <- sum(res_adj_dt$padj < SIG_P, na.rm = TRUE)
  cat("DEGs (padj<0.05 sig. gate) with", covar_name, ":", n_deg_adj, "\n")

  # Also extract covariate coefficient
  res_covar <- tryCatch({
    tt <- topTable(fit_adj, coef = "covar_scaled", number = Inf, sort.by = "none")
    tt$gene <- rownames(tt)
    as.data.table(tt)
  }, error = function(e) {
    cat("  Could not extract covariate coefficient:", conditionMessage(e), "\n")
    NULL
  })
  if (!is.null(res_covar)) {
    if ("adj.P.Val" %in% names(res_covar)) setnames(res_covar, "adj.P.Val", "padj")
    n_covar_sig <- sum(res_covar$padj < 0.05, na.rm = TRUE)
    cat("Genes significantly associated with", covar_name, "(padj<0.05):", n_covar_sig, "\n")
  }

  # --- Comparison metrics ---
  cat("\n--- Concordance Metrics ---\n")
  common_genes <- intersect(res_base_dt$gene, res_adj_dt$gene)
  base_aligned <- res_base_dt[match(common_genes, gene)]
  adj_aligned  <- res_adj_dt[match(common_genes, gene)]

  rho_lfc <- cor(base_aligned$logFC, adj_aligned$logFC, method = "spearman", use = "complete.obs")
  r_lfc   <- cor(base_aligned$logFC, adj_aligned$logFC, method = "pearson", use = "complete.obs")
  rho_t   <- cor(base_aligned$t, adj_aligned$t, method = "spearman", use = "complete.obs")

  dir_concord <- mean(sign(base_aligned$logFC) == sign(adj_aligned$logFC), na.rm = TRUE) * 100

  degs_base <- base_aligned[padj < SIG_P, gene]
  degs_adj  <- adj_aligned[padj < SIG_P, gene]
  jaccard <- if (length(union(degs_base, degs_adj)) > 0) {
    length(intersect(degs_base, degs_adj)) / length(union(degs_base, degs_adj))
  } else 0

  gained <- setdiff(degs_adj, degs_base)
  lost   <- setdiff(degs_base, degs_adj)

  lfc_shift <- adj_aligned$logFC - base_aligned$logFC

  cat("Spearman rho (logFC):", round(rho_lfc, 4), "\n")
  cat("Pearson r (logFC):", round(r_lfc, 4), "\n")
  cat("Spearman rho (t-stat):", round(rho_t, 4), "\n")
  cat("Direction concordance:", round(dir_concord, 1), "%\n")
  cat("Jaccard (padj<0.05 sig. gate):", round(jaccard, 4), "\n")
  cat("DEGs unadjusted:", length(degs_base), " | adjusted:", length(degs_adj), "\n")
  cat("Gained:", length(gained), " | Lost:", length(lost), "\n")
  cat("Mean |logFC shift|:", round(mean(abs(lfc_shift), na.rm = TRUE), 4), "\n")

  # Per-gene comparison table
  comparison_dt <- data.table(
    gene = common_genes,
    logFC_unadjusted = base_aligned$logFC,
    logFC_adjusted = adj_aligned$logFC,
    logFC_shift = lfc_shift,
    padj_unadjusted = base_aligned$padj,
    padj_adjusted = adj_aligned$padj,
    t_unadjusted = base_aligned$t,
    t_adjusted = adj_aligned$t,
    sig_unadjusted = base_aligned$padj < SIG_P,
    sig_adjusted = adj_aligned$padj < SIG_P,
    status = ifelse(
      base_aligned$padj < SIG_P &
      adj_aligned$padj < SIG_P, "shared",
      ifelse(base_aligned$padj < SIG_P, "lost",
      ifelse(adj_aligned$padj < SIG_P, "gained", "ns_both")))
  )
  comparison_dt <- comparison_dt[order(abs(logFC_shift), decreasing = TRUE)]

  # Metrics table
  metrics_dt <- data.table(
    analysis = label,
    metric = c("n_samples", "n_datasets", "n_genes_tested",
               "n_deg_unadjusted", "n_deg_adjusted",
               "spearman_rho_logFC", "pearson_r_logFC", "spearman_rho_tstat",
               "direction_concordance_pct", "jaccard_padj005",
               "n_gained", "n_lost",
               "mean_abs_lfc_shift", "median_abs_lfc_shift", "max_abs_lfc_shift",
               "n_covar_sig_genes"),
    value = c(ncol(dge_sub), n_datasets, length(common_genes),
              length(degs_base), length(degs_adj),
              round(rho_lfc, 4), round(r_lfc, 4), round(rho_t, 4),
              round(dir_concord, 1), round(jaccard, 4),
              length(gained), length(lost),
              round(mean(abs(lfc_shift), na.rm = TRUE), 4),
              round(median(abs(lfc_shift), na.rm = TRUE), 4),
              round(max(abs(lfc_shift), na.rm = TRUE), 4),
              if (!is.null(res_covar)) sum(res_covar$padj < 0.05, na.rm = TRUE) else NA)
  )

  list(
    comparison = comparison_dt,
    metrics = metrics_dt,
    res_base = res_base_dt,
    res_adj = res_adj_dt,
    res_covar = res_covar,
    info = info_sub,
    n_deg_base = n_deg_base,
    n_deg_adj = n_deg_adj
  )
}


# ============================================================
# ANALYSIS A: Steatosis Grade (N=76, GSE130970 only)
# ============================================================
steat_meta <- meta_staging[steatosis_grade != -1 & sample_id %in% pass_samples]
steat_ids <- steat_meta$sample_id
steat_vals <- steat_meta$steatosis_grade

results_steat <- run_dream_comparison(
  dge = dge, meta_matched = meta_matched,
  sample_ids = steat_ids, covar_name = "steatosis_grade",
  covar_values = steat_vals, dataset_col = steat_meta$dataset,
  label = "steatosis_grade_N76", param = param
)


# ============================================================
# ANALYSIS B: NAS Score (N=660, 5 datasets)
# ============================================================
# NAS includes steatosis (0-3), lobular inflammation (0-3), ballooning (0-2)
# NAS correlates with BMI (r~0.4-0.6); steatosis is the largest component
nas_meta <- meta_staging_qc[!is.na(nas_score) & nas_score >= 0]
nas_ids <- nas_meta$sample_id
nas_vals <- nas_meta$nas_score

results_nas <- run_dream_comparison(
  dge = dge, meta_matched = meta_matched,
  sample_ids = nas_ids, covar_name = "nas_score",
  covar_values = nas_vals, dataset_col = nas_meta$dataset,
  label = "nas_score_N660", param = param
)


# ============================================================
# SAVE RESULTS
# ============================================================
cat("\n=== Saving Results ===\n")

# Combined comparison
combined_comparison <- rbind(
  results_steat$comparison[, analysis := "steatosis_grade"],
  results_nas$comparison[, analysis := "nas_score"]
)
fwrite(combined_comparison,
  file.path(OUTDIR, "steatosis_sensitivity_comparison.csv"))

# Combined metrics
combined_metrics <- rbind(results_steat$metrics, results_nas$metrics)
fwrite(combined_metrics, file.path(OUTDIR, "steatosis_sensitivity_metrics.csv"))

# Individual dream results
fwrite(results_steat$res_base, file.path(OUTDIR, "dream_steatosis_unadjusted.csv"))
fwrite(results_steat$res_adj,  file.path(OUTDIR, "dream_steatosis_adjusted.csv"))
fwrite(results_nas$res_base,   file.path(OUTDIR, "dream_nas_unadjusted.csv"))
fwrite(results_nas$res_adj,    file.path(OUTDIR, "dream_nas_adjusted.csv"))

# Covariate coefficients
if (!is.null(results_steat$res_covar)) {
  fwrite(results_steat$res_covar, file.path(OUTDIR, "steatosis_covariate_effects.csv"))
}
if (!is.null(results_nas$res_covar)) {
  fwrite(results_nas$res_covar, file.path(OUTDIR, "nas_covariate_effects.csv"))
}

# Top shifted genes
top_shift_steat <- results_steat$comparison[order(abs(logFC_shift), decreasing = TRUE)][1:min(50, .N)]
top_shift_nas   <- results_nas$comparison[order(abs(logFC_shift), decreasing = TRUE)][1:min(50, .N)]
fwrite(top_shift_steat, file.path(OUTDIR, "top_shifted_genes_steatosis.csv"))
fwrite(top_shift_nas,   file.path(OUTDIR, "top_shifted_genes_nas.csv"))

# Full results CSV (for paper reference)
full_results <- data.table(
  gene = results_nas$comparison$gene,
  logFC_unadjusted = results_nas$comparison$logFC_unadjusted,
  logFC_nas_adjusted = results_nas$comparison$logFC_adjusted,
  logFC_shift_nas = results_nas$comparison$logFC_shift,
  padj_unadjusted = results_nas$comparison$padj_unadjusted,
  padj_nas_adjusted = results_nas$comparison$padj_adjusted,
  sig_unadjusted = results_nas$comparison$sig_unadjusted,
  sig_nas_adjusted = results_nas$comparison$sig_adjusted,
  status_nas = results_nas$comparison$status
)
fwrite(full_results, file.path(BASE, "RNA-seq/results/audit_sensitivity/steatosis_sensitivity_results.csv"))
fwrite(combined_comparison, file.path(BASE, "RNA-seq/results/audit_sensitivity/steatosis_sensitivity_comparison.csv"))


# ============================================================
# FIGURES
# ============================================================
cat("\n=== Generating Figures ===\n")

# --- Panel A: Sample availability ---
# Bar chart: N samples with steatosis per dataset + steatosis grade distribution
avail_dt <- meta_staging[, .(
  n_steatosis = sum(steatosis_grade != -1),
  n_nas = sum(!is.na(nas_score) & nas_score >= 0),
  n_total = .N
), by = dataset][order(-n_nas)]

avail_long <- melt(avail_dt, id.vars = "dataset",
  measure.vars = c("n_steatosis", "n_nas"),
  variable.name = "covariate", value.name = "n")
avail_long[, covariate := factor(covariate,
  levels = c("n_steatosis", "n_nas"),
  labels = c("Steatosis grade", "NAS score"))]
avail_long[, dataset := factor(dataset, levels = avail_dt$dataset)]

p_avail <- ggplot(avail_long, aes(x = dataset, y = n, fill = covariate)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = c(masld_colors$up, masld_colors$down)) +
  labs(x = NULL, y = "Samples with data", fill = NULL,
       title = "BMI proxy availability by cohort") +
  theme_masld(base_size = 7) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        legend.position = "top")

# Steatosis grade distribution (inset-like)
steat_dist <- meta_staging[steatosis_grade != -1, .N, by = steatosis_grade][order(steatosis_grade)]
steat_dist[, steatosis_grade := factor(steatosis_grade)]

p_steat_dist <- ggplot(steat_dist, aes(x = steatosis_grade, y = N)) +
  geom_col(fill = masld_colors$up, width = 0.6) +
  labs(x = "Steatosis grade", y = "N samples",
       title = "GSE130970 steatosis (N=76)") +
  theme_masld(base_size = 7)


# --- Panel B: logFC scatter (NAS-adjusted, primary analysis) ---
nas_comp <- results_nas$comparison
nas_comp[, sig_cat := fifelse(
  sig_unadjusted & sig_adjusted, "Both significant",
  fifelse(sig_unadjusted & !sig_adjusted, "Lost after adjustment",
  fifelse(!sig_unadjusted & sig_adjusted, "Gained after adjustment",
  "Neither significant")))]

rho_nas <- combined_metrics[analysis == "nas_score_N660" & metric == "spearman_rho_logFC", as.numeric(value)]
r_nas   <- combined_metrics[analysis == "nas_score_N660" & metric == "pearson_r_logFC", as.numeric(value)]
dir_nas <- combined_metrics[analysis == "nas_score_N660" & metric == "direction_concordance_pct", as.numeric(value)]

p_scatter_nas <- ggplot(nas_comp, aes(x = logFC_unadjusted, y = logFC_adjusted, color = sig_cat)) +
  geom_point(size = 0.3, alpha = 0.4) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = c(
    "Both significant" = masld_colors$up,
    "Lost after adjustment" = "#FF9800",
    "Gained after adjustment" = "#4CAF50",
    "Neither significant" = masld_colors$ns
  )) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, size = 2.2,
    label = sprintf("rho = %.3f\nr = %.3f\nDir. = %.1f%%", rho_nas, r_nas, dir_nas)) +
  labs(x = "logFC (unadjusted)", y = "logFC (NAS-adjusted)",
       title = sprintf("NAS-adjusted limma-voom QW (N=%s, %s datasets)",
         combined_metrics[analysis == "nas_score_N660" & metric == "n_samples", value],
         combined_metrics[analysis == "nas_score_N660" & metric == "n_datasets", value]),
       color = NULL) +
  theme_masld(base_size = 7) +
  theme(legend.position = "bottom", legend.key.size = unit(3, "mm"),
        legend.text = element_text(size = 5)) +
  guides(color = guide_legend(override.aes = list(size = 1.5)))


# --- Panel C: DEG overlap bar chart ---
status_counts_nas <- nas_comp[, .N, by = status]
status_counts_nas[, status := factor(status, levels = c("shared", "lost", "gained", "ns_both"))]

jac_nas <- combined_metrics[analysis == "nas_score_N660" & metric == "jaccard_padj005", as.numeric(value)]

p_overlap <- ggplot(status_counts_nas[status != "ns_both"], aes(x = status, y = N, fill = status)) +
  geom_col(width = 0.6) +
  scale_fill_manual(values = c(
    shared = masld_colors$up,
    lost = "#FF9800",
    gained = "#4CAF50"
  )) +
  labs(x = NULL, y = "Number of DEGs",
       title = sprintf("DEG overlap (Jaccard = %.3f)", jac_nas)) +
  theme_masld(base_size = 7) +
  theme(legend.position = "none")


# --- Panel D: Top 20 genes with largest |logFC shift| (NAS) ---
top20 <- nas_comp[order(abs(logFC_shift), decreasing = TRUE)][1:min(20, .N)]

# Try to get gene symbols
gene_map <- dream_primary[, .(gene, gene_symbol = gene)]  # fallback
if ("gene_symbol" %in% names(dream_primary)) {
  gene_map <- dream_primary[, .(gene, gene_symbol)]
}
top20 <- merge(top20, gene_map, by = "gene", all.x = TRUE)
top20[is.na(gene_symbol), gene_symbol := sub("\\..*", "", gene)]
top20[, gene_label := gene_symbol]
top20 <- top20[order(abs(logFC_shift))]
top20[, gene_label := factor(gene_label, levels = gene_label)]

p_shifts <- ggplot(top20, aes(x = logFC_shift, y = gene_label,
                               fill = ifelse(logFC_shift > 0, "Increased", "Decreased"))) +
  geom_col() +
  scale_fill_manual(values = c(Increased = masld_colors$up, Decreased = masld_colors$down)) +
  labs(x = "logFC shift (NAS-adjusted - unadjusted)",
       y = NULL, fill = NULL,
       title = "Top 20 genes: largest logFC shift") +
  theme_masld(base_size = 7) +
  theme(legend.position = "bottom", legend.key.size = unit(3, "mm"))


# --- Panel E: NAS/steatosis coefficient histogram ---
if (!is.null(results_nas$res_covar)) {
  covar_dt <- results_nas$res_covar
  covar_dt[, sig := padj < 0.05]
  n_sig_covar <- sum(covar_dt$sig, na.rm = TRUE)

  p_covar <- ggplot(covar_dt, aes(x = logFC, fill = sig)) +
    geom_histogram(bins = 80, alpha = 0.8) +
    scale_fill_manual(values = c("FALSE" = masld_colors$ns, "TRUE" = masld_colors$up),
                      labels = c("NS", "padj < 0.05")) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    labs(x = "NAS coefficient (logFC per scaled unit)",
         y = "Number of genes", fill = NULL,
         title = sprintf("NAS covariate effect (%d sig. genes)", n_sig_covar)) +
    theme_masld(base_size = 7) +
    theme(legend.position = "top")
} else {
  p_covar <- ggplot() + theme_void() + labs(title = "NAS covariate: not available")
}


# --- Compose figure ---
layout <- "
AABB
CCDD
EEEE
"
p_combined <- (p_avail + p_steat_dist) / (p_scatter_nas + p_overlap) / (p_shifts + p_covar) +
  plot_annotation(
    title = "BMI Sensitivity Analysis: Steatosis/NAS as BMI Proxy",
    subtitle = sprintf(
      "NAS-adjusted (N=%s, %s datasets): rho=%.3f, Jaccard=%.3f | Steatosis (N=76, 1 dataset): see metrics",
      combined_metrics[analysis == "nas_score_N660" & metric == "n_samples", value],
      combined_metrics[analysis == "nas_score_N660" & metric == "n_datasets", value],
      rho_nas, jac_nas
    ),
    theme = theme(
      plot.title = element_text(size = 9, face = "plain"),
      plot.subtitle = element_text(size = 6)
    )
  ) +
  plot_layout(heights = c(1, 1, 1))

ggsave(file.path(FIGDIR, "figS_steatosis_bmi_sensitivity.pdf"),
  p_combined, width = 7.5, height = 9, device = cairo_pdf)
cat("Figure saved to", file.path(FIGDIR, "figS_steatosis_bmi_sensitivity.pdf"), "\n")


# ============================================================
# SUMMARY
# ============================================================
cat("\n")
cat("============================================================\n")
cat("SUMMARY: Steatosis/BMI Sensitivity Analysis\n")
cat("============================================================\n")

cat("\nA) Steatosis Grade (N=76, GSE130970, 1 dataset):\n")
sm <- results_steat$metrics
cat("   Spearman rho (logFC):", sm[metric == "spearman_rho_logFC", value], "\n")
cat("   Pearson r (logFC):", sm[metric == "pearson_r_logFC", value], "\n")
cat("   Direction concordance:", sm[metric == "direction_concordance_pct", value], "%\n")
cat("   Jaccard (padj<0.05 sig. gate):", sm[metric == "jaccard_padj005", value], "\n")
cat("   DEGs unadjusted:", sm[metric == "n_deg_unadjusted", value],
    " | adjusted:", sm[metric == "n_deg_adjusted", value], "\n")
cat("   Gained:", sm[metric == "n_gained", value],
    " | Lost:", sm[metric == "n_lost", value], "\n")

cat("\nB) NAS Score (N=660, 5 datasets):\n")
nm <- results_nas$metrics
cat("   Spearman rho (logFC):", nm[metric == "spearman_rho_logFC", value], "\n")
cat("   Pearson r (logFC):", nm[metric == "pearson_r_logFC", value], "\n")
cat("   Direction concordance:", nm[metric == "direction_concordance_pct", value], "%\n")
cat("   Jaccard (padj<0.05 sig. gate):", nm[metric == "jaccard_padj005", value], "\n")
cat("   DEGs unadjusted:", nm[metric == "n_deg_unadjusted", value],
    " | adjusted:", nm[metric == "n_deg_adjusted", value], "\n")
cat("   Gained:", nm[metric == "n_gained", value],
    " | Lost:", nm[metric == "n_lost", value], "\n")
cat("   Genes with sig NAS effect:", nm[metric == "n_covar_sig_genes", value], "\n")

# Interpretation
rho_val <- as.numeric(nm[metric == "spearman_rho_logFC", value])
jac_val <- as.numeric(nm[metric == "jaccard_padj005", value])

cat("\nINTERPRETATION (NAS-adjusted, primary):\n")
if (rho_val > 0.95 && jac_val > 0.80) {
  cat("  NEGLIGIBLE BMI/steatosis confounding (rho>0.95, Jaccard>0.80)\n")
  cat("  Adding steatosis/NAS as covariate does not materially alter dream DEG results.\n")
} else if (rho_val > 0.85) {
  cat("  MODERATE BMI/steatosis effect (rho 0.85-0.95)\n")
  cat("  Some genes affected -- check top_shifted_genes files for details.\n")
} else {
  cat("  SUBSTANTIAL BMI/steatosis effect (rho<0.85)\n")
  cat("  NAS adjustment meaningfully changes results. Investigate top shifted genes.\n")
}

cat("\nLIMITATIONS:\n")
cat("  1. Steatosis grade available only for GSE130970 (N=76, single dataset).\n")
cat("  2. NAS is an imperfect BMI proxy: NAS includes inflammation + ballooning\n")
cat("     (not purely BMI-driven), so adjustment may over-correct.\n")
cat("  3. NAS and disease status are correlated (Controls have NAS=0), so\n")
cat("     adjusting for NAS partially adjusts for disease -- conservative test.\n")
cat("  4. True BMI data would be a stronger test of this confound.\n")

cat("\nOutput files:\n")
cat("  ", file.path(OUTDIR, "steatosis_sensitivity_comparison.csv"), "\n")
cat("  ", file.path(OUTDIR, "steatosis_sensitivity_metrics.csv"), "\n")
cat("  ", file.path(OUTDIR, "top_shifted_genes_steatosis.csv"), "\n")
cat("  ", file.path(OUTDIR, "top_shifted_genes_nas.csv"), "\n")
cat("  ", file.path(FIGDIR, "figS_steatosis_bmi_sensitivity.pdf"), "\n")

cat("\n=== Steatosis/BMI sensitivity analysis completed:", as.character(Sys.time()), "===\n")
