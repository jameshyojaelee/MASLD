#!/usr/bin/env Rscript
# 06b_fixed_effect_sensitivity.R
# ---------------------------------------------------------------------------
# R2 #2: Fixed-effect model sensitivity
# k=5 datasets is arguably under-identified for (1|dataset) random intercept.
# Fit: ~ group_binary + inferred_sex + dataset  (dataset as FIXED effect)
# via limma-voom (no RE -> dream() is not needed; eBayes() IS required).
# Compare: DEG count, t-stat ratio (FE/RE) per gene, SE distributions.
# Report: estimated sigma2_dataset from the RE model and whether it's near zero.
#
# Input: RNA-seq/results/kallisto/all_cohorts_gene_counts.tsv.gz
# Reference: RNA-seq/results/kallisto/dream_results_kallisto.csv
# Output: RNA-seq/results/audit_sensitivity/fixed_effect_sensitivity/
# ---------------------------------------------------------------------------

set.seed(42)
Sys.setenv(R_PARALLEL_SEED = "42")

suppressPackageStartupMessages({
  library(reformulas); library(lme4); library(data.table); library(edgeR); library(yaml)
})
ns_lme4 <- asNamespace("lme4")
for (fn in c("findbars", "nobars", "subbars", "rebuildFormula")) {
  if (exists(fn, envir = ns_lme4)) try({
    unlockBinding(fn, ns_lme4)
    assign(fn, get(fn, asNamespace("reformulas")), envir = ns_lme4)
    lockBinding(fn, ns_lme4)
  }, silent = TRUE)
}
suppressPackageStartupMessages({ library(variancePartition); library(BiocParallel); library(limma) })

PROJECT <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
RDIR    <- file.path(PROJECT, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
KALL    <- file.path(PROJECT, "RNA-seq/results/kallisto")
OUTDIR  <- file.path(PROJECT, "RNA-seq/results/audit_sensitivity/fixed_effect_sensitivity")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ---- Load kallisto counts ---------------------------------------------------
cat("[", as.character(Sys.time()), "] Loading kallisto counts\n")
kc <- fread(file.path(KALL, "all_cohorts_gene_counts.tsv.gz"))
gene_ids <- kc$gene_id
kc[, gene_id := NULL]
cnts <- as.matrix(kc)
rownames(cnts) <- gene_ids
storage.mode(cnts) <- "double"
cnts[is.na(cnts)] <- 0
cnts <- round(cnts)
rownames(cnts) <- sub("\\..*$", "", rownames(cnts))
cat("  kallisto matrix:", nrow(cnts), "genes x", ncol(cnts), "samples\n")

# ---- Canonical metadata -----------------------------------------------------
cat("[", as.character(Sys.time()), "] Loading metadata\n")
dge <- readRDS(file.path(RDIR, "merged_dge.rds"))
meta_new <- readRDS(file.path(RDIR, "meta_matched.rds"))

yaml_path <- file.path(PROJECT, "config/human_datasets.yaml")
ycfg <- yaml::read_yaml(yaml_path)$datasets
mega_cohorts <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))
cat("  mega cohorts:", paste(mega_cohorts, collapse = ", "), "\n")

ds_samples <- dge$samples
ds_samples$sample_id <- rownames(ds_samples)
keep_mega <- ds_samples$dataset %in% mega_cohorts
canonical_ids <- ds_samples$sample_id[keep_mega]
common_ids <- intersect(canonical_ids, colnames(cnts))
cat("  common samples:", length(common_ids), "\n")

cnts <- cnts[, common_ids]
ds_sub <- ds_samples[match(common_ids, ds_samples$sample_id), ]

# ---- Build DGEList -----------------------------------------------------------
group <- factor(ds_sub$group_binary, levels = c("Control", "Disease"))
y <- DGEList(counts = cnts, samples = data.frame(
  sample_id    = common_ids,
  dataset      = ds_sub$dataset,
  group_binary = ds_sub$group_binary,
  stringsAsFactors = FALSE))
keep_g <- filterByExpr(y, group = group)
y <- y[keep_g, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y, method = "TMM")
cat("  After filterByExpr+TMM:", nrow(y), "genes x", ncol(y), "samples\n")

matched_sex <- meta_new$inferred_sex[match(common_ids, meta_new$sample_id)]
info <- data.frame(
  group_binary = factor(ds_sub$group_binary, levels = c("Control", "Disease")),
  dataset      = droplevels(factor(ds_sub$dataset)),
  inferred_sex = factor(matched_sex),
  stringsAsFactors = FALSE)
rownames(info) <- common_ids

cat("\nGroup x Dataset:\n"); print(table(info$group_binary, info$dataset))

# ---- FE model: limma-voom with dataset as fixed effect -----------------------
form_fe <- ~ group_binary + inferred_sex + dataset
cat("\n[", as.character(Sys.time()), "] Fitting FE model:", deparse(form_fe), "\n")

ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
cat("Using", ncpus, "CPU cores\n")
param <- if (ncpus > 1) MulticoreParam(ncpus) else SerialParam()

design_fe <- model.matrix(form_fe, data = info)
cat("  FE design matrix rank:", qr(design_fe)$rank, "/ ncol:", ncol(design_fe), "\n")

v_fe <- voom(y, design_fe, plot = FALSE)
fit_fe <- lmFit(v_fe, design_fe)
fit_fe <- eBayes(fit_fe)

res_fe <- topTable(fit_fe, coef = "group_binaryDisease", number = Inf, sort.by = "none")
res_fe$gene <- rownames(res_fe)
res_fe_dt <- as.data.table(res_fe)
setnames(res_fe_dt, "adj.P.Val", "padj")

cat("\n===== FE MODEL RESULTS =====\n")
n_fe_03 <- nrow(res_fe_dt[padj < 0.05 & abs(logFC) >= 0.3])
n_fe_05 <- nrow(res_fe_dt[padj < 0.05 & abs(logFC) >= 0.5])
cat("DEGs padj<0.05 & |logFC|>=0.3:", n_fe_03, "\n")
cat("DEGs padj<0.05 & |logFC|>=0.5:", n_fe_05, "\n")

fwrite(res_fe_dt, file.path(OUTDIR, "dream_5cohort_fe_kallisto.csv"))

# ---- RE model: estimate sigma2_dataset from variancePartition ----------------
cat("\n[", as.character(Sys.time()), "] Estimating dataset RE variance (sigma2_dataset)\n")
form_re <- ~ group_binary + inferred_sex + (1|dataset)

# Fit RE voom weights for variance partition
v_re <- suppressWarnings(voomWithDreamWeights(y, form_re, info, BPPARAM = param))

# Extract variance partition to get dataset RE fraction
vp <- tryCatch({
  cat("  Running fitExtractVarPartModel (this may take a while)...\n")
  as.data.frame(fitExtractVarPartModel(v_re, form_re, info, BPPARAM = param))
}, error = function(e) {
  cat("  VarPart failed:", conditionMessage(e), "\n")
  NULL
})

if (!is.null(vp)) {
  vp_dataset_med  <- median(vp$dataset, na.rm = TRUE)
  vp_dataset_mean <- mean(vp$dataset, na.rm = TRUE)
  vp_dataset_q25  <- quantile(vp$dataset, 0.25, na.rm = TRUE)
  vp_dataset_q75  <- quantile(vp$dataset, 0.75, na.rm = TRUE)
  pct_near_zero   <- mean(vp$dataset < 0.01, na.rm = TRUE)
  pct_above_20    <- mean(vp$dataset > 0.20, na.rm = TRUE)
  cat(sprintf("  sigma2_dataset: median=%.4f, mean=%.4f, Q25=%.4f, Q75=%.4f\n",
              vp_dataset_med, vp_dataset_mean, vp_dataset_q25, vp_dataset_q75))
  cat(sprintf("  %% genes with dataset var < 1%%: %.1f%%\n", pct_near_zero * 100))
  cat(sprintf("  %% genes with dataset var > 20%%: %.1f%%\n", pct_above_20 * 100))

  # Save variance partition table
  vp$gene <- rownames(vp)
  fwrite(as.data.table(vp), file.path(OUTDIR, "variance_partition_kallisto.csv"))
} else {
  vp_dataset_med <- vp_dataset_mean <- pct_near_zero <- pct_above_20 <- NA_real_
  vp_dataset_q25 <- vp_dataset_q75 <- NA_real_
}

# ---- Also fit RE dream for direct t-stat comparison --------------------------
cat("\n[", as.character(Sys.time()), "] Fitting RE dream for t-stat comparison...\n")
fit_re <- suppressWarnings(dream(v_re, form_re, info, BPPARAM = param))
res_re <- topTable(fit_re, coef = "group_binaryDisease", number = Inf, sort.by = "none")
res_re$gene <- rownames(res_re)
res_re_dt <- as.data.table(res_re)
setnames(res_re_dt, "adj.P.Val", "padj")

# ---- Load 5-cohort canonical reference (if different from just-fitted RE) ----
ref <- fread(file.path(KALL, "dream_results_kallisto.csv"))
setnames(ref, "adj.P.Val", "padj", skip_absent = TRUE)

# ---- Per-gene comparison (FE vs RE) -----------------------------------------
re_sub <- res_re_dt[, .(gene, logFC_RE = logFC, t_RE = t, padj_RE = padj)]
fe_sub <- res_fe_dt[, .(gene, logFC_FE = logFC, t_FE = t, padj_FE = padj)]
comp <- merge(re_sub, fe_sub, by = "gene")
comp[, t_ratio := t_FE / t_RE]
comp[, se_RE := abs(logFC_RE / t_RE)]
comp[, se_FE := abs(logFC_FE / t_FE)]

fwrite(comp, file.path(OUTDIR, "per_gene_re_vs_fe.csv"))

# DEG sets
deg_re <- comp[padj_RE < 0.05 & abs(logFC_RE) > 0.3, gene]
deg_fe <- comp[padj_FE < 0.05 & abs(logFC_FE) > 0.3, gene]
deg_re_t1 <- comp[padj_RE < 0.05 & abs(logFC_RE) > 0.5, gene]
deg_fe_t1 <- comp[padj_FE < 0.05 & abs(logFC_FE) > 0.5, gene]

jacc <- function(a, b) length(intersect(a, b)) / max(length(union(a, b)), 1)
rho <- cor(comp$logFC_RE, comp$logFC_FE, method = "spearman")
rho_t <- cor(comp$t_RE, comp$t_FE, method = "spearman", use = "complete.obs")

# Direction concordance among shared DEGs
shared <- intersect(deg_re, deg_fe)
dir_conc <- if (length(shared) > 0) {
  sub <- comp[gene %in% shared]
  mean(sign(sub$logFC_RE) == sign(sub$logFC_FE))
} else NA_real_

# t-stat ratio summary
t_ratio_med <- median(comp$t_ratio, na.rm = TRUE)
t_ratio_q25 <- quantile(comp$t_ratio, 0.25, na.rm = TRUE)
t_ratio_q75 <- quantile(comp$t_ratio, 0.75, na.rm = TRUE)

# SE ratio summary
se_ratio <- comp$se_FE / comp$se_RE
se_ratio_med <- median(se_ratio, na.rm = TRUE)

# Check for singular fits in RE model
n_singular <- tryCatch({
  # Count genes where the RE model is singular (sigma2_dataset = 0)
  if (!is.null(vp)) sum(vp$dataset < 1e-6, na.rm = TRUE) else NA_integer_
}, error = function(e) NA_integer_)

# ---- Write REPORT.md --------------------------------------------------------
report_path <- file.path(OUTDIR, "REPORT.md")
sink(report_path)
cat("# Fixed-Effect Model Sensitivity [R2 #2]\n\n")
cat("**Date:**", as.character(Sys.time()), "\n\n")

cat("## Design\n\n")
cat("- **Reviewer concern**: k=5 datasets for `(1|dataset)` random intercept is under-identified\n")
cat("- **Quantifier**: kallisto (tximport gene-level counts)\n")
cat("- **RE model**: `~ group_binary + inferred_sex + (1|dataset)` (dream, Satterthwaite)\n")
cat("- **FE model**: `~ group_binary + inferred_sex + dataset` (limma-voom + eBayes)\n")
cat("- **N samples**:", ncol(y), "across k=", nlevels(info$dataset), "datasets\n\n")

cat("## Dataset RE variance (sigma2_dataset)\n\n")
if (!is.na(vp_dataset_med)) {
  cat(sprintf("- **Median dataset variance fraction**: %.4f (%.1f%% of total)\n", vp_dataset_med, vp_dataset_med * 100))
  cat(sprintf("- **Mean dataset variance fraction**: %.4f\n", vp_dataset_mean))
  cat(sprintf("- **IQR**: [%.4f, %.4f]\n", vp_dataset_q25, vp_dataset_q75))
  cat(sprintf("- **%% genes with dataset var < 1%%**: %.1f%%\n", pct_near_zero * 100))
  cat(sprintf("- **%% genes with dataset var > 20%%**: %.1f%%\n", pct_above_20 * 100))
  if (!is.na(n_singular)) cat(sprintf("- **Genes with effectively zero RE (< 1e-6)**: %d\n", n_singular))
  cat("\n")
  if (vp_dataset_med < 0.05) {
    cat("**Interpretation**: Median dataset variance is < 5% of total — the RE is near-singular\n")
    cat("for most genes. Both RE and FE models should give nearly identical results.\n")
  } else if (vp_dataset_med < 0.15) {
    cat("**Interpretation**: Moderate dataset variance (5-15%). RE is absorbing real batch\n")
    cat("variation; k=5 is low but functional.\n")
  } else {
    cat("**Interpretation**: Substantial dataset variance (> 15%). The RE term is doing\n")
    cat("meaningful work, but k=5 may estimate it imprecisely.\n")
  }
} else {
  cat("Variance partition estimation failed; see log for details.\n")
}

cat("\n## DEG counts (padj < 0.05)\n\n")
cat("| Threshold | RE (dream) | FE (limma-voom) |\n")
cat("|-----------|------------|------------------|\n")
cat(sprintf("| \\|logFC\\| > 0.3 | %d | %d |\n", length(deg_re), length(deg_fe)))
cat(sprintf("| \\|logFC\\| > 0.5 | %d | %d |\n", length(deg_re_t1), length(deg_fe_t1)))

cat("\n## Concordance metrics\n\n")
cat(sprintf("- **LFC Spearman rho (all genes)**: %.4f\n", rho))
cat(sprintf("- **t-stat Spearman rho**: %.4f\n", rho_t))
cat(sprintf("- **Direction concordance (shared DEGs)**: %.4f\n", dir_conc))
cat(sprintf("- **Jaccard (padj<0.05, |LFC|>0.3)**: %.4f (intersect %d / union %d)\n",
            jacc(deg_re, deg_fe), length(intersect(deg_re, deg_fe)), length(union(deg_re, deg_fe))))
cat(sprintf("- **Jaccard (padj<0.05, |LFC|>0.5)**: %.4f (intersect %d / union %d)\n",
            jacc(deg_re_t1, deg_fe_t1), length(intersect(deg_re_t1, deg_fe_t1)), length(union(deg_re_t1, deg_fe_t1))))

cat("\n## t-statistic ratio (FE / RE)\n\n")
cat(sprintf("- **Median t-ratio**: %.4f\n", t_ratio_med))
cat(sprintf("- **IQR**: [%.4f, %.4f]\n", t_ratio_q25, t_ratio_q75))
cat(sprintf("- **SE ratio (FE/RE) median**: %.4f\n", se_ratio_med))
cat("\n")
if (abs(t_ratio_med - 1.0) < 0.05) {
  cat("**Interpretation**: t-ratio is near 1.0 — RE and FE give equivalent inference.\n")
  cat("The k=5 concern is empirically moot for this dataset.\n")
} else if (t_ratio_med > 1.0) {
  cat("**Interpretation**: FE t-stats systematically larger — RE is inflating SE\n")
  cat("(classic symptom of under-identified RE). Consider reporting FE as primary.\n")
} else {
  cat("**Interpretation**: RE t-stats systematically larger — RE is shrinking SE.\n")
  cat("This is expected if dataset effects are real; report both.\n")
}

cat("\n## Verdict\n\n")
if (rho > 0.99 && jacc(deg_re, deg_fe) > 0.95) {
  cat("**RE and FE are functionally equivalent** (rho > 0.99, Jaccard > 0.95).\n")
  cat("The k=5 random intercept is NOT under-identified in practice.\n")
  cat("Report RE as primary (accounts for correlation structure); FE as supplementary.\n")
} else if (rho > 0.95 && jacc(deg_re, deg_fe) > 0.85) {
  cat("**RE and FE are broadly concordant** but with some divergence.\n")
  cat("Report both in supplementary; discuss in Methods.\n")
} else {
  cat("**Material divergence between RE and FE**. Recommend reporting FE as primary\n")
  cat("given the k=5 limitation, with RE as sensitivity.\n")
}

cat("\n## Files\n\n")
cat("- `dream_5cohort_fe_kallisto.csv` — full FE per-gene results\n")
cat("- `per_gene_re_vs_fe.csv` — matched RE/FE comparison per gene\n")
cat("- `variance_partition_kallisto.csv` — per-gene variance partition (if computed)\n")

sink()
cat("\nWrote:", report_path, "\n")
cat("[", as.character(Sys.time()), "] Done.\n")
