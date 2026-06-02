#!/usr/bin/env Rscript
# 05_dream_mega_analysis.R
# ---------------------------------------------------------------------------
# Multi-cohort mega-analysis using dream() from variancePartition.
# Accounts for dataset as a random effect.
# Excludes GSE167523 (no controls) from Disease vs Control contrast.
# Output: results/integration/dream_results.csv
# ---------------------------------------------------------------------------

# ---- Seed pinning (T2.4, 2026-04-22) -----
set.seed(42)
Sys.setenv(R_PARALLEL_SEED = "42")

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

# --- Load merged DGE ---
dge <- readRDS(file.path(RDIR, "merged_dge.rds"))

# --- Cohort selection: enforce config/human_datasets.yaml `include_in_mega` ---
# A cohort enters the Disease-vs-Control mega-analysis only if it has both
# arms (healthy controls and disease samples). Cohorts marked
# `include_in_mega: false` in the yaml are excluded because they have no
# healthy controls (GSE167523, GSE174478, GSE193066, GSE240729).
# (PRJNA512027 — L0/S0 library-prep batches perfectly confounded with disease
# severity — was permanently removed from the pipeline 2026-05-15; matches the
# exclusion in MegaMASLD 2024, Tang & Borlak 2024 Hepatology, Govaere
# 2023 Nat Med/Metab.)
# Without this filter, control-less cohorts contribute Disease samples that
# inflate voom weights / dataset-RE variance without informing the
# group_binary fixed-effect estimate (rho = 0.9992 between 8-cohort and
# 5-cohort dream fits — see results/integration/diagnostic_5cohort/).
yaml_path <- file.path(
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
  "config/human_datasets.yaml")
ycfg <- yaml::read_yaml(yaml_path)$datasets
mega_cohorts <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))
cat("yaml-declared mega cohorts (k =", length(mega_cohorts), "):",
    paste(mega_cohorts, collapse = ", "), "\n")

excluded_by_yaml <- setdiff(unique(dge$samples$dataset), mega_cohorts)
cat("Excluded by yaml include_in_mega=false:",
    paste(excluded_by_yaml, collapse = ", "), "\n")

keep_samples <- dge$samples$dataset %in% mega_cohorts
excl_extra <- Sys.getenv("EXCL_EXTRA", "")
if (nchar(excl_extra) > 0) {
  excl_extra_vec <- trimws(strsplit(excl_extra, ",")[[1]])
  keep_samples <- keep_samples & !dge$samples$dataset %in% excl_extra_vec
  cat("Additional excluded datasets (EXCL_EXTRA):",
      paste(excl_extra_vec, collapse = ", "), "\n")
}
dge_mega <- dge[, keep_samples]
cat("Samples for mega-analysis:", ncol(dge_mega), "\n")
# --- Update metadata with inferred sex ---
# dge$samples might be old, so reload meta_matched to get inferred_sex
meta_new <- readRDS(file.path(RDIR, "meta_matched.rds"))

# Match inferred_sex to dge samples
matched_sex <- meta_new$inferred_sex[match(colnames(dge_mega), meta_new$sample_id)]

# NOTE: Cell-type deconvolution covariates (Hepatocytes, Macrophages) are intentionally
# EXCLUDED from this primary model. Including them absorbs composition-driven disease
# signal (MASLD causes hepatocyte loss + immune infiltration), collapsing DEGs to
# near-zero. Deconvolution-adjusted analysis for classifying hepatocyte-intrinsic vs
# composition-driven DEGs is performed separately in 25_deconv_attribution.R.
info <- data.frame(
  group_binary = factor(dge_mega$samples$group_binary, levels = c("Control", "Disease")),
  dataset = factor(dge_mega$samples$dataset),
  inferred_sex = factor(matched_sex),
  stringsAsFactors = FALSE
)
rownames(info) <- colnames(dge_mega)

cat("Group distribution:\n")
print(table(info$group_binary, info$dataset))
cat("Sex distribution:\n")
print(table(info$inferred_sex, useNA="always"))

# --- dream() ---
# ---------- Model selection ----------
model_choice <- Sys.getenv("DREAM_MODEL", "base")
cat("  Dream model variant:", model_choice, "\n")

form <- switch(model_choice,
  "base"             = ~ group_binary + inferred_sex + (1|dataset),
  "randslope"        = ~ group_binary + inferred_sex + (1 + group_binary | dataset),
  "robust"           = ~ group_binary + inferred_sex + (1|dataset),
  "randslope_robust" = ~ group_binary + inferred_sex + (1 + group_binary | dataset),
  stop("Unknown DREAM_MODEL: ", model_choice)
)
use_robust <- model_choice %in% c("robust", "randslope_robust")
cat("  use_robust:", use_robust, "\n")
cat("\nFormula:", deparse(form), "\n")

# Use parallel processing — detect SLURM_CPUS_PER_TASK or default to 1
ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
cat("Using", ncpus, "CPU cores\n")
param <- if (ncpus > 1) MulticoreParam(ncpus) else SerialParam()
# NOTE: pass param explicitly to dream() and voomWithDreamWeights() rather than
# relying on register(), which sets a global default that dream() may not honour.

if (use_robust) cat("  NOTE: robust estimation not available in this variancePartition version (1.36.3); running standard lmer.\n")

cat("Running voomWithDreamWeights...\n")
v <- suppressWarnings(voomWithDreamWeights(dge_mega, form, info, BPPARAM = param))

cat("Running dream()...\n")
fit <- suppressWarnings(dream(v, form, info, BPPARAM = param))
# NOTE: do NOT call eBayes() after dream(). dream() already computes moderated
# t-statistics via the Satterthwaite approximation; a second eBayes() call
# double-shrinks the variance and produces anti-conservative p-values.

# Extract results for Disease vs Control
coef_name <- "group_binaryDisease"
res <- topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
res$gene <- rownames(res)

res_dt <- as.data.table(res)
setnames(res_dt, "adj.P.Val", "padj")

coef_idx <- which(colnames(fit$coefficients) == coef_name)
gene_order <- res_dt$gene
res_dt[, SE := abs(logFC / t)]
if (!is.null(fit$df.total)) {
  res_dt[, df_total := fit$df.total[match(gene_order, rownames(fit$coefficients))]]
}
if (!is.null(fit$sigma)) {
  stdev_unsc <- fit$stdev.unscaled[match(gene_order, rownames(fit$coefficients)), coef_idx]
  res_dt[, SE_unmoderated := stdev_unsc * fit$sigma[match(gene_order, rownames(fit$coefficients))]]
}

# Summary
sig <- res_dt[padj < 0.1]
sig_lfc <- res_dt[padj < 0.1 & abs(logFC) >= 0.5]
cat("\n===== DREAM MEGA-ANALYSIS RESULTS =====\n")
cat("Total genes tested:", nrow(res_dt), "\n")
cat("DEGs (padj < 0.1):", nrow(sig), "\n")
cat("DEGs (padj < 0.1, |logFC| >= 0.5):", nrow(sig_lfc), "\n")
cat("  Up:", sum(sig_lfc$logFC > 0), "\n")
cat("  Down:", sum(sig_lfc$logFC < 0), "\n")

# NOTE: The "base" model writes to dream_results.csv (no suffix) to maintain
# downstream compatibility. Running this script with DREAM_MODEL=base will
# overwrite the production dream_results.csv. Deliberately choosing a different
# DREAM_MODEL (e.g., "randslope") after benchmarking is how you update production.
# Save — filename includes model variant and any excluded datasets
out_suffix <- if (model_choice == "base") "" else paste0("_", model_choice)
if (nchar(excl_extra) > 0) {
  excl_tag <- paste0("_excl", paste(trimws(strsplit(excl_extra, ",")[[1]]), collapse = "_"))
  out_suffix <- paste0(out_suffix, excl_tag)
}
out_fname <- paste0("dream_results", out_suffix, ".csv")
fwrite(res_dt, file.path(RDIR, out_fname))
cat("\nSaved:", out_fname, "\n")
