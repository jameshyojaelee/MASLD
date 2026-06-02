#!/usr/bin/env Rscript
# dream_loo_cv.R
# ---------------------------------------------------------------------------
# Leave-one-study-out cross-validation for dream mega-analysis.
# Parameterized by HELD_OUT env var: excludes one cohort and re-runs dream.
# Modeled on 05_dream_mega_analysis.R.
# Output: results/integration/loo_cv/dream_loo_{HELD_OUT}.csv
# ---------------------------------------------------------------------------

held_out <- Sys.getenv("HELD_OUT", "")
if (nchar(held_out) == 0) stop("HELD_OUT env var must be set (e.g., GSE213621)")
cat("=== Dream LOO-CV: holding out", held_out, "===\n")

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

# --- Output directory ---
loo_dir <- file.path(RDIR, "loo_cv")
dir.create(loo_dir, recursive = TRUE, showWarnings = FALSE)

# --- Load merged DGE ---
cat("Loading merged DGE...\n")
dge <- readRDS(file.path(RDIR, "merged_dge.rds"))

# --- Cohort selection: enforce config/human_datasets.yaml `include_in_mega` ---
yaml_path <- file.path(
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
  "config/human_datasets.yaml")
ycfg <- yaml::read_yaml(yaml_path)$datasets
mega_cohorts <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))

if (!held_out %in% mega_cohorts) {
  stop(held_out, " is not in the mega-analysis cohort set (yaml include_in_mega). ",
       "Choose from: ", paste(mega_cohorts, collapse = ", "))
}

# NOTE: We do NOT re-run filterByExpr() or calcNormFactors() after subsetting,
# matching the approach in 05_dream_mega_analysis.R. This preserves the same gene
# universe and normalization factors as the full model, making LOO results directly
# comparable. voomWithDreamWeights() will recalculate weights for the subset.
keep_samples <- dge$samples$dataset %in% setdiff(mega_cohorts, held_out)
dge_loo <- dge[, keep_samples]
cat("Mega cohorts (k =", length(mega_cohorts), "):",
    paste(mega_cohorts, collapse = ", "), "\n")
cat("Held out:", held_out, "\n")
cat("Samples remaining:", ncol(dge_loo), "\n")
cat("Genes (same as full model):", nrow(dge_loo), "\n")
cat("Datasets remaining:", paste(unique(dge_loo$samples$dataset), collapse = ", "), "\n")

# --- Update metadata with inferred sex ---
meta_new <- readRDS(file.path(RDIR, "meta_matched.rds"))
stopifnot("All DGE samples must be in meta_matched" =
            all(colnames(dge_loo) %in% meta_new$sample_id))
matched_sex <- meta_new$inferred_sex[match(colnames(dge_loo), meta_new$sample_id)]
n_na_sex <- sum(is.na(matched_sex))
if (n_na_sex > 0) {
  warning(n_na_sex, " samples have NA inferred_sex; these will be dropped by dream().")
}

info <- data.frame(
  group_binary = factor(dge_loo$samples$group_binary, levels = c("Control", "Disease")),
  dataset = droplevels(factor(dge_loo$samples$dataset)),
  inferred_sex = factor(matched_sex),
  stringsAsFactors = FALSE
)
rownames(info) <- colnames(dge_loo)

cat("\nGroup distribution (LOO, excluding", held_out, "):\n")
print(table(info$group_binary, info$dataset))
cat("Sex distribution:\n")
print(table(info$inferred_sex, useNA = "always"))

# --- Check we have >= 2 datasets remaining ---
n_datasets <- length(unique(info$dataset))
if (n_datasets < 2) stop("Need >= 2 datasets for random effect; only have ", n_datasets)

# --- dream() ---
form <- ~ group_binary + inferred_sex + (1|dataset)
cat("\nFormula:", deparse(form), "\n")

ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
cat("Using", ncpus, "CPU cores\n")
param <- if (ncpus > 1) MulticoreParam(ncpus) else SerialParam()

cat("Running voomWithDreamWeights...\n")
v <- suppressWarnings(voomWithDreamWeights(dge_loo, form, info, BPPARAM = param))

cat("Running dream()...\n")
fit <- suppressWarnings(dream(v, form, info, BPPARAM = param))
# NOTE: do NOT call eBayes() after dream()

# Extract results
res <- topTable(fit, coef = "group_binaryDisease", number = Inf, sort.by = "none")
res$gene <- rownames(res)
res_dt <- as.data.table(res)
setnames(res_dt, "adj.P.Val", "padj")

# Summary
sig <- res_dt[padj < 0.1]
cat("\n===== LOO-CV RESULTS (excluding", held_out, ") =====\n")
cat("Total genes tested:", nrow(res_dt), "\n")
cat("DEGs (padj < 0.1):", nrow(sig), "\n")
cat("  Up:", sum(sig$logFC > 0), "\n")
cat("  Down:", sum(sig$logFC < 0), "\n")

# Save
out_file <- file.path(loo_dir, paste0("dream_loo_", held_out, ".csv"))
fwrite(res_dt, out_file)
cat("\nSaved:", out_file, "\n")
cat("Done!\n")
