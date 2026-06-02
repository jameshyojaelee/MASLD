#!/usr/bin/env Rscript
# 05b_dream_nested_re.R
# ---------------------------------------------------------------------------
# Reviewer-defense (Agent B2, 2026-05-20): Dream mega-analysis with optional
# nested random effect (1 | dataset/library_prep_batch) and expanded output
# columns (df.total, SE, isSingular) alongside the canonical logFC, AveExpr,
# t, P.Value, padj. This script is a SENSITIVITY variant of
# 05_dream_mega_analysis.R; it writes to dream_results_nested_re.csv and
# does NOT overwrite the canonical dream_results.csv used by downstream
# scripts.
#
# Inputs:
#   - results/integration/merged_dge.rds (canonical)
#   - agent_status/B5_metadata_with_batch.csv (optional;
#     if present, supplies library_prep_batch per sample)
#
# Output:
#   - results/integration/dream_results_nested_re.csv
#   - results/integration/dream_results_nested_re_meta.json (model metadata)
# ---------------------------------------------------------------------------

set.seed(42)
Sys.setenv(R_PARALLEL_SEED = "42")

suppressPackageStartupMessages({
  library(reformulas)
  library(lme4)
  library(data.table)
  library(edgeR)
  library(yaml)
  library(jsonlite)
})

# Inject reformulas into lme4 namespace BEFORE loading variancePartition
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

# Reviewer-defense (worktree): read inputs from main project; write outputs
# to the worktree's results dir to keep main repo clean.
INPUT_ROOT  <- Sys.getenv("MASLD_INPUT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
OUTPUT_ROOT <- Sys.getenv("MASLD_OUTPUT_ROOT", INPUT_ROOT)

IN_RDIR  <- file.path(INPUT_ROOT,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
OUT_RDIR <- file.path(OUTPUT_ROOT,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
dir.create(OUT_RDIR, recursive = TRUE, showWarnings = FALSE)
RDIR <- OUT_RDIR  # writes go here
B5_BATCH_FILE <- file.path(OUTPUT_ROOT, "agent_status/B5_metadata_with_batch.csv")
if (!file.exists(B5_BATCH_FILE)) {
  B5_BATCH_FILE <- file.path(INPUT_ROOT, "agent_status/B5_metadata_with_batch.csv")
}

# ---------------------------------------------------------------------------
# Load DGE + metadata
# ---------------------------------------------------------------------------
dge <- readRDS(file.path(IN_RDIR, "merged_dge.rds"))

ycfg <- yaml::read_yaml(file.path(INPUT_ROOT, "config/human_datasets.yaml"))$datasets
mega_cohorts <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))
cat("[B2 nested RE] yaml-declared mega cohorts (k =", length(mega_cohorts), "):",
    paste(mega_cohorts, collapse = ", "), "\n")

keep_samples <- dge$samples$dataset %in% mega_cohorts
dge_mega <- dge[, keep_samples]
cat("Samples for mega-analysis:", ncol(dge_mega), "\n")

meta_new <- readRDS(file.path(IN_RDIR, "meta_matched.rds"))
matched_sex <- meta_new$inferred_sex[match(colnames(dge_mega), meta_new$sample_id)]

# ---------------------------------------------------------------------------
# Load library_prep_batch if B5 has committed it; else fall back
# ---------------------------------------------------------------------------
use_nested <- FALSE
batch_source <- "none"
batch_vec <- rep(NA_character_, ncol(dge_mega))
if (file.exists(B5_BATCH_FILE)) {
  cat("[B2 nested RE] Found B5 batch metadata:", B5_BATCH_FILE, "\n")
  b5 <- fread(B5_BATCH_FILE)
  # Expected columns: sample_id, library_prep_batch
  if (all(c("sample_id", "library_prep_batch") %in% names(b5))) {
    batch_vec <- b5$library_prep_batch[match(colnames(dge_mega), b5$sample_id)]
    n_have <- sum(!is.na(batch_vec) & nzchar(batch_vec))
    cat(sprintf("  library_prep_batch coverage: %d / %d samples (%.1f%%)\n",
                n_have, length(batch_vec), 100 * n_have / length(batch_vec)))
    if (n_have / length(batch_vec) >= 0.8) {
      use_nested <- TRUE
      batch_source <- "B5_metadata_with_batch.csv"
      # Fill NA batch with dataset name (so RE collapses to dataset-only for those)
      batch_vec[is.na(batch_vec) | !nzchar(batch_vec)] <-
        as.character(dge_mega$samples$dataset[is.na(batch_vec) | !nzchar(batch_vec)])
    } else {
      cat("  [WARN] coverage <80%; falling back to dataset-only RE\n")
    }
  } else {
    cat("  [WARN] B5 file missing expected columns; falling back\n")
  }
} else {
  cat("[B2 nested RE] B5 batch file not found at", B5_BATCH_FILE, "\n")
  cat("  Falling back to dataset-only RE (documented in commit).\n")
}

info <- data.frame(
  group_binary  = factor(dge_mega$samples$group_binary, levels = c("Control", "Disease")),
  dataset       = factor(dge_mega$samples$dataset),
  inferred_sex  = factor(matched_sex),
  stringsAsFactors = FALSE
)
if (use_nested) {
  info$library_prep_batch <- factor(batch_vec)
}
rownames(info) <- colnames(dge_mega)

cat("Group x dataset:\n"); print(table(info$group_binary, info$dataset))
cat("Sex:\n"); print(table(info$inferred_sex, useNA = "always"))
if (use_nested) {
  cat("library_prep_batch levels:", nlevels(info$library_prep_batch), "\n")
}

# ---------------------------------------------------------------------------
# Formula
# ---------------------------------------------------------------------------
if (use_nested) {
  form <- ~ group_binary + inferred_sex + (1 | dataset / library_prep_batch)
} else {
  form <- ~ group_binary + inferred_sex + (1 | dataset)
}
cat("\nFormula:", deparse(form), "\n")

ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
cat("Using", ncpus, "CPU cores\n")
param <- if (ncpus > 1) MulticoreParam(ncpus) else SerialParam()

cat("Running voomWithDreamWeights...\n")
v <- suppressWarnings(voomWithDreamWeights(dge_mega, form, info, BPPARAM = param))

cat("Running dream()...\n")
fit <- suppressWarnings(dream(v, form, info, BPPARAM = param))

# ---------------------------------------------------------------------------
# Extract results with df.total, SE, isSingular
# ---------------------------------------------------------------------------
res <- topTable(fit, coef = "group_binaryDisease", number = Inf, sort.by = "none")
res$gene <- rownames(res)
res_dt <- as.data.table(res)
setnames(res_dt, "adj.P.Val", "padj")

# df.total: variancePartition stores per-gene moderated df in fit$df.total
if (!is.null(fit$df.total)) {
  res_dt$df.total <- fit$df.total[match(res_dt$gene, rownames(fit))]
} else if (!is.null(fit$df.residual)) {
  res_dt$df.total <- fit$df.residual[match(res_dt$gene, rownames(fit))]
} else {
  res_dt$df.total <- NA_real_
}

# SE = logFC / t (only where t != 0 and finite)
res_dt[, SE := ifelse(is.finite(t) & t != 0, logFC / t, NA_real_)]

# isSingular: extract per-gene if available
is_sing_vec <- rep(NA, nrow(res_dt))
if (!is.null(fit$isSingular)) {
  is_sing_vec <- fit$isSingular[match(res_dt$gene, rownames(fit))]
} else {
  attr_is <- attr(fit, "isSingular")
  if (!is.null(attr_is)) {
    if (length(attr_is) == nrow(res_dt)) is_sing_vec <- attr_is
    else is_sing_vec <- rep(attr_is[1], nrow(res_dt))
  }
}
res_dt$isSingular <- is_sing_vec

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
sig     <- res_dt[padj < 0.1]
sig_lfc <- res_dt[padj < 0.1 & abs(logFC) >= 0.5]
cat("\n===== B2 NESTED-RE DREAM RESULTS =====\n")
cat("Total genes tested:", nrow(res_dt), "\n")
cat("DEGs (padj<0.1):", nrow(sig), "\n")
cat("DEGs (padj<0.1, |logFC|>=0.5):", nrow(sig_lfc),
    "(up:", sum(sig_lfc$logFC > 0), ", down:", sum(sig_lfc$logFC < 0), ")\n")
if (any(!is.na(res_dt$isSingular))) {
  cat("Singular fits:", sum(res_dt$isSingular, na.rm = TRUE), "/",
      sum(!is.na(res_dt$isSingular)), "\n")
}

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
out_csv  <- file.path(RDIR, "dream_results_nested_re.csv")
out_meta <- file.path(RDIR, "dream_results_nested_re_meta.json")
fwrite(res_dt, out_csv)
jsonlite::write_json(list(
  formula        = deparse(form),
  use_nested     = use_nested,
  batch_source   = batch_source,
  n_samples      = ncol(dge_mega),
  n_genes        = nrow(res_dt),
  n_singular     = if (any(!is.na(res_dt$isSingular))) sum(res_dt$isSingular, na.rm = TRUE) else NA,
  n_deg_padj0.1  = nrow(sig),
  n_deg_lfc_0.5  = nrow(sig_lfc),
  timestamp      = as.character(Sys.time())
), out_meta, auto_unbox = TRUE, pretty = TRUE)
cat("\nSaved:", out_csv, "\n")
cat("Saved meta:", out_meta, "\n")
