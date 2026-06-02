#!/usr/bin/env Rscript
# 130b_c13_categorical_sensitivity.R
# ---------------------------------------------------------------------------
# C13 CATEGORICAL-FIBROSIS SENSITIVITY (Pre-registered HOLD condition C1)
#
# Purpose:
#   The canonical C13 NAFL-vs-NASH contrast (Script 130, line 461) models
#   fibrosis as a LINEAR numeric covariate: `~ nafl_nash + fib_numeric + ...`.
#   The pre-registered HOLD condition C1 (see FIGURE_PLAN_REVISED.md:133 and
#   claim_C1_mediation.md §7) requires a CATEGORICAL sensitivity fit:
#     `~ nafl_nash + factor(fibrosis_stage) + sex_covar + (1 | dataset)`
#   This rules out linear-misspecification-driven over-correction and tests
#   whether the 268-gene "true NASH core" survives non-linear fibrosis
#   adjustment.
#
# Reframed baseline (per mediation §5.5, 2026-04-14):
#   Pure-fibrosis reduction on matched 615 samples:
#     - padj<0.01: 1,346 -> 52 = 96.1%
#     - padj<0.05: 2,724 -> 268 = 90.2% (central estimate)
#     - padj<0.10: 3,904 -> 497 = 87.3%
#   Envelope: 87.3-96.1% across thresholds. THIS RETIRES the stale "95.5%".
#
# Sample set (615):
#   Subset meeting all of:
#     diagnosis_harmonized in {NAFL, NASH, Borderline}
#     !is.na(fibrosis_stage) & fibrosis_stage in 0:4
#     dataset has >=2 NAFL AND >=2 NASH samples
#   This reproduces the C13 subset exactly (5 datasets: GSE130970, GSE135251,
#   GSE162694, GSE174478, GSE193066; NAFL=91, NASH=524).
#
# Upstream:
#   - results/integration/merged_counts_raw.rds (counts matrix)
#   - results/integration/meta_matched.rds (unified metadata)
#   - qc/sample_qc_report.csv (pass_technical flag)
#   - results/staging_classifier/modeling_metadata.csv (fib/NAS/dx)
#   - results/progression/c13_nash_vs_nafl_fib_adj_dream.csv (linear baseline)
#
# Downstream:
#   - results/c13_categorical_sensitivity.csv (per-gene table)
#   - results/c13_categorical_sensitivity_summary.csv (summary row)
#
# Output schema (per-gene):
#   gene_id, symbol, logFC_categorical, padj_categorical,
#   logFC_linear, padj_linear, sign_concordant (TRUE/FALSE), in_linear_core (TRUE/FALSE)
#
# Invocation: sbatch docs/archive/code_review_2026-04-21/fix_progress/run_c13_categorical.sh
# (code_review tree was archived on 2026-05-24; the .sh script remains executable at the archived path)
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(reformulas)
  library(lme4)
  library(data.table)
  library(edgeR)
})

# Force lme4 namespace to use reformulas helpers (same pattern as Scripts 05, 130)
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
ODIR <- file.path(INT, "results")
dir.create(ODIR, recursive = TRUE, showWarnings = FALSE)

# --- Parallel setup ---
ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
cat("Using", ncpus, "CPU cores\n")
BPPARAM <- if (ncpus > 1) MulticoreParam(ncpus, progressbar = FALSE) else SerialParam()

# --- Load data (matches Scripts 13 / 130 / c1_same_sample_set.R) ---
cat("Loading data...\n")
counts <- readRDS(file.path(INT, "results/integration/merged_counts_raw.rds"))
meta   <- readRDS(file.path(INT, "results/integration/meta_matched.rds"))
qc     <- fread(file.path(INT, "qc/sample_qc_report.csv"))
meta   <- meta[sample_id %in% qc[pass_technical == TRUE, sample_id]]

# Annotation
annot_file <- file.path(INT, "results/disease_signatures/unified_disease_signatures.csv")
if (file.exists(annot_file)) {
  gene_annot <- fread(annot_file,
    select = c("gene", "symbol", "gene_type", "mouse_gene_id", "mouse_symbol"))
  gene_annot <- unique(gene_annot, by = "gene")
} else {
  cat("WARNING: Gene annotation file not found.\n")
  gene_annot <- data.table(gene = character(), symbol = character(),
    gene_type = character(), mouse_gene_id = character(), mouse_symbol = character())
}

# Merge NAS/fibrosis/dx
modeling_meta <- fread(file.path(INT, "results/staging_classifier/modeling_metadata.csv"))
meta <- merge(meta, modeling_meta[, .(sample_id, fibrosis_stage, nas_score,
  diagnosis_harmonized, nas_group, fib_ge3, nas_ge5)],
  by = "sample_id", all.x = TRUE, suffixes = c("", ".mm"))
for (col in c("fibrosis_stage", "nas_score", "diagnosis_harmonized")) {
  mm_col <- paste0(col, ".mm")
  if (mm_col %in% names(meta)) {
    na_mask <- is.na(meta[[col]]) | meta[[col]] == ""
    if (any(na_mask)) meta[[col]][na_mask] <- meta[[mm_col]][na_mask]
    meta[, (mm_col) := NULL]
  }
}
meta[, sex_covar := inferred_sex]
na_sex <- is.na(meta$sex_covar) | meta$sex_covar == ""
if (any(na_sex) && "sex" %in% names(meta)) {
  meta$sex_covar[na_sex] <- meta$sex[na_sex]
}
meta[, sex_covar := factor(sex_covar)]

# --- Subset to the 615-sample matched C13 cohort ---
meta_c13 <- meta[diagnosis_harmonized %in% c("NAFL", "NASH", "Borderline") &
                 !is.na(fibrosis_stage) & fibrosis_stage %in% 0:4]
meta_c13[, nafl_nash := fifelse(diagnosis_harmonized == "NAFL", "NAFL", "NASH")]
meta_c13[, nafl_nash := factor(nafl_nash, levels = c("NAFL", "NASH"))]

ds_counts <- meta_c13[, .(n_nafl = sum(nafl_nash == "NAFL"),
                          n_nash = sum(nafl_nash == "NASH")), by = dataset]
valid_ds <- ds_counts[n_nafl >= 2 & n_nash >= 2, dataset]
meta_c13 <- meta_c13[dataset %in% valid_ds]

# Factor version of fibrosis_stage (the sensitivity change)
meta_c13[, fib_factor := factor(fibrosis_stage, levels = 0:4)]

cat(sprintf("\n=== Subset summary ===\n"))
cat(sprintf("N samples = %d (NAFL=%d, NASH=%d)\n",
  nrow(meta_c13),
  sum(meta_c13$nafl_nash == "NAFL"),
  sum(meta_c13$nafl_nash == "NASH")))
cat(sprintf("Valid datasets (%d): %s\n",
  length(valid_ds), paste(sort(valid_ds), collapse = ", ")))

# Per-dataset distribution
cat("\nPer-dataset x diagnosis:\n")
for (ds in sort(unique(meta_c13$dataset))) {
  sub <- meta_c13[dataset == ds]
  cat(sprintf("  %s: NAFL=%d, NASH=%d, F-range=%d-%d\n", ds,
    sum(sub$nafl_nash == "NAFL"), sum(sub$nafl_nash == "NASH"),
    min(sub$fibrosis_stage), max(sub$fibrosis_stage)))
}

# Per-stage distribution (diagnostics)
cat("\nPer-stage x diagnosis:\n")
print(table(meta_c13$fibrosis_stage, meta_c13$nafl_nash, dnn = c("F-stage", "dx")))

cat("\nStages represented:\n")
print(table(meta_c13$fib_factor, useNA = "ifany"))

# --- Build DGEList ---
keep_samples <- intersect(meta_c13$sample_id, colnames(counts))
meta_c13 <- meta_c13[sample_id %in% keep_samples]
dge <- DGEList(counts = counts[, keep_samples])
m_ordered <- meta_c13[match(colnames(dge), meta_c13$sample_id)]
dge$samples <- cbind(dge$samples, m_ordered[, .(dataset, sex_covar, nafl_nash, fib_factor)])
dge$samples$dataset   <- factor(dge$samples$dataset)
dge$samples$sex_covar <- factor(dge$samples$sex_covar)
dge$samples$nafl_nash <- factor(dge$samples$nafl_nash, levels = c("NAFL", "NASH"))
dge$samples$fib_factor <- factor(dge$samples$fib_factor, levels = 0:4)

dge <- calcNormFactors(dge, method = "TMM")
keep_genes <- filterByExpr(dge, group = dge$samples$nafl_nash)
dge <- dge[keep_genes, , keep.lib.sizes = FALSE]
cat(sprintf("\nGenes after filterByExpr: %d\n", nrow(dge)))

# --- Dream fit with CATEGORICAL fibrosis ---
form_cat <- ~ nafl_nash + fib_factor + sex_covar + (1 | dataset)
cat(sprintf("\n=== Formula (categorical): %s ===\n", deparse(form_cat)))

t0 <- Sys.time()
vobj <- suppressWarnings(voomWithDreamWeights(dge, form_cat, dge$samples, BPPARAM = BPPARAM))
fit  <- suppressWarnings(dream(vobj, form_cat, dge$samples, BPPARAM = BPPARAM))
cat(sprintf("Dream fit completed in %.1f min\n",
  as.numeric(Sys.time() - t0, units = "mins")))

res_cat <- topTable(fit, coef = "nafl_nashNASH", number = Inf, sort.by = "none")
res_cat$gene <- rownames(res_cat)
res_cat <- as.data.table(res_cat)
setnames(res_cat, "adj.P.Val", "padj_categorical", skip_absent = TRUE)
setnames(res_cat, "logFC", "logFC_categorical", skip_absent = TRUE)
res_cat <- res_cat[, .(gene, logFC_categorical, padj_categorical)]

# --- Load linear baseline (C13 canonical) ---
linear_file <- file.path(INT, "results/progression/c13_nash_vs_nafl_fib_adj_dream.csv")
cat(sprintf("\nLoading linear-fibrosis baseline: %s\n", linear_file))
res_lin <- fread(linear_file)
res_lin <- res_lin[, .(gene, symbol, gene_type,
                       logFC_linear = logFC,
                       padj_linear  = padj)]

# --- Merge per-gene ---
merged <- merge(res_cat, res_lin, by = "gene", all = TRUE)

# Fill symbol gaps from gene_annot
if (nrow(gene_annot) > 0) {
  missing_sym <- is.na(merged$symbol) | merged$symbol == ""
  if (any(missing_sym)) {
    merged <- merge(merged, gene_annot[, .(gene, symbol_fill = symbol,
                                           gene_type_fill = gene_type)],
                    by = "gene", all.x = TRUE)
    merged[missing_sym & !is.na(symbol_fill), symbol := symbol_fill]
    if ("gene_type" %in% names(merged)) {
      na_gt <- is.na(merged$gene_type) | merged$gene_type == ""
      merged[na_gt & !is.na(gene_type_fill), gene_type := gene_type_fill]
    }
    merged[, c("symbol_fill", "gene_type_fill") := NULL]
  }
}

# --- Derived columns ---
# Sign concordance: only meaningful where BOTH tests have non-missing logFC
merged[, sign_concordant := ifelse(
  is.na(logFC_categorical) | is.na(logFC_linear),
  NA,
  sign(logFC_categorical) == sign(logFC_linear)
)]

# Linear-core membership: 268 genes at padj<0.05 in the linear C13 fit
merged[, in_linear_core := !is.na(padj_linear) & padj_linear < 0.05]

# Categorical significance flag (for summary)
merged[, in_categorical_sig := !is.na(padj_categorical) & padj_categorical < 0.05]

setnames(merged, "gene", "gene_id")

# Final column order
out_cols <- c("gene_id", "symbol", "logFC_categorical", "padj_categorical",
              "logFC_linear", "padj_linear", "sign_concordant", "in_linear_core")
for (cc in out_cols) if (!cc %in% names(merged)) merged[[cc]] <- NA
out_dt <- merged[, ..out_cols]

fwrite(out_dt, file.path(ODIR, "c13_categorical_sensitivity.csv"))
cat(sprintf("\nSaved per-gene table: %s\n",
  file.path(ODIR, "c13_categorical_sensitivity.csv")))

# --- Summary row ---
n_linear_core   <- sum(merged$in_linear_core, na.rm = TRUE)
n_categorical   <- sum(merged$in_categorical_sig, na.rm = TRUE)
n_overlap       <- sum(merged$in_linear_core & merged$in_categorical_sig, na.rm = TRUE)
overlap_pct     <- ifelse(n_linear_core > 0, 100 * n_overlap / n_linear_core, NA)

# Sign concordance across genes significant in EITHER model (inclusive of retained overlap set)
either_sig <- merged$in_linear_core | merged$in_categorical_sig
sign_pct <- ifelse(sum(either_sig, na.rm = TRUE) > 0,
                   100 * sum(merged$sign_concordant[either_sig], na.rm = TRUE) /
                     sum(either_sig & !is.na(merged$sign_concordant)),
                   NA)

# LogFC ratio across linear core (robust measure: median |cat/lin|)
core_ratio <- merged[in_linear_core == TRUE & !is.na(logFC_linear) & logFC_linear != 0,
                     abs(logFC_categorical / logFC_linear)]
mean_logFC_ratio <- mean(core_ratio, na.rm = TRUE)
median_logFC_ratio <- median(core_ratio, na.rm = TRUE)

summary_dt <- data.table(
  n_linear_core        = n_linear_core,
  n_categorical_sig    = n_categorical,
  n_overlap            = n_overlap,
  overlap_pct          = overlap_pct,
  sign_concordance_pct = sign_pct,
  mean_logFC_ratio     = mean_logFC_ratio,
  median_logFC_ratio   = median_logFC_ratio
)

fwrite(summary_dt, file.path(ODIR, "c13_categorical_sensitivity_summary.csv"))
cat(sprintf("Saved summary: %s\n",
  file.path(ODIR, "c13_categorical_sensitivity_summary.csv")))

cat("\n=== SENSITIVITY SUMMARY ===\n")
print(summary_dt)

cat("\n=== Script complete ===\n")
