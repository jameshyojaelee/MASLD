#!/usr/bin/env Rscript
# 130_progression_contrasts_dream.R
# ---------------------------------------------------------------------------
# Binary progression contrasts: dream mega-analysis for 6 new contrasts.
#
# Contrasts computed:
#   C3: F3-F4 vs F0-F2 (advanced vs early fibrosis, including F2 in early)
#   C4: NAFL vs Control (disease initiation)
#   C5: NAS >= 5 vs NAS < 5 (clinical NASH threshold)
#   C6: NASH + F3-F4 vs NAFL + F0-F1 (extreme endpoints)
#   C8: F4 vs F0-F3 (cirrhosis binary)
#   C9: F2-F4 vs F0-F1 (F2 inflection point)
#
# Existing contrasts NOT re-run here:
#   C1: MASLD vs Control → Script 05
#   C2: NAFL vs NASH → Script 13
#
# Uses dream() from variancePartition with dataset as random intercept.
# Follows patterns from Scripts 05 and 15b.
#
# NOTE: PRJNA512027 was permanently removed from the pipeline 2026-05-15
# (L0/S0 library batch perfectly confounded with disease severity).
#
# NOTE: This script correctly avoids calling eBayes() after dream().
# Script 13 (NAFL vs NASH) does call eBayes() after dream() — this is
# a latent bug in Script 13 that produces slightly anti-conservative
# p-values. Results from this script are statistically correct.
#
# Output: results/progression/{contrast_id}_dream.csv
#
# DEG THRESHOLD SYSTEM (Two-Tier):
#   Tier 1 (Primary): padj < 0.05, |logFC| > 0.5 — main dream DEGs (Script 05b; migrated 0.3 -> 0.5)
#   Tier 2 (Progression): padj < 0.05, no LFC filter — binary and adjacent contrasts
#   Rationale: Binary contrasts pool multiple stages, diluting per-gene fold changes.
#              Adjacent transitions have lower N (200-400 vs 1,444), making LFC estimates
#              noisier. Cross-contrast comparisons use rank-based enrichment (fgsea)
#              to avoid confounding power with biology.
#   See: figures/supplementary/sensitivity/figS_deg_threshold_landscape.pdf
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(reformulas)
  library(lme4)
  library(data.table)
  library(yaml)
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
})

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts"
INT  <- file.path(BASE, "analysis/integration")
ODIR <- file.path(INT, "results/progression")
dir.create(ODIR, recursive = TRUE, showWarnings = FALSE)

# --- Parallel setup ---
ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
cat("Using", ncpus, "CPU cores\n")
BPPARAM <- if (ncpus > 1) MulticoreParam(ncpus, progressbar = TRUE) else SerialParam()

# --- Load data ---
cat("Loading data...\n")
counts   <- readRDS(file.path(INT, "results/integration/merged_counts_raw.rds"))
meta     <- readRDS(file.path(INT, "results/integration/meta_matched.rds"))
qc       <- fread(file.path(INT, "qc/sample_qc_report.csv"))
meta     <- meta[sample_id %in% qc[pass_technical == TRUE, sample_id]]

# Gene annotations from unified disease signatures (dream_results.csv lacks annotations)
annot_file <- file.path(INT, "results/disease_signatures/unified_disease_signatures.csv")
if (file.exists(annot_file)) {
  gene_annot <- fread(annot_file,
    select = c("gene", "symbol", "gene_type", "mouse_gene_id", "mouse_symbol"))
  gene_annot <- unique(gene_annot, by = "gene")
} else {
  cat("WARNING: Gene annotation file not found — outputs will lack symbol/gene_type columns\n")
  gene_annot <- data.table(gene = character(), symbol = character(),
    gene_type = character(), mouse_gene_id = character(), mouse_symbol = character())
}

# Load modeling metadata for NAS/fibrosis columns
modeling_meta <- fread(file.path(INT, "results/staging_classifier/modeling_metadata.csv"))

# Merge NAS/fibrosis annotations into meta
meta <- merge(meta, modeling_meta[, .(sample_id, fibrosis_stage, nas_score,
  diagnosis_harmonized, nas_group, fib_ge3, nas_ge5)],
  by = "sample_id", all.x = TRUE, suffixes = c("", ".mm"))

# Resolve duplicate columns — prefer modeling_metadata values
for (col in c("fibrosis_stage", "nas_score", "diagnosis_harmonized")) {
  mm_col <- paste0(col, ".mm")
  if (mm_col %in% names(meta)) {
    # Use modeling_metadata values where available
    na_mask <- is.na(meta[[col]]) | meta[[col]] == ""
    if (any(na_mask)) meta[[col]][na_mask] <- meta[[mm_col]][na_mask]
    meta[, (mm_col) := NULL]
  }
}

# Prepare sex covariate
meta[, sex_covar := inferred_sex]
na_sex <- is.na(meta$sex_covar) | meta$sex_covar == ""
if (any(na_sex) && "sex" %in% names(meta)) {
  meta$sex_covar[na_sex] <- meta$sex[na_sex]
}
meta[, sex_covar := factor(sex_covar)]

cat(sprintf("Total QC-passing samples: %d\n", nrow(meta)))
cat(sprintf("Samples with fibrosis staging: %d\n", sum(!is.na(meta$fibrosis_stage))))
cat(sprintf("Samples with NAS score: %d\n", sum(!is.na(meta$nas_score))))
cat(sprintf("Samples with diagnosis_harmonized: %d\n",
  sum(!is.na(meta$diagnosis_harmonized) & meta$diagnosis_harmonized != "")))

# ============================================================
# Helper: run dream for a binary contrast
# ============================================================
run_dream_contrast <- function(contrast_id, meta_sub, group_col, group_levels,
                               formula_str, coef_name, description) {
  cat(sprintf("\n%s\n", strrep("=", 60)))
  cat(sprintf("  CONTRAST %s: %s\n", contrast_id, description))
  cat(sprintf("%s\n", strrep("=", 60)))

  # Ensure factor with correct levels
  meta_sub[[group_col]] <- factor(meta_sub[[group_col]], levels = group_levels)

  # Per-dataset sample counts
  cat("  Per-dataset distribution:\n")
  for (ds in sort(unique(meta_sub$dataset))) {
    sub <- meta_sub[dataset == ds]
    n_per_group <- table(sub[[group_col]])
    cat(sprintf("    %s: %s\n", ds,
      paste(paste0(names(n_per_group), "=", n_per_group), collapse = ", ")))
  }

  # Filter datasets with < 2 samples in either group
  ds_counts <- meta_sub[, .(
    n1 = sum(.SD[[group_col]] == group_levels[1]),
    n2 = sum(.SD[[group_col]] == group_levels[2])
  ), by = dataset, .SDcols = group_col]
  valid_ds <- ds_counts[n1 >= 2 & n2 >= 2, dataset]
  if (length(valid_ds) < 2) {
    cat(sprintf("  WARNING: Only %d valid datasets — need >= 2 for random effect.\n",
      length(valid_ds)))
    if (length(valid_ds) == 0) {
      cat("  SKIPPED: no valid datasets.\n")
      return(NULL)
    }
  }
  meta_sub <- meta_sub[dataset %in% valid_ds]
  cat(sprintf("  Valid datasets: %d (%s)\n", length(valid_ds), paste(valid_ds, collapse = ", ")))

  # Build DGEList from raw counts
  keep_samples <- intersect(meta_sub$sample_id, colnames(counts))
  meta_sub <- meta_sub[sample_id %in% keep_samples]
  dge <- DGEList(counts = counts[, keep_samples])

  # Attach metadata
  m_ordered <- meta_sub[match(colnames(dge), meta_sub$sample_id)]
  dge$samples <- cbind(dge$samples, m_ordered[, .(
    dataset, sex_covar)])
  dge$samples[[group_col]] <- factor(m_ordered[[group_col]], levels = group_levels)
  dge$samples$dataset <- factor(dge$samples$dataset)
  dge$samples$sex_covar <- factor(dge$samples$sex_covar)

  # TMM normalization
  dge <- calcNormFactors(dge, method = "TMM")

  # Filter low-expression genes
  keep_genes <- filterByExpr(dge, group = dge$samples[[group_col]])
  dge <- dge[keep_genes, , keep.lib.sizes = FALSE]

  n1 <- sum(dge$samples[[group_col]] == group_levels[1])
  n2 <- sum(dge$samples[[group_col]] == group_levels[2])
  cat(sprintf("  Samples: %d (%s=%d, %s=%d)\n",
    ncol(dge), group_levels[1], n1, group_levels[2], n2))
  cat(sprintf("  Genes: %d\n", nrow(dge)))

  # Parse formula
  form <- as.formula(formula_str)
  cat(sprintf("  Formula: %s\n", formula_str))

  # Determine if we need random effects
  has_random <- grepl("\\|", formula_str)

  if (has_random && length(valid_ds) < 2) {
    # Single dataset: drop dataset term entirely (adding a constant-valued
    # fixed effect would create a singular design matrix)
    formula_str <- gsub(" \\+ \\(1 \\| dataset\\)", "", formula_str)
    form <- as.formula(formula_str)
    has_random <- FALSE
    cat(sprintf("  Fallback (single dataset, no dataset term): %s\n", formula_str))
  }

  # Run dream (or limma for fixed-effects-only)
  if (has_random) {
    vobj <- suppressWarnings(voomWithDreamWeights(dge, form, dge$samples, BPPARAM = BPPARAM))
    fit <- suppressWarnings(dream(vobj, form, dge$samples, BPPARAM = BPPARAM))
    # NOTE: do NOT call eBayes() after dream() — dream() already computes
    # moderated t-statistics via Satterthwaite approximation.
  } else {
    design <- model.matrix(form, data = dge$samples)
    vobj <- voom(dge, design, plot = FALSE)
    fit <- lmFit(vobj, design)
    fit <- eBayes(fit)
  }

  # Extract results
  res <- topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
  res$gene <- rownames(res)
  res_dt <- as.data.table(res)

  # Add gene annotations
  res_dt <- merge(res_dt, gene_annot, by = "gene", all.x = TRUE)

  # Rename for consistency with Script 05 output convention
  setnames(res_dt, "adj.P.Val", "padj", skip_absent = TRUE)

  # Add contrast metadata
  res_dt[, contrast_id := contrast_id]

  # Summary (Tier 2: padj<0.05 primary, no LFC filter)
  n_sig_05 <- sum(res_dt$padj < 0.05, na.rm = TRUE)
  n_sig_01 <- sum(res_dt$padj < 0.1, na.rm = TRUE)
  n_up <- sum(res_dt$padj < 0.05 & res_dt$logFC > 0, na.rm = TRUE)
  n_down <- sum(res_dt$padj < 0.05 & res_dt$logFC < 0, na.rm = TRUE)
  n_sig_05_lfc02 <- sum(res_dt$padj < 0.05 & abs(res_dt$logFC) > 0.2, na.rm = TRUE)
  n_sig_05_lfc03 <- sum(res_dt$padj < 0.05 & abs(res_dt$logFC) > 0.3, na.rm = TRUE)
  n_sig_05_lfc05 <- sum(res_dt$padj < 0.05 & abs(res_dt$logFC) > 0.5, na.rm = TRUE)
  cat(sprintf("  DEGs (padj<0.05): %d (Up: %d, Down: %d)\n", n_sig_05, n_up, n_down))
  cat(sprintf("  DEGs (padj<0.05, |logFC|>0.2/0.3/0.5): %d / %d / %d\n",
              n_sig_05_lfc02, n_sig_05_lfc03, n_sig_05_lfc05))
  cat(sprintf("  DEGs (padj<0.1): %d\n", n_sig_01))

  # Save
  out_file <- file.path(ODIR, sprintf("%s_dream.csv", contrast_id))
  fwrite(res_dt, out_file)
  cat(sprintf("  Saved: %s\n", basename(out_file)))

  gc()
  return(res_dt)
}

# ============================================================
# C3: F3-F4 vs F0-F2 (Advanced vs Early Fibrosis)
# Note: Script 15b does F3-F4 vs F0-F1 (drops F2). This includes F2 in early.
# ============================================================
meta_c3 <- meta[!is.na(fibrosis_stage) & fibrosis_stage %in% 0:4]
meta_c3[, fib_group := fifelse(fibrosis_stage >= 3, "Advanced", "Early")]

c3_res <- run_dream_contrast(
  contrast_id = "c3_adv_vs_early_fib",
  meta_sub = meta_c3,
  group_col = "fib_group",
  group_levels = c("Early", "Advanced"),
  formula_str = "~ fib_group + sex_covar + (1 | dataset)",
  coef_name = "fib_groupAdvanced",
  description = "F3-F4 vs F0-F2 (Advanced vs Early Fibrosis, F2 included in Early)"
)

# ============================================================
# C4: NAFL vs Control (Disease Initiation)
# ============================================================
# Mega-analysis cohort selection: enforce config/human_datasets.yaml include_in_mega.
# This is the canonical 5-cohort set for any contrast involving Control samples.
ycfg <- yaml::read_yaml(file.path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design", "config/human_datasets.yaml"))$datasets
mega_cohorts <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))
excl_dvc <- setdiff(unique(meta$dataset), mega_cohorts)

stopifnot("group_binary" %in% names(meta))
meta_c4 <- meta[!is.na(diagnosis_harmonized) & diagnosis_harmonized != ""]
meta_c4 <- meta_c4[diagnosis_harmonized == "NAFL" | group_binary == "Control"]
meta_c4 <- meta_c4[!dataset %in% excl_dvc]
meta_c4[, nafl_ctrl := fifelse(group_binary == "Control", "Control", "NAFL")]

c4_res <- run_dream_contrast(
  contrast_id = "c4_nafl_vs_ctrl",
  meta_sub = meta_c4,
  group_col = "nafl_ctrl",
  group_levels = c("Control", "NAFL"),
  formula_str = "~ nafl_ctrl + sex_covar + (1 | dataset)",
  coef_name = "nafl_ctrlNAFL",
  description = "NAFL vs Control (Disease Initiation)"
)

# ============================================================
# C5: NAS >= 5 vs NAS < 5 (Clinical NASH Threshold)
# ============================================================
meta_c5 <- meta[!is.na(nas_score)]
meta_c5[, nas_binary := fifelse(nas_score >= 5, "NAS_high", "NAS_low")]

c5_res <- run_dream_contrast(
  contrast_id = "c5_nas_ge5_vs_lt5",
  meta_sub = meta_c5,
  group_col = "nas_binary",
  group_levels = c("NAS_low", "NAS_high"),
  formula_str = "~ nas_binary + sex_covar + (1 | dataset)",
  coef_name = "nas_binaryNAS_high",
  description = "NAS >= 5 vs NAS < 5 (Clinical NASH Threshold)"
)

# ============================================================
# C6: NASH + F3-F4 vs NAFL + F0-F1 (Extreme Endpoints)
# Borderline (NAS 3-4) grouped with NASH per Kleiner et al. 2005,
# consistent with Script 13 classification.
# ============================================================
meta_c6 <- meta[!is.na(diagnosis_harmonized) & !is.na(fibrosis_stage)]
meta_c6 <- meta_c6[
  (diagnosis_harmonized %in% c("NASH", "Borderline") & fibrosis_stage >= 3) |
  (diagnosis_harmonized == "NAFL" & fibrosis_stage <= 1)
]
meta_c6[, extreme_group := fifelse(
  diagnosis_harmonized == "NAFL" & fibrosis_stage <= 1, "Early_NAFL", "Advanced_NASH")]

c6_res <- run_dream_contrast(
  contrast_id = "c6_extreme_endpoints",
  meta_sub = meta_c6,
  group_col = "extreme_group",
  group_levels = c("Early_NAFL", "Advanced_NASH"),
  formula_str = "~ extreme_group + sex_covar + (1 | dataset)",
  coef_name = "extreme_groupAdvanced_NASH",
  description = "NASH + F3-F4 vs NAFL + F0-F1 (Extreme Endpoints)"
)

# ============================================================
# C8: F4 vs F0-F3 (Cirrhosis Binary)
# ============================================================
meta_c8 <- meta[!is.na(fibrosis_stage) & fibrosis_stage %in% 0:4]
meta_c8[, cirrhosis := fifelse(fibrosis_stage == 4, "Cirrhotic", "Non_cirrhotic")]

c8_res <- run_dream_contrast(
  contrast_id = "c8_cirrhosis",
  meta_sub = meta_c8,
  group_col = "cirrhosis",
  group_levels = c("Non_cirrhotic", "Cirrhotic"),
  formula_str = "~ cirrhosis + sex_covar + (1 | dataset)",
  coef_name = "cirrhosisCirrhotic",
  description = "F4 vs F0-F3 (Cirrhosis Binary)"
)

# ============================================================
# C9: F2-F4 vs F0-F1 (F2 Inflection Point)
# ============================================================
meta_c9 <- meta[!is.na(fibrosis_stage) & fibrosis_stage %in% 0:4]
meta_c9[, f2_group := fifelse(fibrosis_stage >= 2, "Significant_fib", "Minimal_fib")]

c9_res <- run_dream_contrast(
  contrast_id = "c9_f2_inflection",
  meta_sub = meta_c9,
  group_col = "f2_group",
  group_levels = c("Minimal_fib", "Significant_fib"),
  formula_str = "~ f2_group + sex_covar + (1 | dataset)",
  coef_name = "f2_groupSignificant_fib",
  description = "F2-F4 vs F0-F1 (F2 Inflection Point)"
)

# ============================================================
# C11: NASH vs Control (Pure NASH initiation signature)
# Excludes GSE167523 (same as C4).
# Distinct from C1 (which pools NAFL+NASH as "Disease").
# ============================================================
stopifnot("group_binary" %in% names(meta))
meta_c11 <- meta[!is.na(diagnosis_harmonized) & diagnosis_harmonized != ""]
meta_c11 <- meta_c11[diagnosis_harmonized %in% c("NASH", "Borderline") | group_binary == "Control"]
meta_c11 <- meta_c11[!dataset %in% excl_dvc]
meta_c11[, nash_ctrl := fifelse(group_binary == "Control", "Control", "NASH")]

c11_res <- run_dream_contrast(
  contrast_id = "c11_nash_vs_ctrl",
  meta_sub = meta_c11,
  group_col = "nash_ctrl",
  group_levels = c("Control", "NASH"),
  formula_str = "~ nash_ctrl + sex_covar + (1 | dataset)",
  coef_name = "nash_ctrlNASH",
  description = "NASH (including Borderline) vs Control"
)

# ============================================================
# C12: Early NASH (F0-F2) vs Late NASH (F3-F4)
# Within-NASH fibrosis progression. F2 included in Early
# (more power, clinically defensible). Distinct from C3
# which includes NAFL patients.
# ============================================================
meta_c12 <- meta[diagnosis_harmonized %in% c("NASH", "Borderline") &
                 !is.na(fibrosis_stage) & fibrosis_stage %in% 0:4]
meta_c12[, nash_fib_group := fifelse(fibrosis_stage <= 2, "Early_NASH", "Late_NASH")]

c12_res <- run_dream_contrast(
  contrast_id = "c12_early_vs_late_nash",
  meta_sub = meta_c12,
  group_col = "nash_fib_group",
  group_levels = c("Early_NASH", "Late_NASH"),
  formula_str = "~ nash_fib_group + sex_covar + (1 | dataset)",
  coef_name = "nash_fib_groupLate_NASH",
  description = "Early NASH (F0-F2) vs Late NASH (F3-F4) — within-NASH fibrosis progression"
)

# ============================================================
# C13: NASH vs NAFL (fibrosis-adjusted)
# Pure steatohepatitis effect AFTER adjusting for fibrosis stage.
# Cannot use run_dream_contrast() because formula has extra covariate.
# This is the highest-novelty contrast: no competitor has this.
# ============================================================
cat(sprintf("\n%s\n", strrep("=", 60)))
cat("  CONTRAST c13_nash_vs_nafl_fib_adj: NASH vs NAFL (fibrosis-adjusted)\n")
cat(sprintf("%s\n", strrep("=", 60)))

meta_c13 <- meta[diagnosis_harmonized %in% c("NAFL", "NASH", "Borderline") &
                 !is.na(fibrosis_stage) & fibrosis_stage %in% 0:4]
meta_c13[, nafl_nash := fifelse(diagnosis_harmonized == "NAFL", "NAFL", "NASH")]
meta_c13[, nafl_nash := factor(nafl_nash, levels = c("NAFL", "NASH"))]
meta_c13[, fib_numeric := as.numeric(fibrosis_stage)]

# Per-dataset sample counts
cat("  Per-dataset distribution:\n")
for (ds in sort(unique(meta_c13$dataset))) {
  sub <- meta_c13[dataset == ds]
  cat(sprintf("    %s: NAFL=%d, NASH=%d, F-range=%d-%d\n", ds,
    sum(sub$nafl_nash == "NAFL"), sum(sub$nafl_nash == "NASH"),
    min(sub$fibrosis_stage), max(sub$fibrosis_stage)))
}

# Filter datasets with >= 2 samples in BOTH diagnosis groups
ds_counts <- meta_c13[, .(n_nafl = sum(nafl_nash == "NAFL"),
                          n_nash = sum(nafl_nash == "NASH")), by = dataset]
valid_ds <- ds_counts[n_nafl >= 2 & n_nash >= 2, dataset]
meta_c13 <- meta_c13[dataset %in% valid_ds]
cat(sprintf("  Valid datasets: %d (%s)\n", length(valid_ds), paste(valid_ds, collapse = ", ")))

# Build DGEList
keep_samples <- intersect(meta_c13$sample_id, colnames(counts))
meta_c13 <- meta_c13[sample_id %in% keep_samples]
dge_c13 <- DGEList(counts = counts[, keep_samples])
m_c13 <- meta_c13[match(colnames(dge_c13), meta_c13$sample_id)]
dge_c13$samples <- cbind(dge_c13$samples, m_c13[, .(dataset, sex_covar, nafl_nash, fib_numeric)])
dge_c13$samples$dataset <- factor(dge_c13$samples$dataset)
dge_c13$samples$sex_covar <- factor(dge_c13$samples$sex_covar)
dge_c13$samples$nafl_nash <- factor(dge_c13$samples$nafl_nash, levels = c("NAFL", "NASH"))

dge_c13 <- calcNormFactors(dge_c13, method = "TMM")
keep_genes <- filterByExpr(dge_c13, group = dge_c13$samples$nafl_nash)
dge_c13 <- dge_c13[keep_genes, , keep.lib.sizes = FALSE]

n_nafl <- sum(dge_c13$samples$nafl_nash == "NAFL")
n_nash <- sum(dge_c13$samples$nafl_nash == "NASH")
cat(sprintf("  Samples: %d (NAFL=%d, NASH=%d)\n", ncol(dge_c13), n_nafl, n_nash))
cat(sprintf("  Genes: %d\n", nrow(dge_c13)))

# Dream with fibrosis as additional covariate
form_c13 <- ~ nafl_nash + fib_numeric + sex_covar + (1 | dataset)
cat(sprintf("  Formula: %s\n", deparse(form_c13)))

vobj_c13 <- suppressWarnings(voomWithDreamWeights(dge_c13, form_c13, dge_c13$samples, BPPARAM = BPPARAM))
fit_c13 <- suppressWarnings(dream(vobj_c13, form_c13, dge_c13$samples, BPPARAM = BPPARAM))

res_c13 <- topTable(fit_c13, coef = "nafl_nashNASH", number = Inf, sort.by = "none")
res_c13$gene <- rownames(res_c13)
c13_dt <- as.data.table(res_c13)
setnames(c13_dt, "adj.P.Val", "padj", skip_absent = TRUE)
c13_dt <- merge(c13_dt, gene_annot, by = "gene", all.x = TRUE)
c13_dt[, contrast_id := "c13_nash_vs_nafl_fib_adj"]

# Summary (Tier 2: padj<0.05 primary, no LFC filter)
n_sig_05 <- sum(c13_dt$padj < 0.05, na.rm = TRUE)
n_sig_01 <- sum(c13_dt$padj < 0.1, na.rm = TRUE)
n_up <- sum(c13_dt$padj < 0.05 & c13_dt$logFC > 0, na.rm = TRUE)
n_down <- sum(c13_dt$padj < 0.05 & c13_dt$logFC < 0, na.rm = TRUE)
n_sig_05_lfc02 <- sum(c13_dt$padj < 0.05 & abs(c13_dt$logFC) > 0.2, na.rm = TRUE)
n_sig_05_lfc03 <- sum(c13_dt$padj < 0.05 & abs(c13_dt$logFC) > 0.3, na.rm = TRUE)
n_sig_05_lfc05 <- sum(c13_dt$padj < 0.05 & abs(c13_dt$logFC) > 0.5, na.rm = TRUE)
cat(sprintf("  DEGs (padj<0.05): %d (Up: %d, Down: %d)\n", n_sig_05, n_up, n_down))
cat(sprintf("  DEGs (padj<0.05, |logFC|>0.2/0.3/0.5): %d / %d / %d\n",
            n_sig_05_lfc02, n_sig_05_lfc03, n_sig_05_lfc05))
cat(sprintf("  DEGs (padj<0.1): %d\n", n_sig_01))

fwrite(c13_dt, file.path(ODIR, "c13_nash_vs_nafl_fib_adj_dream.csv"))
cat("  Saved: c13_nash_vs_nafl_fib_adj_dream.csv\n")
c13_res <- c13_dt
gc()

# ============================================================
# Summary across all contrasts
# ============================================================
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("  PROGRESSION CONTRAST SUMMARY\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

all_results <- list(c3_res, c4_res, c5_res, c6_res, c8_res, c9_res,
                    c11_res, c12_res, c13_res)
names(all_results) <- c("C3_AdvFib", "C4_NAFLvsCtrl", "C5_NAS5", "C6_Extreme",
                         "C8_Cirrhosis", "C9_F2Inflection",
                         "C11_NASHvsCtrl", "C12_EarlyLateNASH", "C13_FibAdj")

summary_table <- rbindlist(lapply(names(all_results), function(nm) {
  r <- all_results[[nm]]
  if (is.null(r)) {
    cat(sprintf("  %s: SKIPPED\n", nm))
    return(data.table(contrast = nm, n_genes = 0, n_deg_05 = 0,
      n_deg_05_lfc02 = 0, n_deg_05_lfc03 = 0, n_deg_05_lfc05 = 0,
      n_deg_01 = 0, n_up = 0, n_down = 0))
  }
  data.table(
    contrast = nm,
    n_genes = nrow(r),
    n_deg_05 = sum(r$padj < 0.05, na.rm = TRUE),
    n_deg_05_lfc02 = sum(r$padj < 0.05 & abs(r$logFC) > 0.2, na.rm = TRUE),
    n_deg_05_lfc03 = sum(r$padj < 0.05 & abs(r$logFC) > 0.3, na.rm = TRUE),
    n_deg_05_lfc05 = sum(r$padj < 0.05 & abs(r$logFC) > 0.5, na.rm = TRUE),
    n_deg_01 = sum(r$padj < 0.1, na.rm = TRUE),
    n_up = sum(r$padj < 0.05 & r$logFC > 0, na.rm = TRUE),
    n_down = sum(r$padj < 0.05 & r$logFC < 0, na.rm = TRUE)
  )
}))

print(summary_table)

# Save summary
fwrite(summary_table, file.path(ODIR, "progression_contrast_summary.csv"))
cat("\nSaved: progression_contrast_summary.csv\n")

cat("\n=== Script 130 complete ===\n")
