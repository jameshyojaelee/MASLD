#!/usr/bin/env Rscript
# 132_sex_stratified_progression.R
# ---------------------------------------------------------------------------
# Sex-stratified NAFL-vs-NASH analysis (C10 contrast).
#
# Runs THREE dream models to characterise sex differences in the
# NAFL -> NASH transition:
#
#   1. Interaction model (all samples):
#      ~ nafl_nash * sex_covar + (1 | dataset)
#      Coefficient of interest: nafl_nashNASH:sex_covarM
#      Tests whether the NAFL -> NASH transition differs by sex.
#
#   2. Female-only model (females only):
#      ~ nafl_nash + (1 | dataset)
#      Coefficient of interest: nafl_nashNASH
#
#   3. Male-only model (males only):
#      ~ nafl_nash + (1 | dataset)
#      Coefficient of interest: nafl_nashNASH
#
# After all three models, classifies each gene into one of:
#   - female_progression:        significant in female only (padj<0.1)
#   - male_progression:          significant in male only (padj<0.1)
#   - shared_progression:        significant in both, same direction
#   - sex_divergent_progression: significant in both, opposite direction
#   - sex_differential:          interaction padj<0.1 AND |logFC_F - logFC_M| > 0.5
#   - NS:                        not significant in either
#
# Data setup follows Script 13:
#   - Subsets to NAFL + NASH + Borderline samples (diagnosis_harmonized)
#   - Maps to binary nafl_nash: NAFL vs NASH (Borderline grouped with NASH)
# NOTE: This script correctly avoids calling eBayes() after dream().
# dream() already computes moderated t-statistics via the Satterthwaite
# approximation.
#
# Output (to results/progression/):
#   c10_sex_interaction_nafl_nash.csv    — interaction model results
#   c10_female_nafl_nash.csv             — female-only results
#   c10_male_nafl_nash.csv               — male-only results
#   c10_sex_progression_classification.csv — gene classification
#   c10_sex_progression_summary.csv      — summary table
# ---------------------------------------------------------------------------

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
})

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts"
INT  <- file.path(BASE, "analysis/integration")
ODIR <- file.path(INT, "results/progression")
dir.create(ODIR, recursive = TRUE, showWarnings = FALSE)

# --- Parallel setup ---
ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
cat("Using", ncpus, "CPU cores\n")
BPPARAM <- if (ncpus > 1) MulticoreParam(ncpus, progressbar = TRUE) else SerialParam()

# ============================================================
# Load data
# ============================================================
cat("Loading data...\n")
counts <- readRDS(file.path(INT, "results/integration/merged_counts_raw.rds"))
meta   <- readRDS(file.path(INT, "results/integration/meta_matched.rds"))
qc     <- fread(file.path(INT, "qc/sample_qc_report.csv"))
meta   <- meta[sample_id %in% qc[pass_technical == TRUE, sample_id]]

# Apply pass_sex filter (consistent with Script 26 for sex-stratified analysis)
if ("pass_sex" %in% names(qc)) {
  sex_pass_ids <- qc[pass_sex == TRUE, sample_id]
  n_before <- nrow(meta)
  meta <- meta[sample_id %in% sex_pass_ids]
  cat(sprintf("pass_sex filter: %d -> %d samples (removed %d)\n",
    n_before, nrow(meta), n_before - nrow(meta)))
}

# Gene annotations (dream_results.csv lacks annotations; use unified_disease_signatures)
annot_file <- file.path(INT, "results/disease_signatures/unified_disease_signatures.csv")
if (file.exists(annot_file)) {
  gene_annot <- fread(annot_file,
    select = c("gene", "symbol", "gene_type", "mouse_gene_id", "mouse_symbol"))
  gene_annot <- unique(gene_annot, by = "gene")
} else {
  gene_annot <- data.table(gene = character(), symbol = character(),
    gene_type = character(), mouse_gene_id = character(), mouse_symbol = character())
}

# --- Prepare sex covariate (same as Scripts 130/131) ---
meta[, sex_covar := inferred_sex]
na_sex <- is.na(meta$sex_covar) | meta$sex_covar == ""
if (any(na_sex) && "sex" %in% names(meta)) {
  meta$sex_covar[na_sex] <- meta$sex[na_sex]
}
meta[, sex_covar := factor(sex_covar)]

# --- Subset to NAFL + NASH + Borderline (same logic as Script 13) ---
meta_nn <- meta[diagnosis_harmonized %in% c("NAFL", "NASH", "Borderline")]
meta_nn[, nafl_nash := fifelse(diagnosis_harmonized == "NAFL", "NAFL", "NASH")]
meta_nn[, nafl_nash := factor(nafl_nash, levels = c("NAFL", "NASH"))]

cat(sprintf("Total NAFL+NASH+Borderline samples: %d\n", nrow(meta_nn)))
cat("  Binary grouping (NAFL vs NASH+Borderline):\n")
print(meta_nn[, .N, by = .(dataset, nafl_nash)][order(dataset, nafl_nash)])

# --- dataset_subbatch for random effect (PRJNA512027 L0/S0 subbatch logic
# permanently removed 2026-05-15 with the cohort itself) ---
meta_nn[, dataset_subbatch := as.character(dataset)]

# --- Drop samples with NA sex ---
na_sex_nn <- is.na(meta_nn$sex_covar) | as.character(meta_nn$sex_covar) == ""
if (any(na_sex_nn)) {
  cat(sprintf("Dropping %d samples with missing sex annotation\n", sum(na_sex_nn)))
  meta_nn <- meta_nn[!na_sex_nn]
}

cat("Sex x nafl_nash distribution:\n")
print(table(meta_nn$nafl_nash, meta_nn$sex_covar))
cat("Sex x dataset distribution:\n")
print(table(meta_nn$sex_covar, meta_nn$dataset))

# Detect sex labels from data (could be "M"/"F" or "Male"/"Female")
sex_levels <- levels(meta_nn$sex_covar)
cat("Sex factor levels:", paste(sex_levels, collapse = ", "), "\n")
male_label   <- sex_levels[grepl("^M", sex_levels)][1]
female_label <- sex_levels[grepl("^F", sex_levels)][1]
cat("Using male label:", male_label, ", female label:", female_label, "\n")

# ============================================================
# PART 1: Interaction model — nafl_nash * sex_covar
# ============================================================
cat(sprintf("\n%s\n", strrep("=", 60)))
cat("  PART 1: Interaction model (nafl_nash * sex_covar)\n")
cat(sprintf("%s\n", strrep("=", 60)))

# Filter datasets with >= 2 samples per nafl_nash group (pooled across sex)
ds_counts_int <- meta_nn[, .(
  n_nafl = sum(nafl_nash == "NAFL"),
  n_nash = sum(nafl_nash == "NASH")
), by = dataset_subbatch]
valid_ds_int <- ds_counts_int[n_nafl >= 2 & n_nash >= 2, dataset_subbatch]
meta_int <- meta_nn[dataset_subbatch %in% valid_ds_int]

cat(sprintf("Valid datasets (interaction): %d (%s)\n",
  length(valid_ds_int), paste(valid_ds_int, collapse = ", ")))

# Build DGEList
keep_samples_int <- intersect(meta_int$sample_id, colnames(counts))
meta_int <- meta_int[sample_id %in% keep_samples_int]
dge_int <- DGEList(counts = counts[, keep_samples_int])

m_ordered <- meta_int[match(colnames(dge_int), meta_int$sample_id)]
dge_int$samples <- cbind(dge_int$samples, m_ordered[, .(
  nafl_nash, dataset_subbatch, sex_covar)])
dge_int$samples$nafl_nash <- factor(dge_int$samples$nafl_nash, levels = c("NAFL", "NASH"))
dge_int$samples$dataset_subbatch <- factor(dge_int$samples$dataset_subbatch)
dge_int$samples$sex_covar <- factor(dge_int$samples$sex_covar)

dge_int <- calcNormFactors(dge_int, method = "TMM")
keep_genes_int <- filterByExpr(dge_int, group = dge_int$samples$nafl_nash)
dge_int <- dge_int[keep_genes_int, , keep.lib.sizes = FALSE]

cat(sprintf("  Samples: %d (NAFL=%d, NASH=%d)\n", ncol(dge_int),
  sum(dge_int$samples$nafl_nash == "NAFL"),
  sum(dge_int$samples$nafl_nash == "NASH")))
cat(sprintf("  Genes: %d\n", nrow(dge_int)))

# Determine formula — random vs fixed effect
has_random_int <- length(valid_ds_int) >= 2
if (has_random_int) {
  form_int <- ~ nafl_nash * sex_covar + (1 | dataset_subbatch)
} else {
  form_int <- ~ nafl_nash * sex_covar + dataset_subbatch
}
cat(sprintf("  Formula: %s\n", deparse(form_int)))

if (has_random_int) {
  vobj_int <- suppressWarnings(voomWithDreamWeights(dge_int, form_int,
    dge_int$samples, BPPARAM = BPPARAM))
  fit_int <- suppressWarnings(dream(vobj_int, form_int,
    dge_int$samples, BPPARAM = BPPARAM))
  # NOTE: do NOT call eBayes() after dream()
} else {
  design_int <- model.matrix(form_int, data = dge_int$samples)
  vobj_int <- voom(dge_int, design_int, plot = FALSE)
  fit_int <- lmFit(vobj_int, design_int)
  fit_int <- eBayes(fit_int)
}

# Find the interaction coefficient
all_coefs <- colnames(fit_int$coefficients)
cat("  Available coefficients:", paste(all_coefs, collapse = ", "), "\n")
int_coef <- grep("nafl_nash.*sex_covar|sex_covar.*nafl_nash", all_coefs, value = TRUE)
if (length(int_coef) == 0) {
  stop("No interaction coefficient found among: ", paste(all_coefs, collapse = ", "))
}
if (length(int_coef) > 1) int_coef <- int_coef[1]
cat("  Using interaction coefficient:", int_coef, "\n")

res_int <- topTable(fit_int, coef = int_coef, number = Inf, sort.by = "none")
res_int$gene <- rownames(res_int)
res_int_dt <- as.data.table(res_int)
res_int_dt <- merge(res_int_dt, gene_annot, by = "gene", all.x = TRUE)
setnames(res_int_dt, "adj.P.Val", "padj", skip_absent = TRUE)
res_int_dt[, contrast_id := "c10_interaction"]

n_sig_int <- sum(res_int_dt$padj < 0.1, na.rm = TRUE)
cat(sprintf("  Sex-differential DEGs (interaction padj<0.1): %d\n", n_sig_int))
cat(sprintf("    Stronger in males (logFC > 0): %d\n",
  sum(res_int_dt$padj < 0.1 & res_int_dt$logFC > 0, na.rm = TRUE)))
cat(sprintf("    Stronger in females (logFC < 0): %d\n",
  sum(res_int_dt$padj < 0.1 & res_int_dt$logFC < 0, na.rm = TRUE)))

fwrite(res_int_dt, file.path(ODIR, "c10_sex_interaction_nafl_nash.csv"))
cat("  Saved: c10_sex_interaction_nafl_nash.csv\n")

gc()

# ============================================================
# PART 2: Sex-stratified dream runs (female-only and male-only)
# ============================================================
cat(sprintf("\n%s\n", strrep("=", 60)))
cat("  PART 2: Sex-stratified NAFL vs NASH dream runs\n")
cat(sprintf("%s\n", strrep("=", 60)))

run_stratum_nafl_nash <- function(sex_label, sex_name, meta_all, counts_all,
                                   gene_annot, BPPARAM) {
  cat(sprintf("\n--- Stratum: %s (%s) ---\n", sex_name, sex_label))

  meta_s <- meta_all[sex_covar == sex_label]
  cat(sprintf("  Samples: %d\n", nrow(meta_s)))

  if (nrow(meta_s) < 5) {
    cat("  SKIPPED: Too few samples for meaningful analysis.\n")
    return(NULL)
  }

  # Check per-dataset distribution
  cat("  nafl_nash x dataset distribution:\n")
  print(table(meta_s$nafl_nash, meta_s$dataset_subbatch))

  # Require both groups present
  n_nafl <- sum(meta_s$nafl_nash == "NAFL")
  n_nash <- sum(meta_s$nafl_nash == "NASH")
  if (n_nafl == 0 || n_nash == 0) {
    cat(sprintf("  SKIPPED: Missing group (NAFL=%d, NASH=%d)\n", n_nafl, n_nash))
    return(NULL)
  }

  # Filter datasets with >= 2 samples in each group
  ds_counts <- meta_s[, .(
    n_nafl = sum(nafl_nash == "NAFL"),
    n_nash = sum(nafl_nash == "NASH")
  ), by = dataset_subbatch]
  valid_ds <- ds_counts[n_nafl >= 2 & n_nash >= 2, dataset_subbatch]

  if (length(valid_ds) == 0) {
    cat("  SKIPPED: No datasets with >= 2 samples per group.\n")
    return(NULL)
  }
  meta_s <- meta_s[dataset_subbatch %in% valid_ds]
  meta_s[, dataset_subbatch := factor(droplevels(factor(dataset_subbatch)))]

  cat(sprintf("  Valid datasets: %d (%s)\n", length(valid_ds),
    paste(valid_ds, collapse = ", ")))
  cat(sprintf("  Samples after filtering: %d (NAFL=%d, NASH=%d)\n",
    nrow(meta_s), sum(meta_s$nafl_nash == "NAFL"), sum(meta_s$nafl_nash == "NASH")))

  # Build DGEList
  keep_samples <- intersect(meta_s$sample_id, colnames(counts_all))
  meta_s <- meta_s[sample_id %in% keep_samples]
  dge <- DGEList(counts = counts_all[, keep_samples])

  m_ordered <- meta_s[match(colnames(dge), meta_s$sample_id)]
  dge$samples <- cbind(dge$samples, m_ordered[, .(nafl_nash, dataset_subbatch)])
  dge$samples$nafl_nash <- factor(dge$samples$nafl_nash, levels = c("NAFL", "NASH"))
  dge$samples$dataset_subbatch <- factor(dge$samples$dataset_subbatch)

  dge <- calcNormFactors(dge, method = "TMM")
  keep_genes <- filterByExpr(dge, group = dge$samples$nafl_nash)
  dge <- dge[keep_genes, , keep.lib.sizes = FALSE]

  cat(sprintf("  Genes after filterByExpr: %d\n", nrow(dge)))

  # Determine formula — random effect requires >= 2 datasets
  n_ds <- nlevels(dge$samples$dataset_subbatch)
  has_random <- n_ds >= 2

  if (has_random) {
    form_s <- ~ nafl_nash + (1 | dataset_subbatch)
    cat(sprintf("  Formula: ~ nafl_nash + (1 | dataset_subbatch) [%d datasets]\n", n_ds))
  } else {
    form_s <- ~ nafl_nash
    cat(sprintf("  Formula: ~ nafl_nash [single dataset, fixed effect only]\n"))
  }

  if (has_random) {
    vobj <- suppressWarnings(voomWithDreamWeights(dge, form_s, dge$samples, BPPARAM = BPPARAM))
    fit  <- suppressWarnings(dream(vobj, form_s, dge$samples, BPPARAM = BPPARAM))
    # NOTE: do NOT call eBayes() after dream()
  } else {
    design <- model.matrix(form_s, data = dge$samples)
    vobj <- voom(dge, design, plot = FALSE)
    fit  <- lmFit(vobj, design)
    fit  <- eBayes(fit)
  }

  res <- topTable(fit, coef = "nafl_nashNASH", number = Inf, sort.by = "none")
  res$gene <- rownames(res)
  res_dt <- as.data.table(res)
  res_dt <- merge(res_dt, gene_annot, by = "gene", all.x = TRUE)
  setnames(res_dt, "adj.P.Val", "padj", skip_absent = TRUE)
  res_dt[, contrast_id := paste0("c10_", tolower(sex_name))]
  res_dt[, n_samples := ncol(dge)]
  res_dt[, n_datasets := n_ds]

  n_sig <- sum(res_dt$padj < 0.1, na.rm = TRUE)
  n_up  <- sum(res_dt$padj < 0.1 & res_dt$logFC > 0, na.rm = TRUE)
  n_down <- sum(res_dt$padj < 0.1 & res_dt$logFC < 0, na.rm = TRUE)
  cat(sprintf("  DEGs (padj<0.1): %d (Up: %d, Down: %d)\n", n_sig, n_up, n_down))
  cat(sprintf("  DEGs (padj<0.05): %d\n", sum(res_dt$padj < 0.05, na.rm = TRUE)))

  # Save
  out_file <- file.path(ODIR, sprintf("c10_%s_nafl_nash.csv", tolower(sex_name)))
  fwrite(res_dt, out_file)
  cat(sprintf("  Saved: %s\n", basename(out_file)))

  gc()
  return(res_dt)
}

res_female <- run_stratum_nafl_nash(female_label, "female", meta_nn, counts,
  gene_annot, BPPARAM)
res_male   <- run_stratum_nafl_nash(male_label, "male", meta_nn, counts,
  gene_annot, BPPARAM)

# ============================================================
# PART 3: Gene classification
# ============================================================
cat(sprintf("\n%s\n", strrep("=", 60)))
cat("  PART 3: Sex-progression gene classification\n")
cat(sprintf("%s\n", strrep("=", 60)))

PADJ_THRESH    <- 0.1
LFC_DIFF_THRESH <- 0.5

if (!is.null(res_female) && !is.null(res_male)) {
  # Join female and male results on gene
  cls_dt <- merge(
    res_female[, .(gene, logFC_F = logFC, padj_F = padj)],
    res_male[,   .(gene, logFC_M = logFC, padj_M = padj)],
    by = "gene", all = TRUE
  )

  # Mark genes with significant interaction term
  int_sig_genes <- res_int_dt[padj < PADJ_THRESH, gene]

  cls_dt[, sig_F := (!is.na(padj_F) & padj_F < PADJ_THRESH)]
  cls_dt[, sig_M := (!is.na(padj_M) & padj_M < PADJ_THRESH)]
  cls_dt[, same_sign := !is.na(logFC_F) & !is.na(logFC_M) &
    sign(logFC_F) == sign(logFC_M)]
  cls_dt[, lfc_diff := abs(logFC_F - logFC_M)]

  cls_dt[, sex_progression_class := fcase(
    sig_F & sig_M & !same_sign,  "sex_divergent_progression",
    sig_F & sig_M & same_sign,   "shared_progression",
    sig_F & !sig_M,              "female_progression",
    !sig_F & sig_M,              "male_progression",
    default = "NS"
  )]

  # Flag sex_differential: interaction padj<0.1 AND |logFC_F - logFC_M| > 0.5
  cls_dt[, sex_differential := (gene %in% int_sig_genes) &
    !is.na(lfc_diff) & (lfc_diff > LFC_DIFF_THRESH)]

  # Add gene annotations
  cls_dt <- merge(cls_dt, gene_annot, by = "gene", all.x = TRUE)

  # Add interaction padj for reference
  cls_dt <- merge(cls_dt, res_int_dt[, .(gene, padj_interaction = padj,
    logFC_interaction = logFC)],
    by = "gene", all.x = TRUE)

  cat("Sex-progression class distribution:\n")
  print(cls_dt[, .N, by = sex_progression_class][order(-N)])
  cat(sprintf("Sex_differential genes (interaction + |LFC diff| > %.1f): %d\n",
    LFC_DIFF_THRESH, sum(cls_dt$sex_differential, na.rm = TRUE)))

  fwrite(cls_dt, file.path(ODIR, "c10_sex_progression_classification.csv"))
  cat("Saved: c10_sex_progression_classification.csv\n")

  # Top sex-differential genes
  sex_diff_genes <- cls_dt[sex_differential == TRUE][order(padj_interaction)]
  if (nrow(sex_diff_genes) > 0) {
    cat(sprintf("\nTop sex-differential genes (up to 10):\n"))
    top_sd <- head(sex_diff_genes, 10)
    for (i in seq_len(nrow(top_sd))) {
      cat(sprintf("  %s: LFC_F=%.3f, LFC_M=%.3f, diff=%.3f, int_padj=%.2e\n",
        top_sd$symbol[i], top_sd$logFC_F[i], top_sd$logFC_M[i],
        top_sd$lfc_diff[i], top_sd$padj_interaction[i]))
    }
  }
} else {
  cat("WARNING: One or both sex strata failed — classification skipped.\n")
  cls_dt <- NULL
}

# ============================================================
# PART 4: Summary table
# ============================================================
cat(sprintf("\n%s\n", strrep("=", 60)))
cat("  C10 SEX-STRATIFIED NAFL-vs-NASH SUMMARY\n")
cat(sprintf("%s\n", strrep("=", 60)))

summary_rows <- list()

# Interaction model summary
summary_rows[["Interaction"]] <- data.table(
  model = "Interaction",
  coefficient = int_coef,
  n_samples = ncol(dge_int),
  n_genes = nrow(res_int_dt),
  n_deg_01 = sum(res_int_dt$padj < 0.1, na.rm = TRUE),
  n_deg_05 = sum(res_int_dt$padj < 0.05, na.rm = TRUE),
  n_up = sum(res_int_dt$padj < 0.1 & res_int_dt$logFC > 0, na.rm = TRUE),
  n_down = sum(res_int_dt$padj < 0.1 & res_int_dt$logFC < 0, na.rm = TRUE)
)

# Female model summary
if (!is.null(res_female)) {
  summary_rows[["Female"]] <- data.table(
    model = "Female_only",
    coefficient = "nafl_nashNASH",
    n_samples = res_female$n_samples[1],
    n_genes = nrow(res_female),
    n_deg_01 = sum(res_female$padj < 0.1, na.rm = TRUE),
    n_deg_05 = sum(res_female$padj < 0.05, na.rm = TRUE),
    n_up = sum(res_female$padj < 0.1 & res_female$logFC > 0, na.rm = TRUE),
    n_down = sum(res_female$padj < 0.1 & res_female$logFC < 0, na.rm = TRUE)
  )
} else {
  summary_rows[["Female"]] <- data.table(
    model = "Female_only", coefficient = "nafl_nashNASH",
    n_samples = 0, n_genes = 0, n_deg_01 = 0, n_deg_05 = 0, n_up = 0, n_down = 0)
}

# Male model summary
if (!is.null(res_male)) {
  summary_rows[["Male"]] <- data.table(
    model = "Male_only",
    coefficient = "nafl_nashNASH",
    n_samples = res_male$n_samples[1],
    n_genes = nrow(res_male),
    n_deg_01 = sum(res_male$padj < 0.1, na.rm = TRUE),
    n_deg_05 = sum(res_male$padj < 0.05, na.rm = TRUE),
    n_up = sum(res_male$padj < 0.1 & res_male$logFC > 0, na.rm = TRUE),
    n_down = sum(res_male$padj < 0.1 & res_male$logFC < 0, na.rm = TRUE)
  )
} else {
  summary_rows[["Male"]] <- data.table(
    model = "Male_only", coefficient = "nafl_nashNASH",
    n_samples = 0, n_genes = 0, n_deg_01 = 0, n_deg_05 = 0, n_up = 0, n_down = 0)
}

# Classification summary
if (!is.null(cls_dt)) {
  class_dist <- cls_dt[, .N, by = sex_progression_class][order(-N)]
  n_sex_diff <- sum(cls_dt$sex_differential, na.rm = TRUE)
  summary_rows[["Classification"]] <- data.table(
    model = "Classification",
    coefficient = paste(class_dist$sex_progression_class, class_dist$N,
      sep = "=", collapse = "; "),
    n_samples = NA_integer_,
    n_genes = nrow(cls_dt),
    n_deg_01 = n_sex_diff,
    n_deg_05 = NA_integer_,
    n_up = sum(cls_dt$sex_progression_class == "female_progression", na.rm = TRUE),
    n_down = sum(cls_dt$sex_progression_class == "male_progression", na.rm = TRUE)
  )
}

summary_dt <- rbindlist(summary_rows, fill = TRUE)
print(summary_dt)

fwrite(summary_dt, file.path(ODIR, "c10_sex_progression_summary.csv"))
cat("Saved: c10_sex_progression_summary.csv\n")

# Correlation of F vs M logFC (for reporting)
if (!is.null(res_female) && !is.null(res_male)) {
  merged_lfc <- merge(
    res_female[, .(gene, logFC_F = logFC)],
    res_male[,   .(gene, logFC_M = logFC)],
    by = "gene"
  )
  rho <- cor(merged_lfc$logFC_F, merged_lfc$logFC_M, method = "spearman",
    use = "complete.obs")
  cat(sprintf("\nFemale vs Male logFC Spearman rho: %.3f (%d genes)\n",
    rho, nrow(merged_lfc)))
}

cat("\n=== Script 132 complete ===\n")
