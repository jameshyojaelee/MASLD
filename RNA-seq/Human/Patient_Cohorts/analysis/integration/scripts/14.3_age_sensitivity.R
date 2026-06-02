#!/usr/bin/env Rscript
# 14.3_age_sensitivity.R
# Age Confounding Sensitivity Analysis
#
# Problem: Dream formula does not adjust for age. Disease samples are ~5 years
# older than Controls. This script tests whether age confounding materially
# affects dream DEG results.
#
# Strategy: Matched-subset comparison — run dream WITH and WITHOUT age on the
# SAME subset of age-available datasets, isolating the age adjustment effect
# from sample size effects.
#
# Three analyses:
#   A) Age-Disease Diagnostic: Wilcoxon test, Cohen's d, per-dataset correlation
#   B) Matched-Subset Dream Comparison: same 4 datasets, ± age covariate
#   C) Age-Correlated Gene Identification: age-associated genes vs dream DEGs
#
# Usage: Rscript 14.3_age_sensitivity.R
# SLURM: cpu partition, 16 CPUs, 64GB RAM, ~4h

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

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR <- file.path(INT, "results/integration")
OUTDIR <- file.path(BASE, "RNA-seq/results/audit_sensitivity/age_confounding")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

cat("=== 14.3: Age Confounding Sensitivity Analysis ===\n")
cat("Started:", as.character(Sys.time()), "\n\n")

# --- Setup parallel backend ---
ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
cat("Using", ncpus, "CPU cores\n")
param <- if (ncpus > 1) MulticoreParam(ncpus) else SerialParam()

# --- Load data ---
cat("Loading merged DGE...\n")
dge <- readRDS(file.path(RDIR, "merged_dge.rds"))
cat("  DGE dimensions:", nrow(dge), "genes x", ncol(dge), "samples\n")

cat("Loading matched metadata...\n")
meta_matched <- readRDS(file.path(RDIR, "meta_matched.rds"))

cat("Loading unified metadata...\n")
meta_unified <- fread(file.path(INT, "metadata/unified_metadata.csv"))

cat("Loading QC report...\n")
qc <- fread(file.path(INT, "qc/sample_qc_report.csv"))
pass_samples <- qc[pass_technical == TRUE, sample_id]
cat("  QC-passing samples:", length(pass_samples), "\n")

# Load primary dream results for comparison
cat("Loading primary dream results...\n")
dream_file <- file.path(RDIR, "dream_results.csv")
if (!file.exists(dream_file)) stop("dream_results.csv not found — run 05 first")
dream_primary <- fread(dream_file)
if ("adj.P.Val" %in% names(dream_primary) && !"padj" %in% names(dream_primary))
  setnames(dream_primary, "adj.P.Val", "padj")
cat("  Primary dream DEGs (padj<0.1):", sum(dream_primary$padj < 0.1, na.rm = TRUE), "\n\n")

# ============================================================
# ANALYSIS A: Age-Disease Diagnostic
# ============================================================
cat("=== ANALYSIS A: Age-Disease Diagnostic ===\n")

# All samples with age data (including permanently excluded datasets)
meta_with_age <- meta_unified[
  !is.na(age) & age != "" &
  sample_id %in% pass_samples
]
meta_with_age[, age := as.numeric(age)]
meta_with_age <- meta_with_age[!is.na(age)]

cat("Samples with age + QC-passing:", nrow(meta_with_age), "\n")
cat("Datasets with age:", paste(unique(meta_with_age$dataset), collapse = ", "), "\n")

# Overall age difference
ctrl_ages <- meta_with_age[group_binary == "Control", age]
dis_ages  <- meta_with_age[group_binary == "Disease", age]
cat("\nAge summary:\n")
cat("  Control: N=", length(ctrl_ages), ", mean=", round(mean(ctrl_ages), 1),
    " ± ", round(sd(ctrl_ages), 1), "\n")
cat("  Disease: N=", length(dis_ages), ", mean=", round(mean(dis_ages), 1),
    " ± ", round(sd(dis_ages), 1), "\n")

# Wilcoxon test
wilcox_res <- wilcox.test(ctrl_ages, dis_ages)
cat("  Wilcoxon p =", signif(wilcox_res$p.value, 4), "\n")

# Cohen's d
pooled_sd <- sqrt(((length(ctrl_ages) - 1) * sd(ctrl_ages)^2 +
                    (length(dis_ages) - 1) * sd(dis_ages)^2) /
                   (length(ctrl_ages) + length(dis_ages) - 2))
cohens_d <- (mean(dis_ages) - mean(ctrl_ages)) / pooled_sd
cat("  Cohen's d =", round(cohens_d, 3), "\n")

# Per-dataset age × disease
per_dataset_age <- meta_with_age[, {
  if (length(unique(group_binary)) >= 2 && .N >= 5) {
    ctrl <- age[group_binary == "Control"]
    dis  <- age[group_binary == "Disease"]
    if (length(ctrl) >= 2 && length(dis) >= 2) {
      wt <- wilcox.test(ctrl, dis)
      d_val <- (mean(dis) - mean(ctrl)) / sqrt(((length(ctrl)-1)*sd(ctrl)^2 + (length(dis)-1)*sd(dis)^2) / (length(ctrl)+length(dis)-2))
      list(
        n_control = length(ctrl),
        n_disease = length(dis),
        mean_age_control = round(mean(ctrl), 1),
        mean_age_disease = round(mean(dis), 1),
        age_diff = round(mean(dis) - mean(ctrl), 1),
        wilcox_p = signif(wt$p.value, 4),
        cohens_d = round(d_val, 3)
      )
    } else {
      list(n_control = length(ctrl), n_disease = length(dis),
           mean_age_control = NA_real_, mean_age_disease = NA_real_,
           age_diff = NA_real_, wilcox_p = NA_real_, cohens_d = NA_real_)
    }
  } else {
    list(n_control = sum(group_binary == "Control"),
         n_disease = sum(group_binary == "Disease"),
         mean_age_control = NA_real_, mean_age_disease = NA_real_,
         age_diff = NA_real_, wilcox_p = NA_real_, cohens_d = NA_real_)
  }
}, by = dataset]

cat("\nPer-dataset age diagnostics:\n")
print(per_dataset_age)

# Summary output
diag_summary <- data.table(
  metric = c("n_samples_with_age", "n_control", "n_disease",
             "mean_age_control", "sd_age_control",
             "mean_age_disease", "sd_age_disease",
             "age_difference", "wilcoxon_p", "cohens_d"),
  value = c(nrow(meta_with_age), length(ctrl_ages), length(dis_ages),
            round(mean(ctrl_ages), 1), round(sd(ctrl_ages), 1),
            round(mean(dis_ages), 1), round(sd(dis_ages), 1),
            round(mean(dis_ages) - mean(ctrl_ages), 1),
            signif(wilcox_res$p.value, 4), round(cohens_d, 3))
)

fwrite(diag_summary, file.path(OUTDIR, "age_disease_diagnostic.csv"))
fwrite(per_dataset_age, file.path(OUTDIR, "age_disease_per_dataset.csv"))
cat("Diagnostic saved to", file.path(OUTDIR, "age_disease_diagnostic.csv"), "\n\n")

# ============================================================
# ANALYSIS B: Matched-Subset Dream Comparison
# ============================================================
cat("=== ANALYSIS B: Matched-Subset Dream Comparison ===\n")

# Age-annotated cohorts intersected with config/human_datasets.yaml include_in_mega
# (strict consistency with primary 5-cohort mega-analysis: drops GSE174478/193066
# which have no controls — sample size for age sensitivity reduces accordingly).
ycfg <- yaml::read_yaml(file.path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design", "config/human_datasets.yaml"))$datasets
mega_cohorts <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))
EXCL_PERMANENT <- setdiff(names(ycfg), mega_cohorts)
AGE_DREAM_DATASETS <- intersect(c("GSE130970", "GSE162694", "GSE174478", "GSE193066"), mega_cohorts)
cat("Age-eligible mega cohorts (k =", length(AGE_DREAM_DATASETS), "):", paste(AGE_DREAM_DATASETS, collapse=", "), "\n")

# Subset to age-available dream-eligible samples
meta_subset <- meta_unified[
  dataset %in% AGE_DREAM_DATASETS &
  !is.na(age) & age != "" &
  sample_id %in% pass_samples
]
meta_subset[, age := as.numeric(age)]
meta_subset <- meta_subset[!is.na(age)]

subset_samples <- intersect(meta_subset$sample_id, colnames(dge))
cat("Matched subset samples:", length(subset_samples), "\n")
cat("Dataset breakdown:\n")
print(meta_subset[sample_id %in% subset_samples, .N, by = .(dataset, group_binary)])

# Subset DGE
dge_sub <- dge[, subset_samples]

# Match metadata
meta_sub_matched <- meta_subset[match(colnames(dge_sub), meta_subset$sample_id)]
sex_sub <- meta_matched$inferred_sex[match(colnames(dge_sub), meta_matched$sample_id)]

# Build info dataframe (shared between both dream runs)
info_sub <- data.frame(
  group_binary = factor(meta_sub_matched$group_binary, levels = c("Control", "Disease")),
  dataset      = factor(meta_sub_matched$dataset),
  inferred_sex = factor(sex_sub),
  age          = as.numeric(meta_sub_matched$age),
  row.names    = colnames(dge_sub),
  stringsAsFactors = FALSE
)

# Scale age for numerical stability
info_sub$age_scaled <- scale(info_sub$age)[, 1]
cat("Age range:", range(info_sub$age), "\n")
cat("Scaled age range:", round(range(info_sub$age_scaled), 2), "\n")

# Check dataset has multiple levels for random effect
n_datasets <- nlevels(droplevels(info_sub$dataset))
cat("Number of datasets in subset:", n_datasets, "\n")

# Check for disease-only datasets (no controls)
ds_groups <- meta_sub_matched[, .(
  n_control = sum(group_binary == "Control"),
  n_disease = sum(group_binary == "Disease")
), by = dataset]
cat("Dataset group counts:\n")
print(ds_groups)

# Filter by expression
keep_sub <- filterByExpr(dge_sub, group = info_sub$group_binary)
dge_sub <- dge_sub[keep_sub, , keep.lib.sizes = FALSE]
dge_sub <- calcNormFactors(dge_sub, method = "RLE")
cat("After expression filter:", nrow(dge_sub), "genes\n")

# --- Dream Run 1: WITHOUT age (matches primary formula) ---
cat("\n--- Dream Run 1: WITHOUT age ---\n")
if (n_datasets > 1) {
  form_noage <- ~ group_binary + inferred_sex + (1|dataset)
} else {
  form_noage <- ~ group_binary + inferred_sex
}
cat("Formula:", deparse(form_noage), "\n")

v_noage <- suppressWarnings(voomWithDreamWeights(dge_sub, form_noage, info_sub, BPPARAM = param))
fit_noage <- suppressWarnings(dream(v_noage, form_noage, info_sub, BPPARAM = param))
res_noage <- topTable(fit_noage, coef = "group_binaryDisease", number = Inf, sort.by = "none")
res_noage$gene <- rownames(res_noage)
res_noage_dt <- as.data.table(res_noage)
if ("adj.P.Val" %in% names(res_noage_dt)) setnames(res_noage_dt, "adj.P.Val", "padj")

n_deg_noage <- sum(res_noage_dt$padj < 0.1, na.rm = TRUE)
cat("DEGs (padj<0.1) without age:", n_deg_noage, "\n")

# --- Dream Run 2: WITH age ---
cat("\n--- Dream Run 2: WITH age ---\n")
if (n_datasets > 1) {
  form_age <- ~ group_binary + inferred_sex + age_scaled + (1|dataset)
} else {
  form_age <- ~ group_binary + inferred_sex + age_scaled
}
cat("Formula:", deparse(form_age), "\n")

v_age <- suppressWarnings(voomWithDreamWeights(dge_sub, form_age, info_sub, BPPARAM = param))
fit_age <- suppressWarnings(dream(v_age, form_age, info_sub, BPPARAM = param))
res_age <- topTable(fit_age, coef = "group_binaryDisease", number = Inf, sort.by = "none")
res_age$gene <- rownames(res_age)
res_age_dt <- as.data.table(res_age)
if ("adj.P.Val" %in% names(res_age_dt)) setnames(res_age_dt, "adj.P.Val", "padj")

n_deg_age <- sum(res_age_dt$padj < 0.1, na.rm = TRUE)
cat("DEGs (padj<0.1) with age:", n_deg_age, "\n")

# --- Comparison metrics ---
cat("\n--- Concordance Metrics ---\n")

# Align genes
common_genes <- intersect(res_noage_dt$gene, res_age_dt$gene)
cat("Common genes:", length(common_genes), "\n")

noage_aligned <- res_noage_dt[match(common_genes, gene)]
age_aligned   <- res_age_dt[match(common_genes, gene)]

# Spearman correlation of logFC
rho_lfc <- cor(noage_aligned$logFC, age_aligned$logFC, method = "spearman", use = "complete.obs")
r_lfc   <- cor(noage_aligned$logFC, age_aligned$logFC, method = "pearson", use = "complete.obs")
cat("Spearman rho (logFC):", round(rho_lfc, 4), "\n")
cat("Pearson r (logFC):", round(r_lfc, 4), "\n")

# Spearman correlation of t-statistics
rho_t <- cor(noage_aligned$t, age_aligned$t, method = "spearman", use = "complete.obs")
cat("Spearman rho (t-stat):", round(rho_t, 4), "\n")

# Direction concordance
dir_noage <- sign(noage_aligned$logFC)
dir_age   <- sign(age_aligned$logFC)
dir_concord <- mean(dir_noage == dir_age, na.rm = TRUE) * 100
cat("Direction concordance:", round(dir_concord, 1), "%\n")

# Jaccard similarity of DEG sets (padj < 0.1)
degs_noage <- noage_aligned[padj < 0.1, gene]
degs_age   <- age_aligned[padj < 0.1, gene]
jaccard <- length(intersect(degs_noage, degs_age)) /
           length(union(degs_noage, degs_age))
cat("Jaccard similarity (padj<0.1):", round(jaccard, 4), "\n")

# Genes that change status
gained <- setdiff(degs_age, degs_noage)
lost   <- setdiff(degs_noage, degs_age)
cat("Genes gained with age adjustment:", length(gained), "\n")
cat("Genes lost with age adjustment:", length(lost), "\n")

# Per-gene LFC shift
lfc_shift <- age_aligned$logFC - noage_aligned$logFC
cat("Mean |LFC shift|:", round(mean(abs(lfc_shift), na.rm = TRUE), 4), "\n")
cat("Median |LFC shift|:", round(median(abs(lfc_shift), na.rm = TRUE), 4), "\n")
cat("Max |LFC shift|:", round(max(abs(lfc_shift), na.rm = TRUE), 4), "\n")

# Build per-gene comparison table
comparison_dt <- data.table(
  gene = common_genes,
  logFC_noage = noage_aligned$logFC,
  logFC_age = age_aligned$logFC,
  logFC_shift = lfc_shift,
  padj_noage = noage_aligned$padj,
  padj_age = age_aligned$padj,
  t_noage = noage_aligned$t,
  t_age = age_aligned$t,
  sig_noage = noage_aligned$padj < 0.1,
  sig_age = age_aligned$padj < 0.1,
  status_change = ifelse(
    noage_aligned$padj < 0.1 & age_aligned$padj < 0.1, "stable_sig",
    ifelse(noage_aligned$padj >= 0.1 & age_aligned$padj >= 0.1, "stable_ns",
    ifelse(noage_aligned$padj < 0.1 & age_aligned$padj >= 0.1, "lost_with_age",
    "gained_with_age")))
)
comparison_dt <- comparison_dt[order(abs(logFC_shift), decreasing = TRUE)]

fwrite(comparison_dt, file.path(OUTDIR, "age_sensitivity_comparison.csv"))

# Status change summary
status_summary <- comparison_dt[, .N, by = status_change]
cat("\nStatus change summary:\n")
print(status_summary)

# Summary metrics
metrics_dt <- data.table(
  metric = c("n_subset_samples", "n_datasets", "n_genes_tested",
             "n_deg_noage", "n_deg_age",
             "spearman_rho_logFC", "pearson_r_logFC", "spearman_rho_tstat",
             "direction_concordance_pct", "jaccard_padj01",
             "n_gained_with_age", "n_lost_with_age",
             "mean_abs_lfc_shift", "median_abs_lfc_shift", "max_abs_lfc_shift"),
  value = c(length(subset_samples), n_datasets, length(common_genes),
            n_deg_noage, n_deg_age,
            round(rho_lfc, 4), round(r_lfc, 4), round(rho_t, 4),
            round(dir_concord, 1), round(jaccard, 4),
            length(gained), length(lost),
            round(mean(abs(lfc_shift), na.rm = TRUE), 4),
            round(median(abs(lfc_shift), na.rm = TRUE), 4),
            round(max(abs(lfc_shift), na.rm = TRUE), 4))
)
fwrite(metrics_dt, file.path(OUTDIR, "age_sensitivity_metrics.csv"))
cat("\nMetrics saved to", file.path(OUTDIR, "age_sensitivity_metrics.csv"), "\n")

# Save individual dream results for reference
fwrite(res_noage_dt, file.path(OUTDIR, "dream_subset_noage.csv"))
fwrite(res_age_dt, file.path(OUTDIR, "dream_subset_withage.csv"))

# ============================================================
# ANALYSIS C: Age-Correlated Gene Identification
# ============================================================
cat("\n=== ANALYSIS C: Age-Correlated Gene Identification ===\n")

# Use the age-available subset DGE (already filtered and normalized)
# Fit per-gene linear model: expression ~ age + dataset + sex
cat("Computing log2(CPM) for age correlation...\n")
lcpm <- cpm(dge_sub, log = TRUE)

# Prepare covariates
age_vec     <- info_sub$age_scaled
dataset_vec <- info_sub$dataset
sex_vec     <- info_sub$inferred_sex

# Fit per-gene lm
cat("Fitting per-gene linear models (", nrow(lcpm), " genes)...\n")
age_results <- rbindlist(lapply(seq_len(nrow(lcpm)), function(i) {
  expr <- lcpm[i, ]
  df <- data.frame(expr = expr, age = age_vec, dataset = dataset_vec, sex = sex_vec)

  # Handle single-dataset case
  if (n_datasets > 1) {
    fit <- tryCatch(lm(expr ~ age + dataset + sex, data = df), error = function(e) NULL)
  } else {
    fit <- tryCatch(lm(expr ~ age + sex, data = df), error = function(e) NULL)
  }

  if (is.null(fit)) {
    return(data.table(gene = rownames(lcpm)[i], age_coef = NA_real_,
                      age_se = NA_real_, age_t = NA_real_, age_p = NA_real_))
  }

  s <- summary(fit)$coefficients
  if ("age" %in% rownames(s)) {
    data.table(
      gene    = rownames(lcpm)[i],
      age_coef = s["age", "Estimate"],
      age_se   = s["age", "Std. Error"],
      age_t    = s["age", "t value"],
      age_p    = s["age", "Pr(>|t|)"]
    )
  } else {
    data.table(gene = rownames(lcpm)[i], age_coef = NA_real_,
               age_se = NA_real_, age_t = NA_real_, age_p = NA_real_)
  }
}))

# BH correction
age_results[, age_padj := p.adjust(age_p, method = "BH")]

# Summary
n_age_genes <- sum(age_results$age_padj < 0.05, na.rm = TRUE)
cat("Age-associated genes (padj<0.05):", n_age_genes, "\n")
n_age_genes_01 <- sum(age_results$age_padj < 0.1, na.rm = TRUE)
cat("Age-associated genes (padj<0.1):", n_age_genes_01, "\n")

# Cross-reference with primary dream DEGs
dream_degs <- dream_primary[padj < 0.1, gene]
age_genes  <- age_results[age_padj < 0.05, gene]
age_genes_01 <- age_results[age_padj < 0.1, gene]

overlap_005 <- length(intersect(dream_degs, age_genes))
overlap_01  <- length(intersect(dream_degs, age_genes_01))

cat("\nCross-reference with primary dream DEGs (padj<0.1):\n")
cat("  Dream DEGs:", length(dream_degs), "\n")
cat("  Age genes (padj<0.05):", length(age_genes), "\n")
cat("  Overlap:", overlap_005, "(", round(overlap_005 / length(dream_degs) * 100, 1), "% of dream DEGs)\n")
cat("  Age genes (padj<0.1):", length(age_genes_01), "\n")
cat("  Overlap:", overlap_01, "(", round(overlap_01 / length(dream_degs) * 100, 1), "% of dream DEGs)\n")

# Fisher's exact test for enrichment
n_total <- nrow(age_results)
n_dream <- sum(age_results$gene %in% dream_degs)
n_age   <- sum(age_results$gene %in% age_genes)
n_both  <- sum(age_results$gene %in% intersect(dream_degs, age_genes))

fisher_mat <- matrix(c(
  n_both, n_dream - n_both,
  n_age - n_both, n_total - n_dream - n_age + n_both
), nrow = 2)
fisher_res <- fisher.test(fisher_mat)
cat("  Fisher's exact test (enrichment): OR=", round(fisher_res$estimate, 2),
    ", p=", signif(fisher_res$p.value, 4), "\n")

# Add dream overlap flag
age_results[, is_dream_deg := gene %in% dream_degs]

fwrite(age_results, file.path(OUTDIR, "age_associated_genes.csv"))
cat("Age-associated genes saved to", file.path(OUTDIR, "age_associated_genes.csv"), "\n")

# Cross-reference summary
xref_summary <- data.table(
  metric = c("n_dream_degs_padj01", "n_age_genes_padj005", "n_age_genes_padj01",
             "n_overlap_padj005", "pct_dream_degs_age_005",
             "n_overlap_padj01", "pct_dream_degs_age_01",
             "fisher_OR", "fisher_p"),
  value = c(length(dream_degs), length(age_genes), length(age_genes_01),
            overlap_005, round(overlap_005 / length(dream_degs) * 100, 1),
            overlap_01, round(overlap_01 / length(dream_degs) * 100, 1),
            round(fisher_res$estimate, 3), signif(fisher_res$p.value, 4))
)
fwrite(xref_summary, file.path(OUTDIR, "age_dream_crossref_summary.csv"))

# ============================================================
# SUMMARY
# ============================================================
cat("\n")
cat("============================================================\n")
cat("SUMMARY: Age Confounding Sensitivity Analysis\n")
cat("============================================================\n")
cat("\nA) Age-Disease Diagnostic:\n")
cat("   Disease vs Control age diff: ", round(mean(dis_ages) - mean(ctrl_ages), 1),
    " years (p=", signif(wilcox_res$p.value, 4), ", d=", round(cohens_d, 3), ")\n")
cat("\nB) Matched-Subset Dream Comparison (", length(subset_samples), " samples, ",
    n_datasets, " datasets):\n")
cat("   Spearman rho (logFC): ", round(rho_lfc, 4), "\n")
cat("   Pearson r (logFC): ", round(r_lfc, 4), "\n")
cat("   Direction concordance: ", round(dir_concord, 1), "%\n")
cat("   Jaccard (padj<0.1): ", round(jaccard, 4), "\n")
cat("   DEGs without age: ", n_deg_noage, " | DEGs with age: ", n_deg_age, "\n")
cat("   Gained: ", length(gained), " | Lost: ", length(lost), "\n")
cat("   Mean |LFC shift|: ", round(mean(abs(lfc_shift), na.rm = TRUE), 4), "\n")
cat("\nC) Age-Correlated Genes:\n")
cat("   Age-associated (padj<0.05): ", n_age_genes, "\n")
cat("   Overlap with dream DEGs: ", overlap_005, " (",
    round(overlap_005 / length(dream_degs) * 100, 1), "%)\n")
cat("   Fisher OR: ", round(fisher_res$estimate, 2), " (p=",
    signif(fisher_res$p.value, 4), ")\n")

# Interpretation
cat("\nINTERPRETATION:\n")
if (rho_lfc > 0.95 && jaccard > 0.80) {
  cat("  NEGLIGIBLE age confounding (rho>0.95, Jaccard>0.80)\n")
  cat("  Age does not substantially affect dream DEG results.\n")
} else if (rho_lfc > 0.85) {
  cat("  MODERATE age effect (rho 0.85-0.95)\n")
  cat("  Some genes affected — check status_change column in comparison CSV.\n")
} else {
  cat("  SUBSTANTIAL age confounding (rho<0.85)\n")
  cat("  Consider age imputation or re-framing results.\n")
}

cat("\nOutput files:\n")
cat("  ", file.path(OUTDIR, "age_disease_diagnostic.csv"), "\n")
cat("  ", file.path(OUTDIR, "age_disease_per_dataset.csv"), "\n")
cat("  ", file.path(OUTDIR, "age_sensitivity_comparison.csv"), "\n")
cat("  ", file.path(OUTDIR, "age_sensitivity_metrics.csv"), "\n")
cat("  ", file.path(OUTDIR, "dream_subset_noage.csv"), "\n")
cat("  ", file.path(OUTDIR, "dream_subset_withage.csv"), "\n")
cat("  ", file.path(OUTDIR, "age_associated_genes.csv"), "\n")
cat("  ", file.path(OUTDIR, "age_dream_crossref_summary.csv"), "\n")

cat("\n=== 14.3 completed:", as.character(Sys.time()), "===\n")
