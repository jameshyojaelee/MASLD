#!/usr/bin/env Rscript
# 05f_mash_vs_healthy_dream.R
# ---------------------------------------------------------------------------
# MASH-vs-Healthy pooled binary dream() mega-analysis.
#
# Eligible cohorts (control-bearing yaml mega cohorts): GSE126848, GSE130970,
# GSE135251, GSE162694.
#
# MASH_DEF env var selects the contrast definition:
#   - "borderline_grouped" (default): Control vs (NASH + Borderline)
#       writes results/disease_signatures/mash_vs_healthy_dream.csv
#   - "strict": Control vs NASH only (NAS >= 5)
#       writes results/disease_signatures/mash_vs_healthy_dream_strict.csv
#
# Structure mirrors 05_dream_mega_analysis.R (yaml + reformulas boilerplate +
# dream formula) with the sex_for_model machinery from 13_nafl_vs_nash_de.R
# for handling sex-annotation gaps.
# PRJNA512027 permanently removed from the pipeline 2026-05-15 (L0/S0 library
# batch perfectly confounded with disease severity).
# ---------------------------------------------------------------------------

set.seed(42)
Sys.setenv(R_PARALLEL_SEED = "42")

suppressPackageStartupMessages({
  library(reformulas)
  library(lme4)
  library(data.table)
  library(edgeR)
  library(yaml)
})

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
ODIR <- file.path(INT, "results/disease_signatures")
dir.create(ODIR, recursive = TRUE, showWarnings = FALSE)

mash_def <- Sys.getenv("MASH_DEF", "borderline_grouped")
if (!mash_def %in% c("borderline_grouped", "strict")) {
  stop("MASH_DEF must be 'borderline_grouped' or 'strict'; got: ", mash_def)
}
cat("MASH_DEF:", mash_def, "\n")

counts <- readRDS(file.path(RDIR, "merged_counts_raw.rds"))
dge_in <- readRDS(file.path(RDIR, "merged_dge.rds"))
meta   <- readRDS(file.path(RDIR, "meta_matched.rds"))
qc     <- fread(file.path(INT, "qc/sample_qc_report.csv"))
meta   <- meta[sample_id %in% qc[pass_technical == TRUE, sample_id]]

if (mash_def == "borderline_grouped") {
  meta_mh <- meta[diagnosis_harmonized %in% c("Control", "NASH", "Borderline")]
} else {
  meta_mh <- meta[diagnosis_harmonized %in% c("Control", "NASH")]
}
# Defensive PRJNA512027 filter (cohort removed from yaml 2026-05-15; guards
# against stale meta_matched.rds inputs).
meta_mh <- meta_mh[dataset != "PRJNA512027"]
meta_mh[, mash_status := fifelse(diagnosis_harmonized == "Control", "Control", "MASH")]
meta_mh[, mash_status := factor(mash_status, levels = c("Control", "MASH"))]

cat("Pre-cohort-filter sample counts:\n")
print(meta_mh[, .N, by = .(dataset, mash_status)][order(dataset, mash_status)])

cohort_counts <- meta_mh[, .(
  n_control = sum(mash_status == "Control"),
  n_mash = sum(mash_status == "MASH")
), by = dataset]
eligible_cohorts <- cohort_counts[n_control >= 3 & n_mash >= 3, dataset]
cat("\nEligible cohorts (>= 3 in BOTH arms):", paste(sort(eligible_cohorts), collapse = ", "), "\n")
dropped <- setdiff(cohort_counts$dataset, eligible_cohorts)
if (length(dropped) > 0) {
  cat("Dropped cohorts:", paste(sort(dropped), collapse = ", "), "\n")
}

meta_mh <- meta_mh[dataset %in% eligible_cohorts]
meta_mh <- meta_mh[order(sample_id)]
cat("Final samples:", nrow(meta_mh), "\n")

idx_all <- colnames(counts) %in% meta_mh$sample_id
dge_all <- DGEList(counts = counts[, idx_all])
dge_all$samples <- cbind(dge_all$samples,
  meta_mh[match(colnames(dge_all), meta_mh$sample_id),
    .(mash_status, dataset, sex, inferred_sex, age)])

# Dataset is the random-effect grouping (PRJNA512027 L0/S0 subbatch logic
# permanently removed 2026-05-15 with the cohort itself).
dge_all$samples$dataset_subbatch <- as.character(dge_all$samples$dataset)

dge_all$samples$sex_for_model <- dge_all$samples$sex
na_sex <- is.na(dge_all$samples$sex_for_model) | dge_all$samples$sex_for_model == ""
dge_all$samples$sex_for_model[na_sex] <- as.character(dge_all$samples$inferred_sex[na_sex])
dge_all$samples$sex_for_model <- factor(dge_all$samples$sex_for_model)

dge_all <- calcNormFactors(dge_all, method = "TMM")
keep_all <- filterByExpr(dge_all, group = dge_all$samples$mash_status)
dge_all <- dge_all[keep_all, , keep.lib.sizes = FALSE]

cat("Dream samples:", ncol(dge_all), "\n")
cat("Dream genes:", nrow(dge_all), "\n")
cat("Batch levels:", paste(sort(unique(dge_all$samples$dataset_subbatch)), collapse = ", "), "\n")
cat("Sex levels:", paste(levels(dge_all$samples$sex_for_model), collapse = ", "), "\n")
cat("Group distribution:\n")
print(table(dge_all$samples$mash_status, dge_all$samples$dataset))

form <- ~ mash_status + sex_for_model + (1 | dataset_subbatch)
cat("Formula:", deparse(form), "\n")

n_cores <- min(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4")), 32)
cat("Using", n_cores, "CPU cores\n")
BPPARAM <- MulticoreParam(n_cores, progressbar = TRUE)

vobjDream <- voomWithDreamWeights(dge_all, form, dge_all$samples, BPPARAM = BPPARAM)
cat("Running dream...\n")
fitDream <- dream(vobjDream, form, dge_all$samples, BPPARAM = BPPARAM)

dream_tt <- topTable(fitDream, coef = "mash_statusMASH", number = Inf, sort.by = "none")
dream_tt$gene <- rownames(dream_tt)

sig <- sum(dream_tt$adj.P.Val < 0.05)
sig_up <- sum(dream_tt$adj.P.Val < 0.05 & dream_tt$logFC > 0)
sig_down <- sum(dream_tt$adj.P.Val < 0.05 & dream_tt$logFC < 0)
sig_lfc <- sum(dream_tt$adj.P.Val < 0.05 & abs(dream_tt$logFC) > 0.5)
cat("\n===== MASH vs Healthy DREAM RESULTS =====\n")
cat("Total genes tested:", nrow(dream_tt), "\n")
cat(sprintf("DEGs (padj<0.05): %d (Up: %d, Down: %d)\n", sig, sig_up, sig_down))
cat(sprintf("DEGs (padj<0.05, |logFC|>0.5): %d\n", sig_lfc))

out_fname <- if (mash_def == "borderline_grouped") {
  "mash_vs_healthy_dream.csv"
} else {
  "mash_vs_healthy_dream_strict.csv"
}
fwrite(as.data.table(dream_tt), file.path(ODIR, out_fname))
cat("\nSaved:", file.path(ODIR, out_fname), "\n")
