#!/usr/bin/env Rscript
# 05g_masl_vs_healthy_dream.R
# ---------------------------------------------------------------------------
# MASL-vs-Healthy pooled binary dream() mega-analysis.
#
# MASL = NAFL (simple steatosis, NAS<3 or original NAFL label).
# No Borderline ambiguity — diagnosis_harmonized == "NAFL" is the contrast arm.
#
# Eligible cohorts (control-bearing yaml mega cohorts): GSE126848, GSE130970,
# GSE135251, GSE162694.
#
# PRJNA512027 permanently removed from the pipeline 2026-05-15 (L0/S0 library
# batch perfectly confounded with disease severity).
#
# Output: results/disease_signatures/masl_vs_healthy_dream.csv
#
# Structure mirrors 05f_mash_vs_healthy_dream.R.
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

counts <- readRDS(file.path(RDIR, "merged_counts_raw.rds"))
meta   <- readRDS(file.path(RDIR, "meta_matched.rds"))
qc     <- fread(file.path(INT, "qc/sample_qc_report.csv"))
meta   <- meta[sample_id %in% qc[pass_technical == TRUE, sample_id]]

meta_ml <- meta[diagnosis_harmonized %in% c("Control", "NAFL")]
# Defensive PRJNA512027 filter (cohort removed from yaml 2026-05-15; guards
# against stale meta_matched.rds inputs).
meta_ml <- meta_ml[dataset != "PRJNA512027"]
meta_ml[, masl_status := factor(fifelse(diagnosis_harmonized == "Control", "Control", "MASL"),
                                levels = c("Control", "MASL"))]

cat("Pre-cohort-filter sample counts:\n")
print(meta_ml[, .N, by = .(dataset, masl_status)][order(dataset, masl_status)])

cohort_counts <- meta_ml[, .(
  n_control = sum(masl_status == "Control"),
  n_masl    = sum(masl_status == "MASL")
), by = dataset]
eligible_cohorts <- cohort_counts[n_control >= 3 & n_masl >= 3, dataset]
cat("\nEligible cohorts (>= 3 in BOTH arms):", paste(sort(eligible_cohorts), collapse = ", "), "\n")
dropped <- setdiff(cohort_counts$dataset, eligible_cohorts)
if (length(dropped) > 0) {
  cat("Dropped cohorts:", paste(sort(dropped), collapse = ", "), "\n")
}

meta_ml <- meta_ml[dataset %in% eligible_cohorts]
meta_ml <- meta_ml[order(sample_id)]
cat("Final samples:", nrow(meta_ml), "\n")

idx_all <- colnames(counts) %in% meta_ml$sample_id
dge_all <- DGEList(counts = counts[, idx_all])
dge_all$samples <- cbind(dge_all$samples,
  meta_ml[match(colnames(dge_all), meta_ml$sample_id),
    .(masl_status, dataset, sex, inferred_sex, age)])

# Dataset is the random-effect grouping (PRJNA512027 L0/S0 subbatch logic
# permanently removed 2026-05-15 with the cohort itself).
dge_all$samples$dataset_subbatch <- as.character(dge_all$samples$dataset)

dge_all$samples$sex_for_model <- dge_all$samples$sex
na_sex <- is.na(dge_all$samples$sex_for_model) | dge_all$samples$sex_for_model == ""
dge_all$samples$sex_for_model[na_sex] <- as.character(dge_all$samples$inferred_sex[na_sex])
dge_all$samples$sex_for_model <- factor(dge_all$samples$sex_for_model)

dge_all <- calcNormFactors(dge_all, method = "TMM")
keep_all <- filterByExpr(dge_all, group = dge_all$samples$masl_status)
dge_all <- dge_all[keep_all, , keep.lib.sizes = FALSE]

cat("Dream samples:", ncol(dge_all), "\n")
cat("Dream genes:", nrow(dge_all), "\n")
cat("Batch levels:", paste(sort(unique(dge_all$samples$dataset_subbatch)), collapse = ", "), "\n")
cat("Sex levels:", paste(levels(dge_all$samples$sex_for_model), collapse = ", "), "\n")
cat("Group distribution:\n")
print(table(dge_all$samples$masl_status, dge_all$samples$dataset))

form <- ~ masl_status + sex_for_model + (1 | dataset_subbatch)
cat("Formula:", deparse(form), "\n")

n_cores <- min(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4")), 32)
cat("Using", n_cores, "CPU cores\n")
BPPARAM <- MulticoreParam(n_cores, progressbar = TRUE)

vobjDream <- voomWithDreamWeights(dge_all, form, dge_all$samples, BPPARAM = BPPARAM)
cat("Running dream...\n")
fitDream <- dream(vobjDream, form, dge_all$samples, BPPARAM = BPPARAM)

dream_tt <- topTable(fitDream, coef = "masl_statusMASL", number = Inf, sort.by = "none")
dream_tt$gene <- rownames(dream_tt)

sig <- sum(dream_tt$adj.P.Val < 0.05)
sig_up <- sum(dream_tt$adj.P.Val < 0.05 & dream_tt$logFC > 0)
sig_down <- sum(dream_tt$adj.P.Val < 0.05 & dream_tt$logFC < 0)
sig_lfc <- sum(dream_tt$adj.P.Val < 0.05 & abs(dream_tt$logFC) > 0.5)
cat("\n===== MASL vs Healthy DREAM RESULTS =====\n")
cat("Total genes tested:", nrow(dream_tt), "\n")
cat(sprintf("DEGs (padj<0.05): %d (Up: %d, Down: %d)\n", sig, sig_up, sig_down))
cat(sprintf("DEGs (padj<0.05, |logFC|>0.5): %d\n", sig_lfc))

fwrite(as.data.table(dream_tt), file.path(ODIR, "masl_vs_healthy_dream.csv"))
cat("\nSaved:", file.path(ODIR, "masl_vs_healthy_dream.csv"), "\n")
