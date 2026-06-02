#!/usr/bin/env Rscript
# 05e_mash_vs_masl_dream.R
# ---------------------------------------------------------------------------
# MASH-vs-MASL pooled binary dream() mega-analysis.
#
# PRJNA512027 permanently removed from pipeline 2026-05-15 (L0/S0 library batch
# perfectly confounded with disease severity, same exclusion rationale as
# 05_dream_mega_analysis.R).
#
# MASH_DEF env var selects the contrast definition:
#   - "borderline_grouped" (default): NAFL vs (NASH + Borderline NAS 3-4)
#       writes results/disease_signatures/mash_vs_masl_dream.csv
#       (replaces nafl_vs_nash_dream.csv from 13_nafl_vs_nash_de.R as the
#        canonical primary MASH-vs-MASL output)
#   - "strict": NAFL vs NASH only (Borderline excluded, NAS >= 5 = NASH)
#       writes results/disease_signatures/mash_vs_masl_dream_strict.csv
#
# Output goes to results/disease_signatures/. Same dream formula as 13:
#   ~ nafl_nash + sex_for_model + (1 | dataset)
# ---------------------------------------------------------------------------

set.seed(42)
Sys.setenv(R_PARALLEL_SEED = "42")

suppressPackageStartupMessages({
  library(reformulas)
  library(lme4)
  library(data.table)
  library(edgeR)
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
  library(limma)
  library(variancePartition)
  library(BiocParallel)
})

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts"
INT  <- file.path(BASE, "analysis/integration")
RDIR <- file.path(INT, "results/disease_signatures")
dir.create(RDIR, recursive = TRUE, showWarnings = FALSE)

mash_def <- Sys.getenv("MASH_DEF", "borderline_grouped")
if (!mash_def %in% c("borderline_grouped", "strict")) {
  stop("MASH_DEF must be 'borderline_grouped' or 'strict'; got: ", mash_def)
}
cat("MASH_DEF:", mash_def, "\n")

counts <- readRDS(file.path(INT, "results/integration/merged_counts_raw.rds"))
meta   <- readRDS(file.path(INT, "results/integration/meta_matched.rds"))
qc     <- fread(file.path(INT, "qc/sample_qc_report.csv"))
meta   <- meta[sample_id %in% qc[pass_technical == TRUE, sample_id]]

# Defensive filter (PRJNA512027 was removed from yaml 2026-05-15; this guards
# against stale meta_matched.rds inputs).
meta <- meta[dataset != "PRJNA512027"]

if (mash_def == "borderline_grouped") {
  meta_nn <- meta[diagnosis_harmonized %in% c("NAFL", "NASH", "Borderline")]
  meta_nn[, nafl_nash := factor(fifelse(diagnosis_harmonized == "NAFL", "NAFL", "NASH"),
                                levels = c("NAFL", "NASH"))]
} else {
  meta_nn <- meta[diagnosis_harmonized %in% c("NAFL", "NASH")]
  meta_nn[, nafl_nash := factor(diagnosis_harmonized, levels = c("NAFL", "NASH"))]
}

cat("Pre-cohort-filter sample counts:\n")
print(meta_nn[, .N, by = .(dataset, nafl_nash)][order(dataset, nafl_nash)])

ds_counts <- dcast(meta_nn[, .N, by = .(dataset, nafl_nash)],
                   dataset ~ nafl_nash, value.var = "N", fill = 0L)
keep_ds <- ds_counts[NAFL >= 3 & NASH >= 3, dataset]
cat("\nDatasets retained (>=3 in each arm):", paste(keep_ds, collapse = ", "), "\n")
dropped <- setdiff(ds_counts$dataset, keep_ds)
if (length(dropped) > 0) cat("Datasets dropped:", paste(dropped, collapse = ", "), "\n")
meta_nn <- meta_nn[dataset %in% keep_ds]
cat("Final samples:", nrow(meta_nn), "\n")

cat("\n===== DREAM MEGA-ANALYSIS (MASH vs MASL,", mash_def, ") =====\n")

all_nn <- meta_nn[order(sample_id)]
idx_all <- colnames(counts) %in% all_nn$sample_id
dge_all <- DGEList(counts = counts[, idx_all])
dge_all$samples <- cbind(dge_all$samples,
  all_nn[match(colnames(dge_all), all_nn$sample_id),
    .(nafl_nash, dataset, sex, inferred_sex, age)])

dge_all$samples$sex_for_model <- dge_all$samples$sex
na_sex <- is.na(dge_all$samples$sex_for_model) | dge_all$samples$sex_for_model == ""
dge_all$samples$sex_for_model[na_sex] <- as.character(dge_all$samples$inferred_sex[na_sex])
dge_all$samples$sex_for_model <- factor(dge_all$samples$sex_for_model)

dge_all <- calcNormFactors(dge_all, method = "TMM")
keep_all <- filterByExpr(dge_all, group = dge_all$samples$nafl_nash)
dge_all <- dge_all[keep_all, , keep.lib.sizes = FALSE]

cat("Dream samples:", ncol(dge_all), "\n")
cat("Dream genes:", nrow(dge_all), "\n")
cat("Datasets:", paste(sort(unique(dge_all$samples$dataset)), collapse = ", "), "\n")
cat("Sex levels:", paste(levels(dge_all$samples$sex_for_model), collapse = ", "), "\n")

form <- ~ nafl_nash + sex_for_model + (1 | dataset)
cat("Formula:", deparse(form), "\n")

n_cores <- min(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4")), 32)
cat("Using", n_cores, "CPU cores\n")
BPPARAM <- MulticoreParam(n_cores, progressbar = TRUE)

vobjDream <- voomWithDreamWeights(dge_all, form, dge_all$samples, BPPARAM = BPPARAM)
cat("Running dream...\n")
fitDream <- dream(vobjDream, form, dge_all$samples, BPPARAM = BPPARAM)

dream_tt <- topTable(fitDream, coef = "nafl_nashNASH", number = Inf, sort.by = "none")
dream_tt$gene <- rownames(dream_tt)
dream_tt <- as.data.table(dream_tt)

sig <- sum(dream_tt$adj.P.Val < 0.05)
sig_up <- sum(dream_tt$adj.P.Val < 0.05 & dream_tt$logFC > 0)
sig_down <- sum(dream_tt$adj.P.Val < 0.05 & dream_tt$logFC < 0)
sig_lfc <- sum(dream_tt$adj.P.Val < 0.05 & abs(dream_tt$logFC) > 0.5)
cat(sprintf("Total genes tested: %d\n", nrow(dream_tt)))
cat(sprintf("DEGs (padj<0.05): %d (Up: %d, Down: %d)\n", sig, sig_up, sig_down))
cat(sprintf("DEGs (padj<0.05 & |logFC|>0.5): %d\n", sig_lfc))

out_fname <- if (mash_def == "borderline_grouped") {
  "mash_vs_masl_dream.csv"
} else {
  "mash_vs_masl_dream_strict.csv"
}
fwrite(dream_tt, file.path(RDIR, out_fname))
cat("\nSaved:", file.path(RDIR, out_fname), "\n")
