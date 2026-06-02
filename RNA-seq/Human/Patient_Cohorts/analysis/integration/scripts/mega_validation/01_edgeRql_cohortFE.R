#!/usr/bin/env Rscript
# 01_edgeRql_cohortFE.R
# ---------------------------------------------------------------------------
# Mega-validation Arm 1: edgeR v4 quasi-likelihood F-test with cohort as
# fixed effect. Chen 2024 v4 uses bias-corrected QL deviances explicitly
# targeting large-n FDR control. With K=5 cohorts the RE-vs-FE distinction
# is statistically negligible (Hodges 2014).
#
# Inputs:  merged_dge.rds + meta_matched.rds + sample_qc_report.csv +
#          config/human_datasets.yaml (include_in_mega filter)
# Outputs: results/mega_validation/edgeRql/{ql_results.csv, qc.csv}
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(edgeR); library(data.table); library(yaml)
})

PROJECT_ROOT <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT <- file.path(PROJECT_ROOT, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR_IN  <- file.path(INT, "results/integration")
RDIR_OUT <- file.path(INT, "results/mega_validation/edgeRql")
dir.create(RDIR_OUT, recursive = TRUE, showWarnings = FALSE)

# normLibSizes is the v4 alias for calcNormFactors. Fall back for older edgeR.
norm_lib_sizes <- function(x, ...) {
  if (exists("normLibSizes", where = "package:edgeR", mode = "function")) {
    edgeR::normLibSizes(x, ...)
  } else {
    edgeR::calcNormFactors(x, ...)
  }
}

# --- Mega-cohort allowlist from yaml ---
ycfg <- read_yaml(file.path(PROJECT_ROOT, "config/human_datasets.yaml"))$datasets
mega_cohorts <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))
cat("Mega cohorts (k =", length(mega_cohorts), "):",
    paste(mega_cohorts, collapse = ", "), "\n")

# --- Load data ---
dge <- readRDS(file.path(RDIR_IN, "merged_dge.rds"))
meta <- readRDS(file.path(RDIR_IN, "meta_matched.rds"))
qc <- fread(file.path(INT, "qc/sample_qc_report.csv"))
pass_ids <- qc[pass_technical == TRUE, sample_id]
keep_sample <- dge$samples$dataset %in% mega_cohorts &
               colnames(dge) %in% pass_ids
dge_mega <- dge[, keep_sample]

matched_sex <- meta$inferred_sex[match(colnames(dge_mega), meta$sample_id)]
dge_mega$samples$inferred_sex <- factor(matched_sex)
dge_mega$samples$group_binary <- factor(dge_mega$samples$group_binary,
                                        levels = c("Control", "Disease"))
dge_mega$samples$dataset <- factor(dge_mega$samples$dataset)

# Drop NA-sex samples (consistent with dream's behaviour)
sex_na <- is.na(dge_mega$samples$inferred_sex)
if (any(sex_na)) {
  cat("Dropping", sum(sex_na), "samples with NA inferred_sex\n")
  dge_mega <- dge_mega[, !sex_na]
  dge_mega$samples$inferred_sex <- droplevels(dge_mega$samples$inferred_sex)
}

cat("Samples after filter:", ncol(dge_mega), "\n")
cat("Group × dataset table:\n"); print(table(dge_mega$samples$group_binary,
                                              dge_mega$samples$dataset))

# --- Design ---
design <- model.matrix(~ dataset + inferred_sex + group_binary,
                       data = dge_mega$samples)
cat("Design columns:\n"); print(colnames(design))

# --- Filter + norm ---
# filterByExpr uses `group` for "min.count in smallest group" sizing.
# Use group_binary (contrast of interest) for fair head-to-head with dream,
# not dataset (which would let smallest cohort drive the threshold and
# retain a more lenient gene set than the canonical analysis).
keep <- filterByExpr(dge_mega, design = design,
                     group = dge_mega$samples$group_binary)
cat("Genes after filterByExpr:", sum(keep), "/", length(keep), "\n")
dge_mega <- dge_mega[keep, , keep.lib.sizes = FALSE]
dge_mega <- norm_lib_sizes(dge_mega, method = "TMM")

# --- Dispersion + QL fit (Chen 2024 v4 bias-corrected) ---
ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
cat("Using", ncpus, "CPUs\n")
dge_mega <- estimateDisp(dge_mega, design, robust = TRUE)

fit <- glmQLFit(dge_mega, design, robust = TRUE)
qlf <- glmQLFTest(fit, coef = "group_binaryDisease")

res <- topTags(qlf, n = Inf, sort.by = "none")$table
res$gene <- rownames(res)
res_dt <- as.data.table(res)
setnames(res_dt, c("PValue", "FDR"), c("P.Value", "padj"))

# --- Summary ---
cat("\n===== edgeR-QL v4 RESULTS =====\n")
cat("Genes tested:", nrow(res_dt), "\n")
cat("DEGs padj<0.05:", res_dt[padj < 0.05, .N], "\n")
cat("DEGs padj<0.05 & |logFC|>0.5:",
    res_dt[padj < 0.05 & abs(logFC) > 0.5, .N], "\n")

fwrite(res_dt, file.path(RDIR_OUT, "ql_results.csv"))

# --- QC summary ---
qc_dt <- data.table(
  metric  = c("n_samples", "n_genes_tested", "n_DEG_padj005",
              "n_DEG_padj005_lfc05", "common_disp", "trended_disp_median"),
  value   = c(ncol(dge_mega), nrow(res_dt),
              res_dt[padj < 0.05, .N],
              res_dt[padj < 0.05 & abs(logFC) > 0.5, .N],
              dge_mega$common.dispersion,
              median(dge_mega$trended.dispersion, na.rm = TRUE))
)
fwrite(qc_dt, file.path(RDIR_OUT, "qc.csv"))

cat("\nSaved:", file.path(RDIR_OUT, "ql_results.csv"), "\n")
cat("Saved:", file.path(RDIR_OUT, "qc.csv"), "\n")
