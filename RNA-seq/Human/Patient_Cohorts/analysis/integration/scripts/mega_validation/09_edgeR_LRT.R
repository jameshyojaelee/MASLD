#!/usr/bin/env Rscript
# 09_edgeR_LRT.R — Arm 8: edgeR likelihood-ratio test.
# Same design as edgeR-QL (Arm 1) but using glmFit + glmLRT instead of
# glmQLFit + glmQLFTest. At n=847 the two should give nearly identical
# results; included as a robustness check on the QL choice.

suppressPackageStartupMessages({
  library(edgeR); library(data.table); library(yaml)
})

PROJECT_ROOT <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT <- file.path(PROJECT_ROOT, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR_IN  <- file.path(INT, "results/integration")
RDIR_OUT <- file.path(INT, "results/mega_validation/edgeR_LRT")
dir.create(RDIR_OUT, recursive = TRUE, showWarnings = FALSE)

norm_lib_sizes <- function(x, ...) {
  if (exists("normLibSizes", where = asNamespace("edgeR"), mode = "function"))
    edgeR::normLibSizes(x, ...) else edgeR::calcNormFactors(x, ...)
}

ycfg <- read_yaml(file.path(PROJECT_ROOT, "config/human_datasets.yaml"))$datasets
mega_cohorts <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))

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
sex_na <- is.na(dge_mega$samples$inferred_sex)
if (any(sex_na)) {
  dge_mega <- dge_mega[, !sex_na]
  dge_mega$samples$inferred_sex <- droplevels(dge_mega$samples$inferred_sex)
  dge_mega$samples$dataset      <- droplevels(dge_mega$samples$dataset)
}
cat("Samples:", ncol(dge_mega), "\n")

design <- model.matrix(~ dataset + inferred_sex + group_binary,
                       data = dge_mega$samples)
keep <- filterByExpr(dge_mega, design = design,
                     group = dge_mega$samples$group_binary)
cat("Genes after filterByExpr:", sum(keep), "\n")
dge_mega <- dge_mega[keep, , keep.lib.sizes = FALSE]
dge_mega <- norm_lib_sizes(dge_mega, method = "TMM")

dge_mega <- estimateDisp(dge_mega, design, robust = TRUE)
fit <- glmFit(dge_mega, design)
lrt <- glmLRT(fit, coef = "group_binaryDisease")
res <- topTags(lrt, n = Inf, sort.by = "none")$table
res$gene <- rownames(res)
res_dt <- as.data.table(res)
setnames(res_dt, c("PValue", "FDR"), c("P.Value", "padj"))

cat("\n===== edgeR-LRT RESULTS =====\n")
cat("Genes tested:", nrow(res_dt), "\n")
cat("DEGs padj<0.05:", res_dt[padj < 0.05, .N], "\n")
cat("DEGs Tier 1 (|LFC|>0.5):",
    res_dt[padj < 0.05 & abs(logFC) > 0.5, .N], "\n")

fwrite(res_dt, file.path(RDIR_OUT, "edger_lrt_results.csv"))
fwrite(data.table(metric = c("n_samples","n_genes","n_DEG_padj005",
                              "n_DEG_padj005_lfc05"),
                  value  = c(ncol(dge_mega), nrow(res_dt),
                             res_dt[padj < 0.05, .N],
                             res_dt[padj < 0.05 & abs(logFC) > 0.5, .N])),
       file.path(RDIR_OUT, "qc.csv"))
