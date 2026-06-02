#!/usr/bin/env Rscript
# 07_limma_trend.R — Arm 6: limma-trend (no voom precision weights).
# log2-CPM → lmFit → eBayes(trend=TRUE). Tests whether voom's mean-variance
# weighting is actually doing work versus a simpler trend-only moderation.

suppressPackageStartupMessages({
  library(limma); library(edgeR); library(data.table); library(yaml)
})

PROJECT_ROOT <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT <- file.path(PROJECT_ROOT, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR_IN  <- file.path(INT, "results/integration")
RDIR_OUT <- file.path(INT, "results/mega_validation/limma_trend")
dir.create(RDIR_OUT, recursive = TRUE, showWarnings = FALSE)

norm_lib_sizes <- function(x, ...) {
  if (exists("normLibSizes", where = asNamespace("edgeR"), mode = "function"))
    edgeR::normLibSizes(x, ...) else edgeR::calcNormFactors(x, ...)
}

ycfg <- read_yaml(file.path(PROJECT_ROOT, "config/human_datasets.yaml"))$datasets
mega_cohorts <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))
cat("Mega cohorts:", paste(mega_cohorts, collapse = ", "), "\n")

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

# log2-CPM (no precision weights — that's the point of "trend")
logCPM <- cpm(dge_mega, log = TRUE, prior.count = 3)

fit <- lmFit(logCPM, design)
fit <- eBayes(fit, trend = TRUE, robust = TRUE)
res <- topTable(fit, coef = "group_binaryDisease",
                number = Inf, sort.by = "none")
res$gene <- rownames(res)
res_dt <- as.data.table(res)
setnames(res_dt, "adj.P.Val", "padj")

cat("\n===== limma-trend RESULTS =====\n")
cat("Genes tested:", nrow(res_dt), "\n")
cat("DEGs padj<0.05:", res_dt[padj < 0.05, .N], "\n")
cat("DEGs Tier 1 (|LFC|>0.5):",
    res_dt[padj < 0.05 & abs(logFC) > 0.5, .N], "\n")

fwrite(res_dt, file.path(RDIR_OUT, "trend_results.csv"))
fwrite(data.table(metric = c("n_samples","n_genes","n_DEG_padj005",
                              "n_DEG_padj005_lfc05"),
                  value  = c(ncol(dge_mega), nrow(res_dt),
                             res_dt[padj < 0.05, .N],
                             res_dt[padj < 0.05 & abs(logFC) > 0.5, .N])),
       file.path(RDIR_OUT, "qc.csv"))
