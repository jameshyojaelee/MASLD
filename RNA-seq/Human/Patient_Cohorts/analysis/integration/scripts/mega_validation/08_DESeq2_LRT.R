#!/usr/bin/env Rscript
# 08_DESeq2_LRT.R — Arm 7: DESeq2 likelihood-ratio test.
# Reduced model: ~ dataset + inferred_sex
# Full model:    ~ dataset + inferred_sex + group_binary
# Robustness check against the DESeq2-Wald result (arm 5).

suppressPackageStartupMessages({
  library(DESeq2); library(data.table); library(yaml); library(BiocParallel)
})

PROJECT_ROOT <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT <- file.path(PROJECT_ROOT, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR_IN  <- file.path(INT, "results/integration")
RDIR_OUT <- file.path(INT, "results/mega_validation/DESeq2_LRT")
dir.create(RDIR_OUT, recursive = TRUE, showWarnings = FALSE)

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

counts_mat <- as.matrix(dge_mega$counts)
storage.mode(counts_mat) <- "integer"
keep <- rowSums(counts_mat >= 10) >= 10
counts_mat <- counts_mat[keep, ]
cat("Samples:", ncol(counts_mat), "  Genes:", nrow(counts_mat), "\n")

coldata <- DataFrame(
  dataset      = dge_mega$samples$dataset,
  inferred_sex = dge_mega$samples$inferred_sex,
  group_binary = dge_mega$samples$group_binary
)
rownames(coldata) <- colnames(counts_mat)

dds <- DESeqDataSetFromMatrix(counts_mat, coldata,
  design = ~ dataset + inferred_sex + group_binary)

ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
cat("Using", ncpus, "CPUs\n")
if (ncpus > 1) {
  bp <- MulticoreParam(ncpus)
  dds <- DESeq(dds, test = "LRT", reduced = ~ dataset + inferred_sex,
               parallel = TRUE, BPPARAM = bp)
} else {
  dds <- DESeq(dds, test = "LRT", reduced = ~ dataset + inferred_sex)
}

res <- results(dds)
res_dt <- as.data.table(as.data.frame(res), keep.rownames = "gene")
setnames(res_dt,
  c("log2FoldChange", "lfcSE", "stat",     "pvalue", "padj"),
  c("logFC",          "SE",    "LRT_stat", "P.Value", "padj_BH"))
res_dt[, padj := padj_BH]
# DESeq2-LRT's `stat` is a chi-square statistic (always positive). Convert
# to a signed-sqrt(stat) "t-stat proxy" so it can be compared to signed
# t-statistics from Wald-based methods.
res_dt[, t := sign(logFC) * sqrt(pmax(LRT_stat, 0))]

cat("\n===== DESeq2-LRT RESULTS =====\n")
cat("Genes tested:", nrow(res_dt), "\n")
cat("DEGs padj<0.05:", sum(res_dt$padj < 0.05, na.rm = TRUE), "\n")
cat("DEGs Tier 1 (|LFC|>0.5):",
    sum(res_dt$padj < 0.05 & abs(res_dt$logFC) > 0.5, na.rm = TRUE), "\n")

fwrite(res_dt, file.path(RDIR_OUT, "deseq2_lrt_results.csv"))
fwrite(data.table(metric = c("n_samples","n_genes","n_DEG_padj005",
                              "n_DEG_padj005_lfc05"),
                  value  = c(ncol(counts_mat), nrow(res_dt),
                             sum(res_dt$padj < 0.05, na.rm = TRUE),
                             sum(res_dt$padj < 0.05 & abs(res_dt$logFC) > 0.5,
                                 na.rm = TRUE))),
       file.path(RDIR_OUT, "qc.csv"))
