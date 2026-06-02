#!/usr/bin/env Rscript
# 06_DESeq2_cohortFE.R
# ---------------------------------------------------------------------------
# Mega-validation Arm 5: DESeq2 (Love et al. 2014) Wald test with cohort
# as fixed effect. DESeq2 is the most-cited NB-family DE tool; included as
# a reviewer-armor sensitivity arm. Uses the same input matrix + design as
# Arms 1 (edgeR-QL) and dream.
#
# Inputs:  merged_dge.rds + meta_matched.rds + sample_qc_report.csv +
#          config/human_datasets.yaml (include_in_mega filter)
# Outputs: results/mega_validation/DESeq2/{deseq2_results.csv, qc.csv}
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(DESeq2); library(data.table); library(yaml)
})

PROJECT_ROOT <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT <- file.path(PROJECT_ROOT, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR_IN  <- file.path(INT, "results/integration")
RDIR_OUT <- file.path(INT, "results/mega_validation/DESeq2")
dir.create(RDIR_OUT, recursive = TRUE, showWarnings = FALSE)

ycfg <- read_yaml(file.path(PROJECT_ROOT, "config/human_datasets.yaml"))$datasets
mega_cohorts <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))
cat("Mega cohorts (k =", length(mega_cohorts), "):",
    paste(mega_cohorts, collapse = ", "), "\n")

# --- Load data (same as edgeR-QL arm) ---
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

# Drop NA-sex samples (same as dream + edgeR-QL)
sex_na <- is.na(dge_mega$samples$inferred_sex)
if (any(sex_na)) {
  cat("Dropping", sum(sex_na), "samples with NA inferred_sex\n")
  dge_mega <- dge_mega[, !sex_na]
  dge_mega$samples$inferred_sex <- droplevels(dge_mega$samples$inferred_sex)
  dge_mega$samples$dataset      <- droplevels(dge_mega$samples$dataset)
}
cat("Samples after filter:", ncol(dge_mega), "\n")

# --- Light prefilter (parallels filterByExpr scale; speeds DESeq2) ---
# DESeq2 prefers integer counts + does its own filtering on row sums.
counts_mat <- as.matrix(dge_mega$counts)
storage.mode(counts_mat) <- "integer"
keep <- rowSums(counts_mat >= 10) >= 10   # at least 10 samples with >=10 reads
cat("Genes after prefilter:", sum(keep), "/", length(keep), "\n")
counts_mat <- counts_mat[keep, ]

# --- DESeqDataSet ---
coldata <- DataFrame(
  dataset      = dge_mega$samples$dataset,
  inferred_sex = dge_mega$samples$inferred_sex,
  group_binary = dge_mega$samples$group_binary
)
rownames(coldata) <- colnames(counts_mat)

dds <- DESeqDataSetFromMatrix(
  countData = counts_mat,
  colData   = coldata,
  design    = ~ dataset + inferred_sex + group_binary
)

# --- Fit (Wald) ---
ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
cat("Using", ncpus, "CPUs (BiocParallel)\n")
if (ncpus > 1) {
  suppressPackageStartupMessages(library(BiocParallel))
  bp <- MulticoreParam(ncpus)
  dds <- DESeq(dds, test = "Wald", parallel = TRUE, BPPARAM = bp)
} else {
  dds <- DESeq(dds, test = "Wald")
}

# --- Extract disease vs control ---
res <- results(dds, name = "group_binary_Disease_vs_Control")
cat("DESeq2 result columns:\n"); print(colnames(res))

res_dt <- as.data.table(as.data.frame(res), keep.rownames = "gene")
setnames(res_dt,
  c("log2FoldChange", "lfcSE", "stat", "pvalue", "padj"),
  c("logFC",          "SE",    "t",    "P.Value", "padj_BH"))
# DESeq2 already applies independent filtering; padj_BH is its FDR.
res_dt[, padj := padj_BH]  # canonical name across mega_validation arms

cat("\n===== DESeq2 RESULTS =====\n")
cat("Genes tested:", nrow(res_dt), "\n")
cat("With non-NA padj (post indep. filter):", sum(!is.na(res_dt$padj)), "\n")
cat("DEGs padj<0.05:", sum(res_dt$padj < 0.05, na.rm = TRUE), "\n")
cat("DEGs padj<0.05 & |logFC|>0.5:",
    sum(res_dt$padj < 0.05 & abs(res_dt$logFC) > 0.5, na.rm = TRUE), "\n")

fwrite(res_dt, file.path(RDIR_OUT, "deseq2_results.csv"))

qc_dt <- data.table(
  metric = c("n_samples", "n_genes_tested", "n_genes_padj_non_NA",
             "n_DEG_padj005", "n_DEG_padj005_lfc05"),
  value  = c(ncol(counts_mat), nrow(res_dt),
             sum(!is.na(res_dt$padj)),
             sum(res_dt$padj < 0.05, na.rm = TRUE),
             sum(res_dt$padj < 0.05 & abs(res_dt$logFC) > 0.5, na.rm = TRUE))
)
fwrite(qc_dt, file.path(RDIR_OUT, "qc.csv"))

cat("\nSaved:", file.path(RDIR_OUT, "deseq2_results.csv"), "\n")
