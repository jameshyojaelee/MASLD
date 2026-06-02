#!/usr/bin/env Rscript
# 02_voomLmFit_block.R
# ---------------------------------------------------------------------------
# Mega-validation Arm 2: classic two-round voom + duplicateCorrelation
# (Smyth & Altman 2003 idiom; what limma::voomLmFit wraps in newer versions).
# This is the genome-wide-single-rho counterpart to dream's per-gene tau^2;
# expected to slightly under-perform on genes with above-average cohort RE
# variance (Hoffman & Roussos 2020 TUBB2B example).
#
# Note: voomLmFit() is the convenience wrapper added in limma >= 3.56,
# but it is not exported in our current build. Implementing the underlying
# two-round idiom directly gives the same statistical result.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(limma); library(edgeR); library(data.table); library(yaml)
})

PROJECT_ROOT <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT <- file.path(PROJECT_ROOT, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR_IN  <- file.path(INT, "results/integration")
RDIR_OUT <- file.path(INT, "results/mega_validation/voomLmFit")
dir.create(RDIR_OUT, recursive = TRUE, showWarnings = FALSE)

norm_lib_sizes <- function(x, ...) {
  if (exists("normLibSizes", where = asNamespace("edgeR"), mode = "function")) {
    edgeR::normLibSizes(x, ...)
  } else {
    edgeR::calcNormFactors(x, ...)
  }
}

ycfg <- read_yaml(file.path(PROJECT_ROOT, "config/human_datasets.yaml"))$datasets
mega_cohorts <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))
cat("Mega cohorts:", paste(mega_cohorts, collapse = ", "), "\n")

dge  <- readRDS(file.path(RDIR_IN, "merged_dge.rds"))
meta <- readRDS(file.path(RDIR_IN, "meta_matched.rds"))
qc   <- fread(file.path(INT, "qc/sample_qc_report.csv"))
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

# Design: group_binary first; dataset enters via block in duplicateCorrelation
design <- model.matrix(~ group_binary + inferred_sex, data = dge_mega$samples)
# group_binary is the contrast of interest; filterByExpr uses it for sizing.
keep <- filterByExpr(dge_mega, design = design,
                     group = dge_mega$samples$group_binary)
cat("Genes after filter:", sum(keep), "\n")
dge_mega <- dge_mega[keep, , keep.lib.sizes = FALSE]
dge_mega <- norm_lib_sizes(dge_mega, method = "TMM")

block <- dge_mega$samples$dataset

# --- Two-round voom + duplicateCorrelation ---
cat("Round 1: voom + duplicateCorrelation...\n")
v1   <- voom(dge_mega, design, plot = FALSE)
cor1 <- duplicateCorrelation(v1, design, block = block)
cat("  Round 1 consensus rho:", round(cor1$consensus.correlation, 4), "\n")

cat("Round 2: voom (block-aware) + duplicateCorrelation...\n")
v2   <- voom(dge_mega, design, plot = FALSE,
             block = block, correlation = cor1$consensus.correlation)
cor2 <- duplicateCorrelation(v2, design, block = block)
cat("  Round 2 consensus rho:", round(cor2$consensus.correlation, 4), "\n")

# --- lmFit + eBayes ---
fit <- lmFit(v2, design, block = block,
             correlation = cor2$consensus.correlation)
fit <- eBayes(fit, robust = TRUE)

res <- topTable(fit, coef = "group_binaryDisease",
                number = Inf, sort.by = "none")
res$gene <- rownames(res)
res_dt <- as.data.table(res)
setnames(res_dt, "adj.P.Val", "padj")

intra_rho <- cor2$consensus.correlation
cat("\n===== voomLmFit + block=dataset (manual idiom) RESULTS =====\n")
cat("Genes tested:", nrow(res_dt), "\n")
cat("DEGs padj<0.05:", res_dt[padj < 0.05, .N], "\n")
cat("DEGs padj<0.1: ", res_dt[padj < 0.10, .N], "\n")
cat("Final intra-cohort consensus correlation:", round(intra_rho, 4), "\n")

fwrite(res_dt, file.path(RDIR_OUT, "vlm_results.csv"))
fwrite(data.table(intra_rho_round1 = cor1$consensus.correlation,
                  intra_rho_round2 = intra_rho,
                  n_samples = ncol(dge_mega),
                  n_genes   = nrow(res_dt)),
       file.path(RDIR_OUT, "intra_rho.csv"))
