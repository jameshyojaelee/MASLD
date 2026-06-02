#!/usr/bin/env Rscript
# 10_voom_sva.R — Arm 9: voom + sva (surrogate variable analysis).
# Tests whether unmeasured latent batch variance remains after cohort FE
# is accounted for. If SVs explain little additional variance, that's
# evidence the cohort RE / FE specification was sufficient.

suppressPackageStartupMessages({
  library(limma); library(edgeR); library(sva); library(data.table); library(yaml)
})

PROJECT_ROOT <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT <- file.path(PROJECT_ROOT, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR_IN  <- file.path(INT, "results/integration")
RDIR_OUT <- file.path(INT, "results/mega_validation/voom_sva")
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

design_full <- model.matrix(~ dataset + inferred_sex + group_binary,
                            data = dge_mega$samples)
design_null <- model.matrix(~ dataset + inferred_sex,
                            data = dge_mega$samples)

keep <- filterByExpr(dge_mega, design = design_full,
                     group = dge_mega$samples$group_binary)
dge_mega <- dge_mega[keep, , keep.lib.sizes = FALSE]
dge_mega <- norm_lib_sizes(dge_mega, method = "TMM")
cat("Genes after filter:", nrow(dge_mega), "\n")

# Voom first for the log-CPM matrix that sva expects
v <- voom(dge_mega, design_full, plot = FALSE)
y <- v$E

cat("Estimating SVs (num.sv on protected design)...\n")
n_sv <- num.sv(y, design_full, method = "be")
cat("num.sv (Buja-Eyuboglu) =", n_sv, "\n")
n_sv <- max(1, min(n_sv, 20))   # cap at 20 for stability + speed
cat("Using", n_sv, "SVs\n")

sva_obj <- sva(y, mod = design_full, mod0 = design_null, n.sv = n_sv)
SVs <- sva_obj$sv
colnames(SVs) <- paste0("SV", seq_len(ncol(SVs)))

# Re-fit with SVs added to design
design_sva <- cbind(design_full, SVs)
v2 <- voom(dge_mega, design_sva, plot = FALSE)
fit <- lmFit(v2, design_sva)
fit <- eBayes(fit, robust = TRUE)
res <- topTable(fit, coef = "group_binaryDisease",
                number = Inf, sort.by = "none")
res$gene <- rownames(res)
res_dt <- as.data.table(res)
setnames(res_dt, "adj.P.Val", "padj")

cat("\n===== voom + SVA RESULTS =====\n")
cat("Genes tested:", nrow(res_dt), "\n")
cat("DEGs padj<0.05:", res_dt[padj < 0.05, .N], "\n")
cat("DEGs Tier 1 (|LFC|>0.5):",
    res_dt[padj < 0.05 & abs(logFC) > 0.5, .N], "\n")

fwrite(res_dt, file.path(RDIR_OUT, "voom_sva_results.csv"))
# SVs is a matrix with sample-rows; attach sample IDs before writing
sv_dt <- data.table(sample_id = colnames(dge_mega))
sv_dt <- cbind(sv_dt, as.data.table(SVs))
fwrite(sv_dt, file.path(RDIR_OUT, "sva_factors.csv"))
fwrite(data.table(metric = c("n_samples","n_genes","n_SVs","n_DEG_padj005",
                              "n_DEG_padj005_lfc05"),
                  value  = c(ncol(dge_mega), nrow(res_dt), n_sv,
                             res_dt[padj < 0.05, .N],
                             res_dt[padj < 0.05 & abs(logFC) > 0.5, .N])),
       file.path(RDIR_OUT, "qc.csv"))
