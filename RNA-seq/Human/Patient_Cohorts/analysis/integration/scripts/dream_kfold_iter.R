#!/usr/bin/env Rscript
# dream_kfold_iter.R
# ---------------------------------------------------------------------------
# Pillar A — Cohort-stratified k-fold CV.
# Iters 1..10  use K=10 (each fold = 10% test from each cohort).
# Iters 11..60 use K=50 (each fold = 2%  test from each cohort).
#
# Env: KFOLD_ITER (1..60)
# Out: results/integration/robustness/kfold/iter_{KFOLD_ITER}.csv
#      The header row records K and fold so the aggregator can stratify.
# ---------------------------------------------------------------------------

t0 <- proc.time()
suppressPackageStartupMessages({
  library(reformulas); library(lme4); library(data.table); library(edgeR); library(yaml)
})
ns_lme4 <- asNamespace("lme4")
for (fn in c("findbars", "nobars", "subbars", "rebuildFormula")) {
  if (exists(fn, envir = ns_lme4)) {
    try({ unlockBinding(fn, ns_lme4)
          assign(fn, get(fn, asNamespace("reformulas")), envir = ns_lme4)
          lockBinding(fn, ns_lme4) }, silent = TRUE)
  }
}
suppressPackageStartupMessages({ library(variancePartition); library(BiocParallel) })

ITER <- as.integer(Sys.getenv("KFOLD_ITER", "1"))
if (is.na(ITER) || ITER < 1 || ITER > 60) stop("KFOLD_ITER must be 1..60")

if (ITER <= 10) {
  K <- 10L; fold <- ITER
} else {
  K <- 50L; fold <- ITER - 10L
}
cat("=== K-fold iter", ITER, "(K =", K, ", fold =", fold, ") ===\n")

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts"
RDIR <- file.path(BASE, "analysis/integration/results/integration")
OUT_DIR <- file.path(RDIR, "robustness/kfold")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

dge <- readRDS(file.path(RDIR, "merged_dge.rds"))
ycfg <- yaml::read_yaml(file.path(
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
  "config/human_datasets.yaml"))$datasets
mega_cohorts <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))
dge <- dge[, dge$samples$dataset %in% mega_cohorts]

meta_new <- readRDS(file.path(RDIR, "meta_matched.rds"))
matched_sex <- meta_new$inferred_sex[match(colnames(dge), meta_new$sample_id)]
info <- data.frame(
  group_binary = factor(dge$samples$group_binary, levels = c("Control", "Disease")),
  dataset      = droplevels(factor(dge$samples$dataset)),
  inferred_sex = factor(matched_sex),
  stringsAsFactors = FALSE)
rownames(info) <- colnames(dge)
na_sex <- is.na(info$inferred_sex)
if (any(na_sex)) { dge <- dge[, !na_sex]; info <- info[!na_sex, , drop = FALSE] }

# --- Deterministic cohort-stratified fold assignment per K ---
# Same seed regardless of fold so all 10 (or 50) folds share one partition.
set.seed(42 + K * 1009L)
fold_assign <- integer(nrow(info))
for (ds in unique(info$dataset)) {
  ds_idx <- which(info$dataset == ds)
  fold_assign[ds_idx] <- sample(rep(seq_len(K), length.out = length(ds_idx)))
}
train_idx <- which(fold_assign != fold)
cat("Train n=", length(train_idx), "Test n=", nrow(info) - length(train_idx), "\n")

dge_t <- dge[, train_idx]; info_t <- info[train_idx, , drop = FALSE]
info_t$dataset <- droplevels(info_t$dataset)
ds_ok <- names(which(apply(table(info_t$dataset, info_t$group_binary), 1, function(x) all(x > 0))))
if (length(ds_ok) < 2) {
  cat("WARN: <2 datasets retain both groups in train — writing empty\n")
  fwrite(data.table(gene = character(0), logFC = numeric(0), t = numeric(0), padj = numeric(0),
                    K = K, fold = fold),
         file.path(OUT_DIR, sprintf("iter_%02d.csv", ITER)))
  quit(save = "no", status = 0)
}
keep <- info_t$dataset %in% ds_ok
if (!all(keep)) { dge_t <- dge_t[, keep]; info_t <- info_t[keep, , drop = FALSE]
                   info_t$dataset <- droplevels(info_t$dataset) }

ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
param <- if (ncpus > 1) MulticoreParam(ncpus) else SerialParam()
form <- ~ group_binary + inferred_sex + (1|dataset)
v <- suppressWarnings(voomWithDreamWeights(dge_t, form, info_t, BPPARAM = param))
fit <- suppressWarnings(dream(v, form, info_t, BPPARAM = param))
res <- topTable(fit, coef = "group_binaryDisease", number = Inf, sort.by = "none")
res$gene <- rownames(res); res_dt <- as.data.table(res); setnames(res_dt, "adj.P.Val", "padj")
res_dt[, K := K]; res_dt[, fold := fold]

fwrite(res_dt[, .(gene, logFC, t, padj, K, fold)],
       file.path(OUT_DIR, sprintf("iter_%02d.csv", ITER)))
cat(sprintf("Elapsed %.1f min  -- iter %d done\n", (proc.time()-t0)["elapsed"]/60, ITER))
