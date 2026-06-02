#!/usr/bin/env Rscript
# dream_bootstrap_iter.R
# ---------------------------------------------------------------------------
# Pillar A — Cohort-stratified subsampling-without-replacement (50%, B=1000).
# (Subsampling, not replacement-bootstrap: De Bin 2016 Biometrics shows ties
# from with-replacement bootstrap distort rank-based selection.)
#
# Env: BOOT_ITER (1..1000)
# Out: results/integration/robustness/bootstrap/iter_{BOOT_ITER}.csv
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

ITER <- as.integer(Sys.getenv("BOOT_ITER", "1"))
if (is.na(ITER) || ITER < 1) stop("BOOT_ITER must be >=1")
cat("=== Bootstrap iter", ITER, "===\n")

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts"
RDIR <- file.path(BASE, "analysis/integration/results/integration")
OUT_DIR <- file.path(RDIR, "robustness/bootstrap")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

out_file <- file.path(OUT_DIR, sprintf("iter_%04d.csv", ITER))
if (file.exists(out_file)) {
  cat("[Bootstrap iter", ITER, "] Already done — skipping\n")
  quit(save = "no", status = 0)
}

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

# --- Cohort-stratified 50% subsample (no replacement) ---
set.seed(42 + ITER * 7919)
sub_idx <- integer(0)
for (ds in unique(info$dataset)) {
  ds_idx <- which(info$dataset == ds)
  sub_idx <- c(sub_idx, sample(ds_idx, floor(length(ds_idx) / 2)))
}
cat("Subsample n =", length(sub_idx), "/", nrow(info), "\n")

dge_s <- dge[, sub_idx]; info_s <- info[sub_idx, , drop = FALSE]
info_s$dataset <- droplevels(info_s$dataset)
ds_ok <- names(which(apply(table(info_s$dataset, info_s$group_binary), 1, function(x) all(x > 0))))
if (length(ds_ok) < 2) {
  cat("WARN: <2 datasets retain both groups; writing empty result\n")
  fwrite(data.table(gene = character(0), logFC = numeric(0), t = numeric(0), padj = numeric(0)),
         file.path(OUT_DIR, sprintf("iter_%04d.csv", ITER)))
  quit(save = "no", status = 0)
}
keep <- info_s$dataset %in% ds_ok
if (!all(keep)) { dge_s <- dge_s[, keep]; info_s <- info_s[keep, , drop = FALSE]
                   info_s$dataset <- droplevels(info_s$dataset) }

ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
param <- if (ncpus > 1) MulticoreParam(ncpus) else SerialParam()
form <- ~ group_binary + inferred_sex + (1|dataset)
cat("Fitting on n=", ncol(dge_s), "across", nlevels(info_s$dataset), "datasets\n")
v <- suppressWarnings(voomWithDreamWeights(dge_s, form, info_s, BPPARAM = param))
fit <- suppressWarnings(dream(v, form, info_s, BPPARAM = param))
res <- topTable(fit, coef = "group_binaryDisease", number = Inf, sort.by = "none")
res$gene <- rownames(res); res_dt <- as.data.table(res); setnames(res_dt, "adj.P.Val", "padj")
fwrite(res_dt[, .(gene, logFC, t, padj)],
       file.path(OUT_DIR, sprintf("iter_%04d.csv", ITER)))
cat(sprintf("Elapsed %.1f min  -- iter %d done\n", (proc.time()-t0)["elapsed"]/60, ITER))
