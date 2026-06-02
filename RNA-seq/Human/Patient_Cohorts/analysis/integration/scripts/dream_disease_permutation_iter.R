#!/usr/bin/env Rscript
# dream_disease_permutation_iter.R
# ---------------------------------------------------------------------------
# Pillar C — Empirical null calibration via within-cohort label shuffling.
# Permutes group_binary WITHIN each dataset (preserves per-cohort case/control
# counts; tests the cohort-stratified null per Phipson & Smyth 2010).
#
# Env: PERM_ITER (1..1000)
# Out: results/audit_sensitivity/permutation_null/perm_iter_{N}.csv
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

PERM <- as.integer(Sys.getenv("PERM_ITER", "1"))
if (is.na(PERM) || PERM < 1) stop("PERM_ITER must be >=1")
cat("=== Disease-label permutation iter", PERM, "===\n")

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
INT  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR <- file.path(INT, "results/integration")
OUT_DIR <- file.path(BASE, "RNA-seq/results/audit_sensitivity/permutation_null")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

out_file <- file.path(OUT_DIR, sprintf("perm_iter_%04d.csv", PERM))
if (file.exists(out_file)) {
  cat("[Perm iter", PERM, "] Already done — skipping\n")
  quit(save = "no", status = 0)
}

dge <- readRDS(file.path(RDIR, "merged_dge.rds"))
ycfg <- yaml::read_yaml(file.path(BASE, "config/human_datasets.yaml"))$datasets
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

cat("Pre-shuffle group x dataset:\n"); print(table(info$group_binary, info$dataset))

# --- Within-cohort permutation of group_binary (preserves per-cohort counts) ---
set.seed(42 + PERM * 7919)
for (ds in unique(info$dataset)) {
  idx <- which(info$dataset == ds)
  info$group_binary[idx] <- sample(info$group_binary[idx])
}
cat("Post-shuffle group x dataset:\n"); print(table(info$group_binary, info$dataset))

ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
param <- if (ncpus > 1) MulticoreParam(ncpus) else SerialParam()
form <- ~ group_binary + inferred_sex + (1|dataset)
v <- suppressWarnings(voomWithDreamWeights(dge, form, info, BPPARAM = param))
fit <- suppressWarnings(dream(v, form, info, BPPARAM = param))
res <- topTable(fit, coef = "group_binaryDisease", number = Inf, sort.by = "none")
res$gene <- rownames(res); res_dt <- as.data.table(res); setnames(res_dt, "adj.P.Val", "padj")

# Compact output for permutation: gene/logFC/t/padj — aggregator turns this
# into the empirical null distribution per gene + global DEG-count distribution.
fwrite(res_dt[, .(gene, logFC, t, padj)],
       file.path(OUT_DIR, sprintf("perm_iter_%04d.csv", PERM)))
cat(sprintf("Elapsed %.1f min -- perm %d done\n", (proc.time()-t0)["elapsed"]/60, PERM))
