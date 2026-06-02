#!/usr/bin/env Rscript
# dream_cpss_iter.R
# ---------------------------------------------------------------------------
# Pillar A — Complementary-Pairs Stability Selection (Shah & Samworth 2013).
# ONE pair = a 50% half-sample I plus its complement I^c, both fitted with dream.
# Cohort-stratified subsampling preserves per-cohort disease/control proportions.
#
# Env: CPSS_PAIR  (1..100) — pair index. B=100 pairs => 200 dream fits total.
# Out: results/integration/robustness/cpss/iter_{PAIR}_half{A,B}.csv
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

PAIR <- as.integer(Sys.getenv("CPSS_PAIR", "1"))
if (is.na(PAIR) || PAIR < 1) stop("CPSS_PAIR must be >=1")
cat("=== CPSS pair", PAIR, "===\n")

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts"
INT  <- file.path(BASE, "analysis/integration")
RDIR <- file.path(INT, "results/integration")
OUT_DIR <- file.path(RDIR, "robustness/cpss")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

label3 <- sprintf("%03d", PAIR)
out_a <- file.path(OUT_DIR, paste0("iter_", label3, "_halfA.csv"))
if (file.exists(out_a)) {
  cat("[CPSS pair", PAIR, "] Already done — skipping\n")
  quit(save = "no", status = 0)
}

# --- Load + filter to mega cohorts ---
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

# Drop NA-sex (matches dream_loo_cv.R behaviour)
na_sex <- is.na(info$inferred_sex)
if (any(na_sex)) { dge <- dge[, !na_sex]; info <- info[!na_sex, , drop = FALSE] }

cat("Total samples:", ncol(dge), "across", nlevels(info$dataset), "cohorts\n")
cat("Group x dataset:\n"); print(table(info$group_binary, info$dataset))

# --- Cohort-stratified 50/50 split, seeded by PAIR ---
set.seed(42 + PAIR * 7919)
half_A <- integer(0)
for (ds in unique(info$dataset)) {
  ds_idx <- which(info$dataset == ds)
  n_half <- floor(length(ds_idx) / 2)
  half_A <- c(half_A, sample(ds_idx, n_half))
}
half_B <- setdiff(seq_len(nrow(info)), half_A)
cat(sprintf("Half-A: n=%d  Half-B: n=%d\n", length(half_A), length(half_B)))

# --- dream fitting helper ---
ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
param <- if (ncpus > 1) MulticoreParam(ncpus) else SerialParam()
form  <- ~ group_binary + inferred_sex + (1|dataset)

fit_one_half <- function(idx, label) {
  cat("\n--- Fitting half", label, "(n=", length(idx), ") ---\n")
  d <- dge[, idx]; i <- info[idx, , drop = FALSE]
  i$dataset <- droplevels(i$dataset)
  # Drop datasets with single group within this half
  ds_ok <- names(which(apply(table(i$dataset, i$group_binary), 1, function(x) all(x > 0))))
  if (length(ds_ok) < 2) {
    cat("  WARN: <2 datasets retain both groups in half ", label, "; writing NA-only output\n")
    return(NULL)
  }
  keep <- i$dataset %in% ds_ok
  if (!all(keep)) { d <- d[, keep]; i <- i[keep, , drop = FALSE]; i$dataset <- droplevels(i$dataset) }
  v <- suppressWarnings(voomWithDreamWeights(d, form, i, BPPARAM = param))
  fit <- tryCatch(
    suppressWarnings(dream(v, form, i, BPPARAM = param)),
    error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL })
  if (is.null(fit)) return(NULL)
  res <- topTable(fit, coef = "group_binaryDisease", number = Inf, sort.by = "none")
  res$gene <- rownames(res)
  res_dt <- as.data.table(res); setnames(res_dt, "adj.P.Val", "padj")
  res_dt[, .(gene, logFC, t, padj)]
}

resA <- fit_one_half(half_A, "A")
resB <- fit_one_half(half_B, "B")

label3 <- sprintf("%03d", PAIR)
if (!is.null(resA)) fwrite(resA, file.path(OUT_DIR, paste0("iter_", label3, "_halfA.csv")))
if (!is.null(resB)) fwrite(resB, file.path(OUT_DIR, paste0("iter_", label3, "_halfB.csv")))

cat(sprintf("\nElapsed %.1f min\n", (proc.time()-t0)["elapsed"]/60))
cat("Done CPSS pair", PAIR, "\n")
