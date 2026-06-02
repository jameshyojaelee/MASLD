#!/usr/bin/env Rscript
# dream_lto_iter.R
# ---------------------------------------------------------------------------
# Leave-TWO-cohorts-out cross-validation for the dream mega-analysis.
# Cohort-level stress test at ~60 % retention (3 of 5 cohorts trained on).
# Complements LOCO (5 folds at 80 %) on a different axis from patient-level
# CPSS / bootstrap (which subsample patients within all cohorts).
#
# 10 pairs total = C(5, 2). Pair index 1..10 selected by LTO_PAIR_IDX env var.
# Output: results/integration/loo_cv/dream_lto_pair{NN}_held_{C1}_{C2}.csv
# ---------------------------------------------------------------------------

PAIR <- as.integer(Sys.getenv("LTO_PAIR_IDX", "1"))
if (is.na(PAIR) || PAIR < 1 || PAIR > 10) stop("LTO_PAIR_IDX must be 1..10")
cat(sprintf("=== Dream LTO: pair %d/10 ===\n", PAIR))

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

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts"
INT  <- file.path(BASE, "analysis/integration")
RDIR <- file.path(INT, "results/integration")
LOO_DIR <- file.path(RDIR, "loo_cv")
dir.create(LOO_DIR, recursive = TRUE, showWarnings = FALSE)

# --- Load + identify mega cohorts ---
dge <- readRDS(file.path(RDIR, "merged_dge.rds"))
ycfg <- yaml::read_yaml(file.path(
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
  "config/human_datasets.yaml"))$datasets
mega_cohorts <- sort(names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg)))
cat("Mega cohorts (k=", length(mega_cohorts), "): ", paste(mega_cohorts, collapse=", "), "\n", sep="")

# Deterministic enumeration of C(5,2)=10 pairs (alphabetically sorted)
pairs_mat <- t(utils::combn(mega_cohorts, 2))
held_out <- as.character(pairs_mat[PAIR, ])
keep <- setdiff(mega_cohorts, held_out)
cat(sprintf("Pair %d: hold out  {%s, %s}\n", PAIR, held_out[1], held_out[2]))
cat(sprintf("        train on  {%s}\n", paste(keep, collapse=", ")))

# --- Subset DGE ---
keep_samples <- dge$samples$dataset %in% keep
dge_lto <- dge[, keep_samples]
cat("Samples remaining:", ncol(dge_lto), "/", ncol(dge), "\n")

# --- Update inferred sex ---
meta_new <- readRDS(file.path(RDIR, "meta_matched.rds"))
matched_sex <- meta_new$inferred_sex[match(colnames(dge_lto), meta_new$sample_id)]
info <- data.frame(
  group_binary = factor(dge_lto$samples$group_binary, levels = c("Control", "Disease")),
  dataset      = droplevels(factor(dge_lto$samples$dataset)),
  inferred_sex = factor(matched_sex),
  stringsAsFactors = FALSE)
rownames(info) <- colnames(dge_lto)
na_sex <- is.na(info$inferred_sex)
if (any(na_sex)) {
  warning(sum(na_sex), " samples have NA inferred_sex; dropping.")
  dge_lto <- dge_lto[, !na_sex]; info <- info[!na_sex, , drop = FALSE]
  info$dataset <- droplevels(info$dataset)
}
cat("Group x dataset:\n"); print(table(info$group_binary, info$dataset))

# Drop datasets with single group present in this LTO subset
ds_ok <- names(which(apply(table(info$dataset, info$group_binary), 1, function(x) all(x > 0))))
if (length(ds_ok) < 2) {
  cat("WARN: <2 datasets retain both groups in LTO pair", PAIR, "\n")
  fwrite(data.table(gene = character(0), logFC = numeric(0), t = numeric(0), padj = numeric(0)),
         file.path(LOO_DIR, sprintf("dream_lto_pair%02d_held_%s_%s.csv", PAIR, held_out[1], held_out[2])))
  quit(save = "no", status = 0)
}
keep_idx <- info$dataset %in% ds_ok
if (!all(keep_idx)) {
  dge_lto <- dge_lto[, keep_idx]; info <- info[keep_idx, , drop = FALSE]
  info$dataset <- droplevels(info$dataset)
}

# --- dream() ---
form <- ~ group_binary + inferred_sex + (1 | dataset)
ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
param <- if (ncpus > 1) MulticoreParam(ncpus) else SerialParam()
v <- suppressWarnings(voomWithDreamWeights(dge_lto, form, info, BPPARAM = param))
fit <- suppressWarnings(dream(v, form, info, BPPARAM = param))
res <- topTable(fit, coef = "group_binaryDisease", number = Inf, sort.by = "none")
res$gene <- rownames(res); res_dt <- as.data.table(res); setnames(res_dt, "adj.P.Val", "padj")

out_file <- file.path(LOO_DIR, sprintf("dream_lto_pair%02d_held_%s_%s.csv",
                                        PAIR, held_out[1], held_out[2]))
fwrite(res_dt[, .(gene, logFC, t, padj)], out_file)
cat(sprintf("Saved: %s\n", basename(out_file)))
cat(sprintf("DEGs (padj<0.05 & |LFC|>0.5): %d\n",
            sum(res_dt$padj < 0.05 & abs(res_dt$logFC) > 0.5, na.rm = TRUE)))
cat("Done LTO pair", PAIR, "\n")
