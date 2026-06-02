#!/usr/bin/env Rscript
# 05b_dream_ordinal.R
# ---------------------------------------------------------------------------
# Reviewer-defense (Agent B2, 2026-05-20): Parallel ordinal dream() model
# fitting fibrosis_stage as a numeric covariate plus per-transition contrasts
# (F0vF1, F1vF2, F2vF3, F3vF4) via a categorical fit.
#
# Two fits:
#   (1) Numeric (linear) trend: ~ as.numeric(fibrosis_stage) + inferred_sex
#                                 + (1|dataset)
#       -> coef "fib_num" = per-stage LFC slope
#
#   (2) Categorical (per-transition) fit: ~ fib_factor + inferred_sex
#                                          + (1|dataset)
#       with fib_factor = factor(fibrosis_stage, levels = c("0","1","2","3","4"))
#       -> 4 successive-stage contrasts via custom contrast matrix
#          (F1-F0, F2-F1, F3-F2, F4-F3)
#
# Output: results/integration/dream_results_ordinal.csv
#         (long format: gene, transition, logFC, t, P.Value, padj, SE)
# ---------------------------------------------------------------------------

set.seed(42)
Sys.setenv(R_PARALLEL_SEED = "42")

suppressPackageStartupMessages({
  library(reformulas)
  library(lme4)
  library(data.table)
  library(edgeR)
  library(yaml)
})

ns_lme4 <- asNamespace("lme4")
for (fn in c("findbars", "nobars", "subbars", "rebuildFormula")) {
  if (exists(fn, envir = ns_lme4)) {
    try({
      unlockBinding(fn, ns_lme4)
      assign(fn, get(fn, asNamespace("reformulas")), envir = ns_lme4)
      lockBinding(fn, ns_lme4)
    }, silent = TRUE)
  }
}

suppressPackageStartupMessages({
  library(variancePartition)
  library(BiocParallel)
})

INPUT_ROOT  <- Sys.getenv("MASLD_INPUT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
OUTPUT_ROOT <- Sys.getenv("MASLD_OUTPUT_ROOT", INPUT_ROOT)
IN_RDIR  <- file.path(INPUT_ROOT,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
OUT_RDIR <- file.path(OUTPUT_ROOT,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
dir.create(OUT_RDIR, recursive = TRUE, showWarnings = FALSE)
RDIR <- OUT_RDIR

cat("=== 05b_dream_ordinal.R (B2 reviewer defense) ===\n")
cat("Started:", as.character(Sys.time()), "\n\n")

# ---------------------------------------------------------------------------
# Load DGE + metadata; restrict to samples with non-missing fibrosis_stage
# ---------------------------------------------------------------------------
dge  <- readRDS(file.path(IN_RDIR, "merged_dge.rds"))
meta <- readRDS(file.path(IN_RDIR, "meta_matched.rds"))

ycfg <- yaml::read_yaml(file.path(INPUT_ROOT, "config/human_datasets.yaml"))$datasets
mega_cohorts <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))

keep_cohort <- dge$samples$dataset %in% mega_cohorts
dge2 <- dge[, keep_cohort]
cat("Samples after cohort filter:", ncol(dge2), "\n")

# Join fibrosis_stage + inferred_sex
m_match <- meta[match(colnames(dge2), meta$sample_id),
                .(sample_id, fibrosis_stage, inferred_sex, dataset)]
has_fib <- !is.na(m_match$fibrosis_stage)
cat("Samples with fibrosis_stage available:",
    sum(has_fib), "/", length(has_fib), "\n")

# Restrict to samples with fibrosis_stage in 0..4
fib_int <- suppressWarnings(as.integer(as.character(m_match$fibrosis_stage)))
ok <- !is.na(fib_int) & fib_int %in% 0:4
cat("Samples with fibrosis_stage in {0..4}:", sum(ok), "\n")

dge2 <- dge2[, ok]
m_match <- m_match[ok]
fib_int <- fib_int[ok]
cat("Fibrosis stage distribution:\n"); print(table(fib_int, useNA = "always"))

info <- data.frame(
  fib_num      = as.numeric(fib_int),
  fib_factor   = factor(as.character(fib_int), levels = c("0", "1", "2", "3", "4")),
  inferred_sex = factor(m_match$inferred_sex),
  dataset      = factor(m_match$dataset),
  stringsAsFactors = FALSE
)
rownames(info) <- colnames(dge2)
cat("Stage x dataset:\n"); print(table(info$fib_factor, info$dataset))

ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
cat("Using", ncpus, "CPU cores\n")
param <- if (ncpus > 1) MulticoreParam(ncpus) else SerialParam()

results_long <- list()

# ---------------------------------------------------------------------------
# Fit 1: Linear/numeric trend
# ---------------------------------------------------------------------------
cat("\n--- Fit 1: numeric fibrosis trend ---\n")
form_num <- ~ fib_num + inferred_sex + (1 | dataset)
cat("Formula:", deparse(form_num), "\n")

v1 <- suppressWarnings(voomWithDreamWeights(dge2, form_num, info, BPPARAM = param))
fit1 <- suppressWarnings(dream(v1, form_num, info, BPPARAM = param))
tt1 <- topTable(fit1, coef = "fib_num", number = Inf, sort.by = "none")
tt1$gene <- rownames(tt1)
dt1 <- as.data.table(tt1)
setnames(dt1, "adj.P.Val", "padj")
dt1[, SE := ifelse(is.finite(t) & t != 0, logFC / t, NA_real_)]
dt1[, transition := "linear_trend"]
results_long[["linear"]] <- dt1[, .(gene, transition, logFC, AveExpr, t,
                                     P.Value, padj, SE)]
cat("Linear trend DEGs (padj<0.1):", sum(dt1$padj < 0.1, na.rm = TRUE), "\n")

# ---------------------------------------------------------------------------
# Fit 2: Categorical, extract successive-stage contrasts
# ---------------------------------------------------------------------------
cat("\n--- Fit 2: categorical per-stage contrasts ---\n")
form_cat <- ~ 0 + fib_factor + inferred_sex + (1 | dataset)
cat("Formula:", deparse(form_cat), "\n")

# Build contrasts for F1-F0, F2-F1, F3-F2, F4-F3
transitions <- list(
  F0vF1 = c("fib_factor1", "fib_factor0"),
  F1vF2 = c("fib_factor2", "fib_factor1"),
  F2vF3 = c("fib_factor3", "fib_factor2"),
  F3vF4 = c("fib_factor4", "fib_factor3")
)

# Only keep transitions where both stages exist
avail_levels <- levels(droplevels(info$fib_factor))
have_levels <- paste0("fib_factor", avail_levels)

v2 <- suppressWarnings(voomWithDreamWeights(dge2, form_cat, info, BPPARAM = param))

for (tname in names(transitions)) {
  hi <- transitions[[tname]][1]; lo <- transitions[[tname]][2]
  if (!all(c(hi, lo) %in% have_levels)) {
    cat("[skip]", tname, "— missing one of", hi, "/", lo, "\n")
    next
  }
  L <- variancePartition::makeContrastsDream(form_cat, info,
         contrasts = setNames(paste0(hi, " - ", lo), tname))
  cat("Running dream() for", tname, "...\n")
  fit2 <- suppressWarnings(dream(v2, form_cat, info, L = L, BPPARAM = param))
  tt2 <- topTable(fit2, coef = tname, number = Inf, sort.by = "none")
  tt2$gene <- rownames(tt2)
  dt2 <- as.data.table(tt2)
  setnames(dt2, "adj.P.Val", "padj")
  dt2[, SE := ifelse(is.finite(t) & t != 0, logFC / t, NA_real_)]
  dt2[, transition := tname]
  results_long[[tname]] <- dt2[, .(gene, transition, logFC, AveExpr, t,
                                    P.Value, padj, SE)]
  cat("  DEGs (padj<0.1):", sum(dt2$padj < 0.1, na.rm = TRUE), "\n")
}

# ---------------------------------------------------------------------------
# Combine + save
# ---------------------------------------------------------------------------
out <- rbindlist(results_long, use.names = TRUE, fill = TRUE)
out_csv <- file.path(RDIR, "dream_results_ordinal.csv")
fwrite(out, out_csv)
cat("\nSaved:", out_csv, "\n")
cat("Rows:", nrow(out), "(", length(unique(out$gene)), "unique genes x",
    length(unique(out$transition)), "contrasts)\n")
cat("Finished:", as.character(Sys.time()), "\n")
