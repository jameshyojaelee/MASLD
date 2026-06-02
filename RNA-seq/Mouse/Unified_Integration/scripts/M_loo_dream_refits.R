#!/usr/bin/env Rscript
# M_loo_dream_refits.R  (new 2026-04-27)
# ---------------------------------------------------------------------------
# Mouse leave-one-diet-model-out (LOO) refits of the pooled dream mega-analysis.
# For each of the 5 diet models (MCD, HFD, CDAHFD, FPC, LIDPAD), refit
# voomWithDreamWeights + dream() on the remaining 4 diets and save the full
# DE table. Used by M_loo_lfc_stability.R to derive the empirically defensible
# |log2FC| cutoff (mouse analogue of human Fig 1G stability analysis).
#
# Input:  results/merged_counts_raw.rds, results/meta_matched.rds,
#         qc/sample_qc_report.csv
# Output: results/integration/loo_cv/dream_loo_<diet>.csv      (5 files)
#
# Mirrors the human LOO refit logic but uses diet_model as the held-out unit
# (5 biological replicates of MASLD induction; cleanest unit for mouse).
# ---------------------------------------------------------------------------

# ---- Seed pinning ----------------------------------------------------------
set.seed(42)
Sys.setenv(R_PARALLEL_SEED = "42")

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(variancePartition)
  library(BiocParallel)
  library(yaml)
})

BASE   <- Sys.getenv("MASLD_PROJECT_ROOT",
                     "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
MOUSE  <- file.path(BASE, "RNA-seq/Mouse")
INT    <- file.path(MOUSE, "Unified_Integration")
RDIR   <- file.path(INT, "results")
LOO    <- file.path(RDIR, "integration", "loo_cv")
dir.create(LOO, recursive = TRUE, showWarnings = FALSE)

# Auto-discover diet models from YAML config (authoritative registry)
cfg_path <- file.path(BASE, "config/mouse_datasets.yaml")
if (!file.exists(cfg_path)) {
  stop("Missing config: ", cfg_path)
}
mouse_cfg <- yaml::read_yaml(cfg_path)
DIET_MODELS <- mouse_cfg$diet_models
if (is.null(DIET_MODELS) || length(DIET_MODELS) == 0L) {
  stop("No diet_models list found in ", cfg_path)
}
cat("Diet models from YAML config:", paste(DIET_MODELS, collapse = ", "), "\n")

# Optional: focus on a single LOO fold via env var (for sbatch parallelism)
loo_target <- Sys.getenv("LOO_DIET", "")
if (nzchar(loo_target)) {
  stopifnot(loo_target %in% DIET_MODELS)
  DIET_MODELS <- loo_target
  cat(sprintf("LOO_DIET set -> running single fold: %s\n", loo_target))
}

cat("=== M_loo_dream_refits: leave-one-diet-out dream refits ===\n\n")

# ---------------------------------------------------------------------------
# Load data once (shared across folds)
# ---------------------------------------------------------------------------
merged <- readRDS(file.path(RDIR, "merged_counts_raw.rds"))
meta   <- readRDS(file.path(RDIR, "meta_matched.rds"))
qc     <- fread(file.path(INT, "qc/sample_qc_report.csv"))

pass <- qc[pass_qc == TRUE, sample_id]
merged <- merged[, colnames(merged) %in% pass]
meta   <- meta[sample_id %in% pass]

cat(sprintf("Total QC-passing samples: %d across %d datasets, %d diet models\n",
            ncol(merged),
            length(unique(meta$dataset)),
            length(unique(meta$diet_model))))
cat("\nDiet model distribution (Disease vs Control):\n")
print(table(meta$diet_model, meta$group_binary))
cat("\n")

# Mouse meta uses the same field names as human dream
n_cores <- min(parallelly::availableCores(), 16)
cat(sprintf("Using %d CPU cores\n\n", n_cores))

run_dream_fold <- function(held_out) {
  cat(sprintf("============================================================\n"))
  cat(sprintf("  LOO fold: holding out diet model = %s\n", held_out))
  cat(sprintf("============================================================\n"))

  # Hold out: drop ALL samples from the held-out diet (Disease + Control of
  # that diet, since controls were drawn from the same datasets).
  # Keep all samples whose diet_model != held_out.
  keep_idx <- meta$diet_model != held_out
  fold_meta   <- meta[keep_idx]
  fold_counts <- merged[, fold_meta$sample_id]

  cat(sprintf("  Samples retained: %d (held out %d %s samples)\n",
              nrow(fold_meta), sum(!keep_idx), held_out))
  cat(sprintf("    Disease: %d, Control: %d\n",
              sum(fold_meta$group_binary == "Disease"),
              sum(fold_meta$group_binary == "Control")))
  cat(sprintf("    Datasets remaining: %d\n",
              length(unique(fold_meta$dataset))))

  if (sum(fold_meta$group_binary == "Disease") < 2 ||
      sum(fold_meta$group_binary == "Control") < 2) {
    stop("Insufficient samples in fold")
  }

  fold_meta$group_binary <- factor(fold_meta$group_binary,
                                   levels = c("Control", "Disease"))
  fold_meta$dataset <- factor(fold_meta$dataset)

  # filterByExpr on the fold (matches M03 pooled behavior)
  dge <- DGEList(counts = fold_counts)
  keep_g <- filterByExpr(dge, group = fold_meta$group_binary)
  dge <- dge[keep_g, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)
  cat(sprintf("  Genes after filterByExpr: %d\n", nrow(dge)))

  form <- ~ group_binary + (1 | dataset)
  param <- SnowParam(n_cores, "SOCK", progressbar = FALSE)

  cat("  Running voomWithDreamWeights...\n")
  vobj <- voomWithDreamWeights(dge, form, fold_meta, BPPARAM = param)

  cat("  Running dream()...\n")
  fit <- dream(vobj, form, fold_meta, BPPARAM = param)
  fit <- eBayes(fit)

  res <- topTable(fit, coef = "group_binaryDisease",
                  number = Inf, sort.by = "none")
  res$gene <- rownames(res)
  res <- as.data.table(res)
  setcolorder(res, "gene")

  sig <- res[adj.P.Val < 0.05]
  cat(sprintf("  Fold DEGs (padj<0.05): %d\n", nrow(sig)))
  cat(sprintf("    |LFC|>0.5: %d   |LFC|>1.0: %d   |LFC|>1.5: %d\n",
              sum(sig$adj.P.Val < 0.05 & abs(sig$logFC) > 0.5),
              sum(sig$adj.P.Val < 0.05 & abs(sig$logFC) > 1.0),
              sum(sig$adj.P.Val < 0.05 & abs(sig$logFC) > 1.5)))

  out_file <- file.path(LOO, sprintf("dream_loo_%s.csv", held_out))
  fwrite(res[order(adj.P.Val)], out_file)
  cat(sprintf("  Saved: %s\n\n", out_file))

  invisible(NULL)
}

for (dm in DIET_MODELS) {
  run_dream_fold(dm)
}

cat("M_loo_dream_refits complete.\n")
