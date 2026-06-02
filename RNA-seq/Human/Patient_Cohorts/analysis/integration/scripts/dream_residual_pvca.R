#!/usr/bin/env Rscript
# dream_residual_pvca.R
# ---------------------------------------------------------------------------
# Pillar D — variancePartition on dream-equivalent mixed-effects residuals.
# variancePartition 1.36.3 does NOT auto-compute fit$residuals for mixed
# models, so we refit lmer per gene in parallel and extract residuals
# directly. Then re-run variancePartition on the residual matrix.
# Pass: cohort variance in residuals < 5% (vs raw 76.5% in 04_variance_partition).
#
# Out: results/audit_sensitivity/pillar_D_residual_varpart.csv         (per-gene)
#      results/audit_sensitivity/pillar_D_residual_varpart_summary.csv (1 row per component)
# ---------------------------------------------------------------------------

t0 <- proc.time()
suppressPackageStartupMessages({
  library(reformulas); library(lme4); library(data.table); library(edgeR); library(yaml); library(lmerTest)
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

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
INT  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR <- file.path(INT, "results/integration")
OUT_DIR <- file.path(BASE, "RNA-seq/results/audit_sensitivity")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

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

ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
param <- if (ncpus > 1) MulticoreParam(ncpus) else SerialParam()
register(param)

cat("Computing voom weights with dream model...\n")
form_dream <- ~ group_binary + inferred_sex + (1|dataset)
v <- suppressWarnings(voomWithDreamWeights(dge, form_dream, info, BPPARAM = param))

# --- Per-gene mixed-effects residuals via lme4 (parallel) ---
# 34k genes ÷ 16 CPU × ~1s/fit ≈ 35 min. Each gene fit is independent.
cat("Computing mixed-effects residuals per gene (this is the slow step)...\n")
gene_idx <- seq_len(nrow(v$E))
compute_resid_one <- function(i) {
  df <- info
  df$y <- v$E[i, ]
  df$w <- v$weights[i, ]
  fit <- tryCatch(
    suppressMessages(suppressWarnings(
      lme4::lmer(y ~ group_binary + inferred_sex + (1 | dataset),
                 data = df, weights = w, REML = TRUE,
                 control = lme4::lmerControl(
                   check.conv.singular = "ignore",
                   check.conv.grad     = "ignore",
                   check.nobs.vs.nlev  = "ignore",
                   check.nobs.vs.nRE   = "ignore")))),
    error = function(e) NULL)
  if (is.null(fit)) return(rep(NA_real_, ncol(v$E)))
  as.numeric(residuals(fit))
}
# Chunk-wise reporting via per-chunk timing
chunk_size <- 500L
n_genes <- length(gene_idx)
n_chunks <- ceiling(n_genes / chunk_size)
resid_mat <- matrix(NA_real_, nrow = n_genes, ncol = ncol(v$E),
                    dimnames = list(rownames(v$E), colnames(v$E)))
for (c_i in seq_len(n_chunks)) {
  cstart <- (c_i - 1) * chunk_size + 1L
  cend   <- min(c_i * chunk_size, n_genes)
  ix <- gene_idx[cstart:cend]
  t_chunk <- proc.time()
  out <- bplapply(ix, compute_resid_one, BPPARAM = param)
  out <- do.call(rbind, out)
  resid_mat[cstart:cend, ] <- out
  el <- (proc.time() - t_chunk)["elapsed"]
  cat(sprintf("  chunk %d/%d (%d–%d): %.1f s\n", c_i, n_chunks, cstart, cend, el))
}

n_failed <- sum(is.na(resid_mat[, 1]))
cat(sprintf("Genes with failed mixed fit (NA-only residuals): %d / %d\n", n_failed, n_genes))
if (n_failed > 0) {
  cat("Dropping failed genes for variancePartition...\n")
  resid_mat <- resid_mat[!is.na(resid_mat[, 1]), , drop = FALSE]
}

cat("Variance partition on residuals...\n")
form_vp <- ~ (1|group_binary) + (1|inferred_sex) + (1|dataset)
vp <- suppressWarnings(fitExtractVarPartModel(resid_mat, form_vp, info))
vp_dt <- as.data.table(vp, keep.rownames = "gene")
fwrite(vp_dt, file.path(OUT_DIR, "pillar_D_residual_varpart.csv"))

med <- as.list(apply(vp, 2, median))
summary_dt <- data.table(component = names(med), median_var = unlist(med))
summary_dt[, mean_var := apply(vp, 2, mean)[component]]
summary_dt[, q90_var  := apply(vp, 2, function(x) quantile(x, .9, na.rm = TRUE))[component]]
summary_dt[, max_var  := apply(vp, 2, max)[component]]
summary_dt[, n_genes  := nrow(vp)]
fwrite(summary_dt, file.path(OUT_DIR, "pillar_D_residual_varpart_summary.csv"))
cat("\nMedian residual variance fractions:\n"); print(summary_dt)
ds_med <- summary_dt[component == "dataset", median_var]
cat(sprintf("\nPass (cohort residual var < 0.05): %s   (median = %.4f)\n",
            ifelse(length(ds_med) && !is.na(ds_med) && ds_med < 0.05, "PASS", "FAIL"), ds_med))
cat(sprintf("Elapsed %.1f min  -- residual PVCA done\n", (proc.time() - t0)["elapsed"]/60))
