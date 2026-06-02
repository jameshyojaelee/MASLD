#!/usr/bin/env Rscript
# ============================================================================
# 07b_chain_v3_hardened.R
#
# Hardened v3 of 07_chain_346_349_v2.R. Addresses critique #1, #2, #4, #7:
#
#   #1 (continuous axis empty): the v2 chain joined `macrophage_pseudotime_mean`
#      from donor_metadata_v2.tsv (which does NOT contain the column). The
#      column lives in donor_metadata_extended.tsv at:
#      results_gpu_v2/ccc/stage_trajectory/donor_metadata_extended.tsv.
#      Fix: join from BOTH files; the pseudotime column comes from extended.
#
#   #2 (donor imbalance / micro-cohorts): GSE174748 (n=4) and GSE189600 (n=2)
#      cannot estimate a variance component in the (1|dataset) random effect.
#      Fix: collapse both into `dataset_pooled = "Other_small"`; six dataset
#      levels remain in the random effect, each with n >= 20.
#
#   #4 (F-stage circularity): F_stage_augmented_v2 is kNN-predicted from
#      obsm['X_scVI']. Run BOTH:
#        - PRIMARY F-stage LMM on F_stage_documented (n=58 Andrews, circularity-
#          free) -> stage_lr_lmm_fstage_documented_v3.tsv
#        - SENSITIVITY F-stage LMM on F_stage_augmented_v2 (n up to 269) ->
#          stage_lr_lmm_fstage_augmented_v3.tsv
#      Coarse axis (disease_stage_coarse) remains primary cross-axis claim --
#      it is NOT derived from the scVI latent.
#
#   #7 (hierarchical Bonferroni not applied): add `q_bonferroni_family` =
#      `q_within_ct` * 72 (capped at 1) per Meinshausen 2008 hierarchical FDR.
#
# Inputs (read-only):
#   - results_gpu_v2_phase05/ccc/stage_trajectory_v2/per_donor_lr_v2/*.parquet
#     (merged through MERGED_TSV cache)
#   - results_gpu_v2_phase05/mcp/inputs/donor_metadata_v2.tsv
#   - results_gpu_v2/ccc/stage_trajectory/donor_metadata_extended.tsv
#   - RNA-seq/Human/.../dream_results_ashr.csv (bulk anchor)
#   - RNA-seq/Human/.../dream_results_stage_*.csv (per-stage bulk)
#
# Outputs (all under results_gpu_v2_phase05/ccc/stage_trajectory_v3):
#   - stage_lr_lmm_coarse_v3.tsv
#   - stage_lr_lmm_fstage_documented_v3.tsv   (primary, circularity-free)
#   - stage_lr_lmm_fstage_augmented_v3.tsv    (sensitivity)
#   - stage_lr_lmm_continuous_v3.tsv          (no longer empty)
#   - loo_replication_rate_per_lr_v3.tsv
#   - lr_bulk_concordance_v3.tsv
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(lme4)
  library(lmerTest)
  library(parallel)
  library(ggplot2)
  library(patchwork)
})

BASE <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)

N_CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "16"))
options(mc.cores = N_CORES)
cat(sprintf("[parallel] using %d cores\n", N_CORES))

V2_STAGE_DIR <- file.path(
  BASE, "Analysis/SingleCell/results_gpu_v2_phase05/ccc/stage_trajectory_v2"
)
V3_STAGE_DIR <- file.path(
  BASE, "Analysis/SingleCell/results_gpu_v2_phase05/ccc/stage_trajectory_v3"
)
dir.create(V3_STAGE_DIR, recursive = TRUE, showWarnings = FALSE)

PER_DONOR_DIR <- file.path(V2_STAGE_DIR, "per_donor_lr_v2")
META_V2 <- file.path(
  BASE,
  "Analysis/SingleCell/results_gpu_v2_phase05/mcp/inputs/donor_metadata_v2.tsv"
)
META_EXTENDED_V1 <- file.path(
  BASE,
  "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/donor_metadata_extended.tsv"
)
MERGED_TSV <- file.path(V2_STAGE_DIR, "all_donor_lr_scores_v2.tsv.gz")

INT_DIR <- file.path(
  BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration"
)
BULK_ALL <- file.path(INT_DIR, "dream_results_ashr.csv")
BULK_STEATOSIS <- file.path(INT_DIR, "dream_results_stage_steatosis.csv")
BULK_SH <- file.path(INT_DIR, "dream_results_stage_sh.csv")
BULK_CIRRHOSIS <- file.path(INT_DIR, "dream_results_stage_cirrhosis.csv")

MIN_DONORS_PER_LR <- 30
MIN_DONORS_PER_STRATUM <- 10
TOP_N <- 100

# Critique #2: micro-cohort pooling. GSE174748 n=4 and GSE189600 n=2 cannot
# support a variance component; collapse them into a single "Other_small"
# dataset level for the random effect.
SMALL_COHORTS <- c("GSE174748", "GSE189600")

# Critique #7 / Master review M-P0-4: hierarchical Bonferroni family count is
# per-axis ADAPTIVE (number of unique ct_pair groups in the fitted result),
# not a hardcoded 9*8=72. Actual axis counts: coarse=34, documented=21,
# augmented=55, continuous=63 — hardcoded 72 over-corrected by 1.1x-3.4x.
# N_CT_PAIR_FAMILIES retained as fallback default but should not be used.
N_CT_PAIR_FAMILIES_FALLBACK <- 72L

# ----------------------------------------------------------------------------
# 0. Consolidate per-donor parquets -> merged TSV (avoids arrow R dep)
# ----------------------------------------------------------------------------
if (!file.exists(MERGED_TSV)) {
  cat(sprintf("[merge] consolidating parquets in %s\n", PER_DONOR_DIR))
  parquets <- list.files(PER_DONOR_DIR, pattern = "_lr_scores\\.parquet$",
                         full.names = TRUE)
  if (length(parquets) == 0) {
    stop("No per-donor parquets found in ", PER_DONOR_DIR)
  }
  merge_py <- tempfile(fileext = ".py")
  writeLines(c(
    "import os, sys",
    "import pandas as pd",
    "from pathlib import Path",
    sprintf("PER_DONOR = Path(%s)", shQuote(PER_DONOR_DIR)),
    sprintf("OUT = Path(%s)",        shQuote(MERGED_TSV)),
    "parquets = sorted(PER_DONOR.glob('*_lr_scores.parquet'))",
    "dfs = [pd.read_parquet(p) for p in parquets]",
    "merged = pd.concat(dfs, ignore_index=True)",
    "print(f'shape={merged.shape}')",
    "merged.to_csv(OUT, sep='\\t', index=False, compression='gzip')"
  ), merge_py)
  cmd <- sprintf("micromamba run -n rapids_singlecell python -u %s",
                 shQuote(merge_py))
  res <- system(cmd, intern = FALSE)
  if (res != 0 || !file.exists(MERGED_TSV)) {
    stop("Merge step failed (code=", res, ")")
  }
}

lr_long <- fread(MERGED_TSV)
cat(sprintf("[input] %d donor x LR rows from merged v2 TSV\n", nrow(lr_long)))

# ----------------------------------------------------------------------------
# Critique #1 FIX: join from BOTH donor metadata files
#   - donor_metadata_v2.tsv: F_stage_{documented,augmented_v2,inferred_v2},
#       disease_stage_coarse, age, sex_numeric
#   - donor_metadata_extended.tsv: macrophage_pseudotime_mean (col 34),
#       hepatocyte_pseudotime_mean (col 36)
# ----------------------------------------------------------------------------
meta_v2 <- fread(META_V2)
cat(sprintf("[input] %d donors in donor_metadata_v2.tsv (%d cols)\n",
            nrow(meta_v2), ncol(meta_v2)))

if (file.exists(META_EXTENDED_V1)) {
  meta_ext <- fread(META_EXTENDED_V1)
  cat(sprintf("[input] %d donors in donor_metadata_extended.tsv (%d cols)\n",
              nrow(meta_ext), ncol(meta_ext)))
  # Pull just sample + pseudotime columns; join on sample
  pt_cols <- intersect(
    c("sample", "macrophage_pseudotime_mean", "hepatocyte_pseudotime_mean"),
    names(meta_ext)
  )
  if (length(pt_cols) >= 2) {
    meta_pt <- meta_ext[, ..pt_cols]
    meta_v2 <- merge(meta_v2, meta_pt, by = "sample", all.x = TRUE)
    cat(sprintf("[meta] joined pseudotime columns: %s\n",
                paste(setdiff(pt_cols, "sample"), collapse = ", ")))
    cat(sprintf("[meta] non-NA macrophage_pseudotime_mean: %d / %d donors\n",
                sum(!is.na(meta_v2$macrophage_pseudotime_mean)), nrow(meta_v2)))
  } else {
    cat("[warn] donor_metadata_extended.tsv missing pseudotime cols\n")
  }
} else {
  cat(sprintf("[warn] donor_metadata_extended.tsv not found at %s\n",
              META_EXTENDED_V1))
}

# Master review M-P0-17: age + sex covariates are recoverable. Sex from
# snrna_donor_sex.csv (XIST/DDX3Y k-means inference). Override / fill
# whatever upstream metadata wrote in.
SNRNA_SEX <- file.path(BASE,
  "Analysis/SingleCell/metadata/snrna_donor_sex.csv")
if (file.exists(SNRNA_SEX)) {
  sex_overlay <- fread(SNRNA_SEX)
  if (all(c("sample", "inferred_sex") %in% names(sex_overlay))) {
    sex_overlay[, sex_numeric_overlay := fifelse(inferred_sex == "Female", 1L,
                                          fifelse(inferred_sex == "Male", 0L, NA_integer_))]
    meta_v2 <- merge(meta_v2,
                     sex_overlay[, .(sample, sex_numeric_overlay)],
                     by = "sample", all.x = TRUE)
    if (!"sex_numeric" %in% names(meta_v2)) meta_v2[, sex_numeric := NA_integer_]
    n_before <- sum(!is.na(meta_v2$sex_numeric))
    meta_v2[is.na(sex_numeric), sex_numeric := sex_numeric_overlay]
    meta_v2[, sex_numeric_overlay := NULL]
    n_after <- sum(!is.na(meta_v2$sex_numeric))
    cat(sprintf("[meta] sex_numeric filled from snrna_donor_sex.csv: %d -> %d donors\n",
                n_before, n_after))
  }
} else {
  cat(sprintf("[warn] %s missing -> sex covariate will be 100%% NA\n", SNRNA_SEX))
}

# Critique #2 FIX: pool micro-cohorts into dataset_pooled for random effect.
# Keep original `dataset` for diagnostics / cross-axis breakdown.
meta_v2[, dataset_pooled := ifelse(dataset %in% SMALL_COHORTS,
                                   "Other_small",
                                   as.character(dataset))]
cat("[meta] dataset_pooled breakdown:\n")
print(meta_v2[, .N, by = dataset_pooled])

covar_cols <- intersect(
  c("sample", "dataset", "dataset_pooled",
    "disease_stage_coarse", "disease_stage_numeric",
    "F_stage_documented", "F_stage_inferred_v2", "F_stage_augmented_v2",
    "F_stage_source_v2", "macrophage_pseudotime_mean",
    "hepatocyte_pseudotime_mean", "age", "sex_numeric"),
  names(meta_v2)
)
lr_long <- merge(lr_long, meta_v2[, ..covar_cols], by = "sample", all.x = TRUE)

lr_long[, score := -log10(pmax(magnitude_rank, 1e-4))]
stage_levels <- c("Healthy", "Steatosis", "Steatohepatitis", "Cirrhosis")
lr_long[, disease_stage_coarse := factor(disease_stage_coarse, levels = stage_levels)]
lr_long <- lr_long[!is.na(disease_stage_coarse)]
lr_long[, lr_pair := paste(ligand_complex, receptor_complex, sep = "__")]
lr_long[, ct_pair := paste(source, target, sep = "->")]

# ----------------------------------------------------------------------------
# Helper: fit one model. Uses dataset_pooled for the random effect (#2 fix).
# Returns NULL if model fails or axis term has no coefficients.
# ----------------------------------------------------------------------------
fit_one <- function(d, axis_term) {
  if (nrow(d) < MIN_DONORS_PER_LR) return(NULL)
  n_pools <- length(unique(d$dataset_pooled))
  use_random <- n_pools >= 3
  fixed <- c(axis_term)
  # Master review M-P0-3: dataset_pooled must NOT be both fixed and random.
  # If random effect is used, dataset_pooled is absorbed there.
  if (!use_random && n_pools > 1) fixed <- c(fixed, "dataset_pooled")
  fixed <- c(fixed,
             "log10(pmax(n_source_cells,1))",
             "log10(pmax(n_target_cells,1))")
  if ("age" %in% names(d) &&
      length(unique(stats::na.omit(d$age))) > 1 &&
      sum(!is.na(d$age)) >= 10) {
    fixed <- c(fixed, "age")
  }
  if ("sex_numeric" %in% names(d) &&
      length(unique(stats::na.omit(d$sex_numeric))) > 1 &&
      sum(!is.na(d$sex_numeric)) >= 10) {
    fixed <- c(fixed, "sex_numeric")
  }
  rhs <- paste(fixed, collapse = " + ")
  if (use_random) {
    rhs <- paste(rhs, "+ (1 | dataset_pooled)")
    fit <- try(lmerTest::lmer(stats::as.formula(paste("score ~", rhs)),
                              data = d, REML = TRUE), silent = TRUE)
  } else {
    fit <- try(stats::lm(stats::as.formula(paste("score ~", rhs)),
                         data = d), silent = TRUE)
  }
  if (inherits(fit, "try-error")) return(NULL)
  co <- try(summary(fit)$coefficients, silent = TRUE)
  if (inherits(co, "try-error")) return(NULL)
  rn <- rownames(co)
  hits <- rn[startsWith(rn, axis_term)]
  if (!length(hits)) return(NULL)
  out <- as.data.table(co[hits, , drop = FALSE], keep.rownames = "term")
  ncol_stats <- ncol(out) - 1L
  if (ncol_stats == 5L) {
    setnames(out, c("term", "Estimate", "StdErr", "df", "tval", "pval"))
  } else if (ncol_stats == 4L) {
    setnames(out, c("term", "Estimate", "StdErr", "tval", "pval"))
    out[, df := NA_real_]
  } else {
    return(NULL)
  }
  out[, n_donors := nrow(d)]
  out[, model_class := ifelse(use_random, "lmer", "lm")]
  out[, is_singular := if (use_random) isTRUE(lme4::isSingular(fit)) else NA]
  out
}

fit_axis <- function(axis_term, only_with_value = NULL,
                     min_donors_per_stratum = MIN_DONORS_PER_STRATUM) {
  cat(sprintf("\n[axis] %s\n", axis_term))
  d_in <- lr_long
  if (!is.null(only_with_value)) {
    d_in <- d_in[!is.na(get(only_with_value))]
  }
  if (nrow(d_in) == 0) return(data.table())
  splits <- split(d_in, by = c("ct_pair", "lr_pair"), drop = TRUE)
  cat(sprintf("[axis] %d (ct_pair, lr_pair) groups\n", length(splits)))

  wrapper <- function(d) {
    if (axis_term == "disease_stage_coarse") {
      tab <- table(d$disease_stage_coarse)
      if (any(tab < min_donors_per_stratum)) return(NULL)
    }
    fit_res <- fit_one(d, axis_term)
    if (is.null(fit_res)) return(NULL)
    fit_res[, ct_pair := d$ct_pair[1]]
    fit_res[, lr_pair := d$lr_pair[1]]
    fit_res[, source := d$source[1]]
    fit_res[, target := d$target[1]]
    fit_res[, ligand_complex := d$ligand_complex[1]]
    fit_res[, receptor_complex := d$receptor_complex[1]]
    fit_res
  }

  t0 <- Sys.time()
  if (N_CORES > 1 && .Platform$OS.type != "windows") {
    res_list <- parallel::mclapply(splits, wrapper, mc.cores = N_CORES)
  } else {
    res_list <- lapply(splits, wrapper)
  }
  cat(sprintf("[axis] %s took %.1f min\n",
              axis_term, as.numeric(Sys.time() - t0, units = "mins")))
  out <- rbindlist(res_list[!sapply(res_list, is.null)], fill = TRUE)
  if (nrow(out) == 0) return(out)
  # Critique #7 FIX / Master review M-P0-4: TWO-LEVEL correction.
  # Step 1: BH within (ct_pair, term) -> q_within_ct
  out[, q_within_ct := stats::p.adjust(pval, method = "BH"),
      by = .(ct_pair, term)]
  # Step 2: Bonferroni at the family of ct_pair groups actually FIT on this axis
  #         (not a hardcoded 72; previous v3 over-corrected by up to 3.4x).
  n_families <- length(unique(out$ct_pair))
  cat(sprintf("[axis %s] q_bonferroni_family uses n_families=%d (actual fitted ct_pair count)\n",
              axis_term, n_families))
  out[, n_families := n_families]
  out[, q_bonferroni_family := pmin(q_within_ct * n_families, 1)]
  out
}

# ============================================================================
# (a) 346 equivalent: stage-LR LMM
#   - Coarse axis (PRIMARY cross-axis claim, n=269, NOT scVI-derived)
#   - F-stage documented-only (PRIMARY F-stage, n=58 Andrews, circularity-free)
#   - F-stage augmented (SENSITIVITY, n up to 269)
#   - Continuous: macrophage_pseudotime_mean (now joined correctly)
# ============================================================================
cat("\n========== (346 v3): mixed-effects fits ==========\n")
res_coarse <- fit_axis("disease_stage_coarse")
fwrite(res_coarse, file.path(V3_STAGE_DIR, "stage_lr_lmm_coarse_v3.tsv"),
       sep = "\t")
cat(sprintf("[output] %d rows -> stage_lr_lmm_coarse_v3.tsv\n",
            nrow(res_coarse)))

# Critique #4 FIX: PRIMARY F-stage = documented-only (Andrews n=58)
res_fstage_doc <- data.table()
if ("F_stage_documented" %in% names(lr_long)) {
  n_doc <- sum(!is.na(lr_long$F_stage_documented))
  cat(sprintf("[axis] F_stage_documented: %d non-NA donor x LR rows\n", n_doc))
  if (n_doc > 0) {
    lr_long[, F_stage_doc_numeric := as.numeric(F_stage_documented)]
    res_fstage_doc <- fit_axis("F_stage_doc_numeric",
                               only_with_value = "F_stage_doc_numeric")
  }
}
fwrite(res_fstage_doc,
       file.path(V3_STAGE_DIR, "stage_lr_lmm_fstage_documented_v3.tsv"),
       sep = "\t")
cat(sprintf("[output] %d rows -> stage_lr_lmm_fstage_documented_v3.tsv (PRIMARY)\n",
            nrow(res_fstage_doc)))

# SENSITIVITY: F_stage_augmented_v2 (circularity-susceptible but powered)
res_fstage_aug <- data.table()
if ("F_stage_augmented_v2" %in% names(lr_long)) {
  n_aug <- sum(!is.na(lr_long$F_stage_augmented_v2))
  cat(sprintf("[axis] F_stage_augmented_v2: %d non-NA donor x LR rows\n", n_aug))
  if (n_aug > 0) {
    lr_long[, F_stage_aug_numeric := as.numeric(F_stage_augmented_v2)]
    res_fstage_aug <- fit_axis("F_stage_aug_numeric",
                               only_with_value = "F_stage_aug_numeric")
  }
}
fwrite(res_fstage_aug,
       file.path(V3_STAGE_DIR, "stage_lr_lmm_fstage_augmented_v3.tsv"),
       sep = "\t")
cat(sprintf("[output] %d rows -> stage_lr_lmm_fstage_augmented_v3.tsv (sensitivity)\n",
            nrow(res_fstage_aug)))

# Critique #1 FIX: continuous axis now has data
res_cont <- data.table()
if ("macrophage_pseudotime_mean" %in% names(lr_long) &&
    sum(!is.na(lr_long$macrophage_pseudotime_mean)) >= 30) {
  cat(sprintf("[axis] macrophage_pseudotime_mean: %d non-NA donor x LR rows\n",
              sum(!is.na(lr_long$macrophage_pseudotime_mean))))
  res_cont <- fit_axis("macrophage_pseudotime_mean",
                       only_with_value = "macrophage_pseudotime_mean")
} else {
  cat("[axis] macrophage_pseudotime_mean missing or insufficient\n")
}
fwrite(res_cont,
       file.path(V3_STAGE_DIR, "stage_lr_lmm_continuous_v3.tsv"),
       sep = "\t")
cat(sprintf("[output] %d rows -> stage_lr_lmm_continuous_v3.tsv\n",
            nrow(res_cont)))

# ============================================================================
# (b) 347 equivalent: LOO-cohort replication on coarse axis SH term
# ============================================================================
cat("\n========== (347 v3): LOO dataset replication ==========\n")
sh_term <- "disease_stage_coarseSteatohepatitis"
top <- res_coarse[term == sh_term][order(pval)][, head(.SD, TOP_N)]
cat(sprintf("[loo] selected top %d SH-vs-Healthy hits\n", nrow(top)))

# Use ORIGINAL `dataset` for hold-out (not dataset_pooled) so we can
# stress-test against the dominant Wang dataset specifically.
datasets <- sort(unique(lr_long$dataset))
cat(sprintf("[loo] %d datasets: %s\n", length(datasets),
            paste(datasets, collapse = ", ")))

fit_one_holdout <- function(d, ds_held_out) {
  d <- d[dataset != ds_held_out]
  if (nrow(d) < 20) return(NULL)
  if (length(unique(d$disease_stage_coarse)) < 2) return(NULL)
  use_random <- length(unique(d$dataset_pooled)) >= 3
  fixed <- c("disease_stage_coarse",
             "log10(pmax(n_source_cells,1))",
             "log10(pmax(n_target_cells,1))")
  if (length(unique(d$dataset_pooled)) > 1) fixed <- c(fixed, "dataset_pooled")
  rhs <- paste(fixed, collapse = " + ")
  if (use_random) {
    fit <- try(lmerTest::lmer(stats::as.formula(
      paste("score ~", rhs, "+ (1|dataset_pooled)")),
      data = d, REML = TRUE), silent = TRUE)
  } else {
    fit <- try(stats::lm(stats::as.formula(paste("score ~", rhs)), data = d),
               silent = TRUE)
  }
  if (inherits(fit, "try-error")) return(NULL)
  co <- try(summary(fit)$coefficients, silent = TRUE)
  if (inherits(co, "try-error")) return(NULL)
  rn <- rownames(co)
  hit_row <- rn[startsWith(rn, sh_term)]
  if (!length(hit_row)) return(NULL)
  list(estimate = co[hit_row, 1], pval = co[hit_row, ncol(co)])
}

run_loo_rows <- function() {
  out_rows <- list()
  for (ix in seq_len(nrow(top))) {
    hit <- top[ix]
    d_sub <- lr_long[ct_pair == hit$ct_pair & lr_pair == hit$lr_pair]
    if (nrow(d_sub) == 0) next
    for (ds in datasets) {
      f <- fit_one_holdout(d_sub, ds)
      if (is.null(f)) {
        replicates <- NA
        estimate <- NA_real_
        pval <- NA_real_
      } else {
        replicates <- (sign(f$estimate) == sign(hit$Estimate)) && (f$pval < 0.10)
        estimate <- f$estimate
        pval <- f$pval
      }
      out_rows[[length(out_rows) + 1]] <- data.table(
        ct_pair = hit$ct_pair, lr_pair = hit$lr_pair,
        full_estimate = hit$Estimate,
        full_q_within_ct = hit$q_within_ct,
        full_q_bonferroni_family = hit$q_bonferroni_family,
        held_out = ds, holdout_estimate = estimate,
        holdout_pval = pval, replication_concordant = replicates
      )
    }
    if (ix %% 10 == 0) cat(sprintf("[loo] %d / %d\n", ix, nrow(top)))
  }
  rbindlist(out_rows)
}

if (nrow(top) > 0 && length(datasets) > 1) {
  loo_dt <- run_loo_rows()
  rep_rate <- loo_dt[, .(
    n_holdouts_total = length(datasets),
    n_holdouts_run = sum(!is.na(replication_concordant)),
    n_replicating = sum(replication_concordant, na.rm = TRUE),
    n_failed_fits = sum(is.na(replication_concordant))
  ), by = .(ct_pair, lr_pair)]
  rep_rate[, replication_rate_strict := n_replicating / n_holdouts_total]
  rep_rate[, replication_rate_among_run := ifelse(
    n_holdouts_run > 0, n_replicating / n_holdouts_run, NA_real_
  )]
  rep_rate[, replication_rate := replication_rate_among_run]
  fwrite(rep_rate,
         file.path(V3_STAGE_DIR, "loo_replication_rate_per_lr_v3.tsv"),
         sep = "\t")
  cat(sprintf("[output] %d LR pairs summarized -> loo_replication_rate_per_lr_v3.tsv\n",
              nrow(rep_rate)))
} else {
  rep_rate <- data.table()
  cat("[loo] skipped (insufficient top hits or single cohort)\n")
}

# ============================================================================
# (c) 348 equivalent: bulk concordance
# ============================================================================
cat("\n========== (348 v3): bulk concordance ==========\n")

load_bulk <- function(path) {
  if (!file.exists(path)) return(NULL)
  bulk <- fread(path)
  gene_col <- intersect(
    c("symbol", "gene_symbol", "hgnc_symbol", "gene", "Gene", "ID"),
    names(bulk)
  )[1]
  lfc_col <- intersect(c("logFC", "log2FoldChange", "lfc"), names(bulk))[1]
  padj_col <- intersect(c("padj", "adj.P.Val", "BH", "FDR"), names(bulk))[1]
  if (is.null(gene_col) || is.null(lfc_col) || is.null(padj_col)) return(NULL)
  out <- bulk[, .(gene = get(gene_col),
                  bulk_lfc = get(lfc_col),
                  bulk_padj = get(padj_col))]
  out <- out[!is.na(gene) & gene != ""]
  out <- out[order(bulk_padj)][, .SD[1], by = gene]
  out
}

bulk_all <- load_bulk(BULK_ALL)
bulk_st <- load_bulk(BULK_STEATOSIS)
bulk_sh <- load_bulk(BULK_SH)
bulk_cirr <- load_bulk(BULK_CIRRHOSIS)

score_concord <- function(lmm, anchor) {
  if (nrow(lmm) == 0 || is.null(anchor)) return(data.table())
  lmm <- copy(lmm)
  lmm[, ligand := sub("_.*", "", ligand_complex)]
  lmm[, receptor := sub("_.*", "", receptor_complex)]
  m_lig <- merge(lmm, anchor, by.x = "ligand", by.y = "gene", all.x = TRUE)
  setnames(m_lig, c("bulk_lfc", "bulk_padj"),
                  c("lig_bulk_lfc", "lig_bulk_padj"))
  m <- merge(m_lig, anchor, by.x = "receptor", by.y = "gene", all.x = TRUE)
  setnames(m, c("bulk_lfc", "bulk_padj"),
              c("rec_bulk_lfc", "rec_bulk_padj"))
  m[, sc_sign := sign(Estimate)]
  m[, lig_concordant := (sign(lig_bulk_lfc) == sc_sign) & (lig_bulk_padj < 0.1)]
  m[, rec_concordant := (sign(rec_bulk_lfc) == sc_sign) & (rec_bulk_padj < 0.1)]
  # Master review M-P0-9: R's three-valued logic treats FALSE & NA as FALSE,
  # which silently counted untested-receptor rows (153 per axis) in the
  # concordance denominator. Require BOTH sides to be testable (non-NA).
  m[, both_concordant := lig_concordant & rec_concordant &
                         !is.na(lig_bulk_padj) & !is.na(rec_bulk_padj)]
  m
}

score_save_axis <- function(term_label, anchor, axis_name) {
  slice <- res_coarse[term == term_label]
  if (nrow(slice) == 0) return(data.table())
  if (is.null(anchor)) {
    anchor <- bulk_all
    axis_tag <- paste0(axis_name, "_fallback_all_MASLD")
  } else {
    axis_tag <- axis_name
  }
  conc <- score_concord(slice, anchor)
  if (nrow(conc) > 0) conc[, axis := axis_tag]
  conc
}

steat_conc <- score_save_axis(
  "disease_stage_coarseSteatosis", bulk_st, "Steatosis_vs_Healthy"
)
sh_conc <- score_save_axis(
  "disease_stage_coarseSteatohepatitis", bulk_sh, "Steatohepatitis_vs_Healthy"
)
cirr_conc <- score_save_axis(
  "disease_stage_coarseCirrhosis", bulk_cirr, "Cirrhosis_vs_Healthy"
)

cont_conc <- data.table()
if (nrow(res_cont) > 0 && !is.null(bulk_all)) {
  cont_conc <- score_concord(res_cont, bulk_all)
  if (nrow(cont_conc) > 0) {
    cont_conc[, axis := "macrophage_pseudotime_all_MASLD"]
  }
}

fstage_doc_conc <- data.table()
if (nrow(res_fstage_doc) > 0 && !is.null(bulk_all)) {
  fstage_doc_slice <- res_fstage_doc[term == "F_stage_doc_numeric"]
  fstage_doc_conc <- score_concord(fstage_doc_slice, bulk_all)
  if (nrow(fstage_doc_conc) > 0) {
    fstage_doc_conc[, axis := "F_stage_documented_all_MASLD"]
  }
}

fstage_aug_conc <- data.table()
if (nrow(res_fstage_aug) > 0 && !is.null(bulk_all)) {
  fstage_aug_slice <- res_fstage_aug[term == "F_stage_aug_numeric"]
  fstage_aug_conc <- score_concord(fstage_aug_slice, bulk_all)
  if (nrow(fstage_aug_conc) > 0) {
    fstage_aug_conc[, axis := "F_stage_augmented_all_MASLD"]
  }
}

all_conc <- rbindlist(
  list(steat_conc, sh_conc, cirr_conc, cont_conc,
       fstage_doc_conc, fstage_aug_conc),
  fill = TRUE
)
fwrite(all_conc, file.path(V3_STAGE_DIR, "lr_bulk_concordance_v3.tsv"),
       sep = "\t")
cat(sprintf("[output] %d total rows -> lr_bulk_concordance_v3.tsv\n",
            nrow(all_conc)))

if (nrow(all_conc) > 0) {
  hdr <- all_conc[!is.na(both_concordant),
                  .(n_testable = .N,
                    n_both = sum(both_concordant, na.rm = TRUE),
                    pct_both = 100 * mean(both_concordant, na.rm = TRUE),
                    pct_lig = 100 * mean(lig_concordant, na.rm = TRUE),
                    pct_rec = 100 * mean(rec_concordant, na.rm = TRUE)),
                  by = axis]
  cat("\n[headline] both_concordant per axis:\n")
  print(hdr)
}

# ============================================================================
# (d) sanity counts: log how many LR pairs pass each correction level
# ============================================================================
cat("\n========== sanity: q_within_ct vs q_bonferroni_family ==========\n")
for (axis_name in c("coarse", "fstage_documented", "fstage_augmented",
                    "continuous")) {
  obj_name <- switch(axis_name,
                     coarse = "res_coarse",
                     fstage_documented = "res_fstage_doc",
                     fstage_augmented = "res_fstage_aug",
                     continuous = "res_cont")
  res <- get(obj_name)
  if (nrow(res) == 0) {
    cat(sprintf("  %-22s: 0 rows\n", axis_name))
    next
  }
  n_w005 <- sum(res$q_within_ct < 0.05, na.rm = TRUE)
  n_f005 <- sum(res$q_bonferroni_family < 0.05, na.rm = TRUE)
  cat(sprintf("  %-22s: %d rows | q_within_ct<0.05: %d | q_bonferroni_family<0.05: %d\n",
              axis_name, nrow(res), n_w005, n_f005))
}

cat("\n[done] 07b_chain_v3_hardened finished\n")
