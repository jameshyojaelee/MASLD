#!/usr/bin/env Rscript
# ============================================================================
# 347_stage_ccc_loo_dataset.R
#
# Leave-one-dataset-out replication: hold out each of the 7 scRNA datasets,
# refit the disease_stage_coarse mixed model on the remaining 6, and score
# the top-100 stage-progressive LR pairs in the held-out dataset.
#
# Outputs:
#   loo_dataset_replication.tsv
#   Columns: ct_pair, lr_pair, term, full_estimate, full_padj_within_ct,
#            held_out, holdout_estimate, holdout_pval, replication_concordant
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(lme4)
  library(lmerTest)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
OUT_DIR <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory")
PER_DONOR_DIR <- file.path(OUT_DIR, "per_donor_lr")
META_EXT <- file.path(OUT_DIR, "donor_metadata_extended.tsv")
COARSE_TSV <- file.path(OUT_DIR, "stage_lr_lmm_coarse.tsv")
FSTAGE_TSV <- file.path(OUT_DIR, "stage_lr_lmm_fstage.tsv")

TOP_N <- 100      # top stage-progressive LR pairs to replicate
MIN_DONORS_PER_LR <- 20  # relaxed for held-out subset

stopifnot(file.exists(COARSE_TSV))
stopifnot(file.exists(META_EXT))

meta <- fread(META_EXT)
# Backwards-compat: protocol-contamination flag and clean F_stage column may
# be missing if Phase 1 atlas refresh hasn't completed yet.
if (!"exclude_stage_analysis" %in% names(meta))
  meta[, exclude_stage_analysis := FALSE]
if (!"F_stage_augmented_clean" %in% names(meta))
  meta[, F_stage_augmented_clean := if ("F_stage_augmented" %in% names(meta))
       F_stage_augmented else NA_real_]
res_full <- fread(COARSE_TSV)
res_fstage <- if (file.exists(FSTAGE_TSV)) fread(FSTAGE_TSV) else NULL

# Pick the top-N stage-progressive LR pairs from the full fit ---------------
# Definition: most-significant Steatohepatitis-vs-Healthy contrast (the SH
# stratum is the largest non-Healthy bin and the most reliable stage signal).
sh_term <- "disease_stage_coarseSteatohepatitis"
top <- res_full[term == sh_term][order(pval)][, head(.SD, TOP_N)]
cat(sprintf("[top] selected top %d SH-vs-Healthy hits from full fit\n", nrow(top)))

# Load merged TSV (produced by Script 345b) ----------------------------------
MERGED_TSV <- file.path(OUT_DIR, "all_donor_lr_scores.tsv.gz")
if (!file.exists(MERGED_TSV)) {
  stop("Merged LR scores TSV not found - run Script 345b first")
}
lr_long <- fread(MERGED_TSV)
meta_cols <- meta[, .(sample, dataset, disease_stage_coarse,
                      exclude_stage_analysis,
                      age = if ("age" %in% names(meta)) age else NA_real_,
                      sex_numeric = if ("sex_numeric" %in% names(meta)) sex_numeric else NA_integer_,
                      F_stage_inferred = if ("F_stage_inferred" %in% names(meta)) F_stage_inferred else NA_real_,
                      F_stage_augmented_clean = F_stage_augmented_clean)]
lr_long <- merge(lr_long, meta_cols, by = "sample", all.x = TRUE)
# Drop protocol-contaminated donors (GSE136103, Liver_Atlas) before any LMM.
n_before_excl <- nrow(lr_long)
n_donors_before_excl <- uniqueN(lr_long$sample)
lr_long <- lr_long[is.na(exclude_stage_analysis) | exclude_stage_analysis == FALSE]
cat(sprintf("[filter] excluded %d rows (%d donors) flagged exclude_stage_analysis; %d rows (%d donors) remain\n",
            n_before_excl - nrow(lr_long),
            n_donors_before_excl - uniqueN(lr_long$sample),
            nrow(lr_long), uniqueN(lr_long$sample)))
lr_long[, score := -log10(pmax(magnitude_rank, 1e-4))]
lr_long[, disease_stage_coarse := factor(disease_stage_coarse,
        levels = c("Healthy", "Steatosis", "Steatohepatitis", "Cirrhosis"))]
# Prefer the contamination-clean augmented F-stage; fall back to inferred.
if ("F_stage_augmented_clean" %in% names(lr_long) &&
    sum(!is.na(lr_long$F_stage_augmented_clean)) > 0) {
  lr_long[, F_stage_numeric := as.numeric(F_stage_augmented_clean)]
} else {
  lr_long[, F_stage_numeric := as.numeric(F_stage_inferred)]
}
lr_long_all <- copy(lr_long)
lr_long <- lr_long[!is.na(disease_stage_coarse)]
lr_long[, ct_pair := paste(source, target, sep = "->")]
lr_long[, lr_pair := paste(ligand_complex, receptor_complex, sep = "__")]
lr_long_all[, ct_pair := paste(source, target, sep = "->")]
lr_long_all[, lr_pair := paste(ligand_complex, receptor_complex, sep = "__")]

datasets <- sort(unique(lr_long$dataset))
cat(sprintf("[loo] %d datasets (coarse axis): %s\n", length(datasets),
            paste(datasets, collapse = ", ")))

# Generic holdout fit. axis = "coarse" (disease_stage_coarse) or "fstage"
# (F_stage_numeric continuous). For fstage the term name returned is
# `F_stage_numeric` (a single slope term).
fit_one_holdout <- function(d, ds_held_out, axis = "coarse") {
  d <- d[dataset != ds_held_out]
  if (nrow(d) < MIN_DONORS_PER_LR) return(NULL)
  if (axis == "coarse") {
    if (length(unique(d$disease_stage_coarse)) < 2) return(NULL)
    axis_term <- "disease_stage_coarse"
    hit_prefix <- "disease_stage_coarseSteatohepatitis"
  } else if (axis == "fstage") {
    d <- d[!is.na(F_stage_numeric)]
    if (nrow(d) < MIN_DONORS_PER_LR) return(NULL)
    if (length(unique(d$F_stage_numeric)) < 2) return(NULL)
    axis_term <- "F_stage_numeric"
    hit_prefix <- "F_stage_numeric"
  } else stop("unknown axis")
  use_random <- length(unique(d$dataset)) >= 3
  fixed <- c(axis_term,
             "log10(pmax(n_source_cells,1))",
             "log10(pmax(n_target_cells,1))")
  # Master review M-P0-3: dataset must NOT be both fixed and random.
  if (!use_random && length(unique(d$dataset)) > 1) fixed <- c(fixed, "dataset")
  rhs <- paste(fixed, collapse = " + ")
  if (use_random) {
    fit <- try(lmerTest::lmer(stats::as.formula(paste("score ~", rhs, "+ (1|dataset)")),
                              data = d, REML = TRUE), silent = TRUE)
  } else {
    fit <- try(stats::lm(stats::as.formula(paste("score ~", rhs)), data = d), silent = TRUE)
  }
  if (inherits(fit, "try-error")) return(NULL)
  co <- try(summary(fit)$coefficients, silent = TRUE)
  if (inherits(co, "try-error")) return(NULL)
  rn <- rownames(co)
  hit_row <- rn[startsWith(rn, hit_prefix)]
  if (!length(hit_row)) return(NULL)
  list(estimate = co[hit_row, 1], pval = co[hit_row, ncol(co)])
}

run_loo <- function(top_hits, lr_src, term_label, axis, datasets_axis) {
  out_rows <- list()
  for (ix in seq_len(nrow(top_hits))) {
    hit <- top_hits[ix]
    d_sub <- lr_src[ct_pair == hit$ct_pair & lr_pair == hit$lr_pair]
    if (nrow(d_sub) == 0) next
    for (ds in datasets_axis) {
      f <- fit_one_holdout(d_sub, ds, axis = axis)
      if (is.null(f)) {
        replicates <- NA
        estimate   <- NA_real_
        pval       <- NA_real_
      } else {
        replicates <- (sign(f$estimate) == sign(hit$Estimate)) && (f$pval < 0.10)
        estimate   <- f$estimate
        pval       <- f$pval
      }
      out_rows[[length(out_rows) + 1]] <- data.table(
        ct_pair = hit$ct_pair,
        lr_pair = hit$lr_pair,
        term    = term_label,
        full_estimate = hit$Estimate,
        full_padj_within_ct = hit$padj_within_ct,
        held_out = ds,
        holdout_estimate = estimate,
        holdout_pval = pval,
        replication_concordant = replicates
      )
    }
    if (ix %% 10 == 0) cat(sprintf("[loo:%s] %d / %d hits processed\n",
                                   axis, ix, nrow(top_hits)))
  }
  rbindlist(out_rows)
}

# --- Coarse SH-vs-Healthy LOO (legacy) -------------------------------------
out_dt <- run_loo(top, lr_long, sh_term, axis = "coarse", datasets_axis = datasets)
fwrite(out_dt, file.path(OUT_DIR, "loo_dataset_replication.tsv"), sep = "\t")

# Aggregate replication rate per hit ---------------------------------------
summarize_rep <- function(out_dt, n_holdouts_total, datasets_used, axis_label) {
  if (nrow(out_dt) == 0) return(data.table())
  rep_rate <- out_dt[, .(n_holdouts_total = n_holdouts_total,
                         n_holdouts_run   = sum(!is.na(replication_concordant)),
                         n_replicating    = sum(replication_concordant, na.rm = TRUE),
                         n_failed_fits    = sum(is.na(replication_concordant))),
                     by = .(ct_pair, lr_pair)]
  rep_rate[, replication_rate_strict    := n_replicating / n_holdouts_total]
  rep_rate[, replication_rate_among_run := ifelse(n_holdouts_run > 0,
                                                  n_replicating / n_holdouts_run,
                                                  NA_real_)]
  rep_rate[, replication_rate := replication_rate_among_run]
  rep_rate[, axis := axis_label]

  n_lr <- nrow(rep_rate)
  cat(sprintf("\n[summary:%s] LOO scope: %d datasets (%s)\n",
              axis_label, n_holdouts_total, paste(datasets_used, collapse = ", ")))
  cat(sprintf("[summary:%s] failed holdout fits across all LR x dataset cells: %d / %d\n",
              axis_label, sum(rep_rate$n_failed_fits), n_holdouts_total * n_lr))
  cat(sprintf("[summary:%s] median replication_rate_strict     = %.3f (denom = %d)\n",
              axis_label, median(rep_rate$replication_rate_strict), n_holdouts_total))
  cat(sprintf("[summary:%s] median replication_rate_among_run  = %.3f (run-only denom)\n",
              axis_label, median(rep_rate$replication_rate_among_run, na.rm = TRUE)))
  for (k in seq(n_holdouts_total, 1)) {
    cat(sprintf("[summary:%s]   %d / %d strict: %d LR pairs (%.1f%%)\n",
                axis_label, k, n_holdouts_total,
                sum(rep_rate$n_replicating >= k),
                100 * mean(rep_rate$n_replicating >= k)))
  }
  rep_rate
}

n_holdouts_total <- length(datasets)
rep_rate <- summarize_rep(out_dt, n_holdouts_total, datasets, "coarse_SH_vs_Healthy")
fwrite(rep_rate, file.path(OUT_DIR, "loo_replication_rate_per_lr.tsv"), sep = "\t")
cat(sprintf("\n[output] %d holdout rows -> loo_dataset_replication.tsv\n",
            nrow(out_dt)))
cat(sprintf("[output] %d LR pairs summarized -> loo_replication_rate_per_lr.tsv\n",
            nrow(rep_rate)))
strict_thresh <- 4 / n_holdouts_total
cat(sprintf("[headline:coarse] %d / %d LR pairs replicate in >= %d of %d holdouts\n",
            sum(rep_rate$replication_rate_strict >= strict_thresh),
            nrow(rep_rate), ceiling(strict_thresh * n_holdouts_total),
            n_holdouts_total))

# --- F-stage continuous axis LOO -------------------------------------------
rep_rate_fstage <- data.table()
if (!is.null(res_fstage) && nrow(res_fstage) > 0) {
  fstage_term <- "F_stage_numeric"
  lr_long_fs <- lr_long_all[!is.na(F_stage_numeric)]
  lr_long_fs[, ct_pair := paste(source, target, sep = "->")]
  lr_long_fs[, lr_pair := paste(ligand_complex, receptor_complex, sep = "__")]
  datasets_fs <- sort(unique(lr_long_fs$dataset))
  cat(sprintf("\n[loo:fstage] %d datasets (F-stage gated): %s\n",
              length(datasets_fs), paste(datasets_fs, collapse = ", ")))

  # Top-N by |Estimate| * -log10(p)
  res_fstage[, rank_score := abs(Estimate) * (-log10(pmax(pval, 1e-300)))]
  top_fs <- res_fstage[term == fstage_term][order(-rank_score)][, head(.SD, TOP_N)]
  cat(sprintf("[top:fstage] selected top %d F-stage hits from full fit\n",
              nrow(top_fs)))

  out_dt_fs <- run_loo(top_fs, lr_long_fs, fstage_term,
                       axis = "fstage", datasets_axis = datasets_fs)
  fwrite(out_dt_fs, file.path(OUT_DIR, "loo_dataset_replication_fstage.tsv"),
         sep = "\t")
  n_holdouts_fs <- length(datasets_fs)
  rep_rate_fstage <- summarize_rep(out_dt_fs, n_holdouts_fs, datasets_fs,
                                   "F_stage_continuous")
  fwrite(rep_rate_fstage,
         file.path(OUT_DIR, "loo_replication_rate_per_lr_fstage.tsv"), sep = "\t")
  cat(sprintf("\n[output] %d holdout rows -> loo_dataset_replication_fstage.tsv\n",
              nrow(out_dt_fs)))
  cat(sprintf("[output] %d LR pairs summarized -> loo_replication_rate_per_lr_fstage.tsv\n",
              nrow(rep_rate_fstage)))
  strict_thresh_fs <- if (n_holdouts_fs > 0) (n_holdouts_fs - 1) / n_holdouts_fs else NA
  if (n_holdouts_fs >= 6) strict_thresh_fs <- 5 / n_holdouts_fs
  cat(sprintf("[headline:fstage] %d / %d LR pairs replicate in all %d holdouts (perfect)\n",
              sum(rep_rate_fstage$n_replicating == n_holdouts_fs),
              nrow(rep_rate_fstage), n_holdouts_fs))
  if (n_holdouts_fs >= 6) {
    cat(sprintf("[headline:fstage] %d / %d LR pairs replicate in >= 5/%d holdouts\n",
                sum(rep_rate_fstage$n_replicating >= 5),
                nrow(rep_rate_fstage), n_holdouts_fs))
  }
} else {
  cat("\n[loo:fstage] stage_lr_lmm_fstage.tsv not found; F-stage LOO skipped.\n")
}

# --- Combined output -------------------------------------------------------
rep_rate$axis <- "coarse_SH_vs_Healthy"
combined <- rbindlist(list(rep_rate, rep_rate_fstage), fill = TRUE)
fwrite(combined,
       file.path(OUT_DIR, "loo_replication_rate_per_lr_combined.tsv"), sep = "\t")
cat(sprintf("\n[combined] %d total LR rows -> loo_replication_rate_per_lr_combined.tsv\n",
            nrow(combined)))

# --- Jaccard top-50 overlap between axes -----------------------------------
if (nrow(rep_rate_fstage) > 0) {
  top50_coarse <- rep_rate[order(-replication_rate_strict, -n_replicating)
                          ][, head(.SD, 50)
                          ][, paste(ct_pair, lr_pair, sep = "::")]
  top50_fs <- rep_rate_fstage[order(-replication_rate_strict, -n_replicating)
                             ][, head(.SD, 50)
                             ][, paste(ct_pair, lr_pair, sep = "::")]
  inter <- length(intersect(top50_coarse, top50_fs))
  uni   <- length(union(top50_coarse, top50_fs))
  jacc  <- if (uni > 0) inter / uni else NA
  cat(sprintf("\n[compare] top-50 Jaccard (coarse_SH vs F-stage) = %.3f (inter=%d, union=%d)\n",
              jacc, inter, uni))
}
