#!/usr/bin/env Rscript
# ============================================================================
# 347b_gse136103_holdout_sensitivity.R
#
# Dissociation-artifact sensitivity test for stage-stratified CCC.
#
# GSE136103 (Ramachandran 2019) used enzymatic dissociation, which depletes
# fragile hepatocytes and over-represents macrophages / lymphoid cells. A
# disproportionate share of our Cirrhosis donors come from GSE136103 (9 of 28).
# This script asks: which Cirrhosis-emergent LR pairs survive when GSE136103
# is held out, vs which are dissociation-driven?
#
# Strategy:
#   1. Take the top-100 Cirrhosis-vs-Healthy hits (term ==
#      "disease_stage_coarseCirrhosis", Estimate > 0, sorted by pval) from
#      the full mega-fit (stage_lr_lmm_coarse.tsv).
#   2. Refit the SAME lmer (mirroring Script 346/347) on the full atlas AND
#      on the GSE136103-excluded atlas, per LR pair.
#   3. Classify each pair as survives_dissoc / weakened / lost / reversed /
#      untestable.
#
# Reuses the per-donor merged LR scores; does NOT re-run LIANA.
#
# Output: cirrhosis_dissoc_sensitivity.tsv
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(lme4)
  library(lmerTest)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
OUT_DIR <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory")
META_EXT <- file.path(OUT_DIR, "donor_metadata_extended.tsv")
COARSE_TSV <- file.path(OUT_DIR, "stage_lr_lmm_coarse.tsv")
MERGED_TSV <- file.path(OUT_DIR, "all_donor_lr_scores.tsv.gz")

TOP_N <- 100
HELD_OUT_DATASET <- "GSE136103"
TARGET_TERM <- "disease_stage_coarseCirrhosis"
MIN_CIRRHOSIS_FOR_TEST <- 3
PVAL_CUTOFF <- 0.10

stopifnot(file.exists(COARSE_TSV), file.exists(META_EXT), file.exists(MERGED_TSV))

# ---- Load full-fit Cirrhosis-vs-Healthy hits ---------------------------------
res_full <- fread(COARSE_TSV)
top <- res_full[term == TARGET_TERM & Estimate > 0][order(pval)][, head(.SD, TOP_N)]
cat(sprintf("[top] selected top %d Cirrhosis-emergent (positive-effect) hits\n",
            nrow(top)))

# ---- Load metadata + merged LR scores (mirror 346/347 setup) -----------------
meta <- fread(META_EXT)
lr_long <- fread(MERGED_TSV)
cat(sprintf("[input] %d donor x LR-pair rows loaded\n", nrow(lr_long)))

covar_cols <- c("sample", "dataset", "disease_stage_coarse",
                "age", "sex_numeric")
covar_cols <- intersect(covar_cols, names(meta))
lr_long <- merge(lr_long, meta[, ..covar_cols], by = "sample", all.x = TRUE)
lr_long[, score := -log10(pmax(magnitude_rank, 1e-4))]
lr_long[, disease_stage_coarse := factor(disease_stage_coarse,
        levels = c("Healthy", "Steatosis", "Steatohepatitis", "Cirrhosis"))]
lr_long <- lr_long[!is.na(disease_stage_coarse)]
lr_long[, ct_pair := paste(source, target, sep = "->")]
lr_long[, lr_pair := paste(ligand_complex, receptor_complex, sep = "__")]

# Sanity-check holdout sample sizes
unique_donors_full <- unique(meta[, .(sample, dataset, disease_stage_coarse)])
unique_donors_held_in <- unique_donors_full[dataset != HELD_OUT_DATASET]
n_cirrhosis_full <- sum(unique_donors_full$disease_stage_coarse == "Cirrhosis",
                        na.rm = TRUE)
n_cirrhosis_held_in_atlas <- sum(
  unique_donors_held_in$disease_stage_coarse == "Cirrhosis", na.rm = TRUE)
cat(sprintf("[atlas] full=%d donors, %d Cirrhosis; held-in=%d donors, %d Cirrhosis\n",
            nrow(unique_donors_full), n_cirrhosis_full,
            nrow(unique_donors_held_in), n_cirrhosis_held_in_atlas))

# ---- Fit helper (mirrors fit_one_holdout in Script 347) ----------------------
fit_cirrhosis <- function(d) {
  if (nrow(d) < 20) return(NULL)
  tab <- table(d$disease_stage_coarse)
  # Require at least one Healthy AND Cirrhosis donor in the subset
  if (!("Healthy" %in% names(tab)) || tab[["Healthy"]] < 3) return(NULL)
  if (!("Cirrhosis" %in% names(tab)) || tab[["Cirrhosis"]] < 3) return(NULL)
  n_ds <- length(unique(d$dataset))
  use_random <- n_ds >= 3
  fixed <- c("disease_stage_coarse",
             "log10(pmax(n_source_cells,1))",
             "log10(pmax(n_target_cells,1))")
  # Master review M-P0-3: dataset must NOT be both fixed and random.
  if (!use_random && n_ds > 1) fixed <- c(fixed, "dataset")
  rhs <- paste(fixed, collapse = " + ")
  if (use_random) {
    fit <- try(lmerTest::lmer(stats::as.formula(paste("score ~", rhs, "+ (1|dataset)")),
                              data = d, REML = TRUE), silent = TRUE)
  } else {
    fit <- try(stats::lm(stats::as.formula(paste("score ~", rhs)),
                         data = d), silent = TRUE)
  }
  if (inherits(fit, "try-error")) return(NULL)
  co <- try(summary(fit)$coefficients, silent = TRUE)
  if (inherits(co, "try-error")) return(NULL)
  rn <- rownames(co)
  hit_row <- rn[startsWith(rn, TARGET_TERM)]
  if (!length(hit_row)) return(NULL)
  list(estimate = co[hit_row, 1],
       pval     = co[hit_row, ncol(co)],
       n_total  = nrow(d),
       n_donors_unique = length(unique(d$sample)),
       n_cirrhosis = as.integer(tab[["Cirrhosis"]]))
}

# ---- Iterate over top-N hits -------------------------------------------------
out_rows <- vector("list", nrow(top))
for (ix in seq_len(nrow(top))) {
  hit <- top[ix]
  d_sub <- lr_long[ct_pair == hit$ct_pair & lr_pair == hit$lr_pair]
  if (nrow(d_sub) == 0) next

  # Full atlas
  f_full <- fit_cirrhosis(d_sub)
  # Holdout: drop GSE136103
  d_hold <- d_sub[dataset != HELD_OUT_DATASET]
  f_hold <- fit_cirrhosis(d_hold)

  full_estimate <- if (!is.null(f_full)) f_full$estimate else NA_real_
  full_pval     <- if (!is.null(f_full)) f_full$pval     else NA_real_
  hold_estimate <- if (!is.null(f_hold)) f_hold$estimate else NA_real_
  hold_pval     <- if (!is.null(f_hold)) f_hold$pval     else NA_real_
  n_donors_held_in <- if (!is.null(f_hold)) f_hold$n_donors_unique else
                      length(unique(d_hold$sample))
  n_cirrhosis_held_in <- if (!is.null(f_hold)) f_hold$n_cirrhosis else
    sum(unique(d_hold[, .(sample, disease_stage_coarse)])$disease_stage_coarse ==
        "Cirrhosis", na.rm = TRUE)

  # Classification
  if (is.na(hold_estimate) || is.na(hold_pval) ||
      n_cirrhosis_held_in < MIN_CIRRHOSIS_FOR_TEST) {
    cls <- "untestable"
  } else if (sign(hold_estimate) != sign(full_estimate)) {
    cls <- "reversed"
  } else if (!is.na(full_pval) && full_pval < PVAL_CUTOFF &&
             hold_pval < PVAL_CUTOFF) {
    cls <- "survives_dissoc"
  } else if (!is.na(full_pval) && full_pval < PVAL_CUTOFF &&
             hold_pval >= PVAL_CUTOFF) {
    cls <- "weakened"
  } else {
    cls <- "lost"
  }

  out_rows[[ix]] <- data.table(
    ct_pair = hit$ct_pair,
    lr_pair = hit$lr_pair,
    source  = hit$source,
    target  = hit$target,
    ligand_complex   = hit$ligand_complex,
    receptor_complex = hit$receptor_complex,
    term = TARGET_TERM,
    full_estimate = full_estimate,
    full_pval     = full_pval,
    holdout_estimate = hold_estimate,
    holdout_pval     = hold_pval,
    n_donors_held_in = n_donors_held_in,
    n_cirrhosis_held_in = n_cirrhosis_held_in,
    classification = cls
  )

  if (ix %% 20 == 0) cat(sprintf("[loop] %d / %d hits processed\n", ix, nrow(top)))
}

out_dt <- rbindlist(out_rows, fill = TRUE)
out_path <- file.path(OUT_DIR, "cirrhosis_dissoc_sensitivity.tsv")
fwrite(out_dt, out_path, sep = "\t")
cat(sprintf("\n[output] %d rows -> %s\n", nrow(out_dt), out_path))

# ---- Summary -----------------------------------------------------------------
cat("\n[summary] classification table:\n")
print(table(out_dt$classification, useNA = "ifany"))

cat("\n[summary] sample sizes (median across top-100):\n")
cat(sprintf("  median n_donors_held_in    = %.0f\n",
            stats::median(out_dt$n_donors_held_in, na.rm = TRUE)))
cat(sprintf("  median n_cirrhosis_held_in = %.0f\n",
            stats::median(out_dt$n_cirrhosis_held_in, na.rm = TRUE)))

cat("\n[summary] non-surviving hits (weakened / lost / reversed):\n")
non_surv <- out_dt[classification %in% c("weakened", "lost", "reversed")]
print(non_surv[, .(ct_pair, lr_pair, classification,
                   full_estimate, full_pval,
                   holdout_estimate, holdout_pval,
                   n_cirrhosis_held_in)])

cat("\n[done] dissociation sensitivity test complete\n")
