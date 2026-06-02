#!/usr/bin/env Rscript
# ============================================================================
# 347c_gse136103_holdout_sensitivity_sh.R
#
# Sister script to 347b. Extends the GSE136103-holdout dissociation-artifact
# sensitivity test from the Cirrhosis-vs-Healthy contrast (347b) to the
# Steatohepatitis-vs-Healthy contrast, which is the contrast actually used
# by Script 349b/349c for the paralog-dedup top-30 figure.
#
# Identical methodology to 347b (top-100 by pval among positive-effect
# coarse-stage LMM hits; refit on full atlas vs GSE136103-excluded atlas;
# survives / weakened / lost / reversed / untestable classification), but
# parameterised on:
#   TARGET_TERM      = "disease_stage_coarseSteatohepatitis"
#   MIN_STAGE_FOR_TEST = require >=3 SH donors in the held-in subset
#
# Then merges 347b (Cirrhosis) + 347c (SH) outputs into one combined
# `all_dissoc_sensitivity.tsv` with a `contrast` tag, so downstream
# annotation (349c) can cover both the SH-emergent dedup top-30 figure and
# the cirrhosis-emergent sensitivity figures.
#
# Outputs:
#   results_gpu_v2/ccc/stage_trajectory/sh_dissoc_sensitivity.tsv
#   results_gpu_v2/ccc/stage_trajectory/all_dissoc_sensitivity.tsv
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(lme4)
  library(lmerTest)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
OUT_DIR  <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory")
META_EXT <- file.path(OUT_DIR, "donor_metadata_extended.tsv")
COARSE_TSV <- file.path(OUT_DIR, "stage_lr_lmm_coarse.tsv")
MERGED_TSV <- file.path(OUT_DIR, "all_donor_lr_scores.tsv.gz")
CIRR_TSV   <- file.path(OUT_DIR, "cirrhosis_dissoc_sensitivity.tsv")

TOP_N <- 100
HELD_OUT_DATASET <- "GSE136103"
TARGET_TERM <- "disease_stage_coarseSteatohepatitis"
TARGET_LEVEL <- "Steatohepatitis"
MIN_STAGE_FOR_TEST <- 3
PVAL_CUTOFF <- 0.10

stopifnot(file.exists(COARSE_TSV), file.exists(META_EXT), file.exists(MERGED_TSV))

# ---- Load full-fit SH-vs-Healthy hits ---------------------------------------
res_full <- fread(COARSE_TSV)
top <- res_full[term == TARGET_TERM & Estimate > 0][order(pval)][, head(.SD, TOP_N)]
cat(sprintf("[top] selected top %d SH-emergent (positive-effect) hits\n",
            nrow(top)))

# ---- Load metadata + merged LR scores (mirror 346/347/347b) ------------------
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
n_sh_full <- sum(unique_donors_full$disease_stage_coarse == TARGET_LEVEL,
                 na.rm = TRUE)
n_sh_held_in <- sum(
  unique_donors_held_in$disease_stage_coarse == TARGET_LEVEL, na.rm = TRUE)
cat(sprintf("[atlas] full=%d donors, %d %s; held-in=%d donors, %d %s\n",
            nrow(unique_donors_full), n_sh_full, TARGET_LEVEL,
            nrow(unique_donors_held_in), n_sh_held_in, TARGET_LEVEL))

# ---- Fit helper -------------------------------------------------------------
fit_target <- function(d) {
  if (nrow(d) < 20) return(NULL)
  tab <- table(d$disease_stage_coarse)
  if (!("Healthy" %in% names(tab)) || tab[["Healthy"]] < 3) return(NULL)
  if (!(TARGET_LEVEL %in% names(tab)) || tab[[TARGET_LEVEL]] < 3) return(NULL)
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
       n_target = as.integer(tab[[TARGET_LEVEL]]))
}

# ---- Iterate over top-N hits ------------------------------------------------
out_rows <- vector("list", nrow(top))
for (ix in seq_len(nrow(top))) {
  hit <- top[ix]
  d_sub <- lr_long[ct_pair == hit$ct_pair & lr_pair == hit$lr_pair]
  if (nrow(d_sub) == 0) next

  f_full <- fit_target(d_sub)
  d_hold <- d_sub[dataset != HELD_OUT_DATASET]
  f_hold <- fit_target(d_hold)

  full_estimate <- if (!is.null(f_full)) f_full$estimate else NA_real_
  full_pval     <- if (!is.null(f_full)) f_full$pval     else NA_real_
  hold_estimate <- if (!is.null(f_hold)) f_hold$estimate else NA_real_
  hold_pval     <- if (!is.null(f_hold)) f_hold$pval     else NA_real_
  n_donors_held_in <- if (!is.null(f_hold)) f_hold$n_donors_unique else
                      length(unique(d_hold$sample))
  n_target_held_in <- if (!is.null(f_hold)) f_hold$n_target else
    sum(unique(d_hold[, .(sample, disease_stage_coarse)])$disease_stage_coarse ==
        TARGET_LEVEL, na.rm = TRUE)

  if (is.na(hold_estimate) || is.na(hold_pval) ||
      n_target_held_in < MIN_STAGE_FOR_TEST) {
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
    n_target_held_in = n_target_held_in,
    classification = cls
  )

  if (ix %% 20 == 0) cat(sprintf("[loop] %d / %d hits processed\n", ix, nrow(top)))
}

out_dt <- rbindlist(out_rows, fill = TRUE)
# Mirror cirrhosis schema: rename n_target_held_in -> n_cirrhosis_held_in slot
# is intentionally NOT done -- keep the stage-specific column name so the
# combined table is unambiguous.
out_path <- file.path(OUT_DIR, "sh_dissoc_sensitivity.tsv")
fwrite(out_dt, out_path, sep = "\t")
cat(sprintf("\n[output] %d rows -> %s\n", nrow(out_dt), out_path))

cat("\n[summary] SH classification table:\n")
print(table(out_dt$classification, useNA = "ifany"))
cat(sprintf("  median n_donors_held_in = %.0f\n",
            stats::median(out_dt$n_donors_held_in, na.rm = TRUE)))
cat(sprintf("  median n_%s_held_in   = %.0f\n", tolower(TARGET_LEVEL),
            stats::median(out_dt$n_target_held_in, na.rm = TRUE)))

# ============================================================================
# Combine 347b (Cirrhosis) + 347c (SH) into one table
# ============================================================================
cat("\n[combine] merging 347b (Cirrhosis) + 347c (SH) into all_dissoc_sensitivity.tsv\n")

# Normalise the per-stage donor-count column into a single n_stage_held_in
# column so the two tables stack cleanly.
sh_out <- copy(out_dt)
sh_out[, contrast := "SH_vs_Healthy"]
setnames(sh_out, "n_target_held_in", "n_stage_held_in")

if (file.exists(CIRR_TSV)) {
  cirr_out <- fread(CIRR_TSV)
  cirr_out[, contrast := "Cirrhosis_vs_Healthy"]
  if ("n_cirrhosis_held_in" %in% names(cirr_out)) {
    setnames(cirr_out, "n_cirrhosis_held_in", "n_stage_held_in")
  }
  combined <- rbindlist(list(cirr_out, sh_out), use.names = TRUE, fill = TRUE)
} else {
  warning(sprintf("[combine] %s not found; writing SH-only combined table",
                  CIRR_TSV))
  combined <- sh_out
}

combined_path <- file.path(OUT_DIR, "all_dissoc_sensitivity.tsv")
fwrite(combined, combined_path, sep = "\t")
cat(sprintf("[output] %d rows -> %s\n", nrow(combined), combined_path))

cat("\n[summary] combined classification x contrast table:\n")
print(table(combined$classification, combined$contrast, useNA = "ifany"))

cat("\n[done] SH dissociation sensitivity + combined table written\n")
