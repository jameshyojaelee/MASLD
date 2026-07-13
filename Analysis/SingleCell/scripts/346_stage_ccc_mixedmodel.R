#!/usr/bin/env Rscript
# ============================================================================
# 346_stage_ccc_mixedmodel.R
#
# Mixed-effects regression of per-donor LIANA magnitude_rank against three
# stage axes:
#   (a) disease_stage_coarse  (4-level categorical: Healthy, Steatosis, SH, Cirrhosis)
#   (b) F_stage_documented    (ordered, NA where unavailable)
#   (c) macrophage_pseudotime_mean (continuous)
#
# Model per (source, target, ligand_complex, receptor_complex):
#   -log10(magnitude_rank) ~ stage_axis + log10(n_source_cells)
#                          + log10(n_target_cells) + dataset
#                          + (1 | dataset)
# (age, sex added when bulk-overlay populated them; else dropped.)
#
# Multiple testing:
#   - BH within each (source, target) pair across LR pairs
#   - Bonferroni at the family of cell-type pairs at alpha = 0.05
#
# Outputs:
#   stage_lr_lmm_coarse.tsv
#   stage_lr_lmm_fstage.tsv
#   stage_lr_lmm_continuous.tsv
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(lme4)
  library(lmerTest)
  library(parallel)
})

N_CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8"))
options(mc.cores = N_CORES)
cat(sprintf("[parallel] using %d cores for lmer batch\n", N_CORES))

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
OUT_DIR <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory")

# ---------------------------------------------------------------------------
# Pseudoreplication fix (2026-07-12): the integrated atlas obs["sample"] is a
# sequencing RUN, not a biological donor, for GSE244832 / GSE185477 / GSE202379.
# When DONOR_COLLAPSE=TRUE (default) we (i) read the donor-collapsed LR scores
# produced by the patched 345 (one LIANA run per TRUE donor, cells pooled across
# runs), (ii) collapse the run-level donor_metadata_extended to true-donor level,
# and (iii) route the F-stage axis through F_stage_documented (Andrews SAF
# scores) rather than the leaked scVI F_stage_inferred. Outputs go to a separate
# subdir so the run-level (buggy) outputs are preserved for before/after diff.
# Set DONOR_COLLAPSE=FALSE to exactly reproduce the pre-fix run-level behavior.
# ---------------------------------------------------------------------------
DONOR_COLLAPSE <- toupper(Sys.getenv("DONOR_COLLAPSE", "TRUE")) == "TRUE"
OUT_WRITE_DIR <- Sys.getenv("CCC_OUT_DIR",
  if (DONOR_COLLAPSE) file.path(OUT_DIR, "donor_collapsed") else OUT_DIR)
dir.create(OUT_WRITE_DIR, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("[config] DONOR_COLLAPSE=%s | writing outputs to %s\n",
            DONOR_COLLAPSE, OUT_WRITE_DIR))

# Master review M-P0-7: enforce the pre-registered bootstrap gate. If 343q
# wrote BOOTSTRAP_FALLBACK_REQUIRED.flag, F_stage_augmented_v2 failed
# stability acceptance and must NOT be used as a primary axis.
#
# Updated 2026-05-22 (protocol-contamination remediation): the gate now
# WARNS rather than HALTs, because the F-stage column selector downstream
# respects the flag (uses F_stage_inferred when flag exists, unless explicit
# ALLOW_F_STAGE_AUGMENTED_DESPITE_FLAG=TRUE override). This lets the LMM run
# on F_stage_inferred while keeping the bootstrap-gate semantics intact.
.boot_flag <- file.path(OUT_DIR, "BOOTSTRAP_FALLBACK_REQUIRED.flag")
.override <- toupper(Sys.getenv("ALLOW_F_STAGE_AUGMENTED_DESPITE_FLAG", "FALSE"))
if (file.exists(.boot_flag) && .override != "TRUE") {
  msg <- paste(readLines(.boot_flag), collapse = " | ")
  warning("[bootstrap gate] Pre-registered bootstrap gate FAILED: ", msg,
          "\n  Falling back to F_stage_inferred (343b ordinal) for F-stage axis.",
          "\n  To force F_stage_augmented despite the gate, set",
          " ALLOW_F_STAGE_AUGMENTED_DESPITE_FLAG=TRUE.")
}
if (file.exists(.boot_flag) && .override == "TRUE") {
  warning("[bootstrap gate] Override active: F_stage_augmented retained despite ",
          "pre-reg failure. Document this in methods.")
}
PER_DONOR_DIR <- file.path(OUT_DIR, "per_donor_lr")
META_EXT <- file.path(OUT_DIR, "donor_metadata_extended.tsv")

# Pre-registered gates (env-overridable for sensitivity runs after donor-collapse
# shrinks the effective N; defaults preserve the original thresholds).
MIN_DONORS_PER_LR <- as.integer(Sys.getenv("MIN_DONORS_PER_LR", "30"))
MIN_DONORS_PER_STRATUM <- as.integer(Sys.getenv("MIN_DONORS_PER_STRATUM", "10"))
cat(sprintf("[gates] MIN_DONORS_PER_LR=%d MIN_DONORS_PER_STRATUM=%d\n",
            MIN_DONORS_PER_LR, MIN_DONORS_PER_STRATUM))

cat(sprintf("[input] reading donor metadata from %s\n", META_EXT))
meta <- fread(META_EXT)
cat(sprintf("[input] %d run-level rows in donor_metadata_extended\n", nrow(meta)))

# --- Collapse run-level metadata to TRUE-donor level ------------------------
# Phenotype/technical-invariant columns (dataset, stage, F-stage, age, sex,
# exclude flag) are constant within a biological donor -> take first non-NA
# (with a consistency check for the load-bearing stage + documented-F-stage).
# Cell-count-weighted columns (frac_Hepatocytes, pseudotime, progressor_frac)
# are re-derived as an n_cells-weighted mean, which for frac_Hepatocytes equals
# the pooled hepatocyte fraction the donor-pooled LIANA run actually used.
if (DONOR_COLLAPSE) {
  source(file.path(BASE, "Analysis/SingleCell/scripts/lib_donor_collapse.R"))
  srr2donor <- build_srr_to_donor_map(BASE)
  first_non_na <- function(x) {
    idx <- which(!is.na(x))
    if (length(idx)) x[idx[1]] else x[NA_integer_]  # typed NA preserves class
  }
  wmean <- function(x, w) {
    ok <- !is.na(x) & !is.na(w) & w > 0
    if (!any(ok)) return(NA_real_)
    sum(x[ok] * w[ok]) / sum(w[ok])
  }
  # Restrict to columns 346 actually consumes (+ n_cells weight) before collapse
  # to avoid coercion edge cases on the 70-column extended table.
  wt_cols  <- c("frac_Hepatocytes", "macrophage_pseudotime_mean",
                "hepatocyte_pseudotime_mean", "progressor_frac")
  inv_cols <- c("dataset", "disease_stage_coarse", "disease_stage_numeric",
                "F_stage_documented", "F_stage_inferred", "F_stage_augmented",
                "F_stage_augmented_clean", "F_stage_source",
                "exclude_stage_analysis", "age", "sex_numeric")
  inv_cols <- intersect(inv_cols, names(meta))
  wt_cols  <- intersect(wt_cols,  names(meta))
  weight_col <- if ("n_cells" %in% names(meta)) "n_cells" else NULL
  meta <- meta[, c("sample", inv_cols, wt_cols,
                   if (!is.null(weight_col)) weight_col), with = FALSE]
  meta[, biological_donor := ifelse(sample %in% names(srr2donor),
                                    srr2donor[sample], sample)]
  # Consistency check on the two load-bearing donor-invariant phenotypes
  for (cc in intersect(c("disease_stage_coarse", "F_stage_documented"), inv_cols)) {
    nbad <- meta[, uniqueN(stats::na.omit(get(cc))), by = biological_donor][V1 > 1, .N]
    if (nbad > 0)
      warning(sprintf("[donor-collapse] %d donor(s) disagree on '%s' across runs",
                      nbad, cc))
  }
  meta_inv <- meta[, lapply(.SD, first_non_na), by = biological_donor,
                   .SDcols = inv_cols]
  if (length(wt_cols) && !is.null(weight_col)) {
    meta_wt <- unique(meta[, .(biological_donor)])
    for (cn in wt_cols) {
      tmp <- meta[, .(v = wmean(get(cn), get(weight_col))), by = biological_donor]
      setnames(tmp, "v", cn)
      meta_wt <- merge(meta_wt, tmp, by = "biological_donor")
    }
    meta <- merge(meta_inv, meta_wt, by = "biological_donor")
  } else {
    meta <- meta_inv
  }
  setnames(meta, "biological_donor", "sample")
  cat(sprintf("[donor-collapse] collapsed metadata to %d true donors\n",
              nrow(meta)))
}
cat(sprintf("[input] %d donors\n", nrow(meta)))

# Load merged TSV (produced by Script 345b from per-donor parquets) --------
MERGED_TSV <- Sys.getenv("LR_SCORES_TSV",
  file.path(OUT_DIR, if (DONOR_COLLAPSE) "all_donor_lr_scores_donorcollapsed.tsv.gz"
                     else "all_donor_lr_scores.tsv.gz"))
if (!file.exists(MERGED_TSV)) {
  stop(paste("Merged LR scores TSV not found at", MERGED_TSV,
             "- run Script 345b first"))
}
lr_long <- fread(MERGED_TSV)
cat(sprintf("[input] %d donor x LR-pair rows loaded from %s\n",
            nrow(lr_long), MERGED_TSV))
# The donor-collapsed LR table already carries `dataset` + `biological_donor`
# (written by 345). Drop them so the authoritative copies come from the
# collapsed metadata join below (avoids a dataset.x/.y collision).
if (DONOR_COLLAPSE) {
  drop_dup <- intersect(c("dataset", "biological_donor"), names(lr_long))
  if (length(drop_dup)) lr_long[, (drop_dup) := NULL]
}

# Backwards-compat: protocol-contamination flag and clean F_stage column may
# be missing if Phase 1 atlas refresh hasn't completed yet.
if (!"exclude_stage_analysis" %in% names(meta))
  meta[, exclude_stage_analysis := FALSE]
if (!"F_stage_augmented_clean" %in% names(meta))
  meta[, F_stage_augmented_clean := F_stage_augmented]

# Join donor metadata into the long table -----------------------------------
covar_cols <- c("sample", "dataset", "disease_stage_coarse",
                "disease_stage_numeric", "F_stage_documented",
                "F_stage_inferred", "F_stage_augmented",
                "F_stage_augmented_clean", "F_stage_source",
                "macrophage_pseudotime_mean", "progressor_frac",
                "frac_Hepatocytes", "exclude_stage_analysis",
                "age", "sex_numeric")
covar_cols <- intersect(covar_cols, names(meta))
lr_long <- merge(lr_long, meta[, ..covar_cols],
                 by = "sample", all.x = TRUE)

# Numeric transform of the response (LIANA magnitude_rank: lower=stronger)
lr_long[, score := -log10(pmax(magnitude_rank, 1e-4))]

# Fix factor levels for disease_stage_coarse with Healthy as reference -----
stage_levels <- c("Healthy", "Steatosis", "Steatohepatitis", "Cirrhosis")
lr_long[, disease_stage_coarse := factor(disease_stage_coarse,
                                         levels = stage_levels)]

# Drop HCC + rows where coarse stage is NA
lr_long <- lr_long[!is.na(disease_stage_coarse)]

# Drop protocol-contaminated donors (GSE136103, Liver_Atlas) flagged in
# Phase 1 atlas refresh. Belt-and-suspenders: the merged parquet should
# already exclude them, but a script-level guard is cheap.
n_before_excl <- nrow(lr_long)
n_donors_before_excl <- uniqueN(lr_long$sample)
lr_long <- lr_long[is.na(exclude_stage_analysis) | exclude_stage_analysis == FALSE]
cat(sprintf("[filter] excluded %d rows (%d donors) flagged exclude_stage_analysis; %d rows (%d donors) remain\n",
            n_before_excl - nrow(lr_long),
            n_donors_before_excl - uniqueN(lr_long$sample),
            nrow(lr_long), uniqueN(lr_long$sample)))

# Define LR pair key
lr_long[, lr_pair := paste(ligand_complex, receptor_complex, sep = "__")]
lr_long[, ct_pair := paste(source, target, sep = "->")]

# Pseudoreplication guard: after donor-collapse every (donor, ct_pair, lr_pair)
# must be unique -- i.e. each true donor contributes at most ONE observation to
# each LR-pair regression. If this fails, run-level rows leaked through and the
# mixed model would be pseudoreplicated again.
if (DONOR_COLLAPSE) {
  dup <- lr_long[, .N, by = .(sample, ct_pair, lr_pair)][N > 1]
  if (nrow(dup) > 0) {
    stop(sprintf("[pseudorep-guard] FAIL: %d (donor,ct_pair,lr_pair) keys are duplicated -- donor-collapse did not take", nrow(dup)))
  }
  cat(sprintf("[pseudorep-guard] PASS: %d donors, one row per (donor, ct_pair, lr_pair)\n",
              uniqueN(lr_long$sample)))
}

# Helper: fit model for one (ct_pair, lr_pair) on one stage axis ------------
# Reviewer-response (2026-05-22): for LR pairs where source OR target is
# "Hepatocytes", add frac_Hepatocytes as a fixed-effect covariate to absorb
# the residual 0%-vs-82% compositional variation that a dataset random
# intercept alone cannot fully account for.
fit_one <- function(d, axis_term, axis_levels = NULL, is_hep_pair = FALSE) {
  if (nrow(d) < MIN_DONORS_PER_LR) return(NULL)
  # Build formula dynamically based on what covariates have variance
  # Compute use_random FIRST so dataset is not added to both fixed and random effects
  use_random <- length(unique(d$dataset)) > 1 && length(unique(d$dataset)) >= 3
  fixed <- c(axis_term)
  if (!use_random && length(unique(d$dataset)) == 2L) fixed <- c(fixed, "dataset")
  fixed <- c(fixed,
             "log10(pmax(n_source_cells,1))",
             "log10(pmax(n_target_cells,1))")
  if (is_hep_pair && "frac_Hepatocytes" %in% names(d) &&
      length(unique(stats::na.omit(d$frac_Hepatocytes))) > 1 &&
      sum(!is.na(d$frac_Hepatocytes)) >= 10) {
    fixed <- c(fixed, "frac_Hepatocytes")
  }
  if (length(unique(stats::na.omit(d$age))) > 1 && sum(!is.na(d$age)) >= 10) {
    fixed <- c(fixed, "age")
  }
  if (length(unique(stats::na.omit(d$sex_numeric))) > 1 &&
      sum(!is.na(d$sex_numeric)) >= 10) {
    fixed <- c(fixed, "sex_numeric")
  }
  rhs <- paste(fixed, collapse = " + ")
  if (use_random) {
    rhs <- paste(rhs, "+ (1 | dataset)")
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
  # Column naming depends on whether lm (4 stat cols: Est, SE, t, p) or
  # lmerTest::lmer (5 stat cols: Est, SE, df, t, p) was used.
  ncol_stats <- ncol(out) - 1L  # subtract the rowname/term column
  if (ncol_stats == 5L) {
    setnames(out, c("term", "Estimate", "StdErr", "df", "tval", "pval"))
  } else if (ncol_stats == 4L) {
    setnames(out, c("term", "Estimate", "StdErr", "tval", "pval"))
    out[, df := NA_real_]
  } else {
    return(NULL)
  }
  out[, n_donors := nrow(d)]
  out
}

# Iterate over (ct_pair, lr_pair) for each axis -----------------------------
fit_axis <- function(axis_term, only_with_value = NULL) {
  cat(sprintf("\n[axis] fitting %s\n", axis_term))
  d_in <- lr_long
  if (!is.null(only_with_value)) {
    d_in <- d_in[!is.na(get(only_with_value))]
  }
  if (nrow(d_in) == 0) return(data.table())

  splits <- split(d_in, by = c("ct_pair", "lr_pair"), drop = TRUE)
  cat(sprintf("[axis] %d (ct_pair, lr_pair) groups to fit (parallel mc=%d)\n",
              length(splits), N_CORES))

  fit_wrapper <- function(d) {
    if (axis_term == "disease_stage_coarse") {
      tab <- table(d$disease_stage_coarse)
      if (any(tab < MIN_DONORS_PER_STRATUM)) return(NULL)
    }
    # Reviewer-response: add frac_Hepatocytes covariate when source or target
    # is Hepatocytes to absorb between-dataset compositional variance.
    is_hep_pair <- isTRUE(d$source[1] == "Hepatocytes") ||
                   isTRUE(d$target[1] == "Hepatocytes")
    fit_res <- fit_one(d, axis_term, is_hep_pair = is_hep_pair)
    if (is.null(fit_res)) return(NULL)
    fit_res[, ct_pair := d$ct_pair[1]]
    fit_res[, lr_pair := d$lr_pair[1]]
    fit_res[, source  := d$source[1]]
    fit_res[, target  := d$target[1]]
    fit_res[, ligand_complex   := d$ligand_complex[1]]
    fit_res[, receptor_complex := d$receptor_complex[1]]
    fit_res[, hep_pair_adjusted := is_hep_pair]
    fit_res
  }

  t0 <- Sys.time()
  if (N_CORES > 1 && .Platform$OS.type != "windows") {
    res_list <- parallel::mclapply(splits, fit_wrapper, mc.cores = N_CORES)
  } else {
    res_list <- lapply(splits, fit_wrapper)
  }
  cat(sprintf("[axis] %s lmer batch took %.1f min\n", axis_term,
              as.numeric(Sys.time() - t0, units = "mins")))

  out <- rbindlist(res_list[!sapply(res_list, is.null)], fill = TRUE)
  if (nrow(out) == 0) return(out)
  # BH within each (ct_pair) for each axis-term
  out[, padj_within_ct := stats::p.adjust(pval, method = "BH"),
      by = .(ct_pair, term)]
  # Family-Bonferroni across ct_pair families
  n_families <- length(unique(out$ct_pair))
  out[, family_bonferroni := pmin(pval * n_families, 1)]
  out
}

# Axis (a): disease_stage_coarse (3 contrasts vs Healthy)
res_coarse <- fit_axis("disease_stage_coarse")
fwrite(res_coarse, file.path(OUT_WRITE_DIR, "stage_lr_lmm_coarse.tsv"), sep = "\t")
cat(sprintf("[output] %d rows -> %s/stage_lr_lmm_coarse.tsv\n",
            nrow(res_coarse), OUT_WRITE_DIR))

# Axis (b): F-stage (ordered, treated as numeric for slope).
# Mega-review A7.2 (2026-06-13): F_stage_augmented is LEAKED — its scVI QWK of
# 0.74-0.76 collapses to jackknife 0.286 / held-out Andrews 0.0 because the
# augmented-anchor training set was screened with leak. The WS3 stage-CCC
# cascade therefore routes the F-stage axis through F_stage_inferred (343b
# ordinal, behind its pre-registered bootstrap gate from 343q), NOT the
# augmented column. F_stage_augmented(_clean) is no longer eligible as the
# primary axis here, regardless of BOOTSTRAP_FALLBACK_REQUIRED.flag presence or
# the ALLOW_F_STAGE_AUGMENTED_DESPITE_FLAG override (both retained upstream only
# for the producer's own provenance; they no longer feed this cascade).
#   - Donor-level protocol contamination (GSE136103 + Liver_Atlas) is handled
#     separately by the exclude_stage_analysis filter upstream of this block.
# DONOR_COLLAPSE mode routes the F-stage axis through F_stage_DOCUMENTED
# (Andrews GSE202379 SAF scores) ONLY -- never the leaked scVI F_stage_inferred,
# per the project operational rule (restrict F-stage to the documented Andrews
# donors; scVI transfer is unlabeled-untestable cross-cohort). This is a
# single-cohort axis, so it fits a plain lm on true donors (no dataset random
# effect). It is deliberately underpowered-honest: N is small after collapse.
if (DONOR_COLLAPSE) {
  fstage_col <- if ("F_stage_documented" %in% names(lr_long) &&
                    sum(!is.na(lr_long$F_stage_documented)) > 0)
    "F_stage_documented" else NULL
  fstage_out <- "stage_lr_lmm_fstage_documented.tsv"
} else if ("F_stage_inferred" %in% names(lr_long) &&
    sum(!is.na(lr_long$F_stage_inferred)) >=
      sum(!is.na(lr_long$F_stage_documented))) {
  fstage_col <- "F_stage_inferred"
  fstage_out <- "stage_lr_lmm_fstage.tsv"
} else if ("F_stage_documented" %in% names(lr_long) &&
           sum(!is.na(lr_long$F_stage_documented)) > 0) {
  fstage_col <- "F_stage_documented"
  fstage_out <- "stage_lr_lmm_fstage.tsv"
} else { fstage_col <- NULL; fstage_out <- NULL }
if (!is.null(fstage_col)) {
  n_fstage_donors <- uniqueN(lr_long[!is.na(get(fstage_col))]$sample)
  cat(sprintf("[axis] F-stage column = %s (%d non-NA rows, %d true donors)\n",
              fstage_col, sum(!is.na(lr_long[[fstage_col]])), n_fstage_donors))
  lr_long[, F_stage_numeric := as.numeric(get(fstage_col))]
  res_fstage <- fit_axis("F_stage_numeric", only_with_value = "F_stage_numeric")
  fwrite(res_fstage, file.path(OUT_WRITE_DIR, fstage_out), sep = "\t")
  cat(sprintf("[output] %d rows -> %s/%s\n",
              nrow(res_fstage), OUT_WRITE_DIR, fstage_out))
} else {
  cat("[axis] no documented F-stage; skipping F-stage axis\n")
}

# Axis (c): macrophage_pseudotime_mean (continuous)
res_cont <- fit_axis("macrophage_pseudotime_mean",
                     only_with_value = "macrophage_pseudotime_mean")
fwrite(res_cont, file.path(OUT_WRITE_DIR, "stage_lr_lmm_continuous.tsv"), sep = "\t")
cat(sprintf("[output] %d rows -> %s/stage_lr_lmm_continuous.tsv\n",
            nrow(res_cont), OUT_WRITE_DIR))

cat("\n[done] mixed-effects fits complete\n")
