#!/usr/bin/env Rscript
# ============================================================================
# 07d_chain_v3_final.R
#
# Assemble the final hardened headline LR list. Joins outputs from:
#   - 07b_chain_v3_hardened: stage_lr_lmm_{coarse,fstage_documented,
#       fstage_augmented,continuous}_v3.tsv + lr_bulk_concordance_v3.tsv +
#       loo_replication_rate_per_lr_v3.tsv
#   - 07c_endocrine_lr_filter: endocrine_lr_flags.tsv
#   - 350_v3_bootstrap: bootstrap_ci_v3.tsv
#   - 350b_permutation_null_v3: permutation_fdr_v3.tsv
#   - 351_leverage_diagnostics_v3: leverage_diagnostics_v3.tsv
#   - 353_three_way_concordance: three_way_concordance_v3.tsv (optional)
#
# Headline gates (must satisfy ALL):
#   - q_bonferroni_family < 0.05 in at least one of the 3 primary axes
#       (coarse, documented F-stage, continuous)
#   - bootstrap CI excludes 0 in the coarse axis (where bootstrap was run)
#   - empirical_q < 0.10 from permutation null (coarse axis)
#   - max_leverage_pct < 0.20 (no single donor dominates beta)
#   - replicates in >=2 of 3 method-rankings (three-way) if available
#
# Paracrine-only sub-list additionally requires endocrine_suspect == FALSE.
#
# Output: stage_trajectory_v3/stage_lr_headline_v3.tsv (all-LR ranked table)
#         stage_trajectory_v3/stage_lr_paracrine_headline_v3.tsv (filtered)
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
V3_DIR <- file.path(BASE,
  "Analysis/SingleCell/results_gpu_v2_phase05/ccc/stage_trajectory_v3")
TARGET_TERM <- "disease_stage_coarseSteatohepatitis"
LEVERAGE_THRESHOLD <- 0.20

read_safe <- function(p) {
  if (file.exists(p)) fread(p) else NULL
}

coarse <- read_safe(file.path(V3_DIR, "stage_lr_lmm_coarse_v3.tsv"))
doc <- read_safe(file.path(V3_DIR, "stage_lr_lmm_fstage_documented_v3.tsv"))
aug <- read_safe(file.path(V3_DIR, "stage_lr_lmm_fstage_augmented_v3.tsv"))
cont <- read_safe(file.path(V3_DIR, "stage_lr_lmm_continuous_v3.tsv"))
endo <- read_safe(file.path(V3_DIR, "endocrine_lr_flags.tsv"))
bulk <- read_safe(file.path(V3_DIR, "lr_bulk_concordance_v3.tsv"))
loo <- read_safe(file.path(V3_DIR, "loo_replication_rate_per_lr_v3.tsv"))
boot <- read_safe(file.path(V3_DIR, "bootstrap_ci_v3.tsv"))
perm <- read_safe(file.path(V3_DIR, "permutation_fdr_v3.tsv"))
lev <- read_safe(file.path(V3_DIR, "leverage_diagnostics_v3.tsv"))
threeway <- read_safe(file.path(V3_DIR, "three_way_concordance_v3.tsv"))

cat(sprintf("[load] coarse=%s | doc=%s | aug=%s | cont=%s\n",
            ifelse(is.null(coarse), 0, nrow(coarse)),
            ifelse(is.null(doc), 0, nrow(doc)),
            ifelse(is.null(aug), 0, nrow(aug)),
            ifelse(is.null(cont), 0, nrow(cont))))
cat(sprintf("[load] endo=%s | bulk=%s | loo=%s | boot=%s | perm=%s | lev=%s | 3way=%s\n",
            ifelse(is.null(endo), 0, nrow(endo)),
            ifelse(is.null(bulk), 0, nrow(bulk)),
            ifelse(is.null(loo), 0, nrow(loo)),
            ifelse(is.null(boot), 0, nrow(boot)),
            ifelse(is.null(perm), 0, nrow(perm)),
            ifelse(is.null(lev), 0, nrow(lev)),
            ifelse(is.null(threeway), 0, nrow(threeway))))

if (is.null(coarse) || nrow(coarse) == 0) stop("coarse axis LMM is empty")

# --- Base: coarse SH term, ranked --------------------------------------
sh <- coarse[term == TARGET_TERM]
sh[, lr_key := paste(ct_pair, lr_pair, sep = "||")]
sh[, rank_score := abs(Estimate) * (-log10(pmax(pval, 1e-30)))]
setorder(sh, -rank_score)
sh[, rank_coarse := seq_len(.N)]

# Rename so we know which axis
setnames(sh, c("Estimate", "StdErr", "tval", "pval",
               "q_within_ct", "q_bonferroni_family", "n_donors", "df"),
         c("coarse_Estimate_SH", "coarse_StdErr_SH", "coarse_tval_SH",
           "coarse_pval_SH", "coarse_q_within_ct_SH",
           "coarse_q_bonferroni_family_SH", "coarse_n_donors_SH", "coarse_df_SH"),
         skip_absent = TRUE)

# Strip 'term' (no longer needed) but keep ligand/receptor
sh[, term := NULL]

final <- sh[, .(lr_key, ct_pair, lr_pair, ligand_complex, receptor_complex,
                rank_coarse, rank_score,
                coarse_Estimate_SH, coarse_pval_SH,
                coarse_q_within_ct_SH, coarse_q_bonferroni_family_SH,
                coarse_n_donors_SH)]

# --- Documented F-stage axis -------------------------------------------
if (!is.null(doc) && nrow(doc) > 0) {
  doc_use <- doc[term == "F_stage_doc_numeric",
                 .(ct_pair, lr_pair, Estimate, pval,
                   q_within_ct, q_bonferroni_family, n_donors)]
  setnames(doc_use,
           c("Estimate", "pval", "q_within_ct",
             "q_bonferroni_family", "n_donors"),
           c("doc_Estimate", "doc_pval", "doc_q_within_ct",
             "doc_q_bonferroni_family", "doc_n_donors"))
  final <- merge(final, doc_use, by = c("ct_pair", "lr_pair"), all.x = TRUE)
}

# --- Augmented F-stage axis (sensitivity) ------------------------------
if (!is.null(aug) && nrow(aug) > 0) {
  aug_use <- aug[term == "F_stage_aug_numeric",
                 .(ct_pair, lr_pair, Estimate, pval,
                   q_within_ct, q_bonferroni_family, n_donors)]
  setnames(aug_use,
           c("Estimate", "pval", "q_within_ct",
             "q_bonferroni_family", "n_donors"),
           c("aug_Estimate", "aug_pval", "aug_q_within_ct",
             "aug_q_bonferroni_family", "aug_n_donors"))
  final <- merge(final, aug_use, by = c("ct_pair", "lr_pair"), all.x = TRUE)
}

# --- Continuous (macrophage pseudotime) --------------------------------
if (!is.null(cont) && nrow(cont) > 0) {
  cont_use <- cont[term == "macrophage_pseudotime_mean",
                   .(ct_pair, lr_pair, Estimate, pval,
                     q_within_ct, q_bonferroni_family, n_donors)]
  setnames(cont_use,
           c("Estimate", "pval", "q_within_ct",
             "q_bonferroni_family", "n_donors"),
           c("cont_Estimate", "cont_pval", "cont_q_within_ct",
             "cont_q_bonferroni_family", "cont_n_donors"))
  final <- merge(final, cont_use, by = c("ct_pair", "lr_pair"), all.x = TRUE)
}

# --- Endocrine flag -----------------------------------------------------
if (!is.null(endo) && nrow(endo) > 0) {
  endo_use <- endo[, .(ligand_complex, receptor_complex,
                       endocrine_suspect, ligand_is_hepatokine_curated)]
  final <- merge(final, endo_use,
                 by = c("ligand_complex", "receptor_complex"), all.x = TRUE)
} else {
  final[, endocrine_suspect := FALSE]
  final[, ligand_is_hepatokine_curated := FALSE]
}

# --- Bulk concordance ---------------------------------------------------
if (!is.null(bulk) && nrow(bulk) > 0) {
  bk <- bulk[axis == "Steatohepatitis_vs_Healthy",
             .(ct_pair, lr_pair, both_concordant, lig_concordant, rec_concordant)]
  setnames(bk,
           c("both_concordant", "lig_concordant", "rec_concordant"),
           c("bulk_sh_both_concordant", "bulk_sh_lig_concordant",
             "bulk_sh_rec_concordant"))
  final <- merge(final, bk, by = c("ct_pair", "lr_pair"), all.x = TRUE)
}

# --- LOO replication ----------------------------------------------------
if (!is.null(loo) && nrow(loo) > 0) {
  loo_use <- loo[, .(ct_pair, lr_pair, replication_rate,
                     n_holdouts_run, n_replicating)]
  setnames(loo_use, "replication_rate", "loo_replication_rate")
  final <- merge(final, loo_use, by = c("ct_pair", "lr_pair"), all.x = TRUE)
}

# --- Bootstrap CIs ------------------------------------------------------
if (!is.null(boot) && nrow(boot) > 0) {
  boot_use <- boot[, .(ct_pair, lr_pair,
                       boot_median, boot_q025, boot_q975, ci_excludes_zero)]
  final <- merge(final, boot_use, by = c("ct_pair", "lr_pair"), all.x = TRUE)
}

# --- Permutation null ---------------------------------------------------
if (!is.null(perm) && nrow(perm) > 0) {
  perm_use <- perm[, .(ct_pair, lr_pair, empirical_pval, empirical_q,
                       n_perms_more_extreme, n_perms_valid)]
  final <- merge(final, perm_use, by = c("ct_pair", "lr_pair"), all.x = TRUE)
}

# --- Leverage -----------------------------------------------------------
if (!is.null(lev) && nrow(lev) > 0) {
  lev_use <- lev[, .(ct_pair, lr_pair, max_leverage_pct,
                     max_leverage_donor, max_leverage_dataset,
                     leverage_concerning)]
  final <- merge(final, lev_use, by = c("ct_pair", "lr_pair"), all.x = TRUE)
}

# --- Three-way method concordance --------------------------------------
if (!is.null(threeway) && nrow(threeway) > 0) {
  threeway_use <- threeway[, .(ct_pair, lr_pair,
                               scvi_rank = get("scvi_rank"),
                               harmony_rank = get("harmony_rank"),
                               scanorama_rank = if ("scanorama_rank" %in% names(threeway))
                                 get("scanorama_rank") else NA_integer_,
                               n_methods_in_topN = get("n_methods_in_topN"))]
  final <- merge(final, threeway_use, by = c("ct_pair", "lr_pair"), all.x = TRUE)
}

# --- Headline gates -----------------------------------------------------
# Master review M-P0-5: pmin across k axes is the Tippett min-P statistic.
# Under independence the null distribution is Beta(1, k); applying no
# correction inflates significance ~k-fold. Use Tippett correction per row,
# where the effective k is the number of NON-NA axes contributing to the
# pmin (a pair tested only on 1 of 3 axes is not penalised k=3).
#
# Master review M-P0-11: the DOCUMENTED F-stage axis is fit with lm() (single
# dataset, Andrews-only) while coarse / augmented / continuous are fit with
# lmer(). Combining their p-values in pmin equates incommensurable estimands.
# Default: EXCLUDE documented F-stage from q_min for the headline gate; report
# it separately for support. Override via env var INCLUDE_DOC_IN_QMIN=TRUE
# (only with peer sign-off + matching lmer re-fit).
.include_doc_qmin <- toupper(Sys.getenv("INCLUDE_DOC_IN_QMIN", "FALSE")) == "TRUE"
axis_q_cols <- c(
  "coarse_q_bonferroni_family_SH",
  if (.include_doc_qmin && "doc_q_bonferroni_family" %in% names(final))
    "doc_q_bonferroni_family" else NULL,
  if ("cont_q_bonferroni_family" %in% names(final)) "cont_q_bonferroni_family" else NULL
)
cat(sprintf("[gate] q_min computed over %d axes: %s\n",
            length(axis_q_cols), paste(axis_q_cols, collapse = ", ")))
final[, q_min_bonferroni := do.call(pmin,
                                    c(lapply(axis_q_cols,
                                             function(c) final[[c]]),
                                      list(na.rm = TRUE)))]
final[, n_axes_tested := rowSums(!is.na(final[, ..axis_q_cols]))]
# Tippett: q_tippett = 1 - (1 - q_min)^k. Equivalent to multiplying min-p by
# effective number of axes when q_min is small.
final[, q_min_tippett := pmin(1 - (1 - pmin(q_min_bonferroni, 1))^pmax(n_axes_tested, 1L),
                              1)]
# Documented axis q-value retained as standalone support column (not in q_min).
if ("doc_q_bonferroni_family" %in% names(final)) {
  final[, supports_doc_fstage := !is.na(doc_q_bonferroni_family) &
                                 doc_q_bonferroni_family < 0.05]
}

passes_significance <- function(dt) {
  # Use Tippett-corrected q_min (M-P0-5).
  !is.na(dt$q_min_tippett) & dt$q_min_tippett < 0.05
}
# Master review M-P0-6: gates must return NA (not FALSE) when the underlying
# evidence is missing. The downstream headline check `is.na(gate) | gate`
# requires NA to mean "untested" rather than "failed". `passes_leverage` was
# already correctly NA-tolerant; mirror that pattern for bootstrap / perm /
# threeway. Fix is expected to lift the headline from 5 to ~12 LR pairs.
passes_bootstrap <- function(dt) {
  if (!"ci_excludes_zero" %in% names(dt)) return(rep(NA, nrow(dt)))
  ifelse(is.na(dt$ci_excludes_zero), NA, dt$ci_excludes_zero == TRUE)
}
passes_perm <- function(dt) {
  if (!"empirical_q" %in% names(dt)) return(rep(NA, nrow(dt)))
  ifelse(is.na(dt$empirical_q), NA, dt$empirical_q < 0.10)
}
passes_leverage <- function(dt) {
  if (!"max_leverage_pct" %in% names(dt)) return(rep(NA, nrow(dt)))
  # Leverage: NA means not measured, so do not exclude.
  ifelse(is.na(dt$max_leverage_pct), NA, dt$max_leverage_pct <= LEVERAGE_THRESHOLD)
}
passes_threeway <- function(dt) {
  if (!"n_methods_in_topN" %in% names(dt)) return(rep(NA, nrow(dt)))
  ifelse(is.na(dt$n_methods_in_topN), NA, dt$n_methods_in_topN >= 2)
}

final[, gate_significance := passes_significance(final)]
final[, gate_bootstrap := passes_bootstrap(final)]
final[, gate_permutation := passes_perm(final)]
final[, gate_leverage := passes_leverage(final)]
final[, gate_threeway := passes_threeway(final)]

# Headline = significance AND (bootstrap OR not yet computed) AND
# (permutation OR not yet computed) AND (leverage OR not yet computed) AND
# (threeway OR not yet computed). Missing gates do not exclude.
final[, headline := gate_significance &
                    (is.na(gate_bootstrap) | gate_bootstrap) &
                    (is.na(gate_permutation) | gate_permutation) &
                    (is.na(gate_leverage) | gate_leverage) &
                    (is.na(gate_threeway) | gate_threeway)]
final[is.na(headline), headline := FALSE]
final[, paracrine_headline := headline &
                              (is.na(endocrine_suspect) | endocrine_suspect == FALSE)]

setorder(final, -headline, q_min_bonferroni, coarse_pval_SH)
final[, lr_key := NULL]

OUT_ALL <- file.path(V3_DIR, "stage_lr_headline_v3.tsv")
fwrite(final, OUT_ALL, sep = "\t")
cat(sprintf("[output] %d total LR pairs -> %s\n", nrow(final), OUT_ALL))

OUT_PARA <- file.path(V3_DIR, "stage_lr_paracrine_headline_v3.tsv")
fwrite(final[paracrine_headline == TRUE], OUT_PARA, sep = "\t")
cat(sprintf("[output] %d paracrine-headline LR pairs -> %s\n",
            sum(final$paracrine_headline, na.rm = TRUE), OUT_PARA))

cat("\n[summary]\n")
cat(sprintf("  total LR pairs ranked: %d\n", nrow(final)))
cat(sprintf("  passes_significance: %d\n",
            sum(final$gate_significance, na.rm = TRUE)))
cat(sprintf("  passes all gates (headline): %d\n",
            sum(final$headline, na.rm = TRUE)))
cat(sprintf("  paracrine_headline: %d\n",
            sum(final$paracrine_headline, na.rm = TRUE)))
cat(sprintf("  endocrine_suspect among headline: %d\n",
            sum(final$headline & final$endocrine_suspect, na.rm = TRUE)))
cat("[done]\n")
