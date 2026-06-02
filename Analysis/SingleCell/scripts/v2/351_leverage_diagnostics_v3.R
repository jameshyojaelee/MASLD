#!/usr/bin/env Rscript
# ============================================================================
# 351_leverage_diagnostics_v3.R
#
# Addresses critique #8: per-donor LMMs at n=30 with no leverage diagnostics.
# For each headline LR pair (top-100 SH from coarse axis), refit the lmer
# and compute per-donor diagnostics:
#   - hatvalues equivalent (donor-level influence)
#   - dfbeta-like contribution to disease_stage_coarseSteatohepatitis estimate
#   - Cook's-distance equivalent
#
# Use a leave-one-donor-out refit strategy as the principal diagnostic because
# `cooks.distance` is not directly defined for lmerMod. For each donor d in
# each LR pair shard:
#   - refit the LMM without donor d
#   - delta_beta_d = beta_full - beta_{-d}
#   - leverage_pct = abs(delta_beta_d) / sum_{d'} abs(delta_beta_{d'})
# Then per LR pair: max_leverage_pct, donor_with_max_leverage.
#
# Flag LR pairs where any single donor contributes >20% of the total
# |delta_beta|.
#
# Output: stage_trajectory_v3/leverage_diagnostics_v3.tsv
#   Columns: ct_pair, lr_pair, n_donors, max_leverage_pct, max_leverage_donor,
#            max_leverage_dataset, leverage_concerning (boolean)
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(lme4)
  library(parallel)
})

N_CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "16"))
options(mc.cores = N_CORES)
TOP_N <- 100L
TARGET_TERM <- "disease_stage_coarseSteatohepatitis"
LEVERAGE_THRESHOLD <- 0.20
SMALL_COHORTS <- c("GSE174748", "GSE189600")

cat(sprintf("[config] TOP_N=%d  target=%s  cores=%d  threshold=%.2f\n",
            TOP_N, TARGET_TERM, N_CORES, LEVERAGE_THRESHOLD))

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
V2_DIR <- file.path(BASE,
  "Analysis/SingleCell/results_gpu_v2_phase05/ccc/stage_trajectory_v2")
V3_DIR <- file.path(BASE,
  "Analysis/SingleCell/results_gpu_v2_phase05/ccc/stage_trajectory_v3")
META_V2 <- file.path(BASE,
  "Analysis/SingleCell/results_gpu_v2_phase05/mcp/inputs/donor_metadata_v2.tsv")
MERGED_TSV <- file.path(V2_DIR, "all_donor_lr_scores_v2.tsv.gz")
LMM_COARSE <- file.path(V3_DIR, "stage_lr_lmm_coarse_v3.tsv")

if (!file.exists(LMM_COARSE)) stop("Missing v3 LMM at ", LMM_COARSE)

lmm <- fread(LMM_COARSE)
sh_hits <- lmm[term == TARGET_TERM]
sh_hits[, neg_abs_est := -abs(Estimate)]
setorder(sh_hits, pval, neg_abs_est)
sh_hits[, neg_abs_est := NULL]
top_hits <- sh_hits[seq_len(min(TOP_N, nrow(sh_hits)))]
top_hits[, lr_key := paste(ct_pair, lr_pair, sep = "||")]
cat(sprintf("[input] top-%d LR pairs\n", nrow(top_hits)))

meta <- fread(META_V2)
meta[, dataset_pooled := ifelse(dataset %in% SMALL_COHORTS, "Other_small",
                                as.character(dataset))]
meta_use <- meta[, .(sample, dataset, dataset_pooled, disease_stage_coarse)]
stage_levels <- c("Healthy", "Steatosis", "Steatohepatitis", "Cirrhosis")
meta_use <- meta_use[!is.na(disease_stage_coarse) & disease_stage_coarse != ""]
meta_use[, disease_stage_coarse := factor(disease_stage_coarse, levels = stage_levels)]

lr_long <- fread(MERGED_TSV)
lr_long[, lr_pair := paste(ligand_complex, receptor_complex, sep = "__")]
lr_long[, ct_pair := paste(source, target, sep = "->")]
lr_long[, lr_key := paste(ct_pair, lr_pair, sep = "||")]
lr_long <- lr_long[lr_key %in% top_hits$lr_key]
lr_long[, score := -log10(pmax(magnitude_rank, 1e-4))]
lr_long <- merge(lr_long, meta_use, by = "sample")
shards <- split(lr_long, by = "lr_key", drop = TRUE)
shards <- shards[top_hits$lr_key]
shards <- shards[!sapply(shards, is.null)]
cat(sprintf("[shards] %d LR pairs ready\n", length(shards)))

fit_lmm <- function(d) {
  use_random <- length(unique(d$dataset_pooled)) >= 3L
  formula_str <- "score ~ disease_stage_coarse + log10(pmax(n_source_cells,1)) + log10(pmax(n_target_cells,1)) + dataset_pooled"
  if (use_random) formula_str <- paste(formula_str, "+ (1|dataset_pooled)")
  fit <- tryCatch(
    if (use_random) {
      lme4::lmer(stats::as.formula(formula_str), data = d, REML = TRUE,
                 control = lmerControl(check.conv.singular = .makeCC("ignore", tol = 1e-4),
                                       calc.derivs = FALSE))
    } else {
      stats::lm(stats::as.formula(formula_str), data = d)
    },
    error = function(e) NULL,
    warning = function(w) NULL)
  if (is.null(fit)) return(NA_real_)
  fe <- if (use_random) tryCatch(lme4::fixef(fit), error = function(e) NULL)
        else tryCatch(stats::coef(fit), error = function(e) NULL)
  if (is.null(fe) || !(TARGET_TERM %in% names(fe))) return(NA_real_)
  unname(fe[TARGET_TERM])
}

leverage_one_shard <- function(d) {
  beta_full <- fit_lmm(d)
  if (is.na(beta_full)) {
    return(data.table(lr_key = d$lr_key[1],
                      n_donors = length(unique(d$sample)),
                      beta_full = NA_real_,
                      max_leverage_pct = NA_real_,
                      max_leverage_donor = NA_character_,
                      max_leverage_dataset = NA_character_))
  }
  donors <- unique(d$sample)
  deltas <- numeric(length(donors))
  names(deltas) <- donors
  for (s in donors) {
    d_loo <- d[sample != s]
    if (length(unique(d_loo$disease_stage_coarse)) < 2L) {
      deltas[s] <- NA_real_; next
    }
    beta_loo <- fit_lmm(d_loo)
    if (is.na(beta_loo)) { deltas[s] <- NA_real_; next }
    deltas[s] <- beta_full - beta_loo
  }
  abs_total <- sum(abs(deltas), na.rm = TRUE)
  if (abs_total == 0 || is.na(abs_total)) {
    return(data.table(lr_key = d$lr_key[1],
                      n_donors = length(donors),
                      beta_full = beta_full,
                      max_leverage_pct = 0,
                      max_leverage_donor = NA_character_,
                      max_leverage_dataset = NA_character_))
  }
  leverage_pct <- abs(deltas) / abs_total
  i_max <- which.max(leverage_pct)
  donor_max <- donors[i_max]
  ds_max <- unique(d[sample == donor_max, dataset])[1]
  data.table(lr_key = d$lr_key[1],
             n_donors = length(donors),
             beta_full = beta_full,
             max_leverage_pct = leverage_pct[i_max],
             max_leverage_donor = donor_max,
             max_leverage_dataset = as.character(ds_max))
}

cat(sprintf("[leverage] running LOO per shard on %d cores\n", N_CORES))
t0 <- Sys.time()
res_list <- parallel::mclapply(shards, leverage_one_shard, mc.cores = N_CORES)
cat(sprintf("[leverage] done in %.1f min\n",
            as.numeric(Sys.time() - t0, units = "mins")))

out_dt <- rbindlist(res_list)
out_dt[, ct_pair := sub("\\|\\|.*", "", lr_key)]
out_dt[, lr_pair := sub("^[^|]+\\|\\|", "", lr_key)]
out_dt[, leverage_concerning := !is.na(max_leverage_pct) &
                                 max_leverage_pct > LEVERAGE_THRESHOLD]
setcolorder(out_dt,
  c("ct_pair", "lr_pair", "n_donors", "beta_full",
    "max_leverage_pct", "max_leverage_donor", "max_leverage_dataset",
    "leverage_concerning"))
out_dt[, lr_key := NULL]
setorder(out_dt, -max_leverage_pct)
OUT <- file.path(V3_DIR, "leverage_diagnostics_v3.tsv")
fwrite(out_dt, OUT, sep = "\t")
cat(sprintf("[output] %d rows -> %s\n", nrow(out_dt), OUT))
n_concern <- sum(out_dt$leverage_concerning, na.rm = TRUE)
cat(sprintf("[summary] %d / %d (%.1f%%) LR pairs have a donor contributing >%.0f%% of |delta_beta|\n",
            n_concern, nrow(out_dt), 100 * n_concern / nrow(out_dt),
            100 * LEVERAGE_THRESHOLD))
cat("[done] leverage diagnostics complete\n")
