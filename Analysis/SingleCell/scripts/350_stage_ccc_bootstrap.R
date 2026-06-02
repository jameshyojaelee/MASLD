#!/usr/bin/env Rscript
# ============================================================================
# 350_stage_ccc_bootstrap.R
#
# Donor-stratified bootstrap confidence intervals for the top-100
# stage-progressive LR pairs (term = disease_stage_coarseSteatohepatitis)
# from Script 346.
#
# For each of B=200 replicates:
#   - resample donors WITH replacement WITHIN each disease_stage_coarse
#     stratum (preserves Healthy ~104, Steatosis ~36, SH ~101, Cirrhosis ~28).
#   - refit the same lmer used by Script 346:
#       score ~ disease_stage_coarse + log10(n_source_cells)
#             + log10(n_target_cells) + dataset + (1|dataset)
#   - store the SH-vs-Healthy estimate.
# Then per LR pair: 2.5% / 50% / 97.5% quantiles + ci_excludes_zero flag.
#
# Output:
#   Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/bootstrap_ci.tsv
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(lme4)
  library(lmerTest)
  library(parallel)
})

# --- Config ----------------------------------------------------------------
N_CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8"))
options(mc.cores = N_CORES)
B <- 200L
TOP_N <- 100L
TARGET_TERM <- "disease_stage_coarseSteatohepatitis"

cat(sprintf("[config] B=%d  TOP_N=%d  target=%s  cores=%d\n",
            B, TOP_N, TARGET_TERM, N_CORES))

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
OUT_DIR <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory")
META_EXT  <- file.path(OUT_DIR, "donor_metadata_extended.tsv")
MERGED_TSV <- file.path(OUT_DIR, "all_donor_lr_scores.tsv.gz")
LMM_COARSE <- file.path(OUT_DIR, "stage_lr_lmm_coarse.tsv")

# --- Pick top-100 LR pairs -------------------------------------------------
cat(sprintf("[input] reading %s\n", LMM_COARSE))
lmm <- fread(LMM_COARSE)
sh_hits <- lmm[term == TARGET_TERM]
cat(sprintf("[input] %d rows for term=%s\n", nrow(sh_hits), TARGET_TERM))

# Order by p-value (ascending) and break ties by descending |Estimate|.
sh_hits[, neg_abs_est := -abs(Estimate)]
setorder(sh_hits, pval, neg_abs_est)
sh_hits[, neg_abs_est := NULL]
top_hits <- sh_hits[seq_len(min(TOP_N, nrow(sh_hits)))]
top_hits[, original_estimate := Estimate]
top_hits[, original_pval := pval]
top_lr_keys <- top_hits[, .(ct_pair, lr_pair, ligand_complex, receptor_complex,
                         original_estimate, original_pval)]
cat(sprintf("[select] picked %d top LR pairs\n", nrow(top_lr_keys)))

# --- Load metadata + LR scores --------------------------------------------
cat(sprintf("[input] reading %s\n", META_EXT))
meta <- fread(META_EXT)
# Backwards-compat: protocol-contamination flag and clean F_stage column may
# be missing if Phase 1 atlas refresh hasn't completed yet.
if (!"exclude_stage_analysis" %in% names(meta))
  meta[, exclude_stage_analysis := FALSE]
if (!"F_stage_augmented_clean" %in% names(meta))
  meta[, F_stage_augmented_clean := if ("F_stage_augmented" %in% names(meta))
       F_stage_augmented else NA_real_]
covar_cols <- intersect(
  c("sample", "dataset", "disease_stage_coarse", "exclude_stage_analysis"),
  names(meta))
meta_use <- meta[, ..covar_cols]
meta_use <- meta_use[!is.na(disease_stage_coarse) &
                     disease_stage_coarse != ""]
stage_levels <- c("Healthy", "Steatosis", "Steatohepatitis", "Cirrhosis")
meta_use[, disease_stage_coarse := factor(disease_stage_coarse,
                                          levels = stage_levels)]
meta_use <- meta_use[!is.na(disease_stage_coarse)]
# Drop protocol-contaminated donors before bootstrap stratification.
n_before_excl <- nrow(meta_use)
meta_use <- meta_use[is.na(exclude_stage_analysis) |
                     exclude_stage_analysis == FALSE]
cat(sprintf("[filter] excluded %d donors flagged exclude_stage_analysis; %d donors remain for bootstrap\n",
            n_before_excl - nrow(meta_use), nrow(meta_use)))
cat(sprintf("[meta] %d donors in scope\n", nrow(meta_use)))
cat("[meta] stage distribution:\n"); print(table(meta_use$disease_stage_coarse))

cat(sprintf("[input] reading %s\n", MERGED_TSV))
lr_long <- fread(MERGED_TSV)
lr_long[, lr_pair := paste(ligand_complex, receptor_complex, sep = "__")]
lr_long[, ct_pair := paste(source, target, sep = "->")]

# Restrict LR rows to the top-100 lr_keys (huge speedup).
lr_long <- lr_long[paste(ct_pair, lr_pair, sep = "||") %in%
                   paste(top_lr_keys$ct_pair, top_lr_keys$lr_pair, sep = "||")]
cat(sprintf("[filter] %d rows kept for top-%d LR pairs\n",
            nrow(lr_long), nrow(top_lr_keys)))

# Score same as 346: -log10(magnitude_rank), floor at 1e-4.
lr_long[, score := -log10(pmax(magnitude_rank, 1e-4))]

# Join only the donor-level covars we need; keep the (already-set) per-row
# n_source_cells / n_target_cells from the per-donor liana output.
lr_long <- merge(lr_long, meta_use, by = "sample", all.x = FALSE)
lr_long[, disease_stage_coarse := factor(disease_stage_coarse,
                                         levels = stage_levels)]
cat(sprintf("[merge] %d donor-LR rows after metadata join\n", nrow(lr_long)))

# --- Build per-LR data shards once -----------------------------------------
lr_long[, lr_key := paste(ct_pair, lr_pair, sep = "||")]
shards <- split(lr_long, by = "lr_key", drop = TRUE)
top_lr_keys[, lr_key := paste(ct_pair, lr_pair, sep = "||")]
shards <- shards[top_lr_keys$lr_key]  # align order
shards <- shards[!sapply(shards, is.null)]
cat(sprintf("[shards] %d LR pairs have donor-LR rows\n", length(shards)))

# --- Donor pool per stratum for resampling ---------------------------------
# Stratified resample is on DONORS, so build a vector of donor IDs per
# stratum once. Each bootstrap replicate samples |Stratum| donors with
# replacement, then we filter each shard down to the resampled donors,
# duplicating rows where a donor was drawn k times.
donors_per_stratum <- split(meta_use$sample, meta_use$disease_stage_coarse)
stratum_sizes <- vapply(donors_per_stratum, length, integer(1))
cat("[bootstrap] stratum sizes: \n"); print(stratum_sizes)

# Per-shard indexer: list of row-indices lr_keyed by donor (sample).
build_donor_index <- function(d) {
  split(seq_len(nrow(d)), d$sample)
}
shard_donor_idx <- lapply(shards, build_donor_index)

# --- Bootstrap fitting function -------------------------------------------
fit_one_replicate <- function(b) {
  set.seed(42L + b)
  # Stratified donor resample
  drawn <- unlist(lapply(donors_per_stratum, function(v) {
    sample(v, size = length(v), replace = TRUE)
  }), use.names = FALSE)
  # Multiplicity table (donor -> count drawn)
  draw_counts <- tabulate(match(drawn, unique(drawn)))
  unique_drawn <- unique(drawn)
  mult <- setNames(draw_counts, unique_drawn)

  out <- numeric(length(shards))
  for (i in seq_along(shards)) {
    d <- shards[[i]]
    idx_by_donor <- shard_donor_idx[[i]]
    # Which donors in this shard appear in the resample?
    keep_donors <- intersect(names(idx_by_donor), unique_drawn)
    if (!length(keep_donors)) { out[i] <- NA_real_; next }
    # Expand rows according to multiplicity
    row_idx <- unlist(lapply(keep_donors, function(s) {
      rep(idx_by_donor[[s]], times = mult[[s]])
    }), use.names = FALSE)
    dd <- d[row_idx]
    # Require all 4 strata present + enough variation
    tab <- table(dd$disease_stage_coarse)
    if (any(tab < 5L)) { out[i] <- NA_real_; next }
    if (length(unique(dd$dataset)) < 2L) { out[i] <- NA_real_; next }
    # Same model as 346 (coarse axis, with random intercept per dataset).
    fit <- tryCatch(
      lmerTest::lmer(
        score ~ disease_stage_coarse +
          log10(pmax(n_source_cells, 1)) +
          log10(pmax(n_target_cells, 1)) +
          dataset + (1 | dataset),
        data = dd, REML = TRUE,
        control = lmerControl(check.conv.singular = .makeCC("ignore", tol = 1e-4),
                              calc.derivs = FALSE)),
      error = function(e) NULL,
      warning = function(w) {
        suppressWarnings(lmerTest::lmer(
          score ~ disease_stage_coarse +
            log10(pmax(n_source_cells, 1)) +
            log10(pmax(n_target_cells, 1)) +
            dataset + (1 | dataset),
          data = dd, REML = TRUE,
          control = lmerControl(check.conv.singular = .makeCC("ignore", tol = 1e-4),
                                calc.derivs = FALSE)))
      })
    if (is.null(fit) || inherits(fit, "try-error")) {
      out[i] <- NA_real_; next
    }
    fe <- tryCatch(lme4::fixef(fit), error = function(e) NULL)
    if (is.null(fe) || !(TARGET_TERM %in% names(fe))) {
      out[i] <- NA_real_; next
    }
    out[i] <- unname(fe[TARGET_TERM])
  }
  out
}

# --- Run B replicates in parallel -----------------------------------------
cat(sprintf("[bootstrap] running B=%d replicates over %d shards on %d cores\n",
            B, length(shards), N_CORES))
t0 <- Sys.time()
if (N_CORES > 1 && .Platform$OS.type != "windows") {
  rep_mat <- parallel::mclapply(seq_len(B), fit_one_replicate,
                                mc.cores = N_CORES,
                                mc.preschedule = TRUE)
} else {
  rep_mat <- lapply(seq_len(B), fit_one_replicate)
}
elapsed_min <- as.numeric(Sys.time() - t0, units = "mins")
cat(sprintf("[bootstrap] done in %.1f min\n", elapsed_min))

# rep_mat is list-of-B; each element a length(shards) numeric vector.
# Stack into a (B x shards) matrix.
boot_mat <- do.call(rbind, lapply(rep_mat, function(v) {
  if (length(v) == length(shards)) v else rep(NA_real_, length(shards))
}))
colnames(boot_mat) <- names(shards)
cat(sprintf("[bootstrap] matrix dim: %d x %d\n",
            nrow(boot_mat), ncol(boot_mat)))

# --- Quantiles + CI flag ---------------------------------------------------
qtable <- data.table(
  lr_key = colnames(boot_mat),
  boot_median = apply(boot_mat, 2L, stats::median, na.rm = TRUE),
  boot_q025   = apply(boot_mat, 2L, stats::quantile, probs = 0.025,
                      na.rm = TRUE, names = FALSE),
  boot_q975   = apply(boot_mat, 2L, stats::quantile, probs = 0.975,
                      na.rm = TRUE, names = FALSE),
  n_replicates = apply(boot_mat, 2L, function(x) sum(!is.na(x)))
)
qtable[, ci_excludes_zero := !is.na(boot_q025) & !is.na(boot_q975) &
         sign(boot_q025) == sign(boot_q975) & sign(boot_q025) != 0]

# Merge back to original lr_keys + estimates
out_dt <- merge(top_lr_keys, qtable, by = "lr_key", all.x = TRUE)
setcolorder(out_dt,
  c("ct_pair", "lr_pair", "ligand_complex", "receptor_complex",
    "original_estimate", "original_pval",
    "boot_median", "boot_q025", "boot_q975",
    "ci_excludes_zero", "n_replicates"))
out_dt[, lr_key := NULL]
setorder(out_dt, original_pval)

OUT_FILE <- file.path(OUT_DIR, "bootstrap_ci.tsv")
fwrite(out_dt, OUT_FILE, sep = "\t")
cat(sprintf("[output] %d rows -> %s\n", nrow(out_dt), OUT_FILE))

# --- Summary ---------------------------------------------------------------
n_exc <- sum(out_dt$ci_excludes_zero, na.rm = TRUE)
ci_width <- out_dt$boot_q975 - out_dt$boot_q025
rel_width <- ci_width / pmax(abs(out_dt$original_estimate), 1e-9)
cat("\n[summary] ------------------------------------------------------\n")
cat(sprintf("  CI-excludes-zero: %d / %d (%.1f%%)\n",
            n_exc, nrow(out_dt), 100 * n_exc / nrow(out_dt)))
cat(sprintf("  Median 95%% CI width: %.4f\n",
            stats::median(ci_width, na.rm = TRUE)))
cat(sprintf("  Median CI width / |estimate|: %.3f\n",
            stats::median(rel_width, na.rm = TRUE)))
cat("\n  Tightest (smallest CI width / |estimate|) top 5:\n")
ord_tight <- order(rel_width, na.last = NA)[seq_len(min(5, length(rel_width)))]
print(out_dt[ord_tight,
  .(ct_pair, lr_pair, original_estimate, boot_q025, boot_q975,
    ci_width = boot_q975 - boot_q025)])
cat("\n  Widest top 5:\n")
ord_wide <- order(rel_width, decreasing = TRUE, na.last = NA)[
  seq_len(min(5, length(rel_width)))]
print(out_dt[ord_wide,
  .(ct_pair, lr_pair, original_estimate, boot_q025, boot_q975,
    ci_width = boot_q975 - boot_q025)])

cat("\n[done] bootstrap CI pass complete\n")
