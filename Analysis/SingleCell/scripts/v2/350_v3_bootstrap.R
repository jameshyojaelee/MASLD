#!/usr/bin/env Rscript
# ============================================================================
# 350_v3_bootstrap.R
#
# v3 of 350_stage_ccc_bootstrap.R. Donor-stratified bootstrap CIs for the
# top-100 stage-progressive LR pairs from the coarse axis SH term, fit with
# the same model as 07b_chain_v3_hardened (with dataset_pooled random effect).
#
# B=500 replicates; output: stage_trajectory_v3/bootstrap_ci_v3.tsv.
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(lme4)
  library(lmerTest)
  library(parallel)
})

N_CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "16"))
options(mc.cores = N_CORES)
B <- as.integer(Sys.getenv("BOOTSTRAP_B", "500"))
TOP_N <- 100L
TARGET_TERM <- "disease_stage_coarseSteatohepatitis"
SMALL_COHORTS <- c("GSE174748", "GSE189600")

cat(sprintf("[config] B=%d  TOP_N=%d  target=%s  cores=%d\n",
            B, TOP_N, TARGET_TERM, N_CORES))

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
V2_DIR <- file.path(BASE,
  "Analysis/SingleCell/results_gpu_v2_phase05/ccc/stage_trajectory_v2")
V3_DIR <- file.path(BASE,
  "Analysis/SingleCell/results_gpu_v2_phase05/ccc/stage_trajectory_v3")
dir.create(V3_DIR, recursive = TRUE, showWarnings = FALSE)

META_V2 <- file.path(BASE,
  "Analysis/SingleCell/results_gpu_v2_phase05/mcp/inputs/donor_metadata_v2.tsv")
MERGED_TSV <- file.path(V2_DIR, "all_donor_lr_scores_v2.tsv.gz")
LMM_COARSE <- file.path(V3_DIR, "stage_lr_lmm_coarse_v3.tsv")

# --- Pick top-100 LR pairs from v3 coarse LMM -----------------------------
if (!file.exists(LMM_COARSE)) {
  stop("v3 LMM coarse not found at ", LMM_COARSE,
       " — run 07b_chain_v3_hardened first")
}
cat(sprintf("[input] reading %s\n", LMM_COARSE))
lmm <- fread(LMM_COARSE)
sh_hits <- lmm[term == TARGET_TERM]
cat(sprintf("[input] %d rows for term=%s\n", nrow(sh_hits), TARGET_TERM))

sh_hits[, neg_abs_est := -abs(Estimate)]
setorder(sh_hits, pval, neg_abs_est)
sh_hits[, neg_abs_est := NULL]
top_hits <- sh_hits[seq_len(min(TOP_N, nrow(sh_hits)))]
top_hits[, original_estimate := Estimate]
top_hits[, original_pval := pval]
top_hits[, original_q_within_ct := q_within_ct]
top_hits[, original_q_bonferroni_family := q_bonferroni_family]
top_lr_keys <- top_hits[, .(ct_pair, lr_pair, ligand_complex, receptor_complex,
                            original_estimate, original_pval,
                            original_q_within_ct, original_q_bonferroni_family)]
cat(sprintf("[select] picked %d top LR pairs\n", nrow(top_lr_keys)))

# --- Load metadata + LR scores (with pooled dataset) ----------------------
meta <- fread(META_V2)
meta[, dataset_pooled := ifelse(dataset %in% SMALL_COHORTS,
                                "Other_small",
                                as.character(dataset))]
covar_cols <- intersect(
  c("sample", "dataset", "dataset_pooled", "disease_stage_coarse"),
  names(meta))
meta_use <- meta[, ..covar_cols]
meta_use <- meta_use[!is.na(disease_stage_coarse) & disease_stage_coarse != ""]
stage_levels <- c("Healthy", "Steatosis", "Steatohepatitis", "Cirrhosis")
meta_use[, disease_stage_coarse := factor(disease_stage_coarse, levels = stage_levels)]
meta_use <- meta_use[!is.na(disease_stage_coarse)]
cat(sprintf("[meta] %d donors in scope\n", nrow(meta_use)))
cat("[meta] stage distribution:\n"); print(table(meta_use$disease_stage_coarse))

cat(sprintf("[input] reading %s\n", MERGED_TSV))
lr_long <- fread(MERGED_TSV)
lr_long[, lr_pair := paste(ligand_complex, receptor_complex, sep = "__")]
lr_long[, ct_pair := paste(source, target, sep = "->")]
lr_long <- lr_long[paste(ct_pair, lr_pair, sep = "||") %in%
                   paste(top_lr_keys$ct_pair, top_lr_keys$lr_pair, sep = "||")]
cat(sprintf("[filter] %d rows kept for top-%d LR pairs\n",
            nrow(lr_long), nrow(top_lr_keys)))
lr_long[, score := -log10(pmax(magnitude_rank, 1e-4))]
lr_long <- merge(lr_long, meta_use, by = "sample", all.x = FALSE)
lr_long[, disease_stage_coarse := factor(disease_stage_coarse, levels = stage_levels)]
cat(sprintf("[merge] %d donor-LR rows after metadata join\n", nrow(lr_long)))

lr_long[, lr_key := paste(ct_pair, lr_pair, sep = "||")]
shards <- split(lr_long, by = "lr_key", drop = TRUE)
top_lr_keys[, lr_key := paste(ct_pair, lr_pair, sep = "||")]
shards <- shards[top_lr_keys$lr_key]
shards <- shards[!sapply(shards, is.null)]
cat(sprintf("[shards] %d LR pairs have donor-LR rows\n", length(shards)))

# Donors per stratum
donors_per_stratum <- split(meta_use$sample, meta_use$disease_stage_coarse)

build_donor_index <- function(d) split(seq_len(nrow(d)), d$sample)
shard_donor_idx <- lapply(shards, build_donor_index)

fit_one_replicate <- function(b) {
  set.seed(42L + b)
  drawn <- unlist(lapply(donors_per_stratum, function(v) {
    sample(v, size = length(v), replace = TRUE)
  }), use.names = FALSE)
  unique_drawn <- unique(drawn)
  draw_counts <- tabulate(match(drawn, unique_drawn))
  mult <- setNames(draw_counts, unique_drawn)

  out <- numeric(length(shards))
  for (i in seq_along(shards)) {
    d <- shards[[i]]
    idx_by_donor <- shard_donor_idx[[i]]
    keep_donors <- intersect(names(idx_by_donor), unique_drawn)
    if (!length(keep_donors)) { out[i] <- NA_real_; next }
    row_idx <- unlist(lapply(keep_donors, function(s) {
      rep(idx_by_donor[[s]], times = mult[[s]])
    }), use.names = FALSE)
    dd <- d[row_idx]
    tab <- table(dd$disease_stage_coarse)
    if (any(tab < 5L)) { out[i] <- NA_real_; next }
    if (length(unique(dd$dataset_pooled)) < 2L) { out[i] <- NA_real_; next }
    use_random <- length(unique(dd$dataset_pooled)) >= 3L
    formula_str <- "score ~ disease_stage_coarse + log10(pmax(n_source_cells,1)) + log10(pmax(n_target_cells,1)) + dataset_pooled"
    if (use_random) formula_str <- paste(formula_str, "+ (1|dataset_pooled)")
    fit <- tryCatch(
      if (use_random) {
        lmerTest::lmer(
          stats::as.formula(formula_str),
          data = dd, REML = TRUE,
          control = lmerControl(check.conv.singular = .makeCC("ignore", tol = 1e-4),
                                calc.derivs = FALSE))
      } else {
        stats::lm(stats::as.formula(formula_str), data = dd)
      },
      error = function(e) NULL,
      warning = function(w) {
        if (use_random) {
          suppressWarnings(lmerTest::lmer(
            stats::as.formula(formula_str),
            data = dd, REML = TRUE,
            control = lmerControl(check.conv.singular = .makeCC("ignore", tol = 1e-4),
                                  calc.derivs = FALSE)))
        } else {
          stats::lm(stats::as.formula(formula_str), data = dd)
        }
      })
    if (is.null(fit) || inherits(fit, "try-error")) {
      out[i] <- NA_real_; next
    }
    fe <- if (use_random) tryCatch(lme4::fixef(fit), error = function(e) NULL)
          else tryCatch(stats::coef(fit), error = function(e) NULL)
    if (is.null(fe) || !(TARGET_TERM %in% names(fe))) {
      out[i] <- NA_real_; next
    }
    out[i] <- unname(fe[TARGET_TERM])
  }
  out
}

cat(sprintf("[bootstrap] B=%d over %d shards on %d cores\n",
            B, length(shards), N_CORES))
t0 <- Sys.time()
rep_mat <- parallel::mclapply(seq_len(B), fit_one_replicate, mc.cores = N_CORES)
cat(sprintf("[bootstrap] done in %.1f min\n",
            as.numeric(Sys.time() - t0, units = "mins")))

boot_mat <- do.call(rbind, lapply(rep_mat, function(v) {
  if (length(v) == length(shards)) v else rep(NA_real_, length(shards))
}))
colnames(boot_mat) <- names(shards)
cat(sprintf("[bootstrap] matrix dim: %d x %d\n",
            nrow(boot_mat), ncol(boot_mat)))

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

out_dt <- merge(top_lr_keys, qtable, by = "lr_key", all.x = TRUE)
setcolorder(out_dt,
  c("ct_pair", "lr_pair", "ligand_complex", "receptor_complex",
    "original_estimate", "original_pval",
    "original_q_within_ct", "original_q_bonferroni_family",
    "boot_median", "boot_q025", "boot_q975",
    "ci_excludes_zero", "n_replicates"))
out_dt[, lr_key := NULL]
setorder(out_dt, original_pval)

OUT_FILE <- file.path(V3_DIR, "bootstrap_ci_v3.tsv")
fwrite(out_dt, OUT_FILE, sep = "\t")
cat(sprintf("[output] %d rows -> %s\n", nrow(out_dt), OUT_FILE))

n_exc <- sum(out_dt$ci_excludes_zero, na.rm = TRUE)
cat(sprintf("\n[summary] CI excludes zero: %d / %d (%.1f%%)\n",
            n_exc, nrow(out_dt), 100 * n_exc / nrow(out_dt)))
cat("[done] bootstrap CI pass complete\n")
