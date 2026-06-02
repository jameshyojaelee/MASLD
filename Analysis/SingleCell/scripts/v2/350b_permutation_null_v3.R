#!/usr/bin/env Rscript
# ============================================================================
# 350b_permutation_null_v3.R
#
# Addresses critique #3: bulk-concordance and LMM hits reported without a
# permutation null. With ~7,000 LMM tests per axis and donor-level
# heterogeneity, BH-FDR could be optimistic if the model is mis-specified.
#
# Design: B=500 donor-stage-label permutations. For each replicate:
#   1. Shuffle the (sample -> disease_stage_coarse) mapping while preserving
#      stratum sizes (i.e. random reassignment among donors).
#   2. Re-fit the SAME coarse-axis LMM as 07b_chain_v3_hardened, retaining
#      the SH-vs-Healthy term estimate / p-value.
# Then: empirical p-value per LR pair = fraction of permutations with
# |estimate_perm| >= |estimate_observed|. Empirical FDR at α = fraction of
# permutations whose smallest p-value <= α, vs BH-FDR-expected α.
#
# Output: stage_trajectory_v3/permutation_fdr_v3.tsv
#   Columns: ct_pair, lr_pair, original_estimate, original_pval,
#            n_perms_more_extreme, empirical_pval, empirical_q (BH)
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(lme4)
  library(lmerTest)
  library(parallel)
})

N_CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "16"))
options(mc.cores = N_CORES)
B <- as.integer(Sys.getenv("PERMUTATION_B", "500"))
TOP_N <- 100L  # bound permutation expense to top-N LR pairs
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
cat(sprintf("[input] top-%d LR pairs selected\n", nrow(top_hits)))

meta <- fread(META_V2)
meta[, dataset_pooled := ifelse(dataset %in% SMALL_COHORTS, "Other_small",
                                as.character(dataset))]
meta_use <- meta[, .(sample, dataset, dataset_pooled, disease_stage_coarse)]
meta_use <- meta_use[!is.na(disease_stage_coarse) & disease_stage_coarse != ""]
stage_levels <- c("Healthy", "Steatosis", "Steatohepatitis", "Cirrhosis")
meta_use[, disease_stage_coarse := factor(disease_stage_coarse, levels = stage_levels)]
meta_use <- meta_use[!is.na(disease_stage_coarse)]

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

donor_stage <- unique(meta_use[, .(sample, disease_stage_coarse)])

fit_perm_replicate <- function(b) {
  set.seed(123L + b)
  # Permute stage labels across donors (preserves stratum sizes)
  perm_stage <- sample(donor_stage$disease_stage_coarse,
                       size = nrow(donor_stage), replace = FALSE)
  perm_map <- data.table(sample = donor_stage$sample,
                         disease_stage_perm = perm_stage)
  out <- numeric(length(shards))
  for (i in seq_along(shards)) {
    d <- copy(shards[[i]])
    d <- merge(d[, !c("disease_stage_coarse"), with = FALSE],
               perm_map, by = "sample")
    setnames(d, "disease_stage_perm", "disease_stage_coarse")
    d[, disease_stage_coarse := factor(disease_stage_coarse, levels = stage_levels)]
    tab <- table(d$disease_stage_coarse)
    if (any(tab < 5L)) { out[i] <- NA_real_; next }
    if (length(unique(d$dataset_pooled)) < 2L) { out[i] <- NA_real_; next }
    use_random <- length(unique(d$dataset_pooled)) >= 3L
    formula_str <- "score ~ disease_stage_coarse + log10(pmax(n_source_cells,1)) + log10(pmax(n_target_cells,1)) + dataset_pooled"
    if (use_random) formula_str <- paste(formula_str, "+ (1|dataset_pooled)")
    fit <- tryCatch(
      if (use_random) {
        lmerTest::lmer(stats::as.formula(formula_str), data = d, REML = TRUE,
          control = lmerControl(check.conv.singular = .makeCC("ignore", tol = 1e-4),
                                calc.derivs = FALSE))
      } else {
        stats::lm(stats::as.formula(formula_str), data = d)
      },
      error = function(e) NULL,
      warning = function(w) NULL)
    if (is.null(fit)) { out[i] <- NA_real_; next }
    fe <- if (use_random) tryCatch(lme4::fixef(fit), error = function(e) NULL)
          else tryCatch(stats::coef(fit), error = function(e) NULL)
    if (is.null(fe) || !(TARGET_TERM %in% names(fe))) {
      out[i] <- NA_real_; next
    }
    out[i] <- unname(fe[TARGET_TERM])
  }
  out
}

cat(sprintf("[perm] B=%d on %d cores\n", B, N_CORES))
t0 <- Sys.time()
perm_list <- parallel::mclapply(seq_len(B), fit_perm_replicate, mc.cores = N_CORES)
cat(sprintf("[perm] done in %.1f min\n",
            as.numeric(Sys.time() - t0, units = "mins")))

perm_mat <- do.call(rbind, lapply(perm_list, function(v) {
  if (length(v) == length(shards)) v else rep(NA_real_, length(shards))
}))
colnames(perm_mat) <- names(shards)

# Empirical p-value: fraction of permutations with |perm_est| >= |original_est|
obs_est <- top_hits$Estimate
names(obs_est) <- top_hits$lr_key
obs_est <- obs_est[colnames(perm_mat)]

n_more_extreme <- sapply(seq_along(obs_est), function(j) {
  v <- perm_mat[, j]
  v <- v[!is.na(v)]
  if (!length(v)) return(NA_integer_)
  sum(abs(v) >= abs(obs_est[j]))
})
n_valid <- apply(perm_mat, 2, function(v) sum(!is.na(v)))
emp_p <- ifelse(n_valid > 0,
                (n_more_extreme + 1) / (n_valid + 1),
                NA_real_)

out_dt <- data.table(
  lr_key = colnames(perm_mat),
  original_estimate = obs_est,
  n_perms_valid = n_valid,
  n_perms_more_extreme = n_more_extreme,
  empirical_pval = emp_p
)
out_dt[, empirical_q := stats::p.adjust(empirical_pval, method = "BH")]
out_dt <- merge(top_hits[, .(lr_key, ct_pair, lr_pair, ligand_complex,
                              receptor_complex, original_pval = pval,
                              original_q_within_ct = q_within_ct,
                              original_q_bonferroni_family = q_bonferroni_family)],
                out_dt, by = "lr_key", all.x = TRUE)
out_dt[, lr_key := NULL]
setorder(out_dt, empirical_pval, original_pval)
OUT <- file.path(V3_DIR, "permutation_fdr_v3.tsv")
fwrite(out_dt, OUT, sep = "\t")
cat(sprintf("[output] %d rows -> %s\n", nrow(out_dt), OUT))
cat(sprintf("[summary] emp_p<0.05: %d ; emp_p<0.10: %d ; emp_q<0.10: %d\n",
            sum(out_dt$empirical_pval < 0.05, na.rm = TRUE),
            sum(out_dt$empirical_pval < 0.10, na.rm = TRUE),
            sum(out_dt$empirical_q < 0.10, na.rm = TRUE)))
cat("[done] permutation null pass complete\n")
