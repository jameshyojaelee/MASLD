#!/usr/bin/env Rscript
# Arm A5: a supervised severity score, presumed dead until it proves otherwise.
#
# This repo already killed a staging classifier twice over. It leaked (feature
# selection and imputation happened inside the CV loop, so a resubstitution AUROC
# of 0.87 was reported as cross-validated against a true leave-one-cohort-out
# ~0.70), and more damningly a valid null beat it: 500 RANDOM genes gave AUROC
# 0.860 against the curated set's 0.853, p = 0.733.
#
# So this arm is built to be disqualified. Everything happens inside the fold:
# residualisation for sex, standardisation, lambda selection. And it runs its own
# null of 200 draws of 500 random genes through the identical pipeline. Per the
# prespecification, A5 may not be nominated the winning axis unless its held-out
# C-index exceeds the 95th percentile of that null.

suppressPackageStartupMessages({
  library(data.table); library(glmnet)
})
source(file.path(Sys.getenv("MASLD_PROJECT_ROOT"),
                 "scripts/analysis/continuous_axis_benchmark/lib_axis_common.R"))

OUT <- Sys.getenv("CAB_OUT_ROOT", "")
assert_true(nzchar(OUT), "CAB_OUT_ROOT is unset")
MATRIX_ID <- Sys.getenv("CAB_MATRIX", "E_qn")   # no ComBat: cohort is the fold

pre <- read_prespec()
N_NULL <- pre$a_priori_predictions$P_a5_presumed_dead$n_null_draws
NULL_SIZE <- pre$a_priori_predictions$P_a5_presumed_dead$null_gene_set_size
set.seed(pre$seeds$master)

meta <- load_manifest()
E <- readRDS(file.path(OUT, "arms", paste0(MATRIX_ID, ".rds")))
popC <- population(meta, "POP-C")
info <- popC[!is.na(fibrosis_stage)]
info[, dataset := droplevels(factor(dataset))]
folds <- levels(info$dataset)
log_step("arm A5: n=", nrow(info), " over ", length(folds), " leave-one-cohort-out folds")

# Somers' D_xy rescaled to a C-index. Rank based, so it tolerates the heavy ties
# in a 5-level ordinal outcome.
c_index <- function(pred, truth) {
  ok <- is.finite(pred) & is.finite(truth)
  if (sum(ok) < 10L) return(NA_real_)
  p <- pred[ok]; t <- truth[ok]
  cc <- 0; dc <- 0
  for (i in seq_along(t)) {
    dt <- t[-i] - t[i]; dp <- p[-i] - p[i]
    cmp <- dt != 0
    cc <- cc + sum(sign(dt[cmp]) == sign(dp[cmp]))
    dc <- dc + sum(sign(dt[cmp]) == -sign(dp[cmp]))
  }
  if ((cc + dc) == 0) return(NA_real_)
  0.5 * (((cc - dc) / (cc + dc)) + 1)
}

# One leave-one-cohort-out pass. Every transform is fitted on training only.
run_loco <- function(genes, seed) {
  set.seed(seed)
  preds <- rep(NA_real_, nrow(info))
  for (f in folds) {
    tr <- which(info$dataset != f); te <- which(info$dataset == f)
    if (length(tr) < 30L || length(te) < 10L) next
    Xtr <- t(E[genes, info$sample_id[tr], drop = FALSE])
    Xte <- t(E[genes, info$sample_id[te], drop = FALSE])
    # Residualise for sex using training coefficients only.
    sx_tr <- model.matrix(~ inferred_sex, data = info[tr])
    sx_te <- model.matrix(~ inferred_sex, data = info[te])
    B <- qr.solve(sx_tr, Xtr)
    Xtr <- Xtr - sx_tr %*% B
    Xte <- Xte - sx_te %*% B
    ctr <- colMeans(Xtr); sdv <- apply(Xtr, 2, sd); sdv[sdv < 1e-8] <- 1
    Xtr <- sweep(sweep(Xtr, 2, ctr), 2, sdv, "/")
    Xte <- sweep(sweep(Xte, 2, ctr), 2, sdv, "/")
    cv <- cv.glmnet(Xtr, info$fibrosis_stage[tr], alpha = 0.5,
                    family = "gaussian", nfolds = 5)
    preds[te] <- as.numeric(predict(cv, newx = Xte, s = "lambda.min"))
  }
  preds
}

log_step("observed model on all ", nrow(E), " genes")
obs_pred <- run_loco(rownames(E), pre$seeds$master)
obs_c <- c_index(obs_pred, info$fibrosis_stage)
per_fold <- rbindlist(lapply(folds, function(f) {
  i <- which(info$dataset == f)
  data.table(fold = f, n = length(i),
             c_index = c_index(obs_pred[i], info$fibrosis_stage[i]))
}))
log_step("observed held-out C-index = ", round(obs_c, 4))

log_step("mandatory null: ", N_NULL, " draws of ", NULL_SIZE, " random genes")
null_c <- numeric(N_NULL)
for (b in seq_len(N_NULL)) {
  set.seed(pre$seeds$sim_base + b)
  g <- sample(rownames(E), NULL_SIZE)
  null_c[b] <- c_index(run_loco(g, pre$seeds$sim_base + b), info$fibrosis_stage)
  if (b %% 20 == 0) log_step("  null ", b, "/", N_NULL,
                             "  running q95=", round(quantile(null_c[1:b], .95, na.rm = TRUE), 4))
}

q95 <- quantile(null_c, 0.95, na.rm = TRUE)
emp_p <- (1 + sum(null_c >= obs_c, na.rm = TRUE)) / (1 + sum(is.finite(null_c)))
passes <- is.finite(obs_c) && obs_c > q95

axis_dt <- data.table(sample_id = info$sample_id, dataset = info$dataset,
                      fibrosis_stage = info$fibrosis_stage, nas_score = info$nas_score,
                      axis_raw = obs_pred)
axis_dt[, axis_rank_within_cohort := rank_within(axis_raw, dataset)]
axis_dt[, `:=`(arm = "A5", matrix_id = MATRIX_ID)]
write_tsv_once(axis_dt, file.path(OUT, "arms", "A5_axis.tsv"))
write_tsv_once(data.table(draw = seq_len(N_NULL), c_index = null_c),
               file.path(OUT, "arms", "A5_random_gene_null.tsv"))
write_tsv_once(per_fold, file.path(OUT, "arms", "A5_per_fold.tsv"))

summary_dt <- data.table(
  arm = "A5", matrix_id = MATRIX_ID, n_donors = nrow(info), n_folds = length(folds),
  observed_c_index = obs_c,
  null_median_c_index = median(null_c, na.rm = TRUE),
  null_q95_c_index = q95,
  empirical_p_vs_random_genes = emp_p,
  disqualified = !passes,
  prior_precedent = "a prior staging classifier gave 0.860 for 500 random genes vs 0.853 curated, p=0.733",
  note = if (passes) "A5 beat its own random-gene null and may be considered"
         else "A5 did not beat 500 random genes; reported as a negative control, not an axis"
)
write_tsv_once(summary_dt, file.path(OUT, "arms", "A5_axis_summary.tsv"))
print(summary_dt)
log_step("ARM_A5_COMPLETE disqualified=", !passes)
