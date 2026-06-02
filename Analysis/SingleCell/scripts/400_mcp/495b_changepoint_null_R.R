#!/usr/bin/env Rscript
# 495b_changepoint_null_R.R
# Stage-shuffle null distribution for segmented-regression breakpoints (using R `segmented` package)
# to address reviewer concern that breakpoints at F2 are mechanical.

suppressPackageStartupMessages({
  library(data.table)
  library(segmented)
})

root <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
meta <- fread(file.path(root, "Analysis/SingleCell/results_gpu_v2/mcp/inputs/donor_metadata.tsv"))
scores_f <- file.path(root, "Analysis/SingleCell/results_gpu_v2/mcp/integration/program_donor_scores.tsv.gz")
stopifnot(file.exists(scores_f))
scores <- fread(scores_f)  # sample, program, score

out_dir <- file.path(root, "Analysis/SingleCell/results_gpu_v2/mcp/reviewer_defense")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(42)
N_PERM <- 200
ALL_PROGS <- unique(scores$program)
rows <- list()
for (p in ALL_PROGS) {
  sub <- merge(scores[program == p, .(sample, score)],
               meta[, .(sample, disease_stage_numeric)], by = "sample")
  sub <- sub[!is.na(disease_stage_numeric) & !is.na(score)]
  if (nrow(sub) < 30) next

  # Observed fit
  lm0 <- try(lm(score ~ disease_stage_numeric, data = sub), silent = TRUE)
  seg0 <- try(segmented(lm0, seg.Z = ~ disease_stage_numeric,
                        psi = list(disease_stage_numeric = 1.5)), silent = TRUE)
  obs_bp <- if (!inherits(seg0, "try-error") && !is.null(seg0$psi)) seg0$psi[, "Est."] else NA_real_
  obs_aic_lin <- AIC(lm0)
  obs_aic_seg <- if (!inherits(seg0, "try-error")) AIC(seg0) else NA_real_

  # Null: shuffle stage labels N times
  null_bps <- numeric(N_PERM); null_bps[] <- NA_real_
  for (i in seq_len(N_PERM)) {
    sh <- sub[, .(score, disease_stage_numeric = sample(disease_stage_numeric))]
    lm_s <- try(lm(score ~ disease_stage_numeric, data = sh), silent = TRUE)
    seg_s <- try(segmented(lm_s, seg.Z = ~ disease_stage_numeric,
                           psi = list(disease_stage_numeric = 1.5)),
                 silent = TRUE, outFile = "/dev/null")
    if (!inherits(seg_s, "try-error") && !is.null(seg_s$psi)) {
      null_bps[i] <- seg_s$psi[, "Est."]
    }
  }

  null_bps <- null_bps[!is.na(null_bps)]
  n_near_2 <- sum(abs(null_bps - 2.0) <= 0.1)
  rows[[length(rows) + 1]] <- data.table(
    program = p,
    obs_breakpoint = obs_bp,
    obs_near_2 = !is.na(obs_bp) && abs(obs_bp - 2.0) <= 0.1,
    n_null_fits = length(null_bps),
    null_frac_near_2 = if (length(null_bps)) n_near_2 / length(null_bps) else NA_real_,
    null_median_bp = if (length(null_bps)) median(null_bps) else NA_real_,
    obs_aic_delta = obs_aic_lin - obs_aic_seg,
    permutation_p = if (!is.na(obs_bp) && length(null_bps)) {
      # p-value: probability null |bp-obs| >= 0 with tighter CI (here: null near 2.0 vs obs near 2.0)
      if (abs(obs_bp - 2.0) <= 0.1) n_near_2 / length(null_bps) else NA_real_
    } else NA_real_
  )
}
out <- rbindlist(rows, fill = TRUE)
fwrite(out, file.path(out_dir, "changepoint_shuffle_null_R.tsv"), sep = "\t")

cat(sprintf("[495b] Observed breakpoints near stage 2.0 (±0.1): %d / %d programs\n",
            sum(out$obs_near_2, na.rm = TRUE), nrow(out)))
cat(sprintf("[495b] Null mean fraction near 2.0 across %d permutations: %.2f%%\n",
            N_PERM, 100 * mean(out$null_frac_near_2, na.rm = TRUE)))
cat(sprintf("[495b] Interpretation: if null fraction is already high (e.g. ~25%%), observed clustering at F2 is mechanical.\n"))
