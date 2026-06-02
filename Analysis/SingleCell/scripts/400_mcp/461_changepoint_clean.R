#!/usr/bin/env Rscript
# 461_changepoint_clean.R — SENSITIVITY VARIANT of 461_changepoint.R
#
# Repeats the segmented + GAM changepoint screen but DROPS donors flagged
# `exclude_stage_analysis == TRUE` in donor_metadata_extended.tsv (protocol-
# contaminated GSE136103 NPC-enriched donors + dubious-Healthy donors).
#
# Output: program_changepoint_clean.tsv  (same schema as program_changepoint.tsv)
# Plus mirror at results_gpu_v2/mcp/crn_changepoints_clean.csv for top-level access.

suppressPackageStartupMessages({
  library(data.table)
  library(segmented)
  library(mgcv)
})

root      <- Sys.getenv("MASLD_PROJECT_ROOT",
                        "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
meta_path <- file.path(root, "Analysis/SingleCell/results_gpu_v2/mcp/inputs/donor_metadata.tsv")
ext_path  <- file.path(root, "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/donor_metadata_extended.tsv")
scores_f  <- file.path(root, "Analysis/SingleCell/results_gpu_v2/mcp/integration/program_donor_scores.tsv.gz")
out_path  <- file.path(root, "Analysis/SingleCell/results_gpu_v2/mcp/integration/program_changepoint_clean.tsv")
mirror_p  <- file.path(root, "Analysis/SingleCell/results_gpu_v2/mcp/crn_changepoints_clean.csv")

stopifnot(file.exists(scores_f), file.exists(meta_path), file.exists(ext_path))

scores <- fread(scores_f)    # long: sample, program, score
meta   <- fread(meta_path)
ext    <- fread(ext_path)

# Merge in exclude flag
ext_sub <- ext[, .(sample, exclude_stage_analysis, protocol_group)]
meta <- merge(meta, ext_sub, by = "sample", all.x = TRUE)
meta[is.na(exclude_stage_analysis), exclude_stage_analysis := FALSE]

n_pre  <- nrow(meta)
n_excl <- sum(meta$exclude_stage_analysis)
meta_clean <- meta[exclude_stage_analysis == FALSE]
cat(sprintf("[461c] CLEAN: excluded %d/%d donors (exclude_stage_analysis==TRUE)\n", n_excl, n_pre))
cat(sprintf("[461c] remaining donors: %d\n", nrow(meta_clean)))

stage <- meta_clean[, .(sample, disease_stage_numeric)]

# Verify stage distribution post-exclusion
cat("[461c] disease_stage_numeric post-exclusion:\n")
print(table(stage$disease_stage_numeric, useNA = "ifany"))

rows <- list()
programs <- unique(scores$program)
cat(sprintf("[461c] programs to test: %d\n", length(programs)))

for (p in programs) {
  sub <- merge(scores[program == p, .(sample, score)], stage, by = "sample")
  sub <- sub[!is.na(disease_stage_numeric)]
  if (nrow(sub) < 15) next

  # Linear
  lm0 <- try(lm(score ~ disease_stage_numeric, data = sub), silent = TRUE)
  if (inherits(lm0, "try-error")) next

  # GAM — degrade k if too few unique stage values (k must be <= unique_x)
  n_unique <- length(unique(sub$disease_stage_numeric))
  k_use <- min(4, max(3, n_unique))
  gam0 <- try(gam(score ~ s(disease_stage_numeric, k = k_use), data = sub, method = "REML"),
              silent = TRUE)

  # Segmented (1 breakpoint, seed in middle)
  # NOTE: after exclusion, stage 3 disappears; seed at midpoint of available range
  seed_psi <- mean(range(sub$disease_stage_numeric, na.rm = TRUE))
  seg1 <- try(segmented(lm0, seg.Z = ~ disease_stage_numeric,
                        psi = list(disease_stage_numeric = seed_psi)),
              silent = TRUE)
  if (inherits(seg1, "try-error")) seg1 <- NULL

  aics <- c(
    linear = AIC(lm0),
    gam    = if (!inherits(gam0, "try-error")) AIC(gam0) else NA,
    seg1   = if (!is.null(seg1)) AIC(seg1) else NA
  )
  best <- names(which.min(aics))
  delta_aic <- aics[best] - aics["linear"]

  bp_stage <- NA_real_
  bp_ci_lo <- NA_real_; bp_ci_hi <- NA_real_
  if (!is.null(seg1) && best == "seg1") {
    bp_stage <- tryCatch(seg1$psi[, "Est."], error = function(e) NA_real_)
    ci <- tryCatch(confint.segmented(seg1, level = 0.95)[[1]], error = function(e) NULL)
    if (!is.null(ci) && is.matrix(ci) && nrow(ci) >= 1) {
      lo_col <- grep("low", colnames(ci), value = TRUE, ignore.case = TRUE)[1]
      hi_col <- grep("up|high", colnames(ci), value = TRUE, ignore.case = TRUE)[1]
      if (!is.na(lo_col)) bp_ci_lo <- ci[1, lo_col]
      if (!is.na(hi_col)) bp_ci_hi <- ci[1, hi_col]
    }
  }

  rows[[length(rows) + 1]] <- data.table(
    program            = p,
    best_model         = best,
    aic_linear         = aics["linear"],
    aic_gam            = aics["gam"],
    aic_seg1           = aics["seg1"],
    delta_aic          = delta_aic,
    breakpoint_stage   = bp_stage,
    breakpoint_ci_low  = bp_ci_lo,
    breakpoint_ci_high = bp_ci_hi,
    n                  = nrow(sub)
  )
}

out <- rbindlist(rows, fill = TRUE)
fwrite(out, out_path, sep = "\t")
fwrite(out, mirror_p, sep = ",")
cat(sprintf("[461c] wrote %d rows; best_model distribution:\n", nrow(out)))
print(table(out$best_model, useNA = "ifany"))

bp_hits <- out[best_model == "seg1" & !is.na(breakpoint_stage)]
if (nrow(bp_hits) > 0) {
  cat("[461c] Breakpoint distribution (stage scale; 0=Healthy, 1=Steatosis, 2=Steatohep, 3=Cirrhosis):\n")
  print(summary(bp_hits$breakpoint_stage))
}
cat("[461c] DONE.\n")
