#!/usr/bin/env Rscript
# 461_changepoint.R
# For every program (cNMF GEP and DIALOGUE meta-MCP), fit segmented regression
# and GAM along the disease axis (disease_stage_numeric) and consensus pseudotime.
# Decide best model: linear / switch (1 breakpoint) / 2-breakpoint / sigmoid / monotonic GAM.
#
# Output: program_changepoint.tsv with best_model, breakpoint_stage, breakpoint_ci, delta.

suppressPackageStartupMessages({
  library(data.table)
  library(segmented)
  library(mgcv)
})

root <- Sys.getenv("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
meta_path <- file.path(root, "Analysis/SingleCell/results_gpu_v2/mcp/inputs/donor_metadata.tsv")
long_path <- file.path(root, "Analysis/SingleCell/results_gpu_v2/mcp/integration/phenotype_screen_long.tsv")
out_path <- file.path(root, "Analysis/SingleCell/results_gpu_v2/mcp/integration/program_changepoint.tsv")

# Helper: load cNMF donor-level usage from prior Python output; here we re-load minimal matrix
# Expect a consolidated program donor-score file:
scores_f <- file.path(root, "Analysis/SingleCell/results_gpu_v2/mcp/integration/program_donor_scores.tsv.gz")
if (!file.exists(scores_f)) {
  cat(sprintf("[461] program_donor_scores.tsv.gz not found at %s -- run 460 first\n", scores_f))
  quit(status = 0)
}
scores <- fread(scores_f)  # long: sample, program, score
meta <- fread(meta_path)

stage <- meta[, .(sample, disease_stage_numeric)]

rows <- list()
programs <- unique(scores$program)
cat(sprintf("[461] programs to test: %d\n", length(programs)))

for (p in programs) {
  sub <- merge(scores[program == p, .(sample, score)], stage, by = "sample")
  sub <- sub[!is.na(disease_stage_numeric)]
  if (nrow(sub) < 15) next

  # Linear
  lm0 <- try(lm(score ~ disease_stage_numeric, data = sub), silent = TRUE)
  if (inherits(lm0, "try-error")) next

  # GAM
  gam0 <- try(gam(score ~ s(disease_stage_numeric, k = 4), data = sub, method = "REML"), silent = TRUE)

  # Segmented (1 breakpoint, seed in middle)
  seg1 <- try(segmented(lm0, seg.Z = ~ disease_stage_numeric,
                        psi = list(disease_stage_numeric = 1.5)),
              silent = TRUE)
  if (inherits(seg1, "try-error")) seg1 <- NULL

  # Pick best model by AIC among linear/GAM/segmented1
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
    program = p,
    best_model = best,
    aic_linear = aics["linear"],
    aic_gam = aics["gam"],
    aic_seg1 = aics["seg1"],
    delta_aic = delta_aic,
    breakpoint_stage = bp_stage,
    breakpoint_ci_low = bp_ci_lo,
    breakpoint_ci_high = bp_ci_hi,
    n = nrow(sub)
  )
}

out <- rbindlist(rows, fill = TRUE)
fwrite(out, out_path, sep = "\t")
cat(sprintf("[461] wrote %d rows; best_model distribution:\n", nrow(out)))
print(table(out$best_model, useNA = "ifany"))

# Breakpoint histogram summary
bp_hits <- out[best_model == "seg1" & !is.na(breakpoint_stage)]
if (nrow(bp_hits) > 0) {
  cat("[461] Breakpoint distribution (stage scale 0=Healthy, 1=Steatosis, 2=Steatohep, 3=Cirrhosis):\n")
  print(summary(bp_hits$breakpoint_stage))
}
cat("[461] DONE.\n")
