#!/usr/bin/env Rscript
# 06: how much of the published concordance is aggregation?
#
# Kamzolas et al. report Pearson R = 0.96-1.0 between histology and trajectory
# position. Their Extended Data Fig. 1 legend says that coefficient is computed
# "between the mean trajectory position of each scoring group" -- so the unit is
# the histology CATEGORY, and each panel has 3 to 8 points. Averaging inside a
# group removes within-group variance from the numerator but not from the
# denominator, so a weak donor-level relationship becomes a near-perfect
# group-level one almost automatically.
#
# This script reports the same statistic three ways: their group-mean version,
# the donor-level version, and a simulation that answers the actual question --
# given the donor-level correlation we observe and the group sizes we have, how
# often does the group-mean estimator alone produce >= 0.96?

suppressPackageStartupMessages({ library(data.table) })
source(file.path(Sys.getenv("MASLD_PROJECT_ROOT"),
                 "scripts/analysis/continuous_axis_benchmark/lib_axis_common.R"))

OUT <- Sys.getenv("CAB_OUT_ROOT", "")
assert_true(nzchar(OUT), "CAB_OUT_ROOT is unset")
dir.create(file.path(OUT, "q1"), recursive = TRUE, showWarnings = FALSE)

pre <- read_prespec()
N_SIM <- pre$aggregation_inflation$n_sim
N_BOOT <- 2000L
set.seed(pre$seeds$sim_base)

meta <- load_manifest()
arms <- c("A1", "A2", "A3", "A4", "A5", "A6")
available <- arms[file.exists(file.path(OUT, "arms", paste0(arms, "_axis.tsv")))]
assert_true(length(available) > 0L, "No arm axis files exist yet")
log_step("arms available: ", paste(available, collapse = ", "))

variables <- list(
  fibrosis_stage = function(d) d$fibrosis_stage,
  nas_score      = function(d) d$nas_score,
  nas_bin        = function(d) fcase(is.na(d$nas_score), NA_integer_,
                                     d$nas_score == 0L, 0L,
                                     d$nas_score <= 2L, 1L,
                                     d$nas_score <= 4L, 2L,
                                     default = 3L)
)

boot_cor <- function(x, y, method) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 20L) return(c(NA_real_, NA_real_, NA_real_))
  obs <- cor(x, y, method = method)
  bs <- replicate(N_BOOT, {
    i <- sample.int(length(x), replace = TRUE)
    if (sd(x[i]) == 0 || sd(y[i]) == 0) NA_real_ else cor(x[i], y[i], method = method)
  })
  c(obs, quantile(bs, 0.025, na.rm = TRUE), quantile(bs, 0.975, na.rm = TRUE))
}

rows <- list(); sims <- list()
for (arm in available) {
  ax <- fread(file.path(OUT, "arms", paste0(arm, "_axis.tsv")))
  d0 <- merge(meta[, .(sample_id, dataset)], ax[, .(sample_id, axis_raw)], by = "sample_id")
  d0 <- merge(d0, meta[, .(sample_id, fibrosis_stage, nas_score)], by = "sample_id")
  for (vn in names(variables)) {
    d <- copy(d0); d[, level := variables[[vn]](d)]
    d <- d[is.finite(level) & is.finite(axis_raw)]
    if (!nrow(d) || uniqueN(d$level) < 3L) next

    grp <- d[, .(mean_axis = mean(axis_raw), n = .N), by = level][order(level)]
    r_group <- cor(grp$level, grp$mean_axis)                       # their estimator
    rp <- boot_cor(d$axis_raw, d$level, "pearson")                 # donor level
    rs <- boot_cor(d$axis_raw, d$level, "spearman")

    # Variance decomposition and ICC of the axis across histology levels.
    gm <- d[, .(m = mean(axis_raw), n = .N, v = var(axis_raw)), by = level]
    between <- sum(gm$n * (gm$m - mean(d$axis_raw))^2) / (nrow(gm) - 1)
    within <- sum((gm$n - 1) * ifelse(is.na(gm$v), 0, gm$v)) / (nrow(d) - nrow(gm))
    nbar <- mean(gm$n)
    icc <- (between - within) / (between + (nbar - 1) * within)

    # The simulation. Hold the donor-level correlation and the group sizes fixed,
    # generate data with NO structure beyond that correlation, and see what the
    # group-mean estimator reports.
    r_true <- rp[1]
    if (is.finite(r_true)) {
      lev <- d$level
      sim_r <- replicate(N_SIM, {
        z <- r_true * scale(lev)[, 1] + sqrt(max(0, 1 - r_true^2)) * rnorm(length(lev))
        g <- tapply(z, lev, mean)
        suppressWarnings(cor(as.numeric(names(g)), as.numeric(g)))
      })
      sims[[paste(arm, vn)]] <- data.table(
        arm = arm, variable = vn, r_patient_used = r_true,
        sim_median_r_group = median(sim_r, na.rm = TRUE),
        sim_p_r_group_ge_0.96 = mean(sim_r >= 0.96, na.rm = TRUE),
        sim_p_r_group_ge_r_observed = mean(sim_r >= r_group, na.rm = TRUE),
        n_sim = N_SIM)
    }

    rows[[paste(arm, vn)]] <- data.table(
      arm = arm, variable = vn, n_donors = nrow(d), n_levels = nrow(grp),
      min_level_n = min(grp$n), max_level_n = max(grp$n),
      r_group_kamzolas_estimator = r_group,
      r_patient_pearson = rp[1], r_patient_pearson_lo = rp[2], r_patient_pearson_hi = rp[3],
      r_patient_spearman = rs[1], r_patient_spearman_lo = rs[2], r_patient_spearman_hi = rs[3],
      inflation_group_minus_patient = r_group - rp[1],
      icc_axis_across_levels = icc,
      between_level_var = between, within_level_var = within)
  }
}

res <- rbindlist(rows); simres <- rbindlist(sims)
setorder(res, arm, variable)
write_tsv_once(res, file.path(OUT, "q1", "aggregation_inflation.tsv"))
write_tsv_once(simres, file.path(OUT, "q1", "aggregation_inflation_simulation.tsv"))
print(res[, .(arm, variable, n_levels, r_group_kamzolas_estimator,
              r_patient_pearson, inflation_group_minus_patient)])
print(simres[, .(arm, variable, r_patient_used, sim_median_r_group, sim_p_r_group_ge_0.96)])
log_step("Q1_AGGREGATION_COMPLETE")
