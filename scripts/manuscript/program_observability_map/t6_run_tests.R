#!/usr/bin/env Rscript

# T6 unit tests, on synthetic data with a KNOWN answer.
#
# The point of this file is that a null result is only worth reporting if the
# machinery that produced it can be shown to detect a real effect first. So the
# central tests plant an interaction of a stated size and require the pipeline
# to find it, then plant none and require it to find none. A third test plants
# no interaction but makes the cohorts differ in sex ratio, stage distribution
# and slope simultaneously, and requires the pooled model to report a false
# interaction while the within-cohort meta path does not -- which is the reason
# the real analysis is fitted within cohort.

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", args[grepl("^--file=", args)])))
source(file.path(script_dir, "analysis_lib.R"))
source(file.path(script_dir, "calibration_lib.R"))
source(file.path(script_dir, "systems_contract.R"))
source(file.path(script_dir, "t6_modifier_lib.R"))

failures <- 0L
check <- function(label, expr) {
  ok <- isTRUE(tryCatch(expr, error = function(e) {
    message("    error: ", conditionMessage(e)); FALSE
  }))
  message(if (ok) "  PASS  " else "  FAIL  ", label)
  if (!ok) failures <<- failures + 1L
  invisible(ok)
}

set.seed(20260812L)

COHORT_N <- c(GSE130970 = 76L, GSE135251 = 214L, GSE162694 = 86L, GSE174478 = 93L)
N_PROGRAMS <- 113L
N_PLANTED <- 10L

# Synthetic substrate shaped like the real discovery design: four cohorts, one
# donor per row, a five-level fibrosis axis, a nine-level NAS axis, and a binary
# modifier. `delta` is the planted sex-by-fibrosis interaction, in program-score
# SD per fibrosis stage per unit of the modifier.
make_synthetic <- function(delta, n_programs = N_PROGRAMS, n_planted = N_PLANTED,
                           female_fraction = rep(0.5, 4L), stage_shift = rep(0, 4L),
                           cohort_slope = rep(0.3, 4L)) {
  meta_rows <- list()
  for (i in seq_along(COHORT_N)) {
    n <- COHORT_N[[i]]
    fib <- pmin(4L, pmax(0L, round(stats::rnorm(n, mean = 2 + stage_shift[[i]], sd = 1.2))))
    meta_rows[[i]] <- data.table(
      sample_id = paste0(names(COHORT_N)[[i]], "_", seq_len(n)),
      dataset = names(COHORT_N)[[i]],
      fibrosis_stage = as.numeric(fib),
      nas_score = as.numeric(pmin(8L, pmax(0L, round(stats::rnorm(n, mean = 4, sd = 1.8))))),
      sex_male = as.numeric(stats::rbinom(n, 1L, 1 - female_fraction[[i]])),
      age = stats::rnorm(n, mean = 50, sd = 13),
      cohort_slope = cohort_slope[[i]]
    )
  }
  meta <- rbindlist(meta_rows)
  programs <- sprintf("program_%03d", seq_len(n_programs))
  scores <- matrix(stats::rnorm(n_programs * nrow(meta)), nrow = n_programs,
                   dimnames = list(programs, meta$sample_id))
  fib_c <- meta[, fibrosis_stage - mean(fibrosis_stage), by = dataset]$V1
  sex_c <- meta[, sex_male - mean(sex_male), by = dataset]$V1
  for (j in seq_len(n_programs)) {
    scores[j, ] <- scores[j, ] + meta$cohort_slope * fib_c
    if (j <= n_planted) scores[j, ] <- scores[j, ] + delta * sex_c * fib_c
  }
  # Rescale within cohort exactly as the real projection does, so the planted
  # effect is expressed in the same program-score SD units the real one is.
  scores <- standardize_scores_by_cohort(scores, meta$dataset)
  list(meta = meta, scores = scores, programs = programs,
       planted = programs[seq_len(n_planted)])
}

# The arm under test, end to end: within-cohort HC3 interaction fits, REML +
# Knapp-Hartung meta-analysis, BH across the program family.
run_arm <- function(synthetic, modifier_column = "sex_male") {
  cohort_effects <- t6_fit_cohort_interactions(
    synthetic$scores, synthetic$meta, names(COHORT_N), modifier_column,
    include_sex_covariate = FALSE, model_kind = "test_arm")
  meta_effects <- meta_analyze(cohort_effects, min_cohorts = 3L)
  meta_effects[estimable == TRUE, q_value := p.adjust(p_value_meta, method = "BH"),
               by = axis]
  meta_effects[]
}

message("\n--- matrix HC3 agrees with the per-feature reference on the T6 design ---")
syn <- make_synthetic(delta = 0.5)
check("hc3_fit_matrix reproduces hc3_fit interaction beta and SE to 1e-10", {
  md <- syn$meta[dataset == "GSE130970"]
  X <- t6_design_matrix(md, "sex_male", include_sex_covariate = FALSE)
  Y <- syn$scores[, md$sample_id, drop = FALSE]
  matrix_fit <- hc3_fit_matrix(Y, X)
  d <- copy(md)
  d[, `:=`(modifier = X[, "modifier"], fibrosis = X[, "fibrosis"], nas = X[, "nas"])]
  worst <- 0
  for (feature in rownames(Y)[1:12]) {
    d[, score := as.numeric(Y[feature, ])]
    ref <- hc3_fit(d$score, d, score ~ modifier + fibrosis + nas +
                     modifier:fibrosis + modifier:nas)
    worst <- max(worst,
                 abs(ref$coefficients[["modifier:fibrosis"]] -
                       matrix_fit$coefficients["modifier:fibrosis", feature]),
                 abs(ref$se[["modifier:fibrosis"]] -
                       matrix_fit$se["modifier:fibrosis", feature]))
  }
  message("    worst absolute discrepancy: ", format(worst, digits = 3))
  worst < 1e-10
})

check("centering the modifier leaves the interaction coefficient unchanged", {
  md <- syn$meta[dataset == "GSE162694"]
  X <- t6_design_matrix(md, "sex_male", include_sex_covariate = FALSE)
  shifted <- copy(md)[, sex_male := sex_male + 7]
  X_shifted <- t6_design_matrix(shifted, "sex_male", include_sex_covariate = FALSE)
  Y <- syn$scores[, md$sample_id, drop = FALSE]
  a <- hc3_fit_matrix(Y, X)$coefficients["modifier:fibrosis", ]
  b <- hc3_fit_matrix(Y, X_shifted)$coefficients["modifier:fibrosis", ]
  max(abs(a - b)) < 1e-10
})

message("\n--- PLANTED INTERACTION: the test must find an effect that is there ---")

# A single planted size proves less than a curve. What matters for the power
# statement is where BH detection turns on, so the whole curve is measured and
# printed: nothing at zero, nothing worth claiming at a modification the size of
# the main effects themselves, and reliable recovery only once the modification
# is several times that. That shape IS the finding this system reports.
recovery <- rbindlist(lapply(c(0, 0.5, 1.0, 1.5, 2.0), function(delta) {
  syn_d <- make_synthetic(delta)
  result <- run_arm(syn_d)[axis == "fibrosis"]
  hits <- result[q_value < 0.05, feature_id]
  data.table(
    planted_delta = delta,
    median_recovered_beta = result[feature_id %in% syn_d$planted, median(abs(beta_meta))],
    median_se_meta = result[, median(se_meta, na.rm = TRUE)],
    n_planted_called = length(intersect(hits, syn_d$planted)),
    n_unplanted_called = length(setdiff(hits, syn_d$planted)))
}))
message("    BH recovery curve on 10 planted of 113 programs:")
print(recovery)

syn_big <- make_synthetic(delta = 2.0)
big_result <- run_arm(syn_big)
check("a planted 2.00 sex-by-fibrosis interaction is recovered in >= 8 of 10 programs", {
  hits <- big_result[axis == "fibrosis" & q_value < 0.05, feature_id]
  true_hits <- intersect(hits, syn_big$planted)
  message("    planted 2.00: ", length(true_hits), "/", N_PLANTED,
          " planted programs called, ", length(setdiff(hits, syn_big$planted)),
          " unplanted programs called")
  length(true_hits) >= 8L && length(setdiff(hits, syn_big$planted)) <= 2L
})

check("BH recovery is nothing at zero and rises with the planted size", {
  recovery[planted_delta == 0, n_planted_called] == 0L &&
    !is.unsorted(recovery$n_planted_called) &&
    recovery[planted_delta == 2.0, n_planted_called] >= 8L
})

check("the planted interaction does NOT leak onto the unplanted NAS axis", {
  nas_hits <- big_result[axis == "nas" & q_value < 0.05, .N]
  message("    NAS-axis calls with no NAS interaction planted: ", nas_hits)
  nas_hits <= 2L
})

# The estimate and the multiplicity-corrected test fail at different sizes, and
# confusing the two is exactly how an underpowered null gets written up as an
# absence. At a planted 0.50 the coefficient comes back essentially unbiased
# while BH calls nothing: the estimator is fine, the corrected test is what runs
# out. This is the quantitative basis for the power statement in the report.
check("a planted 0.50 is estimated without bias but is NOT called by BH", {
  syn_mid <- make_synthetic(delta = 0.5)
  result <- run_arm(syn_mid)[axis == "fibrosis"]
  ratio <- result[feature_id %in% syn_mid$planted, median(abs(beta_meta))] / 0.5
  called <- result[q_value < 0.05, .N]
  message("    planted 0.50: recovered/planted = ", format(ratio, digits = 3),
          "; BH calls = ", called, "; median p on planted = ",
          format(result[feature_id %in% syn_mid$planted, median(p_value_meta)], digits = 3))
  ratio > 0.75 && ratio < 1.05 && called == 0L
})

# The two minimum detectable effects this system reports are meant to bracket
# the size BH actually resolves. This checks that they do, rather than assuming
# the alpha = 0.05 figure the discovery release quotes is the operating one.
check("the reported minimum detectable effects bracket the empirical detection point", {
  # The detection point is the smallest planted size at which BH recovers most of
  # the planted programs. It should sit above the alpha = 0.05 minimum detectable
  # effect the discovery release quotes and below the family-wise one, which is
  # the whole reason this system reports both rather than either alone.
  detected <- recovery[n_planted_called >= 8L][1L]
  mde_nominal <- minimum_detectable_effect(detected$median_se_meta, 4L)
  mde_family <- minimum_detectable_effect(detected$median_se_meta, 4L,
                                          alpha = 0.05 / N_PROGRAMS)
  message("    detection point |beta| = ", format(detected$median_recovered_beta, digits = 3),
          " at median se_meta ", format(detected$median_se_meta, digits = 3),
          " | MDE at alpha 0.05 = ", format(mde_nominal, digits = 3),
          " | MDE at alpha 0.05/113 = ", format(mde_family, digits = 3))
  mde_nominal < detected$median_recovered_beta &&
    detected$median_recovered_beta < mde_family
})

message("\n--- NULL: the test must find nothing when nothing is there ---")
check("zero planted interaction gives at most 1 BH call per 113 programs", {
  syn0 <- make_synthetic(delta = 0)
  result <- run_arm(syn0)
  calls <- result[q_value < 0.05, .N, by = axis]
  message("    calls under the complete null: ",
          paste(sprintf("%s=%d", calls$axis, calls$N), collapse = " "),
          if (!nrow(calls)) "none" else "")
  nrow(calls) == 0L || max(calls$N) <= 1L
})

message("\n--- an interaction BELOW the design's resolution must not be claimed ---")
check("a planted 0.15 interaction is mostly not called, and lands in the CI", {
  syn_small <- make_synthetic(delta = 0.15)
  result <- run_arm(syn_small)
  called <- result[axis == "fibrosis" & feature_id %in% syn_small$planted &
                     q_value < 0.05, .N]
  covered <- result[axis == "fibrosis" & feature_id %in% syn_small$planted,
                    mean(ci_lower <= 0.15 & ci_upper >= 0.15)]
  message("    planted 0.15: called ", called, "/", N_PLANTED,
          "; 95% CI covers the truth in ", round(100 * covered), "% of them")
  called <= 3L && covered >= 0.7
})

message("\n--- POOLING MANUFACTURES INTERACTION, which is why we fit within cohort ---")
check("pooled model invents a sex-by-fibrosis interaction the within-cohort meta does not", {
  syn_conf <- make_synthetic(
    delta = 0,
    female_fraction = c(0.75, 0.30, 0.75, 0.30),
    stage_shift = c(-0.9, 0.9, -0.9, 0.9),
    cohort_slope = c(0.10, 0.55, 0.10, 0.55))
  pooled <- copy(syn_conf$meta)
  Xp <- t6_design_matrix(pooled, "sex_male", include_sex_covariate = FALSE)
  pooled_fit <- hc3_fit_matrix(syn_conf$scores[, pooled$sample_id, drop = FALSE], Xp)
  pooled_stat <- pooled_fit$coefficients["modifier:fibrosis", ] /
    pooled_fit$se["modifier:fibrosis", ]
  pooled_p <- 2 * stats::pt(-abs(pooled_stat), df = pooled_fit$df)
  pooled_calls <- sum(p.adjust(pooled_p, method = "BH") < 0.05)
  within <- run_arm(syn_conf)
  within_calls <- within[axis == "fibrosis" & q_value < 0.05, .N]
  message("    zero interaction planted. pooled model calls ", pooled_calls,
          "/113; within-cohort meta calls ", within_calls, "/113")
  pooled_calls > within_calls && within_calls <= 1L
})

message("\n--- stratified permutation preserves the modifier-axis relationship ---")
check("stratified permutation leaves the modifier-by-axis table exactly unchanged", {
  observed <- syn$meta[, table(dataset, fibrosis_stage, sex_male)]
  ok <- TRUE
  for (i in 1:5) {
    permuted <- t6_permute_modifier(syn$meta, "sex_male",
                                    stratum_column = "fibrosis_stage")
    ok <- ok && identical(as.vector(observed),
                          as.vector(permuted[, table(dataset, fibrosis_stage, sex_male)]))
  }
  ok
})

check("unstratified permutation does NOT preserve it, so the two nulls differ", {
  observed <- syn$meta[, table(dataset, fibrosis_stage, sex_male)]
  differs <- FALSE
  for (i in 1:5) {
    permuted <- t6_permute_modifier(syn$meta, "sex_male", stratum_column = NULL)
    differs <- differs || !identical(
      as.vector(observed), as.vector(permuted[, table(dataset, fibrosis_stage, sex_male)]))
  }
  differs
})

check("stratified permutation still moves most donors on this design", {
  fraction <- t6_permutable_fraction(syn$meta, stratum_column = "fibrosis_stage")
  message("    permutable fraction under fibrosis stratification: ",
          format(fraction, digits = 4))
  fraction > 0.95
})

check("the stratified null is calibrated on synthetic data with no interaction", {
  syn0 <- make_synthetic(delta = 0)
  observed <- run_arm(syn0)[axis == "fibrosis"]
  null_p <- lapply(seq_len(30L), function(i) {
    permuted_meta <- t6_permute_modifier(syn0$meta, "sex_male",
                                         stratum_column = "fibrosis_stage")
    ce <- t6_fit_cohort_interactions(syn0$scores, permuted_meta, names(COHORT_N),
                                     "sex_male", FALSE, "null_arm")
    me <- meta_analyze(ce, min_cohorts = 3L)
    me[axis == "fibrosis" & estimable == TRUE, p_value_meta]
  })
  calibration <- calibrate_view(observed$p_value_meta, function(i) null_p[[i]],
                               n_reps = 30L, n_tests = N_PROGRAMS,
                               label = "synthetic_sex_x_fibrosis")
  message("    mean null BH calls per 113: ", format(calibration$null_call_mean, digits = 3),
          " (max ", calibration$null_call_max, ")")
  isTRUE(calibration$reportable)
})

message("\n--- minimum detectable between-cohort SD behaves ---")
check("MDE tau falls when the per-cohort standard errors fall", {
  grid <- seq(0.02, 1.00, by = 0.02)
  set.seed(1L); tight <- t6_minimum_detectable_tau(rep(0.03, 4L), grid, n_sims = 4000L)
  set.seed(1L); loose <- t6_minimum_detectable_tau(rep(0.30, 4L), grid, n_sims = 4000L)
  message("    MDE tau at SE 0.03: ", tight, "; at SE 0.30: ", loose)
  is.finite(tight) && is.finite(loose) && tight < loose
})

check("MDE tau delivers the power it claims, verified by an independent simulation", {
  grid <- seq(0.02, 1.00, by = 0.02)
  se <- c(0.08, 0.05, 0.09, 0.07)
  set.seed(2L)
  tau <- t6_minimum_detectable_tau(se, grid, n_sims = 8000L)
  set.seed(99L)
  n <- 20000L
  B <- matrix(stats::rnorm(n * 4L, sd = rep(sqrt(se^2 + tau^2), each = n)), nrow = n)
  w <- 1 / se^2
  bw <- as.numeric(B %*% w) / sum(w)
  Q <- rowSums(sweep(B, 1L, bw, "-")^2 * rep(w, each = n))
  achieved <- mean(Q > stats::qchisq(0.95, df = 3L))
  message("    claimed tau ", tau, " achieves power ", format(achieved, digits = 3))
  is.finite(tau) && achieved >= 0.78 && achieved <= 0.92
})

check("Q test holds its nominal size when tau is truly zero", {
  se <- c(0.08, 0.05, 0.09, 0.07)
  set.seed(3L)
  n <- 40000L
  B <- matrix(stats::rnorm(n * 4L, sd = rep(se, each = n)), nrow = n)
  w <- 1 / se^2
  bw <- as.numeric(B %*% w) / sum(w)
  Q <- rowSums(sweep(B, 1L, bw, "-")^2 * rep(w, each = n))
  size <- mean(Q > stats::qchisq(0.95, df = 3L))
  message("    empirical size of the Q test at tau = 0: ", format(size, digits = 3))
  abs(size - 0.05) < 0.005
})

message("\n--- support states and contract assertions ---")
check("support state splits a null into informative and indeterminate", {
  states <- t6_support_state(
    q_value = c(0.001, 0.6, 0.6, NA_real_),
    ci_lower = c(0.30, -0.05, -0.90, NA_real_),
    ci_upper = c(0.70, 0.05, 0.90, NA_real_),
    estimable = c(TRUE, TRUE, TRUE, FALSE),
    threshold = 0.20)
  identical(states, c("modified", "informative_null", "indeterminate", "untestable"))
})

check("assert_negatives_have_mde rejects a negative with no minimum detectable effect", {
  bad <- data.table(support = "informative_null", mde = NA_real_)
  inherits(try(assert_negatives_have_mde(bad, "support", "mde"), silent = TRUE), "try-error")
})

check("assert_counts_calibrated refuses a count from an uncalibrated view", {
  rows <- data.table(view = "sex_x_fibrosis", reportable = FALSE)
  inherits(try(assert_counts_calibrated(1L, rows), silent = TRUE), "try-error")
})

check("the report vocabulary passes the language gate", {
  isTRUE(assert_language(c(
    "cross-sectional association between a frozen program score and recorded fibrosis stage",
    "modified", "informative_null", "indeterminate", "untestable",
    "sex_x_fibrosis", "age_x_nas", "cohort_heterogeneity_fibrosis")))
})

message("")
if (failures) {
  fail(failures, " T6 unit test(s) failed")
}
cat("T6_TESTS_PASS\n")
