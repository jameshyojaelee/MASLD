#!/usr/bin/env Rscript

# Synthetic checks for the T1 machinery. Everything here runs on data whose
# truth is known by construction, so a failure means the code is wrong rather
# than the liver being complicated.

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_dir <- dirname(normalizePath(script_arg))
source(file.path(script_dir, "config.R"))
source(file.path(script_dir, "analysis_lib.R"))
source(file.path(script_dir, "calibration_lib.R"))
source(file.path(script_dir, "t1_cellstate_lib.R"))

passed <- 0L
check <- function(condition, label) {
  if (!isTRUE(condition)) stop("FAIL: ", label, call. = FALSE)
  passed <<- passed + 1L
  message("  ok  ", label)
}

set.seed(20260812L)

message("[1] t1_term_effects reproduces the per-feature reference fit")
n <- 30L
stage <- rep(0:3, length.out = n)
Y <- rbind(f1 = 0.4 * stage + rnorm(n), f2 = rnorm(n), f3 = -0.6 * stage + rnorm(n, sd = 2))
colnames(Y) <- paste0("d", seq_len(n))
X <- cbind(`(Intercept)` = 1, F_stage_documented = stage)
fast <- t1_term_effects(Y, X, "F_stage_documented")
reference <- rbindlist(lapply(rownames(Y), function(f) {
  fit <- hc3_fit(Y[f, ], data.frame(y = Y[f, ], F_stage_documented = stage),
                 y ~ F_stage_documented)
  data.table(feature_id = f, beta = fit$coefficients[["F_stage_documented"]],
             se = fit$se[["F_stage_documented"]], p_value = fit$p_value[["F_stage_documented"]])
}))
merged <- merge(fast, reference, by = "feature_id", suffixes = c("_fast", "_ref"))
check(max(abs(merged$beta_fast - merged$beta_ref)) < 1e-10, "coefficients agree")
check(max(abs(merged$se_fast - merged$se_ref)) < 1e-10, "HC3 standard errors agree")
check(max(abs(merged$p_value_fast - merged$p_value_ref)) < 1e-10, "p-values agree")

message("[2] weighted fit reduces to the unweighted fit at equal weights")
w <- rep(sqrt(7), n)
weighted <- t1_term_effects(sweep(Y, 2L, w, "*"), X * w, "F_stage_documented")
check(max(abs(weighted$beta - fast$beta)) < 1e-10, "equal weights leave the coefficient unchanged")
check(max(abs(weighted$se - fast$se)) < 1e-10, "equal weights leave the SE unchanged")

message("[3] weighted fit recovers truth better when precision varies")
# Half the donors carry a hundred times fewer cells, so their scores are noisy.
truth <- 0.5
cells <- rep(c(1000, 10), length.out = n)
y_noisy <- truth * stage + rnorm(n, sd = 1 / sqrt(cells / 1000))
Yn <- matrix(y_noisy, nrow = 1L, dimnames = list("f", colnames(Y)))
plain <- t1_term_effects(Yn, X, "F_stage_documented")
prec <- t1_term_effects(sweep(Yn, 2L, sqrt(cells), "*"), X * sqrt(cells), "F_stage_documented")
check(prec$se < plain$se, "precision weighting tightens the SE when cell counts differ")

message("[4] CLR drops a structurally absent lineage and centres the rest")
counts <- cbind(a = c(100, 200, 50, 80), b = c(10, 5, 20, 15),
                c = c(0, 0, 0, 3), d = c(40, 60, 30, 25))
rownames(counts) <- paste0("d", 1:4)
res <- clr_cell_counts(counts)
check(identical(res$dropped_lineages, "c"), "the lineage zero in 3 of 4 donors is dropped")
check(all(abs(rowSums(res$clr)) < 1e-12), "CLR rows sum to zero")
dense <- clr_cell_counts(counts[, c("a", "b", "d")])
check(length(dense$dropped_lineages) == 0L, "no dense lineage is dropped")

message("[5] classification recovers activity, abundance and indeterminate by construction")
mk <- function(feature_id, cell_type, q_act, beta_act, se_act, q_abu, beta_abu, se_abu, supported) {
  data.table(feature_id = feature_id, cell_type = cell_type,
             q_activity = q_act, beta_activity = beta_act,
             ci_lower_activity = beta_act - 1.96 * se_act,
             ci_upper_activity = beta_act + 1.96 * se_act,
             mde_activity = 2.8 * se_act,
             q_abundance = q_abu, beta_abundance = beta_abu,
             ci_lower_abundance = beta_abu - 1.96 * se_abu,
             ci_upper_abundance = beta_abu + 1.96 * se_abu,
             mde_abundance = 2.8 * se_abu,
             beta_bulk_fibrosis = 0.3, bulk_supported = supported)
}
synthetic <- rbindlist(list(
  # activity present, abundance absent
  mk("p_activity", "hepatocytes", 0.001, 0.40, 0.05, 0.90, 0.01, 0.02, TRUE),
  # activity ruled out tightly, abundance present
  mk("p_abundance", "fibroblasts", 0.80, 0.01, 0.03, 0.001, 0.50, 0.06, TRUE),
  # both present
  mk("p_both", "macrophages", 0.001, 0.45, 0.05, 0.01, 0.40, 0.08, TRUE),
  # activity null but the interval still admits a 0.20 effect: underpowered
  mk("p_indeterminate", "tcells", 0.60, 0.05, 0.30, 0.70, 0.02, 0.30, TRUE),
  # both arms null and the activity null is informative
  mk("p_neither", "cholangiocytes", 0.70, 0.01, 0.03, 0.60, 0.01, 0.03, TRUE),
  # activity present but pointing the other way from the bulk effect
  mk("p_discordant", "hepatocytes", 0.001, -0.40, 0.05, 0.90, 0.01, 0.02, TRUE),
  # no bulk association to attribute
  mk("p_nobulk", "hepatocytes", 0.001, 0.40, 0.05, 0.90, 0.01, 0.02, FALSE)))
classified <- t1_classify(synthetic, observability_effect_threshold)
expected <- c(p_activity = "activity", p_abundance = "abundance",
              p_both = "activity_and_abundance",
              p_indeterminate = "indeterminate_underpowered",
              p_neither = "neither_arm_resolved",
              p_discordant = "activity_discordant_sign",
              p_nobulk = "no_bulk_association_to_attribute")
for (f in names(expected)) {
  check(classified[feature_id == f, attribution] == expected[[f]],
        paste0(f, " classifies as ", expected[[f]]))
}
check(!any(is.na(classified$attribution)), "every program receives an attribution")

message("[6] an abundance call is never made on an activity null that is not informative")
wide_null <- mk("p_wide", "fibroblasts", 0.90, 0.02, 0.50, 0.001, 0.50, 0.06, TRUE)
check(t1_classify(wide_null, observability_effect_threshold)$attribution ==
        "indeterminate_underpowered",
      "a wide activity interval blocks the abundance call")

message("[7] the procedure is calibrated on pure noise")
# 117 independent noise programs on 37 donors with a real stage vector. The
# permutation null and the observed data are both null here, so both should
# return essentially no BH calls.
n_donors <- 37L
n_programs <- 117L
stage_real <- rep(0:3, length.out = n_donors)
noise <- matrix(rnorm(n_programs * n_donors), nrow = n_programs,
                dimnames = list(paste0("p", seq_len(n_programs)), paste0("d", seq_len(n_donors))))
design_of <- function(s) cbind(`(Intercept)` = 1, F_stage_documented = s)
observed_p <- t1_term_effects(zscore_rows(noise), design_of(stage_real),
                              "F_stage_documented")$p_value
cal <- calibrate_view(
  observed_p = observed_p,
  permuted_p_fn = function(i) {
    set.seed(1000L + i)
    t1_term_effects(zscore_rows(noise), design_of(sample(stage_real)),
                    "F_stage_documented")$p_value
  },
  n_reps = 300L, n_tests = n_programs, label = "synthetic_noise",
  observed_p_is_empirical = FALSE)
check(cal$null_call_mean <= DEFAULT_CALIBRATION_TOLERANCE,
      sprintf("null false-call mean %.3f is within tolerance", cal$null_call_mean))
check(cal$reportable, "a parametric-p view is not blocked by the permutation resolution floor")
check(gate_view(cal, strict = FALSE)$passed, "the gate passes on calibrated noise")

message("[8] the procedure detects a planted within-cell-type effect")
planted <- noise
planted[1:10, ] <- planted[1:10, ] + 0.6 * matrix(rep(stage_real, each = 10L), nrow = 10L)
planted_p <- t1_term_effects(zscore_rows(planted), design_of(stage_real),
                             "F_stage_documented")$p_value
n_called <- sum(p.adjust(planted_p, method = "BH") < 0.05)
check(n_called >= 8L, sprintf("at least 8 of 10 planted programs are recovered (got %d)", n_called))
check(all(which(p.adjust(planted_p, method = "BH") < 0.05) <= 10L),
      "no unplanted program is called")

message("[9] the gate refuses an uncalibrated view")
bad <- calibrate_view(observed_p = runif(117),
                      permuted_p_fn = function(i) rep(1e-8, 117),
                      n_reps = 20L, n_tests = 117L, label = "deliberately_broken")
check(!bad$reportable, "a view that calls everything under the null is not reportable")
check(!gate_view(bad, strict = FALSE)$passed, "gate_view reports the failure")
check(inherits(try(gate_view(bad, strict = TRUE), silent = TRUE), "try-error"),
      "strict gating stops the run")

message("[10] an empirical-p view is gated on the permutation count backing it")
# This guards a real error: passing calibrate_view the cost of the calibration
# draw instead of the number of permutations the observed empirical p was
# computed against, which understates the resolution floor.
empirical_view <- function(n_reps) {
  calibrate_view(observed_p = runif(21), permuted_p_fn = function(i) runif(21),
                 n_reps = n_reps, n_tests = 21L, label = "empirical_aggregate",
                 observed_p_is_empirical = TRUE)
}
check(!empirical_view(200L)$resolvable, "200 permutations cannot resolve BH over 21 tests")
check(empirical_view(2000L)$resolvable, "2000 permutations can")
check(required_permutations(21L) <= 2000L, "required_permutations agrees")

message("\nT1 TESTS PASSED: ", passed, " checks")
