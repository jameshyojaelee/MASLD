# Shared calibration harness for every system in the program framework.
#
# WHY THIS IS SHARED AND NOT REIMPLEMENTED PER SYSTEM
# Three separate results in this workstream were wrong in ways only a
# permutation null revealed, and each failure was a different mechanism:
#
#   - a joint meta test referenced chi-square(2) while its univariate arm
#     referenced Knapp-Hartung t on 3 df, giving 13.0 false calls per 113;
#   - an HC3 Wald on a sparse categorical block gave 37.9 false calls per 113 on
#     the NAS axis while behaving correctly (0.9) on fibrosis;
#   - a donor-level statistic reported zero significant donors that turned out to
#     be arithmetically impossible to reject, because the permutation resolution
#     floor sat above the BH threshold the design needed.
#
# The first two are calibration failures and the third is a resolution failure.
# They are different problems and both have to be checked, so both are checked
# here rather than left to each team to remember.

suppressPackageStartupMessages({
  library(data.table)
})

# Expected number of BH rejections under the complete null is well below one, so
# a view averaging more than one false call per family is not usable for a count.
DEFAULT_CALIBRATION_TOLERANCE <- 1.0

# TWO DISTINCT CHECKS, AND THEY APPLY TO DIFFERENT VIEWS.
#
#  CALIBRATION asks whether the view invents signal under the null. It applies to
#  every view without exception, and it is measured by permutation regardless of
#  where the observed p-values came from.
#
#  RESOLUTION asks whether the observed p-value can be small enough to survive BH
#  at all. It applies ONLY when the observed p is itself empirical, because an
#  empirical p cannot fall below 1/(B+1). A parametric p from a t, F or
#  chi-square reference has no such floor, so the check is meaningless there.
#
# Conflating these was an error in the first version of this file, caught by its
# own test: a perfectly calibrated view with parametric p-values was marked
# unreportable because 300 permutations could not resolve a 113-test family it
# never needed to resolve.
#
# observed_p       : p-values from the real data
# permuted_p_fn(i) : p-values for permutation i, used to measure the null
# n_tests          : size of the multiplicity family the count is drawn from
# observed_p_is_empirical : TRUE when observed_p came from a permutation
calibrate_view <- function(observed_p, permuted_p_fn, n_reps, n_tests = length(observed_p),
                           alpha = 0.05, label = NA_character_,
                           tolerance = DEFAULT_CALIBRATION_TOLERANCE,
                           observed_p_is_empirical = FALSE) {
  stopifnot(is.function(permuted_p_fn), n_reps >= 1L)
  observed_calls <- sum(p.adjust(observed_p, method = "BH") < alpha, na.rm = TRUE)
  null_calls <- vapply(seq_len(n_reps), function(i) {
    p <- permuted_p_fn(i)
    sum(p.adjust(p, method = "BH") < alpha, na.rm = TRUE)
  }, numeric(1))

  resolution_floor <- if (observed_p_is_empirical) 1 / (n_reps + 1) else 0
  bh_threshold_needed <- alpha / n_tests
  calibrated <- mean(null_calls) <= tolerance
  resolvable <- resolution_floor < bh_threshold_needed

  list(
    label = label,
    observed_calls = observed_calls,
    null_call_mean = mean(null_calls),
    null_call_max = max(null_calls),
    null_call_sd = stats::sd(null_calls),
    n_reps = n_reps,
    n_tests = n_tests,
    tolerance = tolerance,
    observed_p_is_empirical = observed_p_is_empirical,
    calibrated = calibrated,
    resolution_floor = resolution_floor,
    bh_threshold_needed = bh_threshold_needed,
    resolvable = resolvable,
    reportable = calibrated && resolvable,
    null_calls = null_calls
  )
}

# The gate. A system calls this before writing any count. Refusing here is the
# intended behaviour, not an error condition to work around.
gate_view <- function(calibration, strict = TRUE) {
  reasons <- c(
    if (!calibration$calibrated)
      sprintf("uncalibrated: %.2f mean false calls per %d tests (tolerance %.2f)",
              calibration$null_call_mean, calibration$n_tests, calibration$tolerance),
    if (!calibration$resolvable)
      sprintf("unresolvable: permutation floor %.3g exceeds the BH threshold %.3g; run at least %d permutations",
              calibration$resolution_floor, calibration$bh_threshold_needed,
              ceiling(calibration$n_tests / 0.05))
  )
  if (length(reasons) && strict) {
    stop(sprintf("CALIBRATION GATE FAILED for '%s': %s",
                 calibration$label, paste(reasons, collapse = "; ")), call. = FALSE)
  }
  invisible(list(passed = !length(reasons), reasons = reasons))
}

# One row per view, written next to every count a system reports. A count
# without its row is not reportable.
calibration_row <- function(calibration) {
  data.table(
    view = calibration$label,
    observed_calls = calibration$observed_calls,
    null_call_mean = calibration$null_call_mean,
    null_call_max = calibration$null_call_max,
    n_permutations = calibration$n_reps,
    n_tests = calibration$n_tests,
    resolution_floor = calibration$resolution_floor,
    bh_threshold_needed = calibration$bh_threshold_needed,
    calibrated = calibration$calibrated,
    resolvable = calibration$resolvable,
    reportable = calibration$reportable
  )
}

# Minimum permutations for a family of n_tests to be resolvable at alpha.
required_permutations <- function(n_tests, alpha = 0.05) ceiling(n_tests / alpha)

# Joint within-cohort permutation of the histologic variables. Permuting the
# columns together preserves their mutual correlation (fibrosis and NAS are
# correlated at r = 0.474 pooled), so the null removes the association with
# expression without also destroying the relationship between the axes.
permute_histology_within_cohort <- function(meta, columns, cohort_column = "dataset") {
  out <- copy(meta)
  out[, .perm := sample.int(.N), by = c(cohort_column)]
  for (column in columns) {
    out[, (column) := get(column)[.perm], by = c(cohort_column)]
  }
  out[, .perm := NULL]
  out[]
}
