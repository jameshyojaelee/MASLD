# T6 library: effect modification of program-axis associations by sex, age and
# cohort.
#
# THE ESTIMAND
# For a frozen program p and a recorded histologic axis a, the discovery release
# reports a single meta-analytic slope beta_pa: how far the program score moves,
# in program-score SD, per one unit of the recorded axis. This system asks
# whether that one number is hiding two: whether the slope differs by sex, by
# age, or between cohorts.
#
# WHY THE MODEL IS FITTED WITHIN COHORT AND ONLY THEN META-ANALYSED
# The four discovery cohorts differ in sex ratio (F fraction 0.42 to 0.67), in
# stage composition, and in every technical respect. A pooled model with a
# sex-by-stage product term cannot distinguish an interaction from the fact that
# the cohorts with more women also have different stage distributions and
# different library preparation. Fitting the product term inside each cohort
# removes cohort composition from the estimate by construction, and the
# meta-analysis then asks whether the within-cohort interactions agree.
#
# WHY THE COHORT ARM IS A HETEROGENEITY QUESTION AND NOT AN INTERACTION FIT
# "Does cohort modify the association" is exactly the question Cochran's Q asks
# of the per-cohort slopes that already exist. Refitting a cohort-by-stage
# product term on pooled data would answer a different and worse-identified
# question, because cohort is confounded with everything. The heterogeneity
# statistics come from the same REML + Knapp-Hartung fits the discovery release
# used.
#
# WHAT A NULL HERE IS AND IS NOT
# An interaction test has far less power than the main effect it modifies, so a
# null result is uninterpretable without the smallest interaction the design
# could have resolved. Every negative in this system therefore carries a minimum
# detectable interaction, and the cohort arm carries a minimum detectable
# between-cohort standard deviation obtained by simulation from each program's
# own observed standard errors.
#
# LANGUAGE
# Fibrosis stage and NAS are cross-sectional recorded histologic variables. A
# slope here is a cross-sectional association, and a modified slope is a
# difference between two cross-sectional associations. Nothing in this file
# implies ordering or movement.

suppressPackageStartupMessages({
  library(data.table)
})

# --------------------------------------------------------------------- design

# Explicit numeric design. Built by hand rather than through model.matrix() so
# the interaction column names are stable and so the same matrix construction is
# reused unchanged inside the permutation loop, where a factor level that
# vanishes under permutation would otherwise silently rename a coefficient.
#
# The modifier is centered within cohort. Centering leaves the interaction
# coefficient exactly unchanged and makes the axis main effects readable as the
# slope at the cohort's average modifier value rather than at modifier zero.
t6_design_matrix <- function(md, modifier_column, include_sex_covariate) {
  fib <- md$fibrosis_stage - mean(md$fibrosis_stage)
  nas <- md$nas_score - mean(md$nas_score)
  mod <- md[[modifier_column]]
  mod <- mod - mean(mod)
  columns <- list(`(Intercept)` = rep(1, nrow(md)))
  if (include_sex_covariate) {
    sex <- md$sex_male - mean(md$sex_male)
    columns[["sex_male"]] <- sex
  }
  columns[["modifier"]] <- mod
  columns[["fibrosis"]] <- fib
  columns[["nas"]] <- nas
  columns[["modifier:fibrosis"]] <- mod * fib
  columns[["modifier:nas"]] <- mod * nas
  X <- do.call(cbind, columns)
  colnames(X) <- names(columns)
  X
}

T6_INTERACTION_TERMS <- c(fibrosis = "modifier:fibrosis", nas = "modifier:nas")

# Per-cohort interaction coefficients for every program at once. Within a cohort
# the design is identical across programs and the score matrix is complete, so
# hc3_fit_matrix() returns all of them from four matrix products. run_tests.R in
# this directory already proves hc3_fit_matrix agrees with the per-feature
# hc3_fit to 1e-10, and t6_run_tests.R re-checks it on this design.
#
# Returns the column set meta_analyze() consumes, so the meta path is literally
# the discovery path and not a reimplementation of it.
t6_fit_cohort_interactions <- function(scores, meta, cohorts, modifier_column,
                                       include_sex_covariate, model_kind) {
  rows <- list()
  k <- 0L
  for (cohort in cohorts) {
    md <- meta[dataset == cohort]
    md <- md[is.finite(md[[modifier_column]])]
    if (!nrow(md)) next
    sample_index <- match(md$sample_id, colnames(scores))
    Y <- scores[, sample_index, drop = FALSE]
    usable <- rownames(Y)[rowSums(is.finite(Y)) == ncol(Y)]
    X <- t6_design_matrix(md, modifier_column, include_sex_covariate)
    fit <- hc3_fit_matrix(Y[usable, , drop = FALSE], X)
    for (axis in names(T6_INTERACTION_TERMS)) {
      term <- T6_INTERACTION_TERMS[[axis]]
      estimable <- isTRUE(fit$estimable) && term %in% rownames(fit$coefficients)
      beta <- if (estimable) fit$coefficients[term, ] else
        setNames(rep(NA_real_, length(usable)), usable)
      se <- if (estimable) fit$se[term, ] else
        setNames(rep(NA_real_, length(usable)), usable)
      k <- k + 1L
      rows[[k]] <- data.table(
        feature_id = usable,
        cohort = cohort,
        model_kind = model_kind,
        axis = axis,
        term = term,
        estimable = estimable,
        failure_reason = if (estimable) NA_character_ else
          if (!is.null(fit$failure_reason)) fit$failure_reason else "term_not_estimable",
        beta = as.numeric(beta[usable]),
        se_hc3 = as.numeric(se[usable]),
        n = nrow(X),
        residual_df = if (!is.null(fit$df)) fit$df else NA_integer_,
        n_model_columns = ncol(X)
      )
    }
    dropped <- setdiff(rownames(scores), usable)
    if (length(dropped)) {
      message("      ", cohort, ": ", length(dropped),
              " program(s) had a non-finite score and were not fitted")
    }
  }
  if (!length(rows)) fail("No cohort produced an interaction fit")
  out <- rbindlist(rows, use.names = TRUE, fill = TRUE)
  out[, critical := stats::qt(0.975, df = residual_df)]
  out[, `:=`(ci_lower = beta - critical * se_hc3, ci_upper = beta + critical * se_hc3,
             statistic_hc3 = beta / se_hc3)]
  out[, p_value := 2 * stats::pt(-abs(statistic_hc3), df = residual_df)]
  out[, critical := NULL]
  out[]
}

# ---------------------------------------------------------------- permutation

# Null for an interaction test.
#
# The naive null permutes the modifier within cohort. That destroys three things
# at once: the modifier's association with expression, its interaction with the
# axis, and its association with the axis itself. Only the second is the null
# hypothesis. Breaking the third changes the collinearity of the design, so the
# permuted design can be better conditioned than the real one and the measured
# false-call rate comes out optimistic.
#
# The stratified null permutes the modifier within cohort AND within level of the
# axis being tested. The modifier's joint distribution with that axis is then
# preserved exactly, so the permuted design has the same collinearity as the
# real one, and the only thing removed is the interaction. This is the null the
# gate is taken from. The unstratified version is computed as a diagnostic so
# the difference between the two remains auditable rather than asserted.
#
# Donors sitting alone in a stratum cannot move. permutable_fraction records how
# many can, because a stratification fine enough to immobilise the sample would
# produce a flattering null by doing nothing.
t6_permute_modifier <- function(meta, modifier_column, cohort_column = "dataset",
                                stratum_column = NULL) {
  out <- copy(meta)
  by_columns <- if (is.null(stratum_column)) cohort_column else
    c(cohort_column, stratum_column)
  out[, .perm := sample.int(.N), by = by_columns]
  out[, (modifier_column) := get(modifier_column)[.perm], by = by_columns]
  out[, .perm := NULL]
  out[]
}

t6_permutable_fraction <- function(meta, cohort_column = "dataset",
                                   stratum_column = NULL) {
  by_columns <- if (is.null(stratum_column)) cohort_column else
    c(cohort_column, stratum_column)
  sizes <- meta[, .N, by = by_columns]
  sum(sizes[N > 1L, N]) / sum(sizes$N)
}

# ------------------------------------------------------- minimum detectable

# Minimum detectable between-cohort standard deviation of a slope, for the
# Cochran Q test the cohort arm reports.
#
# There is no clean closed form: Q under heterogeneity is a weighted sum of
# non-central terms whose weights depend on the individual standard errors. So
# it is simulated from each program's own observed per-cohort standard errors,
# which is the only thing that makes the answer program-specific rather than a
# single design-level number that would hide the programs with the weakest
# cohorts.
#
# Returns the smallest tau on the grid reaching the requested power. NA means no
# tau on the grid reached it, which is itself the answer for that program.
t6_minimum_detectable_tau <- function(se_vector, tau_grid, n_sims = 2000L,
                                      power = 0.80, alpha = 0.05) {
  se_vector <- se_vector[is.finite(se_vector) & se_vector > 0]
  k <- length(se_vector)
  if (k < 3L) return(NA_real_)
  critical <- stats::qchisq(1 - alpha, df = k - 1L)
  w <- 1 / se_vector^2
  for (tau in tau_grid) {
    sd_total <- sqrt(se_vector^2 + tau^2)
    B <- matrix(stats::rnorm(n_sims * k, mean = 0, sd = rep(sd_total, each = n_sims)),
                nrow = n_sims, ncol = k)
    weighted_mean <- as.numeric(B %*% w) / sum(w)
    Q <- rowSums(sweep(B, 1L, weighted_mean, "-")^2 * rep(w, each = n_sims))
    if (mean(Q > critical) >= power) return(tau)
  }
  NA_real_
}

# ------------------------------------------------------------ support states

# Support state for a modifier arm. Deliberately the same three-way split the
# discovery release uses, so a null here means the same thing a null there does:
# either the design excluded an interaction of the stated size, or it could not
# have seen one.
t6_support_state <- function(q_value, ci_lower, ci_upper, estimable, threshold) {
  fcase(
    !estimable, "untestable",
    is.finite(q_value) & q_value < 0.05, "modified",
    informative_null(ci_lower, ci_upper, threshold), "informative_null",
    default = "indeterminate"
  )
}

# Uniformity of the observed p-value family, as a description only. The KS
# p-value here references a theoretical uniform, which is the wrong reference for
# a family of 113 correlated program scores: the programs share genes and share
# donors, so their p-values are dependent and the theoretical KS null is not
# theirs. It is reported because it is the conventional summary, and it is
# reported next to t6_family_signal(), which does the same job against the right
# null and is the one to read.
t6_uniformity <- function(p_values, label) {
  p <- p_values[is.finite(p_values)]
  if (length(p) < 5L) {
    return(data.table(view = label, n = length(p), ks_statistic = NA_real_,
                      ks_p_theoretical = NA_real_, median_p = NA_real_,
                      fraction_below_0.05 = NA_real_))
  }
  ks <- suppressWarnings(stats::ks.test(p, "punif"))
  data.table(view = label, n = length(p), ks_statistic = as.numeric(ks$statistic),
             ks_p_theoretical = as.numeric(ks$p.value), median_p = stats::median(p),
             fraction_below_0.05 = mean(p < 0.05))
}

# Family-level test for ANY modification in a view, against the view's own
# permutation null.
#
# WHY THIS EXISTS
# BH across 113 programs asks whether an individual program is modified. That is
# a demanding question and it is not the only one worth answering: a set of
# genuinely modified programs, each too small to clear multiplicity, would leave
# BH at zero while shifting the whole family's p-value distribution. Reporting
# "zero programs" without checking that would overstate the null.
#
# Two statistics, both referred to the permutation null already computed for the
# calibration gate, so no new assumption enters:
#   excess_small_p : the fraction of the family below 0.05. Sensitive to many
#                    weak modifications spread across programs.
#   min_p          : the single strongest program. Sensitive to one strong
#                    modification that BH still refused.
# Both are empirical, so their resolution floor is 1/(B+1) and they are only
# interpretable when B is large.
t6_family_signal <- function(observed_p, null_matrix, label) {
  finite_observed <- observed_p[is.finite(observed_p)]
  if (!length(finite_observed) || !is.matrix(null_matrix) || !nrow(null_matrix)) {
    fail("Family-level signal test needs observed p-values and a permutation null")
  }
  observed_fraction <- mean(finite_observed < 0.05)
  observed_min <- min(finite_observed)
  null_fraction <- apply(null_matrix, 1L, function(p) mean(p[is.finite(p)] < 0.05))
  null_min <- apply(null_matrix, 1L, function(p) {
    p <- p[is.finite(p)]
    if (!length(p)) NA_real_ else min(p)
  })
  null_fraction <- null_fraction[is.finite(null_fraction)]
  null_min <- null_min[is.finite(null_min)]
  data.table(
    view = label,
    n_permutations = nrow(null_matrix),
    observed_fraction_below_0.05 = observed_fraction,
    null_median_fraction_below_0.05 = stats::median(null_fraction),
    p_excess_small = (1 + sum(null_fraction >= observed_fraction)) / (length(null_fraction) + 1),
    observed_min_p = observed_min,
    null_median_min_p = stats::median(null_min),
    p_strongest_program = (1 + sum(null_min <= observed_min)) / (length(null_min) + 1),
    resolution_floor = 1 / (nrow(null_matrix) + 1)
  )
}
