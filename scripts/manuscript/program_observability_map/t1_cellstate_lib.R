#!/usr/bin/env Rscript

# T1 shared functions: within-cell-type program ACTIVITY vs cell-type ABUNDANCE.
#
# Everything statistical here delegates to the shared contract libraries. This
# file holds only what is specific to T1: the roster construction, the precision
# gate on cells per donor, the compositional transform, and the classification
# rule. Keeping the rule in a function rather than inline is what lets the
# synthetic tests in t1_tests.R exercise the actual production rule.

suppressPackageStartupMessages({
  library(data.table)
})

# ---------------------------------------------------------------------------
# Cell-type label mapping between the Hotspot scoring cell types and the
# phase05 lineage count columns. Only these five carry Hotspot programs.
# ---------------------------------------------------------------------------
T1_CELL_TYPES <- c("cholangiocytes", "fibroblasts", "hepatocytes",
                   "macrophages", "tcells")

T1_LINEAGE_COLUMN <- c(
  cholangiocytes = "n_Cholangiocytes",
  fibroblasts    = "n_Fibroblasts",
  hepatocytes    = "n_Hepatocytes",
  macrophages    = "n_Macrophages",
  tcells         = "n_T_cells"
)

# ---------------------------------------------------------------------------
# Centred log-ratio on cell counts.
#
# Counts are compositional: the five scored lineages plus the four unscored ones
# sum to the cells that were captured, so a rise in one lineage's proportion is
# not separable from a fall in the others on the raw scale. CLR puts every
# lineage in the same unconstrained space against the geometric mean of the
# whole composition.
#
# Two rules, both fixed before any model was fitted:
#   - a lineage that is zero in more than half the roster is dropped before the
#     geometric mean is taken. In GSE202379 that is B cells, zero in 36 of 37
#     donors; leaving it in makes the reference geometric mean a function of the
#     pseudocount rather than of the data.
#   - remaining zeros take a 0.5-cell pseudocount. A pseudocount below one cell
#     is not a fractional cell, it is a claim that the observation is far more
#     extreme than "none were seen", and it turns a single absent lineage into
#     the largest term in the log-ratio.
clr_cell_counts <- function(count_matrix, zero_lineage_max_fraction = 0.5,
                            pseudocount = 0.5) {
  stopifnot(is.matrix(count_matrix), nrow(count_matrix) > 0L)
  zero_fraction <- colMeans(count_matrix == 0)
  keep <- zero_fraction <= zero_lineage_max_fraction
  if (!any(keep)) fail("No lineage survives the CLR sparsity rule")
  retained <- count_matrix[, keep, drop = FALSE]
  augmented <- retained
  augmented[augmented == 0] <- pseudocount
  logged <- log(augmented)
  clr <- logged - rowMeans(logged)
  list(clr = clr, retained_lineages = colnames(retained),
       dropped_lineages = colnames(count_matrix)[!keep])
}

# ---------------------------------------------------------------------------
# One HC3 linear fit of many features on one shared design, returning the
# coefficient, SE, t and p for a single named term. hc3_fit_matrix() is the
# shared implementation; this only unpacks the term of interest and attaches the
# reference distribution. Parametric p by construction, so calibrate_view() is
# told the observed p is not empirical and the resolution floor does not apply.
# ---------------------------------------------------------------------------
t1_term_effects <- function(Y, design, term) {
  fit <- hc3_fit_matrix(Y, design)
  if (!isTRUE(fit$estimable)) {
    return(data.table(feature_id = rownames(Y), beta = NA_real_, se = NA_real_,
                      statistic = NA_real_, p_value = NA_real_,
                      df = NA_integer_, n = NA_integer_,
                      estimable = FALSE, failure_reason = fit$failure_reason))
  }
  beta <- fit$coefficients[term, ]
  se <- fit$se[term, ]
  statistic <- beta / se
  data.table(
    feature_id = names(beta), beta = as.numeric(beta), se = as.numeric(se),
    statistic = as.numeric(statistic),
    p_value = 2 * stats::pt(-abs(as.numeric(statistic)), df = fit$df),
    df = fit$df, n = fit$n, estimable = TRUE, failure_reason = NA_character_
  )
}

# Standardise each feature across the donors actually modelled, so the
# coefficient is in program-score SD per unit documented F stage and is on the
# same nominal scale as the bulk map's beta_meta_fibrosis. This mirrors the
# bulk arm's standardize_scores_by_cohort(); it is a rescaling only.
t1_standardize <- function(Y) {
  out <- zscore_rows(Y)
  keep <- apply(out, 1L, function(v) all(is.finite(v)))
  list(Y = out[keep, , drop = FALSE], dropped = rownames(out)[!keep])
}

# ---------------------------------------------------------------------------
# Attribution rule.
#
# The question is what a program's BULK fibrosis association is made of. Within
# a cell type the program score cannot move because the liver contains more of
# that cell type, so a within-cell-type effect is activity. The abundance arm is
# the host cell type's own CLR abundance on the same donors, which is what a
# bulk composition term would be picking up.
#
# The asymmetry that had to be designed out: activity and abundance must be
# judged on the SAME donors with the SAME gate. Scoring activity on gated donors
# and abundance on every donor would give the abundance arm a larger, cleaner
# sample and manufacture abundance calls.
#
# A null is only allowed to become evidence when it is informative. An
# unsupported arm whose interval still admits an effect of the frozen
# observability threshold says nothing, and the program goes to indeterminate
# with its minimum detectable effect attached.
t1_classify <- function(map, effect_threshold, alpha = 0.05) {
  d <- copy(map)
  d[, activity_called := is.finite(q_activity) & q_activity < alpha]
  d[, abundance_called := is.finite(q_abundance) & q_abundance < alpha]
  d[, activity_informative_null := !activity_called &
      informative_null(ci_lower_activity, ci_upper_activity, effect_threshold)]
  d[, abundance_informative_null := !abundance_called &
      informative_null(ci_lower_abundance, ci_upper_abundance, effect_threshold)]
  d[, activity_sign_agrees := is.finite(beta_activity) & is.finite(beta_bulk_fibrosis) &
      sign(beta_activity) == sign(beta_bulk_fibrosis)]

  d[, attribution := NA_character_]
  # Programs whose bulk fibrosis association is itself absent are a separate
  # stratum. There is no bulk effect to attribute, so attributing one would be
  # answering a question the data never asked.
  d[bulk_supported == FALSE, attribution := "no_bulk_association_to_attribute"]

  d[bulk_supported == TRUE & activity_called & activity_sign_agrees &
      abundance_called, attribution := "activity_and_abundance"]
  d[bulk_supported == TRUE & activity_called & activity_sign_agrees &
      !abundance_called, attribution := "activity"]
  d[bulk_supported == TRUE & activity_called & !activity_sign_agrees,
    attribution := "activity_discordant_sign"]
  d[bulk_supported == TRUE & !activity_called & activity_informative_null &
      abundance_called, attribution := "abundance"]
  d[bulk_supported == TRUE & !activity_called & activity_informative_null &
      !abundance_called, attribution := "neither_arm_resolved"]
  d[bulk_supported == TRUE & !activity_called & !activity_informative_null,
    attribution := "indeterminate_underpowered"]
  d[]
}

# ---------------------------------------------------------------------------
# Aggregate statistics.
#
# WHY THESE EXIST. Per-program BH over 117 tests on 13 to 37 donors is a blunt
# instrument, and a table of 117 non-rejections is compatible both with "no
# within-cell-type effect exists" and with "this design cannot see one". These
# three statistics pool across programs, so they are far better powered than any
# single program's test, and each one has an exact permutation null from the
# same stream that calibrates the per-program view.
#
#   mean_t_squared     : is there any within-cell-type stage signal at all?
#   spearman_vs_bulk   : do the within-cell-type effects line up with the bulk
#                        fibrosis effects they are supposed to explain?
#   sign_agreement     : the same question with no distributional assumption.
#
# All three are one-sided upward. A negative answer to spearman_vs_bulk is not
# evidence that the bulk effect is abundance; it is an absence of evidence that
# it is activity, and the report has to say it that way.
t1_aggregate_statistics <- function(effects, bulk_beta, strata = NULL) {
  d <- merge(effects[, .(feature_id, cell_type, beta, statistic, p_value)],
             bulk_beta, by = "feature_id")
  # The four bulk-untestable programs have no bulk coefficient to align against.
  # Dropping them here rather than letting cor() return NA keeps the observed
  # statistic and its permutation null defined on the same set of programs.
  d <- d[is.finite(beta) & is.finite(beta_bulk_fibrosis)]
  compute <- function(sub, label) {
    if (nrow(sub) < 4L) return(NULL)
    data.table(
      stratum = label, n_programs = nrow(sub),
      mean_t_squared = mean(sub$statistic^2, na.rm = TRUE),
      spearman_vs_bulk = suppressWarnings(stats::cor(
        sub$beta, sub$beta_bulk_fibrosis, method = "spearman")),
      sign_agreement = mean(sign(sub$beta) == sign(sub$beta_bulk_fibrosis), na.rm = TRUE),
      n_nominal = sum(sub$p_value < 0.05, na.rm = TRUE))
  }
  rows <- list(compute(d, "all_programs"),
               compute(d[bulk_supported == TRUE], "bulk_supported"))
  if (!is.null(strata)) {
    for (s in strata) rows[[length(rows) + 1L]] <- compute(d[cell_type == s], paste0("cell_type_", s))
  }
  rbindlist(Filter(Negate(is.null), rows))
}
