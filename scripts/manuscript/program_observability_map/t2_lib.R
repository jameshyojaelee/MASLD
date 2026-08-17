# T2 (liver proteome) system helpers for the program observability map.
#
# The question this system exists to answer: the bulk RNA substrate records NAS
# only as a composite, so nothing there can say whether a frozen Hotspot program
# tracks fat, ballooning or lobular inflammation. PXD051911 is the only substrate
# in this repository carrying donor-level Kleiner fibrosis, NAS AND all three NAS
# components alongside BMI, age and sex, so it is the only place the composite
# can be separated into its parts.
#
# Everything statistical is delegated to analysis_lib.R / calibration_lib.R. What
# lives here is substrate-specific and would be wrong to generalise.

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
})

# ---------------------------------------------------------------------------
# Substrate constants, fixed before the run.
# ---------------------------------------------------------------------------

T2_HISTOLOGY_COMPONENTS <- c("steatosis", "ballooning", "inflammation")
T2_PERMUTED_COLUMNS <- c("steatosis", "ballooning", "inflammation",
                         "nas", "fibrosis", "is_masld", "ballooning_binary")

# ---------------------------------------------------------------------------
# Acquisition batch.
#
# The 58 liver runs come from two acquisitions that differ in date AND in
# instrument/column string: 20191017_FP_Qexn2_1015_22cmr2um_* (n = 33) and
# 11122020_FP_ATLASLiver_vis1_* (n = 25). This is not a nuisance label that can
# be dropped: per-donor detection depth differs between them at p = 4.9e-08, so
# a batch term is mandatory in every model.
# ---------------------------------------------------------------------------
t2_batch_from_filename <- function(filenames) {
  prefix <- sub("_.*$", "", filenames)
  unknown <- setdiff(unique(prefix), c("20191017", "11122020"))
  if (length(unknown)) {
    fail("Unrecognised acquisition prefix in liver filenames: ",
         paste(unknown, collapse = ", "))
  }
  factor(ifelse(prefix == "20191017", "b2019", "b2020"), levels = c("b2020", "b2019"))
}

# ---------------------------------------------------------------------------
# Donor metadata for the liver arm.
#
# The biological unit is patient_name. assert_biological_unit() is satisfied
# non-vacuously here because patient_name is a distinct column from the run
# filename that carries the sample identity, and the check below confirms the
# one-run-per-donor claim rather than assuming it.
# ---------------------------------------------------------------------------
t2_build_liver_meta <- function(meta_path, sample_columns) {
  meta <- fread(meta_path)
  required <- c("patient_name", "sample_group", "gender", "alder", "bmi",
                "kleiner_fibrosis_grade", "steatosis_score",
                "hepatocellular_ballooning_score", "lobular_inflammation_score",
                "nafld_activity_score", "saf_diagnosis",
                "liver_proteomics_filename")
  missing <- setdiff(required, names(meta))
  if (length(missing)) fail("meta_data.txt lacks columns: ", paste(missing, collapse = ", "))

  m <- meta[liver_proteomics_filename %in% sample_columns]
  assert_true(nrow(m) == length(sample_columns),
              sprintf("Liver quant has %d sample columns but %d metadata rows resolve",
                      length(sample_columns), nrow(m)))
  assert_true(uniqueN(m$liver_proteomics_filename) == nrow(m),
              "A liver run filename appears on more than one metadata row")
  assert_true(uniqueN(m$patient_name) == nrow(m),
              "A donor contributes more than one liver run; the model needs a donor term")
  assert_true(all(m$sample_group == "initial_sample"),
              "A liver run is not from the initial_sample visit")

  m[, sample_id := liver_proteomics_filename]
  m[, dataset := "PXD051911"]
  m[, batch := t2_batch_from_filename(liver_proteomics_filename)]
  m[, sex_final := factor(gender)]
  m[, fibrosis := as.integer(sub("^F", "", kleiner_fibrosis_grade))]
  m[, steatosis := as.numeric(steatosis_score)]
  m[, ballooning := as.numeric(hepatocellular_ballooning_score)]
  m[, inflammation := as.numeric(lobular_inflammation_score)]
  m[, nas := as.numeric(nafld_activity_score)]
  m[, is_masld := as.integer(saf_diagnosis != "No_MASLD")]
  m[, ballooning_binary := as.integer(ballooning >= 1)]
  m[, age_z := as.numeric(scale(as.numeric(alder)))]
  m[, bmi_z := as.numeric(scale(as.numeric(bmi)))]

  axis_columns <- c("steatosis", "ballooning", "inflammation", "nas", "fibrosis",
                    "age_z", "bmi_z", "is_masld")
  for (column in axis_columns) {
    assert_true(all(is.finite(m[[column]])),
                paste0("Donor axis '", column, "' has a missing value; complete-case ",
                       "donor selection must be explicit, not silent"))
  }
  assert_true(all(m$nas == m$steatosis + m$ballooning + m$inflammation),
              "NAS is not the sum of its three recorded components")

  m[, sample_order := match(sample_id, sample_columns)]
  setorder(m, sample_order)
  m[]
}

# ---------------------------------------------------------------------------
# Protein matrix.
#
# The log2 / quantile / median-collapse block is preserved from
# Analysis/Multimodal_Program_Projection/scripts/02_proteomics_projection.R
# (lines 118-135) because the deposited values are RAW LINEAR INTENSITY. Skipping
# log2 makes every model a test on a heavy-tailed scale, and skipping the
# between-array quantile step leaves per-run loading differences in the score.
#
# WHY THE SCORING SET IS COMPLETE-CASE PROTEINS ONLY.
# 15.1% of the matrix is missing, and the missingness is not random with respect
# to the very axes this system tests: per-donor detection depth differs by batch
# (p = 4.9e-08) and correlates with ballooning (rho = -0.296, p = 0.024) and NAS
# (rho = -0.283, p = 0.031). An available-case program score - a weighted mean
# over whichever proteins a donor happened to have detected - would therefore
# differ between donors partly because their detection depth differs, and that
# difference tracks histology. It would manufacture an association that a
# permutation null could not catch, because permutation destroys the histology
# link in the scores but not the depth-histology confounding in the observed
# data. Restricting to proteins observed in all 58 donors costs coverage and
# costs it honestly: every donor's score is a weighted mean over an identical
# protein set, so no between-donor difference can come from detection.
# ---------------------------------------------------------------------------
t2_build_protein_matrix <- function(quant_path, sample_columns,
                                    source_max_na_fraction = 0.50) {
  raw <- fread(quant_path)
  assert_true(all(c("ProteinAccessions", "Genes", "ProteinDescriptions") %in% names(raw)),
              "Liver quant lacks its three identifier columns")
  assert_true(all(sample_columns %in% names(raw)),
              "A requested sample column is absent from the liver quant matrix")

  n_groups_total <- nrow(raw)
  r <- raw[!is.na(Genes) & Genes != "" & !grepl(";", Genes, fixed = TRUE)]
  n_groups_unambiguous <- nrow(r)

  expr <- as.matrix(r[, ..sample_columns])
  storage.mode(expr) <- "double"
  assert_true(max(expr, na.rm = TRUE) > 1000,
              paste("Liver intensities look already log-transformed; the log2 step",
                    "below would be applied twice"))
  expr[!is.finite(expr) | expr <= 0] <- NA_real_
  na_fraction_all <- mean(is.na(expr))
  expr <- log2(expr)
  rownames(expr) <- toupper(trimws(r$Genes))

  keep <- rowMeans(is.na(expr)) <= source_max_na_fraction
  expr <- expr[keep, , drop = FALSE]
  expr <- normalizeBetweenArrays(expr, method = "quantile")

  collapse_one <- function(idx) {
    z <- expr[idx, , drop = FALSE]
    if (nrow(z) == 1L) return(as.numeric(z[1L, ]))
    apply(z, 2L, stats::median, na.rm = TRUE)
  }
  by_symbol <- split(seq_len(nrow(expr)), rownames(expr))
  collapsed <- t(vapply(by_symbol, collapse_one, numeric(ncol(expr))))
  colnames(collapsed) <- colnames(expr)
  collapsed[!is.finite(collapsed)] <- NA_real_

  complete <- rowSums(is.na(collapsed)) == 0L
  list(
    source_universe = collapsed,
    scoring_universe = collapsed[complete, , drop = FALSE],
    per_donor_na_fraction = colMeans(is.na(collapsed)),
    stats = data.table(
      n_protein_groups_deposited = n_groups_total,
      n_protein_groups_unambiguous = n_groups_unambiguous,
      matrix_na_fraction = na_fraction_all,
      source_max_na_fraction = source_max_na_fraction,
      n_symbols_source_universe = nrow(collapsed),
      n_symbols_scoring_universe = sum(complete)
    )
  )
}

# ---------------------------------------------------------------------------
# Frozen-membership adapter.
#
# SCHEMA TRAP. In the retired frozen_program_membership.tsv, `mapped_symbol` is a
# BOOLEAN flag and `gene_symbol` carries the symbol. In the frozen 117-program
# program_membership_v2.tsv, `mapped_symbol` IS the symbol string and there is no
# `gene_symbol` column. Any code carried over from the 22-program pipeline that
# filters `mapped_symbol == TRUE` returns zero genes on v2 and silently produces
# an all-NA score matrix rather than an error. This adapter normalises v2 into
# the (feature_id, gene_symbol, original_l1_weight) shape score_programs() wants
# and refuses either schema it does not recognise.
# ---------------------------------------------------------------------------
t2_membership_adapter <- function(membership) {
  m <- copy(membership)
  if (is.logical(m$mapped_symbol)) {
    assert_true("gene_symbol" %in% names(m),
                "Legacy membership schema has a boolean mapped_symbol but no gene_symbol")
    m <- m[mapped_symbol %in% TRUE]
    symbol <- m$gene_symbol
  } else {
    assert_true(is.character(m$mapped_symbol),
                "mapped_symbol is neither a boolean flag nor a symbol string")
    symbol <- m$mapped_symbol
  }
  assert_true("original_l1_weight" %in% names(m),
              "Membership lacks original_l1_weight")
  out <- data.table(program_uid = m$program_uid,
                    mapped_symbol = toupper(trimws(symbol)),
                    original_l1_weight = as.numeric(m$original_l1_weight))
  out <- out[!is.na(mapped_symbol) & mapped_symbol != "" & is.finite(original_l1_weight)]
  assert_true(nrow(out) > 0L,
              "Membership adapter returned zero genes; the schema trap has fired")
  out[]
}

# Retained L1 weight and observed gene count of each program against a symbol
# universe. This is the observability table for all 117 programs, computed
# independently of whether a program is later scored.
t2_program_coverage <- function(membership, symbols, label) {
  m <- membership[, .(weight = sum(original_l1_weight)),
                  by = .(program_uid, mapped_symbol)]
  total <- m[, .(total_weight = sum(weight), n_genes = .N), by = program_uid]
  observed <- m[mapped_symbol %in% symbols,
                .(retained_weight = sum(weight), n_observed = .N), by = program_uid]
  out <- merge(total, observed, by = "program_uid", all.x = TRUE)
  out[is.na(retained_weight), `:=`(retained_weight = 0, n_observed = 0L)]
  out[, retained_l1 := retained_weight / total_weight]
  out[, universe := label]
  out[]
}

# ---------------------------------------------------------------------------
# The declared coverage rule, applied.
#
# A program score is a weighted mean of the z-scored proteins it can see. Below
# some retained fraction of the program's own L1 weight that mean stops being a
# measurement of the program and becomes a measurement of an arbitrary subset of
# it, and no amount of downstream statistics repairs that. The rule is declared
# in the contract before the run and enforced here.
# ---------------------------------------------------------------------------
t2_apply_coverage_rule <- function(coverage, min_retained_l1, min_observed_genes) {
  out <- copy(coverage)
  out[, testable := is.finite(retained_l1) & retained_l1 >= min_retained_l1 &
        n_observed >= min_observed_genes]
  out[]
}

# ---------------------------------------------------------------------------
# Design matrices. One row per view, all sharing the donor covariate block.
#
# Every view carries batch, sex, age and BMI. Batch because detection depth and
# ballooning both differ across it; BMI because this is a bariatric-enriched
# series where steatosis and adiposity are entangled; age and sex as standing
# donor covariates.
# ---------------------------------------------------------------------------
T2_COVARIATE_TERMS <- c("batch", "sex_final", "age_z", "bmi_z")

t2_design <- function(meta, test_terms, covariates = T2_COVARIATE_TERMS) {
  form <- stats::as.formula(paste("~", paste(c(covariates, test_terms), collapse = " + ")))
  X <- stats::model.matrix(form, data = as.data.frame(meta))
  assert_true(nrow(X) == nrow(meta), "Design matrix dropped rows; a covariate is missing")
  X
}

# Which model coefficient carries each test term. Numeric terms keep their name;
# a factor term would not, and this refuses rather than guessing.
t2_term_coefficients <- function(X, test_terms) {
  hits <- vapply(test_terms, function(term) {
    idx <- which(colnames(X) == term)
    assert_true(length(idx) == 1L,
                paste0("Test term '", term, "' does not map to exactly one coefficient; ",
                       "encode it as numeric or name the contrast explicitly"))
    idx
  }, integer(1))
  setNames(hits, test_terms)
}

# ---------------------------------------------------------------------------
# One view: matrix HC3 over all testable programs against one design.
#
# The p-values are parametric (HC3 Wald against a t reference on n - p df), so
# they carry no permutation resolution floor. calibrate_view() is still called on
# every view without exception, because calibration asks a different question
# from resolution: whether the view invents signal under the null.
# ---------------------------------------------------------------------------
t2_fit_view <- function(scores, meta, test_terms, view_label,
                        covariates = T2_COVARIATE_TERMS) {
  X <- t2_design(meta, test_terms, covariates)
  fit <- hc3_fit_matrix(scores, X)
  if (!isTRUE(fit$estimable)) {
    fail("View '", view_label, "' is not estimable: ", fit$failure_reason)
  }
  idx <- t2_term_coefficients(X, test_terms)
  rows <- rbindlist(lapply(names(idx), function(term) {
    beta <- fit$coefficients[idx[[term]], ]
    se <- fit$se[idx[[term]], ]
    statistic <- beta / se
    axis <- meta[[term]]
    data.table(
      view = view_label,
      term = term,
      term_scale = if (all(axis %in% c(0, 1))) "binary" else "grade",
      axis_sd = stats::sd(axis),
      program_uid = colnames(fit$coefficients),
      beta = as.numeric(beta),
      se_hc3 = as.numeric(se),
      statistic_hc3 = as.numeric(statistic),
      p_value = as.numeric(2 * stats::pt(-abs(statistic), df = fit$df)),
      n_donors = fit$n,
      residual_df = fit$df,
      n_model_columns = ncol(X),
      max_leverage = max(fit$leverage)
    )
  }))
  list(rows = rows, fit = fit, design = X, term_index = idx)
}

# Just the p-vector, for the permutation null. Same code path as the observed
# fit so a calibration failure cannot be an artefact of a second implementation.
t2_view_p <- function(scores, meta, test_terms, covariates = T2_COVARIATE_TERMS) {
  X <- t2_design(meta, test_terms, covariates)
  fit <- hc3_fit_matrix(scores, X)
  if (!isTRUE(fit$estimable)) return(rep(NA_real_, nrow(scores) * length(test_terms)))
  idx <- t2_term_coefficients(X, test_terms)
  unlist(lapply(names(idx), function(term) {
    statistic <- fit$coefficients[idx[[term]], ] / fit$se[idx[[term]], ]
    2 * stats::pt(-abs(statistic), df = fit$df)
  }), use.names = FALSE)
}

# ---------------------------------------------------------------------------
# Minimum detectable effect for a single stratum.
#
# minimum_detectable_effect() in analysis_lib.R is written for the meta-analytic
# case and takes the number of cohorts, from which it forms df = n_cohorts - 1.
# There is one cohort here and the relevant reference is the model's own
# residual df, so this passes residual_df + 1 to recover exactly
# (qt(1 - alpha/2, residual_df) + qnorm(power)) * se. It is the same formula and
# the same function, not a second implementation.
# ---------------------------------------------------------------------------
t2_mde <- function(se, residual_df, power = 0.80, alpha = 0.05) {
  minimum_detectable_effect(se_meta = se, n_cohorts = residual_df + 1L,
                            power = power, alpha = alpha)
}

# ---------------------------------------------------------------------------
# Support classification. Every negative carries its MDE, and the two kinds of
# negative are kept apart: a program whose interval excludes an effect worth
# having is evidence of absence at that scale, while a program whose interval
# still admits one is a statement about this design's power and nothing else.
#
# WHY THERE ARE TWO EFFECT SCALES HERE.
# The estimates are reported in assay-native units, program-score SD per
# one-unit increment of the recorded grade. The observability threshold cannot
# be applied on that scale, because a one-unit increment does not mean the same
# thing on every axis: NAS runs 0-8 while ballooning runs 0-2, so a fixed
# threshold in native units silently demands a four-times-larger effect of NAS
# than of ballooning and any "underpowered" verdict would be partly an artefact
# of axis range. The threshold is therefore evaluated per axis SD, which is the
# same quantity on every ordinal axis. Binary terms are left alone: their native
# coefficient already IS a difference between two groups in score SD, and
# dividing that by the prevalence-dependent SD of a 0/1 variable would make it
# less comparable, not more.
# ---------------------------------------------------------------------------
t2_classify <- function(result, effect_threshold, alpha = 0.05) {
  out <- copy(result)
  critical <- stats::qt(1 - alpha / 2, df = out$residual_df)
  out[, ci_lower := beta - critical * se_hc3]
  out[, ci_upper := beta + critical * se_hc3]
  out[, mde := t2_mde(se_hc3, residual_df)]

  out[, threshold_scale_factor := fifelse(term_scale == "binary", 1, axis_sd)]
  out[, beta_per_axis_sd := beta * threshold_scale_factor]
  out[, se_per_axis_sd := se_hc3 * threshold_scale_factor]
  out[, ci_lower_per_axis_sd := ci_lower * threshold_scale_factor]
  out[, ci_upper_per_axis_sd := ci_upper * threshold_scale_factor]
  out[, mde_per_axis_sd := mde * threshold_scale_factor]

  out[, support := fifelse(
    is.finite(q_value) & q_value < alpha, "supported",
    fifelse(informative_null(ci_lower_per_axis_sd, ci_upper_per_axis_sd,
                             effect_threshold),
            "informative_null", "indeterminate"))]
  out[, effect_threshold := effect_threshold]
  out[, effect_threshold_scale := "program-score SD per axis SD"]
  out[]
}

# ---------------------------------------------------------------------------
# Contrast Wald on the HC3 covariance.
#
# wald_test() in analysis_lib.R tests "these coefficients are zero". An
# equal-weighting hypothesis needs "these linear combinations are zero", which
# is the only thing generalised here. The F reference and its df2 come from the
# same place and for the same small-sample reason documented there.
#
# The identity contrast reduces exactly to wald_test(), which is how the unit
# test pins this down rather than trusting the algebra.
# ---------------------------------------------------------------------------
contrast_wald <- function(fit, contrast) {
  if (!isTRUE(fit$estimable)) {
    return(list(estimable = FALSE, failure_reason = fit$failure_reason))
  }
  terms <- colnames(contrast)
  if (!all(terms %in% names(fit$coefficients))) {
    return(list(estimable = FALSE, failure_reason = "contrast_term_absent_from_design"))
  }
  b <- fit$coefficients[terms]
  V <- fit$vcov[terms, terms, drop = FALSE]
  Cb <- contrast %*% b
  CVC <- contrast %*% V %*% t(contrast)
  CVCinv <- tryCatch(solve(CVC), error = function(e) NULL)
  if (is.null(CVCinv)) {
    return(list(estimable = FALSE, failure_reason = "singular_contrast_covariance"))
  }
  df1 <- nrow(contrast)
  statistic <- as.numeric(t(Cb) %*% CVCinv %*% Cb) / df1
  list(estimable = TRUE, failure_reason = NA_character_,
       f_statistic = statistic, df1 = df1, df2 = fit$df,
       p_value = stats::pf(statistic, df1 = df1, df2 = fit$df, lower.tail = FALSE))
}

CONTRAST_OMNIBUS <- local({
  m <- diag(3L); colnames(m) <- T2_HISTOLOGY_COMPONENTS; m
})
# Equal weighting: steatosis - ballooning = 0 and ballooning - inflammation = 0.
# Together these say all three coefficients are equal, which is exactly the claim
# that the NAS composite's own equal-weight sum is adequate for the program.
CONTRAST_EQUAL_WEIGHT <- local({
  m <- rbind(c(1, -1, 0), c(0, 1, -1))
  colnames(m) <- T2_HISTOLOGY_COMPONENTS; m
})

# ---------------------------------------------------------------------------
# How many independent things are these program scores?
#
# A count of supported programs is not a count of findings if the scores are
# near-copies of one another. The frozen programs were built on single-cell
# expression and nothing guarantees they stay distinct once projected onto 3,991
# liver proteins, so the effective dimension of the score matrix is reported
# next to every count that uses it. The participation ratio is used because it
# does not need a variance cutoff chosen after the fact.
# ---------------------------------------------------------------------------
t2_score_dimensionality <- function(scores) {
  X <- t(scores)
  correlation <- stats::cor(X)
  eigenvalues <- stats::prcomp(X, center = TRUE, scale. = TRUE)$sdev^2
  offdiag <- abs(correlation[upper.tri(correlation)])
  data.table(
    n_programs = ncol(X),
    n_donors = nrow(X),
    pc1_variance_explained = eigenvalues[[1]] / sum(eigenvalues),
    pc1_to_pc3_variance_explained = sum(eigenvalues[1:min(3L, length(eigenvalues))]) /
      sum(eigenvalues),
    effective_dimension_participation_ratio =
      sum(eigenvalues)^2 / sum(eigenvalues^2),
    median_abs_pairwise_correlation = stats::median(offdiag),
    max_abs_pairwise_correlation = max(offdiag)
  )
}

# Which of the three components does a program track once they are separated?
# Reported only from the joint (mutually adjusted) view, because the marginal
# view cannot distinguish a component from its correlates.
t2_assign_feature <- function(joint_rows) {
  d <- joint_rows[term %in% T2_HISTOLOGY_COMPONENTS]
  d[, .(
    n_supported = sum(support == "supported"),
    tracks = {
      hit <- term[support == "supported"]
      if (!length(hit)) NA_character_ else paste(sort(hit), collapse = "+")
    },
    n_informative_null = sum(support == "informative_null"),
    n_indeterminate = sum(support == "indeterminate"),
    min_mde = min(mde, na.rm = TRUE),
    max_mde = max(mde, na.rm = TRUE)
  ), by = program_uid]
}
