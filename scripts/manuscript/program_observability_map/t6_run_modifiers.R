#!/usr/bin/env Rscript

# T6: is any frozen program's fibrosis or NAS association modified by sex, age
# or cohort?
#
# WHAT THE DISCOVERY RELEASE ASSUMES AND THIS SYSTEM CHECKS
# The discovery axis map reports one slope per program per axis, adjusted for
# sex, and treats it as a single number that applies to every donor in every
# cohort. That is an assumption, not a finding. If the slope for a program were
# twice as steep in men as in women, or absent in one cohort, the meta-analytic
# average would still be reported and would describe nobody. This system tests
# the assumption on all three modifiers the discovery design can carry.
#
# THE EXPECTATION IS A NULL, AND THAT IS THE DELIVERABLE
# An interaction test is far weaker than the main effect it modifies. This
# repository has already retired one sex-dimorphic result that collapsed from 290
# genes to 8 once cohort was held fixed, and a narrative built on top of it. The
# useful output here is therefore not a count of modified programs; it is the
# size of the modification the design could have resolved, attached to every
# program that returns nothing. A null without that number is not a finding.
#
# WHAT IS AND IS NOT ESTIMABLE ON THIS SUBSTRATE
#   sex    - all four discovery cohorts, 469 donors. GSE135251, the largest
#            cohort at 214 donors, has NO RECORDED SEX; its sex_final comes from
#            inferred XIST/DDX3Y k-means. A sensitivity arm drops it and refits
#            on the 255 donors whose sex is recorded, at the cost of falling to
#            the three-cohort minimum.
#   age    - three cohorts, 255 donors. GSE135251 records no age at all, so the
#            largest cohort cannot enter this arm and the meta-analysis sits at
#            its three-cohort floor with two Knapp-Hartung degrees of freedom.
#   cohort - a heterogeneity question, answered by Cochran's Q on the per-cohort
#            slopes the discovery release already fitted, plus a simulated
#            minimum detectable between-cohort standard deviation.
#
# LANGUAGE
# Fibrosis stage and NAS are cross-sectional recorded histologic variables and
# every slope here is a cross-sectional association. Nothing in this file
# implies ordering, movement, or a subgroup of donors defined by anything other
# than the recorded modifier.

suppressPackageStartupMessages({
  library(edgeR)
  library(data.table)
  library(parallel)
})

args <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", args[grepl("^--file=", args)])))
source(file.path(script_dir, "config.R"))
source(file.path(script_dir, "analysis_lib.R"))
source(file.path(script_dir, "calibration_lib.R"))
source(file.path(script_dir, "systems_contract.R"))
source(file.path(script_dir, "t6_modifier_lib.R"))

set.seed(seed)
n_modifier_null <- as.integer(Sys.getenv("T6_MODIFIER_PERMUTATIONS", "1000"))
n_diagnostic_null <- as.integer(Sys.getenv("T6_DIAGNOSTIC_PERMUTATIONS", "300"))
n_heterogeneity_null <- as.integer(Sys.getenv("T6_HETEROGENEITY_PERMUTATIONS", "1000"))
n_tau_sims <- as.integer(Sys.getenv("T6_TAU_SIMULATIONS", "4000"))
n_cores <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))

system_root <- Sys.getenv("T6_OUTPUT_ROOT",
                          file.path(workstream_root, "systems", "T6_modifiers"))
if (dir.exists(system_root)) fail("Refusing to overwrite: ", system_root)

assert_inclusion_criterion(c("evidence_interpretation", "observability"))

# Prespecified before the run, and identical in all three arms so the three
# nulls mean the same thing. 0.20 program-score SD is the discovery release's
# observability threshold: it sits just below the median absolute main effect of
# the supported set, so it asks whether the design could have resolved a
# modification the size of the associations it did resolve.
#
#   sex    : difference in slope between men and women, per unit of the axis.
#   age    : change in slope per DECADE of age, per unit of the axis. Age is
#            coded in decades for exactly this reason; per-year units would make
#            the threshold unreadable and invite a comparison against a number
#            on a different scale.
#   cohort : between-cohort standard deviation of the slope.
t6_effect_threshold <- observability_effect_threshold
t6_family_size <- expected_testable_programs

# ---------------------------------------------------------------------- inputs

message("[1/12] Loading sealed inputs")
if (!file.exists(nonholdout_dge) || !file.exists(nonholdout_meta)) {
  fail("Sealed non-holdout partition missing")
}
if (!dir.exists(discovery_root)) fail("Discovery release missing: ", discovery_root)

dge <- readRDS(nonholdout_dge)
meta <- as.data.table(readRDS(nonholdout_meta))
assert_true(ncol(dge) == 1097L, "Non-holdout DGE must contain 1,097 samples")
assert_true(identical(colnames(dge), meta$sample_id), "DGE and metadata sample order mismatch")
assert_holdout_sealed(meta, holdout_cohort, gse193066_crosswalk)
frozen <- assert_frozen_programs(program_registry, program_membership, 117L)

discovery_meta <- meta[
  dataset %in% discovery_cohorts & !is.na(fibrosis_stage) & !is.na(nas_score)
]
for (cohort in names(discovery_expected_n)) {
  assert_true(discovery_meta[dataset == cohort, .N] == discovery_expected_n[[cohort]],
              paste0("Discovery count drift for ", cohort))
}
assert_true(nrow(discovery_meta) == 469L, "Joint discovery donor count must be 469")

# The biological unit. sample_id is the unit column and the uniqueness check on
# it would be vacuous on its own, so pairing_root is supplied and the absence of
# any donor-pairing table for these cohorts is verified on disk instead.
assert_biological_unit(discovery_meta, "sample_id", cohorts = discovery_cohorts,
                       pairing_root = file.path(project_root, "data"))
assert_true(uniqueN(discovery_meta$sample_id) == 469L, "Donor identifiers must be unique")

discovery_meta[, sex_male := as.integer(sex_final == "M")]
assert_true(all(!is.na(discovery_meta$sex_male)), "sex_final must be complete on the discovery set")
discovery_meta[, age := as.numeric(age)]
discovery_meta[, age_decades := age / 10]

sex_provenance <- discovery_meta[, .(
  n = .N, n_male = sum(sex_male), female_fraction = mean(1 - sex_male),
  sex_source = paste(sort(unique(sex_source)), collapse = "/"),
  n_age = sum(is.finite(age)),
  age_mean = mean(age, na.rm = TRUE), age_sd = stats::sd(age, na.rm = TRUE)
), by = dataset][order(dataset)]
message("      sex and age provenance by cohort:")
print(sex_provenance)

age_cohorts <- sex_provenance[n_age == n, dataset]
assert_true(length(age_cohorts) >= 3L,
            "The age arm needs at least three cohorts with complete age")
recorded_sex_cohorts <- sex_provenance[sex_source == "annotated", dataset]
assert_true(length(recorded_sex_cohorts) >= 3L,
            "The recorded-sex sensitivity arm needs at least three cohorts")
inferred_sex_cohorts <- setdiff(discovery_cohorts, recorded_sex_cohorts)

# --------------------------------------------------------------------- scoring

message("[2/12] Projecting the frozen programs (identical to the discovery path)")
logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
annotation <- fread(gene_annotation, select = c("gene_id", "gene_name"))
symbol_matrix <- collapse_symbols(logcpm, rownames(logcpm), annotation)
rm(logcpm, dge); gc()

program_scores <- score_programs(symbol_matrix, meta, frozen$membership, frozen$registry,
                                 program_weight_coverage)
assert_true(sum(program_scores$coverage$testable) == expected_testable_programs,
            "Testable program count drifted from the discovery release")
testable_programs <- program_scores$coverage[testable == TRUE, feature_id]
S <- program_scores$primary[testable_programs, discovery_meta$sample_id, drop = FALSE]
assert_true(all(is.finite(S)), "A testable program has a non-finite discovery score")
rm(symbol_matrix); gc()

# ------------------------------------------------- provenance against discovery

message("[3/12] Reproducing the discovery linear fits before modifying them")
primary_cohort <- fit_feature_models(S, discovery_meta, discovery_cohorts, "primary")
primary_meta <- meta_analyze(primary_cohort)
discovery_meta_effects <- fread(file.path(discovery_root, "meta_effects.tsv"))
reference <- discovery_meta_effects[model_kind == "primary" & estimable == TRUE,
                                    .(feature_id, axis, beta_meta, se_meta, tau2, p_heterogeneity)]
reproduced <- merge(primary_meta[estimable == TRUE,
                                 .(feature_id, axis, beta_meta, se_meta, tau2, p_heterogeneity)],
                    reference, by = c("feature_id", "axis"), suffixes = c("", "_discovery"))
assert_true(nrow(reproduced) == 2L * length(testable_programs),
            "Discovery cross-check did not match every program-axis")
reproduction_error <- reproduced[, max(
  abs(beta_meta - beta_meta_discovery), abs(se_meta - se_meta_discovery),
  abs(tau2 - tau2_discovery), abs(p_heterogeneity - p_heterogeneity_discovery))]
message("      worst absolute difference from the sealed discovery release: ",
        format(reproduction_error, digits = 3))
assert_true(reproduction_error < 1e-8,
            paste0("This system does not reproduce the discovery linear fits (worst ",
                   reproduction_error, "); the modifier arms would not be comparable"))

# ------------------------------------------------------------------- the arms

fit_arm <- function(meta_table, cohorts, modifier_column, include_sex_covariate, label) {
  cohort_effects <- t6_fit_cohort_interactions(
    S, meta_table, cohorts, modifier_column, include_sex_covariate, label)
  meta_effects <- meta_analyze(cohort_effects, min_cohorts = 3L)
  meta_effects[, arm := label]
  list(cohort = cohort_effects[, arm := label][], meta = meta_effects)
}

message("[4/12] Fitting the sex arm on all four cohorts")
sex_arm <- fit_arm(discovery_meta, discovery_cohorts, "sex_male", FALSE, "sex")

message("[5/12] Fitting the recorded-sex sensitivity arm (",
        paste(recorded_sex_cohorts, collapse = ", "), ")")
sex_recorded_arm <- fit_arm(discovery_meta[dataset %in% recorded_sex_cohorts],
                            recorded_sex_cohorts, "sex_male", FALSE, "sex_recorded_only")

message("[6/12] Fitting the age arm on ", paste(age_cohorts, collapse = ", "))
age_meta_table <- discovery_meta[dataset %in% age_cohorts & is.finite(age_decades)]
message("      age arm donors: ", nrow(age_meta_table))
age_arm <- fit_arm(age_meta_table, age_cohorts, "age_decades", TRUE, "age")

# --------------------------------------------------------------- permutation

# One permutation replicate of a modifier arm. The stratum is the axis being
# tested, so the modifier's joint distribution with that axis survives the
# permutation and only the interaction is removed.
null_replicate <- function(meta_table, cohorts, modifier_column, include_sex_covariate,
                           axis_name, stratum_column, replicate_seed) {
  set.seed(replicate_seed)
  permuted <- t6_permute_modifier(meta_table, modifier_column,
                                  stratum_column = stratum_column)
  cohort_effects <- t6_fit_cohort_interactions(
    scores = S, meta = permuted, cohorts = cohorts,
    modifier_column = modifier_column,
    include_sex_covariate = include_sex_covariate, model_kind = "null")
  result <- meta_analyze(cohort_effects[axis == axis_name], min_cohorts = 3L)
  result[match(testable_programs, feature_id), p_value_meta]
}

null_campaign <- function(meta_table, cohorts, modifier_column, include_sex_covariate,
                          axis_name, stratum_column, n_reps, seed_offset, label) {
  message("      null campaign: ", label, " (", n_reps, " permutations, ",
          n_cores, " cores)")
  out <- mclapply(seq_len(n_reps), function(i) {
    null_replicate(meta_table, cohorts, modifier_column, include_sex_covariate,
                   axis_name, stratum_column, seed + seed_offset + i)
  }, mc.cores = n_cores)
  bad <- vapply(out, function(x) !is.numeric(x) || length(x) != length(testable_programs),
                logical(1))
  if (any(bad)) fail("Permutation replicate failed in campaign ", label, ": ",
                     paste(utils::head(which(bad), 3L), collapse = ", "))
  out
}

axis_strata <- c(fibrosis = "fibrosis_stage", nas = "nas_score")

message("[7/12] Measuring the permutation null for every modifier view")
calibration_rows <- list()
gated_calibrations <- list()
observed_p_store <- list()
uniformity_rows <- list()
null_store <- list()
seed_offset <- 0L

arms <- list(
  list(label = "sex", arm = sex_arm, meta_table = discovery_meta,
       cohorts = discovery_cohorts, modifier = "sex_male", sex_covariate = FALSE,
       gated = TRUE),
  list(label = "sex_recorded_only", arm = sex_recorded_arm,
       meta_table = discovery_meta[dataset %in% recorded_sex_cohorts],
       cohorts = recorded_sex_cohorts, modifier = "sex_male", sex_covariate = FALSE,
       gated = TRUE),
  list(label = "age", arm = age_arm, meta_table = age_meta_table,
       cohorts = age_cohorts, modifier = "age_decades", sex_covariate = TRUE,
       gated = TRUE)
)

for (spec in arms) {
  for (axis_name in names(axis_strata)) {
    view <- paste0(spec$label, "_x_", axis_name)
    observed <- spec$arm$meta[axis == axis_name][match(testable_programs, feature_id)]
    permutable <- t6_permutable_fraction(spec$meta_table,
                                         stratum_column = axis_strata[[axis_name]])
    seed_offset <- seed_offset + 100000L
    stratified <- null_campaign(spec$meta_table, spec$cohorts, spec$modifier,
                                spec$sex_covariate, axis_name,
                                axis_strata[[axis_name]], n_modifier_null,
                                seed_offset, paste0(view, " [stratified]"))
    seed_offset <- seed_offset + 100000L
    unstratified <- null_campaign(spec$meta_table, spec$cohorts, spec$modifier,
                                  spec$sex_covariate, axis_name, NULL,
                                  n_diagnostic_null, seed_offset,
                                  paste0(view, " [unstratified diagnostic]"))
    calib <- calibrate_view(observed$p_value_meta, function(i) stratified[[i]],
                            n_reps = n_modifier_null, n_tests = t6_family_size,
                            label = view)
    calib_diag <- calibrate_view(observed$p_value_meta, function(i) unstratified[[i]],
                                 n_reps = n_diagnostic_null, n_tests = t6_family_size,
                                 label = paste0(view, "__unstratified_diagnostic"))
    row <- calibration_row(calib)
    row[, `:=`(null_scheme = "modifier permuted within cohort and axis level",
               permutable_fraction = permutable,
               unstratified_null_call_mean = calib_diag$null_call_mean,
               gated = spec$gated)]
    calibration_rows[[length(calibration_rows) + 1L]] <- row
    gated_calibrations[[view]] <- calib
    calibration_rows[[length(calibration_rows) + 1L]] <- calibration_row(calib_diag)[
      , `:=`(null_scheme = "modifier permuted within cohort only (diagnostic)",
             permutable_fraction = t6_permutable_fraction(spec$meta_table),
             unstratified_null_call_mean = NA_real_, gated = FALSE)]
    uniformity_rows[[length(uniformity_rows) + 1L]] <-
      t6_uniformity(observed$p_value_meta, view)
    null_store[[view]] <- do.call(rbind, stratified)
    observed_p_store[[view]] <- observed$p_value_meta
    message("      ", view, ": observed BH calls ", calib$observed_calls,
            " | null mean ", format(calib$null_call_mean, digits = 3),
            " (max ", calib$null_call_max, ") | unstratified null mean ",
            format(calib_diag$null_call_mean, digits = 3),
            " | reportable ", calib$reportable)
  }
}

# ------------------------------------------------------------- cohort arm null

message("[8/12] Measuring the permutation null for the cohort-heterogeneity views")
# The primary model returns both axes from one fit, so one replicate serves both
# views. The null is the joint within-cohort histology permutation from
# calibration_lib: it removes the association with expression while preserving
# the fibrosis-NAS correlation, which leaves Q with its nominal null.
heterogeneity_replicate <- function(replicate_seed) {
  set.seed(replicate_seed)
  permuted <- permute_histology_within_cohort(
    discovery_meta, c("fibrosis_stage", "nas_score"), "dataset")
  cohort_effects <- fit_feature_models(S, permuted, discovery_cohorts, "primary")
  result <- meta_analyze(cohort_effects)
  lapply(c(fibrosis = "fibrosis", nas = "nas"), function(a) {
    result[axis == a][match(testable_programs, feature_id), p_heterogeneity]
  })
}

seed_offset <- seed_offset + 100000L
message("      null campaign: cohort heterogeneity (", n_heterogeneity_null,
        " permutations, ", n_cores, " cores)")
heterogeneity_null <- mclapply(seq_len(n_heterogeneity_null), function(i) {
  heterogeneity_replicate(seed + seed_offset + i)
}, mc.cores = n_cores)

for (axis_name in c("fibrosis", "nas")) {
  view <- paste0("cohort_heterogeneity_", axis_name)
  null_p <- lapply(heterogeneity_null, function(x) x[[axis_name]])
  bad <- vapply(null_p, function(x) !is.numeric(x) || length(x) != length(testable_programs),
                logical(1))
  if (any(bad)) fail("Heterogeneity permutation replicate failed for ", view)
  observed <- primary_meta[axis == axis_name][match(testable_programs, feature_id)]
  calib <- calibrate_view(observed$p_heterogeneity, function(i) null_p[[i]],
                          n_reps = n_heterogeneity_null, n_tests = t6_family_size,
                          label = view)
  calibration_rows[[length(calibration_rows) + 1L]] <- calibration_row(calib)[
    , `:=`(null_scheme = "histologic axes jointly permuted within cohort",
           permutable_fraction = 1, unstratified_null_call_mean = NA_real_,
           gated = TRUE)]
  gated_calibrations[[view]] <- calib
  uniformity_rows[[length(uniformity_rows) + 1L]] <-
    t6_uniformity(observed$p_heterogeneity, view)
  null_store[[view]] <- do.call(rbind, null_p)
  observed_p_store[[view]] <- observed$p_heterogeneity
  message("      ", view, ": observed BH calls ", calib$observed_calls,
          " | null mean ", format(calib$null_call_mean, digits = 3),
          " (max ", calib$null_call_max, ") | reportable ", calib$reportable)
}

calibration_rows <- rbindlist(calibration_rows, use.names = TRUE, fill = TRUE)
uniformity <- rbindlist(uniformity_rows, use.names = TRUE, fill = TRUE)

# Family-level signal, against each view's own permutation null. This is the
# check that stops "zero programs survive BH" from being reported as "there is
# nothing here": a spread of modifications too small to clear multiplicity would
# still move these two statistics.
family_signal <- rbindlist(lapply(names(observed_p_store), function(view) {
  t6_family_signal(observed_p_store[[view]], null_store[[view]], view)
}))
family_signal[, `:=`(
  q_excess_small = p.adjust(p_excess_small, method = "BH"),
  q_strongest_program = p.adjust(p_strongest_program, method = "BH"))]
uniformity <- merge(uniformity, family_signal, by = "view", all.x = TRUE)

# ------------------------------------------------------------ program results

message("[9/12] Assembling per-program results with minimum detectable effects")

# Two minimum detectable effects per program, because the reported count is a BH
# count and the alpha = 0.05 figure the discovery release quotes is not the
# threshold a program in a 113-member family actually has to clear.
#   mde_nominal : alpha 0.05, the discovery release's own definition, comparable
#                 to the main-effect MDE printed in the axis map.
#   mde_family  : alpha 0.05/113, the size a program must reach to be called when
#                 it is the only signal in the family. Conservative.
# The synthetic planted-interaction test in t6_run_tests.R confirms the size BH
# actually resolves falls between these two.
assemble_arm <- function(arm_result, label, modifier_units) {
  out <- copy(arm_result$meta)
  out[, `:=`(
    view = paste0(label, "_x_", axis),
    modifier = label,
    modifier_units = modifier_units,
    mde_nominal = minimum_detectable_effect(se_meta, n_cohorts),
    mde_family = minimum_detectable_effect(se_meta, n_cohorts,
                                           alpha = 0.05 / t6_family_size)
  )]
  out[estimable == TRUE, q_value := p.adjust(p_value_meta, method = "BH",
                                             n = t6_family_size), by = axis]
  out[, support := t6_support_state(q_value, ci_lower, ci_upper, estimable,
                                    t6_effect_threshold)]
  out[, effect_threshold := t6_effect_threshold]
  out[]
}

sex_results <- assemble_arm(sex_arm, "sex", "program-score SD per axis unit, male minus female")
sex_recorded_results <- assemble_arm(sex_recorded_arm, "sex_recorded_only",
                                     "program-score SD per axis unit, male minus female")
age_results <- assemble_arm(age_arm, "age",
                            "program-score SD per axis unit per decade of age")
modifier_results <- rbindlist(list(sex_results, sex_recorded_results, age_results),
                              use.names = TRUE, fill = TRUE)

message("[10/12] Simulating the minimum detectable between-cohort SD for the cohort arm")
tau_grid <- seq(0.01, 1.50, by = 0.01)
cohort_se <- primary_cohort[estimable == TRUE & is.finite(se_hc3) & se_hc3 > 0,
                            .(feature_id, axis, cohort, se_hc3)]
tau_keys <- unique(cohort_se[, .(feature_id, axis)])
tau_values <- unlist(mclapply(seq_len(nrow(tau_keys)), function(i) {
  set.seed(seed + 900000L + i)
  se <- cohort_se[feature_id == tau_keys$feature_id[[i]] & axis == tau_keys$axis[[i]], se_hc3]
  t6_minimum_detectable_tau(se, tau_grid, n_sims = n_tau_sims)
}, mc.cores = n_cores))
tau_keys[, mde_tau := as.numeric(tau_values)]
# A program whose standard errors are so wide that no tau on the grid reaches 80
# percent power still has a minimum detectable effect statement: it is "larger
# than the top of the grid". Recording that as NA would make the negative
# uninterpretable, which is exactly what the contract forbids, so it is recorded
# as a censored lower bound instead.
tau_keys[, `:=`(
  mde_tau_censored = !is.finite(mde_tau),
  mde_tau_lower_bound = fifelse(is.finite(mde_tau), mde_tau, max(tau_grid))
)]

cohort_results <- merge(
  primary_meta[estimable == TRUE, .(
    feature_id, axis, beta_meta, se_meta, ci_lower, ci_upper, p_value_meta,
    tau2, i2, q_heterogeneity, p_heterogeneity, n_cohorts, n_same_direction)],
  tau_keys, by = c("feature_id", "axis"), all.x = TRUE)
cohort_results[, `:=`(view = paste0("cohort_heterogeneity_", axis), modifier = "cohort",
                      modifier_units = "between-cohort SD of the slope, program-score SD per axis unit",
                      tau = sqrt(pmax(tau2, 0)))]
cohort_results[, q_value := p.adjust(p_heterogeneity, method = "BH", n = t6_family_size),
               by = axis]
cohort_results[, support := fcase(
  is.finite(q_value) & q_value < 0.05, "modified",
  is.finite(mde_tau) & tau < t6_effect_threshold & mde_tau <= t6_effect_threshold,
  "informative_null",
  default = "indeterminate")]
cohort_results[, effect_threshold := t6_effect_threshold]

# The contract's negative-carries-an-MDE rule, enforced on both result tables and
# on both MDE definitions, so no negative can be quoted without one.
assert_negatives_have_mde(modifier_results, "support", "mde_nominal")
assert_negatives_have_mde(modifier_results, "support", "mde_family")
assert_negatives_have_mde(cohort_results, "support", "mde_tau_lower_bound")

# ------------------------------------------------------------------- counting

message("[11/12] Gating and counting")
# min() over an all-NA column returns Inf with a warning; a family in which no
# program is estimable should report NA, not a sentinel that looks like a number.
safe_min <- function(x) if (any(is.finite(x))) min(x[is.finite(x)]) else NA_real_
safe_max <- function(x) if (any(is.finite(x))) max(x[is.finite(x)]) else NA_real_
# gate_view() is called non-strict so that a failing view withholds its own count
# instead of aborting every other view's. Refusing to report is the intended
# behaviour of the gate, not an error to work around.
gate_reasons <- rbindlist(lapply(gated_calibrations, function(calib) {
  report <- gate_view(calib, strict = FALSE)
  data.table(view = calib$label, passed = report$passed,
             reason = if (length(report$reasons)) paste(report$reasons, collapse = "; ")
             else NA_character_)
}))
print(gate_reasons)
withheld <- gate_reasons[passed == FALSE, view]
if (length(withheld)) {
  message("      WITHHELD (calibration gate failed): ", paste(withheld, collapse = ", "))
}

reportable_views <- calibration_rows[gated == TRUE & reportable == TRUE, view]
counts <- rbindlist(list(
  modifier_results[view %in% reportable_views, .(
    modifier = unique(modifier), n_tested = sum(estimable), n_cohorts = max(n_cohorts),
    n_modified = sum(support == "modified"),
    n_informative_null = sum(support == "informative_null"),
    n_indeterminate = sum(support == "indeterminate"),
    n_untestable = sum(support == "untestable"),
    median_se_meta = stats::median(se_meta, na.rm = TRUE),
    median_mde_nominal = stats::median(mde_nominal, na.rm = TRUE),
    median_mde_family = stats::median(mde_family, na.rm = TRUE),
    max_abs_beta = safe_max(abs(beta_meta)),
    min_q = safe_min(q_value)), by = view],
  cohort_results[view %in% reportable_views, .(
    modifier = "cohort", n_tested = .N, n_cohorts = max(n_cohorts),
    n_modified = sum(support == "modified"),
    n_informative_null = sum(support == "informative_null"),
    n_indeterminate = sum(support == "indeterminate"),
    n_untestable = 0L,
    median_se_meta = stats::median(se_meta, na.rm = TRUE),
    median_mde_nominal = stats::median(mde_tau, na.rm = TRUE),
    median_mde_family = NA_real_,
    max_abs_beta = safe_max(tau),
    min_q = safe_min(q_value)), by = view]
), use.names = TRUE, fill = TRUE)

assert_counts_calibrated(counts, calibration_rows[view %in% counts$view])

cohort_tau_summary <- cohort_results[, .(
  n = .N, n_tau2_zero = sum(tau2 <= .Machine$double.eps),
  median_tau = stats::median(tau), max_tau = max(tau),
  median_i2 = stats::median(i2, na.rm = TRUE),
  median_mde_tau = stats::median(mde_tau, na.rm = TRUE),
  median_mde_tau_lower_bound = stats::median(mde_tau_lower_bound),
  n_mde_tau_censored = sum(mde_tau_censored),
  n_mde_tau_below_threshold = sum(is.finite(mde_tau) & mde_tau <= t6_effect_threshold)
), by = axis]

# ---------------------------------------------------------------------- seal

message("[12/12] Sealing")
tmp <- atomic_dir(system_root)
on.exit(if (dir.exists(tmp)) unlink(tmp, recursive = TRUE), add = TRUE)

write_tsv(modifier_results, file.path(tmp, "modifier_interaction_meta.tsv"))
write_tsv(rbindlist(list(sex_arm$cohort, sex_recorded_arm$cohort, age_arm$cohort),
                    use.names = TRUE, fill = TRUE),
          file.path(tmp, "modifier_interaction_per_cohort.tsv"))
write_tsv(cohort_results, file.path(tmp, "cohort_heterogeneity.tsv"))
write_tsv(cohort_tau_summary, file.path(tmp, "cohort_heterogeneity_summary.tsv"))
write_tsv(calibration_rows, file.path(tmp, "calibration.tsv"))
write_tsv(uniformity, file.path(tmp, "pvalue_uniformity_and_family_signal.tsv"))
write_tsv(counts, file.path(tmp, "counts.tsv"))
write_tsv(sex_provenance, file.path(tmp, "modifier_provenance_by_cohort.tsv"))
write_tsv(primary_meta, file.path(tmp, "reproduced_discovery_linear_meta.tsv"))
saveRDS(null_store, file.path(tmp, "permutation_null_pvalues.rds"), compress = "xz")
writeLines(capture.output(sessionInfo()), file.path(tmp, "sessionInfo.txt"))

inputs <- c(nonholdout_dge, nonholdout_meta, program_registry, program_membership,
            gene_annotation, file.path(discovery_root, "meta_effects.tsv"),
            file.path(discovery_root, "program_axis_map.tsv"))
write_tsv(data.table(path = inputs, sha256 = vapply(inputs, sha256_file, character(1))),
          file.path(tmp, "input_manifest.tsv"))

estimand_text <- paste(
  "Difference between two cross-sectional associations. For each frozen program and",
  "each recorded histologic axis, the change in the program-score slope per unit of",
  "the modifier: male minus female, per decade of age, or between cohorts. Sex and",
  "age interactions are fitted inside each cohort and then meta-analysed, never",
  "pooled, because the cohorts differ in sex ratio and stage composition at once and",
  "a pooled product term cannot separate an interaction from that difference. The",
  "cohort arm is Cochran heterogeneity of the per-cohort slopes, not a refitted",
  "cohort product term.")
assert_language(c(estimand_text, modifier_results$support, cohort_results$support,
                  counts$view, calibration_rows$view, calibration_rows$null_scheme,
                  uniformity$view, modifier_results$modifier_units,
                  cohort_results$modifier_units))

jsonlite::write_json(list(
  state = "T6_MODIFIERS_COMPLETE_HOLDOUT_UNOPENED",
  release_id = release_id, workstream_id = workstream_id, system = "T6_modifiers",
  estimand = estimand_text,
  n_donors_sex_arm = nrow(discovery_meta),
  n_donors_sex_recorded_arm = nrow(discovery_meta[dataset %in% recorded_sex_cohorts]),
  n_donors_age_arm = nrow(age_meta_table),
  cohorts_sex_arm = discovery_cohorts,
  cohorts_age_arm = age_cohorts,
  cohorts_with_inferred_sex = inferred_sex_cohorts,
  sex_inference_caveat = paste(
    "GSE135251, the largest discovery cohort at 214 donors, has no recorded sex and",
    "no recorded age. Its sex_final is inferred from XIST/DDX3Y k-means, so the",
    "four-cohort sex arm rests on inferred sex for 46 percent of its donors. The",
    "sex_recorded_only arm refits on the 255 donors with recorded sex."),
  n_testable_programs = length(testable_programs),
  effect_threshold = t6_effect_threshold,
  multiplicity_family_size = t6_family_size,
  n_modifier_permutations = n_modifier_null,
  n_diagnostic_permutations = n_diagnostic_null,
  n_heterogeneity_permutations = n_heterogeneity_null,
  n_tau_simulations = n_tau_sims,
  discovery_reproduction_error = reproduction_error,
  withheld_views = withheld,
  holdout_accessed = FALSE
), file.path(tmp, "T6_READY.json"), pretty = TRUE, auto_unbox = TRUE)

artifacts <- setdiff(list.files(tmp, recursive = TRUE), "artifact_manifest.tsv")
write_tsv(data.table(artifact = artifacts,
                     sha256 = vapply(file.path(tmp, artifacts), sha256_file, character(1)),
                     size_bytes = file.info(file.path(tmp, artifacts))$size),
          file.path(tmp, "artifact_manifest.tsv"))
publish_dir(tmp, system_root)
Sys.chmod(list.files(system_root, full.names = TRUE), mode = "0440")

cat("\n=== T6 PROVENANCE ===\n"); print(sex_provenance)
cat("\n=== T6 CALIBRATION ===\n")
print(calibration_rows[, .(view, observed_calls, null_call_mean, null_call_max,
                           unstratified_null_call_mean, permutable_fraction,
                           calibrated, reportable, gated)])
cat("\n=== T6 COUNTS ===\n"); print(counts)
cat("\n=== T6 SUPPORT STATES ===\n")
print(modifier_results[, .N, by = .(view, support)][order(view, support)])
print(cohort_results[, .N, by = .(view, support)][order(view, support)])
cat("\n=== T6 COHORT HETEROGENEITY ===\n"); print(cohort_tau_summary)
cat("\n=== T6 FAMILY-LEVEL SIGNAL AND P-VALUE UNIFORMITY ===\n")
print(uniformity[, .(view, n, median_p, fraction_below_0.05,
                     null_median_fraction_below_0.05, p_excess_small, q_excess_small,
                     observed_min_p, p_strongest_program, q_strongest_program)])
cat("T6_COMPLETE\t", system_root, "\n", sep = "")
