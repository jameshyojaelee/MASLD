#!/usr/bin/env Rscript

# T5: how much of each frozen program's fibrosis association is coupled to
# measured cell-type composition.
#
# THE QUESTION, AND THE ONE IT IS NOT
# The composition proportions are deconvolved from the same logCPM matrix the
# program scores are projected from. A program built out of a lineage's marker
# genes therefore shares variance with that lineage's own estimated abundance
# partly BY CONSTRUCTION, before any biology is involved. Nothing in this system
# can separate that mechanical share from the biological one. Every quantity here
# is therefore an UPPER BOUND on composition attribution and is reported as
# coupling. The causal version of the question is answered independently by a
# sibling system on within-cell-type single-cell scores.
#
# WHAT IS NEW RELATIVE TO THE DISCOVERY RELEASE
# The discovery `composition_pc` arm reports one number per program: the fibrosis
# coefficient after adjustment. That confounds three different things -- real
# shared variance, the variance cost of three extra parameters, and collinearity
# between the composition block and fibrosis. This system separates them:
#
#   1. A commonality partition of R-squared into unique-fibrosis,
#      unique-composition and shared, per program, with the shared term signed so
#      suppression stays visible.
#   2. Bootstrap confidence intervals by resampling donors within cohort, because
#      R-squared fractions have no clean closed-form standard error.
#   3. TWO same-cohort, same-donor control arms with exactly three extra
#      parameters each. A prior attempt at this comparison used an `age` model on
#      three cohorts as its control; that control already returned 0 of 113 at
#      min q 0.078, so it had no power to show anything and the comparison was
#      uninterpretable. Both controls here are fitted on the identical 469 donors
#      and identical four cohorts as the test arm.
#        - composition_permuted: composition columns permuted within cohort.
#          Isolates the pure degrees-of-freedom cost of three noise parameters.
#        - composition_stage_matched: permuted within cohort AND fibrosis stage.
#          Preserves whatever marginal association composition has with stage, so
#          it isolates the degrees-of-freedom cost PLUS the collinearity cost,
#          leaving donor-level coupling as the only thing the real arm adds.
#
# CALIBRATION
# Each countable view gets the permutation null that matches its own claim, not a
# generic one. A view asking whether fibrosis contributes uniquely is nulled by
# permuting histology; a view asking whether composition is coupled is nulled by
# permuting composition. Using one null for both would be answering the wrong
# question in one of the two cases.
#
# LANGUAGE
# Fibrosis stage is a cross-sectional recorded histologic variable. No ordering,
# no movement, and no claim about how a donor got to their stage.

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
source(file.path(script_dir, "t5_commonality_lib.R"))

set.seed(seed)
n_coupling_null <- as.integer(Sys.getenv("T5_COUPLING_PERMUTATIONS", "4000"))
n_histology_null <- as.integer(Sys.getenv("T5_HISTOLOGY_PERMUTATIONS", "500"))
n_control_draws <- as.integer(Sys.getenv("T5_CONTROL_DRAWS", "200"))
n_boot <- as.integer(Sys.getenv("T5_BOOTSTRAPS", "2000"))
n_cores <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))

system_root <- Sys.getenv("T5_OUTPUT_ROOT",
                          file.path(workstream_root, "systems", "T5_composition"))
if (dir.exists(system_root)) fail("Refusing to overwrite: ", system_root)

assert_inclusion_criterion(c("evidence_interpretation", "observability"))

# --------------------------------------------------------------------- inputs

message("[1/11] Loading sealed inputs")
if (!file.exists(composition_accepted_seal)) fail("Composition acceptance seal missing")
seal <- jsonlite::read_json(composition_accepted_seal, simplifyVector = TRUE)
assert_true(identical(seal$accepted_method, "bayesprism"),
            "The accepted composition estimator changed; the seal must be re-read")
assert_true(isFALSE(seal$logratio_specification_available),
            "The seal now permits a log-ratio specification; this system assumes CLR components")

dge <- readRDS(nonholdout_dge)
meta <- as.data.table(readRDS(nonholdout_meta))
composition <- fread(nonholdout_composition, check.names = FALSE)
assert_true(ncol(dge) == 1097L, "Non-holdout DGE must contain 1,097 samples")
assert_true(identical(colnames(dge), meta$sample_id), "DGE and metadata sample order mismatch")
assert_true(identical(composition$sample_id, meta$sample_id),
            "Composition and metadata sample order mismatch")
assert_holdout_sealed(meta, holdout_cohort, gse193066_crosswalk)

registry <- fread(program_registry)
membership <- fread(program_membership)
assert_frozen_programs(program_registry, program_membership, expected_programs = 117L)

message("[2/11] Projecting the frozen programs (identical to the discovery run)")
logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
annotation <- fread(gene_annotation, select = c("gene_id", "gene_name"))
symbol_matrix <- collapse_symbols(logcpm, rownames(logcpm), annotation)
rm(logcpm, dge); gc()
program_scores <- score_programs(symbol_matrix, meta, membership, registry,
                                 program_weight_coverage)
assert_true(sum(program_scores$coverage$testable) == expected_testable_programs,
            "Testable program count drift against the discovery release")
testable_programs <- program_scores$coverage[testable == TRUE, feature_id]
S <- program_scores$primary[testable_programs, , drop = FALSE]
rm(symbol_matrix); gc()

discovery_meta <- meta[dataset %in% discovery_cohorts &
                         !is.na(fibrosis_stage) & !is.na(nas_score)]
assert_true(nrow(discovery_meta) == 469L, "Discovery census drift")

# Biological unit, stated non-vacuously. The sealed metadata carries no
# participant column, so any "one sample per participant" assertion built from
# sample_id compares a column to itself and always passes. The substantive
# evidence is external and is checked instead: no donor-pairing table exists for
# any discovery cohort, and the only bulk cohort with repeat biopsies
# (GSE193066) is the sealed holdout and is not opened here.
assert_true(!"participant_id" %in% names(meta),
            "A participant column now exists; replace the external pairing check with it")
pairing_files <- file.path(project_root, "data", discovery_cohorts,
                           "metadata", "donor_pairing.csv")
assert_true(!any(file.exists(pairing_files)),
            paste0("A discovery cohort has a donor-pairing table, so repeat sampling is ",
                   "possible and every model here needs a donor term: ",
                   paste(pairing_files[file.exists(pairing_files)], collapse = ", ")))
message("      biological unit: 469 donors, one bulk sample each, ",
        "no donor-pairing table in any discovery cohort")

message("[3/11] Rebuilding the accepted composition covariates")
eligible <- if (length(seal$eligible_covariate_lineages))
  as.character(seal$eligible_covariate_lineages) else character(0)
lineages <- setdiff(names(composition)[vapply(composition, is.numeric, logical(1))],
                    c("sample_id", "dataset"))
composition_scores <- build_composition_scores(
  composition, meta, lineages, composition_pseudocount,
  eligible_lineages = eligible, n_components = composition_n_clr_components
)
pc_columns <- composition_scores$pc_columns
assert_true(length(pc_columns) == 3L, "Expected three CLR principal components")

composition_meta <- merge(discovery_meta, composition_scores$pcs,
                          by = "sample_id", all.x = TRUE, sort = FALSE)
setkey(composition_meta, NULL)
assert_true(!anyNA(composition_meta[, ..pc_columns]),
            "A discovery donor has no composition covariate; the decomposition would drop donors")

# --------------------------------------------------------- per-cohort designs

# Blocks are built so that cbind(C, A, B) spans exactly the design the discovery
# `composition_pc` arm fits. Step 5 checks that claim against the sealed output
# rather than asserting it.
build_blocks <- function(md, pcs) {
  md <- copy(md)
  md[, fibrosis_centered := fibrosis_stage - mean(fibrosis_stage)]
  md[, nas_centered := nas_score - mean(nas_score)]
  C <- stats::model.matrix(~ sex_final + nas_centered, data = md)
  A <- matrix(md$fibrosis_centered, ncol = 1L, dimnames = list(NULL, "fibrosis_centered"))
  B <- as.matrix(pcs)
  colnames(B) <- pc_columns
  list(C = C, A = A, B = B, meta = md)
}

cohort_rows <- lapply(discovery_cohorts, function(cohort) which(composition_meta$dataset == cohort))
names(cohort_rows) <- discovery_cohorts
Y_cohort <- lapply(discovery_cohorts, function(cohort) {
  rows <- cohort_rows[[cohort]]
  t(S[, composition_meta$sample_id[rows], drop = FALSE])
})
names(Y_cohort) <- discovery_cohorts
blocks <- lapply(discovery_cohorts, function(cohort) {
  rows <- cohort_rows[[cohort]]
  build_blocks(composition_meta[rows], composition_meta[rows, ..pc_columns])
})
names(blocks) <- discovery_cohorts
baseline_parts <- lapply(discovery_cohorts, function(cohort) {
  t5_baseline_parts(Y_cohort[[cohort]], blocks[[cohort]]$C, blocks[[cohort]]$A)
})
names(baseline_parts) <- discovery_cohorts

message("[4/11] Observed commonality decomposition")
observed_cohort <- lapply(discovery_cohorts, function(cohort) {
  fit <- t5_commonality_cohort(Y_cohort[[cohort]], blocks[[cohort]]$C,
                               blocks[[cohort]]$A, blocks[[cohort]]$B,
                               parts = baseline_parts[[cohort]])
  if (is.null(fit)) fail("Rank-deficient design in cohort ", cohort)
  fit
})
names(observed_cohort) <- discovery_cohorts
observed_agg <- t5_aggregate_cohorts(observed_cohort)
additivity_error <- max(abs(with(observed_agg,
  r2_covariates + unique_fibrosis + unique_composition + common - r2_full)))
assert_true(additivity_error < 1e-12,
            paste0("Commonality components do not sum to the full R-squared: ", additivity_error))
message("      additivity residual: ", format(additivity_error, scientific = TRUE))

message("[5/11] Provenance check against the sealed discovery arm")
sealed_effects <- fread(file.path(discovery_root, "cohort_effects.tsv"))
sealed_pc <- sealed_effects[model_kind == "composition_pc" & axis == "fibrosis"]
recomputed <- rbindlist(lapply(discovery_cohorts, function(cohort) {
  bl <- blocks[[cohort]]
  X <- cbind(bl$C, bl$A, bl$B)
  fit <- hc3_fit_matrix(t(Y_cohort[[cohort]]), X)
  data.table(feature_id = colnames(fit$coefficients), cohort = cohort,
             beta_recomputed = fit$coefficients["fibrosis_centered", ],
             se_recomputed = fit$se["fibrosis_centered", ])
}))
provenance <- merge(sealed_pc[, .(feature_id, cohort, beta, se_hc3)], recomputed,
                    by = c("feature_id", "cohort"))
assert_true(nrow(provenance) == length(testable_programs) * length(discovery_cohorts),
            "Provenance join lost rows")
beta_drift <- max(abs(provenance$beta - provenance$beta_recomputed))
se_drift <- max(abs(provenance$se_hc3 - provenance$se_recomputed))
assert_true(beta_drift < 1e-9 && se_drift < 1e-9,
            paste0("T5 design blocks do not reproduce the sealed composition_pc arm: ",
                   "max |beta| drift ", beta_drift, ", max |se| drift ", se_drift))
message("      reproduces the sealed arm exactly (max beta drift ",
        format(beta_drift, scientific = TRUE), ")")

supported_fibrosis <- fread(file.path(discovery_root, "program_axis_map.tsv"))[
  q_value_fibrosis < 0.05, feature_id]
assert_true(length(supported_fibrosis) == 27L,
            paste0("Expected 27 linear fibrosis-supported programs; found ",
                   length(supported_fibrosis)))

# ------------------------------------------------------------- null: coupling

# Composition permuted within cohort. Under this null the composition block has
# the same marginal distribution and the same correlation among its columns, but
# no link to the donor's own expression. Anything the decomposition still assigns
# to `common` or `unique_composition` is manufactured by the estimator, not
# measured, and that floor is what the observed values are read against.
message("[6/11] Composition-permutation null (", n_coupling_null, " draws)")
coupling_null_one <- function(draw_seed) {
  set.seed(draw_seed)
  common <- unique_comp <- unique_fib <- numeric(length(testable_programs))
  p_comp <- matrix(NA_real_, nrow = length(testable_programs), ncol = length(discovery_cohorts))
  weights <- numeric(length(discovery_cohorts))
  for (i in seq_along(discovery_cohorts)) {
    cohort <- discovery_cohorts[[i]]
    bl <- blocks[[cohort]]
    Bp <- bl$B[sample.int(nrow(bl$B)), , drop = FALSE]
    fit <- t5_commonality_cohort(Y_cohort[[cohort]], bl$C, bl$A, Bp,
                                 parts = baseline_parts[[cohort]])
    weights[[i]] <- fit$n
    common <- common + fit$common * fit$n
    unique_comp <- unique_comp + fit$unique_composition * fit$n
    unique_fib <- unique_fib + fit$unique_fibrosis * fit$n
    p_comp[, i] <- fit$p_unique_composition
  }
  total <- sum(weights)
  list(common = common / total, unique_composition = unique_comp / total,
       unique_fibrosis = unique_fib / total,
       p_unique_composition = t5_fisher_combine(p_comp)$p_value)
}
coupling_null <- mclapply(seq_len(n_coupling_null),
                          function(b) coupling_null_one(seed + 100000L + b),
                          mc.cores = n_cores)
null_common <- do.call(cbind, lapply(coupling_null, `[[`, "common"))
null_unique_composition <- do.call(cbind, lapply(coupling_null, `[[`, "unique_composition"))
null_unique_fibrosis <- do.call(cbind, lapply(coupling_null, `[[`, "unique_fibrosis"))
null_p_composition <- do.call(cbind, lapply(coupling_null, `[[`, "p_unique_composition"))
rownames(null_common) <- rownames(null_unique_composition) <-
  rownames(null_unique_fibrosis) <- rownames(null_p_composition) <- testable_programs
rm(coupling_null); gc()

common_null_mean <- rowMeans(null_common)
null_fraction <- null_common / observed_agg$fibrosis_incremental
observed_fraction_raw <- t5_coupling_fraction(observed_agg$common,
                                              observed_agg$fibrosis_incremental)
observed_fraction_corrected <- t5_coupling_fraction(observed_agg$common,
                                                    observed_agg$fibrosis_incremental,
                                                    common_null_mean = common_null_mean)
coupling_p <- vapply(seq_along(testable_programs), function(j) {
  (1 + sum(null_fraction[j, ] >= observed_fraction_raw[[j]])) / (ncol(null_fraction) + 1)
}, numeric(1))
coupling_q <- p.adjust(coupling_p, method = "BH")

# ------------------------------------------------------------ null: histology

message("[7/11] Histology-permutation null (", n_histology_null, " draws)")
histology_null_one <- function(draw_seed) {
  set.seed(draw_seed)
  permuted <- permute_histology_within_cohort(composition_meta,
                                              c("fibrosis_stage", "nas_score"))
  p_fib <- matrix(NA_real_, nrow = length(testable_programs), ncol = length(discovery_cohorts))
  for (i in seq_along(discovery_cohorts)) {
    cohort <- discovery_cohorts[[i]]
    rows <- cohort_rows[[cohort]]
    bl <- build_blocks(permuted[rows], permuted[rows, ..pc_columns])
    fit <- t5_commonality_cohort(Y_cohort[[cohort]], bl$C, bl$A, bl$B)
    p_fib[, i] <- fit$p_unique_fibrosis
  }
  t5_fisher_combine(p_fib)$p_value
}
histology_null <- mclapply(seq_len(n_histology_null),
                           function(b) histology_null_one(seed + 200000L + b),
                           mc.cores = n_cores)
null_p_fibrosis <- do.call(cbind, histology_null)
rownames(null_p_fibrosis) <- testable_programs
rm(histology_null); gc()

p_unique_fibrosis <- t5_fisher_combine(
  do.call(cbind, lapply(observed_cohort, `[[`, "p_unique_fibrosis")))$p_value
p_unique_composition <- t5_fisher_combine(
  do.call(cbind, lapply(observed_cohort, `[[`, "p_unique_composition")))$p_value

# ---------------------------------------------------------------- calibration

message("[8/11] Calibrating every view before any count is quoted")
calibrations <- list(
  unique_fibrosis = calibrate_view(
    observed_p = p_unique_fibrosis,
    permuted_p_fn = function(i) null_p_fibrosis[, i],
    n_reps = ncol(null_p_fibrosis), n_tests = length(testable_programs),
    label = "unique_fibrosis_given_composition", observed_p_is_empirical = FALSE),
  unique_composition = calibrate_view(
    observed_p = p_unique_composition,
    permuted_p_fn = function(i) null_p_composition[, i],
    n_reps = ncol(null_p_composition), n_tests = length(testable_programs),
    label = "unique_composition_given_fibrosis", observed_p_is_empirical = FALSE),
  # Leave-one-out: each null draw is scored against the other draws, which is the
  # same procedure the observed statistic went through. Computed by rank rather
  # than by dropping a column, which is the same number and far cheaper:
  # #{k != i: T_k >= T_i} = R - rank_min(T_i).
  coupling = calibrate_view(
    observed_p = coupling_p,
    permuted_p_fn = local({
      loo_p <- t(apply(null_fraction, 1L, function(v) {
        (1 + length(v) - rank(v, ties.method = "min")) / length(v)
      }))
      function(i) loo_p[, i]
    }),
    n_reps = ncol(null_fraction), n_tests = length(testable_programs),
    label = "coupling_above_permuted_composition", observed_p_is_empirical = TRUE)
)
calibration_rows <- rbindlist(lapply(calibrations, calibration_row))
for (name in names(calibrations)) {
  gate <- gate_view(calibrations[[name]], strict = FALSE)
  if (!gate$passed) message("      WITHHELD  ", name, ": ", paste(gate$reasons, collapse = "; "))
}

# The majority-coupled threshold count is not a p-value view, so it gets its own
# measured false-call rate: how many programs cross the same threshold when the
# composition block carries no donor-level information at all.
majority_threshold <- 0.5
null_majority_counts <- colSums(
  sweep(null_fraction, 1L, common_null_mean / observed_agg$fibrosis_incremental, "-") >
    majority_threshold, na.rm = TRUE)
observed_majority <- sum(observed_fraction_corrected > majority_threshold, na.rm = TRUE)
observed_majority_supported <- sum(
  observed_fraction_corrected[match(supported_fibrosis, testable_programs)] > majority_threshold,
  na.rm = TRUE)
null_majority_supported <- colSums(
  sweep(null_fraction[match(supported_fibrosis, testable_programs), , drop = FALSE], 1L,
        (common_null_mean / observed_agg$fibrosis_incremental)[
          match(supported_fibrosis, testable_programs)], "-") > majority_threshold,
  na.rm = TRUE)
threshold_rows <- data.table(
  view = c("majority_coupled_all_testable", "majority_coupled_fibrosis_supported"),
  observed_calls = c(observed_majority, observed_majority_supported),
  null_call_mean = c(mean(null_majority_counts), mean(null_majority_supported)),
  null_call_max = c(max(null_majority_counts), max(null_majority_supported)),
  n_permutations = ncol(null_fraction),
  n_tests = c(length(testable_programs), length(supported_fibrosis)),
  resolution_floor = 0, bh_threshold_needed = NA_real_,
  calibrated = c(mean(null_majority_counts), mean(null_majority_supported)) <=
    DEFAULT_CALIBRATION_TOLERANCE,
  resolvable = TRUE
)
threshold_rows[, reportable := calibrated & resolvable]
calibration_rows <- rbind(calibration_rows, threshold_rows, fill = TRUE)

# ------------------------------------------------------------------ bootstrap

message("[9/11] Bootstrapping donors within cohort (", n_boot, " resamples)")
boot_one <- function(draw_seed) {
  set.seed(draw_seed)
  common <- unique_comp <- unique_fib <- increment <- numeric(length(testable_programs))
  total <- 0
  for (cohort in discovery_cohorts) {
    bl <- blocks[[cohort]]
    n <- nrow(bl$B)
    idx <- sample.int(n, n, replace = TRUE)
    fit <- t5_commonality_cohort(Y_cohort[[cohort]][idx, , drop = FALSE],
                                 bl$C[idx, , drop = FALSE],
                                 bl$A[idx, , drop = FALSE],
                                 bl$B[idx, , drop = FALSE])
    if (is.null(fit)) return(NULL)
    total <- total + n
    common <- common + fit$common * n
    unique_comp <- unique_comp + fit$unique_composition * n
    unique_fib <- unique_fib + fit$unique_fibrosis * n
    increment <- increment + fit$fibrosis_incremental * n
  }
  cbind(common / total, unique_comp / total, unique_fib / total, increment / total)
}
boot <- mclapply(seq_len(n_boot), function(b) boot_one(seed + 300000L + b),
                 mc.cores = n_cores)
boot_ok <- !vapply(boot, is.null, logical(1))
message("      usable resamples: ", sum(boot_ok), " / ", n_boot)
boot_common <- do.call(cbind, lapply(boot[boot_ok], function(x) x[, 1L]))
boot_unique_composition <- do.call(cbind, lapply(boot[boot_ok], function(x) x[, 2L]))
boot_unique_fibrosis <- do.call(cbind, lapply(boot[boot_ok], function(x) x[, 3L]))
boot_increment <- do.call(cbind, lapply(boot[boot_ok], function(x) x[, 4L]))
rm(boot); gc()
boot_fraction_raw <- boot_common / boot_increment
boot_fraction_corrected <- (boot_common - common_null_mean) / boot_increment
ci <- function(M) {
  data.table(lower = apply(M, 1L, stats::quantile, 0.025, na.rm = TRUE),
             upper = apply(M, 1L, stats::quantile, 0.975, na.rm = TRUE))
}

# ---------------------------------------------------------------- control arms

# Same donors, same cohorts, same code path, same number of extra parameters.
message("[10/11] Same-cohort control arms (", n_control_draws, " draws each)")
arm_summary <- function(effects, reference) {
  d <- merge(effects[axis == "fibrosis", .(feature_id, beta_meta, se_meta,
                                           statistic_meta, q_value)],
             reference[axis == "fibrosis", .(feature_id, beta_reference = beta_meta,
                                             se_reference = se_meta,
                                             t_reference = statistic_meta)],
             by = "feature_id")
  s <- d[feature_id %in% supported_fibrosis]
  data.table(
    n_supported = sum(d$q_value < 0.05, na.rm = TRUE),
    min_q = min(d$q_value, na.rm = TRUE),
    median_abs_beta_supported = stats::median(abs(s$beta_meta), na.rm = TRUE),
    median_beta_ratio_supported = stats::median(abs(s$beta_meta) / abs(s$beta_reference),
                                                na.rm = TRUE),
    median_se_ratio_supported = stats::median(s$se_meta / s$se_reference, na.rm = TRUE),
    median_abs_t_supported = stats::median(abs(s$statistic_meta), na.rm = TRUE),
    n_sign_preserved = sum(sign(d$beta_meta) == sign(d$beta_reference), na.rm = TRUE)
  )
}
fit_arm <- function(md) {
  meta_analyze(fit_feature_models(S, md, discovery_cohorts, "composition_pc"))
}
# BH family matches the discovery release exactly: one model kind, both axes,
# 226 tests. Adjusting within axis instead would silently change every q value
# relative to the sealed arm this system is compared against.
add_q <- function(effects) {
  effects[, q_value := p.adjust(p_value_meta, method = "BH", n = .N)]
  effects[]
}
reference_effects <- add_q(meta_analyze(fit_feature_models(S, composition_meta,
                                                           discovery_cohorts, "primary")))
real_effects <- add_q(fit_arm(composition_meta))

control_draw <- function(draw_seed, stratified) {
  set.seed(draw_seed)
  md <- copy(composition_meta)
  P <- as.matrix(md[, ..pc_columns])
  P <- t5_permute_block_within(P, md$dataset,
                               strata = if (stratified) md$fibrosis_stage else NULL)
  for (j in seq_along(pc_columns)) md[[pc_columns[[j]]]] <- P[, j]
  arm_summary(add_q(fit_arm(md)), reference_effects)
}
control_permuted <- rbindlist(mclapply(seq_len(n_control_draws),
  function(b) control_draw(seed + 400000L + b, FALSE), mc.cores = n_cores))
control_stage <- rbindlist(mclapply(seq_len(n_control_draws),
  function(b) control_draw(seed + 500000L + b, TRUE), mc.cores = n_cores))

control_table <- rbindlist(list(
  cbind(arm = "primary_no_composition", draw = "observed",
        arm_summary(reference_effects, reference_effects)),
  cbind(arm = "composition_pc_real", draw = "observed",
        arm_summary(real_effects, reference_effects)),
  cbind(arm = "composition_permuted", draw = "control_mean",
        control_permuted[, lapply(.SD, mean, na.rm = TRUE)]),
  cbind(arm = "composition_permuted", draw = "control_p2.5",
        control_permuted[, lapply(.SD, stats::quantile, 0.025, na.rm = TRUE)]),
  cbind(arm = "composition_permuted", draw = "control_p97.5",
        control_permuted[, lapply(.SD, stats::quantile, 0.975, na.rm = TRUE)]),
  cbind(arm = "composition_stage_matched", draw = "control_mean",
        control_stage[, lapply(.SD, mean, na.rm = TRUE)]),
  cbind(arm = "composition_stage_matched", draw = "control_p2.5",
        control_stage[, lapply(.SD, stats::quantile, 0.025, na.rm = TRUE)]),
  cbind(arm = "composition_stage_matched", draw = "control_p97.5",
        control_stage[, lapply(.SD, stats::quantile, 0.975, na.rm = TRUE)])
), use.names = TRUE, fill = TRUE)

# Minimum detectable effect for the negatives, on the same Knapp-Hartung
# reference the arm itself used.
real_effects[, mde := minimum_detectable_effect(se_meta, n_cohorts)]

# And the same thing in the units this system actually reports. A program with no
# supported unique-composition term is only interpretable next to the smallest
# partial R-squared the design could have resolved at that cohort's n.
min_detectable_partial_r2 <- function(n, df1, df_error, alpha, power = 0.80) {
  critical <- stats::qf(1 - alpha, df1, df_error)
  shortfall <- function(f2) {
    stats::pf(critical, df1, df_error, ncp = f2 * n, lower.tail = FALSE) - power
  }
  upper <- 1
  while (shortfall(upper) < 0 && upper < 1e6) upper <- upper * 2
  f2 <- stats::uniroot(shortfall, c(1e-10, upper), tol = 1e-12)$root
  f2 / (1 + f2)
}
mde_r2 <- rbindlist(lapply(discovery_cohorts, function(cohort) {
  fit <- observed_cohort[[cohort]]
  rbindlist(lapply(list(list("fibrosis", 1L), list("composition_block", 3L)), function(spec) {
    data.table(
      cohort = cohort, n_donors = fit$n, block = spec[[1L]], df = spec[[2L]],
      df_error = fit$df_error,
      min_detectable_partial_r2_alpha05 = min_detectable_partial_r2(
        fit$n, spec[[2L]], fit$df_error, 0.05),
      min_detectable_partial_r2_bonferroni = min_detectable_partial_r2(
        fit$n, spec[[2L]], fit$df_error, 0.05 / length(testable_programs)))
  }))
}))
mde_r2 <- rbind(mde_r2, mde_r2[, .(cohort = "donor_weighted_mean", n_donors = sum(n_donors),
                                   df_error = NA_integer_,
                                   min_detectable_partial_r2_alpha05 =
                                     stats::weighted.mean(min_detectable_partial_r2_alpha05, n_donors),
                                   min_detectable_partial_r2_bonferroni =
                                     stats::weighted.mean(min_detectable_partial_r2_bonferroni, n_donors)),
                              by = .(block, df)], use.names = TRUE, fill = TRUE)

# ------------------------------------------------------ construction ceiling

# Evidence for the ceiling, not a correction to it. If program scores track their
# OWN lineage's estimated abundance more than other programs do, part of the
# shared variance is definitional. Correlations are computed within cohort and
# combined on the Fisher z scale, because raw proportions are not comparable
# across cohorts.
lineage_of <- c(cholangiocytes = "Cholangiocytes", fibroblasts = "Fibroblasts",
                hepatocytes = "Hepatocytes", macrophages = "Macrophages",
                tcells = "T cells")
prop <- composition[match(composition_meta$sample_id, sample_id)]
own_lineage <- registry[match(testable_programs, program_uid),
                        lineage_of[cell_type]]
within_cohort_rho <- function(x, y) {
  z <- 0; w <- 0
  for (cohort in discovery_cohorts) {
    rows <- cohort_rows[[cohort]]
    r <- suppressWarnings(stats::cor(x[rows], y[rows], method = "spearman"))
    if (!is.finite(r)) next
    r <- max(min(r, 1 - 1e-12), -1 + 1e-12)
    z <- z + atanh(r) * (length(rows) - 3L)
    w <- w + (length(rows) - 3L)
  }
  if (w == 0) return(NA_real_)
  tanh(z / w)
}
score_by_sample <- S[, composition_meta$sample_id, drop = FALSE]
ceiling_table <- rbindlist(lapply(seq_along(testable_programs), function(j) {
  program <- testable_programs[[j]]
  own <- own_lineage[[j]]
  data.table(
    feature_id = program,
    cell_type = registry[program_uid == program, cell_type],
    own_lineage = own,
    rho_own_lineage_within_cohort = within_cohort_rho(score_by_sample[j, ], prop[[own]]),
    rho_own_lineage_pooled = suppressWarnings(
      stats::cor(score_by_sample[j, ], prop[[own]], method = "spearman")),
    rho_hepatocytes_within_cohort = within_cohort_rho(score_by_sample[j, ],
                                                      prop[["Hepatocytes"]]),
    rho_macrophages_within_cohort = within_cohort_rho(score_by_sample[j, ],
                                                      prop[["Macrophages"]])
  )
}))

# ------------------------------------------------------------------- assembly

message("[11/11] Assembling, gating and sealing")
program_table <- data.table(
  feature_id = testable_programs,
  cell_type = registry[match(testable_programs, program_uid), cell_type],
  module_name = registry[match(testable_programs, program_uid), module_name],
  robust_display = registry[match(testable_programs, program_uid), robust_display],
  fibrosis_supported_linear = testable_programs %in% supported_fibrosis,
  r2_covariates = observed_agg$r2_covariates,
  r2_full = observed_agg$r2_full,
  unique_fibrosis = observed_agg$unique_fibrosis,
  unique_composition = observed_agg$unique_composition,
  common = observed_agg$common,
  fibrosis_incremental = observed_agg$fibrosis_incremental,
  common_null_mean = common_null_mean,
  common_corrected = observed_agg$common - common_null_mean,
  unique_composition_null_mean = rowMeans(null_unique_composition),
  unique_fibrosis_null_mean = rowMeans(null_unique_fibrosis),
  coupling_fraction_raw = observed_fraction_raw,
  coupling_fraction_corrected = observed_fraction_corrected,
  suppression = observed_agg$common < 0,
  p_unique_fibrosis = p_unique_fibrosis,
  q_unique_fibrosis = p.adjust(p_unique_fibrosis, method = "BH"),
  p_unique_composition = p_unique_composition,
  q_unique_composition = p.adjust(p_unique_composition, method = "BH"),
  coupling_empirical_p = coupling_p,
  coupling_q = coupling_q
)
program_table <- cbind(
  program_table,
  setnames(ci(boot_common), c("common_ci_lower", "common_ci_upper")),
  setnames(ci(boot_unique_fibrosis), c("unique_fibrosis_ci_lower", "unique_fibrosis_ci_upper")),
  setnames(ci(boot_unique_composition), c("unique_composition_ci_lower",
                                          "unique_composition_ci_upper")),
  setnames(ci(boot_fraction_raw), c("coupling_fraction_raw_ci_lower",
                                    "coupling_fraction_raw_ci_upper")),
  setnames(ci(boot_fraction_corrected), c("coupling_fraction_corrected_ci_lower",
                                          "coupling_fraction_corrected_ci_upper"))
)
program_table <- merge(program_table, ceiling_table[, !c("cell_type")], by = "feature_id",
                       sort = FALSE)
program_table <- merge(
  program_table,
  real_effects[axis == "fibrosis", .(feature_id, beta_composition_adjusted = beta_meta,
                                     se_composition_adjusted = se_meta,
                                     q_composition_adjusted = q_value,
                                     mde_composition_adjusted = mde)],
  by = "feature_id", sort = FALSE)
program_table <- merge(
  program_table,
  reference_effects[axis == "fibrosis", .(feature_id, beta_unadjusted = beta_meta,
                                          se_unadjusted = se_meta)],
  by = "feature_id", sort = FALSE)
program_table <- cbind(program_table, decompose_effect(
  program_table$beta_unadjusted, program_table$beta_composition_adjusted)[
    , .(beta_retained_fraction = retained_fraction, beta_sign_preserved = sign_preserved)])
program_table[, coupling_class := fcase(
  common < 0, "suppression",
  coupling_q < 0.05 & coupling_fraction_corrected > 0.5, "majority_coupled",
  coupling_q < 0.05, "coupled_minority",
  default = "no_coupling_above_permuted_floor"
)]

per_cohort_table <- rbindlist(lapply(discovery_cohorts, function(cohort) {
  fit <- observed_cohort[[cohort]]
  data.table(cohort = cohort, n_donors = fit$n, feature_id = testable_programs,
             r2_covariates = fit$r2_covariates, r2_full = fit$r2_full,
             unique_fibrosis = fit$unique_fibrosis,
             unique_composition = fit$unique_composition, common = fit$common,
             fibrosis_incremental = fit$fibrosis_incremental,
             p_unique_fibrosis = fit$p_unique_fibrosis,
             p_unique_composition = fit$p_unique_composition)
}))

boot_counts <- data.table(
  n_majority_coupled_corrected = colSums(boot_fraction_corrected > majority_threshold,
                                         na.rm = TRUE),
  n_suppression = colSums(boot_common < 0, na.rm = TRUE)
)

# Every count that will be quoted, next to the null rate of the view it came from.
counts <- data.table(
  quantity = c("n_testable_programs", "n_fibrosis_supported_linear",
               "n_unique_fibrosis_supported", "n_unique_composition_supported",
               "n_coupling_above_permuted_floor", "n_majority_coupled_all_testable",
               "n_majority_coupled_fibrosis_supported", "n_suppression"),
  value = c(length(testable_programs), length(supported_fibrosis),
            sum(program_table$q_unique_fibrosis < 0.05, na.rm = TRUE),
            sum(program_table$q_unique_composition < 0.05, na.rm = TRUE),
            sum(program_table$coupling_q < 0.05, na.rm = TRUE),
            observed_majority, observed_majority_supported,
            sum(program_table$common < 0, na.rm = TRUE)),
  view = c(NA, NA, "unique_fibrosis_given_composition",
           "unique_composition_given_fibrosis", "coupling_above_permuted_composition",
           "majority_coupled_all_testable", "majority_coupled_fibrosis_supported",
           NA)
)
counts <- merge(counts, calibration_rows[, .(view, null_call_mean, null_call_max, reportable)],
                by = "view", all.x = TRUE, sort = FALSE)
counts[is.na(view), reportable := TRUE]
counts[, reportable := as.logical(reportable)]
counts[reportable == FALSE, `:=`(value = NA_real_, withheld_reason = "failed_calibration_gate")]

# Non-vacuous: every count that still carries a value and came from a view must
# name a view whose measured null false-call rate passed the gate.
quoted <- counts[!is.na(value) & !is.na(view)]
assert_true(all(quoted$view %in% calibration_rows[reportable == TRUE, view]),
            paste0("A count is quoted from a view that failed the calibration gate: ",
                   paste(setdiff(quoted$view, calibration_rows[reportable == TRUE, view]),
                         collapse = ", ")))
assert_counts_calibrated(quoted, calibration_rows[view %in% quoted$view])
assert_negatives_have_mde(
  program_table[, .(support = fcase(q_unique_composition < 0.05, "supported",
                                    default = "unsupported"),
                    mde = mde_composition_adjusted)],
  "support", "mde")

# ---------------------------------------------------------------------- write

tmp <- atomic_dir(system_root)
on.exit(if (dir.exists(tmp)) unlink(tmp, recursive = TRUE), add = TRUE)
write_tsv(program_table, file.path(tmp, "program_composition_decomposition.tsv"))
write_tsv(per_cohort_table, file.path(tmp, "per_cohort_decomposition.tsv"))
write_tsv(calibration_rows, file.path(tmp, "calibration.tsv"))
write_tsv(counts, file.path(tmp, "counts.tsv"))
write_tsv(control_table, file.path(tmp, "control_arm_comparison.tsv"))
write_tsv(control_permuted[, arm := "composition_permuted"][],
          file.path(tmp, "control_draws_permuted.tsv"))
write_tsv(control_stage[, arm := "composition_stage_matched"][],
          file.path(tmp, "control_draws_stage_matched.tsv"))
write_tsv(ceiling_table, file.path(tmp, "construction_ceiling_probe.tsv"))
write_tsv(boot_counts, file.path(tmp, "bootstrap_counts.tsv"))
write_tsv(mde_r2, file.path(tmp, "minimum_detectable_partial_r2.tsv"))
saveRDS(list(null_common = null_common, null_fraction = null_fraction,
             null_unique_composition = null_unique_composition,
             null_p_composition = null_p_composition,
             null_p_fibrosis = null_p_fibrosis),
        file.path(tmp, "permutation_nulls.rds"), compress = "xz")
writeLines(capture.output(sessionInfo()), file.path(tmp, "sessionInfo.txt"))

inputs <- c(nonholdout_dge, nonholdout_meta, nonholdout_composition, program_registry,
            program_membership, gene_annotation, composition_accepted_seal,
            file.path(discovery_root, "cohort_effects.tsv"),
            file.path(discovery_root, "program_axis_map.tsv"))
write_tsv(data.table(path = inputs, sha256 = vapply(inputs, sha256_file, character(1))),
          file.path(tmp, "input_manifest.tsv"))

estimand_text <- paste(
  "Coupling between a frozen program's fibrosis-associated variance and BayesPrism",
  "cell-type composition estimated from the same expression matrix. The shared term",
  "is an upper bound on composition attribution and is never attribution: the",
  "proportions are deconvolved from the same measurements as the scores, so a",
  "lineage-marker program shares variance with its own estimated abundance partly by",
  "construction.")
assert_language(c(estimand_text, program_table$coupling_class, counts$quantity,
                  control_table$arm, calibration_rows$view))

jsonlite::write_json(list(
  state = "T5_COMPOSITION_COUPLING_COMPLETE_HOLDOUT_UNOPENED",
  release_id = release_id, workstream_id = workstream_id, system = "T5_composition",
  estimand = estimand_text,
  composition_method = seal$accepted_method,
  composition_specification = "clr_principal_components",
  n_donors = nrow(composition_meta), n_cohorts = length(discovery_cohorts),
  n_testable_programs = length(testable_programs),
  n_coupling_permutations = n_coupling_null,
  n_histology_permutations = n_histology_null,
  n_bootstrap = sum(boot_ok), n_control_draws = n_control_draws,
  additivity_residual = additivity_error,
  holdout_accessed = FALSE
), file.path(tmp, "T5_READY.json"), pretty = TRUE, auto_unbox = TRUE)

artifacts <- setdiff(list.files(tmp, recursive = TRUE), "artifact_manifest.tsv")
write_tsv(data.table(artifact = artifacts,
                     sha256 = vapply(file.path(tmp, artifacts), sha256_file, character(1)),
                     size_bytes = file.info(file.path(tmp, artifacts))$size),
          file.path(tmp, "artifact_manifest.tsv"))
publish_dir(tmp, system_root)
Sys.chmod(list.files(system_root, full.names = TRUE), mode = "0440")

cat("\n=== T5 CALIBRATION ===\n"); print(calibration_rows)
cat("\n=== T5 COUNTS ===\n"); print(counts)
cat("\n=== T5 CONTROL ARMS ===\n"); print(control_table)
cat("\n=== T5 COUPLING CLASSES ===\n"); print(program_table[, .N, by = coupling_class])
cat("\n=== CEILING PROBE (median |rho| to own lineage) ===\n")
print(program_table[, .(median_abs_rho = stats::median(abs(rho_own_lineage_within_cohort),
                                                       na.rm = TRUE), .N),
                    by = fibrosis_supported_linear])
cat("T5_COMPLETE\t", system_root, "\n", sep = "")
