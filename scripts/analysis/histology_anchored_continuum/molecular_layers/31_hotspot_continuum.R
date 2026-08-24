#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(splines)
})

.hotspot_args <- commandArgs(trailingOnly = FALSE)
.hotspot_file <- sub("^--file=", "", .hotspot_args[grepl("^--file=", .hotspot_args)])
.hotspot_dir <- if (length(.hotspot_file)) dirname(normalizePath(.hotspot_file)) else getwd()
.hotspot_dir_override <- Sys.getenv("HAC_ML_SCRIPT_DIR", unset = "")
if (nzchar(.hotspot_dir_override)) .hotspot_dir <- normalizePath(.hotspot_dir_override)
source(file.path(.hotspot_dir, "lib_molecular_layers.R"))

hotspot_fit_continuum <- function(data, stage_column) {
  d <- copy(data)[
    is.finite(outcome_z) & is.finite(axis_raw) &
      !is.na(get(stage_column)) & !is.na(inferred_sex)
  ]
  empty <- data.table(
    estimable = FALSE, failure_reason = "insufficient_complete_design",
    n_participants = nrow(d), beta = NA_real_, se = NA_real_,
    ci_low = NA_real_, ci_high = NA_real_, p_value = NA_real_,
    nonlinear_p_value = NA_real_
  )
  if (nrow(d) < 20L || uniqueN(d[[stage_column]]) < 2L ||
      uniqueN(d$inferred_sex) < 2L) return(empty)
  d[, `:=`(
    axis_z = ml_standardize(axis_raw),
    stage_factor = factor(get(stage_column)),
    sex_factor = factor(inferred_sex)
  )]
  linear <- tryCatch(
    lm(outcome_z ~ axis_z + stage_factor + sex_factor, data = d),
    error = function(e) e
  )
  if (inherits(linear, "error") ||
      !"axis_z" %in% rownames(coef(summary(linear)))) {
    empty[, failure_reason := if (inherits(linear, "error")) {
      conditionMessage(linear)
    } else {
      "axis_not_estimable"
    }]
    return(empty)
  }
  coefficient <- coef(summary(linear))["axis_z", ]
  nonlinear <- tryCatch(
    lm(outcome_z ~ ns(axis_z, df = 3) + stage_factor + sex_factor, data = d),
    error = function(e) NULL
  )
  nonlinear_p <- NA_real_
  if (!is.null(nonlinear)) {
    comparison <- tryCatch(anova(linear, nonlinear), error = function(e) NULL)
    if (!is.null(comparison) && nrow(comparison) == 2L) {
      nonlinear_p <- comparison$`Pr(>F)`[[2L]]
    }
  }
  data.table(
    estimable = TRUE, failure_reason = NA_character_, n_participants = nrow(d),
    beta = unname(coefficient[["Estimate"]]),
    se = unname(coefficient[["Std. Error"]]),
    ci_low = unname(coefficient[["Estimate"]] -
                      qt(0.975, df.residual(linear)) * coefficient[["Std. Error"]]),
    ci_high = unname(coefficient[["Estimate"]] +
                       qt(0.975, df.residual(linear)) * coefficient[["Std. Error"]]),
    p_value = unname(coefficient[["Pr(>|t|)"]]),
    nonlinear_p_value = nonlinear_p
  )
}

hotspot_adjusted_shape <- function(data, stage_column) {
  grid <- seq(0, 1, length.out = 101L)
  d <- copy(data)[
    is.finite(outcome_z) & is.finite(axis_raw) &
      !is.na(get(stage_column)) & !is.na(inferred_sex)
  ]
  blank <- data.table(percentile = grid, predicted = NA_real_, predicted_z = NA_real_)
  if (nrow(d) < 20L || uniqueN(d[[stage_column]]) < 2L ||
      uniqueN(d$inferred_sex) < 2L) return(blank)
  d[, `:=`(
    axis_z = ml_standardize(axis_raw),
    stage_factor = factor(get(stage_column)),
    sex_factor = factor(inferred_sex)
  )]
  fit <- tryCatch(
    lm(outcome_z ~ ns(axis_z, df = 3) + stage_factor + sex_factor, data = d),
    error = function(e) NULL
  )
  if (is.null(fit)) return(blank)
  modal_stage <- names(which.max(table(d$stage_factor)))[[1L]]
  modal_sex <- names(which.max(table(d$sex_factor)))[[1L]]
  new_data <- data.table(
    axis_z = as.numeric(quantile(d$axis_z, grid, names = FALSE, type = 8)),
    stage_factor = factor(modal_stage, levels = levels(d$stage_factor)),
    sex_factor = factor(modal_sex, levels = levels(d$sex_factor))
  )
  prediction <- tryCatch(as.numeric(predict(fit, newdata = new_data)),
                         error = function(e) rep(NA_real_, length(grid)))
  data.table(
    percentile = grid,
    predicted = prediction,
    predicted_z = ml_standardize(prediction)
  )
}

hotspot_within_stage_permutation <- function(data, replicates, seed) {
  d <- copy(data)[
    is.finite(outcome_z) & is.finite(axis_raw) &
      !is.na(fibrosis_stage) & !is.na(inferred_sex)
  ]
  if (nrow(d) < 20L || uniqueN(d$fibrosis_stage) < 2L ||
      uniqueN(d$inferred_sex) < 2L) {
    return(data.table(
      estimable = FALSE, n_participants = nrow(d), residual_pearson = NA_real_,
      permutation_p_value = NA_real_, permutation_replicates = replicates,
      n_finite_permutations = 0L
    ))
  }
  d[, `:=`(
    axis_z = ml_standardize(axis_raw),
    stage_factor = factor(fibrosis_stage),
    sex_factor = factor(inferred_sex)
  )]
  axis_residual <- residuals(lm(axis_z ~ stage_factor + sex_factor, data = d))
  outcome_residual <- residuals(lm(outcome_z ~ stage_factor + sex_factor, data = d))
  observed <- suppressWarnings(cor(axis_residual, outcome_residual))
  strata <- interaction(d$stage_factor, d$sex_factor, drop = TRUE)
  groups <- split(seq_along(axis_residual), strata)
  set.seed(seed)
  null <- replicate(replicates, {
    permuted <- axis_residual
    for (group in groups) {
      permuted[group] <- sample(axis_residual[group], length(group), replace = FALSE)
    }
    suppressWarnings(cor(permuted, outcome_residual))
  })
  n_finite <- sum(is.finite(null))
  empirical <- (1 + sum(abs(null) >= abs(observed), na.rm = TRUE)) / (1 + n_finite)
  data.table(
    estimable = is.finite(observed), n_participants = nrow(d),
    residual_pearson = observed, permutation_p_value = empirical,
    permutation_replicates = replicates, n_finite_permutations = n_finite
  )
}

hotspot_run <- function() {
  contract <- ml_read_contract()
  output_dir <- file.path(ml_out_root(must_exist = TRUE), "programs", "hotspot")
  ml_assert(!dir.exists(output_dir), paste0("Refusing to overwrite: ", output_dir))
  ml_ensure_dir(output_dir)
  source_root <- ml_source_root(contract)

  registry <- fread(ml_resolve(contract$program_registry))
  membership_raw <- fread(ml_resolve(contract$program_membership))
  source_scores <- fread(file.path(source_root, "programs", "program_scores.tsv.gz"))
  source_testability <- fread(file.path(source_root, "programs", "program_testability.tsv"))
  axes <- ml_load_axes(contract)
  ml_assert(nrow(registry) == contract$program_family_size &&
              uniqueN(registry$program_uid) == contract$program_family_size,
            "Frozen Hotspot registry is not the complete 117-program family")
  ml_assert(uniqueN(source_scores$program_uid) == contract$program_family_size,
            "Frozen signature-excluded score family is incomplete")
  ml_assert(uniqueN(source_testability$program_uid) == contract$program_family_size,
            "Frozen Hotspot testability family is incomplete")
  ml_assert(!anyDuplicated(source_scores[, .(sample_id, program_uid)]),
            "Hotspot scores contain duplicate participant-program rows")

  signature_symbols <- ml_signature_symbols(contract)
  membership <- membership_raw[
    !is.na(mapped_symbol) & mapped_symbol != "",
    .(original_l1_weight = sum(as.numeric(original_l1_weight))),
    by = .(program_uid, gene_symbol = toupper(mapped_symbol))
  ]
  total <- membership[, .(total_l1_weight = sum(original_l1_weight)), by = program_uid]
  membership[, excluded_signature_gene := gene_symbol %in% toupper(signature_symbols)]
  retained <- membership[excluded_signature_gene == FALSE]
  retained <- merge(retained, total, by = "program_uid", all.x = TRUE)
  retained[, `:=`(
    retained_program_l1_weight = sum(original_l1_weight),
    n_retained_genes = uniqueN(gene_symbol)
  ), by = program_uid]
  retained[, `:=`(
    retained_l1_fraction = retained_program_l1_weight / total_l1_weight,
    retained_relative_weight = original_l1_weight / retained_program_l1_weight
  )]
  retained[, testable_by_signature_exclusion :=
             n_retained_genes >= contract$minimum_set_genes &
             retained_l1_fraction >= contract$minimum_retained_fraction]
  retained <- merge(
    retained,
    registry[, .(program_uid, cell_type, module, module_name, robust_display)],
    by = "program_uid", all.x = TRUE
  )
  ml_write_tsv_once(retained,
                    file.path(output_dir, "hotspot_signature_excluded_membership.tsv.gz"))

  evaluation <- as.character(contract$evaluation_cohorts)
  testability <- source_testability[dataset %in% evaluation]
  ml_assert(nrow(testability) == length(evaluation) * contract$program_family_size,
            "Evaluation-cohort Hotspot testability rows are incomplete")
  activated <- testability[grepl("Activated stellate.*PDGFRA", module_name)]
  ml_assert(nrow(activated) == length(evaluation),
            "Activated stellate (PDGFRA) testability row is absent")
  ml_assert(all(!activated$testable) &&
              all(abs(activated$retained_l1_fraction - 0.720) < 0.002),
            "Activated stellate (PDGFRA) must remain untestable at 72.0% retained L1")
  testability[, `:=`(
    family_size = contract$program_family_size,
    minimum_genes = contract$minimum_set_genes,
    minimum_retained_l1_fraction = contract$minimum_retained_fraction,
    signature_exclusion_frozen = TRUE
  )]
  ml_write_tsv_once(testability, file.path(output_dir, "hotspot_testability.tsv"))

  scores <- merge(
    source_scores,
    registry[, .(program_uid, cell_type, module, module_name, robust_display)],
    by = "program_uid", all.x = TRUE
  )
  setnames(scores, "program_score", "score_raw")
  scores[, outcome_z := score_raw]
  scores[, `:=`(
    molecular_layer = "frozen_hotspot_program",
    signature_excluded = TRUE,
    same_expression_substrate = TRUE
  )]
  ml_write_tsv_once(scores, file.path(output_dir, "hotspot_participant_scores.tsv.gz"))

  model_rows <- list()
  nas_rows <- list()
  shape_rows <- list()
  registry_ids <- registry$program_uid
  for (axis_name in as.character(contract$co_primary_axes)) {
    for (cohort in evaluation) {
      joined <- merge(
        scores[dataset == cohort],
        axes[axis_id == axis_name & dataset == cohort,
             .(sample_id, axis_raw, axis_percentile)],
        by = "sample_id", all = FALSE
      )
      for (program in registry_ids) {
        block <- joined[program_uid == program]
        fit <- hotspot_fit_continuum(block, "fibrosis_stage")
        fit[, `:=`(
          program_uid = program, axis_id = axis_name, dataset = cohort,
          adjustment = "factor(fibrosis_stage)+inferred_sex"
        )]
        model_rows[[length(model_rows) + 1L]] <- fit
        shape <- hotspot_adjusted_shape(block, "fibrosis_stage")
        shape[, `:=`(
          program_uid = program, axis_id = axis_name, dataset = cohort,
          adjustment = "factor(fibrosis_stage)+inferred_sex"
        )]
        shape_rows[[length(shape_rows) + 1L]] <- shape
        if (cohort == "GSE162694") {
          nas <- hotspot_fit_continuum(block, "nas_score")
          nas[, `:=`(
            program_uid = program, axis_id = axis_name, dataset = cohort,
            adjustment = "factor(nas_score)+inferred_sex"
          )]
          nas_rows[[length(nas_rows) + 1L]] <- nas
        }
      }
    }
  }
  models <- rbindlist(model_rows, fill = TRUE)
  models[, `:=`(
    q_value = ml_complete_bh(p_value, contract$program_family_size),
    nonlinear_q_value = ml_complete_bh(
      nonlinear_p_value, contract$program_family_size
    )
  ), by = .(axis_id, dataset)]
  models <- merge(models,
                  registry[, .(program_uid, cell_type, module, module_name, robust_display)],
                  by = "program_uid", all.x = TRUE)
  models <- merge(
    models,
    testability[, .(program_uid, dataset, testable, retained_l1_fraction,
                    n_observed_genes)],
    by = c("program_uid", "dataset"), all.x = TRUE
  )
  models[testable != TRUE | is.na(testable), `:=`(
    estimable = FALSE,
    failure_reason = "signature_exclusion_testability_gate_failed",
    beta = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
    p_value = NA_real_, q_value = NA_real_, nonlinear_p_value = NA_real_,
    nonlinear_q_value = NA_real_
  )]
  ml_write_tsv_once(models, file.path(output_dir, "hotspot_fibrosis_models.tsv"))

  nas <- rbindlist(nas_rows, fill = TRUE)
  nas[, `:=`(
    q_value = ml_complete_bh(p_value, contract$program_family_size),
    nonlinear_q_value = ml_complete_bh(
      nonlinear_p_value, contract$program_family_size
    )
  ), by = axis_id]
  nas <- merge(nas,
               registry[, .(program_uid, cell_type, module, module_name, robust_display)],
               by = "program_uid", all.x = TRUE)
  ml_write_tsv_once(nas, file.path(output_dir, "hotspot_nas_sensitivity.tsv"))

  shapes <- rbindlist(shape_rows, fill = TRUE)
  shapes <- merge(shapes, registry[, .(program_uid, module_name, robust_display)],
                  by = "program_uid", all.x = TRUE)
  ml_write_tsv_once(shapes,
                    file.path(output_dir, "hotspot_adjusted_spline_shapes.tsv.gz"))
  shape_concordance <- shapes[, {
    wide <- dcast(.SD, percentile ~ dataset, value.var = "predicted_z")
    rho <- NA_real_
    if (all(evaluation %in% names(wide))) {
      x <- wide[[evaluation[[1L]]]]
      y <- wide[[evaluation[[2L]]]]
      ok <- is.finite(x) & is.finite(y)
      if (sum(ok) >= 3L) rho <- suppressWarnings(cor(x[ok], y[ok], method = "spearman"))
    }
    list(shape_spearman_rho = rho)
  }, by = .(program_uid, axis_id)]
  nonlinear_gate <- models[, .(
    nonlinear_bh_both_cohorts = .N == length(evaluation) &&
      all(is.finite(nonlinear_q_value) & nonlinear_q_value < 0.05)
  ), by = .(program_uid, axis_id)]
  shape_concordance <- merge(shape_concordance, nonlinear_gate,
                             by = c("program_uid", "axis_id"))
  shape_concordance[, timing_gate_pass := nonlinear_bh_both_cohorts &
                      is.finite(shape_spearman_rho) & shape_spearman_rho >= 0.70]
  ml_write_tsv_once(shape_concordance,
                    file.path(output_dir, "hotspot_nonlinear_timing_gates.tsv"))

  meta_rows <- list()
  for (axis_name in as.character(contract$co_primary_axes)) {
    for (program in registry_ids) {
      block <- models[
        axis_id == axis_name & program_uid == program & dataset %in% evaluation &
          estimable == TRUE & testable == TRUE
      ]
      fixed_meta <- ml_fixed_meta(block$beta, block$se)
      random_meta <- ml_random_meta(block$beta, block$se)
      row <- fixed_meta
      setnames(row, c("beta", "se", "z", "p_value", "ci_low", "ci_high"),
               paste0("fixed_", c("beta", "se", "z", "p_value", "ci_low", "ci_high")))
      row[, `:=`(
        random_beta = random_meta$beta, random_se = random_meta$se,
        random_p_value = random_meta$p_value, random_ci_low = random_meta$ci_low,
        random_ci_high = random_meta$ci_high, random_tau2 = random_meta$tau2,
        program_uid = program, axis_id = axis_name,
        direction_concordant = nrow(block) == length(evaluation) &&
          length(unique(sign(block$beta))) == 1L,
        cohort_betas = paste(block$dataset, signif(block$beta, 6),
                             sep = "=", collapse = ";")
      )]
      meta_rows[[length(meta_rows) + 1L]] <- row
    }
  }
  meta <- rbindlist(meta_rows, fill = TRUE)
  meta[, fixed_q_value := ml_complete_bh(fixed_p_value, contract$program_family_size),
       by = axis_id]
  meta <- merge(meta,
                registry[, .(program_uid, cell_type, module, module_name, robust_display)],
                by = "program_uid", all.x = TRUE)
  ml_write_tsv_once(meta, file.path(output_dir, "hotspot_meta_analysis.tsv"))

  membership_state <- meta[, .(
    both_axes_meta_q = all(as.character(contract$co_primary_axes) %in% axis_id) &&
      all(is.finite(fixed_q_value) & fixed_q_value < 0.05),
    all_four_cohort_axis_directions = all(direction_concordant)
  ), by = program_uid]
  membership_state[, continuum_associated :=
                     both_axes_meta_q & all_four_cohort_axis_directions]
  membership_state[, evidence_label := fifelse(
    continuum_associated, "continuum_associated", "not_continuum_associated"
  )]
  membership_state <- merge(
    registry[, .(program_uid, cell_type, module, module_name, robust_display)],
    membership_state, by = "program_uid", all.x = TRUE
  )
  ml_write_tsv_once(membership_state,
                    file.path(output_dir, "hotspot_continuum_membership.tsv"))

  requested_replicates <- contract$permutation_replicates
  if (identical(Sys.getenv("HAC_ML_TEST_MODE", unset = "0"), "1")) {
    requested_replicates <- as.integer(Sys.getenv(
      "HAC_ML_PERMUTATIONS", unset = as.character(requested_replicates)
    ))
  }
  ml_assert(requested_replicates == 10000L ||
              identical(Sys.getenv("HAC_ML_TEST_MODE", unset = "0"), "1"),
            "Production Hotspot permutations must use 10,000 replicates")
  permutation_rows <- list()
  for (axis_index in seq_along(contract$co_primary_axes)) {
    axis_name <- contract$co_primary_axes[[axis_index]]
    for (cohort_index in seq_along(evaluation)) {
      cohort <- evaluation[[cohort_index]]
      joined <- merge(
        scores[dataset == cohort],
        axes[axis_id == axis_name & dataset == cohort, .(sample_id, axis_raw)],
        by = "sample_id", all = FALSE
      )
      for (program_index in seq_along(registry_ids)) {
        program <- registry_ids[[program_index]]
        result <- hotspot_within_stage_permutation(
          joined[program_uid == program], requested_replicates,
          contract$seed + axis_index * 100000L + cohort_index * 10000L + program_index
        )
        result[, `:=`(
          program_uid = program, axis_id = axis_name, dataset = cohort,
          permutation_strata = "fibrosis_stage_by_inferred_sex",
          two_sided = TRUE
        )]
        permutation_rows[[length(permutation_rows) + 1L]] <- result
      }
    }
  }
  permutations <- rbindlist(permutation_rows, fill = TRUE)
  permutations[, permutation_q_value := ml_complete_bh(
    permutation_p_value, contract$program_family_size
  ), by = .(axis_id, dataset)]
  permutations <- merge(
    permutations,
    registry[, .(program_uid, cell_type, module, module_name, robust_display)],
    by = "program_uid", all.x = TRUE
  )
  ml_write_tsv_once(permutations,
                    file.path(output_dir, "hotspot_all117_within_stage_permutations.tsv"))

  fixed_axis <- axes[axis_id == "fixed_projection"]
  window_rows <- list()
  for (cohort in evaluation) {
    joined <- merge(
      scores[dataset == cohort],
      fixed_axis[dataset == cohort, .(sample_id, axis_percentile)],
      by = "sample_id", all = FALSE
    )
    for (window_index in seq_along(contract$window_centers)) {
      center <- contract$window_centers[[window_index]]
      lower <- center - contract$window_width / 2
      upper <- center + contract$window_width / 2
      final <- window_index == length(contract$window_centers)
      block <- joined[
        is.finite(axis_percentile) & axis_percentile >= lower &
          (axis_percentile < upper | (final & axis_percentile <= upper))
      ]
      summary <- block[, .(
        n_participants = sum(is.finite(outcome_z)),
        mean_score = if (any(is.finite(outcome_z))) mean(outcome_z, na.rm = TRUE) else NA_real_,
        se_score = if (sum(is.finite(outcome_z)) > 1L) {
          sd(outcome_z, na.rm = TRUE) / sqrt(sum(is.finite(outcome_z)))
        } else {
          NA_real_
        }
      ), by = program_uid]
      complete <- data.table(program_uid = registry_ids)
      summary <- merge(complete, summary, by = "program_uid", all.x = TRUE)
      summary[is.na(n_participants), n_participants := 0L]
      summary[, `:=`(
        dataset = cohort, axis_id = "fixed_projection", window_id = window_index,
        window_center = center, window_lower = lower, window_upper = upper,
        lower_inclusive = TRUE, upper_inclusive = final, visualization_only = TRUE
      )]
      window_rows[[length(window_rows) + 1L]] <- summary
    }
  }
  windows <- rbindlist(window_rows, fill = TRUE)
  windows <- merge(windows, registry[, .(program_uid, module_name, robust_display)],
                   by = "program_uid", all.x = TRUE)
  ml_write_tsv_once(windows,
                    file.path(output_dir, "hotspot_fixed_windows.tsv.gz"))

  focal_ids <- unname(unlist(contract$focal_programs, use.names = FALSE))
  focal_gate <- meta[program_uid %in% focal_ids, .(
    both_axes_meta_q = all(as.character(contract$co_primary_axes) %in% axis_id) &&
      all(is.finite(fixed_q_value) & fixed_q_value < 0.05),
    both_axes_direction_concordant = all(direction_concordant)
  ), by = program_uid]
  focal_gate[, figure4f_cross_sectional_gate_pass :=
               both_axes_meta_q & both_axes_direction_concordant]
  focal_gate <- merge(
    registry[, .(program_uid, module_name)][program_uid %in% focal_ids],
    focal_gate, by = "program_uid", all.x = TRUE
  )
  ml_assert(nrow(focal_gate) == length(focal_ids), "Focal Hotspot gate rows are incomplete")
  ml_write_tsv_once(focal_gate, file.path(output_dir, "hotspot_focal_figure4f_gate.tsv"))

  ml_write_session_info(file.path(output_dir, "sessionInfo.txt"))
  message("HOTSPOT_CONTINUUM_COMPLETE: ", output_dir)
}

if (!identical(Sys.getenv("HAC_ML_TEST_MODE", unset = "0"), "1")) hotspot_run()
