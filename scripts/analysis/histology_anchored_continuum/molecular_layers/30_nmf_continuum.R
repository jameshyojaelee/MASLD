#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(splines)
})

.nmf_args <- commandArgs(trailingOnly = FALSE)
.nmf_file <- sub("^--file=", "", .nmf_args[grepl("^--file=", .nmf_args)])
.nmf_dir <- if (length(.nmf_file)) dirname(normalizePath(.nmf_file)) else getwd()
.nmf_dir_override <- Sys.getenv("HAC_ML_SCRIPT_DIR", unset = "")
if (nzchar(.nmf_dir_override)) .nmf_dir <- normalizePath(.nmf_dir_override)
source(file.path(.nmf_dir, "lib_molecular_layers.R"))

nmf_fit_continuum <- function(data, stage_column) {
  d <- copy(data)
  d <- d[
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

nmf_adjusted_shape <- function(data, stage_column) {
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

nmf_nnls_coefficients <- function(basis, expression) {
  ml_assert(requireNamespace("nnls", quietly = TRUE),
            "The existing nnls package is required for fixed-W reconstruction")
  ml_assert(nrow(basis) == nrow(expression), "Fixed-W expression row mismatch")
  answer <- vapply(seq_len(ncol(expression)), function(index) {
    nnls::nnls(basis, expression[, index])$x
  }, numeric(ncol(basis)))
  rownames(answer) <- colnames(basis)
  colnames(answer) <- colnames(expression)
  answer
}

nmf_build_fixed_basis <- function(cache_path, annotation_path, signature_symbols,
                                  evaluation_cohorts, manifest) {
  ml_assert(requireNamespace("NMF", quietly = TRUE),
            "The frozen NMF cache requires the existing NMF package")
  cache <- readRDS(cache_path)
  ml_assert(all(c("mat_nn", "nmf_results") %in% names(cache)),
            "Frozen NMF cache schema drift")
  annotation <- fread(annotation_path, select = c("gene_id", "gene_name"))
  annotation[, gene_id_base := ml_base_gene_id(gene_id)]
  annotation <- unique(annotation, by = "gene_id_base")
  cache_gene_base <- ml_base_gene_id(rownames(cache$mat_nn))
  cache_symbol <- toupper(annotation$gene_name[match(cache_gene_base, annotation$gene_id_base)])
  keep_signature_excluded <- is.na(cache_symbol) |
    !cache_symbol %in% toupper(signature_symbols)

  score_rows <- list()
  gate_rows <- list()
  for (k_value in c(4L, 6L)) {
    fit <- cache$nmf_results[[as.character(k_value)]]
    ml_assert(!is.null(fit), paste0("Frozen NMF k=", k_value, " fit is absent"))
    basis <- NMF::basis(fit)
    released_h <- NMF::coef(fit)
    ml_assert(identical(rownames(basis), rownames(cache$mat_nn)),
              paste0("Frozen NMF k=", k_value, " basis rows do not match mat_nn"))
    ml_assert(identical(colnames(released_h), colnames(cache$mat_nn)),
              paste0("Frozen NMF k=", k_value, " coefficient samples drift"))
    colnames(basis) <- paste0("P", seq_len(k_value))
    rownames(released_h) <- colnames(basis)

    full_h <- nmf_nnls_coefficients(basis, cache$mat_nn)
    excluded_basis <- basis[keep_signature_excluded, , drop = FALSE]
    excluded_expression <- cache$mat_nn[keep_signature_excluded, , drop = FALSE]
    excluded_h <- nmf_nnls_coefficients(excluded_basis, excluded_expression)

    for (program_index in seq_len(k_value)) {
      program_code <- paste0("P", program_index)
      outcome_id <- paste0("k", k_value, "_", program_code)
      retained_l1 <- sum(excluded_basis[, program_index]) / sum(basis[, program_index])
      full_rho <- suppressWarnings(cor(
        full_h[program_index, ], released_h[program_index, ], method = "spearman"
      ))
      correlation_by_cohort <- setNames(rep(NA_real_, length(evaluation_cohorts)),
                                        evaluation_cohorts)
      for (cohort in evaluation_cohorts) {
        samples <- manifest[dataset == cohort, sample_id]
        samples <- intersect(samples, colnames(released_h))
        if (length(samples) >= 3L) {
          correlation_by_cohort[[cohort]] <- suppressWarnings(cor(
            excluded_h[program_index, samples], released_h[program_index, samples],
            method = "spearman"
          ))
        }
      }
      gate_pass <- is.finite(full_rho) && full_rho >= 0.99 &&
        is.finite(retained_l1) && retained_l1 >= 0.80 &&
        all(is.finite(correlation_by_cohort) & correlation_by_cohort >= 0.90)
      gate_rows[[length(gate_rows) + 1L]] <- data.table(
        outcome_id = outcome_id, k = k_value, program_code = program_code,
        n_basis_genes = nrow(basis),
        n_signature_excluded_basis_genes = sum(!keep_signature_excluded),
        n_retained_basis_genes = nrow(excluded_basis),
        retained_l1_fraction = retained_l1,
        full_fixed_w_vs_released_h_rho = full_rho,
        evaluation_cohort = names(correlation_by_cohort),
        signature_excluded_vs_native_rho = as.numeric(correlation_by_cohort),
        reconstruction_gate_pass = gate_pass
      )
      score_rows[[length(score_rows) + 1L]] <- data.table(
        sample_id = colnames(excluded_h), outcome_id = outcome_id,
        k = k_value, program_code = program_code,
        score_raw = as.numeric(excluded_h[program_index, ]),
        score_variant = "signature_excluded_fixed_w"
      )
    }
  }
  list(scores = rbindlist(score_rows), gates = rbindlist(gate_rows))
}

nmf_run <- function() {
  contract <- ml_read_contract()
  out_root <- ml_out_root(must_exist = TRUE)
  output_dir <- file.path(out_root, "programs", "nmf")
  ml_assert(!dir.exists(output_dir), paste0("Refusing to overwrite: ", output_dir))
  ml_ensure_dir(output_dir)

  manifest <- fread(ml_resolve(contract$sample_manifest))
  ml_assert(!anyDuplicated(manifest$analysis_unit_id),
            "NMF inference requires one row per biological participant")
  axes <- ml_load_axes(contract)
  native_source <- fread(ml_resolve(contract$nmf_loadings))
  ml_assert(all(c("sample_id", "k", "program_code", "program_label",
                  "continuous_loading") %in% names(native_source)),
            "Frozen ten-axis NMF loading schema drift")
  native_source[, outcome_id := paste0("k", k, "_", program_code)]
  ml_assert(uniqueN(native_source$outcome_id) == contract$nmf_family_size,
            "Frozen NMF family is not ten axes")
  ml_assert(all(native_source[, uniqueN(sample_id), by = outcome_id]$V1 == 1104L),
            "Frozen NMF loading participant census drift")

  native <- merge(
    native_source[, .(
      sample_id, outcome_id, k, program_code, program_label,
      score_raw = continuous_loading
    )],
    manifest[, .(
      sample_id, analysis_unit_id, dataset, inferred_sex, fibrosis_stage, nas_score
    )],
    by = "sample_id", all = FALSE
  )
  native[, `:=`(
    score_variant = "native_released_loading",
    outcome_z = ml_standardize(score_raw)
  ), by = .(dataset, outcome_id)]

  fixed_state <- data.table(
    analysis = "signature_excluded_fixed_w",
    state = "blocked_source_or_dependency",
    reason = NA_character_
  )
  fixed <- NULL
  fixed_gates <- data.table()
  fixed_result <- tryCatch(
    nmf_build_fixed_basis(
      ml_resolve(contract$nmf_cache), ml_resolve(contract$gene_annotation),
      ml_signature_symbols(contract), as.character(contract$evaluation_cohorts), manifest
    ),
    error = function(e) e
  )
  if (inherits(fixed_result, "error")) {
    fixed_state[, reason := conditionMessage(fixed_result)]
  } else {
    fixed_gates <- fixed_result$gates
    fixed <- merge(
      fixed_result$scores,
      unique(native_source[, .(outcome_id, program_label)]),
      by = "outcome_id", all.x = TRUE
    )
    fixed <- merge(
      fixed,
      manifest[, .(
        sample_id, analysis_unit_id, dataset, inferred_sex, fibrosis_stage, nas_score
      )],
      by = "sample_id", all = FALSE
    )
    fixed[, outcome_z := ml_standardize(score_raw), by = .(dataset, outcome_id)]
    fixed_state[, `:=`(
      state = if (all(fixed_gates$reconstruction_gate_pass)) {
        "all_ten_axes_inference_eligible"
      } else {
        "partial_or_failed_gates_native_descriptive_only"
      },
      reason = paste0(
        uniqueN(fixed_gates[reconstruction_gate_pass == TRUE, outcome_id]),
        "/", contract$nmf_family_size, " axes pass all fixed-W gates"
      )
    )]
  }
  ml_write_tsv_once(fixed_state, file.path(output_dir, "nmf_fixed_basis_state.tsv"))
  if (nrow(fixed_gates)) {
    ml_write_tsv_once(fixed_gates,
                      file.path(output_dir, "nmf_fixed_basis_reconstruction.tsv"))
  } else {
    ml_write_tsv_once(data.table(
      outcome_id = character(), reconstruction_gate_pass = logical()
    ), file.path(output_dir, "nmf_fixed_basis_reconstruction.tsv"))
  }

  scores <- rbindlist(list(native, fixed), use.names = TRUE, fill = TRUE)
  gate_summary <- if (nrow(fixed_gates)) {
    unique(fixed_gates[, .(outcome_id, reconstruction_gate_pass)])
  } else {
    data.table(outcome_id = unique(scores$outcome_id), reconstruction_gate_pass = FALSE)
  }
  scores <- merge(scores, gate_summary, by = "outcome_id", all.x = TRUE)
  scores[, inference_eligible :=
           score_variant == "signature_excluded_fixed_w" & reconstruction_gate_pass == TRUE]
  ml_write_tsv_once(scores[, .(
    sample_id, analysis_unit_id, dataset, inferred_sex, fibrosis_stage, nas_score,
    outcome_id, k, program_code, program_label, score_variant, score_raw, outcome_z,
    reconstruction_gate_pass, inference_eligible
  )], file.path(output_dir, "nmf_participant_scores.tsv.gz"))

  evaluation <- as.character(contract$evaluation_cohorts)
  model_rows <- list()
  nas_rows <- list()
  shape_rows <- list()
  for (variant in unique(scores$score_variant)) {
    for (axis_name in as.character(contract$co_primary_axes)) {
      for (cohort in evaluation) {
        joined <- merge(
          scores[score_variant == variant & dataset == cohort],
          axes[axis_id == axis_name & dataset == cohort,
               .(sample_id, axis_raw, axis_percentile)],
          by = "sample_id", all = FALSE
        )
        for (outcome in unique(native$outcome_id)) {
          block <- joined[outcome_id == outcome]
          fit <- nmf_fit_continuum(block, "fibrosis_stage")
          fit[, `:=`(
            outcome_id = outcome, score_variant = variant, axis_id = axis_name,
            dataset = cohort, adjustment = "factor(fibrosis_stage)+inferred_sex"
          )]
          model_rows[[length(model_rows) + 1L]] <- fit
          shape <- nmf_adjusted_shape(block, "fibrosis_stage")
          shape[, `:=`(
            outcome_id = outcome, score_variant = variant, axis_id = axis_name,
            dataset = cohort, adjustment = "factor(fibrosis_stage)+inferred_sex"
          )]
          shape_rows[[length(shape_rows) + 1L]] <- shape
          if (cohort == "GSE162694") {
            nas <- nmf_fit_continuum(block, "nas_score")
            nas[, `:=`(
              outcome_id = outcome, score_variant = variant, axis_id = axis_name,
              dataset = cohort, adjustment = "factor(nas_score)+inferred_sex"
            )]
            nas_rows[[length(nas_rows) + 1L]] <- nas
          }
        }
      }
    }
  }
  models <- rbindlist(model_rows, fill = TRUE)
  models[, `:=`(
    q_value = ml_complete_bh(p_value, contract$nmf_family_size),
    nonlinear_q_value = ml_complete_bh(
      nonlinear_p_value, contract$nmf_family_size
    )
  ), by = .(score_variant, axis_id, dataset)]
  models <- merge(models, unique(native[, .(outcome_id, k, program_code, program_label)]),
                  by = "outcome_id", all.x = TRUE)
  models <- merge(models, gate_summary, by = "outcome_id", all.x = TRUE)
  models[, inference_eligible :=
           score_variant == "signature_excluded_fixed_w" & reconstruction_gate_pass == TRUE]
  ml_write_tsv_once(models, file.path(output_dir, "nmf_fibrosis_models.tsv"))

  nas_models <- rbindlist(nas_rows, fill = TRUE)
  nas_models[, q_value := ml_complete_bh(p_value, contract$nmf_family_size),
             by = .(score_variant, axis_id)]
  nas_models <- merge(
    nas_models, unique(native[, .(outcome_id, k, program_code, program_label)]),
    by = "outcome_id", all.x = TRUE
  )
  ml_write_tsv_once(nas_models, file.path(output_dir, "nmf_nas_sensitivity.tsv"))

  shapes <- rbindlist(shape_rows, fill = TRUE)
  ml_write_tsv_once(shapes, file.path(output_dir, "nmf_adjusted_spline_shapes.tsv.gz"))
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
  }, by = .(outcome_id, score_variant, axis_id)]
  nonlinear_gate <- models[, .(
    nonlinear_bh_both_cohorts = .N == length(evaluation) &&
      all(is.finite(nonlinear_q_value) & nonlinear_q_value < 0.05)
  ), by = .(outcome_id, score_variant, axis_id)]
  shape_concordance <- merge(shape_concordance, nonlinear_gate,
                             by = c("outcome_id", "score_variant", "axis_id"))
  shape_concordance[, timing_gate_pass := nonlinear_bh_both_cohorts &
                      is.finite(shape_spearman_rho) & shape_spearman_rho >= 0.70]
  ml_write_tsv_once(shape_concordance,
                    file.path(output_dir, "nmf_nonlinear_timing_gates.tsv"))

  meta_rows <- list()
  for (variant in unique(models$score_variant)) {
    for (axis_name in as.character(contract$co_primary_axes)) {
      for (outcome in unique(models$outcome_id)) {
        block <- models[
          score_variant == variant & axis_id == axis_name & outcome_id == outcome &
            dataset %in% evaluation & estimable == TRUE
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
          outcome_id = outcome, score_variant = variant, axis_id = axis_name,
          direction_concordant = nrow(block) == length(evaluation) &&
            length(unique(sign(block$beta))) == 1L,
          cohort_betas = paste(block$dataset, signif(block$beta, 6),
                               sep = "=", collapse = ";")
        )]
        meta_rows[[length(meta_rows) + 1L]] <- row
      }
    }
  }
  meta <- rbindlist(meta_rows, fill = TRUE)
  meta[, fixed_q_value := ml_complete_bh(fixed_p_value, contract$nmf_family_size),
       by = .(score_variant, axis_id)]
  meta <- merge(meta, unique(native[, .(outcome_id, k, program_code, program_label)]),
                by = "outcome_id", all.x = TRUE)
  meta <- merge(meta, gate_summary, by = "outcome_id", all.x = TRUE)
  meta[, inference_eligible :=
         score_variant == "signature_excluded_fixed_w" & reconstruction_gate_pass == TRUE]
  ml_write_tsv_once(meta, file.path(output_dir, "nmf_meta_analysis.tsv"))

  membership <- meta[score_variant == "signature_excluded_fixed_w", .(
    both_axes_meta_q = all(as.character(contract$co_primary_axes) %in% axis_id) &&
      all(is.finite(fixed_q_value) & fixed_q_value < 0.05),
    all_four_cohort_axis_directions = all(direction_concordant),
    fixed_w_gate_pass = all(reconstruction_gate_pass)
  ), by = outcome_id]
  membership <- merge(
    data.table(outcome_id = unique(native$outcome_id)), membership,
    by = "outcome_id", all.x = TRUE
  )
  membership[is.na(both_axes_meta_q), `:=`(
    both_axes_meta_q = FALSE,
    all_four_cohort_axis_directions = FALSE,
    fixed_w_gate_pass = FALSE
  )]
  membership[, continuum_associated :=
               both_axes_meta_q & all_four_cohort_axis_directions & fixed_w_gate_pass]
  membership[, evidence_label := fifelse(
    continuum_associated, "continuum_associated", "same_substrate_alignment_descriptive"
  )]
  ml_write_tsv_once(membership,
                    file.path(output_dir, "nmf_continuum_membership.tsv"))

  fixed_axis <- axes[axis_id == "fixed_projection"]
  window_rows <- list()
  for (variant in unique(scores$score_variant)) {
    for (cohort in evaluation) {
      joined <- merge(
        scores[score_variant == variant & dataset == cohort],
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
          mean_score = mean(outcome_z, na.rm = TRUE),
          se_score = sd(outcome_z, na.rm = TRUE) / sqrt(sum(is.finite(outcome_z)))
        ), by = outcome_id]
        summary[, `:=`(
          score_variant = variant, dataset = cohort, axis_id = "fixed_projection",
          window_id = window_index, window_center = center,
          window_lower = lower, window_upper = upper,
          lower_inclusive = TRUE, upper_inclusive = final,
          visualization_only = TRUE
        )]
        window_rows[[length(window_rows) + 1L]] <- summary
      }
    }
  }
  windows <- rbindlist(window_rows, fill = TRUE)
  windows <- merge(windows,
                   unique(native[, .(outcome_id, k, program_code, program_label)]),
                   by = "outcome_id", all.x = TRUE)
  ml_write_tsv_once(windows, file.path(output_dir, "nmf_fixed_windows.tsv.gz"))

  ml_write_session_info(file.path(output_dir, "sessionInfo.txt"))
  message("NMF_CONTINUUM_COMPLETE: ", output_dir)
}

if (!identical(Sys.getenv("HAC_ML_TEST_MODE", unset = "0"), "1")) nmf_run()
