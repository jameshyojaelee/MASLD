suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
})

bulk_required_columns <- function(x, required, label) {
  absent <- setdiff(required, names(x))
  ml_assert(!length(absent), paste0(
    label, " lacks required columns: ", paste(absent, collapse = ", ")
  ))
  invisible(TRUE)
}

bulk_validate_design <- function(design, coefficient, label) {
  ml_assert(is.matrix(design) && nrow(design) > ncol(design),
            paste0(label, " has insufficient residual degrees of freedom"))
  ml_assert(qr(design)$rank == ncol(design), paste0(label, " is rank deficient"))
  ml_assert(coefficient %in% colnames(design),
            paste0(label, " lacks coefficient ", coefficient))
  invisible(TRUE)
}

bulk_voom_fit <- function(dge, sample_ids, design, coefficient, gene_annotation,
                          cohort, analysis, axis_id = NA_character_) {
  ml_assert(setequal(sample_ids, rownames(design)),
            paste0(analysis, " design/sample identity mismatch"))
  design <- design[sample_ids, , drop = FALSE]
  bulk_validate_design(design, coefficient, paste(cohort, analysis, axis_id))
  cohort_dge <- dge[, sample_ids]
  ml_assert(identical(colnames(cohort_dge), sample_ids),
            paste0(cohort, " DGE could not be ordered to the model design"))

  voom <- limma::voomWithQualityWeights(cohort_dge, design, plot = FALSE)
  fit <- limma::eBayes(limma::lmFit(voom, design))
  index <- match(coefficient, colnames(fit$coefficients))
  se <- fit$stdev.unscaled[, index] * sqrt(fit$s2.post)
  result <- data.table(
    gene_id_versioned = rownames(fit$coefficients),
    gene_id_base = ml_base_gene_id(rownames(fit$coefficients)),
    gene_name = gene_annotation$gene_name[
      match(ml_base_gene_id(rownames(fit$coefficients)), gene_annotation$gene_id_base)
    ],
    beta = as.numeric(fit$coefficients[, index]),
    se = as.numeric(se),
    t_value = as.numeric(fit$t[, index]),
    raw_p_value = as.numeric(fit$p.value[, index]),
    ave_expr = as.numeric(fit$Amean),
    n_participants = length(sample_ids),
    residual_df = as.numeric(fit$df.total),
    dataset = cohort,
    analysis = analysis,
    axis_id = axis_id,
    coefficient = coefficient,
    normalization = "frozen_F_five_TMM_factors_plus_cohort_specific_voomWithQualityWeights"
  )
  ml_assert(nrow(result) == nrow(dge) && !anyDuplicated(result$gene_id_base),
            paste0(cohort, " ", analysis, " did not return the complete gene family"))
  ml_assert(all(is.finite(result$beta)) && all(is.finite(result$se)) &&
              all(result$se > 0) && all(is.finite(result$raw_p_value)),
            paste0(cohort, " ", analysis, " returned invalid coefficients"))
  list(table = result, voom = voom, fit = fit)
}

bulk_nested_f_test <- function(linear_sse, nonlinear_sse, df_difference,
                               nonlinear_residual_df) {
  ml_assert(df_difference > 0L, "Nested-model numerator degrees of freedom must be positive")
  ml_assert(all(nonlinear_residual_df > 0),
            "Nested-model residual degrees of freedom must be positive")
  numerator <- pmax(0, linear_sse - nonlinear_sse) / df_difference
  denominator <- nonlinear_sse / nonlinear_residual_df
  f_value <- numerator / denominator
  p_value <- pf(
    f_value, df1 = df_difference, df2 = nonlinear_residual_df,
    lower.tail = FALSE
  )
  list(f_value = f_value, p_value = p_value)
}

bulk_spline_timing <- function(voom, linear_fit, linear_design, axis_z,
                               gene_annotation, cohort, axis_id, signature_ids,
                               family_size, grid_probabilities = seq(0, 1, 0.01)) {
  ml_assert(identical(colnames(voom$E), rownames(linear_design)),
            "Spline voom/design sample order differs")
  ml_assert(length(axis_z) == nrow(linear_design) && all(is.finite(axis_z)),
            "Spline axis is incomplete")
  axis_column <- match("axis_z", colnames(linear_design))
  ml_assert(!is.na(axis_column), "Linear design lacks axis_z")
  base_design <- linear_design[, -axis_column, drop = FALSE]
  spline_basis <- splines::ns(axis_z, df = 3L)
  colnames(spline_basis) <- paste0("axis_spline_", seq_len(ncol(spline_basis)))
  nonlinear_design <- cbind(base_design, spline_basis)
  ml_assert(qr(nonlinear_design)$rank == ncol(nonlinear_design),
            paste0(cohort, " / ", axis_id, " nonlinear design is rank deficient"))
  augmented_rank <- qr(cbind(nonlinear_design, linear_design))$rank
  ml_assert(augmented_rank == ncol(nonlinear_design),
            paste0(cohort, " / ", axis_id, " linear model is not nested in the spline model"))
  df_difference <- ncol(nonlinear_design) - ncol(linear_design)
  ml_assert(df_difference == 2L,
            "A three-df natural spline must add two degrees of freedom over linear timing")

  nonlinear_fit <- limma::lmFit(voom, nonlinear_design)
  linear_sse <- linear_fit$sigma^2 * linear_fit$df.residual
  nonlinear_sse <- nonlinear_fit$sigma^2 * nonlinear_fit$df.residual
  test <- bulk_nested_f_test(
    linear_sse, nonlinear_sse, df_difference, nonlinear_fit$df.residual
  )
  genes <- rownames(nonlinear_fit$coefficients)
  result <- data.table(
    gene_id_versioned = genes,
    gene_id_base = ml_base_gene_id(genes),
    gene_name = gene_annotation$gene_name[
      match(ml_base_gene_id(genes), gene_annotation$gene_id_base)
    ],
    dataset = cohort,
    axis_id = axis_id,
    spline_df = 3L,
    df_difference = df_difference,
    f_value = as.numeric(test$f_value),
    raw_p_value = as.numeric(test$p_value),
    linear_sse = as.numeric(linear_sse),
    nonlinear_sse = as.numeric(nonlinear_sse),
    nonlinear_residual_df = as.numeric(nonlinear_fit$df.residual),
    n_participants = ncol(voom$E)
  )
  result[, axis_component := gene_id_base %in% signature_ids]
  result[, inference_eligible := !axis_component]
  result[, p_value := fifelse(inference_eligible, raw_p_value, NA_real_)]
  result[, bh_q_value := ml_complete_bh(p_value, family_size)]
  result[, `:=`(
    bh_family_size = family_size,
    comparison = "natural_spline_df3_vs_linear",
    model = "expression ~ ns(continuum_z, df=3) + factor(fibrosis_stage) + inferred_sex",
    evidence_role = fifelse(
      axis_component,
      "axis_constituent_not_independently_testable",
      "continuum_nonlinear_timing_test"
    )
  )]

  grid_axis <- as.numeric(quantile(
    axis_z, probs = grid_probabilities, names = FALSE, type = 8
  ))
  grid_basis <- predict(spline_basis, newx = grid_axis)
  spline_columns <- match(colnames(spline_basis), colnames(nonlinear_fit$coefficients))
  ml_assert(!anyNA(spline_columns), "Spline coefficients are absent from nonlinear fit")
  predictions <- nonlinear_fit$coefficients[, spline_columns, drop = FALSE] %*%
    t(grid_basis)
  rownames(predictions) <- genes
  colnames(predictions) <- sprintf("p%03d", round(100 * grid_probabilities))
  list(table = result, predictions = predictions, grid_axis = grid_axis)
}

bulk_row_spearman <- function(x, y) {
  ml_assert(is.matrix(x) && is.matrix(y) && identical(dim(x), dim(y)),
            "Shape matrices have incompatible dimensions")
  vapply(seq_len(nrow(x)), function(index) {
    suppressWarnings(cor(x[index, ], y[index, ], method = "spearman"))
  }, numeric(1))
}

bulk_meta_table <- function(cohort_table, family_size, signature_ids = character()) {
  required <- c("gene_id_versioned", "gene_id_base", "gene_name", "axis_id",
                "dataset", "beta", "se")
  bulk_required_columns(cohort_table, required, "cohort coefficient table")
  result <- cohort_table[, {
    fixed <- ml_fixed_meta(beta, se)
    random <- ml_random_meta(beta, se)
    list(
      estimable = fixed$estimable,
      beta = fixed$beta,
      se = fixed$se,
      z_value = fixed$z,
      raw_p_value = fixed$p_value,
      ci_low = fixed$ci_low,
      ci_high = fixed$ci_high,
      random_beta = random$beta,
      random_se = random$se,
      raw_random_p_value = random$p_value,
      random_ci_low = random$ci_low,
      random_ci_high = random$ci_high,
      random_tau2 = random$tau2,
      n_cohorts = fixed$n_cohorts,
      direction_concordant = uniqueN(sign(beta[is.finite(beta) & beta != 0])) == 1L,
      cohort_directions = paste(dataset, sign(beta), sep = ":", collapse = ";")
    )
  }, by = .(gene_id_versioned, gene_id_base, gene_name, axis_id)]
  result[, axis_component := gene_id_base %in% signature_ids]
  result[, inference_eligible := !axis_component]
  result[, p_value := fifelse(inference_eligible, raw_p_value, NA_real_)]
  result[, random_p_value := fifelse(inference_eligible, raw_random_p_value, NA_real_)]
  result[, bh_q_value := ml_complete_bh(p_value, family_size), by = axis_id]
  result[, `:=`(
    bh_family_size = family_size,
    evidence_role = fifelse(
      axis_component,
      "axis_constituent_not_independently_testable",
      "continuum_association_stage_sex_adjusted"
    )
  )]
  result
}

bulk_continuum_membership <- function(meta_table, cohort_table, axes,
                                      family_size, signature_ids = character()) {
  ml_assert(setequal(unique(meta_table$axis_id), axes),
            "Continuum meta-analysis does not contain exactly the co-primary axes")
  meta_wide <- dcast(
    meta_table,
    gene_id_versioned + gene_id_base + gene_name + axis_component ~ axis_id,
    value.var = c("beta", "bh_q_value", "direction_concordant")
  )
  expected_q <- paste0("bh_q_value_", axes)
  expected_beta <- paste0("beta_", axes)
  expected_direction <- paste0("direction_concordant_", axes)
  bulk_required_columns(meta_wide, c(expected_q, expected_beta, expected_direction),
                        "wide continuum meta-analysis")

  cohort_direction <- cohort_table[, .(
    four_model_direction_concordant = .N == length(axes) * 2L &&
      uniqueN(sign(beta[is.finite(beta) & beta != 0])) == 1L,
    n_finite_cohort_axis_betas = sum(is.finite(beta)),
    cohort_axis_directions = paste(dataset, axis_id, sign(beta), sep = ":", collapse = ";")
  ), by = .(gene_id_versioned, gene_id_base)]
  result <- merge(meta_wide, cohort_direction,
                  by = c("gene_id_versioned", "gene_id_base"), all.x = TRUE, sort = FALSE)
  q_pass <- rowSums(as.data.frame(result[, ..expected_q]) < 0.05, na.rm = FALSE) == length(expected_q)
  direction_pass <- rowSums(as.data.frame(result[, ..expected_direction]), na.rm = FALSE) ==
    length(expected_direction)
  result[, continuum_associated := !axis_component & q_pass & direction_pass &
           four_model_direction_concordant]
  result[, `:=`(
    membership_class = fifelse(
      axis_component, "axis_constituent_not_independently_testable",
      fifelse(continuum_associated, "continuum_associated", "not_continuum_associated")
    ),
    bh_family_size = family_size,
    required_axes = paste(axes, collapse = ";"),
    required_evaluation_cohorts = "GSE162694;GSE213621"
  )]
  setorder(result, gene_id_versioned)
  result
}

bulk_weighted_coefficient <- function(y, design, weights, coefficient) {
  keep <- is.finite(y) & is.finite(weights) & weights > 0 &
    rowSums(!is.finite(design)) == 0L
  y <- y[keep]
  weights <- weights[keep]
  design <- design[keep, , drop = FALSE]
  if (length(y) <= ncol(design) || qr(design)$rank != ncol(design) ||
      !coefficient %in% colnames(design)) {
    return(data.table(
      estimable = FALSE, beta = NA_real_, se = NA_real_, t_value = NA_real_,
      p_value = NA_real_, residual_df = NA_integer_, n_participants = length(y)
    ))
  }
  fit <- lm.wfit(design, y, weights)
  residual_df <- length(y) - fit$rank
  sigma2 <- sum(weights * fit$residuals^2) / residual_df
  coefficient_index <- match(coefficient, colnames(design))
  unscaled <- chol2inv(qr.R(fit$qr))[coefficient_index, coefficient_index]
  se <- sqrt(sigma2 * unscaled)
  beta <- unname(fit$coefficients[[coefficient]])
  t_value <- beta / se
  data.table(
    estimable = is.finite(beta) && is.finite(se) && se > 0,
    beta = beta,
    se = se,
    t_value = t_value,
    p_value = 2 * pt(abs(t_value), residual_df, lower.tail = FALSE),
    residual_df = residual_df,
    n_participants = length(y)
  )
}

bulk_orient_pc1 <- function(expression, frozen_loading) {
  fit <- prcomp(t(expression), center = TRUE, scale. = FALSE, rank. = 1L)
  loading <- fit$rotation[, 1L]
  shared <- intersect(names(loading), names(frozen_loading))
  ml_assert(length(shared) >= 100L, "LOO signature PC1 has too few frozen loading genes")
  cosine <- sum(loading[shared] * frozen_loading[shared]) /
    sqrt(sum(loading[shared]^2) * sum(frozen_loading[shared]^2))
  ml_assert(is.finite(cosine) && abs(cosine) > sqrt(.Machine$double.eps),
            "LOO signature PC1 orientation is undefined")
  as.numeric(fit$x[, 1L]) * if (cosine < 0) -1 else 1
}

bulk_loo_axis <- function(axis_id, target_gene, expression, cohort_loadings,
                          projection_loadings) {
  genes <- rownames(expression)
  ml_assert(target_gene %in% genes, paste0("LOO target is absent: ", target_gene))
  if (axis_id == "signature_pc1") {
    loading <- cohort_loadings[axis_id == "signature_pc1"]
    frozen <- loading$frozen_loading
    names(frozen) <- loading$gene_id_base
    use <- setdiff(intersect(genes, names(frozen)), target_gene)
    ml_assert(length(use) >= 100L, "LOO signature PC1 retains too few genes")
    score <- bulk_orient_pc1(expression[use, , drop = FALSE], frozen)
  } else if (axis_id == "fixed_projection") {
    loading <- projection_loadings[used_for_fixed_projection == TRUE]
    weights <- loading$discovery_loading_oriented
    centers <- loading$discovery_center
    names(weights) <- loading$gene_id_base
    names(centers) <- loading$gene_id_base
    use <- setdiff(intersect(genes, names(weights)), target_gene)
    ml_assert(length(use) >= 100L, "LOO fixed projection retains too few genes")
    score <- as.numeric(crossprod(
      weights[use], sweep(expression[use, , drop = FALSE], 1L, centers[use], "-")
    ))
  } else {
    ml_fail("Unsupported co-primary LOO axis: ", axis_id)
  }
  names(score) <- colnames(expression)
  score
}

bulk_rank_cor <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3L) return(NA_real_)
  suppressWarnings(cor(x[keep], y[keep], method = "spearman"))
}
