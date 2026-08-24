#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(Matrix)
  library(splines)
})

ml_pathway_script_dir <- Sys.getenv("HAC_ML_SCRIPT_DIR", unset = "")
if (!nzchar(ml_pathway_script_dir)) {
  ml_pathway_args <- commandArgs(trailingOnly = FALSE)
  ml_pathway_file <- sub("^--file=", "", ml_pathway_args[grepl("^--file=", ml_pathway_args)])
  stopifnot(length(ml_pathway_file) == 1L)
  ml_pathway_script_dir <- dirname(normalizePath(ml_pathway_file))
}
source(file.path(ml_pathway_script_dir, "lib_molecular_layers.R"))

pathway_pick_input <- function(env_name, candidates, description) {
  override <- Sys.getenv(env_name, unset = "")
  if (nzchar(override)) return(ml_resolve(override))
  existing <- candidates[file.exists(candidates)]
  ml_assert(
    length(existing) == 1L,
    paste0(
      description, " requires exactly one upstream input; found ", length(existing),
      ". Set ", env_name, ". Candidates: ", paste(candidates, collapse = ";")
    )
  )
  normalizePath(existing, mustWork = TRUE)
}

pathway_read_gmt <- function(path, collection, expected_size) {
  fields <- strsplit(readLines(path, warn = FALSE), "\t", fixed = TRUE)
  ml_assert(length(fields) == expected_size,
            paste0(collection, " GMT family-size drift: ", length(fields)))
  ids <- vapply(fields, `[[`, character(1), 1L)
  ml_assert(!anyDuplicated(ids), paste0(collection, " GMT has duplicated set IDs"))
  members <- rbindlist(lapply(seq_along(fields), function(index) {
    genes <- unique(toupper(trimws(fields[[index]][-(1:2)])))
    genes <- genes[!is.na(genes) & nzchar(genes)]
    data.table(set_id = ids[[index]], gene_symbol = genes)
  }))
  ml_assert(uniqueN(members$set_id) == expected_size,
            paste0(collection, " contains an empty gene set"))
  list(ids = ids, membership = members)
}

pathway_symbol_expression <- function(dge_path, manifest_path, annotation_path) {
  dge <- readRDS(dge_path)
  meta <- fread(manifest_path, na.strings = c("", "NA"))
  required_meta <- c("sample_id", "analysis_unit_id", "dataset", "inferred_sex",
                     "fibrosis_stage", "nas_score")
  ml_assert(all(required_meta %in% names(meta)), "Sample manifest schema drift")
  ml_assert(!anyDuplicated(meta$sample_id), "Sample manifest contains duplicate samples")
  ml_assert(!anyDuplicated(meta$analysis_unit_id),
            "Cells or technical rows cannot be donor replicates")
  counts <- if (inherits(dge, "DGEList")) dge$counts else dge$counts
  ml_assert(!is.null(counts), "DGE object has no count matrix")
  ml_assert(setequal(colnames(counts), meta$sample_id), "DGE/manifest sample sets differ")
  log_cpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
  log_cpm <- log_cpm[, meta$sample_id, drop = FALSE]

  annotation <- fread(annotation_path, select = c("gene_id", "gene_name"))
  annotation[, gene_id_base := ml_base_gene_id(gene_id)]
  annotation[, gene_symbol := toupper(trimws(gene_name))]
  annotation <- unique(annotation[
    !is.na(gene_symbol) & nzchar(gene_symbol), .(gene_id_base, gene_symbol)
  ])
  ml_assert(annotation[, uniqueN(gene_symbol), by = gene_id_base][V1 > 1L, .N] == 0L,
            "One Ensembl base ID maps to multiple symbols")
  map <- annotation[match(ml_base_gene_id(rownames(log_cpm)), gene_id_base)]
  keep <- !is.na(map$gene_symbol) & nzchar(map$gene_symbol)
  expression <- rowsum(log_cpm[keep, , drop = FALSE], map$gene_symbol[keep],
                       reorder = FALSE, na.rm = TRUE)
  divisor <- as.numeric(table(map$gene_symbol[keep])[rownames(expression)])
  expression <- expression / divisor
  storage.mode(expression) <- "double"
  list(expression = expression, metadata = meta)
}

pathway_gene_annotation <- function(gene_ids, annotation_path) {
  annotation <- fread(annotation_path, select = c("gene_id", "gene_name"))
  annotation[, `:=`(
    gene_id_base = ml_base_gene_id(gene_id),
    gene_symbol = toupper(trimws(gene_name))
  )]
  annotation <- unique(annotation[, .(gene_id_base, gene_symbol)])
  ml_assert(annotation[, uniqueN(gene_symbol), by = gene_id_base][V1 > 1L, .N] == 0L,
            "One Ensembl base ID maps to multiple symbols")
  result <- annotation[match(ml_base_gene_id(gene_ids), gene_id_base)]
  data.table(
    gene_id = as.character(gene_ids), gene_id_base = ml_base_gene_id(gene_ids),
    gene_symbol = result$gene_symbol
  )
}

pathway_voom_statistics <- function(dge, samples, design, coefficient, annotation,
                                    dataset, axis_id, model_formula) {
  ml_assert(length(samples) == nrow(design), "Voom sample/design dimensions differ")
  subset <- dge[, samples, keep.lib.sizes = FALSE]
  subset <- edgeR::calcNormFactors(subset)
  voom <- limma::voomWithQualityWeights(subset, design = design, plot = FALSE)
  fit <- limma::eBayes(limma::lmFit(voom, design), robust = TRUE)
  coefficient_index <- match(coefficient, colnames(design))
  ml_assert(!is.na(coefficient_index), paste0("Missing model coefficient: ", coefficient))
  estimate <- fit$coefficients[, coefficient_index]
  moderated_se <- fit$stdev.unscaled[, coefficient_index] * sqrt(fit$s2.post)
  result <- copy(annotation)
  result[, `:=`(
    dataset = dataset, axis_id = axis_id,
    beta = as.numeric(estimate), se = as.numeric(moderated_se),
    moderated_t = as.numeric(fit$t[, coefficient_index]),
    p_value = as.numeric(fit$p.value[, coefficient_index]),
    average_expression = as.numeric(fit$Amean),
    n_participants = length(samples), model_formula = model_formula,
    inference_method = "limma_voomWithQualityWeights_eBayes_robust"
  )]
  result
}

pathway_generate_rank_statistics <- function(dge_path, manifest_path, annotation_path,
                                             axes, contract) {
  dge <- readRDS(dge_path)
  metadata <- fread(manifest_path, na.strings = c("", "NA"))
  ml_assert(nrow(dge$counts) == contract$gene_family_size,
            "Expression universe is not the frozen 23,370-gene family")
  ml_assert(setequal(colnames(dge$counts), metadata$sample_id),
            "DGE/manifest sample sets differ during rank-statistic generation")
  annotation <- pathway_gene_annotation(rownames(dge$counts), annotation_path)
  continuum_rows <- list()
  for (cohort in contract$evaluation_cohorts) {
    for (axis_name in contract$co_primary_axes) {
      d <- merge(
        metadata[dataset == cohort],
        axes[dataset == cohort & axis_id == axis_name],
        by = c("sample_id", "dataset"), all = FALSE, sort = FALSE
      )
      d <- d[is.finite(axis_raw) & !is.na(fibrosis_stage) & !is.na(inferred_sex)]
      d[, `:=`(
        axis_z = ml_standardize(axis_raw),
        stage_factor = factor(fibrosis_stage),
        sex_factor = factor(inferred_sex)
      )]
      design <- model.matrix(~ axis_z + stage_factor + sex_factor, data = d)
      ml_assert(qr(design)$rank == ncol(design),
                paste0("Rank-deficient continuum gene model: ", cohort, "/", axis_name))
      continuum_rows[[length(continuum_rows) + 1L]] <- pathway_voom_statistics(
        dge, d$sample_id, design, "axis_z", annotation, cohort, axis_name,
        "expression ~ continuum_z + factor(fibrosis_stage) + inferred_sex"
      )
    }
  }
  continuum <- rbindlist(continuum_rows, use.names = TRUE, fill = TRUE)
  continuum[, q_value := ml_complete_bh(p_value, contract$gene_family_size),
            by = .(dataset, axis_id)]
  continuum[, bh_family_size := as.integer(contract$gene_family_size)]

  disease <- metadata[
    dataset %in% contract$source_overlap_cohorts & !is.na(group_binary) &
      !is.na(inferred_sex)
  ]
  disease[, `:=`(
    disease_factor = relevel(factor(group_binary), ref = "Control"),
    dataset_factor = factor(dataset),
    sex_factor = factor(inferred_sex)
  )]
  design <- model.matrix(~ disease_factor + dataset_factor + sex_factor, data = disease)
  ml_assert(qr(design)$rank == ncol(design),
            "Rank-deficient source-overlap disease-control gene model")
  disease_stats <- pathway_voom_statistics(
    dge, disease$sample_id, design, "disease_factorDisease", annotation,
    "GSE126848+GSE130970+GSE135251", "source_overlap_disease_vs_control",
    "expression ~ disease_control + cohort + inferred_sex"
  )
  disease_stats[, q_value := ml_complete_bh(p_value, contract$gene_family_size)]
  disease_stats[, bh_family_size := as.integer(contract$gene_family_size)]
  list(continuum = continuum, disease = disease_stats)
}

pathway_zscore_rows <- function(matrix) {
  center <- rowMeans(matrix, na.rm = TRUE)
  spread <- apply(matrix, 1L, sd, na.rm = TRUE)
  answer <- sweep(matrix, 1L, center, "-")
  variable <- is.finite(spread) & spread > 0
  answer[variable, ] <- sweep(answer[variable, , drop = FALSE], 1L,
                             spread[variable], "/")
  answer[!variable, ] <- NA_real_
  answer
}

pathway_score_collection <- function(expression, metadata, set_ids, membership,
                                     signature_symbols, minimum_genes,
                                     minimum_fraction) {
  membership <- unique(membership[, .(
    set_id = as.character(set_id), gene_symbol = toupper(gene_symbol)
  )])
  original <- membership[, .(n_original_genes = uniqueN(gene_symbol)), by = set_id]
  excluded <- membership[gene_symbol %in% signature_symbols, .(
    n_signature_genes = uniqueN(gene_symbol),
    signature_genes = paste(sort(unique(gene_symbol)), collapse = ";")
  ), by = set_id]
  retained <- membership[!gene_symbol %in% signature_symbols]
  retained_count <- retained[, .(n_signature_retained = uniqueN(gene_symbol)), by = set_id]

  coverage_rows <- list()
  score_blocks <- list()
  for (cohort in unique(metadata$dataset)) {
    sample_ids <- metadata[dataset == cohort, sample_id]
    z <- pathway_zscore_rows(expression[, sample_ids, drop = FALSE])
    observed_genes <- rownames(z)[rowSums(is.finite(z)) == ncol(z)]
    observed <- retained[gene_symbol %in% observed_genes]
    observed_count <- observed[, .(n_observed_genes = uniqueN(gene_symbol)), by = set_id]
    coverage <- Reduce(function(x, y) merge(x, y, by = "set_id", all.x = TRUE), list(
      data.table(set_id = set_ids), original, retained_count, observed_count, excluded
    ))
    for (column in c("n_signature_retained", "n_observed_genes", "n_signature_genes")) {
      set(coverage, which(is.na(coverage[[column]])), column, 0L)
    }
    coverage[is.na(signature_genes), signature_genes := ""]
    coverage[, `:=`(
      dataset = cohort,
      signature_retained_fraction = n_signature_retained / n_original_genes,
      observed_retained_fraction = n_observed_genes / n_original_genes
    )]
    coverage[, testable := n_observed_genes >= minimum_genes &
               observed_retained_fraction >= minimum_fraction]
    coverage[, testability_reason := fifelse(
      testable, "testable",
      fifelse(n_observed_genes < minimum_genes, "fewer_than_minimum_observed_genes",
              "less_than_80_percent_original_membership_observed_after_signature_exclusion")
    )]
    coverage_rows[[cohort]] <- coverage

    testable_ids <- coverage[testable == TRUE, set_id]
    if (!length(testable_ids)) next
    scored_membership <- observed[set_id %in% testable_ids]
    gene_levels <- observed_genes
    set_index <- match(scored_membership$set_id, testable_ids)
    gene_index <- match(scored_membership$gene_symbol, gene_levels)
    set_sizes <- tabulate(set_index, nbins = length(testable_ids))
    weights <- Matrix::sparseMatrix(
      i = set_index, j = gene_index, x = 1 / set_sizes[set_index],
      dims = c(length(testable_ids), length(gene_levels)),
      dimnames = list(testable_ids, gene_levels)
    )
    scores <- as.matrix(weights %*% z[gene_levels, , drop = FALSE])
    scores <- t(apply(scores, 1L, ml_standardize))
    if (length(testable_ids) == 1L) {
      scores <- matrix(scores, nrow = 1L,
                       dimnames = list(testable_ids, sample_ids))
    } else {
      rownames(scores) <- testable_ids
      colnames(scores) <- sample_ids
    }
    score_blocks[[cohort]] <- scores
  }
  list(coverage = rbindlist(coverage_rows, use.names = TRUE, fill = TRUE),
       scores = score_blocks)
}

pathway_fit_matrix <- function(outcomes, axis_raw, stage, sex) {
  keep <- is.finite(axis_raw) & !is.na(stage) & !is.na(sex)
  outcomes <- outcomes[, keep, drop = FALSE]
  axis_z <- ml_standardize(axis_raw[keep])
  stage <- factor(stage[keep])
  sex <- factor(sex[keep])
  n <- length(axis_z)
  if (n < 20L || nlevels(stage) < 2L || nlevels(sex) < 2L) return(NULL)
  outcome_z <- t(apply(outcomes, 1L, ml_standardize))
  complete_rows <- rowSums(is.finite(outcome_z)) == n
  base <- model.matrix(~ stage + sex)
  linear <- cbind(base, axis_z = axis_z)
  nonlinear <- cbind(base, splines::ns(axis_z, df = 3))
  ml_assert(qr(linear)$rank == ncol(linear), "Linear pathway design is rank deficient")
  ml_assert(qr(nonlinear)$rank == ncol(nonlinear),
            "Spline pathway design is rank deficient")
  fit_linear <- lm.fit(linear, t(outcome_z[complete_rows, , drop = FALSE]))
  fit_spline <- lm.fit(nonlinear, t(outcome_z[complete_rows, , drop = FALSE]))
  coefficient <- fit_linear$coefficients[ncol(linear), ]
  residual_df <- n - ncol(linear)
  sigma2 <- colSums(fit_linear$residuals^2) / residual_df
  covariance_multiplier <- solve(crossprod(linear))[ncol(linear), ncol(linear)]
  se <- sqrt(sigma2 * covariance_multiplier)
  t_value <- coefficient / se
  linear_p <- 2 * pt(abs(t_value), df = residual_df, lower.tail = FALSE)
  spline_df <- n - ncol(nonlinear)
  linear_sse <- colSums(fit_linear$residuals^2)
  spline_sse <- colSums(fit_spline$residuals^2)
  df_difference <- ncol(nonlinear) - ncol(linear)
  f_value <- pmax(0, (linear_sse - spline_sse) / df_difference) /
    (spline_sse / spline_df)
  nonlinear_p <- pf(f_value, df_difference, spline_df, lower.tail = FALSE)
  data.table(
    set_id = rownames(outcome_z)[complete_rows], estimable = TRUE, n = n,
    beta = coefficient, se = se, t_value = t_value, p_value = linear_p,
    ci_low = coefficient - qt(0.975, residual_df) * se,
    ci_high = coefficient + qt(0.975, residual_df) * se,
    nonlinear_f = f_value, nonlinear_p_value = nonlinear_p,
    spline_df = 3L
  )
}

pathway_fit_all <- function(scores, coverage, metadata, axes, set_ids,
                            evaluation_cohorts, co_primary_axes, stage_column,
                            family_size) {
  rows <- list()
  for (cohort in evaluation_cohorts) {
    score_matrix <- scores[[cohort]]
    for (axis_name in co_primary_axes) {
      base <- data.table(set_id = set_ids, dataset = cohort, axis_id = axis_name)
      axis <- axes[dataset == cohort & get("axis_id") == axis_name]
      meta <- merge(
        metadata[dataset == cohort],
        axis[, .(sample_id, axis_raw, axis_percentile)],
        by = "sample_id", all = FALSE, sort = FALSE
      )
      fit <- NULL
      if (!is.null(score_matrix)) {
        common <- intersect(meta$sample_id, colnames(score_matrix))
        meta <- meta[match(common, sample_id)]
        fit <- pathway_fit_matrix(
          score_matrix[, common, drop = FALSE], meta$axis_raw,
          meta[[stage_column]], meta$inferred_sex
        )
      }
      if (is.null(fit)) {
        fit <- data.table(
          set_id = character(), estimable = logical(), n = integer(), beta = numeric(),
          se = numeric(), t_value = numeric(), p_value = numeric(), ci_low = numeric(),
          ci_high = numeric(), nonlinear_f = numeric(), nonlinear_p_value = numeric(),
          spline_df = integer()
        )
      }
      result <- merge(base, fit, by = "set_id", all.x = TRUE, sort = FALSE)
      result <- merge(
        result,
        coverage[dataset == cohort, .(set_id, testable, testability_reason)],
        by = "set_id", all.x = TRUE, sort = FALSE
      )
      result[is.na(estimable), estimable := FALSE]
      result[, `:=`(
        stage_adjustment = paste0("factor(", stage_column, ")"),
        q_value = ml_complete_bh(p_value, family_size),
        nonlinear_q_value = ml_complete_bh(nonlinear_p_value, family_size),
        bh_family_size = family_size
      )]
      rows[[length(rows) + 1L]] <- result
    }
  }
  rbindlist(rows, use.names = TRUE, fill = TRUE)
}

pathway_meta_models <- function(models, set_ids, axes, family_size) {
  result <- rbindlist(lapply(axes, function(axis) {
    rbindlist(lapply(set_ids, function(set_name) {
      d <- models[models$axis_id == axis & models$set_id == set_name]
      fixed <- ml_fixed_meta(d$beta, d$se)
      random <- ml_random_meta(d$beta, d$se)
      data.table(
        set_id = set_name, axis_id = axis,
        fixed_beta = fixed$beta, fixed_se = fixed$se,
        fixed_p_value = fixed$p_value, fixed_ci_low = fixed$ci_low,
        fixed_ci_high = fixed$ci_high,
        random_beta = random$beta, random_se = random$se,
        random_p_value = random$p_value, random_tau2 = random$tau2,
        direction_concordant = nrow(d[estimable == TRUE]) == 2L &&
          (all(d[estimable == TRUE]$beta > 0) || all(d[estimable == TRUE]$beta < 0)),
        n_cohorts = fixed$n_cohorts
      )
    }))
  }))
  result[, fixed_q_value := ml_complete_bh(fixed_p_value, family_size), by = axis_id]
  result[, random_q_value := ml_complete_bh(random_p_value, family_size), by = axis_id]
  result[, bh_family_size := family_size]
  result
}

pathway_spline_shapes <- function(models, scores, metadata, axes, set_ids,
                                  cohorts, axis_names, grid = seq(0, 1, length.out = 101L)) {
  summary_rows <- list()
  prediction_rows <- list()
  for (axis_name in axis_names) {
    for (set_name in set_ids) {
      model_rows <- models[models$axis_id == axis_name & models$set_id == set_name]
      eligible <- nrow(model_rows[estimable == TRUE]) == length(cohorts) &&
        all(is.finite(model_rows[estimable == TRUE]$nonlinear_q_value)) &&
        all(model_rows[estimable == TRUE]$nonlinear_q_value < 0.05)
      if (!eligible) {
        summary_rows[[length(summary_rows) + 1L]] <- data.table(
          set_id = set_name, axis_id = axis_name, shape_tested = FALSE,
          shape_spearman = NA_real_, shape_concordant = FALSE,
          reason = "nonlinear_complete_family_gate_not_met",
          n_grid_points = length(grid)
        )
        next
      }
      cohort_predictions <- list()
      complete <- TRUE
      for (cohort in cohorts) {
        matrix <- scores[[cohort]]
        if (is.null(matrix) || !set_name %in% rownames(matrix)) {
          complete <- FALSE
          next
        }
        axis <- axes[dataset == cohort & axis_id == axis_name]
        d <- merge(metadata[dataset == cohort], axis,
                   by = c("sample_id", "dataset"), all = FALSE, sort = FALSE)
        d <- d[sample_id %in% colnames(matrix) & is.finite(axis_raw) &
                 !is.na(fibrosis_stage) & !is.na(inferred_sex)]
        d[, `:=`(
          outcome_z = as.numeric(matrix[set_name, sample_id]),
          axis_z = ml_standardize(axis_raw),
          stage_factor = factor(fibrosis_stage),
          sex_factor = factor(inferred_sex)
        )]
        fit <- tryCatch(
          lm(outcome_z ~ splines::ns(axis_z, df = 3) + stage_factor + sex_factor,
             data = d),
          error = function(e) NULL
        )
        if (is.null(fit)) {
          complete <- FALSE
          next
        }
        axis_grid <- as.numeric(quantile(d$axis_z, probs = grid, names = FALSE, type = 8))
        stage_reference <- names(sort(table(d$stage_factor), decreasing = TRUE))[[1L]]
        sex_reference <- names(sort(table(d$sex_factor), decreasing = TRUE))[[1L]]
        newdata <- data.frame(
          axis_z = axis_grid,
          stage_factor = factor(stage_reference, levels = levels(d$stage_factor)),
          sex_factor = factor(sex_reference, levels = levels(d$sex_factor))
        )
        prediction <- as.numeric(predict(fit, newdata = newdata))
        cohort_predictions[[cohort]] <- prediction
        prediction_rows[[length(prediction_rows) + 1L]] <- data.table(
          set_id = set_name, axis_id = axis_name, dataset = cohort,
          percentile = grid, predicted_score = prediction,
          reference_fibrosis_stage = stage_reference,
          reference_sex = sex_reference
        )
      }
      rho <- if (complete && length(cohort_predictions) == length(cohorts)) {
        suppressWarnings(cor(cohort_predictions[[cohorts[[1L]]]],
                             cohort_predictions[[cohorts[[2L]]]], method = "spearman"))
      } else {
        NA_real_
      }
      summary_rows[[length(summary_rows) + 1L]] <- data.table(
        set_id = set_name, axis_id = axis_name,
        shape_tested = is.finite(rho), shape_spearman = rho,
        shape_concordant = is.finite(rho) && rho >= 0.70,
        reason = if (is.finite(rho)) "tested" else "spline_prediction_failure",
        n_grid_points = length(grid)
      )
    }
  }
  predictions <- rbindlist(prediction_rows, use.names = TRUE, fill = TRUE)
  if (!nrow(predictions)) {
    predictions <- data.table(
      set_id = character(), axis_id = character(), dataset = character(),
      percentile = numeric(), predicted_score = numeric(),
      reference_fibrosis_stage = character(), reference_sex = character()
    )
  }
  list(summary = rbindlist(summary_rows, use.names = TRUE, fill = TRUE),
       predictions = predictions)
}

pathway_normalize_gene_stats <- function(path, co_primary_axes,
                                         evaluation_cohorts = NULL) {
  x <- fread(path)
  symbol_column <- intersect(c("gene_symbol", "symbol", "gene_name"), names(x))[1L]
  axis_column <- intersect(c("axis_id", "continuum_axis", "score_id"), names(x))[1L]
  ml_assert(!is.na(symbol_column) && !is.na(axis_column),
            "Continuum gene-statistic file needs gene-symbol and axis columns")
  setnames(x, c(symbol_column, axis_column), c("gene_symbol", "axis_id"), skip_absent = TRUE)
  if (!is.null(evaluation_cohorts)) {
    dataset_column <- intersect(c("dataset", "cohort"), names(x))[1L]
    ml_assert(!is.na(dataset_column), "Cohort gene statistics need a dataset column")
    if (dataset_column != "dataset") setnames(x, dataset_column, "dataset")
    x <- x[dataset %in% evaluation_cohorts]
  }
  x <- x[axis_id %in% co_primary_axes]
  statistic_column <- intersect(
    c("moderated_t", "t_value", "t", "z_value", "z", "statistic"), names(x)
  )[1L]
  if (!is.na(statistic_column)) {
    x[, rank_statistic := as.numeric(get(statistic_column))]
  } else {
    ml_assert(all(c("beta", "se") %in% names(x)),
              "Gene statistics require a moderated statistic or beta/se")
    x[, rank_statistic := beta / se]
  }
  x[, gene_symbol := toupper(trimws(gene_symbol))]
  x[is.finite(rank_statistic) & !is.na(gene_symbol) & nzchar(gene_symbol)]
}

pathway_meta_ranks <- function(cohort_statistics) {
  ml_assert(all(c("dataset", "beta", "se") %in% names(cohort_statistics)),
            "Cohort gene statistics need beta and se for ranked meta-analysis")
  cohort_statistics[, {
    meta <- ml_fixed_meta(beta, se)
    .(rank_statistic = meta$z, meta_beta = meta$beta, meta_se = meta$se,
      n_cohorts = meta$n_cohorts)
  }, by = .(axis_id, gene_symbol)]
}

pathway_ranked_enrichment <- function(set_ids, membership, ranks,
                                      signature_symbols, family_size, seed) {
  membership <- unique(membership[, .(
    set_id = as.character(set_id), gene_symbol = toupper(gene_symbol)
  )])
  original <- membership[, .(n_original_genes = uniqueN(gene_symbol)), by = set_id]
  retained_membership <- membership[!gene_symbol %in% signature_symbols]
  retained <- retained_membership[, .(
    n_signature_retained = uniqueN(gene_symbol)
  ), by = set_id]
  available_symbols <- unique(ranks[is.finite(rank_statistic), gene_symbol])
  available <- retained_membership[gene_symbol %in% available_symbols, .(
    n_ranked_genes = uniqueN(gene_symbol)
  ), by = set_id]
  base <- Reduce(function(x, y) merge(x, y, by = "set_id", all.x = TRUE), list(
    data.table(set_id = set_ids), original, retained, available
  ))
  for (column in c("n_signature_retained", "n_ranked_genes")) {
    set(base, which(is.na(base[[column]])), column, 0L)
  }
  base[, `:=`(
    retained_fraction = n_signature_retained / n_original_genes,
    ranked_fraction = n_ranked_genes / n_original_genes
  )]
  base[, testable := n_ranked_genes >= 10L & ranked_fraction >= 0.8]
  if (!requireNamespace("fgsea", quietly = TRUE)) {
    base[, `:=`(
      pathway_size = NA_integer_, enrichment_score = NA_real_, nes = NA_real_,
      p_value = NA_real_, q_value = NA_real_, leading_edge = "",
      state = "not_run_dependency_unavailable", bh_family_size = family_size
    )]
    return(base)
  }
  pathways <- split(
    retained_membership[set_id %in% base[testable == TRUE, set_id]]$gene_symbol,
    retained_membership[set_id %in% base[testable == TRUE, set_id]]$set_id
  )
  if (!length(pathways)) {
    base[, `:=`(
      pathway_size = NA_integer_, enrichment_score = NA_real_, nes = NA_real_,
      p_value = NA_real_, q_value = NA_real_, leading_edge = "",
      state = "untestable_lt10_or_lt80_percent_original_membership_ranked",
      bh_family_size = family_size
    )]
    return(base)
  }
  ranks <- ranks[order(-abs(rank_statistic))][!duplicated(gene_symbol)]
  statistic <- setNames(ranks$rank_statistic, ranks$gene_symbol)
  statistic <- sort(statistic[is.finite(statistic)], decreasing = TRUE)
  set.seed(seed)
  result <- as.data.table(fgsea::fgseaMultilevel(
    pathways = pathways, stats = statistic, minSize = 10L,
    maxSize = Inf, eps = 0
  ))
  if (nrow(result)) {
    result[, leadingEdge := vapply(leadingEdge, paste, collapse = ";", character(1))]
    setnames(result,
             c("pathway", "size", "ES", "NES", "pval", "padj", "leadingEdge"),
             c("set_id", "pathway_size", "enrichment_score", "nes", "p_value",
               "fgsea_padj", "leading_edge"))
  }
  result <- merge(base, result, by = "set_id", all.x = TRUE, sort = FALSE)
  result[, `:=`(
    q_value = ml_complete_bh(p_value, family_size),
    state = fifelse(
      is.finite(p_value), "tested",
      fifelse(!testable,
              "untestable_lt10_or_lt80_percent_original_membership_ranked",
              "untestable_ranked_enrichment")
    ),
    bh_family_size = family_size
  )]
  result
}

pathway_window_summary <- function(scores, metadata, axes, display_sets,
                                   cohorts, centers, width) {
  rows <- list()
  for (cohort in cohorts) {
    matrix <- scores[[cohort]]
    if (is.null(matrix)) next
    chosen <- intersect(display_sets, rownames(matrix))
    axis <- axes[dataset == cohort & axis_id == "fixed_projection"]
    meta <- merge(metadata[dataset == cohort], axis, by = c("sample_id", "dataset"),
                  all = FALSE, sort = FALSE)
    meta <- meta[match(intersect(meta$sample_id, colnames(matrix)), sample_id)]
    for (set_id in chosen) {
      values <- matrix[set_id, meta$sample_id]
      for (index in seq_along(centers)) {
        lower <- max(0, centers[[index]] - width / 2)
        upper <- min(1, centers[[index]] + width / 2)
        keep <- is.finite(meta$axis_percentile) & meta$axis_percentile >= lower &
          if (index == length(centers)) meta$axis_percentile <= upper else
            meta$axis_percentile < upper
        rows[[length(rows) + 1L]] <- data.table(
          dataset = cohort, set_id = set_id, window_id = index,
          center = centers[[index]], lower = lower, upper = upper,
          n_participants = sum(keep), mean_score = mean(values[keep], na.rm = TRUE),
          se_score = sd(values[keep], na.rm = TRUE) / sqrt(sum(keep)),
          visualization_only = TRUE
        )
      }
    }
  }
  rbindlist(rows, use.names = TRUE, fill = TRUE)
}

pathway_main <- function() {
  contract <- ml_read_contract()
  out_root <- ml_out_root(must_exist = TRUE)
  output <- file.path(out_root, "pathway_tf", "pathways")
  ml_assert(!dir.exists(output), paste0("Refusing to overwrite pathway output: ", output))
  ml_ensure_dir(output)

  dge_path <- ml_resolve(contract$dge_rds)
  manifest_path <- ml_resolve(contract$sample_manifest)
  annotation_path <- ml_resolve(contract$gene_annotation)
  pathway_dir <- ml_resolve(contract$pathway_dir)
  axes <- ml_load_axes(contract)
  gene_model_override <- Sys.getenv("HAC_ML_CONTINUUM_GENE_MODELS", unset = "")
  disease_override <- Sys.getenv("HAC_ML_DISEASE_GENE_STATS", unset = "")
  if (nzchar(gene_model_override) || nzchar(disease_override)) {
    ml_assert(nzchar(gene_model_override) && nzchar(disease_override),
              paste("Continuum and disease rank-statistic overrides must be supplied",
                    "together to avoid mixed estimands"))
    gene_model_path <- ml_resolve(gene_model_override)
    disease_path <- ml_resolve(disease_override)
    rank_statistic_source <- "explicit_environment_overrides"
  } else {
    rank_statistics <- pathway_generate_rank_statistics(
      dge_path, manifest_path, annotation_path, axes, contract
    )
    gene_model_path <- file.path(output, "continuum_gene_cohort_models.tsv.gz")
    disease_path <- file.path(output, "source_overlap_disease_gene_statistics.tsv.gz")
    ml_write_tsv_once(rank_statistics$continuum, gene_model_path)
    ml_write_tsv_once(rank_statistics$disease, disease_path)
    rank_statistic_source <- "generated_within_parallel_pathway_tf_workstream"
    rm(rank_statistics)
    gc()
  }
  input_manifest <- data.table(
    input_id = c("dge", "manifest", "annotation", "continuum_gene_models",
                 "disease_gene_statistics"),
    path = c(dge_path, manifest_path, annotation_path, gene_model_path, disease_path)
  )
  input_manifest[, sha256 := vapply(path, ml_sha256, character(1))]
  input_manifest[, rank_statistic_source := rank_statistic_source]
  ml_write_tsv_once(input_manifest, file.path(output, "input_manifest.tsv"))

  expression_data <- pathway_symbol_expression(dge_path, manifest_path, annotation_path)
  expression <- expression_data$expression
  metadata <- expression_data$metadata
  signature <- ml_signature_symbols(contract)
  continuum_cohort_stats <- pathway_normalize_gene_stats(
    gene_model_path, contract$co_primary_axes, contract$evaluation_cohorts
  )
  continuum_cohort_stats <- continuum_cohort_stats[!gene_symbol %in% signature]
  continuum_meta_stats <- pathway_meta_ranks(continuum_cohort_stats)
  disease_stats <- fread(disease_path)
  disease_symbol <- intersect(c("gene_symbol", "symbol", "gene_name"), names(disease_stats))[1L]
  disease_stat <- intersect(
    c("moderated_t", "t_value", "t", "z_value", "z", "statistic"),
    names(disease_stats)
  )[1L]
  ml_assert(!is.na(disease_symbol) && !is.na(disease_stat),
            "Disease-control statistics need gene-symbol and ranked-statistic columns")
  disease_ranks <- disease_stats[, .(
    gene_symbol = toupper(trimws(get(disease_symbol))),
    rank_statistic = as.numeric(get(disease_stat))
  )][!gene_symbol %in% signature & is.finite(rank_statistic)]

  collection_index <- list()
  for (collection in names(contract$pathway_collections)) {
    family_size <- as.integer(contract$pathway_collections[[collection]])
    collection_out <- file.path(output, collection)
    ml_ensure_dir(collection_out)
    gmt_path <- file.path(pathway_dir, paste0(collection, ".gmt"))
    gmt <- pathway_read_gmt(gmt_path, collection, family_size)
    scored <- pathway_score_collection(
      expression, metadata, gmt$ids, gmt$membership, signature,
      as.integer(contract$minimum_set_genes), contract$minimum_retained_fraction
    )
    scored$coverage[, `:=`(collection = collection, family_size = family_size)]
    ml_write_tsv_once(scored$coverage, file.path(collection_out, "testability.tsv"))

    score_long <- rbindlist(lapply(names(scored$scores), function(cohort) {
      matrix <- scored$scores[[cohort]]
      if (is.null(matrix)) return(NULL)
      answer <- as.data.table(as.table(matrix))
      setnames(answer, c("set_id", "sample_id", "pathway_score"))
      answer[, `:=`(set_id = as.character(set_id), sample_id = as.character(sample_id),
                    pathway_score = as.numeric(pathway_score), dataset = cohort)]
      answer
    }), use.names = TRUE, fill = TRUE)
    ml_write_tsv_once(score_long, file.path(collection_out, "donor_scores.tsv.gz"))

    fibrosis_models <- pathway_fit_all(
      scored$scores, scored$coverage, metadata, axes, gmt$ids,
      contract$evaluation_cohorts, contract$co_primary_axes,
      "fibrosis_stage", family_size
    )
    ml_write_tsv_once(fibrosis_models, file.path(collection_out, "fibrosis_models.tsv"))
    meta <- pathway_meta_models(
      fibrosis_models, gmt$ids, contract$co_primary_axes, family_size
    )
    ml_write_tsv_once(meta, file.path(collection_out, "meta_analysis.tsv"))
    shapes <- pathway_spline_shapes(
      fibrosis_models, scored$scores, metadata, axes, gmt$ids,
      contract$evaluation_cohorts, contract$co_primary_axes
    )
    ml_write_tsv_once(
      shapes$summary, file.path(collection_out, "nonlinear_shape_concordance.tsv")
    )
    ml_write_tsv_once(
      shapes$predictions,
      file.path(collection_out, "nonlinear_shape_predictions.tsv.gz")
    )

    nas_models <- pathway_fit_all(
      scored$scores, scored$coverage, metadata, axes, gmt$ids,
      "GSE162694", contract$co_primary_axes, "nas_score", family_size
    )
    ml_write_tsv_once(nas_models, file.path(collection_out, "nas_sensitivity.tsv"))

    continuum_enrichment <- rbindlist(lapply(contract$co_primary_axes, function(axis) {
      result <- pathway_ranked_enrichment(
        gmt$ids, gmt$membership,
        continuum_meta_stats[axis_id == axis], signature, family_size,
        as.integer(contract$seed)
      )
      result[, axis_id := axis]
      result
    }), use.names = TRUE, fill = TRUE)
    ml_write_tsv_once(
      continuum_enrichment,
      file.path(collection_out, "ranked_continuum_enrichment.tsv")
    )
    disease_enrichment <- pathway_ranked_enrichment(
      gmt$ids, gmt$membership, disease_ranks, signature, family_size,
      as.integer(contract$seed)
    )
    disease_enrichment[, axis_id := "source_overlap_disease_vs_control"]
    ml_write_tsv_once(
      disease_enrichment,
      file.path(collection_out, "ranked_disease_baseline_enrichment.tsv")
    )
    concordance <- merge(
      continuum_enrichment[, .(set_id, axis_id, continuum_nes = nes,
                                continuum_q = q_value)],
      disease_enrichment[, .(set_id, disease_nes = nes, disease_q = q_value)],
      by = "set_id", all.x = TRUE
    )
    concordance[, direction_concordant := is.finite(continuum_nes) &
                  is.finite(disease_nes) & sign(continuum_nes) == sign(disease_nes)]
    ml_write_tsv_once(concordance, file.path(collection_out, "disease_concordance.tsv"))

    if (collection == "hallmark") {
      missing_display <- setdiff(contract$display_hallmarks, gmt$ids)
      ml_assert(!length(missing_display), paste0(
        "Fixed Hallmark display roster missing: ", paste(missing_display, collapse = ";")
      ))
      windows <- pathway_window_summary(
        scored$scores, metadata, axes, contract$display_hallmarks,
        contract$evaluation_cohorts, as.numeric(contract$window_centers),
        as.numeric(contract$window_width)
      )
      ml_write_tsv_once(windows, file.path(collection_out, "fixed_display_windows.tsv"))
    }
    collection_index[[collection]] <- data.table(
      collection = collection, family_size = family_size,
      n_testable_gse162694 = scored$coverage[
        dataset == "GSE162694" & testable == TRUE, .N
      ],
      n_testable_gse213621 = scored$coverage[
        dataset == "GSE213621" & testable == TRUE, .N
      ],
      confirmatory = collection == "hallmark",
      interpretation = if (collection == "hallmark") "confirmatory" else "exploratory"
    )
    rm(scored, score_long, fibrosis_models, meta, shapes, nas_models,
       continuum_enrichment, disease_enrichment, concordance)
    gc()
  }
  ml_write_tsv_once(rbindlist(collection_index),
                    file.path(output, "collection_index.tsv"))
  ml_write_session_info(file.path(output, "sessionInfo.txt"))
  message("PATHWAY_CONTINUUM_COMPLETE: ", output)
}

if (!identical(Sys.getenv("HAC_ML_LIBRARY_ONLY", unset = "0"), "1")) {
  pathway_main()
}
