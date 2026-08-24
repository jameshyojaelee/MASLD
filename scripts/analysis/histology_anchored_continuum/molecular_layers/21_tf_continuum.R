#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(Matrix)
})

ml_tf_script_dir <- Sys.getenv("HAC_ML_SCRIPT_DIR", unset = "")
if (!nzchar(ml_tf_script_dir)) {
  ml_tf_args <- commandArgs(trailingOnly = FALSE)
  ml_tf_file <- sub("^--file=", "", ml_tf_args[grepl("^--file=", ml_tf_args)])
  stopifnot(length(ml_tf_file) == 1L)
  ml_tf_script_dir <- dirname(normalizePath(ml_tf_file))
}
source(file.path(ml_tf_script_dir, "lib_molecular_layers.R"))

tf_pick_input <- function(env_name, candidates, description) {
  override <- Sys.getenv(env_name, unset = "")
  if (nzchar(override)) return(ml_resolve(override))
  existing <- candidates[file.exists(candidates)]
  ml_assert(
    length(existing) == 1L,
    paste0(description, " requires exactly one upstream input; found ", length(existing),
           ". Set ", env_name, ". Candidates: ", paste(candidates, collapse = ";"))
  )
  normalizePath(existing, mustWork = TRUE)
}

tf_fixed_roster <- function(path, expected_size) {
  source <- fread(path)
  tf_column <- intersect(c("tf", "source"), names(source))[1L]
  ml_assert(!is.na(tf_column), "Frozen TF roster source has no TF column")
  roster <- sort(unique(as.character(source[[tf_column]])))
  ml_assert(length(roster) == expected_size,
            paste0("Frozen TF family-size drift: ", length(roster)))
  roster
}

tf_load_network <- function(roster) {
  ml_assert(requireNamespace("dorothea", quietly = TRUE),
            "dorothea is unavailable; no runtime network download is permitted")
  suppressPackageStartupMessages(library(dorothea))
  data("dorothea_hs", package = "dorothea", envir = environment())
  network <- as.data.table(get("dorothea_hs", envir = environment()))[
    confidence %in% c("A", "B", "C"),
    .(tf = as.character(tf), target = toupper(as.character(target)),
      mor = as.numeric(mor), confidence = as.character(confidence))
  ]
  network <- network[tf %in% roster & !is.na(target) & nzchar(target) & is.finite(mor)]
  network <- network[, .(
    mor = mean(mor), confidence = paste(sort(unique(confidence)), collapse = "/")
  ), by = .(tf, target)]
  ml_assert(setequal(unique(network$tf), roster),
            "The installed DoRothEA A/B/C network does not contain the frozen TF roster")
  network
}

tf_normalize_gene_stats <- function(path, axes, cohorts) {
  x <- fread(path)
  symbol_column <- intersect(c("gene_symbol", "symbol", "gene_name"), names(x))[1L]
  axis_column <- intersect(c("axis_id", "continuum_axis", "score_id"), names(x))[1L]
  dataset_column <- intersect(c("dataset", "cohort"), names(x))[1L]
  ml_assert(!is.na(symbol_column) && !is.na(axis_column) && !is.na(dataset_column),
            "TF input requires gene-symbol, continuum-axis, and cohort columns")
  if (symbol_column != "gene_symbol") setnames(x, symbol_column, "gene_symbol")
  if (axis_column != "axis_id") setnames(x, axis_column, "axis_id")
  if (dataset_column != "dataset") setnames(x, dataset_column, "dataset")
  statistic_column <- intersect(
    c("moderated_t", "t_value", "t", "z_value", "z", "statistic"), names(x)
  )[1L]
  if (!is.na(statistic_column)) {
    x[, rank_statistic := as.numeric(get(statistic_column))]
  } else {
    ml_assert(all(c("beta", "se") %in% names(x)),
              "TF input requires a moderated statistic or beta/se")
    x[, rank_statistic := beta / se]
  }
  x[, gene_symbol := toupper(trimws(gene_symbol))]
  x[axis_id %in% axes & dataset %in% cohorts & is.finite(rank_statistic) &
      !is.na(gene_symbol) & nzchar(gene_symbol)]
}

tf_run_wmean_family <- function(statistics, network, roster, signature,
                                family_size, permutations, seed) {
  ml_assert(requireNamespace("decoupleR", quietly = TRUE),
            "decoupleR is unavailable; frozen run_wmean analysis cannot run")
  statistics <- statistics[!gene_symbol %in% signature]
  statistics <- statistics[order(-abs(rank_statistic))][!duplicated(gene_symbol)]
  matrix <- matrix(
    statistics$rank_statistic, ncol = 1L,
    dimnames = list(statistics$gene_symbol, "continuum")
  )
  available_network <- network[!target %in% signature & target %in% rownames(matrix)]
  set.seed(seed)
  result <- as.data.table(decoupleR::run_wmean(
    mat = matrix, net = as.data.frame(available_network),
    .source = "tf", .target = "target", .mor = "mor",
    times = permutations, minsize = 5L
  ))
  result <- result[statistic == "norm_wmean"]
  if ("source" %in% names(result)) setnames(result, "source", "tf")
  base <- data.table(tf = roster)
  result <- merge(base, result, by = "tf", all.x = TRUE, sort = FALSE)
  result[, `:=`(
    q_value = ml_complete_bh(p_value, family_size),
    testable = is.finite(score) & is.finite(p_value),
    testability_reason = fifelse(is.finite(score) & is.finite(p_value), "testable",
                                 "run_wmean_minsize_or_statistic_coverage_failure"),
    bh_family_size = family_size,
    seed = seed,
    permutations = permutations,
    confidence_levels = "A/B/C",
    signature_genes_excluded = TRUE
  )]
  result
}

tf_symbol_expression <- function(dge_path, manifest_path, annotation_path) {
  dge <- readRDS(dge_path)
  metadata <- fread(manifest_path, na.strings = c("", "NA"))
  ml_assert(all(c("sample_id", "analysis_unit_id", "dataset", "inferred_sex",
                  "fibrosis_stage", "nas_score") %in% names(metadata)),
            "Sample manifest schema drift")
  ml_assert(!anyDuplicated(metadata$analysis_unit_id),
            "Cells or technical rows cannot be donor replicates")
  counts <- dge$counts
  ml_assert(!is.null(counts) && setequal(colnames(counts), metadata$sample_id),
            "DGE/manifest sample sets differ")
  log_cpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
  log_cpm <- log_cpm[, metadata$sample_id, drop = FALSE]
  annotation <- fread(annotation_path, select = c("gene_id", "gene_name"))
  annotation[, `:=`(
    gene_id_base = ml_base_gene_id(gene_id),
    gene_symbol = toupper(trimws(gene_name))
  )]
  annotation <- unique(annotation[
    !is.na(gene_symbol) & nzchar(gene_symbol), .(gene_id_base, gene_symbol)
  ])
  ml_assert(annotation[, uniqueN(gene_symbol), by = gene_id_base][V1 > 1L, .N] == 0L,
            "One Ensembl base ID maps to multiple symbols")
  map <- annotation[match(ml_base_gene_id(rownames(log_cpm)), gene_id_base)]
  keep <- !is.na(map$gene_symbol) & nzchar(map$gene_symbol)
  expression <- rowsum(log_cpm[keep, , drop = FALSE], map$gene_symbol[keep],
                       reorder = FALSE, na.rm = TRUE)
  expression <- expression /
    as.numeric(table(map$gene_symbol[keep])[rownames(expression)])
  storage.mode(expression) <- "double"
  list(expression = expression, metadata = metadata)
}

tf_zscore_rows <- function(matrix) {
  center <- rowMeans(matrix, na.rm = TRUE)
  spread <- apply(matrix, 1L, sd, na.rm = TRUE)
  answer <- sweep(matrix, 1L, center, "-")
  variable <- is.finite(spread) & spread > 0
  answer[variable, ] <- sweep(answer[variable, , drop = FALSE], 1L,
                             spread[variable], "/")
  answer[!variable, ] <- NA_real_
  answer
}

tf_score_regulons <- function(expression, metadata, network, roster, signature,
                              minimum_fraction = 0.8, minimum_targets = 5L) {
  original <- network[, .(
    original_l1 = sum(abs(mor)), n_original_targets = uniqueN(target)
  ), by = tf]
  excluded <- network[target %in% signature, .(
    excluded_l1 = sum(abs(mor)), n_signature_targets = uniqueN(target),
    signature_targets = paste(sort(unique(target)), collapse = ";")
  ), by = tf]
  retained <- network[!target %in% signature]
  coverage_rows <- list()
  score_blocks <- list()
  for (cohort in unique(metadata$dataset)) {
    samples <- metadata[dataset == cohort, sample_id]
    z <- tf_zscore_rows(expression[, samples, drop = FALSE])
    observed_genes <- rownames(z)[rowSums(is.finite(z)) == ncol(z)]
    available <- retained[target %in% observed_genes]
    observed <- available[, .(
      observed_l1 = sum(abs(mor)), n_observed_targets = uniqueN(target)
    ), by = tf]
    coverage <- Reduce(function(x, y) merge(x, y, by = "tf", all.x = TRUE), list(
      data.table(tf = roster), original, excluded, observed
    ))
    for (column in c("excluded_l1", "n_signature_targets", "observed_l1",
                     "n_observed_targets")) {
      set(coverage, which(is.na(coverage[[column]])), column, 0)
    }
    coverage[is.na(signature_targets), signature_targets := ""]
    coverage[, `:=`(
      dataset = cohort,
      retained_l1_fraction = (original_l1 - excluded_l1) / original_l1,
      observed_l1_fraction = observed_l1 / original_l1
    )]
    coverage[, testable := n_observed_targets >= minimum_targets &
               observed_l1_fraction >= minimum_fraction]
    coverage[, testability_reason := fifelse(
      testable, "testable",
      fifelse(n_observed_targets < minimum_targets, "fewer_than_five_observed_targets",
              "less_than_80_percent_original_l1_observed_after_signature_exclusion")
    )]
    coverage_rows[[cohort]] <- coverage
    testable_tfs <- coverage[testable == TRUE, tf]
    if (!length(testable_tfs)) next
    selected <- available[tf %in% testable_tfs]
    gene_levels <- observed_genes
    tf_index <- match(selected$tf, testable_tfs)
    target_index <- match(selected$target, gene_levels)
    denominator <- coverage[match(testable_tfs, tf), observed_l1]
    ml_assert(all(is.finite(denominator) & denominator > 0),
              paste0(cohort, " regulon denominators are invalid"))
    weights <- Matrix::sparseMatrix(
      i = tf_index, j = target_index, x = selected$mor / denominator[tf_index],
      dims = c(length(testable_tfs), length(gene_levels)),
      dimnames = list(testable_tfs, gene_levels)
    )
    scores <- as.matrix(weights %*% z[gene_levels, , drop = FALSE])
    scores <- t(apply(scores, 1L, ml_standardize))
    if (length(testable_tfs) == 1L) {
      scores <- matrix(scores, nrow = 1L, dimnames = list(testable_tfs, samples))
    } else {
      rownames(scores) <- testable_tfs
      colnames(scores) <- samples
    }
    ml_assert(any(is.finite(scores)), paste0(cohort, " regulon scores are all non-finite"))
    score_blocks[[cohort]] <- scores
  }
  list(coverage = rbindlist(coverage_rows, use.names = TRUE, fill = TRUE),
       scores = score_blocks)
}

tf_fit_outcome_nonlinear <- function(data, stage_column) {
  d <- copy(data)
  d <- d[is.finite(outcome_z) & is.finite(axis_raw) &
           !is.na(get(stage_column)) & !is.na(inferred_sex)]
  if (nrow(d) < 20L || uniqueN(d[[stage_column]]) < 2L ||
      uniqueN(d$inferred_sex) < 2L) {
    return(data.table(
      estimable = FALSE, n = nrow(d), beta = NA_real_, se = NA_real_,
      p_value = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
      nonlinear_p_value = NA_real_, spline_df = 3L
    ))
  }
  d[, `:=`(
    axis_z = ml_standardize(axis_raw),
    stage_factor = factor(get(stage_column)),
    sex_factor = factor(inferred_sex)
  )]
  linear <- lm(outcome_z ~ axis_z + stage_factor + sex_factor, data = d)
  nonlinear <- lm(
    outcome_z ~ splines::ns(axis_z, df = 3) + stage_factor + sex_factor,
    data = d
  )
  coefficient <- coef(summary(linear))["axis_z", ]
  comparison <- anova(linear, nonlinear)
  data.table(
    estimable = TRUE, n = nrow(d), beta = coefficient[["Estimate"]],
    se = coefficient[["Std. Error"]], p_value = coefficient[["Pr(>|t|)"]],
    ci_low = coefficient[["Estimate"]] - qt(0.975, df.residual(linear)) *
      coefficient[["Std. Error"]],
    ci_high = coefficient[["Estimate"]] + qt(0.975, df.residual(linear)) *
      coefficient[["Std. Error"]],
    nonlinear_p_value = comparison$`Pr(>F)`[[2L]], spline_df = 3L
  )
}

tf_fit_regulon_models <- function(scores, coverage, metadata, axes, roster,
                                  cohorts, axis_names, stage_column, family_size) {
  rows <- list()
  for (cohort in cohorts) {
    matrix <- scores[[cohort]]
    for (axis_name in axis_names) {
      axis <- axes[dataset == cohort & axis_id == axis_name]
      joined <- merge(metadata[dataset == cohort], axis, by = c("sample_id", "dataset"),
                      all = FALSE, sort = FALSE)
      base <- merge(
        data.table(tf = roster),
        coverage[dataset == cohort, .(tf, testable, testability_reason)],
        by = "tf", all.x = TRUE, sort = FALSE
      )
      fits <- list()
      if (!is.null(matrix)) {
        common <- intersect(joined$sample_id, colnames(matrix))
        joined <- joined[match(common, sample_id)]
        for (tf_name in intersect(rownames(matrix), roster)) {
          d <- copy(joined)
          d[, outcome_z := as.numeric(matrix[tf_name, sample_id])]
          fit <- tf_fit_outcome_nonlinear(d, stage_column)
          fit[, tf := tf_name]
          fits[[tf_name]] <- fit
        }
      }
      fit_table <- rbindlist(fits, use.names = TRUE, fill = TRUE)
      if (!nrow(fit_table)) {
        fit_table <- data.table(
          tf = character(), estimable = logical(), n = integer(), beta = numeric(),
          se = numeric(), p_value = numeric(), ci_low = numeric(), ci_high = numeric(),
          nonlinear_p_value = numeric(), spline_df = integer()
        )
      }
      result <- merge(base, fit_table, by = "tf", all.x = TRUE, sort = FALSE)
      result[is.na(estimable), estimable := FALSE]
      result[, `:=`(
        dataset = cohort, axis_id = axis_name,
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

tf_meta_models <- function(models, roster, axes, family_size) {
  result <- rbindlist(lapply(axes, function(axis_name) {
    rbindlist(lapply(roster, function(tf_name) {
      d <- models[models$axis_id == axis_name & models$tf == tf_name]
      fixed <- ml_fixed_meta(d$beta, d$se)
      random <- ml_random_meta(d$beta, d$se)
      data.table(
        tf = tf_name, axis_id = axis_name,
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

tf_fixed_windows <- function(scores, metadata, axes, display_tfs, cohorts,
                             centers, width) {
  rows <- list()
  for (cohort in cohorts) {
    matrix <- scores[[cohort]]
    if (is.null(matrix)) next
    selected <- intersect(display_tfs, rownames(matrix))
    axis <- axes[dataset == cohort & axis_id == "fixed_projection"]
    joined <- merge(metadata[dataset == cohort], axis, by = c("sample_id", "dataset"),
                    all = FALSE, sort = FALSE)
    joined <- joined[match(intersect(joined$sample_id, colnames(matrix)), sample_id)]
    for (tf_name in selected) {
      values <- matrix[tf_name, joined$sample_id]
      for (index in seq_along(centers)) {
        lower <- max(0, centers[[index]] - width / 2)
        upper <- min(1, centers[[index]] + width / 2)
        keep <- is.finite(joined$axis_percentile) & joined$axis_percentile >= lower &
          if (index == length(centers)) joined$axis_percentile <= upper else
            joined$axis_percentile < upper
        rows[[length(rows) + 1L]] <- data.table(
          dataset = cohort, tf = tf_name, window_id = index,
          center = centers[[index]], lower = lower, upper = upper,
          n_participants = sum(keep), mean_score = mean(values[keep], na.rm = TRUE),
          se_score = sd(values[keep], na.rm = TRUE) / sqrt(sum(keep)),
          score_type = "signed_regulon_score", visualization_only = TRUE
        )
      }
    }
  }
  rbindlist(rows, use.names = TRUE, fill = TRUE)
}

tf_main <- function() {
  contract <- ml_read_contract()
  out_root <- ml_out_root(must_exist = TRUE)
  output <- file.path(out_root, "pathway_tf", "tf")
  ml_assert(!dir.exists(output), paste0("Refusing to overwrite TF output: ", output))
  ml_ensure_dir(output)

  gene_model_override <- Sys.getenv("HAC_ML_CONTINUUM_GENE_MODELS", unset = "")
  pathway_local_statistics <- file.path(
    out_root, "pathway_tf", "pathways", "continuum_gene_cohort_models.tsv.gz"
  )
  gene_model_path <- if (nzchar(gene_model_override)) {
    ml_resolve(gene_model_override)
  } else if (file.exists(pathway_local_statistics)) {
    normalizePath(pathway_local_statistics, mustWork = TRUE)
  } else {
    tf_pick_input(
      "HAC_ML_CONTINUUM_GENE_MODELS",
      file.path(out_root, c(
        "bulk/continuum_gene_cohort_models.tsv.gz",
        "bulk/transcript_cohort_models.tsv.gz",
        "bulk/continuum_gene_models.tsv.gz"
      )),
      "Continuum TF activity"
    )
  }
  roster_path <- tf_pick_input(
    "HAC_ML_TF_ROSTER_SOURCE",
    file.path(
      ml_project_root(),
      "figures/candidates/pi-figure-redesign-2026-08-17-v13/source_tables/fig4f_all_tf_activity.tsv"
    ),
    "Frozen 268-TF roster"
  )
  dge_path <- ml_resolve(contract$dge_rds)
  manifest_path <- ml_resolve(contract$sample_manifest)
  annotation_path <- ml_resolve(contract$gene_annotation)
  input_manifest <- data.table(
    input_id = c("continuum_gene_models", "tf_roster_source", "dge", "manifest",
                 "annotation"),
    path = c(gene_model_path, roster_path, dge_path, manifest_path, annotation_path)
  )
  input_manifest[, sha256 := vapply(path, ml_sha256, character(1))]
  ml_write_tsv_once(input_manifest, file.path(output, "input_manifest.tsv"))

  family_size <- as.integer(contract$tf_family_size)
  roster <- tf_fixed_roster(roster_path, family_size)
  network <- tf_load_network(roster)
  ml_write_tsv_once(network, file.path(output, "frozen_dorothea_abc_network.tsv.gz"))
  signature <- ml_signature_symbols(contract)
  network_audit <- network[, .(
    n_targets = uniqueN(target), original_l1 = sum(abs(mor)),
    n_signature_targets = uniqueN(target[target %in% signature]),
    signature_excluded_l1 = sum(abs(mor)[target %in% signature]),
    retained_l1_fraction = sum(abs(mor)[!target %in% signature]) / sum(abs(mor))
  ), by = tf]
  network_audit[, `:=`(
    confidence_levels = "A/B/C", family_size = family_size,
    network_source = "installed_dorothea_hs_no_runtime_download"
  )]
  ml_write_tsv_once(network_audit, file.path(output, "network_audit.tsv"))

  statistics <- tf_normalize_gene_stats(
    gene_model_path, contract$co_primary_axes, contract$evaluation_cohorts
  )
  activity <- rbindlist(lapply(contract$evaluation_cohorts, function(cohort) {
    rbindlist(lapply(contract$co_primary_axes, function(axis_name) {
      result <- tf_run_wmean_family(
        statistics[dataset == cohort & axis_id == axis_name], network, roster,
        signature, family_size, 1000L, as.integer(contract$tf_seed)
      )
      result[, `:=`(dataset = cohort, axis_id = axis_name)]
      result
    }), use.names = TRUE, fill = TRUE)
  }), use.names = TRUE, fill = TRUE)
  ml_assert(all(activity[, .N, by = .(dataset, axis_id)]$N == family_size),
            "TF activity output does not preserve all 268 rows per cohort/axis")
  ml_write_tsv_once(activity, file.path(output, "continuum_tf_activity.tsv"))
  cross_cohort <- activity[, .(
    n_testable = sum(testable),
    direction_concordant = sum(testable) == length(contract$evaluation_cohorts) &&
      (all(score[testable] > 0) || all(score[testable] < 0)),
    maximum_q = if (all(testable)) max(q_value) else NA_real_
  ), by = .(axis_id, tf)]
  expected_tf_tests <- length(contract$evaluation_cohorts) *
    length(contract$co_primary_axes)
  cross_score <- activity[, .(
    all_four_direction_concordant = sum(testable) == expected_tf_tests &&
      (all(score[testable] > 0) || all(score[testable] < 0)),
    all_four_q_below_005 = sum(testable) == expected_tf_tests &&
      all(q_value[testable] < 0.05)
  ), by = tf]
  cross_cohort <- merge(cross_cohort, cross_score, by = "tf", all.x = TRUE)
  cross_cohort[, display_tf := tf %in% contract$display_tfs]
  ml_write_tsv_once(cross_cohort, file.path(output, "tf_cross_cohort_concordance.tsv"))

  expression_data <- tf_symbol_expression(dge_path, manifest_path, annotation_path)
  regulons <- tf_score_regulons(
    expression_data$expression, expression_data$metadata, network, roster, signature,
    minimum_fraction = contract$minimum_retained_fraction, minimum_targets = 5L
  )
  regulons$coverage[, family_size := family_size]
  ml_write_tsv_once(regulons$coverage, file.path(output, "regulon_testability.tsv"))
  score_long <- rbindlist(lapply(names(regulons$scores), function(cohort) {
    matrix <- regulons$scores[[cohort]]
    if (is.null(matrix)) return(NULL)
    answer <- as.data.table(as.table(matrix))
    setnames(answer, c("tf", "sample_id", "regulon_score"))
    answer[, `:=`(tf = as.character(tf), sample_id = as.character(sample_id),
                  regulon_score = as.numeric(regulon_score), dataset = cohort,
                  score_type = "signed_regulon_score")]
    answer
  }), use.names = TRUE, fill = TRUE)
  ml_write_tsv_once(score_long, file.path(output, "donor_regulon_scores.tsv.gz"))

  axes <- ml_load_axes(contract)
  models <- tf_fit_regulon_models(
    regulons$scores, regulons$coverage, expression_data$metadata, axes, roster,
    contract$evaluation_cohorts, contract$co_primary_axes, "fibrosis_stage", family_size
  )
  ml_write_tsv_once(models, file.path(output, "regulon_fibrosis_models.tsv"))
  meta <- tf_meta_models(models, roster, contract$co_primary_axes, family_size)
  ml_write_tsv_once(meta, file.path(output, "regulon_meta_analysis.tsv"))
  nas <- tf_fit_regulon_models(
    regulons$scores, regulons$coverage, expression_data$metadata, axes, roster,
    "GSE162694", contract$co_primary_axes, "nas_score", family_size
  )
  ml_write_tsv_once(nas, file.path(output, "regulon_nas_sensitivity.tsv"))
  windows <- tf_fixed_windows(
    regulons$scores, expression_data$metadata, axes, contract$display_tfs,
    contract$evaluation_cohorts, as.numeric(contract$window_centers),
    as.numeric(contract$window_width)
  )
  ml_write_tsv_once(windows, file.path(output, "fixed_display_regulon_windows.tsv"))
  ml_write_session_info(file.path(output, "sessionInfo.txt"))
  message("TF_CONTINUUM_COMPLETE: ", output)
}

if (!identical(Sys.getenv("HAC_ML_LIBRARY_ONLY", unset = "0"), "1")) {
  tf_main()
}
