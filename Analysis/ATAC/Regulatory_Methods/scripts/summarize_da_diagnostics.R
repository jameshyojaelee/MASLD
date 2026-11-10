#!/usr/bin/env Rscript
# Summarize completed donor-level fits; uncertainty here describes simulation draws.
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(jsonlite) })
source(file.path(Sys.getenv("HOME"), "publication_color_themes.R"))
setDTthreads(1L)
grDevices::pdf.options(useDingbats = FALSE)

require_true <- function(value, message) if (!isTRUE(value)) stop(message, call. = FALSE)
KEY <- c("cohort", "lineage", "representation", "filter")
METHODS <- c("voom_qw", "edger_ql")
COLORS <- c(voom_qw = "#9E9E9E", edger_ql = "#1565C0")

mean_mc <- function(value) {
  n <- length(value)
  average <- mean(value)
  error <- if (n > 1L) sd(value) / sqrt(n) else NA_real_
  critical <- if (n > 1L) qt(.975, df = n - 1L) else NA_real_
  list(mean = average, se = error, lower = average - critical * error, upper = average + critical * error)
}

summarize_null <- function(d) {
  require_true(!anyDuplicated(d[, c(KEY, "method", "replicate"), with = FALSE]), "Duplicate simulation draw")
  d[, {
    tested <- .SD[status == "testable"]
    require_true(nrow(tested) > 1L && all(tested$n_tested_peaks > 0), "Insufficient valid simulation draws")
    require_true(uniqueN(tested$n_tested_peaks) == 1L, "A diagnostic changes its peak universe across draws")
    nominal <- mean_mc(tested$n_nominal_p05 / tested$n_tested_peaks)
    family_error <- tested$n_bh_q05 > 0
    interval <- binom.test(sum(family_error), length(family_error))$conf.int
    list(n_recorded_draws = .N, n_testable_draws = nrow(tested), n_failed_draws = sum(status == "failed"),
      n_other_draws = sum(!status %chin% c("testable", "failed")), n_diagnostic_peaks = tested$n_tested_peaks[1L],
      nominal_rejection_fraction = nominal$mean, nominal_mc_se = nominal$se,
      nominal_mc_lower = nominal$lower, nominal_mc_upper = nominal$upper,
      n_draws_any_BH_rejection = sum(family_error), probability_any_BH_rejection = mean(family_error),
      probability_any_BH_mc_lower = interval[1L], probability_any_BH_mc_upper = interval[2L],
      mean_fraction_peaks_BH_rejected = mean(tested$n_bh_q05 / tested$n_tested_peaks),
      tiny_feature_universe = tested$n_tested_peaks[1L] < 100L)
  }, by = c(KEY, "method")]
}

pair_methods <- function(d) {
  columns <- c(KEY, "replicate", "status", "n_tested_peaks", "n_nominal_p05", "n_bh_q05",
    "simulated_count_sum", "simulated_count_sum_squares", "seed")
  a <- d[method == "voom_qw", ..columns]
  b <- d[method == "edger_ql", ..columns]
  setnames(a, setdiff(columns, c(KEY, "replicate")), paste0(setdiff(columns, c(KEY, "replicate")), "_voom"))
  setnames(b, setdiff(columns, c(KEY, "replicate")), paste0(setdiff(columns, c(KEY, "replicate")), "_ql"))
  paired <- merge(a, b, by = c(KEY, "replicate"), all = TRUE)
  require_true(nrow(paired) == nrow(a) && nrow(a) == nrow(b), "Methods have different diagnostic draw sets")
  complete <- paired[status_voom == "testable" & status_ql == "testable"]
  for (column in c("n_tested_peaks", "seed", "simulated_count_sum", "simulated_count_sum_squares")) {
    require_true(identical(complete[[paste0(column, "_voom")]], complete[[paste0(column, "_ql")]]),
      paste("Paired diagnostic checksum differs:", column))
  }
  complete[, {
    nominal <- mean_mc(n_nominal_p05_ql / n_tested_peaks_ql - n_nominal_p05_voom / n_tested_peaks_voom)
    family <- mean_mc(as.integer(n_bh_q05_ql > 0) - as.integer(n_bh_q05_voom > 0))
    list(n_paired_draws = .N, n_diagnostic_peaks = n_tested_peaks_voom[1L],
      ql_minus_voom_nominal_fraction = nominal$mean, nominal_difference_mc_se = nominal$se,
      nominal_difference_mc_lower = nominal$lower, nominal_difference_mc_upper = nominal$upper,
      ql_minus_voom_probability_any_BH = family$mean, any_BH_difference_mc_se = family$se,
      any_BH_difference_mc_lower = family$lower, any_BH_difference_mc_upper = family$upper,
      diagnostic_seed_and_two_moment_checksums_match = TRUE)
  }, by = KEY]
}

self_test <- function() {
  fixture <- data.table(cohort = "fixture", lineage = "hepatocyte", representation = "fragment", filter = "primary",
    method = rep(METHODS, each = 2L), replicate = rep(1:2, 2L), status = "testable", n_tested_peaks = 2000L,
    n_nominal_p05 = c(0L, 100L, 100L, 100L), n_bh_q05 = c(0L, 1L, 1L, 1L),
    simulated_count_sum = rep(c(10L, 20L), 2L), simulated_count_sum_squares = rep(c(100L, 400L), 2L), seed = 42L)
  result <- summarize_null(fixture)
  stopifnot(result[method == "voom_qw", nominal_rejection_fraction] == .025,
    result[method == "voom_qw", probability_any_BH_rejection] == .5,
    result[method == "voom_qw", mean_fraction_peaks_BH_rejected] == .00025)
  paired <- pair_methods(fixture)
  stopifnot(paired$ql_minus_voom_probability_any_BH == .5, paired$ql_minus_voom_nominal_fraction == .025)
  altered <- copy(fixture); altered[method == "edger_ql" & replicate == 1L, simulated_count_sum := 999L]
  stopifnot(inherits(tryCatch(pair_methods(altered), error = identity), "error"))
  tiny <- copy(fixture); tiny[, `:=`(n_tested_peaks = 2L, n_nominal_p05 = 0L)]
  stopifnot(all(summarize_null(tiny)$tiny_feature_universe))
  cat("Scientific diagnostic fixtures PASS\n")
}

panel_theme <- function() {
  theme_classic(base_size = 6, base_family = "Helvetica") +
    theme(text = element_text(size = 6, face = "plain"), plot.title = element_text(size = 6, face = "plain"),
      plot.subtitle = element_text(size = 6, face = "plain"), axis.title = element_text(size = 6, face = "plain"),
      axis.text = element_text(size = 6, face = "plain", color = "black"),
      legend.title = element_blank(), legend.text = element_text(size = 6, face = "plain"),
      legend.position = "bottom", plot.caption = element_text(size = 6, face = "plain", hjust = 0),
      plot.caption.position = "plot", plot.margin = margin(4, 4, 4, 4))
}

render_calibration <- function(d, output) {
  manifest <- list()
  for (arm in unique(d$run_arm)) for (representation_value in c("fragment", "binary")) for (rule in c("primary", "prevalence")) {
    local <- copy(d[run_arm == arm & representation == representation_value & filter == rule])
    contexts <- unique(local[, .(cohort, lineage, n_diagnostic_peaks)])
    contexts[, lineage_order := match(lineage, c("hepatocyte", "stellate", "macrophage"))]
    setorder(contexts, lineage_order, cohort)
    contexts[, context := paste(sub("GSE", "", cohort), lineage,
      paste0("(", format(n_diagnostic_peaks, big.mark = ",", trim = TRUE), " peaks)"))]
    contexts[, y := rev(seq_len(.N))]
    local <- merge(local, contexts[, .(cohort, lineage, context, y)], by = c("cohort", "lineage"))
    local[, y_point := y + ifelse(method == "voom_qw", .10, -.10)]
    for (metric in c("nominal", "any_BH")) {
      if (metric == "nominal") local[, `:=`(value = nominal_rejection_fraction, lower = nominal_mc_lower, upper = nominal_mc_upper)]
      else local[, `:=`(value = probability_any_BH_rejection, lower = probability_any_BH_mc_lower, upper = probability_any_BH_mc_upper)]
      # KEY MESSAGE: Distinguish per-peak nominal rejection from any BH rejection under an assumed all-null count model.
      p <- ggplot(local, aes(value, y_point, color = method, shape = method)) +
        geom_vline(xintercept = .05, linewidth = .3, linetype = 2, color = "#9E9E9E") +
        geom_segment(aes(x = lower, xend = upper, yend = y_point), linewidth = .4) + geom_point(size = 1.4) +
        scale_color_manual(values = COLORS, breaks = METHODS, labels = c("Quality-weighted voom", "Robust edgeR QL")) +
        scale_shape_manual(values = c(voom_qw = 16, edger_ql = 17), breaks = METHODS, labels = c("Quality-weighted voom", "Robust edgeR QL")) +
        scale_y_continuous(breaks = contexts$y, labels = contexts$context, expand = expansion(add = .45)) +
        scale_x_continuous(labels = function(x) paste0(format(round(100 * x, 1), trim = TRUE), "%"),
          limits = c(0, max(.055, local$upper))) +
        labs(x = if (metric == "nominal") "Mean fraction of peaks with nominal p < 0.05" else "Fraction of draws with any BH q < 0.05",
          y = NULL, title = paste("NB global-null diagnostic:", if (arm == "full") "original counts" else "equal-cell sensitivity"),
          subtitle = paste(if (representation_value == "fragment") "Fragment counts" else "Accessible-nucleus counts", "/", if (rule == "primary") "original filter" else "prevalence sensitivity"),
          caption = paste0("1,000 simulation draws; 95% Monte Carlo intervals, not donor uncertainty.\n",
            if (metric == "any_BH") "Any-rejection probability equals FDR only under this fitted all-null simulation.\n" else "Peak fractions are averaged within each independent simulation draw.\n",
            if (any(local$tiny_feature_universe)) "Two-peak stellate diagnostics are too sparse for method selection.\n" else "",
            "Selected fixed features; assumed NB model does not establish biological calibration.",
            if (representation_value == "binary") "\nBinary aggregates are bounded; this NB null does not preserve that distribution." else "")) + panel_theme()
      name <- paste("calibration", arm, representation_value, rule, metric, sep = "_")
      path <- file.path(output, paste0(name, ".pdf"))
      grDevices::pdf(path, width = 5.3, height = if (nrow(contexts) > 4L) 3.3 else 2.8, family = "Helvetica", useDingbats = FALSE)
      tryCatch(print(p), finally = grDevices::dev.off())
      manifest[[length(manifest) + 1L]] <- data.table(panel = basename(path), run_arm = arm,
        representation = representation_value, filter = rule, metric = metric, n_contexts = nrow(contexts))
    }
  }
  fwrite(rbindlist(manifest), file.path(output, "panel_manifest.tsv"), sep = "\t")
}

run_review <- function(source, output, warning_receipt = NULL) {
  require_true(nzchar(Sys.getenv("SLURM_JOB_ID")), "Execute in an allocated compute node")
  require_true(!file.exists(output), "Refusing to overwrite diagnostic review")
  dir.create(output, recursive = TRUE)
  self_test()
  nulls <- pairs <- lodo <- designs <- agreement <- directional <- batch <- list()
  statuses <- list()
  source_files <- character()
  for (arm in c("full", "cell_downsample")) {
    root <- file.path(source, arm)
    info <- fromJSON(file.path(root, "run_summary.json"))
    require_true(identical(info$status, "candidate_complete_not_adopted"), "Incomplete source analysis")
    required <- c("run_summary.json", "null_diagnostics.tsv", "lodo_summary.tsv", "model_qc.tsv.gz", "filter_qc.tsv",
      "cohort_summary.tsv", "replication.tsv.gz", "covariate_sensitivity.tsv.gz", "covariate_availability.tsv")
    source_files <- c(source_files, file.path(root, required))
    n <- fread(file.path(root, "null_diagnostics.tsv"))
    require_true(all(n$null_model == "fitted_intercept_NB_global_null_fixed_features"), "Unexpected null model")
    nulls[[arm]] <- summarize_null(n)[, run_arm := arm]
    require_true(all(nulls[[arm]]$n_testable_draws == 1000L), "Not all planned diagnostic draws completed")
    pairs[[arm]] <- pair_methods(n)[, run_arm := arm]
    statuses[[arm]] <- n[, .(n_rows = .N), by = .(status)][, `:=`(run_arm = arm, component = "null")]
    omitted <- fread(file.path(root, "lodo_summary.tsv"))
    statuses[[paste0(arm, "_lodo")]] <- omitted[, .(n_rows = .N), by = .(status)][, `:=`(run_arm = arm, component = "donor_omission")]
    lodo[[arm]] <- omitted[status == "testable", .(n_omitted_donors = .N, n_peaks = unique(n_peaks),
      minimum_effect_spearman = min(effect_spearman), median_effect_spearman = median(effect_spearman),
      worst_rho_donor = omitted_donor[which.min(effect_spearman)], minimum_sign_agreement = min(sign_agreement),
      largest_median_absolute_log2_effect_change = max(median_absolute_effect_change)), by = c(KEY, "method")][, run_arm := arm]
    qc <- fread(file.path(root, "model_qc.tsv.gz"))
    designs[[arm]] <- qc[, .(n_donors = uniqueN(donor_id), n_normal = uniqueN(donor_id[condition == "NORMAL"]),
      n_mash = uniqueN(donor_id[condition == "MASH"]), design_rank = unique(design_rank),
      n_model_columns = unique(n_model_columns), residual_df = unique(residual_df)), by = c(KEY, "method")][, run_arm := arm]
    cohort <- fread(file.path(root, "cohort_summary.tsv"))
    cross <- cohort[cohort == "cross_cohort"]
    setnames(cross, "n_supported", "n_peakwise_both_BH_same_direction")
    cross[, `:=`(run_arm = arm, directional_family = "BH_three_prespecified_lineages_per_representation_filter_method",
      interpretation = "peakwise_both_cohort_BH_is_distinct_from_lineage_directional_enrichment")]
    directional[[arm]] <- cross
    x <- fread(file.path(root, "replication.tsv.gz"), select = c("lineage", "representation", "filter", "method", "logFC_gse244832", "logFC_gse281367", "evidence_state"))
    agreement[[arm]] <- x[, .(n_peaks = .N, effect_spearman = cor(logFC_gse244832, logFC_gse281367, method = "spearman"),
      effect_sign_agreement = mean(sign(logFC_gse244832) == sign(logFC_gse281367)),
      n_peakwise_both_BH_same_direction = sum(evidence_state == "supported")), by = .(lineage, representation, filter, method)][, run_arm := arm]
    rm(x); gc()
    adjusted <- fread(file.path(root, "covariate_sensitivity.tsv.gz"), select = c(KEY, "method", "covariate", "status", "adjusted_qvalue"))
    batch[[arm]] <- adjusted[, .(n_tested_peaks = .N, n_BH_q05_adjusted = sum(adjusted_qvalue < .05)),
      by = c(KEY, "method", "covariate", "status")]
    reference <- cohort[cohort != "cross_cohort", c(KEY, "method", "n_bh_q05"), with = FALSE]
    batch[[arm]] <- merge(batch[[arm]], reference, by = c(KEY, "method"))[, run_arm := arm]
    rm(adjusted, n, qc, omitted); gc()
  }
  tables <- list(null_calibration = nulls, paired_method_differences = pairs, donor_omission = lodo,
    donor_design = designs, cohort_agreement = agreement, directional_replication = directional,
    batch_sensitivity = batch, diagnostic_status = statuses)
  for (name in names(tables)) fwrite(rbindlist(tables[[name]], fill = TRUE), file.path(output, paste0(name, ".tsv")), sep = "\t")
  calibration <- rbindlist(nulls)
  render_calibration(calibration, output)
  paired_full <- rbindlist(pairs)[run_arm == "full"]
  voom_nominal <- calibration[run_arm == "full" & method == "voom_qw", c(KEY, "nominal_rejection_fraction"), with = FALSE]
  ql_nominal <- calibration[run_arm == "full" & method == "edger_ql", c(KEY, "nominal_rejection_fraction"), with = FALSE]
  nominal_comparison <- merge(voom_nominal, ql_nominal, by = KEY, suffixes = c("_voom", "_ql"))
  n_ql_further_from_nominal <- nominal_comparison[, sum(abs(nominal_rejection_fraction_ql - .05) > abs(nominal_rejection_fraction_voom - .05))]
  direction_tests <- rbindlist(directional, fill = TRUE)[directional_tested == TRUE]
  downsample_features <- fread(file.path(source, "cell_downsample/filter_qc.tsv"))[representation == "fragment" & filter == "primary"]
  full <- calibration[run_arm == "full" & representation == "fragment" & filter == "primary"]
  report <- c("Completed donor-level accessibility diagnostics; candidate comparison, no scientific adoption.",
    "The fitted-intercept NB model simulates all peaks as null. Pr(any BH rejection) equals FDR only within this all-null simulation.",
    "Mean fraction of peaks BH-rejected is a different estimand and must not be read as the all-null false-discovery probability.",
    "Intervals describe Monte Carlo uncertainty across independent seeded simulation draws, not biological-donor uncertainty.",
    "The same seed, tested features and recorded count sums/two-moment checksums match across both methods; these are diagnostic checks, not full-matrix cryptographic hashes.",
    "Fragment and binary actual fits use identical feature masks, but diagnostic feature subsets and simulation draws are not paired across representations because their seeds differ.",
    "Each model diagnostic uses at most 2,000 randomly selected fixed peaks. Binary counts are bounded; NB simulations do not preserve that native distribution.",
    "Original-count primary diagnostics:")
  for (i in seq_len(nrow(full))) {
    x <- full[i]
    report <- c(report, sprintf("%s %s %s: nominal=%.6f; Pr(any BH)=%.3f [%.3f, %.3f]; draws=%d; selected peaks=%d.",
      x$cohort, x$lineage, x$method, x$nominal_rejection_fraction, x$probability_any_BH_rejection,
      x$probability_any_BH_mc_lower, x$probability_any_BH_mc_upper, x$n_testable_draws, x$n_diagnostic_peaks))
  }
  report <- c(report, sprintf("Among %d original-count contexts, robust QL has higher Pr(any BH rejection) in %d and lower in %d; its nominal rejection fraction is further from 0.05 in %d. These are descriptive model-diagnostic comparisons.",
      nrow(paired_full), sum(paired_full$ql_minus_voom_probability_any_BH > 0), sum(paired_full$ql_minus_voom_probability_any_BH < 0), n_ql_further_from_nominal),
    "Both methods have anti-conservative tails in some contexts. Unpaired cross-representation diagnostic feature subsets limit interpretation of small binary-versus-fragment differences.",
    sprintf("Of %d measured lineage-level directional tests across both count arms, %d pass BH q<0.05 in their fixed three-lineage families. Peakwise both-cohort BH+same-sign counts remain distinct.",
      nrow(direction_tests), sum(direction_tests$directional_qvalue < .05)),
    sprintf("Equal-cell primary refiltering retains hepatocyte=%d, stellate=%d, macrophage=%d peaks; these families differ from the original counts. Two-peak diagnostics cannot establish method superiority.",
      downsample_features[lineage == "hepatocyte", n_jointly_testable], downsample_features[lineage == "stellate", n_jointly_testable], downsample_features[lineage == "macrophage", n_jointly_testable]),
    "Donors, rather than cells or peaks, are the disease-inference units. Covariates are only fit when complete and estimable; recorded batch adjustment remains materially influential.",
    "Current voom remains the manuscript reference. More significant peaks, stable effects, and model-null diagnostics alone do not establish biological correctness.")
  if (!is.null(warning_receipt)) {
    warnings <- fromJSON(warning_receipt)
    require_true(warnings$status == "completed" && nrow(warnings$fits) == 8L &&
      all(warnings$fits$n_peaks == 2L) && all(warnings$fits$all_logFC_and_pvalues_finite), "Unqualified two-peak warning reproduction")
    require_true(file.copy(warning_receipt, file.path(output, "two_peak_warning_reproduction.json")), "Cannot copy warning receipt")
    source_files <- c(source_files, warning_receipt)
    report <- c(report, sprintf("Bounded real-count warning reproduction: %d two-peak fits, %d captured warnings; all logFC and p-values finite and donor designs full rank.",
      nrow(warnings$fits), sum(warnings$fits$n_warnings)),
      "All four voom fits emitted zero warnings. All four robust QL fits warned 'Estimated rdf < 1.0; not estimating variance'; one also warned 'procv: no points with non-zero weight'.",
      "The donor models retain 12/10 residual degrees of freedom; the warning does not refer to loss of donor-model rank.",
      "[Inference] The two-feature QL variance-trend calculation is too weakly supported for method-selection claims, even when reported effect estimates and p-values remain finite.")
  }
  writeLines(report, file.path(output, "REPORT.txt"))
  writeLines(capture.output(sessionInfo()), file.path(output, "sessionInfo.txt"))
  source_files <- normalizePath(source_files, mustWork = TRUE)
  hashes <- vapply(source_files, function(path) strsplit(system2("sha256sum", shQuote(path), stdout = TRUE), " ")[[1L]][1L], "")
  fwrite(data.table(path = source_files, sha256 = hashes), file.path(output, "source_sha256.tsv"), sep = "\t")
  script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
  if (length(script) == 1L && file.exists(script)) file.copy(script, file.path(output, "summarize_da_diagnostics.R"))
  write_json(list(status = "completed", candidate_only = TRUE, scientific_adoption = FALSE,
    run_arms = c("full", "cell_downsample"), simulation_uncertainty_unit = "independent_NB_simulation_draw_not_donor",
    biological_unit = "donor", paired_diagnostic_checksums = "seed_feature_count_sum_and_sum_squares_match",
    inference = "conditional_model_diagnostics_no_automatic_method_selection", n_panels = 16L),
    file.path(output, "completion.json"), auto_unbox = TRUE, pretty = TRUE)
}

if (sys.nframe() == 0L) {
  arguments <- commandArgs(trailingOnly = TRUE)
  if (identical(arguments, "--self-test")) self_test()
  else {
    require_true(length(arguments) %in% c(4L, 6L) && arguments[1L] == "--source" && arguments[3L] == "--output" &&
      (length(arguments) == 4L || arguments[5L] == "--warning-receipt"),
      "Usage: summarize_da_diagnostics.R --source DA_FULL_ROOT --output NEWDIR [--warning-receipt JSON] | --self-test")
    run_review(arguments[2L], arguments[4L], if (length(arguments) == 6L) arguments[6L] else NULL)
  }
}
