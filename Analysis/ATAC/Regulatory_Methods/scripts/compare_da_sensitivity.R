#!/usr/bin/env Rscript
# Descriptive full versus whole-cell-downsample comparison on shared peaks.
suppressPackageStartupMessages({ library(data.table); library(jsonlite) })

sensitivity_require <- function(value, message) if (!isTRUE(value)) stop(message, call. = FALSE)
SENSITIVITY_KEYS <- c("cohort", "lineage", "representation", "filter", "method")

compare_effect_arm <- function(full, downsample, gene_map = NULL) {
  sensitivity_require(!anyDuplicated(full$peak_coordinate) && !anyDuplicated(downsample$peak_coordinate), "Duplicated peak within comparison arm")
  paired <- merge(full[, .(peak_coordinate, full_logFC = logFC, full_qvalue = qvalue)],
    downsample[, .(peak_coordinate, downsample_logFC = logFC, downsample_qvalue = qvalue)], by = "peak_coordinate")
  valid <- is.finite(paired$full_logFC) & is.finite(paired$downsample_logFC)
  effect <- paired[valid]
  delta <- effect$downsample_logFC - effect$full_logFC
  rho <- if (nrow(effect) > 1L && sd(effect$full_logFC) > 0 && sd(effect$downsample_logFC) > 0)
    cor(effect$full_logFC, effect$downsample_logFC, method = "spearman") else NA_real_
  result <- data.table(n_peaks_full = nrow(full), n_peaks_downsample = nrow(downsample),
    n_common_peaks = nrow(paired), n_full_only_peaks = nrow(full) - nrow(paired),
    n_downsample_only_peaks = nrow(downsample) - nrow(paired),
    fraction_full_peaks_retained = if (nrow(full)) nrow(paired) / nrow(full) else NA_real_,
    bh_universe_changed = !setequal(full$peak_coordinate, downsample$peak_coordinate),
    n_finite_common_effects = nrow(effect), effect_spearman = rho,
    effect_sign_agreement = if (nrow(effect)) mean(sign(effect$full_logFC) == sign(effect$downsample_logFC)) else NA_real_,
    median_signed_log2_effect_change = if (length(delta)) median(delta) else NA_real_,
    median_absolute_log2_effect_change = if (length(delta)) median(abs(delta)) else NA_real_,
    p95_absolute_log2_effect_change = if (length(delta)) as.numeric(quantile(abs(delta), .95)) else NA_real_,
    n_bh_q05_full_original_family = sum(full$qvalue < .05, na.rm = TRUE),
    n_bh_q05_downsample_original_family = sum(downsample$qvalue < .05, na.rm = TRUE),
    n_bh_q05_full_on_common_peaks = sum(paired$full_qvalue < .05, na.rm = TRUE),
    n_bh_q05_downsample_on_common_peaks = sum(paired$downsample_qvalue < .05, na.rm = TRUE),
    qvalue_interpretation = "q_values_retain_their_original_tested_peak_families; not_recomputed_on_common_peaks",
    interpretation = "descriptive_sensitivity; significance_count_changes_are_not_method_improvement")
  if (is.null(gene_map)) {
    result[, `:=`(gene_coverage_status = "unavailable_without_explicit_peak_gene_map",
      n_mapped_genes_full = NA_integer_, n_mapped_genes_downsample = NA_integer_,
      n_mapped_genes_common_peaks = NA_integer_, n_genes_shared_across_arms = NA_integer_)]
  } else {
    full_genes <- unique(gene_map[peak_coordinate %chin% full$peak_coordinate, gene_id])
    down_genes <- unique(gene_map[peak_coordinate %chin% downsample$peak_coordinate, gene_id])
    common_genes <- unique(gene_map[peak_coordinate %chin% paired$peak_coordinate, gene_id])
    result[, `:=`(gene_coverage_status = "coverage_under_supplied_mapping; no_new_target_assignment",
      n_mapped_genes_full = length(full_genes), n_mapped_genes_downsample = length(down_genes),
      n_mapped_genes_common_peaks = length(common_genes), n_genes_shared_across_arms = length(intersect(full_genes, down_genes)))]
  }
  result
}

read_sensitivity_inputs <- function(root) {
  needed <- c("peak_coordinate", "logFC", "qvalue", "effect_units", SENSITIVITY_KEYS)
  peaks <- fread(file.path(root, "peak_results.tsv.gz"), select = needed)
  sensitivity_require(!anyDuplicated(peaks[, c(SENSITIVITY_KEYS, "peak_coordinate"), with = FALSE]), "Duplicated comparison rows")
  setindexv(peaks, SENSITIVITY_KEYS)
  qc_path <- file.path(root, "model_qc.tsv.gz")
  qc <- if (file.exists(qc_path)) fread(qc_path, select = c("donor_id", "condition", SENSITIVITY_KEYS)) else NULL
  summary <- fread(file.path(root, "cohort_summary.tsv"))
  summary <- summary[cohort != "cross_cohort"]
  list(peaks = peaks, qc = qc, summary = summary)
}

compare_rosters <- function(full, downsample) {
  if (is.null(full) || is.null(downsample)) return(data.table(donor_roster_status = "unavailable_without_both_model_qc_tables"))
  f <- unique(full[, .(donor_id, condition)])
  d <- unique(downsample[, .(donor_id, condition)])
  setorder(f, donor_id); setorder(d, donor_id)
  data.table(n_donors_full = nrow(f), n_donors_downsample = nrow(d),
    n_shared_donors = length(intersect(f$donor_id, d$donor_id)),
    n_normal_full = sum(f$condition == "NORMAL"), n_mash_full = sum(f$condition == "MASH"),
    n_normal_downsample = sum(d$condition == "NORMAL"), n_mash_downsample = sum(d$condition == "MASH"),
    donor_roster_status = if (!nrow(f) || !nrow(d)) "unavailable_no_fitted_donor_roster" else if (identical(as.data.frame(f), as.data.frame(d))) "identical_donor_conditions" else "changed_donor_roster_or_conditions")
}

run_sensitivity_comparison <- function(args) {
  sensitivity_require(nzchar(Sys.getenv("SLURM_JOB_ID")), "Requires a Slurm compute allocation")
  full <- read_sensitivity_inputs(args$full)
  down <- read_sensitivity_inputs(args$downsample)
  sensitivity_require(!file.exists(args$output), "Refusing to overwrite output")
  gene_map <- NULL
  if (!is.null(args$peak_gene_map)) {
    gene_map <- fread(args$peak_gene_map, select = c("peak_coordinate", "gene_id"))
    sensitivity_require(!anyNA(gene_map) && all(grepl("^ENSG[0-9]+(\\.[0-9]+)?$", gene_map$gene_id)),
      "Gene coverage requires explicitly supplied Ensembl gene identifiers")
    gene_map[, gene_id := sub("\\.[0-9]+$", "", gene_id)]
    gene_map <- unique(gene_map)
  }
  arms <- unique(rbindlist(list(full$summary[, ..SENSITIVITY_KEYS], down$summary[, ..SENSITIVITY_KEYS])))
  rows <- lapply(seq_len(nrow(arms)), function(index) {
    arm <- arms[index]
    f <- full$peaks[arm, on = SENSITIVITY_KEYS, nomatch = 0L]
    d <- down$peaks[arm, on = SENSITIVITY_KEYS, nomatch = 0L]
    units <- unique(c(f$effect_units, d$effect_units))
    sensitivity_require(length(units) <= 1L, "Effect units differ within an arm")
    stats <- compare_effect_arm(f, d, gene_map)
    fq <- if (is.null(full$qc)) NULL else full$qc[arm, on = SENSITIVITY_KEYS, nomatch = 0L]
    dq <- if (is.null(down$qc)) NULL else down$qc[arm, on = SENSITIVITY_KEYS, nomatch = 0L]
    roster <- compare_rosters(fq, dq)
    f_status <- full$summary[arm, on = SENSITIVITY_KEYS, nomatch = 0L]$status
    d_status <- down$summary[arm, on = SENSITIVITY_KEYS, nomatch = 0L]$status
    cbind(arm, data.table(full_status = paste(unique(f_status), collapse = ";"),
      downsample_status = paste(unique(d_status), collapse = ";"),
      effect_units = if (length(units)) units else NA_character_), stats, roster)
  })
  dir.create(args$output, recursive = TRUE)
  fwrite(rbindlist(rows, fill = TRUE), file.path(args$output, "sensitivity_comparison.tsv"), sep = "\t")
  write_json(list(status = "descriptive_comparison_complete", full_root = normalizePath(args$full),
    downsample_root = normalizePath(args$downsample), gene_map = args$peak_gene_map,
    gene_map_label = args$gene_map_label, gene_id_policy = "version_suffix_removed_from_explicit_Ensembl_IDs",
    contrast = "MASH_minus_normal", common_features = "peak_coordinate_intersection_within_each_prespecified_arm",
    bh_universes = "retained_original_families; universe_changes_reported_explicitly",
    biological_unit = "donor", inference = "descriptive; no_method_winner_or_new_p_values"),
    file.path(args$output, "run_summary.json"), auto_unbox = TRUE, pretty = TRUE, null = "null")
  writeLines(capture.output(sessionInfo()), file.path(args$output, "sessionInfo.txt"))
}

parse_sensitivity_args <- function(args) {
  values <- list(peak_gene_map = NULL, gene_map_label = NULL)
  sensitivity_require(length(args) >= 6L && length(args) %% 2L == 0L, "Arguments require --full, --downsample and --output values")
  for (i in seq.int(1L, length(args), by = 2L)) {
    sensitivity_require(args[i] %in% c("--full", "--downsample", "--output", "--peak-gene-map", "--gene-map-label"), "Unknown argument")
    values[[gsub("-", "_", substring(args[i], 3L))]] <- args[i + 1L]
  }
  sensitivity_require(all(vapply(values[c("full", "downsample", "output")], function(x) !is.null(x) && nzchar(x), logical(1))),
    "Usage: compare_da_sensitivity.R --full DIR --downsample DIR --output NEWDIR [--peak-gene-map TSV --gene-map-label LABEL]")
  if (!is.null(values$peak_gene_map)) sensitivity_require(!is.null(values$gene_map_label) && nzchar(values$gene_map_label),
    "Supply --gene-map-label to identify the mapping evidence type")
  values
}

if (sys.nframe() == 0L) run_sensitivity_comparison(parse_sensitivity_args(commandArgs(trailingOnly = TRUE)))
