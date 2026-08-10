#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("usage: 06_analyze_gse175448_human.R <candidate_root>")
candidate <- normalizePath(args[[1L]], mustWork = TRUE)
file_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
stopifnot(length(file_argument) == 1L)
script_dir <- dirname(normalizePath(sub("^--file=", "", file_argument[[1L]])))
source(file.path(script_dir, "lib_plan44_scores.R"))

input_path <- file.path(candidate, "source_expression/GSE175448_Centaur_rawReadCount.txt.gz")
soft_path <- file.path(candidate, "source_metadata/GSE175448_family.soft.gz")
output_dir <- file.path(candidate, "human_response")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(input_path), file.exists(soft_path))

parse_soft_metadata <- function(path) {
  lines <- readLines(gzfile(path), warn = FALSE)
  starts <- grep("^\\^SAMPLE = ", lines)
  ends <- c(starts[-1L] - 1L, length(lines))
  rows <- vector("list", length(starts))
  extract_one <- function(block, prefix) {
    values <- sub(prefix, "", grep(prefix, block, value = TRUE))
    stopifnot(length(values) == 1L)
    values[[1L]]
  }
  for (index in seq_along(starts)) {
    block <- lines[starts[[index]]:ends[[index]]]
    title <- extract_one(block, "^!Sample_title = ")
    description <- extract_one(block, "^!Sample_description = ")
    timing <- extract_one(block, "^!Sample_characteristics_ch1 = timing of biopsy: ")
    improvement <- extract_one(block, "^!Sample_characteristics_ch1 = fibrosis improvement: ")
    description_parts <- strsplit(description, "_", fixed = TRUE)[[1L]]
    source_arm <- description_parts[[2L]]
    rows[[index]] <- data.table(
      sample_title = title,
      participant = description_parts[[1L]],
      source_arm = source_arm,
      arm = if (source_arm == "drugA") "cenicriviroc" else "placebo",
      timing = timing,
      improvement = improvement
    )
  }
  rbindlist(rows)
}

metadata <- parse_soft_metadata(soft_path)
stopifnot(nrow(metadata) == 38L, uniqueN(metadata$participant) == 19L)
participant_response <- metadata[timing == "post-treatment", .(
  participant,
  arm,
  improved = improvement == "improved"
)]
stopifnot(
  participant_response[arm == "cenicriviroc", sum(improved)] == 4L,
  participant_response[arm == "placebo", sum(improved)] == 3L
)

counts_dt <- fread(input_path)
sample_columns <- metadata$sample_title
stopifnot(all(sample_columns %in% names(counts_dt)))
count_matrix <- as.matrix(counts_dt[, ..sample_columns])
storage.mode(count_matrix) <- "numeric"
symbol_counts <- aggregate_counts_by_symbol(count_matrix, counts_dt$gene_name)
colnames(symbol_counts) <- sample_columns
frozen <- read_frozen_plan44(candidate)

fit_response <- function(participant_scores) {
  participant_scores[, arm := factor(arm, levels = c("placebo", "cenicriviroc"))]
  fit <- lm(delta_topology ~ improved + arm, data = participant_scores)
  coefficient <- summary(fit)$coefficients["improvedTRUE", ]
  list(
    estimate = unname(coefficient[[1L]]),
    se = unname(coefficient[[2L]]),
    statistic = unname(coefficient[[3L]]),
    parametric_p = unname(coefficient[[4L]])
  )
}

enumerate_stratified_exact <- function(participant_scores, observed_estimate) {
  participant_scores[, arm := factor(arm, levels = c("placebo", "cenicriviroc"))]
  cvc_indices <- which(participant_scores$arm == "cenicriviroc")
  placebo_indices <- which(participant_scores$arm == "placebo")
  cvc_combinations <- combn(cvc_indices, 4L)
  placebo_combinations <- combn(placebo_indices, 3L)
  null_effects <- numeric(ncol(cvc_combinations) * ncol(placebo_combinations))
  cursor <- 1L
  for (cvc_column in seq_len(ncol(cvc_combinations))) {
    for (placebo_column in seq_len(ncol(placebo_combinations))) {
      labels <- rep(FALSE, nrow(participant_scores))
      labels[cvc_combinations[, cvc_column]] <- TRUE
      labels[placebo_combinations[, placebo_column]] <- TRUE
      design <- model.matrix(~ labels + arm, data = participant_scores)
      null_effects[[cursor]] <- lm.fit(design, participant_scores$delta_topology)$coefficients[["labelsTRUE"]]
      cursor <- cursor + 1L
    }
  }
  stopifnot(length(null_effects) == 15120L, all(is.finite(null_effects)))
  list(
    p_two_sided = mean(abs(null_effects) >= abs(observed_estimate) - 1e-12),
    null_median = median(null_effects),
    null_mean = mean(null_effects),
    n_allocations = length(null_effects)
  )
}

score_one_expression <- function(expression, score_mode, normalization_name) {
  z <- standardize_global_reference(expression, metadata$timing == "pre-treatment")
  scored <- score_frozen_geometry(
    z,
    frozen$membership,
    frozen$geometry,
    mode = score_mode
  )
  required_lineages <- c("hepatocytes", "fibroblasts", "macrophages")
  required_status <- scored$lineage_testability[cell_type %in% required_lineages]
  if (nrow(required_status) != 3L || !all(required_status$testable)) {
    return(list(scored = scored, sample_scores = NULL, participant_scores = NULL, effect = NULL))
  }
  lineage_long <- matrix_to_long(
    scored$lineage_scores,
    "cell_type",
    "sample_title",
    "lineage_score"
  )
  score_wide <- dcast(
    lineage_long[cell_type %in% required_lineages],
    sample_title ~ cell_type,
    value.var = "lineage_score"
  )
  sample_scores <- merge(metadata, score_wide, by = "sample_title", all.x = TRUE, sort = FALSE)
  stopifnot(nrow(sample_scores) == 38L)
  sample_scores[, remodeling := (fibroblasts + macrophages) / 2]
  sample_scores[, topology := remodeling - hepatocytes]
  sample_scores[, `:=`(mode = score_mode, normalization = normalization_name)]

  participant_scores <- dcast(
    sample_scores,
    participant + arm + mode + normalization ~ timing,
    value.var = c("topology", "hepatocytes", "fibroblasts", "macrophages", "remodeling")
  )
  participant_scores <- merge(
    participant_scores,
    participant_response,
    by = c("participant", "arm"),
    all.x = TRUE
  )
  participant_scores[, delta_topology := `topology_post-treatment` - `topology_pre-treatment`]
  participant_scores[, delta_fibroblast_vs_hepatocyte :=
    (`fibroblasts_post-treatment` - `fibroblasts_pre-treatment`) -
      (`hepatocytes_post-treatment` - `hepatocytes_pre-treatment`)]
  participant_scores[, delta_macrophage_vs_hepatocyte :=
    (`macrophages_post-treatment` - `macrophages_pre-treatment`) -
      (`hepatocytes_post-treatment` - `hepatocytes_pre-treatment`)]
  stopifnot(nrow(participant_scores) == 19L, all(!is.na(participant_scores$improved)))
  effect <- fit_response(copy(participant_scores))
  exact <- enumerate_stratified_exact(copy(participant_scores), effect$estimate)
  effect_row <- data.table(
    mode = score_mode,
    normalization = normalization_name,
    contrast = "fibrosis_improved_vs_not_delta_remodeling_minus_hepatocyte",
    estimate = effect$estimate,
    se = effect$se,
    statistic = effect$statistic,
    parametric_p = effect$parametric_p,
    exact_p_two_sided = exact$p_two_sided,
    exact_null_median = exact$null_median,
    exact_null_mean = exact$null_mean,
    n_exact_allocations = exact$n_allocations,
    expected_direction = "negative",
    expected_direction_met = effect$estimate < 0
  )
  list(
    scored = scored,
    sample_scores = sample_scores,
    participant_scores = participant_scores,
    effect = effect_row
  )
}

tmm_expression <- tmm_logcpm(symbol_counts)
modes <- c("primary", "equal_gene", "leave_top_gene", "sign_loading", "leave_top_program")
results <- lapply(modes, function(mode) score_one_expression(tmm_expression, mode, "TMM_logCPM"))
names(results) <- modes
deseq_result <- score_one_expression(
  deseq_sizefactor_log(symbol_counts),
  "primary",
  "DESeq2_sizefactor_log2"
)
results[["deseq_sizefactor"]] <- deseq_result

program_testability <- rbindlist(lapply(names(results), function(name) {
  output <- copy(results[[name]]$scored$program_testability)
  output[, result_id := name]
  output
}))
lineage_testability <- rbindlist(lapply(names(results), function(name) {
  output <- copy(results[[name]]$scored$lineage_testability)
  output[, result_id := name]
  output
}))
sample_scores <- rbindlist(lapply(results, `[[`, "sample_scores"), use.names = TRUE, fill = TRUE)
participant_scores <- rbindlist(lapply(results, `[[`, "participant_scores"), use.names = TRUE, fill = TRUE)
effects <- rbindlist(lapply(results, `[[`, "effect"), use.names = TRUE, fill = TRUE)

primary_participants <- participant_scores[mode == "primary" & normalization == "TMM_logCPM"]
primary_effect <- effects[mode == "primary" & normalization == "TMM_logCPM"]
stopifnot(nrow(primary_effect) == 1L)

leave_one_out <- rbindlist(lapply(primary_participants$participant, function(omitted) {
  subset_scores <- copy(primary_participants[participant != omitted])
  effect <- fit_response(subset_scores)
  data.table(
    omitted_participant = omitted,
    estimate = effect$estimate,
    se = effect$se,
    statistic = effect$statistic,
    parametric_p = effect$parametric_p,
    expected_direction_met = effect$estimate < 0
  )
}))

baseline_adjusted_fit <- lm(
  delta_topology ~ improved + arm + `topology_pre-treatment`,
  data = primary_participants
)
baseline_adjusted_coefficient <- summary(baseline_adjusted_fit)$coefficients["improvedTRUE", ]
lineage_specific <- rbindlist(lapply(
  c("delta_fibroblast_vs_hepatocyte", "delta_macrophage_vs_hepatocyte"),
  function(outcome) {
    model <- lm(reformulate(c("improved", "arm"), response = outcome), data = primary_participants)
    coefficient <- summary(model)$coefficients["improvedTRUE", ]
    data.table(
      sensitivity = outcome,
      estimate = coefficient[[1L]],
      se = coefficient[[2L]],
      statistic = coefficient[[3L]],
      parametric_p = coefficient[[4L]],
      expected_direction_met = coefficient[[1L]] < 0
    )
  }
))
sensitivity <- rbindlist(list(
  effects[, .(
    sensitivity = paste(mode, normalization, sep = "__"),
    estimate,
    se,
    statistic,
    parametric_p,
    expected_direction_met
  )],
  data.table(
    sensitivity = "baseline_axis_adjusted",
    estimate = baseline_adjusted_coefficient[[1L]],
    se = baseline_adjusted_coefficient[[2L]],
    statistic = baseline_adjusted_coefficient[[3L]],
    parametric_p = baseline_adjusted_coefficient[[4L]],
    expected_direction_met = baseline_adjusted_coefficient[[1L]] < 0
  ),
  lineage_specific
), use.names = TRUE, fill = TRUE)

required_mode_direction <- all(
  effects[!(mode == "primary" & normalization == "TMM_logCPM"), expected_direction_met]
)
human_gate_pass <- primary_effect$estimate < 0 && primary_effect$exact_p_two_sided < 0.05 &&
  all(leave_one_out$expected_direction_met) && required_mode_direction

fwrite(metadata, file.path(output_dir, "human_pair_manifest.tsv"), sep = "\t")
fwrite(program_testability, file.path(output_dir, "human_program_testability.tsv"), sep = "\t")
fwrite(lineage_testability, file.path(output_dir, "human_lineage_testability.tsv"), sep = "\t")
fwrite(sample_scores, file.path(output_dir, "human_lineage_scores.tsv"), sep = "\t")
fwrite(participant_scores, file.path(output_dir, "human_participant_changes.tsv"), sep = "\t")
fwrite(effects, file.path(output_dir, "human_response_effects.tsv"), sep = "\t")
fwrite(leave_one_out, file.path(output_dir, "human_leave_one_participant_out.tsv"), sep = "\t")
fwrite(sensitivity, file.path(output_dir, "human_sensitivity.tsv"), sep = "\t")
fwrite(data.table(
  gate = "paired_human_remodeling_disassembly",
  status = if (human_gate_pass) "pass" else "complete_nonconfirmatory",
  primary_estimate = primary_effect$estimate,
  primary_exact_p = primary_effect$exact_p_two_sided,
  primary_direction_met = primary_effect$estimate < 0,
  all_leave_one_participant_direction = all(leave_one_out$expected_direction_met),
  all_required_mode_direction = required_mode_direction,
  n_participants = nrow(primary_participants),
  n_improved = sum(primary_participants$improved),
  n_cenicriviroc = sum(primary_participants$arm == "cenicriviroc"),
  n_placebo = sum(primary_participants$arm == "placebo")
), file.path(output_dir, "human_gate_status.tsv"), sep = "\t")

cat(
  "PLAN44_HUMAN_COMPLETE",
  paste0("gate=", if (human_gate_pass) "pass" else "complete_nonconfirmatory"),
  paste0("estimate=", signif(primary_effect$estimate, 6)),
  paste0("exact_p=", signif(primary_effect$exact_p_two_sided, 6)),
  sep = "\t"
)
cat("\n")
