#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("usage: 05_analyze_gse168285_mps.R <candidate_root>")
candidate <- normalizePath(args[[1L]], mustWork = TRUE)
file_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
stopifnot(length(file_argument) == 1L)
script_dir <- dirname(normalizePath(sub("^--file=", "", file_argument[[1L]])))
source(file.path(script_dir, "lib_plan44_scores.R"))

input_dir <- file.path(candidate, "source_expression")
output_dir <- file.path(candidate, "mps")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

required <- c(
  "GSE168285_raw_counts.txt.gz",
  "GSE168285_gene_annotation.txt.gz",
  "GSE168285_sample_meta_data.txt.gz"
)
stopifnot(all(file.exists(file.path(input_dir, required))))

frozen <- read_frozen_plan44(candidate)
counts_dt <- fread(file.path(input_dir, "GSE168285_raw_counts.txt.gz"))
annotation <- fread(file.path(input_dir, "GSE168285_gene_annotation.txt.gz"))
metadata <- fread(file.path(input_dir, "GSE168285_sample_meta_data.txt.gz"))
stopifnot(nrow(metadata) == 179L, uniqueN(metadata$sampleName) == 179L)
sample_columns <- metadata$sampleName
stopifnot(all(sample_columns %in% names(counts_dt)))

count_matrix <- as.matrix(counts_dt[, ..sample_columns])
storage.mode(count_matrix) <- "numeric"
symbols <- annotation$GeneSymbol[match(counts_dt$rowname, annotation$EnsemblAcc)]
symbol_counts <- aggregate_counts_by_symbol(count_matrix, symbols)
colnames(symbol_counts) <- sample_columns
expression <- tmm_logcpm(symbol_counts)

metadata[, experiment := sub("_.*$", "", sampleName)]
metadata[, high_NPC := as.integer(NPC == "high")]
metadata[, fat := as.integer(grepl("(^|_)Fat($|_)", treatment))]
metadata[, fructose := as.integer(grepl("Fructose", treatment))]
metadata[, cholesterol := as.integer(grepl("Cholesterol", treatment))]
metadata[, LPS := as.integer(grepl("LPS", treatment))]
metadata[, TGF_beta := as.integer(grepl("TGF-B", treatment))]
metadata[, reference_control := treatment == "Reduced_NPC_Lean"]
stopifnot(
  uniqueN(metadata$experiment) == 4L,
  all(metadata[, sum(reference_control), by = experiment]$V1 >= 2L)
)

z <- standardize_by_reference(
  expression,
  strata = metadata$experiment,
  reference = metadata$reference_control
)

modes <- c(
  "primary",
  "equal_gene",
  "leave_top_gene",
  "sign_loading",
  "leave_top_program",
  "shared_gene_removed"
)
all_program_testability <- list()
all_lineage_testability <- list()
all_lineage_scores <- list()
all_condition_scores <- list()
all_effects <- list()
all_lodo <- list()

fit_condition_model <- function(condition_scores) {
  fit <- lm(
    topology ~ experiment + high_NPC + fat + fructose + cholesterol + LPS + TGF_beta,
    data = condition_scores
  )
  coefficients <- summary(fit)$coefficients
  stopifnot(all(c("high_NPC", "fat", "TGF_beta") %in% rownames(coefficients)))
  data.table(
    contrast = c("fat_hepatocyte_specificity", "tgfb_remodeling_specificity", "npc_remodeling_specificity"),
    raw_term = c("fat", "TGF_beta", "high_NPC"),
    orientation = c(1, -1, -1),
    estimate = c(coefficients["fat", 1L], coefficients["TGF_beta", 1L], coefficients["high_NPC", 1L]),
    se = c(coefficients["fat", 2L], coefficients["TGF_beta", 2L], coefficients["high_NPC", 2L]),
    statistic = c(coefficients["fat", 3L], coefficients["TGF_beta", 3L], coefficients["high_NPC", 3L]),
    parametric_p = c(coefficients["fat", 4L], coefficients["TGF_beta", 4L], coefficients["high_NPC", 4L])
  )[, `:=`(
    oriented_estimate = orientation * estimate,
    oriented_statistic = orientation * statistic,
    expected_direction_met = orientation * estimate > 0
  )]
}

for (mode in modes) {
  scored <- score_frozen_geometry(z, frozen$membership, frozen$geometry, mode = mode)
  all_program_testability[[mode]] <- scored$program_testability
  all_lineage_testability[[mode]] <- scored$lineage_testability
  lineage_long <- matrix_to_long(scored$lineage_scores, "cell_type", "sampleName", "lineage_score")
  lineage_long[, mode := mode]
  all_lineage_scores[[mode]] <- lineage_long

  primary_lineages <- c("hepatocytes", "fibroblasts", "macrophages")
  status <- scored$lineage_testability[cell_type %in% primary_lineages]
  if (nrow(status) != 3L || !all(status$testable)) next
  score_wide <- dcast(
    lineage_long[cell_type %in% primary_lineages],
    sampleName + mode ~ cell_type,
    value.var = "lineage_score"
  )
  sample_scores <- merge(
    metadata[, .(
      sampleName,
      experiment,
      treatment,
      batch,
      high_NPC,
      fat,
      fructose,
      cholesterol,
      LPS,
      TGF_beta
    )],
    score_wide,
    by = "sampleName",
    all.x = TRUE,
    sort = FALSE
  )
  stopifnot(nrow(sample_scores) == 179L, complete.cases(sample_scores[, .(hepatocytes, fibroblasts, macrophages)]))
  sample_scores[, remodeling := (fibroblasts + macrophages) / 2]
  sample_scores[, topology := hepatocytes - remodeling]
  condition_scores <- sample_scores[, lapply(.SD, mean), by = .(
    experiment,
    treatment,
    mode,
    high_NPC,
    fat,
    fructose,
    cholesterol,
    LPS,
    TGF_beta
  ), .SDcols = c("hepatocytes", "fibroblasts", "macrophages", "remodeling", "topology")]
  stopifnot(nrow(condition_scores) == 60L)
  all_condition_scores[[mode]] <- condition_scores
  effects <- fit_condition_model(condition_scores)
  effects[, mode := mode]
  all_effects[[mode]] <- effects

  for (omitted in sort(unique(condition_scores$experiment))) {
    subset_scores <- condition_scores[experiment != omitted]
    lodo_effect <- fit_condition_model(subset_scores)
    lodo_effect[, `:=`(mode = mode, omitted_experiment = omitted)]
    all_lodo[[paste(mode, omitted, sep = "_")]] <- lodo_effect
  }
}

program_testability <- rbindlist(all_program_testability, use.names = TRUE)
lineage_testability <- rbindlist(all_lineage_testability, use.names = TRUE)
lineage_scores <- rbindlist(all_lineage_scores, use.names = TRUE)
condition_scores_all <- rbindlist(all_condition_scores, use.names = TRUE)
effects_all <- rbindlist(all_effects, use.names = TRUE)
lodo_all <- rbindlist(all_lodo, use.names = TRUE)

# Outcome-blind Freedman-Lane residual permutation for the frozen primary mode.
primary_scores <- condition_scores_all[mode == "primary"]
full_design <- model.matrix(
  ~ experiment + high_NPC + fat + fructose + cholesterol + LPS + TGF_beta,
  data = primary_scores
)
null_design <- model.matrix(
  ~ experiment + fructose + cholesterol + LPS,
  data = primary_scores
)
stopifnot(qr(full_design)$rank == ncol(full_design))
response <- primary_scores$topology
observed_fit <- lm.fit(full_design, response)
inverse_information <- solve(crossprod(full_design))
residual_df <- nrow(full_design) - qr(full_design)$rank
observed_sigma <- sum(observed_fit$residuals^2) / residual_df
observed_se <- sqrt(diag(inverse_information) * observed_sigma)
observed_t <- observed_fit$coefficients / observed_se
observed_oriented <- c(
  fat = observed_t[["fat"]],
  TGF_beta = -observed_t[["TGF_beta"]],
  high_NPC = -observed_t[["high_NPC"]]
)

null_fit <- lm.fit(null_design, response)
blocks <- split(seq_len(nrow(primary_scores)), primary_scores$experiment)
set.seed(440168285)
n_permutations <- 100000L
null_statistics <- matrix(
  NA_real_,
  nrow = n_permutations,
  ncol = 3L,
  dimnames = list(NULL, names(observed_oriented))
)
for (iteration in seq_len(n_permutations)) {
  permutation <- seq_len(nrow(primary_scores))
  for (indices in blocks) permutation[indices] <- sample(indices, length(indices), replace = FALSE)
  permuted_response <- null_fit$fitted.values + null_fit$residuals[permutation]
  permuted_fit <- lm.fit(full_design, permuted_response)
  sigma_squared <- sum(permuted_fit$residuals^2) / residual_df
  standard_errors <- sqrt(diag(inverse_information) * sigma_squared)
  statistics <- permuted_fit$coefficients / standard_errors
  null_statistics[iteration, ] <- c(
    statistics[["fat"]],
    -statistics[["TGF_beta"]],
    -statistics[["high_NPC"]]
  )
}
max_null <- apply(null_statistics, 1L, max)
min_null <- apply(null_statistics, 1L, min)
permutation_rows <- data.table(
  contrast = c("fat_hepatocyte_specificity", "tgfb_remodeling_specificity", "npc_remodeling_specificity"),
  observed_oriented_statistic = as.numeric(observed_oriented),
  permutation_p_one_sided = vapply(seq_along(observed_oriented), function(index) {
    (1 + sum(null_statistics[, index] >= observed_oriented[[index]])) / (n_permutations + 1)
  }, numeric(1L)),
  maxT_adjusted_p_one_sided = vapply(observed_oriented, function(value) {
    (1 + sum(max_null >= value)) / (n_permutations + 1)
  }, numeric(1L)),
  null_median = apply(null_statistics, 2L, median),
  n_permutations = n_permutations
)
conjunction_p <- (1 + sum(min_null >= min(observed_oriented))) / (n_permutations + 1)

effects_all <- merge(
  effects_all,
  permutation_rows[, .(
    contrast,
    permutation_p_one_sided,
    maxT_adjusted_p_one_sided
  )],
  by = "contrast",
  all.x = TRUE
)
effects_all[mode != "primary", `:=`(
  permutation_p_one_sided = NA_real_,
  maxT_adjusted_p_one_sided = NA_real_
)]

primary_effects <- effects_all[mode == "primary"]
primary_lodo <- lodo_all[mode == "primary"]
primary_sensitivity <- effects_all[mode != "primary"]
required_sensitivity_modes <- c(
  "equal_gene",
  "leave_top_gene",
  "sign_loading",
  "leave_top_program",
  "shared_gene_removed"
)
all_primary_adjusted <- all(primary_effects$maxT_adjusted_p_one_sided < 0.05)
all_primary_direction <- all(primary_effects$expected_direction_met)
lodo_direction <- all(primary_lodo$expected_direction_met)
sensitivity_direction <- setequal(unique(primary_sensitivity$mode), required_sensitivity_modes) &&
  all(primary_sensitivity$expected_direction_met)
mps_gate_pass <- all_primary_adjusted && all_primary_direction && conjunction_p < 0.01 &&
  lodo_direction && sensitivity_direction

fwrite(metadata, file.path(output_dir, "mps_sample_manifest.tsv"), sep = "\t")
fwrite(program_testability, file.path(output_dir, "mps_program_testability.tsv"), sep = "\t")
fwrite(lineage_testability, file.path(output_dir, "mps_lineage_testability.tsv"), sep = "\t")
fwrite(lineage_scores, file.path(output_dir, "mps_lineage_scores.tsv"), sep = "\t")
fwrite(condition_scores_all, file.path(output_dir, "mps_condition_scores.tsv"), sep = "\t")
fwrite(effects_all, file.path(output_dir, "mps_lineage_effects.tsv"), sep = "\t")
fwrite(lodo_all, file.path(output_dir, "mps_leave_one_experiment_out.tsv"), sep = "\t")
fwrite(permutation_rows, file.path(output_dir, "mps_permutation_summary.tsv"), sep = "\t")
fwrite(data.table(
  gate = "controlled_multicellular_assembly",
  status = if (mps_gate_pass) "pass" else "complete_nonconfirmatory",
  all_primary_direction = all_primary_direction,
  all_primary_maxT_p_lt_0_05 = all_primary_adjusted,
  conjunction_p = conjunction_p,
  conjunction_p_lt_0_01 = conjunction_p < 0.01,
  all_lodo_direction = lodo_direction,
  all_sensitivity_direction = sensitivity_direction,
  n_condition_means = nrow(primary_scores),
  n_expression_libraries = nrow(metadata),
  n_human_donors = NA_character_
), file.path(output_dir, "mps_gate_status.tsv"), sep = "\t")

cat(
  "PLAN44_MPS_COMPLETE",
  paste0("gate=", if (mps_gate_pass) "pass" else "complete_nonconfirmatory"),
  paste0("conjunction_p=", signif(conjunction_p, 6)),
  paste0("testable_programs=", sum(program_testability[mode == "primary", testable])),
  sep = "\t"
)
cat("\n")
