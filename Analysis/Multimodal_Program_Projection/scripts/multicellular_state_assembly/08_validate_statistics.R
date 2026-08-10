#!/usr/bin/env Rscript

# Independent, table-level rederivation of the frozen Plan 44 effects.  This
# script intentionally does not source the producer's scoring or model code.

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("usage: 08_validate_statistics.R <candidate_root>")
candidate <- normalizePath(args[[1L]], mustWork = TRUE)

checks <- list()
add_check <- function(check_id, passed, detail) {
  checks[[length(checks) + 1L]] <<- data.table(
    check_id = check_id,
    passed = isTRUE(passed),
    detail = as.character(detail)
  )
}

mps <- fread(file.path(candidate, "mps/mps_condition_scores.tsv"))
mps <- mps[mode == "primary"]
add_check("mps_60_condition_units", nrow(mps) == 60L, nrow(mps))
add_check(
  "mps_topology_fixture",
  max(abs(mps$topology - (mps$hepatocytes - (mps$fibroblasts + mps$macrophages) / 2))) < 1e-12,
  "H-(F+M)/2"
)
mps_fit <- lm(
  topology ~ experiment + high_NPC + fat + fructose + cholesterol + LPS + TGF_beta,
  data = mps
)
mps_coef <- summary(mps_fit)$coefficients
effect_rows <- fread(file.path(candidate, "mps/mps_lineage_effects.tsv"))[mode == "primary"]
term_by_contrast <- c(
  fat_hepatocyte_specificity = "fat",
  npc_remodeling_specificity = "high_NPC",
  tgfb_remodeling_specificity = "TGF_beta"
)
for (contrast_id in names(term_by_contrast)) {
  term <- term_by_contrast[[contrast_id]]
  observed <- effect_rows[effect_rows$contrast == contrast_id]
  add_check(
    paste0("mps_ols_rederived:", contrast_id),
    nrow(observed) == 1L &&
      abs(observed$estimate - mps_coef[term, 1L]) < 1e-12 &&
      abs(observed$se - mps_coef[term, 2L]) < 1e-12 &&
      abs(observed$statistic - mps_coef[term, 3L]) < 1e-12,
    sprintf("beta=%.12g;t=%.12g", mps_coef[term, 1L], mps_coef[term, 3L])
  )
}

# Independently re-execute the sealed blocked Freedman-Lane null from the
# immutable condition table.  Matching the fixed-seed producer is a strict
# reproducibility test; the table-level OLS checks above are algorithmically
# independent of this replay.
full_design <- model.matrix(
  ~ experiment + high_NPC + fat + fructose + cholesterol + LPS + TGF_beta,
  data = mps
)
null_design <- model.matrix(~ experiment + fructose + cholesterol + LPS, data = mps)
full_fit <- lm.fit(full_design, mps$topology)
inverse_information <- solve(crossprod(full_design))
residual_df <- nrow(full_design) - qr(full_design)$rank
sigma_squared <- sum(full_fit$residuals^2) / residual_df
observed_t <- full_fit$coefficients / sqrt(diag(inverse_information) * sigma_squared)
observed_oriented <- c(
  fat = observed_t[["fat"]],
  TGF_beta = -observed_t[["TGF_beta"]],
  high_NPC = -observed_t[["high_NPC"]]
)
null_fit <- lm.fit(null_design, mps$topology)
blocks <- split(seq_len(nrow(mps)), mps$experiment)
set.seed(440168285)
n_permutations <- 100000L
null_statistics <- matrix(NA_real_, nrow = n_permutations, ncol = 3L)
for (iteration in seq_len(n_permutations)) {
  permutation <- seq_len(nrow(mps))
  for (indices in blocks) permutation[indices] <- sample(indices, length(indices))
  response <- null_fit$fitted.values + null_fit$residuals[permutation]
  fit <- lm.fit(full_design, response)
  sigma_squared <- sum(fit$residuals^2) / residual_df
  statistic <- fit$coefficients / sqrt(diag(inverse_information) * sigma_squared)
  null_statistics[iteration, ] <- c(
    statistic[["fat"]],
    -statistic[["TGF_beta"]],
    -statistic[["high_NPC"]]
  )
}
permutation <- fread(file.path(candidate, "mps/mps_permutation_summary.tsv"))
permutation_order <- c(
  fat_hepatocyte_specificity = 1L,
  tgfb_remodeling_specificity = 2L,
  npc_remodeling_specificity = 3L
)
max_null <- apply(null_statistics, 1L, max)
min_null <- apply(null_statistics, 1L, min)
for (contrast_id in names(permutation_order)) {
  index <- permutation_order[[contrast_id]]
  expected_p <- (1 + sum(null_statistics[, index] >= observed_oriented[[index]])) /
    (n_permutations + 1)
  expected_max_t <- (1 + sum(max_null >= observed_oriented[[index]])) /
    (n_permutations + 1)
  observed <- permutation[permutation$contrast == contrast_id]
  add_check(
    paste0("mps_permutation_replayed:", contrast_id),
    nrow(observed) == 1L &&
      abs(observed$permutation_p_one_sided - expected_p) < 1e-15 &&
      abs(observed$maxT_adjusted_p_one_sided - expected_max_t) < 1e-15,
    sprintf("p=%.8g;maxT=%.8g", expected_p, expected_max_t)
  )
}
expected_conjunction <- (1 + sum(min_null >= min(observed_oriented))) /
  (n_permutations + 1)
mps_gate <- fread(file.path(candidate, "mps/mps_gate_status.tsv"))
add_check(
  "mps_conjunction_replayed",
  abs(mps_gate$conjunction_p - expected_conjunction) < 1e-15,
  expected_conjunction
)
add_check(
  "mps_tgfb_opposite_in_every_experiment",
  all(fread(file.path(candidate, "mps/mps_leave_one_experiment_out.tsv"))[
    mode == "primary" & contrast == "tgfb_remodeling_specificity",
    expected_direction_met
  ] == FALSE),
  "4/4 leave-one-experiment-out estimates opposite"
)

human <- fread(file.path(candidate, "human_response/human_participant_changes.tsv"))[
  mode == "primary" & normalization == "TMM_logCPM"
]
add_check("human_19_participant_pairs", nrow(human) == 19L, nrow(human))
add_check(
  "human_delta_topology_fixture",
  max(abs(human$delta_topology -
    (human[["topology_post-treatment"]] - human[["topology_pre-treatment"]]))) < 1e-12,
  "post-pre"
)
human[, arm := factor(arm, levels = c("placebo", "cenicriviroc"))]
human_fit <- lm(delta_topology ~ improved + arm, data = human)
human_coef <- summary(human_fit)$coefficients["improvedTRUE", ]
human_effect <- fread(file.path(candidate, "human_response/human_response_effects.tsv"))[
  mode == "primary" & normalization == "TMM_logCPM"
]
add_check(
  "human_ols_rederived",
  nrow(human_effect) == 1L &&
    abs(human_effect$estimate - human_coef[[1L]]) < 1e-12 &&
    abs(human_effect$se - human_coef[[2L]]) < 1e-12 &&
    abs(human_effect$statistic - human_coef[[3L]]) < 1e-12,
  sprintf("beta=%.12g;t=%.12g", human_coef[[1L]], human_coef[[3L]])
)
cvc_indices <- which(human$arm == "cenicriviroc")
placebo_indices <- which(human$arm == "placebo")
cvc_combinations <- combn(cvc_indices, 4L)
placebo_combinations <- combn(placebo_indices, 3L)
null_effects <- numeric(ncol(cvc_combinations) * ncol(placebo_combinations))
cursor <- 1L
for (cvc_column in seq_len(ncol(cvc_combinations))) {
  for (placebo_column in seq_len(ncol(placebo_combinations))) {
    label <- rep(FALSE, nrow(human))
    label[cvc_combinations[, cvc_column]] <- TRUE
    label[placebo_combinations[, placebo_column]] <- TRUE
    null_effects[[cursor]] <- coef(lm(human$delta_topology ~ label + human$arm))[["labelTRUE"]]
    cursor <- cursor + 1L
  }
}
exact_p <- mean(abs(null_effects) >= abs(human_coef[[1L]]) - 1e-12)
add_check(
  "human_exact_15120_rederived",
  length(null_effects) == 15120L &&
    abs(human_effect$exact_p_two_sided - exact_p) < 1e-15 &&
    abs(mean(null_effects)) < 1e-15,
  sprintf("p=%.12g;null_mean=%.4g", exact_p, mean(null_effects))
)
add_check(
  "human_primary_and_sensitivities_opposite",
  human_coef[[1L]] > 0 && all(
    fread(file.path(candidate, "human_response/human_sensitivity.tsv"))$expected_direction_met == FALSE
  ),
  "positive improved-minus-not effect; expected negative"
)

validation <- rbindlist(checks)
fwrite(
  validation,
  file.path(candidate, "statistical_rederivation.tsv"),
  sep = "\t"
)
if (any(!validation$passed)) {
  print(validation[passed == FALSE])
  stop("Plan 44 statistical validation failed")
}
cat(sprintf("PLAN44_STATISTICS_VALIDATED\tchecks=%d\n", nrow(validation)))
