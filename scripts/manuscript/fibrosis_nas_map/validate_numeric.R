#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_dir <- dirname(normalizePath(script_arg))
source(file.path(script_dir, "config.R"))
source(file.path(script_dir, "analysis_lib.R"))

tol <- 1e-10
close <- function(observed, expected, label) {
  if (length(observed) != length(expected) ||
      any(abs(unname(observed) - unname(expected)) > tol, na.rm = TRUE) ||
      !identical(unname(is.na(observed)), unname(is.na(expected)))) {
    fail("Independent numeric validation failed: ", label)
  }
}

meta <- fread(file.path(discovery_root, "meta_effects.tsv"))
primary <- meta[
  view == "hotspot_program" & score_definition == "weighted_mean_z" &
    model_kind == "primary"
]
assert_true(nrow(primary) == expected_testable_programs * 2L,
            "Primary family size failed independent validation")
close(
  primary$q_value,
  p.adjust(primary$p_value_meta, method = "BH", n = nrow(primary)),
  "primary BH family"
)

cohort <- fread(file.path(discovery_root, "cohort_effects.tsv"))
cohort_primary <- cohort[
  view == "hotspot_program" & score_definition == "weighted_mean_z" &
    model_kind == "primary" & estimable == TRUE
]
close(
  cohort_primary$ci_lower,
  cohort_primary$beta - qt(0.975, cohort_primary$residual_df) * cohort_primary$se_hc3,
  "cohort HC3 lower confidence limits"
)
close(
  cohort_primary$ci_upper,
  cohort_primary$beta + qt(0.975, cohort_primary$residual_df) * cohort_primary$se_hc3,
  "cohort HC3 upper confidence limits"
)

loco <- fread(file.path(discovery_root, "validation_summary.tsv"))
loco_null <- readRDS(file.path(discovery_root, "loco_permutation_null.rds"))
rederived_global <- list()
for (axis_name in c("fibrosis", "nas")) {
  observed <- mean(atanh(pmax(pmin(
    loco[split == "leave_one_cohort_out" & axis == axis_name, estimate],
    0.999999), -0.999999)))
  null <- apply(loco_null[, , axis_name, drop = FALSE], 1L, function(x) {
    mean(atanh(pmax(pmin(x, 0.999999), -0.999999)))
  })
  stored <- loco[split == "leave_one_cohort_out_global" & axis == axis_name]
  close(stored$estimate, observed, paste(axis_name, "global Fisher-z"))
  p <- empirical_p(observed, null, "greater")
  close(stored$empirical_p, p, paste(axis_name, "global empirical p"))
  rederived_global[[axis_name]] <- p
}
global_p <- unlist(rederived_global[c("fibrosis", "nas")])
stored_holm <- loco[split == "leave_one_cohort_out_global"][
  match(names(global_p), axis), p_holm
]
close(stored_holm, p.adjust(global_p, method = "holm"), "global Holm correction")

linearity_contrasts <- fread(file.path(discovery_root, "linearity_contrasts.tsv"))
close(
  linearity_contrasts$q_value,
  p.adjust(
    linearity_contrasts$p_value, method = "BH",
    n = nrow(linearity_contrasts)
  ),
  "categorical adjacent-contrast BH family"
)
estimable_linearity <- linearity_contrasts[estimable == TRUE]
close(
  estimable_linearity$ci_lower,
  estimable_linearity$estimate -
    qt(0.975, estimable_linearity$residual_df) * estimable_linearity$se_hc3,
  "categorical adjacent-contrast HC3 lower confidence limits"
)
close(
  estimable_linearity$ci_upper,
  estimable_linearity$estimate +
    qt(0.975, estimable_linearity$residual_df) * estimable_linearity$se_hc3,
  "categorical adjacent-contrast HC3 upper confidence limits"
)

paired <- fread(file.path(holdout_root, "paired_validation_summary.tsv"))
participants <- fread(file.path(holdout_root, "paired_participant_validation.tsv"))
permutations <- readRDS(file.path(holdout_root, "paired_permutation_nulls.rds"))
bootstraps <- readRDS(file.path(holdout_root, "paired_bootstrap_distributions.rds"))
key <- "hotspot_program::weighted_mean_z"
stored_pair <- paired[view == "hotspot_program" & score_definition == "weighted_mean_z"]
participant_cosines <- participants[
  view == "hotspot_program" & score_definition == "weighted_mean_z" & changed == TRUE,
  cosine
]
observed_pair <- median(participant_cosines, na.rm = TRUE)
close(stored_pair$estimate, observed_pair, "paired median cosine")
close(
  c(stored_pair$ci_lower, stored_pair$ci_upper),
  as.numeric(quantile(bootstraps[[key]], c(0.025, 0.975), na.rm = TRUE)),
  "paired bootstrap confidence interval"
)
close(
  stored_pair$empirical_p,
  empirical_p(observed_pair, permutations[[key]], "greater"),
  "paired empirical p"
)

cat("Independent fibrosis-NAS numeric validation passed\n")
