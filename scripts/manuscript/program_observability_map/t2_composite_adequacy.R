# T2 supplement: is the NAS composite's equal weighting adequate for a program?
#
# PROVENANCE, STATED PLAINLY. This view was NOT in the contract of the sealed
# axis_map run. It was added after that run returned, and it is sealed in its own
# directory and reported as a supplement rather than promoted into the primary
# family. It is not a second attempt at the same hypothesis: the axis_map asks
# "which single component has a nonzero coefficient once the others are held",
# and at n = 58 with the components correlated at Spearman 0.48 to 0.70 that
# question is answerable for almost no program. The question here is different
# and better matched to what the substrate can support.
#
# THE QUESTION. NAS is defined as steatosis + ballooning + inflammation, which is
# an equal-weight sum. So the hypothesis
#
#     H0: beta_steatosis = beta_ballooning = beta_inflammation
#
# is exactly the claim that the composite's own weighting is adequate for that
# program - that nothing is lost by scoring the three together. Rejecting it says
# the program weights the components differently, which is the decomposition
# claim, and it is a 2-df test rather than three 1-df tests competing against
# each other's collinearity.
#
# Alongside it, a 3-df omnibus asks whether the three components jointly carry
# any association at all, without imposing equal weights. Comparing the omnibus
# with the 1-df composite test says what the equal-weighting assumption buys.
#
# Same contract as the axis_map: frozen programs, donors as the unit, the shared
# calibration gate before any count, an MDE on every negative, atomic sealing.

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(jsonlite)
})

script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file <- sub("^--file=", "", args[grepl("^--file=", args)])
  if (length(file)) return(dirname(normalizePath(file[[1]])))
  getwd()
}
here <- script_dir()
source(file.path(here, "config.R"))
source(file.path(here, "analysis_lib.R"))
source(file.path(here, "calibration_lib.R"))
source(file.path(here, "systems_contract.R"))
source(file.path(here, "t2_lib.R"))

system_id <- "T2_proteome"
run_id <- "composite_adequacy"
n_permutations <- as.integer(Sys.getenv("MASLD_T2_PERMUTATIONS", "2000"))
t2_seed <- 20260812L
primary_min_retained_l1 <- 0.50
primary_min_observed_genes <- 8L

liver_quant_path <- file.path(project_root, "data/PXD051911/liver_protein_quant.txt")
liver_meta_path <- file.path(project_root, "data/PXD051911/meta_data.txt")

# contrast_wald(), CONTRAST_OMNIBUS and CONTRAST_EQUAL_WEIGHT live in t2_lib.R
# so that t2_unit_tests.R can pin the identity contrast against wald_test().

joint_formula <- score ~ batch + sex_final + age_z + bmi_z +
  steatosis + ballooning + inflammation

fit_both_contrasts <- function(scores, meta) {
  d <- as.data.frame(meta)
  rbindlist(lapply(rownames(scores), function(program) {
    d$score <- scores[program, ]
    fit <- hc3_fit(d$score, d, joint_formula)
    omnibus <- contrast_wald(fit, CONTRAST_OMNIBUS)
    equal <- contrast_wald(fit, CONTRAST_EQUAL_WEIGHT)
    data.table(
      program_uid = program,
      estimable = isTRUE(fit$estimable),
      beta_steatosis = if (isTRUE(fit$estimable)) fit$coefficients[["steatosis"]] else NA_real_,
      beta_ballooning = if (isTRUE(fit$estimable)) fit$coefficients[["ballooning"]] else NA_real_,
      beta_inflammation = if (isTRUE(fit$estimable)) fit$coefficients[["inflammation"]] else NA_real_,
      omnibus_f = if (isTRUE(omnibus$estimable)) omnibus$f_statistic else NA_real_,
      omnibus_p = if (isTRUE(omnibus$estimable)) omnibus$p_value else NA_real_,
      equal_weight_f = if (isTRUE(equal$estimable)) equal$f_statistic else NA_real_,
      equal_weight_p = if (isTRUE(equal$estimable)) equal$p_value else NA_real_,
      residual_df = if (isTRUE(fit$estimable)) fit$df else NA_integer_,
      max_leverage = if (isTRUE(fit$estimable)) fit$max_leverage else NA_real_)
  }))
}

destination <- Sys.getenv("MASLD_T2_DEST",
                          file.path(workstream_root, "systems", system_id, run_id))
tmp <- atomic_dir(destination)
log_lines <- character(0)
say <- function(...) { line <- paste0(...); log_lines <<- c(log_lines, line); cat(line, "\n", sep = "") }

say("T2 supplement: adequacy of the NAS composite's equal weighting")
say("post-hoc addition to the sealed axis_map run; supplementary, not primary")
assert_inclusion_criterion(c("evidence_interpretation", "observability"))

frozen <- assert_frozen_programs(program_registry, program_membership, 117L)
registry <- frozen$registry
membership <- t2_membership_adapter(frozen$membership)

quant_header <- names(fread(liver_quant_path, nrows = 0L))
sample_columns <- setdiff(quant_header, c("ProteinAccessions", "Genes", "ProteinDescriptions"))
meta <- t2_build_liver_meta(liver_meta_path, sample_columns)
assert_biological_unit(meta, "patient_name")
assert_true(uniqueN(meta$patient_name) == nrow(meta), "A donor contributes more than one liver run")

proteins <- t2_build_protein_matrix(liver_quant_path, meta$sample_id)
coverage <- t2_apply_coverage_rule(
  t2_program_coverage(membership, rownames(proteins$scoring_universe), "scoring_complete_case"),
  primary_min_retained_l1, primary_min_observed_genes)
primary_programs <- coverage[testable == TRUE, program_uid]
scores <- score_programs(proteins$scoring_universe, meta, membership, registry,
                         coverage_threshold = primary_min_retained_l1)$primary[
                           primary_programs, , drop = FALSE]
assert_true(all(is.finite(scores)), "A program score is not finite")
say("programs: ", nrow(scores), "   donors: ", ncol(scores))

observed <- fit_both_contrasts(scores, meta)
assert_true(all(observed$estimable), "A joint contrast model is not estimable")

set.seed(t2_seed)
permuted <- lapply(seq_len(n_permutations), function(i) {
  permute_histology_within_cohort(copy(meta), T2_PERMUTED_COLUMNS, cohort_column = "batch")
})

calibrations <- list(); gates <- list()
for (test in c("omnibus", "equal_weight")) {
  p_column <- paste0(test, "_p")
  cal <- calibrate_view(
    observed_p = observed[[p_column]],
    permuted_p_fn = function(i) fit_both_contrasts(scores, permuted[[i]])[[p_column]],
    n_reps = n_permutations, n_tests = nrow(observed),
    label = test, observed_p_is_empirical = FALSE)
  calibrations[[test]] <- calibration_row(cal)
  gate <- gate_view(cal, strict = FALSE)
  gates[[test]] <- data.table(view = test, passed = gate$passed,
                              reasons = if (length(gate$reasons))
                                paste(gate$reasons, collapse = "; ") else NA_character_)
  say(sprintf("  %-13s tests=%d  null false calls %.2f (max %d)  gate=%s",
              test, nrow(observed), cal$null_call_mean,
              as.integer(cal$null_call_max), if (gate$passed) "PASS" else "FAIL"))
}
calibration <- rbindlist(calibrations)
gate_table <- rbindlist(gates)

observed[, omnibus_q := p.adjust(omnibus_p, "BH")]
observed[, equal_weight_q := p.adjust(equal_weight_p, "BH")]
reportable <- calibration[reportable == TRUE, view]
observed[, omnibus_reportable := "omnibus" %in% reportable]
observed[, equal_weight_reportable := "equal_weight" %in% reportable]

# The interpretation, per program, in one column. The gate check is a scalar
# property of the run, so it guards the assignment rather than sitting inside a
# vectorised branch alongside per-program columns.
both_gates_passed <- all(c("omnibus", "equal_weight") %in% reportable)
if (!both_gates_passed) {
  observed[, composite_verdict := "withheld_calibration_gate"]
} else {
  observed[, composite_verdict := fifelse(
    omnibus_q >= 0.05, "no_component_association_detected",
    fifelse(equal_weight_q < 0.05,
            "components_weighted_unequally_composite_inadequate",
            "equal_weighting_not_rejected_composite_adequate"))]
}
observed <- merge(observed,
                  registry[, .(program_uid, cell_type, module, module_name, module_top_pathway)],
                  by = "program_uid", all.x = TRUE)
observed <- merge(observed, coverage[, .(program_uid, retained_l1, n_observed, n_genes)],
                  by = "program_uid", all.x = TRUE)
setorder(observed, equal_weight_p)

say("verdicts: ", paste(sprintf("%s n=%d", names(table(observed$composite_verdict)),
                                table(observed$composite_verdict)), collapse = "   "))

assert_counts_calibrated(nrow(observed), calibration)

write_json(list(
  system_id = system_id, run_id = run_id, substrate = "PXD051911",
  provenance = paste("added after the sealed axis_map run returned; supplementary,",
                     "not part of the primary family and not a re-test of it"),
  hypotheses = list(
    omnibus = "H0: all three component coefficients are zero (3 df)",
    equal_weight = paste("H0: the three component coefficients are equal (2 df),",
                         "which is the claim that the NAS composite's own",
                         "equal-weight sum is adequate for this program")),
  n_donors = ncol(scores), n_programs = nrow(scores),
  biological_unit = "patient_name",
  coverage_rule = list(min_retained_l1 = primary_min_retained_l1,
                       min_observed_genes = primary_min_observed_genes),
  fdr_families = list(list(family = "omnibus", n_tests = nrow(observed)),
                      list(family = "equal_weight", n_tests = nrow(observed))),
  permutation = list(n_permutations = n_permutations,
                     scheme = "joint histologic permutation within acquisition batch"),
  seed = t2_seed),
  file.path(tmp, "analysis_contract.json"), auto_unbox = TRUE, pretty = TRUE)

write_tsv(observed, file.path(tmp, "composite_adequacy_results.tsv"))
write_tsv(calibration, file.path(tmp, "calibration.tsv"))
write_tsv(gate_table, file.path(tmp, "gate_status.tsv"))
writeLines(log_lines, file.path(tmp, "run_log.txt"))
capture.output(print(sessionInfo()), file = file.path(tmp, "sessionInfo.txt"))
writeLines(sort(system2("sha256sum",
                        c(file.path(here, "t2_lib.R"),
                          file.path(here, "t2_composite_adequacy.R"),
                          file.path(here, "analysis_lib.R"),
                          file.path(here, "calibration_lib.R"),
                          liver_quant_path, liver_meta_path,
                          program_registry, program_membership), stdout = TRUE)),
           file.path(tmp, "input_manifest.sha256"))
writeLines(sort(system2("sha256sum", list.files(tmp, full.names = TRUE), stdout = TRUE)),
           file.path(tmp, "output_manifest.sha256"))
for (f in list.files(tmp, full.names = TRUE)) Sys.chmod(f, "0440")
publish_dir(tmp, destination)
Sys.chmod(destination, "0550")
cat("\nSEALED: ", destination, "\n", sep = "")
