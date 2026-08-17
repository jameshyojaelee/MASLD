# T2: the 117 frozen Hotspot programs on the liver proteome, against steatosis,
# ballooning and lobular inflammation separated.
#
# WHAT THIS SYSTEM IS FOR
# The bulk RNA substrate records the NAF LD activity score only as a composite, so
# on that substrate a program that tracks fat and a program that tracks lobular
# inflammation are indistinguishable. PXD051911 carries donor-level Kleiner
# fibrosis, NAS AND all three NAS components alongside BMI, age and sex on the
# same 58 donors, which makes it the only substrate in this repository where the
# composite can be taken apart. This is a cross-sectional association between a
# program score and a recorded histologic grade. It is not an ordering.
#
# WHAT IT IS NOT FOR
# Fibrosis. The liver arm is Kleiner F0 13 / F1 31 / F2 11 / F3 3 / F4 0. The
# fibrosis view is carried as a descriptive companion with its power stated, and
# it is not part of the primary family.
#
# Inclusion criteria (PAPER.md:211-214): coverage, observability,
# evidence_interpretation.

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

# ===========================================================================
# DECLARED BEFORE THE RUN. Nothing below this block is chosen after seeing a
# result, and the block is written to the contract file before the first model
# is fitted.
# ===========================================================================

system_id <- "T2_proteome"
run_id <- "axis_map"
substrate_id <- "PXD051911"

# COVERAGE RULE.
# A program score is a weighted mean of the z-scored proteins the program can
# see. Below some retained fraction of the program's own L1 weight, that mean
# stops measuring the program and starts measuring an arbitrary subset of it.
# The rule is stated against the SCORING universe - proteins quantified in all
# 58 donors - because that is the set the score is actually built from, not
# against the wider source universe where a protein can be counted as covered
# and still contribute to nobody's score.
primary_min_retained_l1 <- 0.50
primary_min_observed_genes <- 8L
sensitivity_min_retained_l1 <- 0.30
sensitivity_min_observed_genes <- 8L

# EFFECT THRESHOLD. Reused unchanged from config.R
# (observability_effect_threshold = 0.20), where it is defined as 0.20 program-
# score SD per one-unit increment of the recorded axis. A Kleiner component
# grade is such an axis, so the scale carries over without rescaling.
effect_threshold <- observability_effect_threshold

# COVARIATES. Batch is mandatory, not optional: the two acquisitions differ in
# per-donor detection depth at p = 4.9e-08 and in ballooning grade at Wilcoxon
# p = 0.012, so it is a confounder of the very contrast this system reports.
# BMI because the series is bariatric-enriched and adiposity is entangled with
# steatosis; age and sex as standing donor covariates.

# FDR FAMILIES. Declared explicitly because the retired 22-program pipeline used
# the scope string "17_testable_programs_x_5_histology_features" and inheriting
# it at 117 programs would silently misstate every padj. Each view is its own BH
# family, each family is calibrated separately, and no count is quoted from a
# view whose gate fails.
view_spec <- list(
  list(view = "joint_components", tier = "primary", role = "primary",
       terms = T2_HISTOLOGY_COMPONENTS, covariates = T2_COVARIATE_TERMS,
       description = "steatosis, ballooning and lobular inflammation mutually adjusted"),
  list(view = "marginal_steatosis", tier = "primary", role = "marginal",
       terms = "steatosis", covariates = T2_COVARIATE_TERMS,
       description = "steatosis alone, not adjusted for the other two components"),
  list(view = "marginal_ballooning", tier = "primary", role = "marginal",
       terms = "ballooning", covariates = T2_COVARIATE_TERMS,
       description = "ballooning alone, not adjusted for the other two components"),
  list(view = "marginal_inflammation", tier = "primary", role = "marginal",
       terms = "inflammation", covariates = T2_COVARIATE_TERMS,
       description = "lobular inflammation alone, not adjusted for the other two"),
  list(view = "nas_composite", tier = "primary", role = "composite",
       terms = "nas", covariates = T2_COVARIATE_TERMS,
       description = "the composite activity score, for comparison with its parts"),
  list(view = "fibrosis_descriptive", tier = "primary", role = "descriptive",
       terms = "fibrosis", covariates = T2_COVARIATE_TERMS,
       description = "Kleiner fibrosis grade; F4 absent and F3 n = 3, descriptive only"),
  list(view = "masld_vs_control", tier = "primary", role = "diagnosis",
       terms = "is_masld", covariates = T2_COVARIATE_TERMS,
       description = "SAF diagnosis MASL or MASH versus No_MASLD, 46 versus 12"),
  list(view = "joint_components_coverage30", tier = "sensitivity", role = "sensitivity",
       terms = T2_HISTOLOGY_COMPONENTS, covariates = T2_COVARIATE_TERMS,
       description = "the primary joint model on the wider L1 >= 0.30 program set"),
  list(view = "joint_ballooning_binary", tier = "primary", role = "sensitivity",
       terms = c("steatosis", "ballooning_binary", "inflammation"),
       covariates = T2_COVARIATE_TERMS,
       description = "ballooning as present/absent; grade 2 has n = 1 donor"),
  list(view = "joint_components_no_batch", tier = "primary", role = "diagnostic",
       terms = T2_HISTOLOGY_COMPONENTS,
       covariates = setdiff(T2_COVARIATE_TERMS, "batch"),
       description = "batch deliberately omitted, to measure what batch adjustment costs")
)

# The marginal views are one BH family jointly, not three, because they answer
# one question with three tests.
family_of_view <- c(
  joint_components = "joint_components",
  marginal_steatosis = "marginal_components",
  marginal_ballooning = "marginal_components",
  marginal_inflammation = "marginal_components",
  nas_composite = "nas_composite",
  fibrosis_descriptive = "fibrosis_descriptive",
  masld_vs_control = "masld_vs_control",
  joint_components_coverage30 = "joint_components_coverage30",
  joint_ballooning_binary = "joint_ballooning_binary",
  joint_components_no_batch = "joint_components_no_batch"
)

n_calibration_permutations <- as.integer(
  Sys.getenv("MASLD_T2_PERMUTATIONS", "2000"))
t2_seed <- 20260812L

liver_quant_path <- file.path(project_root, "data/PXD051911/liver_protein_quant.txt")
liver_meta_path <- file.path(project_root, "data/PXD051911/meta_data.txt")
plasma_quant_path <- file.path(project_root, "data/PXD051911/plasma_protein_quant.txt")
plasma_bc_path <- file.path(
  project_root, "data/PXD051911/plasma_protein_quant_filtered_and_batch_corrected.txt")

# The sealed destination. MASLD_T2_DEST exists so the smoke run can write to
# scratch; the real run leaves it unset and lands under the workstream.
destination <- Sys.getenv("MASLD_T2_DEST",
                          file.path(workstream_root, "systems", system_id, run_id))
tmp <- atomic_dir(destination)
log_lines <- character(0)
say <- function(...) {
  line <- paste0(...)
  log_lines <<- c(log_lines, line)
  cat(line, "\n", sep = "")
}

say("T2 liver-proteome axis map")
say("destination: ", destination)
say("permutations: ", n_calibration_permutations, "   seed: ", t2_seed)

assert_inclusion_criterion(c("coverage", "observability", "evidence_interpretation"))
for (spec in view_spec) assert_language(spec$description)

# ===========================================================================
# 1. Frozen programs.
# ===========================================================================
frozen <- assert_frozen_programs(program_registry, program_membership, 117L)
registry <- frozen$registry
membership <- t2_membership_adapter(frozen$membership)
say("frozen programs: ", nrow(registry),
    "   membership rows after adapter: ", nrow(membership),
    "   distinct symbols: ", uniqueN(membership$mapped_symbol))
assert_true(uniqueN(membership$program_uid) == 117L,
            "The membership adapter lost a program; the schema trap has fired")

# ===========================================================================
# 2. Substrate.
# ===========================================================================
quant_header <- names(fread(liver_quant_path, nrows = 0L))
sample_columns <- setdiff(quant_header,
                          c("ProteinAccessions", "Genes", "ProteinDescriptions"))
meta <- t2_build_liver_meta(liver_meta_path, sample_columns)

# The biological unit. patient_name is a distinct column from sample_id, so the
# uniqueness check is not a column compared with itself.
assert_biological_unit(meta, "patient_name")
assert_true(!identical(as.character(meta$patient_name), as.character(meta$sample_id)),
            "patient_name was filled from sample_id; the unit check would be vacuous")
assert_true(uniqueN(meta$patient_name) == nrow(meta),
            "A donor contributes more than one liver run")
say("donors: ", uniqueN(meta$patient_name), " (one liver run each, all initial_sample)")
say("batch: ", paste(sprintf("%s n=%d", names(table(meta$batch)), table(meta$batch)),
                     collapse = "  "))

proteins <- t2_build_protein_matrix(liver_quant_path, meta$sample_id)
say("protein groups deposited: ", proteins$stats$n_protein_groups_deposited,
    "   unambiguous single-symbol: ", proteins$stats$n_protein_groups_unambiguous)
say("matrix missingness: ", sprintf("%.2f%%", 100 * proteins$stats$matrix_na_fraction))
say("source universe (<= 50% NA): ", proteins$stats$n_symbols_source_universe, " symbols")
say("scoring universe (complete in all 58): ",
    proteins$stats$n_symbols_scoring_universe, " symbols")

# The evidence for the complete-case decision, recorded rather than asserted.
meta[, detection_na_fraction := proteins$per_donor_na_fraction[sample_id]]
depth_rows <- rbindlist(lapply(
  c("steatosis", "ballooning", "inflammation", "nas", "fibrosis", "bmi_z", "age_z"),
  function(v) {
    ct <- suppressWarnings(stats::cor.test(meta$detection_na_fraction, meta[[v]],
                                           method = "spearman"))
    data.table(axis = v, spearman_rho = unname(ct$estimate), p_value = ct$p.value)
  }))
depth_rows <- rbind(depth_rows, data.table(
  axis = "batch",
  spearman_rho = NA_real_,
  p_value = suppressWarnings(
    stats::wilcox.test(detection_na_fraction ~ batch, data = meta)$p.value)))
say("per-donor detection depth vs batch: p = ",
    signif(depth_rows[axis == "batch", p_value], 3),
    "; vs ballooning: rho = ", signif(depth_rows[axis == "ballooning", spearman_rho], 3),
    ", p = ", signif(depth_rows[axis == "ballooning", p_value], 3))

# ===========================================================================
# 3. Observability of all 117 programs, under both universes.
# ===========================================================================
symbols_source <- rownames(proteins$source_universe)
symbols_scoring <- rownames(proteins$scoring_universe)
coverage_source <- t2_program_coverage(membership, symbols_source, "source_le50pct_na")
coverage_scoring <- t2_program_coverage(membership, symbols_scoring, "scoring_complete_case")

coverage_scoring <- t2_apply_coverage_rule(coverage_scoring,
                                           primary_min_retained_l1,
                                           primary_min_observed_genes)
coverage_scoring[, testable_sensitivity :=
                   is.finite(retained_l1) & retained_l1 >= sensitivity_min_retained_l1 &
                   n_observed >= sensitivity_min_observed_genes]
primary_programs <- coverage_scoring[testable == TRUE, program_uid]
sensitivity_programs <- coverage_scoring[testable_sensitivity == TRUE, program_uid]
say("median retained L1, source universe: ",
    sprintf("%.3f", median(coverage_source$retained_l1)))
say("median retained L1, scoring universe: ",
    sprintf("%.3f", median(coverage_scoring$retained_l1)))
say("PRIMARY testable programs (L1 >= ", primary_min_retained_l1, ", n >= ",
    primary_min_observed_genes, "): ", length(primary_programs), " of 117")
say("SENSITIVITY testable programs (L1 >= ", sensitivity_min_retained_l1, "): ",
    length(sensitivity_programs), " of 117")
assert_true(length(primary_programs) > 0L,
            "No program clears the declared coverage rule; there is nothing to test")

# ===========================================================================
# 4. Contract, written before any model is fitted.
# ===========================================================================
contract <- list(
  system_id = system_id, run_id = run_id, substrate = substrate_id,
  release_id = release_id, workstream_id = workstream_id,
  question = paste("which recorded histologic feature each frozen Hotspot program",
                   "tracks on liver protein when steatosis, ballooning and lobular",
                   "inflammation are separated"),
  inclusion_criteria = c("coverage", "observability", "evidence_interpretation"),
  biological_unit = "patient_name",
  n_donors = nrow(meta),
  frozen_programs = 117L,
  coverage_rule = list(
    stated_against = "scoring universe: proteins quantified in all 58 donors",
    primary = list(min_retained_l1 = primary_min_retained_l1,
                   min_observed_genes = primary_min_observed_genes,
                   n_testable = length(primary_programs)),
    sensitivity = list(min_retained_l1 = sensitivity_min_retained_l1,
                       min_observed_genes = sensitivity_min_observed_genes,
                       n_testable = length(sensitivity_programs))),
  covariates = T2_COVARIATE_TERMS,
  effect_threshold = effect_threshold,
  effect_units = "program-score SD per one-unit increment of the recorded grade",
  fdr_families = lapply(unique(unname(family_of_view)), function(f) {
    views <- names(family_of_view)[family_of_view == f]
    tiers <- vapply(view_spec[vapply(view_spec, function(s) s$view %in% views, logical(1))],
                    function(s) s$tier, character(1))
    n_programs <- if (any(grepl("coverage30", views))) length(sensitivity_programs)
                  else length(primary_programs)
    n_terms <- sum(vapply(view_spec[vapply(view_spec, function(s) s$view %in% views,
                                           logical(1))],
                          function(s) length(s$terms), integer(1)))
    list(family = f, views = views, tier = unname(tiers[[1]]),
         n_programs = n_programs, n_terms = n_terms,
         n_tests = n_programs * n_terms)
  }),
  superseded_scope_string = "17_testable_programs_x_5_histology_features",
  permutation = list(
    n_permutations = n_calibration_permutations,
    scheme = paste("joint permutation of the histologic block within acquisition",
                   "batch, so the batch-histology association and the mutual",
                   "correlation of the axes are both preserved"),
    permuted_columns = T2_PERMUTED_COLUMNS,
    observed_p_is_empirical = FALSE),
  seed = t2_seed,
  fibrosis_note = paste("Kleiner F0 13 / F1 31 / F2 11 / F3 3 / F4 0; the fibrosis",
                        "view is descriptive and underpowered by construction"),
  missingness_decision = paste("program scores are built only from proteins",
                               "quantified in all 58 donors, because per-donor",
                               "detection depth differs by batch and correlates",
                               "with ballooning and NAS")
)
write_json(contract, file.path(tmp, "analysis_contract.json"),
           auto_unbox = TRUE, pretty = TRUE)
say("contract written before fitting: ", file.path(tmp, "analysis_contract.json"))

# ===========================================================================
# 5. Program scores.
# ===========================================================================
score_all <- score_programs(proteins$scoring_universe, meta, membership, registry,
                            coverage_threshold = primary_min_retained_l1)
scores_primary <- score_all$primary[primary_programs, , drop = FALSE]
assert_true(all(is.finite(scores_primary)),
            "A primary program score is not finite after complete-case scoring")
assert_true(all(abs(apply(scores_primary, 1L, stats::sd) - 1) < 1e-8),
            "Program scores are not on the unit-SD scale the effect threshold assumes")

score_sens <- score_programs(proteins$scoring_universe, meta, membership, registry,
                             coverage_threshold = sensitivity_min_retained_l1)
scores_sensitivity <- score_sens$primary[sensitivity_programs, , drop = FALSE]
assert_true(all(is.finite(scores_sensitivity)), "A sensitivity program score is not finite")
say("scored: ", nrow(scores_primary), " primary and ", nrow(scores_sensitivity),
    " sensitivity programs across ", ncol(scores_primary), " donors")

# How many independent things are being counted?
dimensionality <- rbind(
  cbind(set = "primary", t2_score_dimensionality(scores_primary)),
  cbind(set = "sensitivity", t2_score_dimensionality(scores_sensitivity)))
say("primary score matrix: PC1 explains ",
    sprintf("%.1f%%", 100 * dimensionality[set == "primary", pc1_variance_explained]),
    ", effective dimension ",
    sprintf("%.1f", dimensionality[set == "primary",
                                   effective_dimension_participation_ratio]),
    " of ", nrow(scores_primary), " programs")

# ===========================================================================
# 6. Shared permutations, then one fit per view.
# ===========================================================================
set.seed(t2_seed)
permuted_meta <- lapply(seq_len(n_calibration_permutations), function(i) {
  permute_histology_within_cohort(copy(meta), T2_PERMUTED_COLUMNS,
                                  cohort_column = "batch")
})

fit_rows <- list()
calibration_rows <- list()
gate_status <- list()

for (spec in view_spec) {
  Y <- if (identical(spec$tier, "sensitivity")) scores_sensitivity else scores_primary
  fitted <- t2_fit_view(Y, meta, spec$terms, spec$view, spec$covariates)
  rows <- fitted$rows
  rows[, `:=`(tier = spec$tier, role = spec$role,
              family = unname(family_of_view[[spec$view]]),
              n_programs = nrow(Y))]
  fit_rows[[spec$view]] <- rows

  permuted_p_fn <- function(i) {
    t2_view_p(Y, permuted_meta[[i]], spec$terms, spec$covariates)
  }
  cal <- calibrate_view(
    observed_p = rows$p_value, permuted_p_fn = permuted_p_fn,
    n_reps = n_calibration_permutations, n_tests = nrow(rows),
    label = spec$view, observed_p_is_empirical = FALSE)
  calibration_rows[[spec$view]] <- calibration_row(cal)
  gate <- gate_view(cal, strict = FALSE)
  gate_status[[spec$view]] <- data.table(
    view = spec$view, family = unname(family_of_view[[spec$view]]),
    tier = spec$tier, passed = gate$passed,
    reasons = if (length(gate$reasons)) paste(gate$reasons, collapse = "; ") else NA_character_)
  say(sprintf("  %-28s tests=%3d  null false calls %.2f (max %d)  gate=%s",
              spec$view, nrow(rows), cal$null_call_mean,
              as.integer(cal$null_call_max), if (gate$passed) "PASS" else "FAIL"))
}

results <- rbindlist(fit_rows, use.names = TRUE)
calibration <- rbindlist(calibration_rows, use.names = TRUE)
calibration <- merge(calibration,
                     rbindlist(gate_status)[, .(view, family, tier)], by = "view")
gates <- rbindlist(gate_status)

# BH within the declared family, not within a view and not across everything.
results[, q_value := p.adjust(p_value, method = "BH"), by = family]
results <- t2_classify(results, effect_threshold)
assert_negatives_have_mde(results, "support", "mde")

# The gate decides what may be counted. A failing view keeps its estimates, so
# the failure itself is inspectable, but its calls are not reportable.
reportable_views <- calibration[reportable == TRUE, view]
results[, reportable := view %in% reportable_views]
results[reportable == FALSE, support := "withheld_calibration_gate"]

say("views passing the calibration gate: ", length(reportable_views), " of ",
    nrow(calibration))
if (nrow(calibration[reportable == FALSE])) {
  say("WITHHELD: ", paste(calibration[reportable == FALSE, view], collapse = ", "))
}

# ===========================================================================
# 7. Which feature does each program track?
# ===========================================================================
joint <- results[view == "joint_components"]
assignment <- if ("joint_components" %in% reportable_views) {
  t2_assign_feature(joint)
} else {
  data.table(program_uid = character(0))
}
if (nrow(assignment)) {
  assignment <- merge(assignment,
                      registry[, .(program_uid, cell_type, module, module_name,
                                   module_top_pathway)],
                      by = "program_uid", all.x = TRUE)
  assignment <- merge(assignment,
                      coverage_scoring[, .(program_uid, retained_l1, n_observed, n_genes)],
                      by = "program_uid", all.x = TRUE)
  assignment <- merge(
    assignment,
    results[view == "nas_composite", .(program_uid, nas_beta = beta,
                                       nas_q = q_value, nas_support = support)],
    by = "program_uid", all.x = TRUE)
  assignment[, resolution := fifelse(
    !is.na(tracks), "component_resolved",
    fifelse(nas_support == "supported", "composite_only_no_component_resolved",
            fifelse(n_informative_null == 3L, "informative_null_on_all_three",
                    "underpowered_on_at_least_one_component")))]
  setorder(assignment, -n_supported, program_uid)
  say("programs tracking at least one separated component: ",
      assignment[n_supported > 0L, .N], " of ", nrow(assignment))
  say("  by feature: ", paste(sprintf("%s n=%d",
        names(table(assignment$tracks)), table(assignment$tracks)), collapse = "  "))
  say("  composite-supported but no component resolved: ",
      assignment[resolution == "composite_only_no_component_resolved", .N])
  say("  informative null on all three components: ",
      assignment[resolution == "informative_null_on_all_three", .N])
}

# ===========================================================================
# 8. Power, stated per view rather than asserted.
# ===========================================================================
power_summary <- results[, .(
  n_tests = .N,
  n_supported = sum(support == "supported"),
  n_informative_null = sum(support == "informative_null"),
  n_indeterminate = sum(support == "indeterminate"),
  median_mde_native = median(mde, na.rm = TRUE),
  median_mde_per_axis_sd = median(mde_per_axis_sd, na.rm = TRUE),
  min_mde_per_axis_sd = min(mde_per_axis_sd, na.rm = TRUE),
  max_mde_per_axis_sd = max(mde_per_axis_sd, na.rm = TRUE),
  pct_tests_powered_for_threshold =
    100 * mean(mde_per_axis_sd <= effect_threshold, na.rm = TRUE)
), by = .(view, family, tier, role, reportable)]
setorder(power_summary, -reportable, view)

# ===========================================================================
# 9. Liver-and-plasma donor overlap. Two numbers have been quoted; they answer
# different questions and only one of them is the paired count.
# ===========================================================================
meta_all <- fread(liver_meta_path)
plasma_header <- names(fread(plasma_quant_path, nrows = 0L))
plasma_columns <- setdiff(plasma_header,
                          c("ProteinAccessions", "Genes", "ProteinDescriptions"))
plasma_bc_header <- names(fread(plasma_bc_path, nrows = 0L))
plasma_bc_columns <- setdiff(plasma_bc_header,
                             c("ProteinAccessions", "Genes", "ProteinDescriptions"))
liver_donors <- unique(meta$patient_name)
paired_same_visit <- intersect(
  liver_donors,
  meta_all[plasma_proteomics_filename %in% plasma_columns &
             sample_group == "initial_sample", patient_name])
any_visit <- intersect(
  liver_donors, meta_all[plasma_proteomics_filename %in% plasma_columns, patient_name])
later_visit_only <- setdiff(any_visit, paired_same_visit)

plasma_raw <- fread(plasma_quant_path, select = c("ProteinAccessions", "Genes"))
plasma_symbols <- unique(toupper(trimws(
  plasma_raw[!is.na(Genes) & Genes != "" & !grepl(";", Genes, fixed = TRUE), Genes])))
shared_symbols <- intersect(plasma_symbols, symbols_source)
coverage_plasma <- t2_program_coverage(membership, plasma_symbols, "plasma_all_symbols")

overlap <- data.table(
  quantity = c("liver donors",
               "liver donor with baseline plasma from the SAME visit",
               "liver donor with plasma at ANY visit",
               "of those, plasma only from the later bariatric-surgery visit",
               "plasma sample columns (raw)",
               "plasma sample columns (batch-corrected)",
               "plasma donors overall",
               "plasma symbols",
               "symbols shared with the liver source universe",
               "programs at plasma L1 >= 0.50",
               "programs at plasma L1 >= 0.30"),
  value = c(length(liver_donors), length(paired_same_visit), length(any_visit),
            length(later_visit_only), length(plasma_columns), length(plasma_bc_columns),
            uniqueN(meta_all[plasma_proteomics_filename %in% plasma_columns, patient_name]),
            length(plasma_symbols), length(shared_symbols),
            coverage_plasma[retained_l1 >= 0.50, .N],
            coverage_plasma[retained_l1 >= 0.30, .N]))
say("liver-and-plasma donors: ", length(paired_same_visit),
    " paired at the same visit; ", length(any_visit),
    " have plasma at any visit (the extra ", length(later_visit_only),
    " contribute plasma only from the later bariatric-surgery visit, whose",
    " histology is not the liver biopsy's)")
say("plasma program coverage: ", coverage_plasma[retained_l1 >= 0.50, .N],
    " of 117 programs reach L1 >= 0.50; plasma cannot carry the program map")

# ===========================================================================
# 10. Seal.
# ===========================================================================
donor_table <- meta[, .(patient_name, sample_id, batch, sex_final, age_z, bmi_z,
                        steatosis, ballooning, inflammation, nas, fibrosis,
                        is_masld, saf_diagnosis, detection_na_fraction)]
score_table <- as.data.table(scores_primary, keep.rownames = "program_uid")

write_tsv(results, file.path(tmp, "program_histology_results.tsv"))
write_tsv(calibration, file.path(tmp, "calibration.tsv"))
write_tsv(gates, file.path(tmp, "gate_status.tsv"))
write_tsv(power_summary, file.path(tmp, "power_summary.tsv"))
write_tsv(assignment, file.path(tmp, "program_feature_assignment.tsv"))
write_tsv(coverage_scoring, file.path(tmp, "program_coverage_scoring_universe.tsv"))
write_tsv(coverage_source, file.path(tmp, "program_coverage_source_universe.tsv"))
write_tsv(coverage_plasma, file.path(tmp, "program_coverage_plasma.tsv"))
write_tsv(score_table, file.path(tmp, "program_scores_primary.tsv"))
write_tsv(donor_table, file.path(tmp, "donor_table.tsv"))
write_tsv(depth_rows, file.path(tmp, "detection_depth_confounding.tsv"))
write_tsv(dimensionality, file.path(tmp, "score_dimensionality.tsv"))
write_tsv(overlap, file.path(tmp, "liver_plasma_overlap.tsv"))
write_tsv(proteins$stats, file.path(tmp, "substrate_stats.tsv"))
writeLines(log_lines, file.path(tmp, "run_log.txt"))

# The contract's own assertion: no count without its calibration row, and no
# count from a view whose gate failed. Reported counts come only from the
# reportable views, so this is asserted on exactly those.
assert_counts_calibrated(results[reportable == TRUE, .N],
                         calibration[reportable == TRUE])

capture.output(print(sessionInfo()), file = file.path(tmp, "sessionInfo.txt"))
writeLines(sort(system2("sha256sum",
                        c(file.path(here, "t2_lib.R"),
                          file.path(here, "t2_proteome_axis_map.R"),
                          file.path(here, "t2_unit_tests.R"),
                          file.path(here, "analysis_lib.R"),
                          file.path(here, "calibration_lib.R"),
                          file.path(here, "systems_contract.R"),
                          liver_quant_path, liver_meta_path,
                          program_registry, program_membership),
                        stdout = TRUE)),
           file.path(tmp, "input_manifest.sha256"))
outputs <- list.files(tmp, full.names = TRUE)
writeLines(sort(system2("sha256sum", outputs, stdout = TRUE)),
           file.path(tmp, "output_manifest.sha256"))

for (f in list.files(tmp, full.names = TRUE)) Sys.chmod(f, "0440")
publish_dir(tmp, destination)
Sys.chmod(destination, "0550")
cat("\nSEALED: ", destination, "\n", sep = "")
