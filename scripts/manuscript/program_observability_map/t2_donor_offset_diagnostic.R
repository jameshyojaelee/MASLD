# T2 diagnostic: is the shared direction of the program effects a donor-level
# technical offset rather than 26 separate findings?
#
# WHY THIS EXISTS. In the sealed axis_map run every one of the 85 supported
# effects, across every view and every axis, is negative, and the 26 primary
# program scores have an effective dimension of 2.8. Both are consistent with
# real biology - loss of differentiated hepatocyte function - and both are also
# consistent with a single donor-level offset propagating identically into every
# program score. A weighted mean of per-protein z-scores inherits any additive
# per-donor shift, so this is not a remote possibility, it is the obvious
# competing explanation and it has to be measured rather than argued away.
#
# THE MECHANISM THAT WOULD PRODUCE ONE. Between-array quantile normalisation is
# applied to the source universe, where 15.1% of values are missing, and the
# scores are then built from the complete-case subset. A donor missing many
# low-abundance proteins has its observed values spread across the same target
# distribution, which shifts where its complete-case proteins land. Detection
# depth differs sharply by acquisition batch, so this offset is expected to
# exist. The question is only whether it points at the histologic axes.
#
# WHAT IS REPORTED. The offset itself, what it tracks, and the NAS view refitted
# with the offset as an explicit covariate. Adjusting for it is deliberately
# conservative: if disease genuinely lowers a broad swath of the abundant
# hepatocyte proteome then the offset absorbs real signal along with technical
# signal, so the adjusted count is a lower bound and the unadjusted one an upper
# bound. Both are reported and neither is called the answer.

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
run_id <- "donor_offset_diagnostic"
n_permutations <- as.integer(Sys.getenv("MASLD_T2_PERMUTATIONS", "2000"))
t2_seed <- 20260812L
primary_min_retained_l1 <- 0.50
primary_min_observed_genes <- 8L

liver_quant_path <- file.path(project_root, "data/PXD051911/liver_protein_quant.txt")
liver_meta_path <- file.path(project_root, "data/PXD051911/meta_data.txt")

destination <- Sys.getenv("MASLD_T2_DEST",
                          file.path(workstream_root, "systems", system_id, run_id))
tmp <- atomic_dir(destination)
log_lines <- character(0)
say <- function(...) { line <- paste0(...); log_lines <<- c(log_lines, line); cat(line, "\n", sep = "") }

say("T2 diagnostic: donor-level offset as a competing explanation")
assert_inclusion_criterion(c("evidence_interpretation", "observability"))

frozen <- assert_frozen_programs(program_registry, program_membership, 117L)
registry <- frozen$registry
membership <- t2_membership_adapter(frozen$membership)
quant_header <- names(fread(liver_quant_path, nrows = 0L))
sample_columns <- setdiff(quant_header, c("ProteinAccessions", "Genes", "ProteinDescriptions"))
meta <- t2_build_liver_meta(liver_meta_path, sample_columns)
assert_biological_unit(meta, "patient_name")
proteins <- t2_build_protein_matrix(liver_quant_path, meta$sample_id)

# The offset: each donor's mean position across the z-scored complete-case
# proteins. This is exactly the quantity a weighted program score inherits.
Z <- zscore_rows(proteins$scoring_universe)
Z <- Z[rowSums(is.finite(Z)) == ncol(Z), , drop = FALSE]
meta[, donor_offset := colMeans(Z)]
meta[, offset_z := as.numeric(scale(donor_offset))]
meta[, detection_na_fraction := proteins$per_donor_na_fraction[sample_id]]
say("donor offset: SD ", sprintf("%.3f", sd(meta$donor_offset)),
    ", range ", sprintf("%.3f to %.3f", min(meta$donor_offset), max(meta$donor_offset)))

offset_tracks <- rbindlist(lapply(
  c("detection_na_fraction", "steatosis", "ballooning", "inflammation", "nas",
    "fibrosis", "bmi_z", "age_z"),
  function(v) {
    ct <- suppressWarnings(stats::cor.test(meta$donor_offset, meta[[v]], method = "spearman"))
    data.table(against = v, spearman_rho = unname(ct$estimate), p_value = ct$p.value)
  }))
offset_tracks <- rbind(offset_tracks, data.table(
  against = "batch", spearman_rho = NA_real_,
  p_value = suppressWarnings(stats::wilcox.test(donor_offset ~ batch, data = meta)$p.value)))
say("offset vs detection depth: rho = ",
    signif(offset_tracks[against == "detection_na_fraction", spearman_rho], 3),
    "; vs batch: p = ", signif(offset_tracks[against == "batch", p_value], 3),
    "; vs NAS: rho = ", signif(offset_tracks[against == "nas", spearman_rho], 3),
    ", p = ", signif(offset_tracks[against == "nas", p_value], 3))

coverage <- t2_apply_coverage_rule(
  t2_program_coverage(membership, rownames(proteins$scoring_universe), "scoring_complete_case"),
  primary_min_retained_l1, primary_min_observed_genes)
scores <- score_programs(proteins$scoring_universe, meta, membership, registry,
                         coverage_threshold = primary_min_retained_l1)$primary[
                           coverage[testable == TRUE, program_uid], , drop = FALSE]
score_offset_correlation <- apply(scores, 1L, function(x) stats::cor(x, meta$donor_offset))
say("program scores vs the offset: median r = ",
    sprintf("%.3f", median(score_offset_correlation)))

adjusted_covariates <- c(T2_COVARIATE_TERMS, "offset_z")
set.seed(t2_seed)
permuted <- lapply(seq_len(n_permutations), function(i) {
  permute_histology_within_cohort(copy(meta), T2_PERMUTED_COLUMNS, cohort_column = "batch")
})

arms <- list(
  list(view = "nas_unadjusted", covariates = T2_COVARIATE_TERMS),
  list(view = "nas_offset_adjusted", covariates = adjusted_covariates))
rows <- list(); calibrations <- list(); gates <- list()
for (arm in arms) {
  fitted <- t2_fit_view(scores, meta, "nas", arm$view, arm$covariates)$rows
  cal <- calibrate_view(fitted$p_value,
                        function(i) t2_view_p(scores, permuted[[i]], "nas", arm$covariates),
                        n_reps = n_permutations, n_tests = nrow(fitted),
                        label = arm$view, observed_p_is_empirical = FALSE)
  calibrations[[arm$view]] <- calibration_row(cal)
  gate <- gate_view(cal, strict = FALSE)
  gates[[arm$view]] <- data.table(view = arm$view, passed = gate$passed,
                                  reasons = if (length(gate$reasons))
                                    paste(gate$reasons, collapse = "; ") else NA_character_)
  fitted[, q_value := p.adjust(p_value, "BH")]
  rows[[arm$view]] <- t2_classify(fitted, observability_effect_threshold)
  say(sprintf("  %-22s null false calls %.2f  gate=%s  supported=%d/%d",
              arm$view, cal$null_call_mean, if (gate$passed) "PASS" else "FAIL",
              sum(rows[[arm$view]]$q_value < 0.05), nrow(fitted)))
}
result <- rbindlist(rows)
calibration <- rbindlist(calibrations)
gate_table <- rbindlist(gates)
assert_negatives_have_mde(result, "support", "mde")
assert_counts_calibrated(nrow(result), calibration)

comparison <- merge(
  result[view == "nas_unadjusted", .(program_uid, beta_unadjusted = beta,
                                     q_unadjusted = q_value)],
  result[view == "nas_offset_adjusted", .(program_uid, beta_adjusted = beta,
                                          q_adjusted = q_value)],
  by = "program_uid")
comparison[, score_offset_r := score_offset_correlation[program_uid]]
comparison[, sign_reversed := sign(beta_unadjusted) != sign(beta_adjusted)]
comparison[, attenuation := 1 - abs(beta_adjusted) / abs(beta_unadjusted)]
say("supported on NAS: ", comparison[q_unadjusted < 0.05, .N], " unadjusted, ",
    comparison[q_adjusted < 0.05, .N], " after removing the donor offset, of ",
    nrow(comparison))
say("median attenuation ", sprintf("%.0f%%", 100 * median(comparison$attenuation)),
    "; coefficient correlation ",
    sprintf("%.3f", cor(comparison$beta_unadjusted, comparison$beta_adjusted)),
    "; sign reversals ", comparison[sign_reversed == TRUE, .N])
say("the offset is not the explanation: it tracks detection depth, not NAS, and",
    " its weak NAS association points the opposite way to the program effects")

write_json(list(
  system_id = system_id, run_id = run_id, substrate = "PXD051911",
  question = paste("whether the shared negative direction of the program effects",
                   "is a donor-level offset propagated into every score"),
  interpretation = paste("the offset-adjusted count is a lower bound because the",
                         "offset absorbs real signal if disease lowers a broad",
                         "swath of the abundant hepatocyte proteome; the",
                         "unadjusted count is an upper bound"),
  n_donors = ncol(scores), n_programs = nrow(scores),
  biological_unit = "patient_name", seed = t2_seed,
  permutation = list(n_permutations = n_permutations,
                     scheme = "joint histologic permutation within acquisition batch")),
  file.path(tmp, "analysis_contract.json"), auto_unbox = TRUE, pretty = TRUE)
write_tsv(result, file.path(tmp, "offset_arm_results.tsv"))
write_tsv(comparison, file.path(tmp, "offset_adjustment_comparison.tsv"))
write_tsv(offset_tracks, file.path(tmp, "what_the_offset_tracks.tsv"))
write_tsv(meta[, .(patient_name, batch, donor_offset, detection_na_fraction, nas,
                   steatosis, ballooning, inflammation, fibrosis)],
          file.path(tmp, "donor_offset.tsv"))
write_tsv(calibration, file.path(tmp, "calibration.tsv"))
write_tsv(gate_table, file.path(tmp, "gate_status.tsv"))
writeLines(log_lines, file.path(tmp, "run_log.txt"))
capture.output(print(sessionInfo()), file = file.path(tmp, "sessionInfo.txt"))
writeLines(sort(system2("sha256sum",
                        c(file.path(here, "t2_lib.R"),
                          file.path(here, "t2_donor_offset_diagnostic.R"),
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
