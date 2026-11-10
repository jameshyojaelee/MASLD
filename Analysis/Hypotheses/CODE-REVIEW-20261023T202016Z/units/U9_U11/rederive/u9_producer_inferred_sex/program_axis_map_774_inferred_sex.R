#!/usr/bin/env Rscript

# Figure 4E program-axis models on the adopted 774-unit bulk roster.
#
# Reuses, unchanged, the scoring, CLR-component, HC3 and meta-analysis functions
# of scripts/manuscript/program_observability_map/analysis_lib.R and the frozen
# constants in config.R. No program is rediscovered, reweighted or filtered by
# outcome. Four arms are fitted with the same code path:
#
#   old469_reproduction      stored BULK-PROGRAM-MAP-v9 inputs, four cohorts
#                            (positive control; must reproduce the stored map)
#   old_three_cohorts        stored inputs without GSE174478 (isolates the loss
#                            of the fourth cohort from the change of roster)
#   new774_three_cohorts     PRIMARY: adopted 774-unit DGE and rebuilt
#                            composition; cohorts carrying both fibrosis and NAS
#   new774_plus_GSE174478    primary arm plus the unchanged GSE174478 records
#                            (outside the 774 units; old four-cohort design)

suppressPackageStartupMessages({
  library(edgeR)
  library(data.table)
})

repo <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
Sys.setenv(MASLD_PROJECT_ROOT = repo)
lib_dir <- file.path(repo, "scripts/manuscript/program_observability_map")
source(file.path(lib_dir, "config.R"))
source(file.path(lib_dir, "analysis_lib.R"))

out_root <- Sys.getenv("FIG4E_774_OUT")
if (!nzchar(out_root)) fail("FIG4E_774_OUT is required")
analysis_root <- file.path(out_root, "analysis")
if (dir.exists(analysis_root)) fail("Refusing to overwrite: ", analysis_root)
dir.create(analysis_root, recursive = TRUE)

# No step below draws random numbers (OLS, HC3, prcomp, REML). The seed is the
# frozen config seed, set so that any later stochastic addition is reproducible.
set.seed(seed)

new_root <- file.path(
  repo, "RNA-seq/results/manuscript_release/candidates",
  "bulk-corrected-controls-20261008T162308Z/primary_source_controls"
)
new_dge_path <- file.path(new_root, "merged_dge.rds")
new_meta_path <- file.path(new_root, "meta_matched.rds")
new_composition_path <- file.path(
  repo, "figures/candidates/composition-774-20261011T151616Z/composition",
  "accepted_composition_774_plus_unchanged4.tsv.gz"
)
unit_manifest_path <- file.path(
  repo, "figures/candidates/bulk-release-corrected-controls-20261008T162308Z",
  "analysis/stage_extensions/five_cohort_sample_manifest.tsv"
)
stored_map_path <- file.path(discovery_root, "program_axis_map.tsv")
stored_meta_effects_path <- file.path(discovery_root, "meta_effects.tsv")

stromal_ecm <- "hotspot_hepatocytes_f05c535ae5bbc0b9"
ductular <- "hotspot_hepatocytes_48f39dd4d817a10e"
three_cohorts <- c("GSE130970", "GSE135251", "GSE162694")
extra_cohort <- "GSE174478"

registry <- fread(program_registry)
membership <- fread(program_membership)
assert_true(nrow(registry) == 117L, "Frozen program registry must contain 117 programs")
assert_true(all(c(stromal_ecm, ductular) %in% registry$program_uid),
            "Named programs missing from the registry")
annotation <- fread(gene_annotation, select = c("gene_id", "gene_name"))

score_substrate <- function(dge, meta, composition, label) {
  message("[", label, "] samples ", ncol(dge), ", genes ", nrow(dge))
  assert_true(setequal(colnames(dge), meta$sample_id), "DGE and metadata rosters differ")
  meta <- meta[match(colnames(dge), sample_id)]
  logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
  symbol_matrix <- collapse_symbols(logcpm, rownames(logcpm), annotation)
  rm(logcpm)
  scores <- score_programs(symbol_matrix, meta, membership, registry, program_weight_coverage)
  composition <- composition[match(meta$sample_id, sample_id)]
  assert_true(identical(composition$sample_id, meta$sample_id),
              "Composition does not cover every record")
  lineages <- setdiff(names(composition)[vapply(composition, is.numeric, logical(1))],
                      c("sample_id", "dataset"))
  assert_true(length(lineages) == 16L, "Composition must carry 16 lineages")
  fractions <- as.matrix(composition[, ..lineages])
  zero_audit <- data.table(
    substrate = label, n_records = nrow(fractions), n_lineages = ncol(fractions),
    n_exact_zero = sum(fractions == 0), n_below_pseudocount = sum(fractions < composition_pseudocount),
    min_fraction = min(fractions), pseudocount = composition_pseudocount,
    max_abs_rowsum_minus_1 = max(abs(rowSums(fractions) - 1))
  )
  comp <- build_composition_scores(
    composition, meta, lineages, composition_pseudocount,
    eligible_lineages = character(0), n_components = composition_n_clr_components
  )
  list(meta = meta, primary = scores$primary, coverage = scores$coverage,
       cohort_coverage = scores$cohort_coverage, pcs = comp$pcs,
       pc_variance = comp$pc_variance, zero_audit = zero_audit)
}

fit_arm <- function(arm, S, meta, pcs, coverage, cohorts, unit_col) {
  message("[", arm, "] fitting")
  common <- c("sample_id", "dataset", "sex_final", "fibrosis_stage", "nas_score", unit_col)
  discovery <- meta[dataset %in% cohorts & !is.na(fibrosis_stage) & !is.na(nas_score), ..common]
  assert_true(all(!is.na(discovery$sex_final)), "Missing sex_final in the model roster")
  assert_true(uniqueN(discovery[[unit_col]]) == nrow(discovery), "More than one record per unit")
  composition_meta <- merge(discovery, pcs, by = "sample_id", all.x = TRUE, sort = FALSE)
  setkey(composition_meta, NULL)
  testable <- coverage[testable == TRUE, feature_id]
  S <- S[testable, , drop = FALSE]
  linear <- fit_feature_models(S, discovery, cohorts, "primary")
  adjusted <- fit_feature_models(S, composition_meta, cohorts, "composition_pc")
  cohort_effects <- rbindlist(list(linear, adjusted), use.names = TRUE, fill = TRUE)
  meta_effects <- rbindlist(list(meta_analyze(linear), meta_analyze(adjusted)),
                            use.names = TRUE, fill = TRUE)
  # Same multiplicity rule as run_discovery.R: one BH family per model, pooling
  # the fibrosis and NAS terms of every eligible program.
  meta_effects[, `:=`(
    multiplicity_family = paste("hotspot_program::weighted_mean_z", model_kind, sep = "::"),
    family_size = .N,
    q_value = p.adjust(p_value_meta, method = "BH", n = .N)
  ), by = model_kind]
  value_cols <- c("beta_meta", "se_meta", "ci_lower", "ci_upper", "p_value_meta",
                  "q_value", "tau2", "n_cohorts", "n_same_direction")
  linear_wide <- dcast(meta_effects[model_kind == "primary"], feature_id ~ axis,
                       value.var = value_cols)
  comp_wide <- dcast(meta_effects[model_kind == "composition_pc"], feature_id ~ axis,
                     value.var = value_cols)
  setnames(comp_wide, setdiff(names(comp_wide), "feature_id"),
           paste0("composition_pc_", setdiff(names(comp_wide), "feature_id")))
  axis_map <- Reduce(function(x, y) merge(x, y, by = "feature_id", all.x = TRUE),
                     list(registry[, .(feature_id = program_uid, cell_type, module, module_name)],
                          coverage[, .(feature_id, weight_coverage, testable)],
                          linear_wide, comp_wide))
  roster <- discovery[, .(n_units = .N, n_female = sum(sex_final == "F"),
                          n_male = sum(sex_final == "M"),
                          fibrosis_levels = paste(sort(unique(fibrosis_stage)), collapse = ","),
                          nas_levels = paste(sort(unique(nas_score)), collapse = ",")),
                      by = dataset]
  arm_dir <- file.path(analysis_root, arm)
  dir.create(arm_dir)
  write_tsv(cohort_effects, file.path(arm_dir, "cohort_effects.tsv"))
  write_tsv(meta_effects, file.path(arm_dir, "meta_effects.tsv"))
  write_tsv(axis_map, file.path(arm_dir, "program_axis_map.tsv"))
  write_tsv(roster, file.path(arm_dir, "model_roster_by_cohort.tsv"))
  write_tsv(discovery, file.path(arm_dir, "model_roster.tsv"))
  list(axis_map = axis_map, roster = roster, n = nrow(discovery))
}

# ---- old substrate -------------------------------------------------------
old_meta <- as.data.table(readRDS(nonholdout_meta))
old_meta[, unit := sample_id]
old <- score_substrate(readRDS(nonholdout_dge), old_meta,
                       fread(nonholdout_composition, check.names = FALSE), "old_1097_records")
gc()

# ---- new substrate -------------------------------------------------------
new_meta <- as.data.table(readRDS(new_meta_path))
new_meta[, sex_final := inferred_sex]  # CODE-REVIEW: Methods covariate
manifest <- fread(unit_manifest_path)
assert_true(nrow(new_meta) == 774L && uniqueN(new_meta$analysis_unit_id) == 774L,
            "The adopted roster must be 774 units, one record each")
assert_true(setequal(manifest$sample_id, new_meta$sample_id), "Unit manifest roster differs")
mm <- manifest[match(new_meta$sample_id, sample_id)]
assert_true(identical(mm$fibrosis_stage, new_meta$fibrosis_stage) &&
              identical(mm$nas_score, new_meta$nas_score),
            "Manifest stage or NAS differs from the model metadata")
new_meta[, unit := analysis_unit_id]
new <- score_substrate(readRDS(new_dge_path), new_meta,
                       fread(new_composition_path, check.names = FALSE), "new_774_units")
gc()

write_tsv(rbind(old$zero_audit, new$zero_audit), file.path(analysis_root, "clr_zero_audit.tsv"))
write_tsv(rbind(cbind(substrate = "old_1097_records", old$pc_variance),
                cbind(substrate = "new_774_units", new$pc_variance)),
          file.path(analysis_root, "composition_pc_variance.tsv"))
write_tsv(merge(old$coverage[, .(feature_id, old_weight_coverage = weight_coverage,
                                 old_testable = testable)],
                new$coverage[, .(feature_id, new_weight_coverage = weight_coverage,
                                 new_testable = testable)], by = "feature_id"),
          file.path(analysis_root, "program_weight_coverage_old_vs_new.tsv"))

# ---- arms ----------------------------------------------------------------
arms <- list()
arms$old469_reproduction <- fit_arm("old469_reproduction", old$primary, old$meta, old$pcs,
                                    old$coverage, discovery_cohorts, "unit")
arms$old_three_cohorts <- fit_arm("old_three_cohorts", old$primary, old$meta, old$pcs,
                                  old$coverage, three_cohorts, "unit")
arms$new774_three_cohorts <- fit_arm("new774_three_cohorts", new$primary, new$meta, new$pcs,
                                     new$coverage, three_cohorts, "unit")

# Primary scores for the three 774-unit cohorts plus the unchanged GSE174478
# records scored in their own (old) expression matrix. Eligibility is the
# 774-unit weight coverage; a program not scored in GSE174478 is fitted on the
# remaining three cohorts.
extra_ids <- old$meta[dataset == extra_cohort, sample_id]
combined_scores <- cbind(
  new$primary,
  old$primary[rownames(new$primary), extra_ids, drop = FALSE]
)
keep <- c("sample_id", "dataset", "sex_final", "fibrosis_stage", "nas_score", "unit")
combined_meta <- rbind(new$meta[, ..keep], old$meta[dataset == extra_cohort, ..keep])
combined_pcs <- rbind(new$pcs, old$pcs[sample_id %in% extra_ids])
assert_true(identical(colnames(combined_scores), combined_meta$sample_id),
            "Combined score and metadata order mismatch")
arms$new774_plus_GSE174478 <- fit_arm(
  "new774_plus_GSE174478", combined_scores, combined_meta, combined_pcs,
  new$coverage, c(three_cohorts, extra_cohort), "unit"
)

# ---- reproduction of the stored map --------------------------------------
stored <- fread(stored_map_path)
repro <- merge(stored, arms$old469_reproduction$axis_map, by = "feature_id",
               suffixes = c("_stored", "_rerun"))
check_cols <- c("beta_meta_fibrosis", "beta_meta_nas", "q_value_fibrosis", "q_value_nas",
                "composition_pc_beta_meta_fibrosis", "composition_pc_beta_meta_nas",
                "composition_pc_q_value_fibrosis", "composition_pc_q_value_nas")
reproduction <- rbindlist(lapply(check_cols, function(col) {
  a <- repro[[paste0(col, "_stored")]]; b <- repro[[paste0(col, "_rerun")]]
  data.table(quantity = col, n_compared = sum(is.finite(a) & is.finite(b)),
             n_na_mismatch = sum(is.finite(a) != is.finite(b)),
             max_abs_difference = max(abs(a - b), na.rm = TRUE))
}))
reproduction[, reproduced := n_na_mismatch == 0L & max_abs_difference < 1e-8]
write_tsv(reproduction, file.path(analysis_root, "stored_map_reproduction.tsv"))
assert_true(all(reproduction$reproduced) &&
              identical(sort(stored[testable == TRUE, feature_id]),
                        sort(arms$old469_reproduction$axis_map[testable == TRUE, feature_id])),
            "The stored 469-record map was not reproduced")

# ---- summaries -----------------------------------------------------------
count_calls <- function(map, prefix, axis) {
  q <- map[[paste0(prefix, "q_value_", axis)]]
  beta <- map[[paste0(prefix, "beta_meta_", axis)]]
  list(n_estimable = sum(is.finite(q)), n_q05 = sum(q < 0.05, na.rm = TRUE),
       n_q05_positive = sum(q < 0.05 & beta > 0, na.rm = TRUE),
       n_q05_negative = sum(q < 0.05 & beta < 0, na.rm = TRUE),
       n_nominal_p05 = sum(map[[paste0(prefix, "p_value_meta_", axis)]] < 0.05, na.rm = TRUE),
       min_q = min(q, na.rm = TRUE),
       median_abs_beta = median(abs(beta), na.rm = TRUE),
       median_se = median(map[[paste0(prefix, "se_meta_", axis)]], na.rm = TRUE))
}
summary_rows <- rbindlist(lapply(names(arms), function(arm) {
  map <- arms[[arm]]$axis_map
  rbindlist(lapply(c("", "composition_pc_"), function(prefix) {
    rbindlist(lapply(c("fibrosis", "nas"), function(axis) {
      c(list(arm = arm, n_units = arms[[arm]]$n, n_cohorts = nrow(arms[[arm]]$roster),
             model = if (prefix == "") "unadjusted" else "composition_adjusted",
             axis = axis, n_registry = nrow(map), n_eligible = sum(map$testable),
             bh_family_size = 2L * sum(map$testable)),
        count_calls(map, prefix, axis))
    }))
  }))
}))
write_tsv(summary_rows, file.path(analysis_root, "arm_summary.tsv"))
write_tsv(rbindlist(lapply(names(arms), function(arm) cbind(arm = arm, arms[[arm]]$roster))),
          file.path(analysis_root, "arm_rosters.tsv"))

compare <- function(reference, target) {
  a <- arms[[reference]]$axis_map; b <- arms[[target]]$axis_map
  m <- merge(a, b, by = "feature_id", suffixes = c("_ref", "_new"))
  rbindlist(lapply(c("", "composition_pc_"), function(prefix) {
    rbindlist(lapply(c("fibrosis", "nas"), function(axis) {
      x <- m[[paste0(prefix, "beta_meta_", axis, "_ref")]]
      y <- m[[paste0(prefix, "beta_meta_", axis, "_new")]]
      qx <- m[[paste0(prefix, "q_value_", axis, "_ref")]]
      qy <- m[[paste0(prefix, "q_value_", axis, "_new")]]
      ok <- is.finite(x) & is.finite(y)
      data.table(
        reference = reference, target = target,
        model = if (prefix == "") "unadjusted" else "composition_adjusted", axis = axis,
        n_programs_both = sum(ok),
        n_same_sign = sum(sign(x[ok]) == sign(y[ok])),
        sign_agreement = mean(sign(x[ok]) == sign(y[ok])),
        pearson_r = stats::cor(x[ok], y[ok]),
        spearman_rho = stats::cor(x[ok], y[ok], method = "spearman"),
        slope_new_on_ref = unname(stats::coef(stats::lm(y[ok] ~ x[ok]))[2]),
        n_q05_ref = sum(qx[ok] < 0.05), n_q05_new = sum(qy[ok] < 0.05),
        n_q05_both = sum(qx[ok] < 0.05 & qy[ok] < 0.05),
        n_q05_ref_same_sign_in_new = sum(qx[ok] < 0.05 & sign(x[ok]) == sign(y[ok]))
      )
    }))
  }))
}
comparisons <- rbindlist(list(
  compare("old469_reproduction", "new774_three_cohorts"),
  compare("old469_reproduction", "new774_plus_GSE174478"),
  compare("old469_reproduction", "old_three_cohorts"),
  compare("old_three_cohorts", "new774_three_cohorts")
))
write_tsv(comparisons, file.path(analysis_root, "old_vs_new_agreement.tsv"))

named <- rbindlist(lapply(names(arms), function(arm) {
  me <- fread(file.path(analysis_root, arm, "meta_effects.tsv"))
  me[feature_id %in% c(stromal_ecm, ductular),
     .(arm = arm, program = fifelse(feature_id == stromal_ecm, "Stromal ECM", "Ductular injury"),
       feature_id, model = fifelse(model_kind == "primary", "unadjusted", "composition_adjusted"),
       axis, beta_meta, ci_lower, ci_upper, p_value_meta, q_value, family_size,
       n_cohorts, n_same_direction, tau2, i2)]
}))
setorder(named, program, axis, model, arm)
write_tsv(named, file.path(analysis_root, "stromal_ecm_ductular_by_arm.tsv"))
named_cohort <- rbindlist(lapply(names(arms), function(arm) {
  ce <- fread(file.path(analysis_root, arm, "cohort_effects.tsv"))
  ce[feature_id %in% c(stromal_ecm, ductular),
     .(arm = arm, program = fifelse(feature_id == stromal_ecm, "Stromal ECM", "Ductular injury"),
       cohort, model_kind, axis, beta, se_hc3, ci_lower, ci_upper, p_value, n)]
}))
write_tsv(named_cohort, file.path(analysis_root, "stromal_ecm_ductular_cohort_effects.tsv"))

per_program <- Reduce(function(x, y) merge(x, y, by = "feature_id", all = TRUE), lapply(
  names(arms), function(arm) {
    cols <- c("feature_id", grep("^(composition_pc_)?(beta_meta|ci_lower|ci_upper|p_value_meta|q_value|n_cohorts)_",
                                 names(arms[[arm]]$axis_map), value = TRUE), "testable")
    d <- arms[[arm]]$axis_map[, ..cols]
    setnames(d, setdiff(names(d), "feature_id"), paste0(arm, "::", setdiff(names(d), "feature_id")))
    d
  }))
per_program <- merge(registry[, .(feature_id = program_uid, cell_type, module, module_name)],
                     per_program, by = "feature_id")
write_tsv(per_program, file.path(analysis_root, "per_program_all_arms.tsv"))

inputs <- c(old_dge = nonholdout_dge, old_meta = nonholdout_meta,
            old_composition = nonholdout_composition, new_dge = new_dge_path,
            new_meta = new_meta_path, new_composition = new_composition_path,
            unit_manifest = unit_manifest_path, program_registry = program_registry,
            program_membership = program_membership, gene_annotation = gene_annotation,
            stored_program_axis_map = stored_map_path,
            analysis_lib = file.path(lib_dir, "analysis_lib.R"),
            config = file.path(lib_dir, "config.R"))
write_tsv(data.table(input = names(inputs), path = inputs,
                     sha256 = vapply(inputs, sha256_file, character(1))),
          file.path(analysis_root, "input_manifest.tsv"))
writeLines(c(paste("seed:", seed, "(no stochastic step is executed)"),
             capture.output(sessionInfo())), file.path(analysis_root, "sessionInfo.txt"))
print(summary_rows)
print(comparisons)
print(named)
print(reproduction)
cat("FIG4E_774_ANALYSIS_COMPLETE\n")
