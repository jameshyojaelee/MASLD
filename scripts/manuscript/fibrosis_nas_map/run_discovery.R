#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(edgeR)
  library(data.table)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_dir <- dirname(normalizePath(script_arg))
source(file.path(script_dir, "config.R"))
source(file.path(script_dir, "analysis_lib.R"))

set.seed(seed)
assert_true(n_permutations == 10000L && n_bootstrap == 10000L,
            "Locked resampling counts must both equal 10,000")
if (!file.exists(contract_file)) fail("Frozen contract missing: ", contract_file)
if (!file.exists(nonholdout_dge) || !file.exists(nonholdout_meta) ||
    !file.exists(nonholdout_composition)) {
  fail("Blind non-holdout partition missing; run partition_inputs.R first")
}
if (dir.exists(discovery_root)) fail("Refusing to overwrite discovery release: ", discovery_root)

tmp <- atomic_dir(discovery_root)
on.exit(if (dir.exists(tmp)) unlink(tmp, recursive = TRUE), add = TRUE)

message("[1/9] Loading sealed non-holdout input")
dge <- readRDS(nonholdout_dge)
meta <- as.data.table(readRDS(nonholdout_meta))
composition <- fread(nonholdout_composition, check.names = FALSE)
assert_true(ncol(dge) == 1097L, "Non-holdout DGE must contain 1,097 samples")
assert_true(!holdout_cohort %in% as.character(dge$samples$dataset),
            "HOLDOUT LEAK: GSE193066 is present in discovery DGE")
assert_true(!holdout_cohort %in% meta$dataset,
            "HOLDOUT LEAK: GSE193066 is present in discovery metadata")
assert_true(identical(colnames(dge), meta$sample_id), "DGE and metadata sample order mismatch")
assert_true(identical(composition$sample_id, meta$sample_id),
            "Composition and metadata sample order mismatch")
assert_true(all(!is.na(meta$sex_final)), "Discovery metadata contains missing sex_final")
if (!"participant_id" %in% names(meta)) meta[, participant_id := sample_id]
if (!"biopsy" %in% names(meta)) meta[, biopsy := NA_character_]

holdout_participants <- fread(gse193066_crosswalk)[, unique(participant_token)]
assert_true(!any(meta$participant_id %in% holdout_participants),
            "HOLDOUT LEAK: participant identifier overlaps discovery and GSE193066")

sample_manifest <- meta[, .(
  sample_id, participant_id, dataset, biopsy, condition, group_binary, sex_final, age,
  fibrosis_stage, nas_score, diagnosis_harmonized
)]
sample_manifest[, role := fifelse(
  dataset %in% discovery_cohorts & !is.na(fibrosis_stage) & !is.na(nas_score),
  "joint_fibrosis_nas_discovery",
  fifelse(dataset %in% transport_cohorts, "transport", "nonholdout_context")
)]
write_tsv(sample_manifest, file.path(tmp, "sample_manifest.tsv"))

discovery_meta <- meta[
  dataset %in% discovery_cohorts & !is.na(fibrosis_stage) & !is.na(nas_score)
]
observed_n <- discovery_meta[, .N, by = dataset]
for (cohort in names(discovery_expected_n)) {
  n <- observed_n[dataset == cohort, N]
  assert_true(length(n) == 1L && n == discovery_expected_n[[cohort]],
              paste0("Discovery count drift for ", cohort, ": expected ",
                     discovery_expected_n[[cohort]], ", observed ", n))
}
assert_true(nrow(discovery_meta) == 469L, "Joint discovery sample count must be 469")

message("[2/9] Building symbol-level logCPM without outcome-based filtering")
logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
annotation <- fread(gene_annotation, select = c("gene_id", "gene_name"))
symbol_matrix <- collapse_symbols(logcpm, rownames(logcpm), annotation)
rm(logcpm)
gc()

registry <- fread(program_registry)
membership <- fread(program_membership)
assert_true(nrow(registry) == 117L, "Frozen program registry must contain 117 programs")
assert_true(uniqueN(registry$program_uid) == 117L, "Program UIDs must be unique")
assert_true(sum(registry$robust_display) == 2L, "Exactly two programs must be robust_display")

message("[3/9] Scoring frozen programs and fixed comparator views")
program_scores <- score_programs(
  symbol_matrix, meta, membership, registry, program_weight_coverage
)
assert_true(sum(program_scores$coverage$testable) == expected_testable_programs,
            paste0("Expected ", expected_testable_programs, " testable programs; observed ",
                   sum(program_scores$coverage$testable)))
expected_untestable <- c(
  "hotspot_cholangiocytes_fe9b79d294a58bf5",
  "hotspot_hepatocytes_c749bbf208c0ff37",
  "hotspot_hepatocytes_00250b1e5f818798",
  "hotspot_tcells_eaa9c262b58d9323"
)
assert_true(setequal(program_scores$coverage[testable == FALSE, feature_id], expected_untestable),
            "Outcome-blind program testability identities drifted")
discovery_program_coverage <- program_scores$cohort_coverage[
  cohort %in% discovery_cohorts & feature_id %in%
    program_scores$coverage[testable == TRUE, feature_id]
]
assert_true(nrow(discovery_program_coverage) == expected_testable_programs * 4L &&
              all(discovery_program_coverage$cohort_testable),
            "A globally testable program lacks 80% nonconstant weight in a discovery cohort")

hallmark_membership <- read_gmt(hallmark_gmt)
assert_true(uniqueN(hallmark_membership$feature_id) == 50L,
            "Frozen Hallmark GMT must contain 50 gene sets")
hallmark_scores <- score_gene_sets(
  symbol_matrix, meta, hallmark_membership,
  gene_set_min_genes, gene_set_coverage, directional = FALSE
)

published_membership <- read_published_panels(published_panel_root, published_panels)
assert_true(uniqueN(published_membership$feature_id) == 4L,
            "Published comparator family must contain four signatures")
published_scores <- score_gene_sets(
  symbol_matrix, meta, published_membership,
  gene_set_min_genes, gene_set_coverage, directional = TRUE
)
signature_manifest <- merge(
  data.table(
    feature_id = names(published_panels),
    source_path = file.path(published_panel_root, unname(published_panels)),
    source_sha256 = vapply(
      file.path(published_panel_root, unname(published_panels)),
      sha256_file, character(1)
    )
  ),
  published_membership[, .(
    n_members = uniqueN(gene_symbol),
    n_up = uniqueN(gene_symbol[direction == "up"]),
    n_down = uniqueN(gene_symbol[direction == "down"])
  ), by = feature_id],
  by = "feature_id"
)
signature_manifest <- merge(
  signature_manifest, published_scores$coverage,
  by = "feature_id", all.x = TRUE
)
write_tsv(signature_manifest, file.path(tmp, "published_signature_manifest.tsv"))

composition_scores <- build_composition_scores(
  composition, meta, fread(composition_testability), composition_pseudocount
)
assert_true(all(composition_scores$coverage[testable == TRUE, n_observed] ==
                  composition_expected_available - 160L),
            "Non-holdout composition availability drift")

view_specs <- list(
  list(view = "hotspot_program", definition = "weighted_mean_z",
       scores = program_scores$primary, coverage = program_scores$coverage),
  list(view = "hotspot_program", definition = "unweighted_mean_z",
       scores = program_scores$unweighted, coverage = program_scores$coverage),
  list(view = "hotspot_program", definition = "weighted_within_sample_rank",
       scores = program_scores$weighted_rank, coverage = program_scores$coverage),
  list(view = "hallmark", definition = "unweighted_mean_z",
       scores = hallmark_scores$scores, coverage = hallmark_scores$coverage),
  list(view = "published_signature", definition = "directional_mean_z",
       scores = published_scores$scores, coverage = published_scores$coverage),
  list(view = "composition", definition = "clr_1e-6",
       scores = composition_scores$scores, coverage = composition_scores$coverage)
)

composition_status <- data.table(
  view = "composition",
  status = "available",
  reason = "rebuilt_against_corrected_1257_sample_census",
  required_path = composition_candidate
)
write_tsv(composition_status, file.path(tmp, "composition_status.tsv"))

score_long <- rbindlist(lapply(view_specs, function(v) {
  scores_long(v$scores, meta, v$view, v$definition, v$coverage)
}), use.names = TRUE, fill = TRUE)
score_long <- merge(
  score_long,
  registry[, .(feature_id = program_uid, cell_type, module_name,
               robust_display, primary_selected, source_stability)],
  by = "feature_id", all.x = TRUE, sort = FALSE
)
score_long[, cohort := dataset]
score_long <- merge(
  score_long,
  program_scores$cohort_coverage,
  by = c("feature_id", "cohort"), all.x = TRUE, sort = FALSE
)
score_long[, coverage := fcoalesce(cohort_weight_coverage, weight_coverage, gene_coverage)]
score_long[, `:=`(sample = sample_id, participant = participant_id)]
setcolorder(
  score_long,
  c("sample", "participant", "cohort", "biopsy", "view", "feature_id",
    "score_definition", "coverage", "score", "testable")
)
fwrite(score_long, file.path(tmp, "feature_scores.tsv.gz"), sep = "\t", quote = FALSE, na = "NA")

message("[4/9] Fitting cohort models with HC3 standard errors")
cohort_tables <- list()
meta_tables <- list()
table_index <- 0L
for (v in view_specs) {
  eligible <- v$coverage[testable == TRUE, feature_id]
  S <- v$scores[intersect(rownames(v$scores), eligible), , drop = FALSE]
  primary <- fit_feature_models(S, discovery_meta, discovery_cohorts, "primary")
  interaction <- fit_feature_models(S, discovery_meta, discovery_cohorts, "interaction")
  age_meta <- discovery_meta[dataset %in% c("GSE130970", "GSE162694", "GSE174478")]
  age <- fit_feature_models(S, age_meta, unique(age_meta$dataset), "age")
  composition_adjusted <- data.table()
  if (v$view == "hotspot_program" && v$definition == "weighted_mean_z") {
    composition_meta <- merge(
      discovery_meta, composition_scores$logratios,
      by = "sample_id", all.x = TRUE, sort = FALSE
    )
    composition_adjusted <- fit_feature_models(
      S, composition_meta, discovery_cohorts, "composition"
    )
  }
  marginal_f <- fit_feature_models(S, discovery_meta, discovery_cohorts, "marginal_fibrosis")
  disease_meta <- copy(meta[
    dataset %in% discovery_cohorts & group_binary %in% c("Control", "Disease")
  ])
  disease_meta[, group_binary := factor(group_binary, levels = c("Control", "Disease"))]
  disease_control <- fit_feature_models(
    S, disease_meta, unique(disease_meta$dataset), "disease_control"
  )
  nash_meta <- copy(meta[
    dataset %in% discovery_cohorts & diagnosis_harmonized %in% c("NAFL", "NASH")
  ])
  nash_meta[, diagnosis_binary := factor(diagnosis_harmonized, levels = c("NAFL", "NASH"))]
  nash_nafl <- fit_feature_models(S, nash_meta, unique(nash_meta$dataset), "nash_nafl")
  cohort <- rbindlist(list(
    primary, interaction, age, composition_adjusted, marginal_f,
    disease_control, nash_nafl
  ), fill = TRUE)
  cohort[, `:=`(view = v$view, score_definition = v$definition)]
  meta_fit <- rbindlist(list(
    meta_analyze(primary), meta_analyze(interaction),
    meta_analyze(age),
    if (nrow(composition_adjusted)) meta_analyze(composition_adjusted) else data.table(),
    meta_analyze(marginal_f),
    meta_analyze(disease_control), meta_analyze(nash_nafl)
  ), fill = TRUE)
  meta_fit[, `:=`(view = v$view, score_definition = v$definition)]
  table_index <- table_index + 1L
  cohort_tables[[table_index]] <- cohort
  meta_tables[[table_index]] <- meta_fit
}
cohort_effects <- rbindlist(cohort_tables, fill = TRUE)
meta_effects <- rbindlist(meta_tables, fill = TRUE)

meta_effects[, `:=`(
  multiplicity_family = paste(view, score_definition, model_kind, sep = "::"),
  family_size = .N,
  q_value = p.adjust(p_value_meta, method = "BH", n = .N)
), by = .(view, score_definition, model_kind)]

primary_family <- meta_effects[
  view == "hotspot_program" & score_definition == "weighted_mean_z" &
    model_kind == "primary"
]
assert_true(nrow(primary_family) == expected_testable_programs * 2L,
            "Complete primary BH family size drift")
assert_true(all(primary_family$estimable),
            "A testable primary program-axis model is non-estimable")

feature_annotations <- registry[, .(
  feature_id = program_uid, cell_type, module_name, robust_display,
  primary_selected, source_stability
)]
cohort_effects <- merge(
  cohort_effects, feature_annotations, by = "feature_id", all.x = TRUE, sort = FALSE
)
meta_effects <- merge(
  meta_effects, feature_annotations, by = "feature_id", all.x = TRUE, sort = FALSE
)

write_tsv(cohort_effects, file.path(tmp, "cohort_effects.tsv"))
write_tsv(meta_effects, file.path(tmp, "meta_effects.tsv"))

message("[5/9] Building complete program effect map and robustness states")
program_meta <- meta_effects[
  view == "hotspot_program" & score_definition == "weighted_mean_z" &
    model_kind == "primary" & axis %in% c("fibrosis", "nas")
]
wide <- dcast(
  program_meta,
  feature_id ~ axis,
  value.var = c("beta_meta", "se_meta", "ci_lower", "ci_upper", "p_value_meta",
                "q_value", "tau2", "i2", "n_cohorts", "n_same_direction")
)
effect_map <- merge(
  registry[, .(feature_id = program_uid, cell_type, module, module_name,
               robust_display, primary_selected, source_stability)],
  program_scores$coverage[, .(feature_id, weight_coverage, testable)],
  by = "feature_id", all.x = TRUE
)
effect_map <- merge(effect_map, wide, by = "feature_id", all.x = TRUE)

sensitivity <- dcast(
  meta_effects[
    view == "hotspot_program" & model_kind == "primary" & axis %in% c("fibrosis", "nas")
  ],
  feature_id + axis ~ score_definition,
  value.var = "beta_meta"
)
sensitivity[, direction_agree :=
  is.finite(weighted_mean_z) & is.finite(unweighted_mean_z) &
    is.finite(weighted_within_sample_rank) &
    sign(weighted_mean_z) == sign(unweighted_mean_z) &
    sign(weighted_mean_z) == sign(weighted_within_sample_rank)]
sens_wide <- dcast(sensitivity, feature_id ~ axis, value.var = "direction_agree")
setnames(sens_wide, c("fibrosis", "nas"), c("score_direction_agree_fibrosis", "score_direction_agree_nas"))
effect_map <- merge(effect_map, sens_wide, by = "feature_id", all.x = TRUE)
composition_sensitivity <- dcast(
  meta_effects[
    view == "hotspot_program" & score_definition == "weighted_mean_z" &
      model_kind == "composition" & axis %in% c("fibrosis", "nas")
  ],
  feature_id ~ axis, value.var = "beta_meta"
)
setnames(
  composition_sensitivity, c("fibrosis", "nas"),
  c("composition_beta_fibrosis", "composition_beta_nas")
)
effect_map <- merge(effect_map, composition_sensitivity, by = "feature_id", all.x = TRUE)
effect_map[, `:=`(
  composition_direction_agree_fibrosis = is.finite(composition_beta_fibrosis) &
    sign(composition_beta_fibrosis) == sign(beta_meta_fibrosis),
  composition_direction_agree_nas = is.finite(composition_beta_nas) &
    sign(composition_beta_nas) == sign(beta_meta_nas)
)]
effect_map[, evidence_state := fifelse(
  !testable, "untestable",
  fifelse(q_value_fibrosis < 0.05 & q_value_nas < 0.05 & sign(beta_meta_fibrosis) == sign(beta_meta_nas),
          "both_axes_same_direction",
  fifelse(q_value_fibrosis < 0.05 & q_value_nas < 0.05,
          "both_axes_opposing",
  fifelse(q_value_fibrosis < 0.05, "fibrosis_supported_nas_unsupported",
  fifelse(q_value_nas < 0.05, "nas_supported_fibrosis_unsupported", "unsupported")))))
]
effect_map[, robust_fibrosis := testable & q_value_fibrosis < 0.05 &
             score_direction_agree_fibrosis & n_same_direction_fibrosis >= 3]
effect_map[, robust_nas := testable & q_value_nas < 0.05 &
             score_direction_agree_nas & n_same_direction_nas >= 3]
cohort_export <- copy(cohort_effects)
cohort_export[, effect_level := "cohort"]
meta_export <- copy(meta_effects)
meta_export[, `:=`(effect_level = "meta", cohort = "REML_KH")]
effect_map_long <- rbindlist(list(cohort_export, meta_export), fill = TRUE)
effect_map_long[, direction := fifelse(
  effect_level == "cohort",
  fifelse(beta > 0, "positive", fifelse(beta < 0, "negative", NA_character_)),
  fifelse(beta_meta > 0, "positive", fifelse(beta_meta < 0, "negative", NA_character_))
)]
untestable_grid <- CJ(
  feature_id = program_scores$coverage[testable == FALSE, feature_id],
  cohort = c(discovery_cohorts, "REML_KH"),
  axis = c("fibrosis", "nas")
)
untestable_grid[, `:=`(
  view = "hotspot_program", score_definition = "weighted_mean_z",
  model_kind = "primary", effect_level = fifelse(cohort == "REML_KH", "meta", "cohort"),
  estimable = FALSE, failure_reason = "program_weight_coverage_below_0.80"
)]
untestable_grid <- merge(
  untestable_grid, feature_annotations, by = "feature_id", all.x = TRUE, sort = FALSE
)
effect_map_long <- rbindlist(list(effect_map_long, untestable_grid), fill = TRUE)
effect_map_long <- merge(
  effect_map_long,
  effect_map[, .(
    feature_id, weight_coverage, testable, evidence_state,
    robust_fibrosis, robust_nas,
    score_direction_agree_fibrosis, score_direction_agree_nas,
    composition_beta_fibrosis, composition_beta_nas,
    composition_direction_agree_fibrosis, composition_direction_agree_nas
  )],
  by = "feature_id", all.x = TRUE, sort = FALSE
)
setorder(effect_map_long, view, score_definition, model_kind, feature_id, axis,
         effect_level, cohort)
message("[6/9] Auditing ordinal linearity without selecting features")
S_primary <- program_scores$primary
testable_programs <- program_scores$coverage[testable == TRUE, feature_id]
linearity_result <- categorical_marginal_means(
  S_primary[testable_programs, , drop = FALSE],
  discovery_meta, discovery_cohorts
)
linearity <- linearity_result$means
linearity_contrasts <- linearity_result$contrasts
linearity <- merge(
  linearity,
  melt(
    effect_map[, .(feature_id, fibrosis = beta_meta_fibrosis, nas = beta_meta_nas)],
    id.vars = "feature_id", variable.name = "axis", value.name = "linear_beta"
  ),
  by = c("feature_id", "axis"), all.x = TRUE
)
linearity_contrasts <- merge(
  linearity_contrasts,
  melt(
    effect_map[, .(feature_id, fibrosis = beta_meta_fibrosis, nas = beta_meta_nas)],
    id.vars = "feature_id", variable.name = "axis", value.name = "linear_beta"
  ),
  by = c("feature_id", "axis"), all.x = TRUE
)
linearity_contrasts[, `:=`(
  multiplicity_family = "categorical_adjacent_contrasts_all_programs_axes_cohorts",
  family_size = .N,
  q_value = p.adjust(p_value, method = "BH", n = .N)
)]
linearity_contrasts[, `:=`(
  opposite_direction = estimable & is.finite(estimate) & is.finite(linear_beta) &
    sign(estimate) != 0 & sign(linear_beta) != 0 & sign(estimate) != sign(linear_beta),
  supported_reversal = estimable & q_value < 0.05 &
    is.finite(estimate) & is.finite(linear_beta) & sign(estimate) != 0 &
    sign(linear_beta) != 0 & sign(estimate) != sign(linear_beta)
)]

reversal_cohort <- linearity_contrasts[, .(
  linearity_reversal = any(supported_reversal),
  n_supported_reversal_contrasts = sum(supported_reversal)
), by = .(feature_id, cohort, axis)]
linearity <- merge(
  linearity, reversal_cohort,
  by = c("feature_id", "cohort", "axis"), all.x = TRUE, sort = FALSE
)
linearity[is.na(linearity_reversal), `:=`(
  linearity_reversal = FALSE,
  n_supported_reversal_contrasts = 0L
)]
write_tsv(linearity, file.path(tmp, "linearity_audit.tsv"))
write_tsv(linearity_contrasts, file.path(tmp, "linearity_contrasts.tsv"))

linearity_cohort <- unique(linearity[, .(
  feature_id, cohort, axis, estimable, linearity_reversal
)])
linearity_summary <- linearity_cohort[, .(
  linearity_reversal = any(linearity_reversal),
  linearity_estimable_cohorts = sum(estimable)
), by = .(feature_id, axis)]
reversal <- dcast(
  linearity_summary,
  feature_id ~ axis, value.var = "linearity_reversal"
)
setnames(reversal, c("fibrosis", "nas"),
         c("linearity_reversal_fibrosis", "linearity_reversal_nas"))
linearity_n <- dcast(
  linearity_summary,
  feature_id ~ axis, value.var = "linearity_estimable_cohorts"
)
setnames(linearity_n, c("fibrosis", "nas"),
         c("linearity_estimable_cohorts_fibrosis", "linearity_estimable_cohorts_nas"))
reversal <- merge(reversal, linearity_n, by = "feature_id")
effect_map <- merge(effect_map, reversal, by = "feature_id", all.x = TRUE)
effect_map[, `:=`(
  monotonic_summary_allowed_fibrosis = !linearity_reversal_fibrosis &
    linearity_estimable_cohorts_fibrosis == 4L,
  monotonic_summary_allowed_nas = !linearity_reversal_nas &
    linearity_estimable_cohorts_nas == 4L
)]
effect_map_long <- merge(effect_map_long, reversal, by = "feature_id", all.x = TRUE)
effect_map_long <- merge(
  effect_map_long,
  effect_map[, .(
    feature_id, monotonic_summary_allowed_fibrosis,
    monotonic_summary_allowed_nas
  )],
  by = "feature_id", all.x = TRUE, sort = FALSE
)
setorder(effect_map_long, view, score_definition, model_kind, feature_id, axis,
         effect_level, cohort)
write_tsv(effect_map, file.path(tmp, "program_map.tsv"))
write_tsv(effect_map_long, file.path(tmp, "effect_map.tsv"))

message("[7/9] Running four-fold leave-one-cohort-out validation")
primary_cohort <- cohort_effects[
  view == "hotspot_program" & score_definition == "weighted_mean_z" &
    model_kind == "primary" & axis %in% c("fibrosis", "nas") &
    feature_id %in% testable_programs & estimable == TRUE
]
fold_rows <- list()
permutation_rho <- array(
  NA_real_, dim = c(n_permutations, length(discovery_cohorts), 2L),
  dimnames = list(NULL, discovery_cohorts, c("fibrosis", "nas"))
)

for (fold_i in seq_along(discovery_cohorts)) {
  held <- discovery_cohorts[[fold_i]]
  training_effects <- primary_cohort[cohort != held]
  training_meta <- meta_analyze(training_effects, min_cohorts = 3L)
  held_effects <- primary_cohort[cohort == held]
  held_wide <- dcast(held_effects, feature_id ~ axis, value.var = "beta")
  train_wide <- dcast(training_meta, feature_id ~ axis, value.var = "beta_meta")
  setnames(train_wide, c("fibrosis", "nas"), c("train_fibrosis", "train_nas"))
  setnames(held_wide, c("fibrosis", "nas"), c("held_fibrosis", "held_nas"))
  cmp <- merge(train_wide, held_wide, by = "feature_id")
  assert_true(nrow(cmp) == expected_testable_programs,
              paste0("LOCO feature census drift in ", held))
  for (axis in c("fibrosis", "nas")) {
    rho <- cor(cmp[[paste0("train_", axis)]], cmp[[paste0("held_", axis)]],
               method = "spearman", use = "complete.obs")
    fold_rows[[length(fold_rows) + 1L]] <- data.table(
      split = "leave_one_cohort_out", cohort = held, view = "hotspot_program",
      axis = axis, metric = "effect_vector_spearman", estimate = rho,
      n_features = nrow(cmp), biological_n = discovery_expected_n[[held]]
    )
  }

  held_meta <- copy(discovery_meta[dataset == held])
  held_scores <- S_primary[testable_programs, held_meta$sample_id, drop = FALSE]
  for (b in seq_len(n_permutations)) {
    perm <- sample.int(nrow(held_meta))
    perm_meta <- copy(held_meta)
    perm_meta[, `:=`(
      fibrosis_stage = held_meta$fibrosis_stage[perm],
      nas_score = held_meta$nas_score[perm]
    )]
    B <- fit_matrix_coefficients(held_scores, perm_meta, "primary")
    if (is.null(B)) next
    permutation_rho[b, fold_i, "fibrosis"] <- cor(
      cmp$train_fibrosis,
      B["fibrosis_centered", match(cmp$feature_id, colnames(B))],
      method = "spearman", use = "complete.obs"
    )
    permutation_rho[b, fold_i, "nas"] <- cor(
      cmp$train_nas,
      B["nas_centered", match(cmp$feature_id, colnames(B))],
      method = "spearman", use = "complete.obs"
    )
  }
}

fold_dt <- rbindlist(fold_rows)
global_rows <- list()
for (axis_name in c("fibrosis", "nas")) {
  observed <- mean(atanh(pmax(pmin(
    fold_dt[axis == axis_name, estimate], 0.999999
  ), -0.999999)))
  null <- apply(permutation_rho[, , axis_name, drop = FALSE], 1L, function(x) {
    mean(atanh(pmax(pmin(x, 0.999999), -0.999999)), na.rm = TRUE)
  })
  global_rows[[axis_name]] <- data.table(
    split = "leave_one_cohort_out_global", cohort = "four_discovery_cohorts",
    view = "hotspot_program", axis = axis_name,
    metric = "mean_fisher_z_effect_vector_spearman",
    estimate = observed,
    empirical_p = empirical_p(observed, null, "greater"),
    positive_folds = sum(fold_dt[axis == axis_name, estimate] > 0),
    n_features = length(testable_programs), biological_n = nrow(discovery_meta)
  )
}
global_dt <- rbindlist(global_rows)
assert_true(all(is.finite(permutation_rho)), "LOCO permutation null contains non-finite values")
global_dt[, p_holm := p.adjust(empirical_p, method = "holm")]
validation <- rbindlist(list(fold_dt, global_dt), fill = TRUE)
write_tsv(validation, file.path(tmp, "validation_summary.tsv"))
saveRDS(permutation_rho, file.path(tmp, "loco_permutation_null.rds"), compress = "xz")

message("[8/9] Writing input, environment, and artifact manifests")
write_tsv(data.table(
  input = c("contract", "nonholdout_dge", "nonholdout_meta", "program_registry",
            "program_membership", "nonholdout_composition", "composition_testability",
            "hallmark_gmt"),
  path = c(contract_file, nonholdout_dge, nonholdout_meta, program_registry,
           program_membership, nonholdout_composition, composition_testability,
           hallmark_gmt),
  sha256 = vapply(c(contract_file, nonholdout_dge, nonholdout_meta, program_registry,
                    program_membership, nonholdout_composition, composition_testability,
                    hallmark_gmt), sha256_file, character(1))
), file.path(tmp, "input_manifest.tsv"))

session <- capture.output(sessionInfo())
writeLines(session, file.path(tmp, "sessionInfo.txt"))
write_tsv(data.table(
  check = c("holdout_absent", "discovery_n", "program_registry_n", "program_testable_n",
            "robust_display_n", "hallmark_n", "published_signature_n",
            "composition_declared_n", "composition_testable_n", "participant_overlap_n"),
  value = c(TRUE, nrow(discovery_meta), nrow(registry), sum(program_scores$coverage$testable),
            sum(registry$robust_display), uniqueN(hallmark_membership$feature_id),
            uniqueN(published_membership$feature_id), nrow(composition_scores$coverage),
            sum(composition_scores$coverage$testable),
            sum(meta$participant_id %in% holdout_participants)),
  expected = c(TRUE, 469, 117, expected_testable_programs, 2, 50, 4, 22, 16, 0)
), file.path(tmp, "validation_checks.tsv"))

artifacts <- list.files(tmp, recursive = TRUE, full.names = TRUE)
artifacts <- artifacts[file.info(artifacts)$isdir == FALSE]
artifact_manifest <- data.table(
  artifact = substring(artifacts, nchar(tmp) + 2L),
  sha256 = vapply(artifacts, sha256_file, character(1)),
  size_bytes = file.info(artifacts)$size
)
write_tsv(artifact_manifest, file.path(tmp, "artifact_manifest.tsv"))
ready_payload <- list(
  state = "DISCOVERY_COMPLETE_HOLDOUT_UNOPENED",
  release_id = release_id,
  workstream_id = workstream_id,
  contract_sha256 = sha256_file(contract_file),
  n_discovery = nrow(discovery_meta),
  n_programs = nrow(registry),
  n_testable_programs = sum(program_scores$coverage$testable),
  n_permutations = n_permutations,
  holdout_accessed = FALSE
)
jsonlite::write_json(ready_payload, file.path(tmp, "DISCOVERY_READY.json"),
                     pretty = TRUE, auto_unbox = TRUE)

message("[9/9] Publishing immutable discovery release")
publish_dir(tmp, discovery_root)
Sys.chmod(list.files(discovery_root, recursive = TRUE, full.names = TRUE), mode = "0440")
message("Published discovery: ", discovery_root)
