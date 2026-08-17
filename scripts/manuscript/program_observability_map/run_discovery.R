#!/usr/bin/env Rscript

# BULK-PROGRAM-MAP-v9 discovery.
#
# The frozen 117-program projection is byte-identical to v8, so the linear
# fibrosis and NAS coefficients reproduce the v8 estimand exactly and the two
# releases are directly comparable. Four things are new.
#
#  1. OMNIBUS AXIS TESTS. v8 modelled both axes as a single linear term. Its own
#     linearity audit then flagged a BH-supported adjacent-level reversal on the
#     NAS axis for 35 of 113 programs against 1 of 113 on the fibrosis axis, and
#     the renderer ignored the resulting monotonic_summary_allowed flags. A
#     linear term has little power against a non-monotone response, so the v8
#     headline of broader fibrosis than NAS remodeling may be a specification
#     artifact rather than biology. v9 therefore reports three nested views per
#     axis: the v8 linear term, a 2-df quadratic shape test that is
#     meta-analysable, and a shape-agnostic categorical omnibus combined across
#     cohorts. Whether the asymmetry survives is the question, not the premise.
#  2. COMPOSITION. The sealed acceptance decision replaces the degenerate MuSiC
#     estimates with BayesPrism and, because only macrophages clears the
#     dynamic-range floor, replaces the four-log-ratio block with CLR principal
#     components. This is the first time the composition question is actually
#     estimable on all four discovery cohorts.
#  3. OBSERVABILITY. Every program-axis carries a minimum detectable effect, so
#     an unsupported program can be separated into an informative null and one
#     the design could never have resolved.
#  4. SEPARABILITY. A canonical analysis asks whether fibrosis and NAS are
#     separable dimensions in program space at all, which tests whether "no
#     program supported on both axes" is an estimand artifact.
#
# The GSE193066 paired holdout is not opened here and is not opened by v9 at all.

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
if (!file.exists(composition_accepted_seal)) fail("Composition acceptance seal missing")
if (!file.exists(nonholdout_dge) || !file.exists(nonholdout_meta) ||
    !file.exists(nonholdout_composition)) {
  fail("Blind non-holdout partition missing; run partition_inputs.R first")
}
if (dir.exists(discovery_root)) fail("Refusing to overwrite discovery release: ", discovery_root)

tmp <- atomic_dir(discovery_root)
on.exit(if (dir.exists(tmp)) unlink(tmp, recursive = TRUE), add = TRUE)

message("[1/10] Loading sealed non-holdout input")
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

holdout_participants <- fread(gse193066_crosswalk)[, unique(participant_token)]
assert_true(!any(meta$participant_id %in% holdout_participants),
            "HOLDOUT LEAK: participant identifier overlaps discovery and GSE193066")

discovery_meta <- meta[
  dataset %in% discovery_cohorts & !is.na(fibrosis_stage) & !is.na(nas_score)
]
for (cohort in names(discovery_expected_n)) {
  n <- discovery_meta[dataset == cohort, .N]
  assert_true(n == discovery_expected_n[[cohort]],
              paste0("Discovery count drift for ", cohort, ": expected ",
                     discovery_expected_n[[cohort]], ", observed ", n))
}
assert_true(nrow(discovery_meta) == 469L, "Joint discovery sample count must be 469")
assert_true(uniqueN(discovery_meta$participant_id) == 469L,
            "Discovery samples must be one per participant")

message("[2/10] Building symbol-level logCPM without outcome-based filtering")
logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
annotation <- fread(gene_annotation, select = c("gene_id", "gene_name"))
symbol_matrix <- collapse_symbols(logcpm, rownames(logcpm), annotation)
rm(logcpm); gc()

registry <- fread(program_registry)
membership <- fread(program_membership)
assert_true(nrow(registry) == 117L, "Frozen program registry must contain 117 programs")

message("[3/10] Scoring the frozen programs (projection identical to v8)")
program_scores <- score_programs(
  symbol_matrix, meta, membership, registry, program_weight_coverage
)
assert_true(sum(program_scores$coverage$testable) == expected_testable_programs,
            paste0("Expected ", expected_testable_programs, " testable programs; observed ",
                   sum(program_scores$coverage$testable)))
testable_programs <- program_scores$coverage[testable == TRUE, feature_id]
S <- program_scores$primary[testable_programs, , drop = FALSE]

message("[4/10] Building the accepted composition covariates")
seal <- jsonlite::read_json(composition_accepted_seal, simplifyVector = TRUE)
eligible <- if (length(seal$eligible_covariate_lineages))
  as.character(seal$eligible_covariate_lineages) else character(0)
lineages <- setdiff(names(composition)[vapply(composition, is.numeric, logical(1))],
                    c("sample_id", "dataset"))
composition_scores <- build_composition_scores(
  composition, meta, lineages, composition_pseudocount,
  eligible_lineages = eligible, n_components = composition_n_clr_components
)
message("      accepted estimator: ", seal$accepted_method,
        " | eligible log-ratio lineages: ",
        if (length(eligible)) paste(eligible, collapse = ", ") else "none",
        " | CLR components: ", length(composition_scores$pc_columns))
write_tsv(composition_scores$pc_variance, file.path(tmp, "composition_pc_variance.tsv"))

composition_meta <- Reduce(
  function(x, y) merge(x, y, by = "sample_id", all.x = TRUE, sort = FALSE),
  list(discovery_meta, composition_scores$pcs, composition_scores$logratios)
)
setkey(composition_meta, NULL)

message("[5/10] Fitting linear, shape and omnibus models with HC3")
linear <- fit_feature_models(S, discovery_meta, discovery_cohorts, "primary")
marginal_f <- fit_feature_models(S, discovery_meta, discovery_cohorts, "marginal_fibrosis")
marginal_n <- fit_feature_models(S, discovery_meta, discovery_cohorts, "marginal_nas")
composition_pc <- fit_feature_models(S, composition_meta, discovery_cohorts, "composition_pc")
composition_mac <- if (length(eligible))
  fit_feature_models(S, composition_meta, discovery_cohorts, "composition_macrophage") else data.table()

shape <- fit_axis_shape_models(S, discovery_meta, discovery_cohorts)
omnibus <- categorical_omnibus(S, discovery_meta, discovery_cohorts)

cohort_effects <- rbindlist(
  list(linear, marginal_f, marginal_n, composition_pc, composition_mac),
  use.names = TRUE, fill = TRUE
)
cohort_effects[, `:=`(view = "hotspot_program", score_definition = "weighted_mean_z")]
write_tsv(cohort_effects, file.path(tmp, "cohort_effects.tsv"))
write_tsv(shape, file.path(tmp, "shape_cohort_effects.tsv"))
write_tsv(omnibus, file.path(tmp, "omnibus_cohort_effects.tsv"))

message("[6/10] Meta-analysing every arm")
meta_effects <- rbindlist(list(
  meta_analyze(linear), meta_analyze(marginal_f), meta_analyze(marginal_n),
  meta_analyze(composition_pc),
  if (nrow(composition_mac)) meta_analyze(composition_mac) else data.table()
), use.names = TRUE, fill = TRUE)
meta_effects[, `:=`(
  multiplicity_family = paste("hotspot_program::weighted_mean_z", model_kind, sep = "::"),
  family_size = .N,
  q_value = p.adjust(p_value_meta, method = "BH", n = .N)
), by = model_kind]

shape_meta <- meta_analyze_shape(shape)
shape_meta[estimable == TRUE, `:=`(
  q_joint = p.adjust(p_joint_meta, method = "BH", n = .N),
  q_linear = p.adjust(p_linear_meta, method = "BH", n = .N),
  q_quadratic = p.adjust(p_quadratic_meta, method = "BH", n = .N)
)]
omnibus_combined <- combine_omnibus(omnibus)
omnibus_combined[, q_combined := p.adjust(p_combined, method = "BH", n = .N), by = axis]
# Calibration status, measured not assumed. A permutation null (10 joint
# within-cohort permutations, scripts/.../calibrate_null.R) gave a mean of 0.9
# false positives per 113 programs on the fibrosis axis and 37.9 on the NAS axis.
# HC3 Wald on a categorical block is anticonservative when levels are thin, and
# NAS carries up to nine levels per cohort against five for fibrosis. The NAS
# omnibus therefore cannot support a count and is excluded from every support
# state below; it stays in the output with its null rate attached so the reason
# travels with the number.
omnibus_combined[, `:=`(
  null_false_positive_rate = fifelse(axis == "nas", 37.9, 0.9),
  calibrated = axis != "nas"
)]

write_tsv(meta_effects, file.path(tmp, "meta_effects.tsv"))
write_tsv(shape_meta, file.path(tmp, "shape_meta_effects.tsv"))
write_tsv(omnibus_combined, file.path(tmp, "omnibus_combined.tsv"))

message("[7/10] Assembling the three-view axis comparison")
linear_wide <- dcast(
  meta_effects[model_kind == "primary" & axis %in% c("fibrosis", "nas")],
  feature_id ~ axis,
  value.var = c("beta_meta", "se_meta", "ci_lower", "ci_upper", "p_value_meta",
                "q_value", "tau2", "n_cohorts", "n_same_direction")
)
comp_wide <- dcast(
  meta_effects[model_kind == "composition_pc" & axis %in% c("fibrosis", "nas")],
  feature_id ~ axis, value.var = c("beta_meta", "q_value")
)
setnames(comp_wide, setdiff(names(comp_wide), "feature_id"),
         paste0("composition_pc_", setdiff(names(comp_wide), "feature_id")))
shape_wide <- dcast(
  shape_meta[estimable == TRUE], feature_id ~ axis,
  value.var = c("beta_linear_meta", "beta_quadratic_meta", "p_joint_meta", "q_joint",
                "q_quadratic")
)
omni_wide <- dcast(omnibus_combined, feature_id ~ axis,
                   value.var = c("combined_chisq", "combined_df", "p_combined", "q_combined"))

axis_map <- Reduce(function(x, y) merge(x, y, by = "feature_id", all.x = TRUE),
                   list(registry[, .(feature_id = program_uid, cell_type, module,
                                     module_name, robust_display)],
                        program_scores$coverage[, .(feature_id, weight_coverage, testable)],
                        linear_wide, comp_wide, shape_wide, omni_wide))

for (axis in c("fibrosis", "nas")) {
  se <- axis_map[[paste0("se_meta_", axis)]]
  k <- axis_map[[paste0("n_cohorts_", axis)]]
  axis_map[[paste0("mde_", axis)]] <- minimum_detectable_effect(se, k)
  axis_map[[paste0("informative_null_", axis)]] <- informative_null(
    axis_map[[paste0("ci_lower_", axis)]], axis_map[[paste0("ci_upper_", axis)]],
    observability_effect_threshold
  )
  # Support state, built only from views whose null rate has been measured. The
  # NAS categorical omnibus is excluded because its measured null rate is 37.9
  # per 113; the fibrosis omnibus is retained at 0.9. Both counts are still
  # written out, so a reader can see what the excluded view would have said.
  lin <- axis_map[[paste0("q_value_", axis)]] < 0.05
  shp <- axis_map[[paste0("q_joint_", axis)]] < 0.05
  omn <- axis_map[[paste0("q_combined_", axis)]] < 0.05
  omnibus_calibrated <- axis != "nas"
  views <- if (omnibus_calibrated) lin | shp | omn else lin | shp
  all_views <- if (omnibus_calibrated) lin & shp & omn else lin & shp
  axis_map[[paste0("support_", axis)]] <- fcase(
    !axis_map$testable, "untestable",
    all_views, "supported_all_calibrated_views",
    views, paste0(
      "supported_", fifelse(lin, "L", ""), fifelse(shp, "S", ""),
      fifelse(omn & omnibus_calibrated, "O", "")
    ),
    axis_map[[paste0("informative_null_", axis)]], "informative_null",
    default = "indeterminate"
  )
  axis_map[[paste0("omnibus_calibrated_", axis)]] <- omnibus_calibrated
}
write_tsv(axis_map, file.path(tmp, "program_axis_map.tsv"))

message("[8/10] Testing whether fibrosis and NAS are separable dimensions")
sep_rows <- list()
for (cohort in discovery_cohorts) {
  md <- discovery_meta[dataset == cohort]
  result <- canonical_separability(
    S[, md$sample_id, drop = FALSE],
    md[, .(fibrosis_stage, nas_score)],
    n_components = separability_n_components,
    n_permutations = n_permutations,
    seed = seed + which(discovery_cohorts == cohort)
  )
  sep_rows[[length(sep_rows) + 1L]] <- data.table(
    cohort = cohort, n_samples = result$n_samples, n_features = result$n_features,
    n_components = result$n_components,
    component = seq_along(result$canonical_correlations),
    canonical_correlation = result$canonical_correlations,
    empirical_p = result$empirical_p, null_median = result$null_median
  )
}
separability <- rbindlist(sep_rows)
separability[, p_holm_within_component := p.adjust(empirical_p, method = "holm"), by = component]
write_tsv(separability, file.path(tmp, "axis_separability.tsv"))

message("[9/10] Writing manifests and the discovery seal")
# Every count carries the measured false-positive rate of the view that produced
# it. A count without a null rate is not reportable in this workstream.
summary_rows <- rbindlist(lapply(c("fibrosis", "nas"), function(axis) {
  data.table(
    axis = axis,
    n_testable = sum(axis_map$testable),
    n_linear_supported = sum(axis_map[[paste0("q_value_", axis)]] < 0.05, na.rm = TRUE),
    linear_null_rate = if (axis == "fibrosis") 0.1 else 0.2,
    n_shape_supported = sum(axis_map[[paste0("q_joint_", axis)]] < 0.05, na.rm = TRUE),
    shape_null_rate = 0.0,
    shape_reference = "F(2, k-1)",
    n_omnibus_supported = sum(axis_map[[paste0("q_combined_", axis)]] < 0.05, na.rm = TRUE),
    omnibus_null_rate = if (axis == "fibrosis") 0.9 else 37.9,
    omnibus_calibrated = axis != "nas",
    n_quadratic_supported = sum(axis_map[[paste0("q_quadratic_", axis)]] < 0.05, na.rm = TRUE),
    n_composition_pc_supported = sum(
      axis_map[[paste0("composition_pc_q_value_", axis)]] < 0.05, na.rm = TRUE),
    n_informative_null = sum(axis_map[[paste0("informative_null_", axis)]], na.rm = TRUE),
    median_mde = median(axis_map[[paste0("mde_", axis)]], na.rm = TRUE)
  )
}))
write_tsv(summary_rows, file.path(tmp, "axis_support_summary.tsv"))

inputs <- c(nonholdout_dge, nonholdout_meta, nonholdout_composition,
            program_registry, program_membership, gene_annotation, composition_accepted_seal)
write_tsv(data.table(
  input = c("nonholdout_dge", "nonholdout_meta", "nonholdout_composition",
            "program_registry", "program_membership", "gene_annotation",
            "composition_acceptance_seal"),
  path = inputs, sha256 = vapply(inputs, sha256_file, character(1))
), file.path(tmp, "input_manifest.tsv"))
writeLines(capture.output(sessionInfo()), file.path(tmp, "sessionInfo.txt"))

jsonlite::write_json(list(
  state = "DISCOVERY_COMPLETE_HOLDOUT_UNOPENED",
  release_id = release_id, workstream_id = workstream_id,
  composition_method = seal$accepted_method,
  eligible_logratio_lineages = eligible,
  composition_specification = if (length(eligible) >= 2L)
    "log_ratio" else "clr_principal_components",
  n_discovery = nrow(discovery_meta),
  n_discovery_participants = uniqueN(discovery_meta$participant_id),
  n_testable_programs = length(testable_programs),
  n_permutations = n_permutations,
  holdout_accessed = FALSE
), file.path(tmp, "DISCOVERY_READY.json"), pretty = TRUE, auto_unbox = TRUE)

artifacts <- setdiff(list.files(tmp, recursive = TRUE), "artifact_manifest.tsv")
write_tsv(data.table(
  artifact = artifacts,
  sha256 = vapply(file.path(tmp, artifacts), sha256_file, character(1)),
  size_bytes = file.info(file.path(tmp, artifacts))$size
), file.path(tmp, "artifact_manifest.tsv"))

message("[10/10] Publishing discovery")
publish_dir(tmp, discovery_root)
Sys.chmod(list.files(discovery_root, recursive = TRUE, full.names = TRUE), mode = "0440")
cat("V9_DISCOVERY_COMPLETE\t", discovery_root, "\n", sep = "")
