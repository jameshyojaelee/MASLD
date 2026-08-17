#!/usr/bin/env Rscript

# T1: is a program's bulk fibrosis association ACTIVITY or ABUNDANCE?
#
# THE QUESTION
# The bulk arm associates each frozen Hotspot program with fibrosis stage in
# whole-liver RNA. Whole liver mixes two different things: a cell type doing
# something different, and there being more or less of that cell type. The bulk
# arm cannot separate them, because its cell-type proportions are deconvolved
# from the same expression matrix the program scores come from, so the
# composition adjustment and the thing being adjusted share their noise.
#
# WHY THIS CONSTRUCTION
# Scoring a program only within the cells of its own cell type removes
# composition by construction rather than by adjustment. A hepatocyte program
# scored across hepatocytes cannot move because a donor has more hepatocytes.
# The abundance arm is then the host cell type's own compositional abundance in
# the same donors, which is what the bulk composition term is trying to capture.
#
# WHAT THIS CANNOT DO, STATED UP FRONT
#   1. NAS is unavailable at donor level for every single-cell donor, so the
#      only histologic axis here is documented fibrosis stage.
#   2. Documented F stage exists in exactly one single-cell dataset
#      (GSE202379). A dataset random effect is therefore not estimable and none
#      is fitted. This is the same reason all_modules.tsv carries an empty
#      F_stage_documented_only_beta for all 117 programs.
#   3. The canonical frozen donor score table contains no F4 donor. All nine
#      documented F4 donors fall outside the frozen 64-donor registry, so the
#      activity arm spans F0 to F3 only while the bulk association it is being
#      compared against spans F0 to F4.
#   4. "Activity" here means within a broad cell-type label. A shift in the mix
#      of sub-states inside that label is abundance at a finer grain and is not
#      distinguishable from activity by this design.
#
# CALIBRATION
# Every count below carries a permutation null false-call rate measured by
# calibrate_view() and passed through gate_view() before it is written.

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_dir <- dirname(normalizePath(script_arg))
source(file.path(script_dir, "config.R"))
source(file.path(script_dir, "analysis_lib.R"))
source(file.path(script_dir, "calibration_lib.R"))
source(file.path(script_dir, "systems_contract.R"))
source(file.path(script_dir, "t1_cellstate_lib.R"))

set.seed(seed)
n_null <- as.integer(Sys.getenv("MASLD_T1_PERMUTATIONS", "2000"))

# Frozen before the first fit.
min_cells_primary <- 20L
min_cells_sensitivity <- c(0L, 50L)
alpha <- 0.05
inclusion_criteria <- c("evidence_interpretation", "experiment_routing")
assert_inclusion_criterion(inclusion_criteria)

output_root <- Sys.getenv("MASLD_T1_OUTPUT",
                          file.path(workstream_root, "systems", "T1_cellstate"))
if (dir.exists(output_root)) fail("Refusing to overwrite: ", output_root)

score_primary <- file.path(program_root, "donor_program_scores_primary.tsv")
score_equal_run <- file.path(program_root, "donor_program_scores_equal_run.tsv")
fstage_map_path <- file.path(program_root, "documented_fstage_donor_map.tsv")
lineage_counts_path <- file.path(
  project_root,
  "Analysis/SingleCell/results_gpu_v2_phase05/ccc/stage_trajectory_v2/lineage_cell_counts_v2_dc.tsv"
)
# Run-level cell counts extracted from the Hotspot per-cell score tables by
# t1_cell_counts.py. The parquet files themselves are what gets hashed into the
# input manifest; this table is a deterministic reduction of them and is
# republished inside the sealed directory.
run_cell_counts_path <- Sys.getenv("MASLD_T1_CELL_COUNTS", "")
if (!nzchar(run_cell_counts_path)) {
  fail("MASLD_T1_CELL_COUNTS must point at the run-level cell counts from t1_cell_counts.py")
}
cell_score_parquets <- file.path(
  project_root, "Analysis/SingleCell/results_gpu_v2/hotspot_modules",
  T1_CELL_TYPES, "cell_scores.parquet")
donor_pairing_files <- file.path(project_root, c(
  "data/GSE244832/metadata/donor_pairing.csv", "data/GSE185477/metadata/donor_pairing.csv",
  "data/GSE202379/metadata/donor_pairing.csv", "data/GSE136103/metadata/donor_pairing.csv"))
bulk_map_path <- file.path(discovery_root, "program_axis_map.tsv")

for (p in c(score_primary, score_equal_run, fstage_map_path, lineage_counts_path,
            run_cell_counts_path, bulk_map_path, program_registry, program_membership,
            cell_score_parquets, donor_pairing_files)) {
  if (!file.exists(p)) fail("Missing input: ", p)
}

message("[1/8] Contract assertions")
assert_frozen_programs(program_registry, program_membership, 117L)
assert_v8_untouched(project_root, file.path(
  candidate_root, "workstreams/BULK-PROGRAM-MAP-v8/contract/analysis_contract.json"))

message("[2/8] Roster")
scores_long_primary <- fread(score_primary)
fstage <- fread(fstage_map_path)
assert_true(uniqueN(fstage$dataset) == 1L,
            "Documented F stage must come from a single dataset for this design to be described honestly")

# The three available rosters and the deliberate choice between them.
#   85 donors  : donor_collapse/donor_scores_all_weighted.tsv, everything scored
#   64 donors  : the frozen program-context-v2 primary table, the release roster
#   37 donors  : the 64 that also carry a documented F stage
# The frozen 64 is the release's own primary score table and the only one the
# frozen analysis_specification.tsv names, so the roster is the intersection of
# that table with documented F stage. Widening it to the 85 would import donors
# the release did not freeze in order to gain F4 donors, which is exactly the
# kind of post hoc roster change the contract exists to stop.
roster <- intersect(unique(scores_long_primary$donor), fstage$donor)
donor_meta <- merge(
  unique(scores_long_primary[donor %in% roster, .(donor, dataset)]),
  fstage[, .(donor, F_stage_documented)], by = "donor")
setnames(donor_meta, "donor", "sample_id")
donor_meta[, participant_id := sample_id]
# Non-vacuous unit assertion: the donor ids are the output of the shared
# run-to-donor collapse, so the check that matters is that no donor is
# represented twice and that every roster donor came through that collapse.
assert_true(uniqueN(donor_meta$sample_id) == nrow(donor_meta),
            "A donor appears twice in the roster")
assert_true(all(grepl("^GSE[0-9]+_", donor_meta$sample_id)),
            "Roster ids are not dataset-prefixed donor ids from the donor collapse")
assert_true(!any(donor_meta$sample_id %in% unique(scores_long_primary$dataset)),
            "A roster id is a dataset id, so the unit column has been filled from the wrong field")
message("      roster donors: ", nrow(donor_meta), " from ", uniqueN(donor_meta$dataset),
        " dataset; F stages present: ",
        paste(sort(unique(donor_meta$F_stage_documented)), collapse = ","))
missing_documented <- setdiff(fstage$donor, roster)
message("      documented donors outside the frozen score table: ", length(missing_documented),
        " (F stages ",
        paste(sort(unique(fstage[donor %in% missing_documented, F_stage_documented])), collapse = ","),
        ")")

message("[3/8] Cells per donor per cell type")
# The same run-to-donor collapse the score tables went through. Runs are
# sequencing runs, not donors: GSE202379 alone is 67 runs over 46 donors, so
# summing cells per run within donor is the only count that matches the score.
source(file.path(project_root, "Analysis/SingleCell/scripts/lib_donor_collapse.R"))
run_to_donor <- build_srr_to_donor_map(project_root)
run_counts <- fread(run_cell_counts_path)
assert_true(all(c("cell_type", "run", "n_cells") %in% names(run_counts)),
            "Run-level cell count table is missing a required column")
run_counts[, donor := fifelse(run %in% names(run_to_donor), run_to_donor[run], run)]
cell_counts <- run_counts[, .(n_cells = as.numeric(sum(n_cells)), n_runs = .N),
                          by = .(cell_type, sample_id = donor)]
cell_counts <- cell_counts[sample_id %in% donor_meta$sample_id]
scored_pairs <- unique(scores_long_primary[donor %in% roster, .(sample_id = donor, cell_type)])
cell_counts <- merge(scored_pairs, cell_counts[, .(sample_id, cell_type, n_cells, n_runs)],
                     by = c("sample_id", "cell_type"), all.x = TRUE)
assert_true(!anyNA(cell_counts$n_cells),
            "A scored donor-cell-type pair has no cell count; the precision gate cannot be applied")

# Independent check that the exact counts and the phase05 lineage table describe
# the same cells, since the abundance arm uses the lineage table.
lineage <- fread(lineage_counts_path)
setnames(lineage, "sample", "sample_id")
lineage <- lineage[sample_id %in% donor_meta$sample_id]
assert_true(nrow(lineage) == nrow(donor_meta), "A roster donor is absent from the lineage count table")
proxy <- melt(lineage, id.vars = "sample_id", variable.name = "column", value.name = "proxy_cells")
proxy[, cell_type := names(T1_LINEAGE_COLUMN)[match(as.character(column), T1_LINEAGE_COLUMN)]]
proxy <- proxy[!is.na(cell_type)]
agreement <- merge(cell_counts, proxy[, .(sample_id, cell_type, proxy_cells)],
                   by = c("sample_id", "cell_type"))
agreement_by_type <- agreement[, .(
  n = .N, pearson = stats::cor(as.numeric(n_cells), as.numeric(proxy_cells))), by = cell_type]
message("      exact vs lineage-table cell counts, Pearson by cell type: ",
        paste(sprintf("%s %.3f", agreement_by_type$cell_type,
                      agreement_by_type$pearson), collapse = "; "))

message("[4/8] Composition")
count_matrix <- as.matrix(lineage[match(donor_meta$sample_id, lineage$sample_id),
                                  setdiff(names(lineage), "sample_id"), with = FALSE])
rownames(count_matrix) <- donor_meta$sample_id
clr_result <- clr_cell_counts(count_matrix)
message("      lineages retained: ", paste(clr_result$retained_lineages, collapse = ", "))
message("      lineages dropped as structurally absent: ",
        if (length(clr_result$dropped_lineages)) paste(clr_result$dropped_lineages, collapse = ", ") else "none")
assert_true(all(T1_LINEAGE_COLUMN %in% clr_result$retained_lineages),
            "A cell type carrying Hotspot programs was dropped from the composition")

# ---------------------------------------------------------------------------
# The fitting procedure, written once so the observed run and every permutation
# go through identical code.
# ---------------------------------------------------------------------------
score_matrix_by_type <- function(long, cell_type_value, donors) {
  sub <- long[cell_type == cell_type_value & donor %in% donors,
              .(program_uid, donor, score)]
  wide <- dcast(sub, program_uid ~ donor, value.var = "score")
  Y <- as.matrix(wide[, -1L, with = FALSE])
  rownames(Y) <- wide$program_uid
  Y[, donors[donors %in% colnames(Y)], drop = FALSE]
}

fit_one_view <- function(long, meta, min_cells, counts, clr_matrix, weighted = FALSE,
                         adjust_cell_yield = FALSE) {
  activity <- list()
  abundance <- list()
  for (ct in T1_CELL_TYPES) {
    eligible <- counts[cell_type == ct & n_cells >= min_cells, sample_id]
    donors <- meta[sample_id %in% eligible, sample_id]
    if (length(donors) < 6L) next
    Y <- score_matrix_by_type(long, ct, donors)
    if (!nrow(Y)) next
    donors <- colnames(Y)
    md <- meta[match(donors, sample_id)]
    standardized <- t1_standardize(Y)
    Ys <- standardized$Y
    if (!nrow(Ys)) next
    donor_cells <- counts[cell_type == ct][match(donors, sample_id), as.numeric(n_cells)]
    X <- cbind(`(Intercept)` = 1, F_stage_documented = md$F_stage_documented)
    if (adjust_cell_yield) {
      # Sensitivity only. Cell yield is a technical quantity that also carries
      # signal: in this roster T-cell yield falls with F stage, so adjusting for
      # it removes a nuisance for some cell types and part of the exposure for
      # others. The pre-declared primary model does not include it.
      X <- cbind(X, log_cells = log(donor_cells) - mean(log(donor_cells)))
    }
    # The abundance arm is standardised on the same donors for the same reason
    # the activity arm is: the frozen observability threshold is expressed in
    # score SD per unit stage, and comparing a CLR log-ratio coefficient against
    # a threshold defined in SD units would be a unit error. The native CLR
    # effect is recoverable as beta_abundance * clr_sd.
    Ya <- matrix(clr_matrix[donors, T1_LINEAGE_COLUMN[[ct]]], nrow = 1L,
                 dimnames = list(ct, donors))
    clr_sd <- stats::sd(Ya[1L, ])
    Ya <- zscore_rows(Ya)
    if (weighted) {
      w <- sqrt(donor_cells)
      Ys <- sweep(Ys, 2L, w, "*")
      Ya <- sweep(Ya, 2L, w, "*")
      X <- X * w
    }
    act <- t1_term_effects(Ys, X, "F_stage_documented")
    act[, cell_type := ct]
    act[, n_donors := length(donors)]
    abu <- t1_term_effects(Ya, X, "F_stage_documented")
    abu[, cell_type := ct]
    abu[, n_donors := length(donors)]
    abu[, clr_sd := clr_sd]
    abu[, median_cells_gated := stats::median(
      counts[cell_type == ct][match(donors, sample_id), as.numeric(n_cells)])]
    activity[[ct]] <- act
    abundance[[ct]] <- abu
  }
  list(activity = rbindlist(activity), abundance = rbindlist(abundance))
}

annotate <- function(effects) {
  d <- copy(effects)
  d[, ci_lower := beta - stats::qt(1 - alpha / 2, df = df) * se]
  d[, ci_upper := beta + stats::qt(1 - alpha / 2, df = df) * se]
  # minimum_detectable_effect() derives its t reference from n_cohorts - 1. This
  # design has one cohort and a residual-df t test, so the residual df is passed
  # through as n_cohorts + 1 to make the function use exactly that df.
  d[, mde := minimum_detectable_effect(se, n_cohorts = df + 1L)]
  d[]
}

message("[5/8] Observed fit, primary view (>= ", min_cells_primary, " cells)")
observed <- fit_one_view(scores_long_primary, donor_meta, min_cells_primary,
                         cell_counts, clr_result$clr)
activity <- annotate(observed$activity)
abundance <- annotate(observed$abundance)
activity[, q_value := p.adjust(p_value, method = "BH")]
abundance[, q_value := p.adjust(p_value, method = "BH")]
message("      programs modelled: ", nrow(activity), " of 117; donors per cell type: ",
        paste(sprintf("%s %d", activity[, .SD[1], by = cell_type]$cell_type,
                      activity[, .SD[1], by = cell_type]$n_donors), collapse = "; "))

message("[6/8] Permutation calibration (", n_null, " permutations)")
permuted_effects <- function(i) {
  set.seed(seed + i)
  permuted_meta <- permute_histology_within_cohort(
    donor_meta, "F_stage_documented", cohort_column = "dataset")
  fit_one_view(scores_long_primary, permuted_meta, min_cells_primary,
               cell_counts, clr_result$clr)
}
# One permutation stream feeds both views so the two calibrations describe the
# same null, and the stream is cached so the two calibrate_view() calls do not
# refit it twice.
null_cache <- vector("list", n_null)
for (i in seq_len(n_null)) null_cache[[i]] <- permuted_effects(i)

calibration_activity <- calibrate_view(
  observed_p = activity$p_value,
  permuted_p_fn = function(i) null_cache[[i]]$activity$p_value,
  n_reps = n_null, n_tests = nrow(activity), alpha = alpha,
  label = "activity_within_cell_type", observed_p_is_empirical = FALSE)
calibration_abundance <- calibrate_view(
  observed_p = abundance$p_value,
  permuted_p_fn = function(i) null_cache[[i]]$abundance$p_value,
  n_reps = n_null, n_tests = nrow(abundance), alpha = alpha,
  label = "cell_type_abundance", observed_p_is_empirical = FALSE)

gate_activity <- gate_view(calibration_activity, strict = FALSE)
gate_abundance <- gate_view(calibration_abundance, strict = FALSE)
message("      activity: ", calibration_activity$observed_calls, " observed calls, null mean ",
        signif(calibration_activity$null_call_mean, 3), " of ", nrow(activity),
        " | reportable: ", calibration_activity$reportable)
message("      abundance: ", calibration_abundance$observed_calls, " observed calls, null mean ",
        signif(calibration_abundance$null_call_mean, 3), " of ", nrow(abundance),
        " | reportable: ", calibration_abundance$reportable)
calibration_rows <- rbindlist(list(calibration_row(calibration_activity),
                                   calibration_row(calibration_abundance)))

message("[7/8] Attribution, aggregates and sensitivity")
bulk <- fread(bulk_map_path)
bulk_slim <- bulk[, .(feature_id, module, module_name,
                      beta_bulk_fibrosis = beta_meta_fibrosis,
                      q_bulk_fibrosis = q_value_fibrosis,
                      support_bulk_fibrosis = support_fibrosis,
                      mde_bulk_fibrosis = mde_fibrosis,
                      composition_pc_beta_bulk = composition_pc_beta_meta_fibrosis,
                      composition_pc_q_bulk = composition_pc_q_value_fibrosis)]
bulk_slim[, bulk_supported := grepl("^supported", support_bulk_fibrosis)]

map <- merge(
  activity[, .(feature_id, cell_type, n_donors_activity = n_donors,
               df_activity = df, beta_activity = beta, se_activity = se,
               ci_lower_activity = ci_lower, ci_upper_activity = ci_upper,
               p_activity = p_value, q_activity = q_value, mde_activity = mde)],
  abundance[, .(cell_type, beta_abundance = beta, se_abundance = se,
                ci_lower_abundance = ci_lower, ci_upper_abundance = ci_upper,
                p_abundance = p_value, q_abundance = q_value,
                mde_abundance = mde, n_donors_abundance = n_donors,
                clr_sd_abundance = clr_sd, median_cells_gated)],
  by = "cell_type")
map <- merge(map, bulk_slim, by = "feature_id", all.x = TRUE)
assert_true(all(map$n_donors_activity == map$n_donors_abundance),
            "The two arms were fitted on different donor sets, which would manufacture an asymmetry")

attributed <- t1_classify(map, observability_effect_threshold, alpha = alpha)
setorder(attributed, cell_type, feature_id)
# Every attribution that rests on an activity null has to carry the smallest
# activity effect this design could have resolved, or the null is not readable.
assert_negatives_have_mde(
  data.table(
    support = fifelse(
      attributed$attribution %in% c("abundance", "neither_arm_resolved",
                                    "indeterminate_underpowered"),
      "indeterminate", "supported"),
    mde = attributed$mde_activity),
  "support", "mde")

# ---------------------------------------------------------------------------
# Aggregate statistics with an exact permutation null. These carry the weight
# the per-program table cannot at this donor count.
# ---------------------------------------------------------------------------
bulk_beta <- bulk_slim[, .(feature_id, beta_bulk_fibrosis, bulk_supported)]
aggregate_observed <- t1_aggregate_statistics(activity, bulk_beta, T1_CELL_TYPES)
aggregate_null <- lapply(seq_len(n_null), function(i) {
  t1_aggregate_statistics(annotate(null_cache[[i]]$activity), bulk_beta, T1_CELL_TYPES)
})
aggregate_statistics <- c("mean_t_squared", "spearman_vs_bulk", "sign_agreement")
aggregate_table <- rbindlist(lapply(seq_len(nrow(aggregate_observed)), function(r) {
  this_stratum <- aggregate_observed$stratum[[r]]
  rbindlist(lapply(aggregate_statistics, function(s) {
    null_values <- vapply(aggregate_null, function(x) {
      v <- x[["stratum"]] == this_stratum
      if (!any(v)) NA_real_ else as.numeric(x[[s]][which(v)[[1L]]])
    }, numeric(1))
    null_values <- null_values[is.finite(null_values)]
    observed_value <- aggregate_observed[[s]][[r]]
    if (!is.finite(observed_value) || length(null_values) < 10L) return(NULL)
    data.table(
      stratum = this_stratum, statistic = s, n_programs = aggregate_observed$n_programs[[r]],
      observed = observed_value,
      null_mean = mean(null_values), null_p95 = as.numeric(stats::quantile(null_values, 0.95)),
      n_null_finite = length(null_values),
      # One-sided upward: the alternative is that within-cell-type effects are
      # larger, or better aligned with the bulk effects, than chance.
      p_empirical = empirical_p(observed_value, null_values, "greater"))
  }), fill = TRUE)
}), fill = TRUE)
aggregate_table[, q_value := p.adjust(p_empirical, method = "BH")]
calibration_aggregate <- calibrate_view(
  observed_p = aggregate_table$p_empirical,
  permuted_p_fn = function(i) {
    # An empirical p computed against its own permutation null is uniform under
    # the null by construction, so the calibration arm of this view is not
    # informative and the check that matters is resolution: an empirical p
    # cannot fall below 1/(n_null + 1), and BH over this family needs
    # 0.05/n_tests. n_reps is therefore set to the permutation count that
    # actually backs the observed p-values, not to the cost of this draw.
    stats::runif(nrow(aggregate_table))
  },
  n_reps = n_null, n_tests = nrow(aggregate_table), alpha = alpha,
  label = "aggregate_activity_statistics", observed_p_is_empirical = TRUE)
gate_aggregate <- gate_view(calibration_aggregate, strict = FALSE)
calibration_rows <- rbind(calibration_rows, calibration_row(calibration_aggregate))
message("      aggregate statistics: ", sum(aggregate_table$q_value < alpha),
        " of ", nrow(aggregate_table), " reject; resolvable: ", calibration_aggregate$resolvable)

# ---------------------------------------------------------------------------
# Observability control on a different, larger roster.
#
# The fibrosis arm is limited to 37 donors in one dataset. If it returns nothing
# there are two explanations, and they have different consequences: the
# within-cell-type view carries no stage signal, or this roster is too small to
# show one. The frozen analysis_specification.tsv already names a coarse
# disease-stage model, score ~ stage_ordinal + factor(dataset), which runs on
# all 64 frozen donors across five datasets. Coarse disease stage is NOT
# fibrosis stage and this arm attributes nothing; it exists only to say whether
# the measurement is capable of resolving anything at this atlas size.
# ---------------------------------------------------------------------------
coarse_meta <- unique(scores_long_primary[
  exclude_stage_analysis == FALSE, .(sample_id = donor, dataset, stage_ordinal)])
assert_true(uniqueN(coarse_meta$sample_id) == nrow(coarse_meta),
            "A donor carries two coarse stage values")
coarse_counts <- run_counts[, .(n_cells = as.numeric(sum(n_cells))), by = .(cell_type, sample_id = donor)]
coarse_counts <- coarse_counts[sample_id %in% coarse_meta$sample_id]
fit_coarse <- function(meta) {
  out <- list()
  for (ct in T1_CELL_TYPES) {
    eligible <- coarse_counts[cell_type == ct & n_cells >= min_cells_primary, sample_id]
    donors <- meta[sample_id %in% eligible, sample_id]
    if (length(donors) < 8L) next
    Y <- score_matrix_by_type(scores_long_primary, ct, donors)
    if (!nrow(Y)) next
    donors <- colnames(Y)
    md <- meta[match(donors, sample_id)]
    Ys <- t1_standardize(Y)$Y
    if (!nrow(Ys)) next
    # Dataset enters as a fixed effect exactly as the frozen specification says.
    # Levels with a single donor are absorbed and cost a degree of freedom; the
    # design is dropped rather than silently rank-reduced if that leaves none.
    X <- stats::model.matrix(~ stage_ordinal + factor(dataset), data = md)
    colnames(X)[colnames(X) == "stage_ordinal"] <- "stage_ordinal"
    if (qr(X)$rank != ncol(X) || nrow(X) <= ncol(X)) next
    e <- t1_term_effects(Ys, X, "stage_ordinal")
    e[, cell_type := ct][, n_donors := length(donors)][, n_datasets := uniqueN(md$dataset)]
    out[[ct]] <- e
  }
  rbindlist(out)
}
coarse_observed <- annotate(fit_coarse(coarse_meta))
coarse_observed[, q_value := p.adjust(p_value, method = "BH")]
coarse_null_p <- function(i) {
  set.seed(seed + 500000L + i)
  annotate(fit_coarse(permute_histology_within_cohort(
    coarse_meta, "stage_ordinal", cohort_column = "dataset")))$p_value
}
calibration_coarse <- calibrate_view(
  observed_p = coarse_observed$p_value, permuted_p_fn = coarse_null_p,
  n_reps = min(n_null, 500L), n_tests = nrow(coarse_observed), alpha = alpha,
  label = "coarse_disease_stage_observability_control", observed_p_is_empirical = FALSE)
gate_coarse <- gate_view(calibration_coarse, strict = FALSE)
calibration_rows <- rbind(calibration_rows, calibration_row(calibration_coarse))
message("      coarse-stage control: ", calibration_coarse$observed_calls, " of ",
        nrow(coarse_observed), " called on ", max(coarse_observed$n_donors),
        " donors; null mean ", signif(calibration_coarse$null_call_mean, 3),
        "; reportable: ", calibration_coarse$reportable)

# Sensitivity arms, all pre-declared: the two other cell-count gates, the
# equal-run score definition the frozen spec names as its sensitivity arm, and a
# precision-weighted fit that uses every donor's score in proportion to the
# number of cells it was averaged over instead of gating on it.
sensitivity <- list()
run_sensitivity <- function(label, long, min_cells, weighted, adjust_cell_yield = FALSE) {
  fitted <- fit_one_view(long, donor_meta, min_cells, cell_counts, clr_result$clr,
                         weighted, adjust_cell_yield)
  if (!nrow(fitted$activity)) return(NULL)
  a <- annotate(fitted$activity)
  b <- annotate(fitted$abundance)
  a[, q_value := p.adjust(p_value, method = "BH")]
  b[, q_value := p.adjust(p_value, method = "BH")]
  data.table(arm = label, n_programs = nrow(a),
             n_activity_called = sum(a$q_value < alpha, na.rm = TRUE),
             n_abundance_called = sum(b$q_value < alpha, na.rm = TRUE),
             median_n_donors = stats::median(a$n_donors),
             concordance_with_primary = stats::cor(
               a[match(activity$feature_id, feature_id), beta],
               activity$beta, use = "complete.obs"))
}
scores_long_equal <- fread(score_equal_run)
for (g in min_cells_sensitivity) {
  sensitivity[[paste0("gate_", g)]] <- run_sensitivity(
    paste0("cell_gate_", g), scores_long_primary, g, FALSE)
}
sensitivity[["equal_run"]] <- run_sensitivity(
  "equal_run_score_definition", scores_long_equal, min_cells_primary, FALSE)
sensitivity[["weighted"]] <- run_sensitivity(
  "precision_weighted_no_gate", scores_long_primary, 0L, TRUE)
sensitivity[["cell_yield_adjusted"]] <- run_sensitivity(
  "cell_yield_adjusted", scores_long_primary, min_cells_primary, FALSE, TRUE)
sensitivity_table <- rbindlist(sensitivity, fill = TRUE)

# Is the precision gate itself stage-dependent? If a cell type's yield falls as
# stage rises, the donors the gate keeps are not a stage-independent subset and
# the activity arm for that cell type is conditioned on the exposure. This is
# reported whether or not it changes a call.
cell_yield_diagnostic <- rbindlist(lapply(T1_CELL_TYPES, function(ct) {
  eligible <- cell_counts[cell_type == ct & n_cells >= min_cells_primary, sample_id]
  md <- merge(donor_meta[sample_id %in% eligible],
              cell_counts[cell_type == ct, .(sample_id, n_cells)], by = "sample_id")
  ungated <- merge(donor_meta, cell_counts[cell_type == ct, .(sample_id, n_cells)],
                   by = "sample_id")
  data.table(
    cell_type = ct, n_gated = nrow(md), n_before_gate = nrow(ungated),
    spearman_log_cells_vs_stage = suppressWarnings(stats::cor(
      log(md$n_cells), md$F_stage_documented, method = "spearman")),
    p_value = suppressWarnings(stats::cor.test(
      log(md$n_cells), md$F_stage_documented, method = "spearman")$p.value),
    stages_lost_to_gate = paste(sort(setdiff(ungated$F_stage_documented,
                                             md$F_stage_documented)), collapse = ","))
}))

# Supplementary abundance-only view on all 46 documented donors including F4.
# It is NOT used for attribution: its roster is not the roster the activity arm
# can use, and comparing arms across different rosters is the asymmetry this
# design exists to avoid.
lineage_all <- fread(lineage_counts_path)
setnames(lineage_all, "sample", "sample_id")
fstage_all <- fstage[donor %in% lineage_all$sample_id]
count_all <- as.matrix(lineage_all[match(fstage_all$donor, lineage_all$sample_id),
                                   setdiff(names(lineage_all), "sample_id"), with = FALSE])
rownames(count_all) <- fstage_all$donor
clr_all <- clr_cell_counts(count_all)
X_all <- cbind(`(Intercept)` = 1, F_stage_documented = fstage_all$F_stage_documented)
Y_all <- t(clr_all$clr[, T1_LINEAGE_COLUMN, drop = FALSE])
rownames(Y_all) <- names(T1_LINEAGE_COLUMN)
abundance_f4 <- annotate(t1_term_effects(Y_all, X_all, "F_stage_documented"))
abundance_f4[, q_value := p.adjust(p_value, method = "BH")]
abundance_f4[, roster := "documented_46_including_F4"]

message("[8/8] Writing")
tmp <- atomic_dir(output_root)
on.exit(if (dir.exists(tmp)) unlink(tmp, recursive = TRUE), add = TRUE)

summary_counts <- attributed[, .N, by = attribution][order(-N)]
write_tsv(attributed, file.path(tmp, "program_attribution.tsv"))
write_tsv(activity, file.path(tmp, "activity_effects.tsv"))
write_tsv(abundance, file.path(tmp, "abundance_effects.tsv"))
write_tsv(abundance_f4, file.path(tmp, "abundance_supplementary_f4_roster.tsv"))
write_tsv(calibration_rows, file.path(tmp, "calibration.tsv"))
write_tsv(aggregate_table, file.path(tmp, "aggregate_statistics.tsv"))
write_tsv(coarse_observed, file.path(tmp, "coarse_stage_observability_control.tsv"))
# The power ceiling, stated as a number rather than as a caveat.
#
# Two comparisons matter. First, the smallest resolvable activity effect against
# the frozen threshold that separates an informative null from an unobservable
# one. Second, the same quantity against the size of the bulk effects this
# system is meant to explain: a design that can only see effects twice the size
# of the ones in question cannot adjudicate them, and saying so with a number is
# more use than saying "underpowered".
#
# The donor extrapolation holds the residual variance and the F-stage
# distribution fixed and scales the standard error as 1/sqrt(n). It is an
# estimate of the order of magnitude required, not a sample-size calculation for
# a specific program.
median_supported_bulk_effect <- stats::median(
  abs(bulk_slim[bulk_supported == TRUE, beta_bulk_fibrosis]), na.rm = TRUE)
power_ceiling <- activity[, .(
  n_programs = .N, n_donors = n_donors[1L],
  median_se = stats::median(se), median_mde = stats::median(mde),
  min_mde = min(mde), observability_effect_threshold = observability_effect_threshold,
  median_mde_over_threshold = stats::median(mde) / observability_effect_threshold,
  n_programs_with_informative_null_possible = sum(mde < observability_effect_threshold),
  median_supported_bulk_effect = median_supported_bulk_effect,
  median_mde_over_bulk_effect = stats::median(mde) / median_supported_bulk_effect),
  by = cell_type]
power_ceiling[, donors_for_threshold_mde := ceiling(
  n_donors * (median_mde / observability_effect_threshold)^2)]
power_ceiling[, donors_for_bulk_effect_mde := ceiling(
  n_donors * (median_mde / median_supported_bulk_effect)^2)]
write_tsv(power_ceiling, file.path(tmp, "power_ceiling.tsv"))
write_tsv(sensitivity_table, file.path(tmp, "sensitivity.tsv"))
write_tsv(cell_yield_diagnostic, file.path(tmp, "cell_yield_vs_stage.tsv"))
write_tsv(summary_counts, file.path(tmp, "attribution_summary.tsv"))
write_tsv(donor_meta, file.path(tmp, "roster.tsv"))
write_tsv(cell_counts, file.path(tmp, "cells_per_donor_per_cell_type.tsv"))
write_tsv(agreement_by_type, file.path(tmp, "cell_count_source_agreement.tsv"))
write_tsv(data.table(
  limit = c("nas_axis", "dataset_random_effect", "f4_absent_from_activity_arm",
            "sub_state_abundance", "roster_choice", "cell_precision_gate"),
  statement = c(
    "NAS is unavailable at donor level for every single-cell donor. No NAS arm was attempted.",
    "Documented F stage exists in one dataset (GSE202379), so no dataset random effect is estimable and none was fitted.",
    paste0("All ", length(missing_documented), " documented donors outside the frozen 64-donor score table are excluded; ",
           "the activity arm spans F stages ",
           paste(sort(unique(donor_meta$F_stage_documented)), collapse = ","),
           " while the bulk association spans F0 to F4."),
    "Activity is defined within a broad cell-type label. A shift in sub-state composition inside a label is not distinguishable from activity here.",
    "The frozen program-context-v2 primary score table (64 donors) intersected with documented F stage. The 85-donor donor-collapse table was not substituted.",
    paste0("Donors contributing fewer than ", min_cells_primary,
           " cells of the program's own cell type are excluded from both arms."))
), file.path(tmp, "stated_limits.tsv"))

write_tsv(run_counts, file.path(tmp, "cells_per_run_per_cell_type.tsv"))
inputs <- c(score_primary, score_equal_run, fstage_map_path, lineage_counts_path,
            bulk_map_path, program_registry, program_membership,
            cell_score_parquets, donor_pairing_files)
write_tsv(data.table(path = inputs, sha256 = vapply(inputs, sha256_file, character(1))),
          file.path(tmp, "input_manifest.tsv"))
writeLines(capture.output(sessionInfo()), file.path(tmp, "sessionInfo.txt"))

report_lines <- c(
  sprintf("roster: %d donors, %d dataset, F stages %s", nrow(donor_meta),
          uniqueN(donor_meta$dataset),
          paste(sort(unique(donor_meta$F_stage_documented)), collapse = ",")),
  sprintf("activity view reportable: %s (null false-call mean %.3f of %d)",
          calibration_activity$reportable, calibration_activity$null_call_mean, nrow(activity)),
  sprintf("abundance view reportable: %s (null false-call mean %.3f of %d)",
          calibration_abundance$reportable, calibration_abundance$null_call_mean, nrow(abundance)))
assert_language(c(report_lines, attributed$attribution))

reportable <- calibration_activity$reportable && calibration_abundance$reportable &&
  calibration_aggregate$reportable && calibration_coarse$reportable
if (reportable) {
  assert_counts_calibrated(summary_counts, calibration_rows)
} else {
  write_tsv(data.table(
    withheld = TRUE,
    reason = paste(c(gate_activity$reasons, gate_abundance$reasons,
                     gate_aggregate$reasons, gate_coarse$reasons), collapse = "; ")),
    file.path(tmp, "COUNTS_WITHHELD.tsv"))
}

jsonlite::write_json(list(
  state = if (reportable) "T1_COMPLETE" else "T1_COMPLETE_COUNTS_WITHHELD",
  system = "T1_cellstate", workstream_id = workstream_id,
  inclusion_criteria = inclusion_criteria,
  n_donors = nrow(donor_meta), n_datasets = uniqueN(donor_meta$dataset),
  f_stages_present = sort(unique(donor_meta$F_stage_documented)),
  min_cells_gate = min_cells_primary,
  n_permutations = n_null,
  activity = list(observed_calls = calibration_activity$observed_calls,
                  null_false_call_mean = calibration_activity$null_call_mean,
                  null_false_call_max = calibration_activity$null_call_max,
                  n_tests = calibration_activity$n_tests,
                  reportable = calibration_activity$reportable),
  abundance = list(observed_calls = calibration_abundance$observed_calls,
                   null_false_call_mean = calibration_abundance$null_call_mean,
                   null_false_call_max = calibration_abundance$null_call_max,
                   n_tests = calibration_abundance$n_tests,
                   reportable = calibration_abundance$reportable),
  aggregate = list(n_rejected = sum(aggregate_table$q_value < alpha),
                   n_tests = nrow(aggregate_table),
                   reportable = calibration_aggregate$reportable),
  coarse_stage_control = list(observed_calls = calibration_coarse$observed_calls,
                              n_tests = calibration_coarse$n_tests,
                              max_donors = max(coarse_observed$n_donors),
                              null_false_call_mean = calibration_coarse$null_call_mean,
                              reportable = calibration_coarse$reportable),
  median_activity_mde = stats::median(activity$mde),
  observability_effect_threshold = observability_effect_threshold,
  attribution = as.list(setNames(summary_counts$N, summary_counts$attribution)),
  nas_attempted = FALSE, holdout_accessed = FALSE
), file.path(tmp, "T1_READY.json"), pretty = TRUE, auto_unbox = TRUE)

publish_dir(tmp, output_root)
Sys.chmod(list.files(output_root, full.names = TRUE), mode = "0440")

cat("\n=== T1 CELL-STATE ATTRIBUTION ===\n")
print(summary_counts)
cat("\ncalibration:\n")
print(calibration_rows[, .(view, observed_calls, null_call_mean, null_call_max, n_tests, reportable)])
cat("\nsensitivity:\n")
print(sensitivity_table)
cat("\ncell yield vs stage within the gated set:\n")
print(cell_yield_diagnostic)
cat("\nabundance, primary roster:\n")
print(abundance[, .(cell_type, n_donors, beta, se, p_value, q_value, mde)])
cat("\naggregate statistics:\n")
print(aggregate_table[, .(stratum, statistic, n_programs, observed = round(observed, 4),
                          null_mean = round(null_mean, 4), p_empirical, q_value)])
cat("\npower ceiling:\n")
print(power_ceiling[, .(cell_type, n_donors, median_mde = round(median_mde, 3),
                        median_mde_over_threshold = round(median_mde_over_threshold, 2),
                        median_mde_over_bulk_effect = round(median_mde_over_bulk_effect, 2),
                        donors_for_threshold_mde, donors_for_bulk_effect_mde)])
cat("\ncoarse disease-stage observability control (different estimand, attributes nothing):\n")
print(coarse_observed[, .(n_programs = .N, n_donors = n_donors[1L], n_datasets = n_datasets[1L],
                          n_called = sum(q_value < alpha)), by = cell_type])
cat("\nT1_COMPLETE\t", output_root, "\n", sep = "")
