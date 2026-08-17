#!/usr/bin/env Rscript

# T3 vocabulary invariance.
#
# WHAT THIS SYSTEM ADDS. The v9 discovery reports 27 of 113 Hotspot programs
# supported on the fibrosis axis. That count has no denominator, in two senses.
# It is not compared against what an arbitrary gene set of the same size and
# expression level achieves on this substrate, and it is not compared against
# any other way of carving the transcriptome into programs. This system supplies
# both. Six vocabularies go through one identical model, alongside 200
# size- and expression-decile-matched random gene sets per program.
#
# INCLUSION CRITERION (PAPER.md:211-214): evidence_interpretation. A supported
# program means something different once you know the matched-random rate and
# how many independent vocabularies recover the same genes.
#
# WHAT IS AND IS NOT FROZEN. The 117 Hotspot programs and their weights are
# frozen and are read, never refit. The other vocabularies are read from their
# own frozen artifacts. Nothing here rediscovers, reweights or re-selects a
# program. The random sets are generated fresh under a recorded seed.
#
# LANGUAGE. Fibrosis stage and NAS are recorded cross-sectional histologic
# scores entering as regressors. Nothing here orders donors or stages.

suppressPackageStartupMessages({
  library(edgeR)
  library(data.table)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_dir <- dirname(normalizePath(script_arg))
source(file.path(script_dir, "config.R"))
source(file.path(script_dir, "analysis_lib.R"))
source(file.path(script_dir, "calibration_lib.R"))
source(file.path(script_dir, "systems_contract.R"))
source(file.path(script_dir, "t3_vocabulary_lib.R"))

# The contract output path. Overridable only so a reduced smoke run can be sent
# somewhere disposable; the sealed release always uses the default.
system_root <- Sys.getenv("T3_SYSTEM_ROOT",
                          file.path(workstream_root, "systems", "T3_vocabulary"))
n_random <- as.integer(Sys.getenv("T3_RANDOM_DRAWS", "200"))
n_calibration <- as.integer(Sys.getenv("T3_CALIBRATION_REPS", "100"))
mc_cores <- as.integer(Sys.getenv("T3_CORES", "16"))
t3_seed <- 20260812L

assert_inclusion_criterion("evidence_interpretation")
assert_language(readLines(file.path(script_dir, "t3_run_vocabulary.R")))
assert_language(readLines(file.path(script_dir, "t3_vocabulary_lib.R")))
if (dir.exists(system_root)) fail("Refusing to overwrite T3 release: ", system_root)
tmp <- atomic_dir(system_root)
on.exit(if (dir.exists(tmp)) unlink(tmp, recursive = TRUE), add = TRUE)
set.seed(t3_seed)

# Six published panels, two more than config.R's default set. Hoang_2019 and
# Suppli_2019 are added deliberately so the leakage table can show the full
# range: three panels derived from a discovery cohort, one from a transport
# cohort, and two from cohorts this Resource never uses.
t3_published_panels <- c(
  Govaere_2020 = "govaere_2020_panel.tsv",
  Pantano_2021 = "pantano_2021_panel.tsv",
  Hoang_2019   = "hoang_2019_panel.tsv",
  Suppli_2019  = "suppli_2019_panel.tsv",
  Moylan_2014  = "moylan_2014_panel.tsv",
  Arendt_2015  = "arendt_2015_panel.tsv"
)

# Vocabulary provenance, verified from the scripts that built each artifact.
# `substrate_leak` is TRUE when the vocabulary was learned on the same bulk
# matrix this system tests it on, which inflates its hit rate by construction.
vocabulary_provenance <- data.table(
  vocabulary = c("hotspot_117", "hallmark_50", "cnmf_16", "bnmf_6", "wgcna_4",
                 "published_panels"),
  derivation = c(
    "single-cell Hotspot per-cell-type modules (frozen 117, program_membership_v2)",
    "MSigDB Hallmark, SHA-pinned GMT",
    "single-cell cNMF k=16 global programs",
    "bulk RNA-seq NMF k=6 basis, top 120 genes per program",
    "WGCNA modules on the human bulk co-expression network (cross-species preservation arm)",
    "published MASLD/MASH gene panels"),
  substrate_leak = c(FALSE, FALSE, FALSE, TRUE, TRUE, NA),
  substrate_leak_note = c(
    "single-cell derived; no overlap with the bulk discovery matrix",
    "curated; no expression data involved",
    "single-cell derived; no overlap with the bulk discovery matrix",
    "learned on RNA-seq/results/subtypes/nmf_shard_clean_k6.rds, the bulk subtype NMF fit; the script that produced that fit is not in this repository, so the substrate overlap is inferred from its location in the bulk subtype results rather than read off its source",
    "network built on RNA-seq/.../integration/results/integration/merged_dge.rds, which contains all four discovery cohorts (Analysis/Cross_Species_Concordance/scripts/03a_wgcna_module_preservation.R:38)",
    "see panel_leakage.tsv, one row per panel")
)

# Cohort of origin for each panel. Discovery-derived panels are reported
# descriptively and never enter a support count.
panel_leakage <- data.table(
  feature_id = names(t3_published_panels),
  source_cohort = c("GSE135251", "GSE162694", "GSE130970", "GSE126848",
                    "GSE49541", "not_a_public_series"),
  leakage_class = c("discovery_cohort", "discovery_cohort", "discovery_cohort",
                    "transport_cohort", "external_clean", "external_clean")
)
panel_leakage[, reportable_as_support := leakage_class == "external_clean"]

message("[1/12] Loading sealed non-holdout input")
dge <- readRDS(nonholdout_dge)
meta <- as.data.table(readRDS(nonholdout_meta))
assert_true(ncol(dge) == 1097L, "Non-holdout DGE must contain 1,097 samples")
assert_true(identical(colnames(dge), meta$sample_id), "DGE and metadata sample order mismatch")
assert_holdout_sealed(meta, holdout_cohort, gse193066_crosswalk)
if (!"participant_id" %in% names(meta)) meta[, participant_id := sample_id]
# Non-vacuous: the unit column is a copy of sample_id, so uniqueness proves
# nothing on its own and the assertion is forced to verify externally that no
# discovery cohort ships a donor-pairing table.
assert_biological_unit(meta, "participant_id", cohorts = discovery_cohorts,
                       pairing_root = file.path(project_root, "data"))
assert_frozen_programs(program_registry, program_membership, expected_programs = 117L)

discovery_meta <- meta[dataset %in% discovery_cohorts &
                         !is.na(fibrosis_stage) & !is.na(nas_score)]
assert_true(nrow(discovery_meta) == 469L, "Discovery census drift")
assert_true(uniqueN(discovery_meta$participant_id) == 469L,
            "Discovery samples must be one per participant")

message("[2/12] Building symbol-level logCPM without outcome-based filtering")
logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
annotation <- fread(gene_annotation, select = c("gene_id", "gene_name"))
symbol_matrix <- collapse_symbols(logcpm, rownames(logcpm), annotation)
rm(logcpm, dge); gc()
ens2sym <- t3_ensembl_symbol_map(annotation)

# The universe every random set is drawn from and every score is computed on.
# Restricted to genes that are non-constant in all four discovery cohorts, so a
# random draw can never contain a gene that is unscoreable where it is used.
universe <- rownames(symbol_matrix)
for (co in discovery_cohorts) {
  j <- which(meta$dataset == co)
  sds <- matrixStats::rowSds(symbol_matrix[universe, j, drop = FALSE])
  universe <- universe[is.finite(sds) & sds > 0]
}
message("      universe: ", length(universe), " symbols of ", nrow(symbol_matrix))

message("[3/12] Replication anchor: reproducing the canonical weighted Hotspot count")
registry <- fread(program_registry)
membership <- fread(program_membership)
program_scores <- score_programs(symbol_matrix, meta, membership, registry,
                                 program_weight_coverage)
assert_true(sum(program_scores$coverage$testable) == expected_testable_programs,
            "Testable-program drift against the sealed discovery run")
S_weighted <- program_scores$primary[program_scores$coverage[testable == TRUE, feature_id], ,
                                     drop = FALSE]
anchor_cohort <- fit_feature_models(S_weighted, discovery_meta, discovery_cohorts, "primary")
anchor_meta <- meta_analyze(anchor_cohort)
anchor_meta[, q_value := p.adjust(p_value_meta, method = "BH", n = .N), by = model_kind]
anchor_counts <- anchor_meta[, .(n_supported = sum(q_value < 0.05, na.rm = TRUE)), by = axis]
print(anchor_counts)
assert_true(anchor_counts[axis == "fibrosis", n_supported] == 27L,
            paste0("Replication anchor failed: expected 27 fibrosis-supported programs, got ",
                   anchor_counts[axis == "fibrosis", n_supported]))
assert_true(anchor_counts[axis == "nas", n_supported] == 4L,
            "Replication anchor failed on the NAS axis")

designs <- t3_designs(discovery_meta, discovery_cohorts)
anchor_fast <- t3_fit_vocabulary(S_weighted[, discovery_meta$sample_id, drop = FALSE],
                                 designs, mc_cores = mc_cores)
anchor_check <- merge(anchor_fast, anchor_meta, by = c("feature_id", "axis"),
                      suffixes = c("", "_ref"))
assert_true(nrow(anchor_check) == nrow(anchor_fast) &&
              max(abs(anchor_check$p_value_meta - anchor_check$p_value_meta_ref)) < 1e-8 &&
              max(abs(anchor_check$beta_meta - anchor_check$beta_meta_ref)) < 1e-8,
            "Fast path disagrees with the canonical fit_feature_models + meta_analyze path")
write_tsv(anchor_check[, .(feature_id, axis, beta_meta, beta_meta_ref,
                           p_value_meta, p_value_meta_ref, q_value, q_value_ref)],
          file.path(tmp, "replication_anchor.tsv"))

message("[4/12] Assembling vocabularies")
vocab <- t3_build_vocabularies(program_membership, hallmark_gmt,
                               file.path(project_root,
                                         "RNA-seq/results/heterogeneity_program/phase0",
                                         "program_gene_membership.tsv"),
                               published_panel_root, t3_published_panels, ens2sym)
vocabularies <- unique(vocab$vocabulary)
message("      ", paste(sprintf("%s (%d programs)", vocabularies,
                                vapply(vocabularies, function(v)
                                  uniqueN(vocab[vocabulary == v, feature_id]), integer(1))),
                        collapse = " | "))

zlist <- t3_cohort_z(symbol_matrix, meta, discovery_cohorts, universe)
bins <- t3_expression_deciles(symbol_matrix)

message("[5/12] Scoring every vocabulary with the identical unweighted definition")
scored <- list()
coverage_rows <- list()
for (v in vocabularies) {
  mv <- vocab[vocabulary == v, .(feature_id, gene_symbol, direction)]
  sg <- score_gene_sets(symbol_matrix, meta, mv, gene_set_min_genes,
                        gene_set_coverage, directional = FALSE)
  sg$coverage[, vocabulary := v]
  coverage_rows[[v]] <- sg$coverage
  testable <- sg$coverage[testable == TRUE, feature_id]
  message("      ", v, ": ", length(testable), " of ",
          uniqueN(mv$feature_id), " testable")
  if (!length(testable)) {
    message("      ", v, ": no testable feature at coverage ", gene_set_coverage,
            "; vocabulary dropped")
    next
  }
  scored[[v]] <- sg$scores[testable, discovery_meta$sample_id, drop = FALSE]
}
dropped_vocabularies <- setdiff(vocabularies, names(scored))
if (length(dropped_vocabularies)) {
  message("      dropped: ", paste(dropped_vocabularies, collapse = ", "))
}
vocabularies <- names(scored)
coverage <- rbindlist(coverage_rows, use.names = TRUE)
write_tsv(coverage, file.path(tmp, "vocabulary_coverage.tsv"))

message("[6/12] Fitting the identical model in every vocabulary")
fits <- lapply(vocabularies, function(v) {
  f <- t3_fit_vocabulary(scored[[v]], designs, mc_cores = mc_cores)
  f[, vocabulary := v]
  f[, mde := minimum_detectable_effect(se_meta, n_cohorts)]
  f[, informative_null := informative_null(ci_lower, ci_upper,
                                           observability_effect_threshold)]
  f[, support := fcase(
    !is.finite(p_value_meta), "not_estimable",
    q_value < 0.05, "supported",
    informative_null == TRUE, "informative_null",
    default = "indeterminate")]
  f
})
names(fits) <- vocabularies
meta_effects <- rbindlist(fits, use.names = TRUE)
assert_negatives_have_mde(meta_effects, "support", "mde")
write_tsv(meta_effects, file.path(tmp, "vocabulary_meta_effects.tsv"))

# The weighted Hotspot arm, so the weighting contribution is readable directly
# rather than confounded into the vocabulary comparison.
anchor_fast[, `:=`(vocabulary = "hotspot_117_weighted",
                   mde = minimum_detectable_effect(se_meta, n_cohorts))]
anchor_fast[, informative_null := informative_null(ci_lower, ci_upper,
                                                   observability_effect_threshold)]
write_tsv(anchor_fast, file.path(tmp, "hotspot_weighted_meta_effects.tsv"))

message("[7/12] Calibration gate: permutation null per vocabulary")
calibration_rows <- list()
gate_failures <- character(0)
for (v in c(vocabularies, "hotspot_117_weighted")) {
  S <- if (v == "hotspot_117_weighted")
    S_weighted[, discovery_meta$sample_id, drop = FALSE] else scored[[v]]
  observed <- if (v == "hotspot_117_weighted") anchor_fast else fits[[v]]
  permuted_p_fn <- function(i) {
    set.seed(t3_seed + 1000L * which(c(vocabularies, "hotspot_117_weighted") == v) + i)
    md <- permute_histology_within_cohort(discovery_meta,
                                          c("fibrosis_stage", "nas_score"))
    d <- t3_designs(md, discovery_cohorts)
    t3_fit_vocabulary(S, d, mc_cores = mc_cores)$p_value_meta
  }
  cal <- calibrate_view(observed$p_value_meta, permuted_p_fn, n_reps = n_calibration,
                        label = paste0(v, "::linear_primary::both_axes"),
                        observed_p_is_empirical = FALSE)
  g <- gate_view(cal, strict = FALSE)
  if (!g$passed) gate_failures <- c(gate_failures, paste0(v, ": ", paste(g$reasons, collapse = "; ")))
  message("      ", v, ": ", sprintf("%.2f", cal$null_call_mean),
          " mean false calls per ", cal$n_tests, " tests (max ", cal$null_call_max,
          "), reportable = ", cal$reportable)
  calibration_rows[[v]] <- calibration_row(cal)
}
calibration <- rbindlist(calibration_rows, use.names = TRUE)
write_tsv(calibration, file.path(tmp, "calibration.tsv"))

message("[8/12] Matched-random null: the denominator")
random_rows <- list()
random_support_sets <- list()
per_program_rows <- list()
for (v in vocabularies) {
  testable <- rownames(scored[[v]])
  members <- split(vocab[vocabulary == v & feature_id %in% testable, gene_symbol],
                   vocab[vocabulary == v & feature_id %in% testable, feature_id])
  members <- lapply(members, function(g) intersect(g, universe))
  draws <- t3_matched_random_sets(members, universe, bins, n_random,
                                  seed = t3_seed + 100000L * match(v, vocabularies))
  flat <- unlist(draws, recursive = FALSE, use.names = FALSE)
  assert_true(is.list(flat) && length(flat) == length(draws) * n_random,
              "Flattening the matched-random draws did not preserve one set per element")
  names(flat) <- paste(rep(names(draws), each = n_random), rep(seq_len(n_random),
                                                              times = length(draws)),
                       sep = "||")
  message("      ", v, ": scoring ", length(flat), " matched-random sets")
  S_rand <- t3_score_sets_fast(zlist, flat)[, discovery_meta$sample_id, drop = FALSE]
  fit_rand <- t3_fit_vocabulary(S_rand, designs, mc_cores = mc_cores)
  # BH is applied WITHIN each synthetic vocabulary, so the null carries the same
  # multiplicity structure as the observed count it is the denominator for.
  fit_rand[, replicate := as.integer(sub("^.*\\|\\|", "", feature_id))]
  fit_rand[, q_within := p.adjust(p_value_meta, method = "BH", n = .N), by = replicate]
  counts <- fit_rand[, .(n_supported = sum(q_within < 0.05, na.rm = TRUE),
                         n_tested = uniqueN(feature_id)), by = .(replicate, axis)]
  counts[, vocabulary := v]
  random_rows[[v]] <- counts

  # Per-program view of the same null: how often a program's OWN 200 matched
  # draws reach a p-value at least as small as the program did. This is the
  # statistic that says which specific programs carry signal beyond their size
  # and expression profile, rather than which vocabularies do on average.
  #
  # Its p-value is empirical, so it has a hard resolution floor of 1/(R+1). At
  # R = 200 that floor sits above the BH threshold a family of this size needs,
  # which is exactly the resolution failure calibration_lib.R was written to
  # catch. The per-program p-values are therefore reported as they are and NO
  # BH-thresholded count is derived from them; required_permutations() records
  # what it would take.
  fit_rand[, base_feature := sub("\\|\\|.*$", "", feature_id)]
  own <- merge(fit_rand[, .(base_feature, axis, p_rand = p_value_meta)],
               fits[[v]][, .(base_feature = feature_id, axis, p_obs = p_value_meta)],
               by = c("base_feature", "axis"))
  pp <- own[, .(n_random = .N,
                p_vs_own_random = (1 + sum(p_rand <= p_obs, na.rm = TRUE)) / (.N + 1)),
            by = .(base_feature, axis)]
  pp[, `:=`(vocabulary = v,
            resolution_floor = 1 / (n_random + 1),
            bh_threshold_needed = 0.05 / (2 * length(members)),
            permutations_required_for_a_count = required_permutations(2 * length(members)))]
  pp[, count_reportable := resolution_floor < bh_threshold_needed]
  per_program_rows[[v]] <- pp
  random_support_sets[[v]] <- lapply(seq_len(n_random), function(r) {
    ids <- fit_rand[replicate == r & axis == "fibrosis" & q_within < 0.05, feature_id]
    if (!length(ids)) return(character(0))
    unique(universe[unlist(flat[ids], use.names = FALSE)])
  })
  rm(S_rand, fit_rand); gc()
}
random_counts <- rbindlist(random_rows, use.names = TRUE)
write_tsv(random_counts, file.path(tmp, "matched_random_counts.tsv"))
per_program_random <- rbindlist(per_program_rows, use.names = TRUE)
setnames(per_program_random, "base_feature", "feature_id")
write_tsv(per_program_random, file.path(tmp, "per_program_vs_matched_random.tsv"))

message("[9/12] Hit rate against the matched-random rate")
hit_rows <- list()
for (v in vocabularies) {
  K <- nrow(scored[[v]])
  for (a in c("fibrosis", "nas")) {
    observed <- fits[[v]][axis == a, sum(q_value < 0.05, na.rm = TRUE)]
    null_counts <- random_counts[vocabulary == v & axis == a, n_supported]
    ci <- t3_wilson(observed, K)
    hit_rows[[length(hit_rows) + 1L]] <- data.table(
      vocabulary = v, axis = a, n_testable = K,
      n_supported = observed,
      hit_rate = observed / K,
      hit_rate_ci_lower = ci$lower, hit_rate_ci_upper = ci$upper,
      random_hit_rate_mean = mean(null_counts) / K,
      random_hit_rate_q025 = as.numeric(stats::quantile(null_counts, 0.025)) / K,
      random_hit_rate_q975 = as.numeric(stats::quantile(null_counts, 0.975)) / K,
      random_n_supported_mean = mean(null_counts),
      random_n_supported_max = max(null_counts),
      n_random_draws = n_random,
      enrichment_over_random = if (mean(null_counts) > 0) observed / mean(null_counts) else NA_real_,
      p_vs_matched_random = t3_empirical_p_ge(observed, null_counts),
      median_mde = median(fits[[v]][axis == a, mde], na.rm = TRUE),
      n_informative_null = fits[[v]][axis == a, sum(support == "informative_null")],
      n_indeterminate = fits[[v]][axis == a, sum(support == "indeterminate")],
      # Nominal, and deliberately not a support count: the per-program empirical
      # p cannot resolve a BH threshold at 200 draws (see
      # per_program_vs_matched_random.tsv, count_reportable).
      n_nominal_beats_own_random = per_program_random[
        vocabulary == v & axis == a, sum(p_vs_own_random < 0.05)],
      per_program_count_reportable = per_program_random[
        vocabulary == v & axis == a, all(count_reportable)]
    )
  }
}
hit_rates <- rbindlist(hit_rows)
hit_rates <- merge(hit_rates, vocabulary_provenance[, .(vocabulary, substrate_leak)],
                   by = "vocabulary", all.x = TRUE)
# Three of the six published panels were derived from a discovery cohort, so
# their row is descriptive and must never be read as an independent support
# count. The two substrate-leaked vocabularies keep their counts but are barred
# from the headline comparison for the same reason.
hit_rates[, reportable_as_support := vocabulary != "published_panels" &
            !(substrate_leak %in% TRUE)]
hit_rates <- merge(hit_rates, calibration[, .(vocabulary = sub("::.*$", "", view),
                                              null_call_mean, reportable)],
                   by = "vocabulary", all.x = TRUE)
# The gate, applied as a gate rather than as an abort. A vocabulary whose
# permutation null invents calls keeps its row so the reason travels with it,
# but its count is blanked and it is excluded from every comparison below.
hit_rates[reportable == FALSE, `:=`(
  n_supported = NA_integer_, hit_rate = NA_real_,
  hit_rate_ci_lower = NA_real_, hit_rate_ci_upper = NA_real_,
  enrichment_over_random = NA_real_, p_vs_matched_random = NA_real_,
  withheld_reason = "calibration gate failed; see calibration.tsv")]
assert_counts_calibrated(hit_rates[reportable == TRUE],
                         calibration[view %in% paste0(hit_rates[reportable == TRUE, vocabulary],
                                                      "::linear_primary::both_axes")])
write_tsv(hit_rates, file.path(tmp, "hit_rates.tsv"))
print(hit_rates[, .(vocabulary, axis, n_testable, n_supported, hit_rate,
                    random_hit_rate_mean, enrichment_over_random, p_vs_matched_random)])

message("[10/12] Cross-vocabulary recovery")
# Recovery is counted only across the five vocabularies that are genuine
# independent clusterings. The published panels are not a clustering and three
# of the six are derived from a discovery cohort, so they are excluded here and
# reported descriptively in panel_leakage.tsv.
recovery_vocabularies <- setdiff(vocabularies, "published_panels")
supported_genes <- lapply(recovery_vocabularies, function(v) {
  ids <- fits[[v]][axis == "fibrosis" & q_value < 0.05, feature_id]
  unique(intersect(vocab[vocabulary == v & feature_id %in% ids, gene_symbol], universe))
})
names(supported_genes) <- recovery_vocabularies
gene_recovery <- data.table(gene_symbol = universe)
for (v in recovery_vocabularies) gene_recovery[, (v) := gene_symbol %in% supported_genes[[v]]]
gene_recovery[, n_vocabularies := rowSums(.SD), .SDcols = recovery_vocabularies]
external <- setdiff(recovery_vocabularies,
                    vocabulary_provenance[substrate_leak %in% TRUE, vocabulary])
gene_recovery[, n_external_vocabularies := rowSums(.SD), .SDcols = external]
write_tsv(gene_recovery[n_vocabularies > 0], file.path(tmp, "gene_recovery.tsv"))

# The same statistic under the matched-random null, so "recovered by 3 or more"
# has a denominator too.
# random_support_sets[[v]][[r]] is already unique within vocabulary v, so a
# tabulation across vocabularies counts vocabularies per gene, not sets per gene.
null_recovery <- vapply(seq_len(n_random), function(r) {
  tab <- table(unlist(lapply(recovery_vocabularies,
                             function(v) random_support_sets[[v]][[r]])))
  sum(tab >= 3L)
}, integer(1))
recovery_summary <- data.table(
  threshold_vocabularies = 3L,
  n_genes_observed = gene_recovery[n_vocabularies >= 3L, .N],
  n_genes_external_only = gene_recovery[n_external_vocabularies >= 3L, .N],
  n_external_vocabularies = length(external),
  random_mean = mean(null_recovery),
  random_max = max(null_recovery),
  p_vs_matched_random = t3_empirical_p_ge(gene_recovery[n_vocabularies >= 3L, .N],
                                          null_recovery)
)
write_tsv(recovery_summary, file.path(tmp, "recovery_summary.tsv"))
print(recovery_summary)

# Program-level correspondence: which supported programs in different
# vocabularies are the same biology, by gene overlap against the hypergeometric.
supported_members <- rbindlist(lapply(recovery_vocabularies, function(v) {
  ids <- fits[[v]][axis == "fibrosis" & q_value < 0.05, feature_id]
  vocab[vocabulary == v & feature_id %in% ids & gene_symbol %in% universe,
        .(vocabulary, feature_id, gene_symbol)]
}))
sets <- split(supported_members$gene_symbol,
              paste(supported_members$vocabulary, supported_members$feature_id, sep = "::"))
pair_rows <- list()
keys <- names(sets)
if (length(keys) >= 2L) {
  for (i in seq_len(length(keys) - 1L)) for (j in (i + 1L):length(keys)) {
    vi <- sub("::.*$", "", keys[[i]]); vj <- sub("::.*$", "", keys[[j]])
    if (vi == vj) next
    a <- sets[[i]]; b <- sets[[j]]
    ov <- length(intersect(a, b))
    pair_rows[[length(pair_rows) + 1L]] <- data.table(
      program_a = keys[[i]], program_b = keys[[j]],
      vocabulary_a = vi, vocabulary_b = vj,
      n_a = length(a), n_b = length(b), n_overlap = ov,
      jaccard = ov / length(union(a, b)),
      p_hyper = stats::phyper(ov - 1L, length(a), length(universe) - length(a),
                              length(b), lower.tail = FALSE)
    )
  }
}
program_pairs <- if (length(pair_rows)) rbindlist(pair_rows) else
  data.table(program_a = character(0), program_b = character(0),
             vocabulary_a = character(0), vocabulary_b = character(0),
             n_a = integer(0), n_b = integer(0), n_overlap = integer(0),
             jaccard = numeric(0), p_hyper = numeric(0))
program_pairs[, q_hyper := if (.N) p.adjust(p_hyper, method = "BH") else numeric(0)]
if (nrow(program_pairs)) setorder(program_pairs, p_hyper)
write_tsv(program_pairs, file.path(tmp, "supported_program_pairs.tsv"))

message("[11/12] Leakage and provenance tables")
panel_report <- if ("published_panels" %in% vocabularies) {
  panel_fit <- fits[["published_panels"]][axis == "fibrosis"]
  merge(panel_leakage, panel_fit[, .(feature_id, beta_meta, se_meta,
                                     p_value_meta, q_value, mde, support)],
        by = "feature_id", all.x = TRUE)
} else {
  copy(panel_leakage)[, note := "vocabulary dropped: no panel cleared the coverage floor"]
}
setorder(panel_report, leakage_class, feature_id)
write_tsv(panel_report, file.path(tmp, "panel_leakage.tsv"))
write_tsv(vocabulary_provenance, file.path(tmp, "vocabulary_provenance.tsv"))

message("[12/12] Sealing")
inputs <- c(nonholdout_dge, nonholdout_meta, program_registry, program_membership,
            hallmark_gmt, gene_annotation,
            file.path(project_root, "RNA-seq/results/heterogeneity_program/phase0",
                      "program_gene_membership.tsv"),
            file.path(published_panel_root, unname(t3_published_panels)))
write_tsv(data.table(path = inputs, sha256 = vapply(inputs, sha256_file, character(1))),
          file.path(tmp, "input_manifest.tsv"))
writeLines(capture.output(sessionInfo()), file.path(tmp, "sessionInfo.txt"))

jsonlite::write_json(list(
  state = if (!length(gate_failures)) "T3_COMPLETE_HOLDOUT_UNOPENED" else "T3_COMPLETE_WITH_WITHHELD_VIEWS",
  release_id = release_id, workstream_id = workstream_id,
  system = "T3_vocabulary_invariance",
  inclusion_criterion = "evidence_interpretation",
  n_discovery = nrow(discovery_meta),
  n_discovery_participants = uniqueN(discovery_meta$participant_id),
  vocabularies = vocabularies,
  n_matched_random_draws = n_random,
  n_calibration_permutations = n_calibration,
  gate_failures = if (length(gate_failures)) gate_failures else "none",
  replication_anchor_fibrosis = 27L,
  replication_anchor_nas = 4L,
  seed = t3_seed,
  holdout_accessed = FALSE
), file.path(tmp, "T3_READY.json"), pretty = TRUE, auto_unbox = TRUE)

artifacts <- setdiff(list.files(tmp, recursive = TRUE), "artifact_manifest.tsv")
write_tsv(data.table(
  artifact = artifacts,
  sha256 = vapply(file.path(tmp, artifacts), sha256_file, character(1)),
  size_bytes = file.info(file.path(tmp, artifacts))$size
), file.path(tmp, "artifact_manifest.tsv"))

publish_dir(tmp, system_root)
Sys.chmod(list.files(system_root, recursive = TRUE, full.names = TRUE), mode = "0440")
cat("T3_VOCABULARY_COMPLETE\t", system_root, "\n", sep = "")
