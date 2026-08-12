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
discovery_ready <- file.path(discovery_root, "DISCOVERY_READY.json")
if (!file.exists(discovery_ready)) fail("Discovery must be sealed before holdout access")
ready <- jsonlite::read_json(discovery_ready, simplifyVector = TRUE)
assert_true(identical(ready$state, "DISCOVERY_COMPLETE_HOLDOUT_UNOPENED"),
            "Discovery READY state does not authorize holdout")
assert_true(identical(ready$contract_sha256, sha256_file(contract_file)),
            "Discovery contract hash mismatch")
if (!file.exists(holdout_dge) || !file.exists(holdout_meta) ||
    !file.exists(holdout_composition)) fail("Sealed holdout partition missing")
if (dir.exists(holdout_root)) fail("Refusing to overwrite holdout release: ", holdout_root)

tmp <- atomic_dir(holdout_root)
on.exit(if (dir.exists(tmp)) unlink(tmp, recursive = TRUE), add = TRUE)

message("[1/7] Opening locked GSE193066 expression after discovery seal")
dge <- readRDS(holdout_dge)
meta <- as.data.table(readRDS(holdout_meta))
composition <- fread(holdout_composition, check.names = FALSE)
assert_true(ncol(dge) == 160L, "Holdout DGE must contain 160 QC-passing samples")
assert_true(all(as.character(dge$samples$dataset) == holdout_cohort),
            "Holdout DGE contains a non-GSE193066 sample")
assert_true(all(meta$dataset == holdout_cohort), "Holdout metadata contains another cohort")
assert_true(identical(colnames(dge), meta$sample_id), "Holdout DGE/metadata order mismatch")
assert_true(identical(composition$sample_id, meta$sample_id),
            "Holdout composition/metadata order mismatch")

crosswalk <- fread(gse193066_crosswalk)
crosswalk <- crosswalk[pass_technical == TRUE & corrected_dge_present == TRUE]
crosswalk_index <- match(meta$sample_id, crosswalk$run_id)
assert_true(!anyNA(crosswalk_index), "Every holdout sample must map to the locked crosswalk")
meta[, `:=`(
  participant_id = crosswalk$participant_token[crosswalk_index],
  biopsy = crosswalk$biopsy[crosswalk_index],
  is_first_biopsy = crosswalk$is_first_biopsy[crosswalk_index],
  is_second_biopsy = crosswalk$is_second_biopsy[crosswalk_index]
)]

logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
annotation <- fread(gene_annotation, select = c("gene_id", "gene_name"))
symbol_matrix <- collapse_symbols(logcpm, rownames(logcpm), annotation)
rm(logcpm)
gc()

registry <- fread(program_registry)
membership <- fread(program_membership)
program_scores <- score_programs(symbol_matrix, meta, membership, registry, program_weight_coverage)
assert_true(sum(program_scores$coverage$testable) == expected_testable_programs,
            "Holdout program testability drifted from the frozen discovery census")
assert_true(all(program_scores$cohort_coverage[
  feature_id %in% program_scores$coverage[testable == TRUE, feature_id], cohort_testable
]), "A discovery-testable program lacks 80% nonconstant weight in the holdout")

hallmark_membership <- read_gmt(hallmark_gmt)
hallmark_scores <- score_gene_sets(
  symbol_matrix, meta, hallmark_membership,
  gene_set_min_genes, gene_set_coverage, directional = FALSE
)
published_membership <- read_published_panels(published_panel_root, published_panels)
published_scores <- score_gene_sets(
  symbol_matrix, meta, published_membership,
  gene_set_min_genes, gene_set_coverage, directional = TRUE
)
composition_scores <- build_composition_scores(
  composition, meta, fread(composition_testability), composition_pseudocount
)
assert_true(all(composition_scores$coverage[testable == TRUE, n_observed] == 160L),
            "Holdout composition availability drift")
view_specs <- list(
  list(view = "hotspot_program", definition = "weighted_mean_z", scores = program_scores$primary,
       coverage = program_scores$coverage),
  list(view = "hotspot_program", definition = "unweighted_mean_z", scores = program_scores$unweighted,
       coverage = program_scores$coverage),
  list(view = "hotspot_program", definition = "weighted_within_sample_rank", scores = program_scores$weighted_rank,
       coverage = program_scores$coverage),
  list(view = "hallmark", definition = "unweighted_mean_z", scores = hallmark_scores$scores,
       coverage = hallmark_scores$coverage),
  list(view = "published_signature", definition = "directional_mean_z", scores = published_scores$scores,
       coverage = published_scores$coverage),
  list(view = "composition", definition = "clr_1e-6", scores = composition_scores$scores,
       coverage = composition_scores$coverage)
)

message("[2/7] Reconstructing biological participant pairs")
pair_meta <- meta

pair_counts <- pair_meta[, .(
  n_first = sum(is_first_biopsy), n_second = sum(is_second_biopsy)
), by = participant_id]
eligible_ids <- pair_counts[n_first == 1L & n_second == 1L, participant_id]
pairs_long <- pair_meta[participant_id %in% eligible_ids]
first <- pairs_long[is_first_biopsy == TRUE, .(
  participant_id, sample_first = sample_id,
  fibrosis_first = fibrosis_stage, nas_first = nas_score
)]
second <- pairs_long[is_second_biopsy == TRUE, .(
  participant_id, sample_second = sample_id,
  fibrosis_second = fibrosis_stage, nas_second = nas_score
)]
pairs <- merge(first, second, by = "participant_id")
pairs <- pairs[complete.cases(pairs)]
pairs[, `:=`(
  delta_fibrosis = fibrosis_second - fibrosis_first,
  delta_nas = nas_second - nas_first
)]
pairs[, changed := delta_fibrosis != 0 | delta_nas != 0]
pairs[, discordant := delta_fibrosis * delta_nas < 0]
assert_true(nrow(pairs) == 54L, paste0("Expected 54 complete pairs; observed ", nrow(pairs)))
assert_true(sum(pairs$changed) == 46L, paste0("Expected 46 changed pairs; observed ", sum(pairs$changed)))
assert_true(sum(!pairs$changed) == 8L, "Expected eight unchanged pairs")
assert_true(sum(pairs$discordant) == 11L, paste0("Expected 11 discordant pairs; observed ", sum(pairs$discordant)))
write_tsv(pairs, file.path(tmp, "paired_participant_manifest.tsv"))

meta_effects <- fread(file.path(discovery_root, "meta_effects.tsv"))
program_map <- fread(file.path(discovery_root, "program_map.tsv"))

message("[3/7] Running locked paired-change validation")
paired_summaries <- list()
participant_rows <- list()
program_validation <- program_map[testable == TRUE, .(
  feature_id, estimable = FALSE, failure_reason = "not_evaluated",
  n_participants = NA_integer_, spearman = NA_real_, p_value = NA_real_
)]
nulls <- list()
bootstraps <- list()
for (v in view_specs) {
  beta <- meta_effects[
    view == v$view & score_definition == v$definition & model_kind == "primary" &
      axis %in% c("fibrosis", "nas") & estimable == TRUE
  ]
  beta_wide <- dcast(beta, feature_id ~ axis, value.var = "beta_meta")
  features <- intersect(beta_wide$feature_id, rownames(v$scores))
  if (v$view == "hotspot_program" && v$definition == "weighted_mean_z") {
    assert_true(length(features) == expected_testable_programs,
                "Primary paired feature census drift")
  }
  beta_wide <- beta_wide[match(features, feature_id)]
  S <- v$scores[features, , drop = FALSE]
  observed <- t(S[, pairs$sample_second, drop = FALSE] - S[, pairs$sample_first, drop = FALSE])
  predicted <- outer(pairs$delta_fibrosis, beta_wide$fibrosis) +
    outer(pairs$delta_nas, beta_wide$nas)
  colnames(observed) <- colnames(predicted) <- features
  rownames(observed) <- rownames(predicted) <- pairs$participant_id
  cosine <- vapply(seq_len(nrow(pairs)), function(i) {
    cosine_similarity(observed[i, ], predicted[i, ])
  }, numeric(1))
  cosine_fibrosis_only <- rep(NA_real_, nrow(pairs))
  cosine_nas_only <- rep(NA_real_, nrow(pairs))
  if (v$view == "hotspot_program" && v$definition == "weighted_mean_z") {
    predicted_fibrosis_only <- outer(pairs$delta_fibrosis, beta_wide$fibrosis)
    predicted_nas_only <- outer(pairs$delta_nas, beta_wide$nas)
    cosine_fibrosis_only <- vapply(seq_len(nrow(pairs)), function(i) {
      cosine_similarity(observed[i, ], predicted_fibrosis_only[i, ])
    }, numeric(1))
    cosine_nas_only <- vapply(seq_len(nrow(pairs)), function(i) {
      cosine_similarity(observed[i, ], predicted_nas_only[i, ])
    }, numeric(1))
  }
  observed_norm <- sqrt(rowSums(observed^2, na.rm = TRUE))
  predicted_norm <- sqrt(rowSums(predicted^2, na.rm = TRUE))
  observed_norm[rowSums(is.finite(observed)) == 0L] <- NA_real_
  predicted_norm[rowSums(is.finite(predicted)) == 0L] <- NA_real_
  participant_rows[[length(participant_rows) + 1L]] <- data.table(
    participant_id = pairs$participant_id,
    view = v$view,
    score_definition = v$definition,
    delta_fibrosis = pairs$delta_fibrosis,
    delta_nas = pairs$delta_nas,
    changed = pairs$changed,
    discordant = pairs$discordant,
    cosine = cosine,
    cosine_fibrosis_only = cosine_fibrosis_only,
    cosine_nas_only = cosine_nas_only,
    observed_delta_norm = observed_norm,
    predicted_delta_norm = predicted_norm
  )

  changed_idx <- which(pairs$changed & is.finite(cosine))
  observed_stat <- stats::median(cosine[changed_idx])
  boot <- replicate(n_bootstrap, {
    idx <- sample(changed_idx, length(changed_idx), replace = TRUE)
    stats::median(cosine[idx], na.rm = TRUE)
  })
  assert_true(all(is.finite(boot)),
              paste0("Paired bootstrap contains non-finite values for ",
                     v$view, "::", v$definition))
  perm <- replicate(n_permutations, {
    shuffled <- sample(changed_idx, length(changed_idx), replace = FALSE)
    pred_perm <- outer(pairs$delta_fibrosis[shuffled], beta_wide$fibrosis) +
      outer(pairs$delta_nas[shuffled], beta_wide$nas)
    stats::median(vapply(seq_along(changed_idx), function(k) {
      cosine_similarity(observed[changed_idx[k], ], pred_perm[k, ])
    }, numeric(1)), na.rm = TRUE)
  })
  assert_true(all(is.finite(perm)),
              paste0("Paired permutation null contains non-finite values for ",
                     v$view, "::", v$definition))
  nulls[[paste(v$view, v$definition, sep = "::")]] <- perm
  bootstraps[[paste(v$view, v$definition, sep = "::")]] <- boot
  paired_summaries[[length(paired_summaries) + 1L]] <- data.table(
    split = "paired_biopsy_changed",
    cohort = holdout_cohort,
    view = v$view,
    score_definition = v$definition,
    axis = "joint_fibrosis_nas",
    metric = "median_participant_cosine",
    estimate = observed_stat,
    ci_lower = stats::quantile(boot, 0.025, na.rm = TRUE),
    ci_upper = stats::quantile(boot, 0.975, na.rm = TRUE),
    empirical_p = empirical_p(observed_stat, perm, "greater"),
    n_participants = length(changed_idx),
    n_features = length(features)
  )

  if (v$view == "hotspot_program" && v$definition == "weighted_mean_z") {
    for (j in seq_along(features)) {
      x <- predicted[changed_idx, j]
      y <- observed[changed_idx, j]
      ok <- is.finite(x) & is.finite(y)
      row <- which(program_validation$feature_id == features[[j]])
      if (sum(ok) < 10L || stats::sd(x[ok]) == 0 || stats::sd(y[ok]) == 0) {
        program_validation[row, `:=`(
          failure_reason = "fewer_than_10_or_zero_variance",
          n_participants = sum(ok)
        )]
        next
      }
      test <- suppressWarnings(stats::cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
      program_validation[row, `:=`(
        estimable = TRUE, failure_reason = NA_character_,
        n_participants = sum(ok), spearman = unname(test$estimate),
        p_value = test$p.value
      )]
    }
  }
}
paired_summary <- rbindlist(paired_summaries, fill = TRUE)
participant_validation <- rbindlist(participant_rows, fill = TRUE)
assert_true(nrow(program_validation) == expected_testable_programs,
            "Paired program-validation family size drift")
program_validation[, q_value := p.adjust(p_value, method = "BH", n = .N)]
write_tsv(paired_summary, file.path(tmp, "paired_validation_summary.tsv"))
write_tsv(participant_validation, file.path(tmp, "paired_participant_validation.tsv"))
write_tsv(program_validation, file.path(tmp, "paired_program_validation.tsv"))
saveRDS(nulls, file.path(tmp, "paired_permutation_nulls.rds"), compress = "xz")
saveRDS(bootstraps, file.path(tmp, "paired_bootstrap_distributions.rds"), compress = "xz")

message("[4/7] Quantifying unchanged and discordant prespecified diagnostics")
primary_participants <- participant_validation[
  view == "hotspot_program" & score_definition == "weighted_mean_z"
]
diagnostics <- rbindlist(list(
  primary_participants[discordant == TRUE, .(
    split = "paired_discordant", metric = "median_participant_cosine",
    estimate = median(cosine, na.rm = TRUE), n_participants = .N
  )],
  primary_participants[discordant == TRUE, .(
    split = "paired_discordant_fibrosis_only", metric = "median_participant_cosine",
    estimate = median(cosine_fibrosis_only, na.rm = TRUE), n_participants = .N
  )],
  primary_participants[discordant == TRUE, .(
    split = "paired_discordant_nas_only", metric = "median_participant_cosine",
    estimate = median(cosine_nas_only, na.rm = TRUE), n_participants = .N
  )],
  primary_participants[changed == FALSE, .(
    split = "paired_unchanged", metric = "median_observed_delta_norm",
    estimate = median(observed_delta_norm, na.rm = TRUE), n_participants = .N
  )],
  primary_participants[changed == TRUE, .(
    split = "paired_changed", metric = "median_observed_delta_norm",
    estimate = median(observed_delta_norm, na.rm = TRUE), n_participants = .N
  )]
), fill = TRUE)
diagnostics[, `:=`(cohort = holdout_cohort, view = "hotspot_program",
                   score_definition = "weighted_mean_z")]
write_tsv(diagnostics, file.path(tmp, "paired_diagnostics.tsv"))

message("[5/7] Running phenotype-native transport tests")
discovery_scores <- fread(file.path(discovery_root, "feature_scores.tsv.gz"))
transport_specs <- list(
  list(cohort = "GSE240729", kind = "marginal_fibrosis", axis = "fibrosis_marginal",
       filter = function(d) d[!is.na(fibrosis_stage)]),
  list(cohort = "GSE126848", kind = "disease_control", axis = "disease_control",
       filter = function(d) { d <- d[group_binary %in% c("Control", "Disease")]; d[, group_binary := factor(group_binary, levels = c("Control", "Disease"))]; d }),
  list(cohort = "GSE126848", kind = "nash_nafl", axis = "nash_nafl",
       filter = function(d) { d <- d[diagnosis_harmonized %in% c("NAFL", "NASH")]; d[, diagnosis_binary := factor(diagnosis_harmonized, levels = c("NAFL", "NASH"))]; d }),
  list(cohort = "GSE167523", kind = "nash_nafl", axis = "nash_nafl",
       filter = function(d) { d <- d[diagnosis_harmonized %in% c("NAFL", "NASH")]; d[, diagnosis_binary := factor(diagnosis_harmonized, levels = c("NAFL", "NASH"))]; d }),
  list(cohort = "GSE213621", kind = "disease_control", axis = "disease_control",
       filter = function(d) { d <- d[group_binary %in% c("Control", "Disease")]; d[, group_binary := factor(group_binary, levels = c("Control", "Disease"))]; d })
)
transport_rows <- list()
transport_views <- Filter(function(v) {
  (v$view == "hotspot_program" && v$definition == "weighted_mean_z") ||
    v$view %in% c("hallmark", "published_signature", "composition")
}, view_specs)
for (v in transport_views) {
  long <- discovery_scores[view == v$view & score_definition == v$definition]
  S <- dcast(long, feature_id ~ sample_id, value.var = "score")
  feature_ids <- S$feature_id
  Y <- as.matrix(S[, -"feature_id"])
  rownames(Y) <- feature_ids
  for (spec in transport_specs) {
    md <- unique(long[dataset == spec$cohort, .(
      sample_id, dataset, sex_final, fibrosis_stage, nas_score,
      group_binary, diagnosis_harmonized
    )])
    md <- spec$filter(copy(md))
    reference <- meta_effects[
      view == v$view & score_definition == v$definition &
        model_kind == spec$kind & axis == spec$axis & estimable == TRUE,
      .(feature_id, beta_reference = beta_meta)
    ]
    analysis_features <- intersect(reference$feature_id, rownames(Y))
    if (length(analysis_features) < 3L) next
    assay_complete <- complete_feature_vector_samples(
      Y, analysis_features, md$sample_id
    )
    md <- md[assay_complete]
    if (nrow(md) < 10L) next
    target <- fit_feature_models(Y, md, spec$cohort, spec$kind)
    target <- target[axis == spec$axis, .(feature_id, beta_target = beta)]
    cmp <- merge(reference, target, by = "feature_id")
    cmp <- cmp[is.finite(beta_reference) & is.finite(beta_target)]
    if (nrow(cmp) < 3L) next
    observed_rho <- cor(cmp$beta_reference, cmp$beta_target,
                        method = "spearman", use = "complete.obs")
    permuted <- replicate(n_permutations, {
      pmd <- copy(md)
      perm <- sample.int(nrow(pmd))
      if (spec$kind == "marginal_fibrosis") pmd[, fibrosis_stage := fibrosis_stage[perm]]
      if (spec$kind == "disease_control") pmd[, group_binary := group_binary[perm]]
      if (spec$kind == "nash_nafl") pmd[, diagnosis_binary := diagnosis_binary[perm]]
      B <- fit_matrix_coefficients(
        Y[cmp$feature_id, pmd$sample_id, drop = FALSE], pmd, spec$kind
      )
      if (is.null(B)) {
        NA_real_
      } else {
        term <- switch(spec$kind,
                       marginal_fibrosis = "fibrosis_centered",
                       disease_control = "group_binaryDisease",
                       nash_nafl = "diagnosis_binaryNASH")
        cor(cmp$beta_reference, B[term, match(cmp$feature_id, colnames(B))],
            method = "spearman", use = "complete.obs")
      }
    })
    assert_true(all(is.finite(permuted)),
                paste0("Transport permutation null is incomplete for ",
                       v$view, "::", spec$cohort, "::", spec$axis))
    transport_rows[[length(transport_rows) + 1L]] <- data.table(
      split = "native_phenotype_transport", cohort = spec$cohort,
      view = v$view, score_definition = v$definition, axis = spec$axis,
      metric = "effect_vector_spearman", estimate = observed_rho,
      empirical_p = empirical_p(observed_rho, permuted, "greater"),
      n_participants = nrow(md), n_features = nrow(cmp)
    )
  }
}
transport <- rbindlist(transport_rows, fill = TRUE)
transport[, p_holm_within_view := p.adjust(empirical_p, method = "holm"), by = view]
write_tsv(transport, file.path(tmp, "transport_validation.tsv"))

discovery_validation <- fread(file.path(discovery_root, "validation_summary.tsv"))
if (!"score_definition" %in% names(discovery_validation)) {
  discovery_validation[, score_definition := "weighted_mean_z"]
}
discovery_validation[, `:=`(
  multiplicity_adjustment = fifelse(
    split == "leave_one_cohort_out_global", "Holm across two global axes", "none"
  ),
  adjusted_p = p_holm
)]
paired_for_summary <- copy(paired_summary)
paired_for_summary[, `:=`(
  multiplicity_adjustment = "none_primary_or_prespecified_view",
  adjusted_p = empirical_p,
  biological_n = n_participants
)]
transport_for_summary <- copy(transport)
transport_for_summary[, `:=`(
  multiplicity_adjustment = "Holm within feature view across native contrasts",
  adjusted_p = p_holm_within_view,
  biological_n = n_participants
)]
validation_summary <- rbindlist(
  list(discovery_validation, paired_for_summary, transport_for_summary),
  fill = TRUE, use.names = TRUE
)
write_tsv(validation_summary, file.path(tmp, "validation_summary.tsv"))

message("[6/7] Writing holdout scores and manifests")
holdout_score_long <- rbindlist(lapply(view_specs, function(v) {
  scores_long(v$scores, meta, v$view, v$definition, v$coverage)
}), fill = TRUE)
holdout_score_long <- merge(
  holdout_score_long,
  registry[, .(feature_id = program_uid, cell_type, module_name,
               robust_display, primary_selected, source_stability)],
  by = "feature_id", all.x = TRUE, sort = FALSE
)
holdout_score_long[, cohort := dataset]
holdout_score_long <- merge(
  holdout_score_long,
  program_scores$cohort_coverage,
  by = c("feature_id", "cohort"), all.x = TRUE, sort = FALSE
)
holdout_score_long[, coverage := fcoalesce(
  cohort_weight_coverage, weight_coverage, gene_coverage
)]
holdout_score_long[, `:=`(sample = sample_id, participant = participant_id)]
setcolorder(
  holdout_score_long,
  c("sample", "participant", "cohort", "biopsy", "view", "feature_id",
    "score_definition", "coverage", "score", "testable")
)
fwrite(holdout_score_long, file.path(tmp, "holdout_feature_scores.tsv.gz"),
       sep = "\t", quote = FALSE, na = "NA")
all_feature_scores <- rbindlist(
  list(fread(file.path(discovery_root, "feature_scores.tsv.gz")), holdout_score_long),
  use.names = TRUE, fill = TRUE
)
fwrite(all_feature_scores, file.path(tmp, "feature_scores.tsv.gz"),
       sep = "\t", quote = FALSE, na = "NA")

write_tsv(data.table(
  input = c("contract", "discovery_ready", "holdout_dge", "holdout_meta",
            "holdout_composition", "composition_testability"),
  path = c(contract_file, discovery_ready, holdout_dge, holdout_meta,
           holdout_composition, composition_testability),
  sha256 = vapply(c(contract_file, discovery_ready, holdout_dge, holdout_meta,
                    holdout_composition, composition_testability),
                  sha256_file, character(1))
), file.path(tmp, "input_manifest.tsv"))
writeLines(capture.output(sessionInfo()), file.path(tmp, "sessionInfo.txt"))

artifacts <- list.files(tmp, recursive = TRUE, full.names = TRUE)
artifacts <- artifacts[file.info(artifacts)$isdir == FALSE]
write_tsv(data.table(
  artifact = substring(artifacts, nchar(tmp) + 2L),
  sha256 = vapply(artifacts, sha256_file, character(1)),
  size_bytes = file.info(artifacts)$size
), file.path(tmp, "artifact_manifest.tsv"))

primary_result <- paired_summary[
  view == "hotspot_program" & score_definition == "weighted_mean_z"
]
ready_payload <- list(
  state = "HOLDOUT_COMPLETE_NO_REFIT",
  release_id = release_id,
  workstream_id = workstream_id,
  contract_sha256 = sha256_file(contract_file),
  discovery_ready_sha256 = sha256_file(discovery_ready),
  n_pairs = nrow(pairs),
  n_changed_pairs = sum(pairs$changed),
  n_unchanged_pairs = sum(!pairs$changed),
  n_discordant_pairs = sum(pairs$discordant),
  paired_median_cosine = primary_result$estimate,
  paired_empirical_p = primary_result$empirical_p,
  paired_ci_lower = primary_result$ci_lower,
  paired_ci_upper = primary_result$ci_upper,
  discovery_refit = FALSE
)
jsonlite::write_json(ready_payload, file.path(tmp, "HOLDOUT_READY.json"),
                     pretty = TRUE, auto_unbox = TRUE)

message("[7/7] Publishing immutable holdout release")
publish_dir(tmp, holdout_root)
Sys.chmod(list.files(holdout_root, recursive = TRUE, full.names = TRUE), mode = "0440")
message("Published holdout: ", holdout_root)
