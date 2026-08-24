#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(splines)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_dir <- dirname(normalizePath(script_arg))
source(file.path(script_dir, "lib_continuum.R"))

pre <- read_prespec()
assert_true(identical(pre$resampling$anchor_bootstrap_unit, "unstratified_donor"),
            "Anchor-bootstrap prespecification drift")
assert_true(identical(
  pre$windows$interval_rule,
  "lower_closed_upper_open_except_final_upper_closed"
), "Fixed-window interval prespecification drift")
root <- project_root()
out <- out_root()
analysis_out <- file.path(out, "programs")
ensure_new_dir(analysis_out)

path_or_default <- function(env_name, default) read_env_path(env_name, default)
dge_path <- path_or_default(
  "HAC_DGE_PATH",
  file.path(root, "RNA-seq/results/manuscript_release/candidates",
            "resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BG001-DECISION",
            "arms/F_five/results/integration/merged_dge.rds")
)
manifest_path <- path_or_default(
  "HAC_MANIFEST_PATH",
  file.path(root, "figures/candidates/pi-figure-redesign-2026-08-13-v3/analysis",
            "stage_extensions/five_cohort_sample_manifest.tsv")
)
annotation_path <- path_or_default(
  "HAC_GENE_ANNOTATION_PATH",
  file.path(root, "RNA-seq/results/manuscript_release/candidates",
            "resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BULK-F-FIVE",
            "frozen_model_inputs/gencode_v49_gene_metadata.tsv.gz")
)
registry_path <- path_or_default(
  "HAC_PROGRAM_REGISTRY_PATH",
  file.path(root, "Analysis/Multimodal_Program_Projection/candidates",
            "program-context-v2-candidate-2026-08-07/hotspot/program_registry_v2.tsv")
)
membership_path <- path_or_default(
  "HAC_PROGRAM_MEMBERSHIP_PATH",
  file.path(root, "Analysis/Multimodal_Program_Projection/candidates",
            "program-context-v2-candidate-2026-08-07/hotspot/program_membership_v2.tsv")
)

unsupervised_path <- file.path(out, "unsupervised", "participant_scores.tsv")
projection_path <- file.path(out, "projection", "participant_scores.tsv")
signature_path <- file.path(out, "reproduction", "signature_genes.tsv")
for (path in c(unsupervised_path, projection_path, signature_path)) {
  assert_true(file.exists(path), paste0("Required upstream output missing: ", path))
}

message("Loading Resource expression and fixed program registry")
dge <- readRDS(dge_path)
meta <- fread(manifest_path, na.strings = c("", "NA"))
assert_true(setequal(colnames(dge$counts), meta$sample_id),
            "DGE/manifest participant sets differ")
assert_true(!anyDuplicated(meta$analysis_unit_id), "Cells or technical rows cannot be donor replicates")
logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
logcpm <- logcpm[, meta$sample_id, drop = FALSE]
assert_true(identical(colnames(logcpm), meta$sample_id),
            "Expression could not be reindexed to manifest order")
annotation <- fread(annotation_path, select = c("gene_id", "gene_name"))
symbol_matrix <- collapse_symbols(logcpm, rownames(logcpm), annotation)
rm(logcpm)
gc()

registry <- fread(registry_path)
membership_raw <- fread(membership_path)
signature <- fread(signature_path)
assert_true(nrow(registry) == pre$programs$family_size, "Program family is not 117")
assert_true(all(c("gene_id_base", "gene_symbol") %in% names(signature)),
            "Signature table schema drift")
annotation[, gene_id_base := base_gene_id(gene_id)]
assert_true(!anyDuplicated(annotation$gene_id_base),
            "GENCODE v49 annotation has duplicated Ensembl base IDs")
signature[, gene_id_base := base_gene_id(gene_id_base)]
signature_mapping <- merge(
  signature[, .(gene_id_base, released_gene_symbol = gene_symbol, in_resource)],
  annotation[, .(gene_id_base, gencode_v49_gene_symbol = gene_name)],
  by = "gene_id_base", all.x = TRUE, sort = FALSE
)
assert_true(nrow(signature_mapping) == pre$signature$published_gene_count &&
              !anyDuplicated(signature_mapping$gene_id_base),
            "The frozen 145-gene signature did not map one-to-one by Ensembl base ID")
assert_true(!anyNA(signature_mapping$gencode_v49_gene_symbol),
            "A signature Ensembl base ID is absent from the GENCODE v49 annotation")
signature_mapping[, symbol_concordant :=
                    toupper(released_gene_symbol) == toupper(gencode_v49_gene_symbol)]
assert_true(all(signature_mapping$symbol_concordant),
            "Released and GENCODE v49 signature symbols disagree")
signature_symbols <- unique(toupper(signature_mapping$gencode_v49_gene_symbol))
write_tsv_once(signature_mapping,
               file.path(analysis_out, "signature_gencode_v49_mapping.tsv"))

total_weight <- membership_raw[, .(
  total_l1_weight = sum(as.numeric(original_l1_weight)),
  n_source_rows = .N,
  n_mapped_original_genes = uniqueN(mapped_symbol[!is.na(mapped_symbol)]),
  unmapped_l1_weight = sum(as.numeric(original_l1_weight)[is.na(mapped_symbol)])
), by = program_uid]
membership <- membership_raw[
  !is.na(mapped_symbol) & mapped_symbol != "",
  .(weight = sum(as.numeric(original_l1_weight))),
  by = .(program_uid, gene_symbol = toupper(mapped_symbol))
]
membership[, excluded_signature_gene := gene_symbol %in% signature_symbols]
excluded_audit <- membership[excluded_signature_gene == TRUE, .(
  n_excluded_signature_genes = uniqueN(gene_symbol),
  excluded_l1_weight = sum(weight),
  excluded_genes = paste(sort(unique(gene_symbol)), collapse = ";")
), by = program_uid]
membership_use <- membership[excluded_signature_gene == FALSE]

coverage_rows <- list()
score_matrix <- matrix(
  NA_real_, nrow = nrow(registry), ncol = ncol(symbol_matrix),
  dimnames = list(registry$program_uid, colnames(symbol_matrix))
)

for (cohort in unique(meta$dataset)) {
  j <- which(meta$dataset == cohort)
  z <- zscore_rows(symbol_matrix[, j, drop = FALSE])
  variable <- rownames(z)[rowSums(is.finite(z)) == ncol(z)]
  for (program in registry$program_uid) {
    mm <- membership_use[program_uid == program & gene_symbol %in% variable]
    total <- total_weight[program_uid == program, total_l1_weight]
    excluded_signature <- excluded_audit[
      program_uid == program, excluded_l1_weight
    ]
    if (!length(excluded_signature)) excluded_signature <- 0
    signature_retained <- total - excluded_signature
    observed <- sum(mm$weight)
    n_observed <- uniqueN(mm$gene_symbol)
    retained_fraction <- signature_retained / total
    observed_fraction <- observed / total
    testable <- is.finite(retained_fraction) &&
      retained_fraction >= pre$programs$minimum_retained_l1_weight &&
      n_observed >= pre$programs$minimum_observed_genes
    coverage_rows[[length(coverage_rows) + 1L]] <- data.table(
      dataset = cohort,
      program_uid = program,
      n_observed_genes = n_observed,
      retained_l1_weight = signature_retained,
      retained_l1_fraction = retained_fraction,
      observed_l1_weight = observed,
      observed_l1_fraction = observed_fraction,
      testability_l1_basis = "one_minus_signature_excluded_original_l1",
      scoring_weight_rule = "original_relative_weights_renormalized_over_observed_genes",
      testable = testable
    )
    if (!testable) next
    idx <- match(mm$gene_symbol, rownames(z))
    assert_true(observed > 0, paste0("Testable program has zero observed weight: ", program))
    w <- mm$weight / observed
    score_matrix[program, j] <- as.numeric(crossprod(w, z[idx, , drop = FALSE]))
  }
}

coverage <- rbindlist(coverage_rows)
coverage <- merge(coverage, total_weight, by = "program_uid", all.x = TRUE)
coverage <- merge(coverage, excluded_audit, by = "program_uid", all.x = TRUE)
coverage[is.na(n_excluded_signature_genes), `:=`(
  n_excluded_signature_genes = 0L,
  excluded_l1_weight = 0,
  excluded_genes = ""
)]
coverage <- merge(
  coverage,
  registry[, .(program_uid, cell_type, module, module_name, robust_display)],
  by = "program_uid", all.x = TRUE
)
write_tsv_once(coverage, file.path(analysis_out, "program_testability.tsv"))

focal_ids <- unname(unlist(pre$programs$focal_programs, use.names = FALSE))
focal_cov <- coverage[program_uid %in% focal_ids]
assert_true(nrow(focal_cov) == length(focal_ids) * uniqueN(meta$dataset),
            "Focal-program coverage rows are incomplete")
expected_focal_retention <- c(
  hotspot_hepatocytes_48f39dd4d817a10e = 0.957397465312792,
  hotspot_hepatocytes_f05c535ae5bbc0b9 = 0.804964557414518
)
observed_focal_retention <- focal_cov[, unique(retained_l1_fraction), by = program_uid]
assert_true(nrow(observed_focal_retention) == length(focal_ids) &&
              all(abs(observed_focal_retention$V1 -
                      expected_focal_retention[observed_focal_retention$program_uid]) < 1e-12),
            "Focal signature-exclusion L1 retention differs from the frozen 95.7%/80.5% audit")
assert_true(all(focal_cov$testable), "A focal program became untestable after signature exclusion")

for (cohort in unique(meta$dataset)) {
  j <- which(meta$dataset == cohort)
  for (i in seq_len(nrow(score_matrix))) {
    score_matrix[i, j] <- standardize_vector(score_matrix[i, j])
  }
}
program_scores <- as.data.table(as.table(score_matrix))
setnames(program_scores, c("program_uid", "sample_id", "program_score"))
program_scores[, `:=`(
  program_uid = as.character(program_uid),
  sample_id = as.character(sample_id),
  program_score = as.numeric(program_score)
)]
program_scores <- merge(program_scores, meta, by = "sample_id", all.x = TRUE, sort = FALSE)
write_tsv_once(program_scores[, .(
  sample_id, analysis_unit_id, dataset, inferred_sex, fibrosis_stage, nas_score,
  program_uid, program_score
)], file.path(analysis_out, "program_scores.tsv.gz"))

message("Loading axes and enforcing method roles")
axes <- rbindlist(list(
  fread(unsupervised_path),
  fread(projection_path)
), use.names = TRUE, fill = TRUE)
required_axis <- c("sample_id", "dataset", "axis_id", "axis_raw", "axis_percentile")
assert_true(all(required_axis %in% names(axes)), "Participant-axis schema drift")
axes <- unique(axes[, ..required_axis], by = c("sample_id", "axis_id"))
assert_true(!anyDuplicated(axes[, .(sample_id, axis_id)]), "Duplicate participant-axis rows")
expected_axes <- c(pre$axes$co_primary, pre$axes$secondary, pre$axes$display_only)
assert_true(all(expected_axes %in% axes$axis_id),
            paste0("Missing expected axes: ", paste(setdiff(expected_axes, axes$axis_id), collapse = ",")))

independent <- as.character(pre$cohorts$independent_primary)
model_axes <- c(pre$axes$co_primary, pre$axes$secondary, pre$axes$display_only)

fit_one <- function(d, stage_column) {
  d <- copy(d)
  d <- d[is.finite(program_score) & is.finite(axis_raw) &
           !is.na(get(stage_column)) & !is.na(inferred_sex)]
  if (nrow(d) < 20L || uniqueN(d[[stage_column]]) < 2L || uniqueN(d$inferred_sex) < 2L) {
    return(data.table(
      estimable = FALSE, failure_reason = "insufficient_complete_design",
      n_donors = nrow(d), beta = NA_real_, se = NA_real_, ci_low = NA_real_,
      ci_high = NA_real_, t_value = NA_real_, p_value = NA_real_,
      nonlinear_p_value = NA_real_
    ))
  }
  d[, axis_z := standardize_vector(axis_raw)]
  d[, stage_factor := factor(get(stage_column))]
  d[, sex_factor := factor(inferred_sex)]
  linear <- tryCatch(lm(program_score ~ axis_z + stage_factor + sex_factor, data = d),
                     error = function(e) e)
  if (inherits(linear, "error") || "axis_z" %notin% rownames(coef(summary(linear)))) {
    reason <- if (inherits(linear, "error")) conditionMessage(linear) else "axis_not_estimable"
    return(data.table(
      estimable = FALSE, failure_reason = reason, n_donors = nrow(d),
      beta = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
      t_value = NA_real_, p_value = NA_real_, nonlinear_p_value = NA_real_
    ))
  }
  cs <- coef(summary(linear))["axis_z", ]
  nonlinear <- tryCatch(
    lm(program_score ~ ns(axis_z, df = 3) + stage_factor + sex_factor, data = d),
    error = function(e) e
  )
  nonlinear_p <- NA_real_
  if (!inherits(nonlinear, "error")) {
    comparison <- tryCatch(anova(linear, nonlinear), error = function(e) NULL)
    if (!is.null(comparison) && nrow(comparison) == 2L) nonlinear_p <- comparison$`Pr(>F)`[[2L]]
  }
  data.table(
    estimable = TRUE, failure_reason = NA_character_, n_donors = nrow(d),
    beta = unname(cs[["Estimate"]]), se = unname(cs[["Std. Error"]]),
    ci_low = unname(cs[["Estimate"]] - qnorm(0.975) * cs[["Std. Error"]]),
    ci_high = unname(cs[["Estimate"]] + qnorm(0.975) * cs[["Std. Error"]]),
    t_value = unname(cs[["t value"]]), p_value = unname(cs[["Pr(>|t|)"]]),
    nonlinear_p_value = nonlinear_p
  )
}

fit_rows <- list()
nas_rows <- list()
for (axis_name in model_axes) {
  for (cohort in independent) {
    axis_part <- axes[axis_id == axis_name & dataset == cohort]
    joined_axis <- merge(
      program_scores[dataset == cohort],
      axis_part[, .(sample_id, axis_raw, axis_percentile)],
      by = "sample_id", all = FALSE
    )
    assert_true(uniqueN(joined_axis$analysis_unit_id) == uniqueN(axis_part$sample_id),
                paste0("Axis/program donor join drift for ", axis_name, "/", cohort))
    for (program in registry$program_uid) {
      result <- fit_one(joined_axis[program_uid == program], "fibrosis_stage")
      result[, `:=`(axis_id = axis_name, dataset = cohort, program_uid = program,
                    stage_adjustment = "factor(fibrosis_stage)")]
      fit_rows[[length(fit_rows) + 1L]] <- result
    }
    if (cohort == "GSE162694") {
      for (program in registry$program_uid) {
        result <- fit_one(joined_axis[program_uid == program], "nas_score")
        result[, `:=`(axis_id = axis_name, dataset = cohort, program_uid = program,
                      stage_adjustment = "factor(nas_score)")]
        nas_rows[[length(nas_rows) + 1L]] <- result
      }
    }
  }
}
fits <- rbindlist(fit_rows, fill = TRUE)
fits[, nonlinear_q_value := p_adjust_complete_family(
  nonlinear_p_value, pre$programs$family_size
), by = .(axis_id, dataset)]
fits <- merge(fits, registry[, .(program_uid, cell_type, module, module_name, robust_display)],
              by = "program_uid", all.x = TRUE)
write_tsv_once(fits, file.path(analysis_out, "program_fibrosis_models.tsv"))

nas_fits <- rbindlist(nas_rows, fill = TRUE)
nas_fits[, `:=`(
  q_value = p_adjust_complete_family(p_value, pre$programs$family_size),
  nonlinear_q_value = p_adjust_complete_family(
    nonlinear_p_value, pre$programs$family_size
  )
), by = axis_id]
nas_fits <- merge(
  nas_fits,
  registry[, .(program_uid, cell_type, module, module_name, robust_display)],
  by = "program_uid", all.x = TRUE
)
write_tsv_once(nas_fits, file.path(analysis_out, "program_nas_sensitivity.tsv"))

meta_rows <- list()
for (axis_name in model_axes) {
  for (program in registry$program_uid) {
    block <- fits[axis_id == axis_name & program_uid == program & estimable == TRUE]
    meta_fit <- fixed_effect_meta(block$beta, block$se)
    row <- as.data.table(meta_fit)
    row[, `:=`(
      axis_id = axis_name,
      program_uid = program,
      direction_concordant = nrow(block) == length(independent) &&
        length(unique(sign(block$beta))) == 1L,
      cohort_betas = paste(block$dataset, signif(block$beta, 6), sep = "=", collapse = ";")
    )]
    meta_rows[[length(meta_rows) + 1L]] <- row
  }
}
meta_results <- rbindlist(meta_rows, fill = TRUE)
meta_results[, q_value := p_adjust_complete_family(
  p_value, pre$programs$family_size
), by = axis_id]
meta_results <- merge(
  meta_results,
  registry[, .(program_uid, cell_type, module, module_name, robust_display)],
  by = "program_uid", all.x = TRUE
)
write_tsv_once(meta_results, file.path(analysis_out, "program_meta_analysis.tsv"))

message("Calculating independent-cohort anchors and co-primary agreement")
anchor_rows <- list()
for (axis_name in model_axes) {
  for (cohort in independent) {
    d <- merge(
      axes[axis_id == axis_name & dataset == cohort],
      meta[dataset == cohort, .(sample_id, fibrosis_stage, nas_score)],
      by = "sample_id", all.x = TRUE
    )
    ok <- is.finite(d$axis_raw) & !is.na(d$fibrosis_stage)
    rho <- suppressWarnings(cor(d$axis_raw[ok], d$fibrosis_stage[ok], method = "spearman"))
    row <- data.table(
      axis_id = axis_name, dataset = cohort, anchor = "fibrosis_stage",
      n_donors = sum(ok), spearman_rho = rho,
      bootstrap_ci_low = NA_real_, bootstrap_ci_high = NA_real_,
      bootstrap_replicates = 0L, bootstrap_scheme = NA_character_
    )
    if (axis_name %in% pre$axes$co_primary) {
      boot <- bootstrap_spearman(
        d$axis_raw, d$fibrosis_stage, rep("all_donors", nrow(d)),
        pre$resampling$bootstrap_replicates,
        pre$seeds$bootstrap + match(axis_name, pre$axes$co_primary) * 100L +
          match(cohort, independent),
        pre$promotion_gates$anchor_bootstrap_confidence
      )
      row[, `:=`(
        spearman_rho = boot$rho,
        bootstrap_ci_low = boot$ci_low,
        bootstrap_ci_high = boot$ci_high,
        bootstrap_replicates = pre$resampling$bootstrap_replicates,
        bootstrap_scheme = "unstratified_donor"
      )]
    }
    anchor_rows[[length(anchor_rows) + 1L]] <- row
  }
}
anchors <- rbindlist(anchor_rows)
write_tsv_once(anchors, file.path(analysis_out, "axis_histology_anchors.tsv"))

agreement_rows <- lapply(independent, function(cohort) {
  wide <- dcast(
    axes[dataset == cohort & axis_id %in% pre$axes$co_primary,
         .(sample_id, axis_id, axis_raw)],
    sample_id ~ axis_id, value.var = "axis_raw"
  )
  a <- pre$axes$co_primary[[1L]]
  b <- pre$axes$co_primary[[2L]]
  ok <- is.finite(wide[[a]]) & is.finite(wide[[b]])
  data.table(
    dataset = cohort, axis_x = a, axis_y = b, n_donors = sum(ok),
    spearman_rho = suppressWarnings(cor(wide[[a]][ok], wide[[b]][ok], method = "spearman"))
  )
})
agreement <- rbindlist(agreement_rows)
write_tsv_once(agreement, file.path(analysis_out, "co_primary_agreement.tsv"))

all_axis_pairs <- combn(expected_axes, 2L, simplify = FALSE)
all_axis_agreement <- rbindlist(lapply(unique(meta$dataset), function(cohort) {
  rbindlist(lapply(all_axis_pairs, function(pair) {
    x <- axes[dataset == cohort & axis_id == pair[[1L]],
              .(sample_id, x = axis_percentile)]
    y <- axes[dataset == cohort & axis_id == pair[[2L]],
              .(sample_id, y = axis_percentile)]
    d <- merge(x, y, by = "sample_id", all = FALSE)
    data.table(
      dataset = cohort,
      evaluation_role = ifelse(cohort %in% independent,
                               "independent_validation", "source_overlap_replication"),
      axis_x = pair[[1L]],
      axis_y = pair[[2L]],
      n_donors = nrow(d),
      spearman_rho = suppressWarnings(cor(d$x, d$y, method = "spearman")),
      pearson_r = suppressWarnings(cor(d$x, d$y, method = "pearson"))
    )
  }))
}))
assert_true(nrow(all_axis_agreement) == uniqueN(meta$dataset) *
              choose(length(expected_axes), 2L),
            "All-scorer agreement table is incomplete")
write_tsv_once(all_axis_agreement,
               file.path(analysis_out, "all_axis_score_agreement.tsv"))

message("Building fixed visualization windows")
window_rows <- list()
centers <- as.numeric(pre$windows$centers)
half_width <- pre$windows$width / 2
for (axis_name in c(pre$axes$co_primary, pre$axes$display_only)) {
  for (cohort in independent) {
    axis_part <- axes[axis_id == axis_name & dataset == cohort]
    joined <- merge(
      program_scores[dataset == cohort],
      axis_part[, .(sample_id, axis_percentile)],
      by = "sample_id", all = FALSE
    )
    for (center in centers) {
      lower <- max(0, center - half_width)
      upper <- min(1, center + half_width)
      is_final_window <- isTRUE(all.equal(center, max(centers), tolerance = 0))
      block <- joined[
        is.finite(axis_percentile) & axis_percentile >= lower &
          (axis_percentile < upper | (is_final_window & axis_percentile <= upper)) &
          is.finite(program_score)
      ]
      summary <- block[, .(
        n_donors = .N,
        mean_score = mean(program_score),
        se_score = sd(program_score) / sqrt(.N)
      ), by = program_uid]
      summary[, `:=`(
        axis_id = axis_name, dataset = cohort, window_center = center,
        window_lower = lower,
        window_upper = upper,
        lower_inclusive = TRUE,
        upper_inclusive = is_final_window
      )]
      window_rows[[length(window_rows) + 1L]] <- summary
    }
  }
}
windows <- rbindlist(window_rows, fill = TRUE)
windows <- merge(windows, registry[, .(program_uid, module_name, robust_display)],
                 by = "program_uid", all.x = TRUE)
write_tsv_once(windows, file.path(analysis_out, "fixed_program_windows.tsv.gz"))

message("Running within-stage permutation controls for focal programs")
permutation_rows <- list()
for (axis_name in pre$axes$co_primary) {
  for (cohort in independent) {
    for (program in focal_ids) {
      d <- merge(
        program_scores[dataset == cohort & program_uid == program],
        axes[axis_id == axis_name & dataset == cohort,
             .(sample_id, axis_raw)],
        by = "sample_id", all = FALSE
      )
      d <- d[is.finite(program_score) & is.finite(axis_raw) &
               !is.na(fibrosis_stage) & !is.na(inferred_sex)]
      d[, `:=`(
        axis_z = standardize_vector(axis_raw),
        stage_factor = factor(fibrosis_stage),
        sex_factor = factor(inferred_sex)
      )]
      rx <- residuals(lm(axis_z ~ stage_factor + sex_factor, data = d))
      ry <- residuals(lm(program_score ~ stage_factor + sex_factor, data = d))
      observed <- suppressWarnings(cor(rx, ry))
      strata <- interaction(d$stage_factor, d$sex_factor, drop = TRUE)
      groups <- split(seq_along(rx), strata)
      set.seed(pre$seeds$permutation + match(axis_name, pre$axes$co_primary) * 1000L +
                 match(cohort, independent) * 100L + match(program, focal_ids))
      null <- replicate(pre$resampling$permutation_replicates, {
        rx_perm <- rx
        for (g in groups) rx_perm[g] <- sample(rx[g], length(g), replace = FALSE)
        suppressWarnings(cor(rx_perm, ry))
      })
      empirical_p <- (1 + sum(abs(null) >= abs(observed), na.rm = TRUE)) /
        (1 + sum(is.finite(null)))
      permutation_rows[[length(permutation_rows) + 1L]] <- data.table(
        axis_id = axis_name, dataset = cohort, program_uid = program,
        n_donors = nrow(d), residual_pearson = observed,
        permutation_p_value = empirical_p,
        permutation_replicates = pre$resampling$permutation_replicates,
        permutation_strata = "fibrosis_stage_by_sex"
      )
    }
  }
}
write_tsv_once(rbindlist(permutation_rows),
               file.path(analysis_out, "focal_within_stage_permutations.tsv"))

nmf_path <- read_env_path(
  "HAC_NMF_AXIS_PATH",
  file.path(root, "RNA-seq/results/continuous_axis_benchmark",
            "20260814T184338Z/arms/A3_axis.tsv")
)
if (file.exists(nmf_path)) {
  nmf_source <- fread(nmf_path)
  assert_true(nrow(nmf_source) == nrow(meta) &&
                setequal(nmf_source$sample_id, meta$sample_id),
              "Frozen NMF axis participant set drift")
  nmf <- nmf_source[, .(sample_id, dataset, nmf_k6 = axis_raw, nmf_k4 = axis_raw_k4)]
  native_scores <- rbindlist(list(
    nmf[, .(sample_id, dataset, axis_id = "nmf_k6", axis_raw = nmf_k6)],
    nmf[, .(sample_id, dataset, axis_id = "nmf_k4", axis_raw = nmf_k4)]
  ))
  native_scores[, axis_percentile := percent_rank(axis_raw), by = .(dataset, axis_id)]
  write_tsv_once(native_scores,
                 file.path(analysis_out, "native_nmf_participant_scores.tsv"))
  context_rows <- list()
  for (axis_name in c(pre$axes$co_primary, pre$axes$display_only)) {
    d <- merge(axes[axis_id == axis_name], nmf, by = c("sample_id", "dataset"),
               all = FALSE)
    for (native in c("nmf_k6", "nmf_k4")) {
      context_rows[[length(context_rows) + 1L]] <- d[, .(
        n_donors = sum(is.finite(axis_raw) & is.finite(get(native))),
        spearman_rho = suppressWarnings(cor(axis_raw, get(native), method = "spearman",
                                            use = "complete.obs"))
      ), by = dataset][, `:=`(axis_id = axis_name, native_axis = native)]
    }
  }
  write_tsv_once(rbindlist(context_rows),
                 file.path(analysis_out, "native_nmf_context.tsv"))
}

composition_path <- file.path(root, "RNA-seq/results/celltype_attribution",
                              "persample_celltype_proportions.csv")
composition_state <- data.table(
  analysis = "composition_adjustment",
  state = "deferred_pending_synchronized_five_cohort_composition_release",
  source_path = composition_path,
  source_sha256 = if (file.exists(composition_path)) sha256_file(composition_path) else NA_character_,
  reason = paste(
    "Composition is supplementary and cannot be mixed into the selected five-cohort",
    "candidate until the synchronized composition release is promoted."
  )
)
write_tsv_once(composition_state, file.path(analysis_out, "composition_sensitivity_state.tsv"))

write_session_info(file.path(analysis_out, "sessionInfo.txt"))
message("PROGRAM_VALIDATION_COMPLETE: ", analysis_out)
