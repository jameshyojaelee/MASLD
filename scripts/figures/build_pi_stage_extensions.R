#!/usr/bin/env Rscript
# KEY MESSAGE: One synchronized five-cohort fragment substrate supports every bulk stage, NAS, and F-versus-F0 source table used in the revised figures.

# Candidate-only extension of the validated fragment-count stage analysis.
# Fits direct F-versus-F0 and adjacent NAS models, and estimates centered
# prespecified NMF P1-P6 activity across fibrosis and NAS groups.

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
})

options(digits = 17, scipen = 999)
set.seed(42)

fail <- function(...) stop(..., call. = FALSE)
sha256_file <- function(path) {
  value <- system2("sha256sum", normalizePath(path, mustWork = TRUE), stdout = TRUE)
  strsplit(value[[1L]], "[[:space:]]+")[[1L]][[1L]]
}
write_tsv_once <- function(x, path) {
  if (file.exists(path)) fail("Refusing overwrite: ", path)
  fwrite(x, path, sep = "\t", quote = FALSE, na = "NA")
}

project_root <- normalizePath(Sys.getenv(
  "MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
), mustWork = TRUE)
candidate_root <- Sys.getenv("FIGURE_CANDIDATE_ROOT", "")
stage_root <- normalizePath(Sys.getenv("STAGE_RELEASE_ROOT", ""), mustWork = TRUE)
if (!nzchar(candidate_root)) fail("FIGURE_CANDIDATE_ROOT is required")
if (file.exists(candidate_root)) fail("FIGURE_CANDIDATE_ROOT must not exist: ", candidate_root)

stage_validation <- file.path(stage_root, "VALIDATED.json")
model_manifest_path <- file.path(stage_root, "model_input_manifest.tsv")
if (!file.exists(stage_validation) || !file.exists(model_manifest_path) ||
    !grepl('"status"[[:space:]]*:[[:space:]]*"VALIDATED"',
           paste(readLines(stage_validation, warn = FALSE), collapse = ""))) {
  fail("STAGE_RELEASE_ROOT must be the validated five-cohort fragment release")
}
model_manifest <- fread(model_manifest_path)
# A corrected-control refit ships expected_counts.json beside its manifest; the
# adopted 2026-08-10 release has none and keeps its 23,370-gene universe.
expected_counts_path <- file.path(stage_root, "expected_counts.json")
n_genes_expected <- if (file.exists(expected_counts_path)) {
  jsonlite::fromJSON(expected_counts_path)$n_genes
} else {
  23370L
}
unit_qualified <- "identity_qualified" %in% names(model_manifest) &&
  isTRUE(as.logical(model_manifest$identity_qualified[[1L]]))
n_samples_expected <- if (unit_qualified) as.integer(model_manifest$n_samples[[1L]]) else 844L
if (nrow(model_manifest) != 1L || model_manifest$n_cohorts != 5L ||
    model_manifest$n_samples != n_samples_expected || model_manifest$n_genes != n_genes_expected ||
    (unit_qualified && (!file.exists(expected_counts_path) ||
      jsonlite::fromJSON(expected_counts_path)$n_samples != n_samples_expected))) {
  fail("Five-cohort model-input contract drift")
}

release_candidate <- file.path(
  project_root,
  "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10"
)
dge_path <- normalizePath(model_manifest$dge_path, mustWork = TRUE)
meta_path <- normalizePath(model_manifest$meta_path, mustWork = TRUE)
annotation_path <- normalizePath(model_manifest$gene_metadata_path, mustWork = TRUE)
nmf_path <- file.path(project_root, "RNA-seq/results/subtypes/nmf_assignments.csv")
context_roster_path <- Sys.getenv("STAGE_CONTEXT_ROSTER", "")
if (unit_qualified) {
  nmf_path <- Sys.getenv("STAGE_FIXED_NMF_SCORES", "")
  if (!nzchar(context_roster_path) || !nzchar(nmf_path)) {
    fail("Qualified stage models require STAGE_CONTEXT_ROSTER and STAGE_FIXED_NMF_SCORES")
  }
  nmf_path <- normalizePath(nmf_path, mustWork = TRUE)
  context_roster_path <- normalizePath(context_roster_path, mustWork = TRUE)
}
program_label_path <- file.path(project_root, "RNA-seq/results/subtypes/program_labels.csv")

out_root <- file.path(candidate_root, "analysis", "stage_extensions")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(out_root)) fail("Cannot create output root")

dge_all <- readRDS(dge_path)
if (!inherits(dge_all, "DGEList") || !identical(dim(dge_all), c(as.integer(n_genes_expected), n_samples_expected))) {
  fail("Five-cohort corrected fragment DGE identity drift")
}
stage_annotation <- fread(annotation_path)[, .(
  gene_id_versioned = gene_id, gene_id_base = ensembl_base, gene_name,
  gene_type = gene_biotype, chromosome,
  is_canonical_chromosome = chromosome %in% c(paste0("chr", 1:22), "chrX", "chrY"),
  display_biotype = fifelse(gene_biotype == "protein_coding", "protein_coding",
                     fifelse(gene_biotype == "lncRNA", "lncRNA", "other"))
)]
stage_annotation <- stage_annotation[match(rownames(dge_all), gene_id_versioned)]
if (anyNA(stage_annotation$gene_id_versioned) || anyDuplicated(stage_annotation$gene_id_versioned)) {
  fail("Five-cohort gene annotation join failed")
}

eligibility <- as.data.table(readRDS(meta_path))[sample_id %in% colnames(dge_all)]
if (unit_qualified) {
  required_unit_columns <- c("sample_id", "dataset", "analysis_unit_id", "unit_evidence")
  roster <- fread(context_roster_path)
  if (!all(required_unit_columns %in% names(eligibility)) ||
      !all(required_unit_columns %in% names(roster))) fail("Qualified unit evidence is missing")
  selected_roster <- roster[sample_id %in% colnames(dge_all)]
  if (nrow(selected_roster) != n_samples_expected || anyDuplicated(selected_roster$sample_id) ||
      anyDuplicated(selected_roster$analysis_unit_id)) fail("Context roster is incomplete or repeated")
  joined_units <- merge(eligibility[, ..required_unit_columns], selected_roster[, ..required_unit_columns],
                        by = c("sample_id", "dataset"), suffixes = c("_model", "_context"))
  if (nrow(joined_units) != n_samples_expected ||
      anyNA(joined_units) ||
      any(!nzchar(joined_units$unit_evidence_model)) ||
      any(joined_units$analysis_unit_id_model != joined_units$analysis_unit_id_context) ||
      any(joined_units$unit_evidence_model != joined_units$unit_evidence_context)) {
    fail("Model metadata and source-qualified context roster disagree")
  }
} else {
  eligibility[, `:=`(analysis_unit_id = paste(dataset, sample_id, sep = "::"),
                    unit_evidence = "historical_library_record_not_identity_verified")]
}
eligibility <- eligibility[!is.na(inferred_sex) & !is.na(dataset)]
if (nrow(eligibility) != n_samples_expected || anyDuplicated(eligibility$sample_id) ||
    anyDuplicated(eligibility$analysis_unit_id)) fail("Biological-unit identity drift")
nmf <- fread(nmf_path)
programs <- paste0("P", 1:6)
if (!all(c("sample_id", programs) %in% names(nmf))) fail("Fixed P1-P6 scores are absent")
# An axis that failed the fixed reconstruction criteria is NA for every unit and
# stays in the output as unestimable. A partly missing or infinite column is an error.
nmf_scores <- as.matrix(nmf[, ..programs])
nmf_axis_unavailable <- colSums(is.na(nmf_scores)) == nrow(nmf_scores)
if (all(nmf_axis_unavailable) || any(is.infinite(nmf_scores)) ||
    anyNA(nmf_scores[, !nmf_axis_unavailable, drop = FALSE])) fail("Fixed P1-P6 scores are absent or non-finite")
if (unit_qualified && (anyDuplicated(nmf$sample_id) ||
    !setequal(nmf$sample_id, eligibility$sample_id))) fail("Fixed NMF projection roster differs from corrected model")
if (identical(Sys.getenv("STAGE_EXTENSION_PREFLIGHT_ONLY"), "1")) {
  cat("STAGE_EXTENSION_PREFLIGHT_PASS:", nrow(eligibility), "analysis units; no models fitted\n")
  quit(save = "no", status = 0L)
}

fit_contrast <- function(sample_table, group_column, reference, comparison, contrast_id, axis_name) {
  selected <- sample_table[get(group_column) %in% c(reference, comparison)]
  selected[, group := factor(get(group_column), levels = c(reference, comparison))]
  unestimable <- function(reason) list(
    result = data.table(),
    audit = data.table(axis = axis_name, contrast = contrast_id, reference = reference,
      comparison = comparison, n_reference = sum(selected$group == reference),
      n_comparison = sum(selected$group == comparison), n_samples = nrow(selected),
      n_cohorts = uniqueN(selected$dataset), covariates = "dataset+inferred_sex",
      bh_method = "Benjamini-Hochberg", bh_family_size = 0L,
      biological_unit = if (unit_qualified) "source_linked_donor_or_conservative_identity_cluster" else "historical_library_record",
      estimable = FALSE, failure_reason = reason),
    samples = selected[, .(sample_id, analysis_unit_id, unit_evidence, dataset, inferred_sex,
                           axis = axis_name, contrast = contrast_id, group = as.character(group))]
  )
  if (uniqueN(selected$group) != 2L || uniqueN(selected$dataset) < 2L) {
    return(unestimable("comparison_group_or_multiple_cohort_support_missing"))
  }
  if (uniqueN(selected$inferred_sex) < 2L) return(unestimable("inferred_sex_covariate_not_estimable"))
  setorder(selected, sample_id)
  dge <- dge_all[, match(selected$sample_id, colnames(dge_all))]
  info <- data.frame(
    dataset = droplevels(factor(selected$dataset)),
    inferred_sex = droplevels(factor(selected$inferred_sex)),
    group = selected$group,
    row.names = selected$sample_id
  )
  design <- model.matrix(~ dataset + inferred_sex + group, info)
  coef_name <- paste0("group", comparison)
  if (!coef_name %in% colnames(design) || qr(design)$rank != ncol(design)) {
    return(unestimable("rank_deficient_design"))
  }
  transformed <- voomWithQualityWeights(dge, design, plot = FALSE)
  fit <- eBayes(lmFit(transformed, design))
  coef_index <- match(coef_name, colnames(design))
  table <- topTable(fit, coef = coef_index, number = Inf, sort.by = "none")
  se <- fit$stdev.unscaled[, coef_index] * sqrt(fit$s2.post)
  critical <- qt(0.975, df = fit$df.total)
  names(critical) <- rownames(fit)
  genes <- rownames(dge)
  result <- cbind(
    stage_annotation,
    data.table(
      logFC = table[genes, "logFC"],
      SE = se[genes],
      CI_low = table[genes, "logFC"] - critical[genes] * se[genes],
      CI_high = table[genes, "logFC"] + critical[genes] * se[genes],
      t = table[genes, "t"],
      P.Value = table[genes, "P.Value"],
      FDR = p.adjust(table[genes, "P.Value"], method = "BH"),
      AveExpr = table[genes, "AveExpr"],
      axis = axis_name,
      contrast = contrast_id,
      reference = reference,
      comparison = comparison,
      n_reference = sum(selected$group == reference),
      n_comparison = sum(selected$group == comparison),
      n_samples = nrow(selected),
      n_cohorts = uniqueN(selected$dataset),
      covariates = "dataset+inferred_sex",
      bh_family_size = nrow(dge)
    )
  )
  if (nrow(result) != nrow(dge_all) || max(abs(result$FDR - p.adjust(result$P.Value, "BH"))) > 1e-12) {
    fail("Complete BH family check failed: ", contrast_id)
  }
  list(
    result = result,
    audit = data.table(
      axis = axis_name, contrast = contrast_id, reference = reference,
      comparison = comparison, n_reference = sum(selected$group == reference),
      n_comparison = sum(selected$group == comparison), n_samples = nrow(selected),
      n_cohorts = uniqueN(selected$dataset), covariates = "dataset+inferred_sex",
      bh_method = "Benjamini-Hochberg", bh_family_size = nrow(dge),
      biological_unit = if (unit_qualified) "source_linked_donor_or_conservative_identity_cluster" else "historical_library_record",
      estimable = TRUE, failure_reason = NA_character_
    ),
    samples = selected[, .(
      sample_id, analysis_unit_id, unit_evidence, dataset, inferred_sex,
      axis = axis_name, contrast = contrast_id, group = as.character(group)
    )]
  )
}

stage_samples <- eligibility[
  dataset %in% c("GSE130970", "GSE135251", "GSE162694") & fibrosis_stage %in% 0:4,
  .(sample_id, analysis_unit_id, unit_evidence, dataset, fibrosis_stage, inferred_sex)
]
stage_samples[, fib_group := paste0("F", fibrosis_stage)]

adjacent_pairs <- lapply(0:3, function(stage) c(paste0("F", stage), paste0("F", stage + 1L)))
adjacent_results <- lapply(adjacent_pairs, function(pair) {
  fit_contrast(stage_samples, "fib_group", pair[[1L]], pair[[2L]],
               paste0(pair[[1L]], "_to_", pair[[2L]]), "fibrosis_adjacent")
})

f0_results <- lapply(paste0("F", 1:4), function(comparison) {
  fit_contrast(stage_samples, "fib_group", "F0", comparison,
               paste0(comparison, "_vs_F0"), "fibrosis")
})

bin_nas <- function(x) fcase(
  x == 0, "NAS0", x <= 2, "NAS1-2", x <= 4, "NAS3-4",
  x <= 8, "NAS5-8", default = NA_character_
)
nas_samples <- eligibility[!is.na(nas_score)]
nas_samples[, nas_group := bin_nas(as.integer(nas_score))]
nas_samples <- nas_samples[!is.na(nas_group), .(
  sample_id, analysis_unit_id, unit_evidence, dataset, inferred_sex, nas_group
)]
if (anyDuplicated(nas_samples$analysis_unit_id)) fail("NAS biological unit is not unique")
nas_pairs <- list(c("NAS0", "NAS1-2"), c("NAS1-2", "NAS3-4"), c("NAS3-4", "NAS5-8"))
nas_adjacent_results <- lapply(nas_pairs, function(pair) {
  contrast_id <- paste0(gsub("-", "_", pair[[2L]]), "_vs_", gsub("-", "_", pair[[1L]]))
  fit_contrast(nas_samples, "nas_group", pair[[1L]], pair[[2L]], contrast_id, "NAS_adjacent")
})

# Common-reference NAS models against strict NAS0, mirroring the F-versus-F0
# set above so the two staging axes carry the same estimand where the figures
# place them side by side. The adjacent chain is retained on its own axis label
# and stays available to the supplement.
nas0_results <- lapply(c("NAS1-2", "NAS3-4", "NAS5-8"), function(comparison) {
  fit_contrast(nas_samples, "nas_group", "NAS0", comparison,
               paste0(gsub("-", "_", comparison), "_vs_NAS0"), "NAS")
})

all_fits <- c(adjacent_results, f0_results, nas_adjacent_results, nas0_results)
write_tsv_once(rbindlist(lapply(all_fits, `[[`, "result")), file.path(out_root, "stage_extension_all_gene_results.tsv"))
write_tsv_once(rbindlist(lapply(all_fits, `[[`, "audit")), file.path(out_root, "stage_extension_design_audit.tsv"))
write_tsv_once(rbindlist(lapply(all_fits, `[[`, "samples")), file.path(out_root, "stage_extension_sample_manifest.tsv"))
write_tsv_once(rbindlist(lapply(adjacent_results, `[[`, "result")), file.path(out_root, "stage_all_gene_results.tsv"))
write_tsv_once(unique(stage_samples), file.path(out_root, "stage_sample_manifest.tsv"))

# Per-cohort disease models support the Figure 3C replication alluvial. Each
# cohort retains its own complete BH family; no pooled p value is substituted.
cohort_design_audit <- list()
fit_cohort <- function(cohort_id) {
  selected <- eligibility[dataset == cohort_id & group_binary %in% c("Control", "Disease")]
  selected[, group := factor(group_binary, levels = c("Control", "Disease"))]
  varying_sex <- uniqueN(selected$inferred_sex) > 1L
  cohort_covariates <- if (varying_sex) "inferred_sex" else "none_constant_sex_omitted"
  record_cohort <- function(estimable, reason = NA_character_) {
    cohort_design_audit[[cohort_id]] <<- data.table(
      dataset = cohort_id, contrast = "Disease_vs_Control", estimable = estimable,
      failure_reason = reason, n_control = sum(selected$group == "Control"),
      n_disease = sum(selected$group == "Disease"), n_analysis_units = nrow(selected),
      covariates = cohort_covariates, bh_family_size = if (estimable) nrow(dge_all) else 0L)
  }
  if (uniqueN(selected$group) != 2L || anyNA(selected$inferred_sex)) {
    record_cohort(FALSE, "group_not_estimable_or_missing_sex")
    return(data.table())
  }
  setorder(selected, sample_id)
  dge <- dge_all[, match(selected$sample_id, colnames(dge_all))]
  info <- data.frame(inferred_sex = droplevels(factor(selected$inferred_sex)),
                     group = selected$group, row.names = selected$sample_id)
  design <- model.matrix(if (varying_sex) ~ inferred_sex + group else ~ group, info)
  if (!"groupDisease" %in% colnames(design) || qr(design)$rank != ncol(design)) {
    record_cohort(FALSE, "rank_deficient_design")
    return(data.table())
  }
  record_cohort(TRUE)
  transformed <- voomWithQualityWeights(dge, design, plot = FALSE)
  fit <- eBayes(lmFit(transformed, design))
  table <- topTable(fit, coef = "groupDisease", number = Inf, sort.by = "none")
  genes <- rownames(dge)
  cbind(stage_annotation, data.table(
    dataset = cohort_id, contrast = "Disease_vs_Control",
    logFC = table[genes, "logFC"], P.Value = table[genes, "P.Value"],
    FDR = p.adjust(table[genes, "P.Value"], "BH"), t = table[genes, "t"],
    n_control = sum(selected$group == "Control"), n_disease = sum(selected$group == "Disease"),
    n_participants = nrow(selected), covariates = cohort_covariates,
    bh_family_size = nrow(dge)
  ))
}
cohort_results <- rbindlist(lapply(sort(unique(eligibility$dataset)), fit_cohort))
write_tsv_once(cohort_results, file.path(out_root, "cohort_disease_all_gene_results.tsv"))
write_tsv_once(rbindlist(cohort_design_audit), file.path(out_root, "cohort_disease_design_audit.tsv"))
write_tsv_once(eligibility[, .(sample_id, analysis_unit_id, unit_evidence, dataset, group_binary, inferred_sex,
                               sex_source, fibrosis_stage, nas_score, diagnosis_harmonized)],
               file.path(out_root, "five_cohort_sample_manifest.tsv"))
cohort_summary <- eligibility[, .(
  n_participants = .N, n_control = sum(group_binary == "Control"),
  n_disease = sum(group_binary == "Disease"),
  n_fibrosis = sum(!is.na(fibrosis_stage)), n_nas = sum(!is.na(nas_score)),
  n_inferred_sex = sum(sex_source != "reported", na.rm = TRUE)
), by = dataset]
write_tsv_once(cohort_summary, file.path(out_root, "five_cohort_summary.tsv"))

# Continuous prespecified NMF axes. Scores are centered within program before
# cohort- and inferred-sex-adjusted group models. No dominant classes are used.
program_labels <- fread(program_label_path)
if (!all(programs %in% names(nmf))) fail("Expected fixed P1-P6 scores are absent")
nmf_long <- melt(nmf[, c("sample_id", programs), with = FALSE], id.vars = "sample_id",
                 variable.name = "program", value.name = "score")
nmf_long[, centered_score := as.numeric(scale(score)), by = program]
nmf_long[, program_label := program_labels$biological_label[match(program, program_labels$program_code)]]

program_group_summary <- function(samples, group_column, groups, axis_name) {
  joined <- merge(nmf_long, samples, by = "sample_id")
  joined <- joined[get(group_column) %in% groups]
  joined[, group := factor(get(group_column), levels = groups)]
  rbindlist(lapply(programs, function(program_id) {
    d <- joined[program == program_id]
    if (any(!is.finite(d$centered_score))) {
      return(data.table(axis = axis_name, program = program_id,
        program_label = program_labels$biological_label[match(program_id, program_labels$program_code)],
        group = groups, centered_activity = NA_real_, SE = NA_real_,
        n_samples = vapply(groups, function(value) sum(d$group == value), integer(1)),
        n_cohorts = vapply(groups, function(value) uniqueN(d[group == value, dataset]), integer(1)),
        covariates = "dataset+inferred_sex", prespecified = TRUE,
        estimable = FALSE, failure_reason = if (all(is.na(d$score))) "fixed_nmf_axis_unavailable_failed_reconstruction" else
          "program_score_has_no_finite_standardized_variation"))
    }
    rbindlist(lapply(groups, function(group_id) {
      x <- d[group == group_id]
      design_group <- tryCatch(model.matrix(~ dataset + inferred_sex + group, data = d),
                               error = function(e) NULL)
      if (!nrow(x) || !any(d$group == groups[[1L]]) || is.null(design_group) ||
          qr(design_group)$rank != ncol(design_group)) {
        return(data.table(axis = axis_name, program = program_id,
          program_label = program_labels$biological_label[match(program_id, program_labels$program_code)],
          group = group_id, centered_activity = NA_real_, SE = NA_real_,
          n_samples = nrow(x), n_cohorts = uniqueN(x$dataset),
          covariates = "dataset+inferred_sex", prespecified = TRUE,
          estimable = FALSE, failure_reason = "group_or_adjusted_design_not_estimable"))
      }
      fit <- lm(centered_score ~ dataset + inferred_sex + group, data = d)
      reference_row <- data.frame(
        dataset = factor(d$dataset[[1L]], levels = levels(factor(d$dataset))),
        inferred_sex = factor(d$inferred_sex[[1L]], levels = levels(factor(d$inferred_sex))),
        group = factor(group_id, levels = groups)
      )
      design_group <- model.matrix(~ dataset + inferred_sex + group, data = d)
      group_columns <- grep("^group", colnames(design_group), value = TRUE)
      estimate <- if (group_id == groups[[1L]]) 0 else unname(coef(fit)[paste0("group", group_id)])
      se <- if (group_id == groups[[1L]]) 0 else sqrt(vcov(fit)[paste0("group", group_id), paste0("group", group_id)])
      data.table(
        axis = axis_name, program = program_id,
        program_label = unique(d$program_label), group = group_id,
        centered_activity = estimate, SE = se,
        n_samples = nrow(x), n_cohorts = uniqueN(x$dataset),
        covariates = "dataset+inferred_sex", prespecified = TRUE,
        estimable = TRUE, failure_reason = NA_character_
      )
    }))
  }))
}

fib_programs <- program_group_summary(stage_samples, "fib_group", paste0("F", 0:4), "fibrosis")
# NAS activity reference group. The default starts at NAS1-2 because the legacy
# scores cover 1 of 46 NAS0 units. STAGE_NAS_ACTIVITY_REFERENCE=NAS0 adds NAS0 as
# the zero reference and requires a finite score on every available axis for every
# NAS0 unit; the fixed weights and scores are not changed.
nas_activity_reference <- Sys.getenv("STAGE_NAS_ACTIVITY_REFERENCE", "NAS1-2")
nas_activity_groups <- c("NAS0", "NAS1-2", "NAS3-4", "NAS5-8")
if (!nas_activity_reference %in% nas_activity_groups[1:2]) {
  fail("STAGE_NAS_ACTIVITY_REFERENCE must be NAS0 or NAS1-2")
}
nas_activity_groups <- nas_activity_groups[match(nas_activity_reference, nas_activity_groups):4]
if (nas_activity_reference == "NAS0") {
  nas0_ids <- nas_samples[nas_group == "NAS0", sample_id]
  nas0_scores <- nmf_scores[match(nas0_ids, nmf$sample_id), !nmf_axis_unavailable, drop = FALSE]
  if (!length(nas0_ids) || any(!is.finite(nas0_scores))) {
    fail("NAS0 reference requires finite fixed NMF scores for every NAS0 unit")
  }
  cat("NAS0_ACTIVITY_REFERENCE:", length(nas0_ids), "NAS0 units with finite scores\n")
}
nas_programs <- program_group_summary(nas_samples, "nas_group", nas_activity_groups, "NAS")
write_tsv_once(rbind(fib_programs, nas_programs), file.path(out_root, "prespecified_nmf_program_activity.tsv"))

source_manifest <- data.table(
  source_id = c("validated_five_cohort_release", "five_cohort_model_manifest", "corrected_fragment_dge",
                "matched_metadata", "gene_annotation", "nmf_scores", "program_labels"),
  path = c(stage_validation, model_manifest_path, dge_path, meta_path, annotation_path, nmf_path, program_label_path)
)
source_manifest[, `:=`(size_bytes = file.info(path)$size, sha256 = vapply(path, sha256_file, character(1)))]
if (unit_qualified) source_manifest <- rbind(source_manifest, data.table(
  source_id = "source_qualified_context_roster", path = context_roster_path,
  size_bytes = file.info(context_roster_path)$size, sha256 = sha256_file(context_roster_path)))
write_tsv_once(source_manifest, file.path(out_root, "input_manifest.tsv"))

writeLines(capture.output(sessionInfo()), file.path(out_root, "sessionInfo.txt"))
outputs <- setdiff(list.files(out_root, full.names = TRUE), file.path(out_root, "output_manifest.tsv"))
output_manifest <- data.table(
  relative_path = basename(outputs), size_bytes = file.info(outputs)$size,
  sha256 = vapply(outputs, sha256_file, character(1))
)
write_tsv_once(output_manifest, file.path(out_root, "output_manifest.tsv"))
cat("STAGE_EXTENSION_COMPLETE:", out_root, "\n")
