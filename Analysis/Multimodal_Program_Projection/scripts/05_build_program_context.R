#!/usr/bin/env Rscript

# Assemble assay-native Figure 4 results without summing modalities or ranking
# programs. Each column retains its own unit, testability, multiplicity control,
# and robustness flag.

suppressPackageStartupMessages(library(data.table))

BASE <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
ROOT <- file.path(BASE, "Analysis/Multimodal_Program_Projection")
OUT <- file.path(ROOT, "results")

paths <- list(
  registry = file.path(OUT, "frozen_programs.tsv"),
  protein = file.path(OUT, "proteomics/module_protein_results.tsv"),
  atac_dynamic = file.path(OUT, "atac/dynamic_program_accessibility.tsv"),
  atac_static = file.path(OUT, "atac/static_gwas_atac_summary.tsv"),
  spatial = file.path(OUT, "spatial/spatial_program_results.tsv")
)
missing <- names(paths)[!file.exists(unlist(paths))]
if (length(missing)) stop("Missing context input(s): ", paste(missing, collapse = ", "))

registry <- fread(paths$registry)
protein <- fread(paths$protein)
atac_dynamic <- fread(paths$atac_dynamic)
atac_static <- fread(paths$atac_static)
spatial <- fread(paths$spatial)

if (nrow(registry) != 22L || anyDuplicated(registry$program_id)) {
  stop("Frozen registry must contain 22 unique programs")
}

wide <- copy(registry[, .(
  program_id, display_order, lineage_order, cell_type, module, program_name,
  fig3_beta, fig3_q, fig3_direction, stability_score, stability_fail,
  gse244832_loo_jaccard
)])

static_wide <- dcast(
  atac_static,
  program_id ~ trait_scope,
  value.var = "n_genes",
  fill = 0
)
setnames(
  static_wide,
  c("direct_disease_PDFF", "liver_enzyme"),
  c("gwas_atac_direct_genes", "gwas_atac_enzyme_genes"),
  skip_absent = TRUE
)
wide <- merge(wide, static_wide, by = "program_id", all.x = TRUE, sort = FALSE)

for (cohort in c("GSE281367", "GSE244832")) {
  cohort_id <- cohort
  prefix <- fifelse(cohort == "GSE281367", "atac_external", "atac_internal")
  x <- atac_dynamic[cohort == cohort_id, .(
    program_id, effect, padj, testable, robust, sensitivity_sign_agree,
    n_measured, retained_l1_weight
  )]
  setnames(
    x,
    setdiff(names(x), "program_id"),
    paste0(prefix, "_", setdiff(names(x), "program_id"))
  )
  wide <- merge(wide, x, by = "program_id", all.x = TRUE, sort = FALSE)
}

protein_x <- protein[, .(
  program_id,
  protein_effect = effect,
  protein_padj = padj,
  protein_testable = testable,
  protein_robust = robust,
  protein_sensitivity_sign_agree = sensitivity_sign_agree,
  protein_n_measured = n_measured,
  protein_retained_l1_weight = retained_l1_weight
)]
wide <- merge(wide, protein_x, by = "program_id", all.x = TRUE, sort = FALSE)

for (dataset in c("GSE192741", "Vu_et_al_2025")) {
  dataset_id <- dataset
  prefix <- fifelse(dataset == "GSE192741", "spatial_gse", "spatial_vu")
  x <- spatial[dataset == dataset_id, .(
    program_id,
    raw_moran_z,
    residual_moran_z,
    residual_padj,
    zonation_moran_z,
    testable,
    robust,
    sensitivity_sign_agree,
    disease_delta_descriptive
  )]
  setnames(
    x,
    setdiff(names(x), "program_id"),
    paste0(prefix, "_", setdiff(names(x), "program_id"))
  )
  wide <- merge(wide, x, by = "program_id", all.x = TRUE, sort = FALSE)
}

setorder(wide, display_order)
if (nrow(wide) != 22L || anyDuplicated(wide$program_id)) {
  stop("Context integration changed the 22-program row contract")
}

long_specs <- list(
  list("Fig. 3 disease trajectory", "fig3_beta", "fig3_q", NULL, NULL, "donor-level beta", "defining"),
  list("Direct disease/PDFF genes", "gwas_atac_direct_genes", NULL, NULL, NULL, "program genes", "descriptive"),
  list("Liver-enzyme genes", "gwas_atac_enzyme_genes", NULL, NULL, NULL, "program genes", "descriptive"),
  list("External ATAC", "atac_external_effect", "atac_external_padj", "atac_external_testable", "atac_external_robust", "score difference", NULL),
  list("Internal ATAC", "atac_internal_effect", "atac_internal_padj", "atac_internal_testable", "atac_internal_robust", "score difference", NULL),
  list("Frozen-module liver DIA-MS", "protein_effect", "protein_padj", "protein_testable", "protein_robust", "adjusted score difference", NULL),
  list("GSE spatial organization", "spatial_gse_residual_moran_z", "spatial_gse_residual_padj", "spatial_gse_testable", "spatial_gse_robust", "matched-null Z", NULL),
  list("Vu spatial organization", "spatial_vu_residual_moran_z", "spatial_vu_residual_padj", "spatial_vu_testable", "spatial_vu_robust", "matched-null Z", NULL),
  list("GSE disease delta", "spatial_gse_disease_delta_descriptive", NULL, "spatial_gse_testable", NULL, "descriptive score difference", "descriptive")
)

long <- rbindlist(lapply(seq_along(long_specs), function(i) {
  s <- long_specs[[i]]
  d <- wide[, .(
    program_id, display_order, cell_type, module, program_name,
    assay_order = i,
    assay = s[[1]],
    value = get(s[[2]]),
    qvalue = if (is.null(s[[3]]) || !s[[3]] %in% names(wide)) NA_real_ else get(s[[3]]),
    testable = if (is.null(s[[4]]) || !s[[4]] %in% names(wide)) !is.na(get(s[[2]])) else as.logical(get(s[[4]])),
    robust = if (is.null(s[[5]]) || !s[[5]] %in% names(wide)) FALSE else as.logical(get(s[[5]])),
    unit = s[[6]]
  )]
  # A valid test whose robust rule does not pass is INDETERMINATE. It is never
  # tested_negative: that state requires a prespecified adequate-negative rule
  # that passed, and no such rule is defined for these assays.
  d[, status := fifelse(!testable, "untested", fifelse(robust, "robust", "indeterminate"))]
  if (!is.null(s[[7]])) d[testable == TRUE, status := s[[7]]]
  if (s[[1]] %in% c("Direct disease/PDFF genes", "Liver-enzyme genes")) {
    d[testable == TRUE, status := fifelse(value > 0, "descriptive", "no_overlap")]
  }
  d
}))

fwrite(wide, file.path(OUT, "program_context_wide.tsv"), sep = "\t", quote = FALSE)
fwrite(long, file.path(OUT, "program_context_long.tsv"), sep = "\t", quote = FALSE)

coverage <- long[, .(
  n_programs = .N,
  n_testable = sum(testable),
  n_robust = sum(robust),
  n_defining = sum(status == "defining"),
  n_descriptive = sum(status == "descriptive"),
  n_no_overlap = sum(status == "no_overlap"),
  n_indeterminate = sum(status == "indeterminate"),
  n_untested = sum(status == "untested")
), by = .(assay_order, assay, unit)]
fwrite(coverage, file.path(OUT, "program_context_coverage.tsv"), sep = "\t", quote = FALSE)

cat("[context] assembled 22 programs without a cross-modality sum or rank\n")
print(coverage)
