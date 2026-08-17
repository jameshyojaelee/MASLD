#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
})

BASE <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
RESULTS <- file.path(BASE, "Analysis/Multimodal_Program_Projection/results")
FIG5 <- file.path(BASE, "figures/main/fig5_molecular_context")

assert <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
}

same_numeric <- function(x, y, tolerance = 1e-12) {
  same_na <- is.na(x) == is.na(y)
  all(same_na) && all(abs(x[!is.na(x)] - y[!is.na(y)]) <= tolerance)
}

assert_grouped_bh <- function(x, p_col, q_col, by, label) {
  if (length(by)) {
    check <- x[, .(
      ok = same_numeric(get(q_col), p.adjust(get(p_col), method = "BH"))
    ), by = by]
  } else {
    check <- data.table(ok = same_numeric(x[[q_col]], p.adjust(x[[p_col]], method = "BH")))
  }
  assert(all(check$ok), paste(label, "BH values do not rederive from primary P values"))
}

programs <- fread(file.path(RESULTS, "frozen_programs.tsv"))
manifest <- fread(file.path(RESULTS, "freeze_manifest.tsv"))
assert(
  all(manifest[role == "input" & !is.na(expected_sha256), contract_match] == TRUE),
  "One or more Figure 3 freeze inputs drifted from the versioned contract"
)
assert(nrow(programs) == 22L, "Frozen registry must contain 22 programs")
assert(uniqueN(programs$program_id) == 22L, "Frozen program IDs must be unique")
expected_lineages <- c(
  hepatocytes = 7L,
  fibroblasts = 6L,
  macrophages = 6L,
  cholangiocytes = 3L
)
observed_lineages <- programs[, .N, by = cell_type]
observed_lineages <- setNames(observed_lineages$N, observed_lineages$cell_type)
assert(
  identical(unname(observed_lineages[names(expected_lineages)]), unname(expected_lineages)),
  "Frozen lineage counts differ from the pinned 7/6/6/3 registry"
)
assert(sum(programs$stability_fail) == 8L, "Pinned registry must contain 8 stability-flagged programs")

atac <- fread(file.path(RESULTS, "atac/dynamic_program_accessibility.tsv"))
atac_counts <- atac[, .(
  rows = .N,
  testable = sum(testable),
  robust = sum(robust)
), by = cohort]
assert(nrow(atac_counts) == 2L, "Expected two dynamic ATAC cohorts")
assert(all(atac_counts$rows == 22L), "Each dynamic ATAC cohort must contain 22 rows")
assert(all(atac_counts$testable == 22L), "Each dynamic ATAC cohort must test all 22 programs")
assert_grouped_bh(atac, "pvalue", "padj", "cohort", "Dynamic ATAC")
atac[, robust_rederived := testable & is.finite(padj) & padj < 0.05 &
  sign(effect) == sign(equal_effect) & sign(effect) == sign(leave_top_effect)]
assert(all(atac$robust == atac$robust_rederived), "Dynamic ATAC robust flags do not rederive")

static <- fread(file.path(RESULTS, "atac/static_gwas_atac_summary.tsv"))
coverage <- fread(file.path(RESULTS, "atac/static_study_coverage.tsv"))
assert(nrow(coverage) == 35L, "Expected 35 registered primary GWAS studies")
hits <- fread(file.path(RESULTS, "atac/static_gwas_atac_hits.tsv"))
if (nrow(hits)) {
  static_rederived <- hits[, .(
    n_source_loci_check = uniqueN(source_locus),
    n_study_loci_check = uniqueN(paste(study, source_locus, sep = "::")),
    n_variants_check = uniqueN(variant_id),
    n_genes_check = uniqueN(gene_symbol),
    n_studies_check = uniqueN(study),
    n_ancestries_check = uniqueN(ancestry)
  ), by = .(program_id, trait_scope)]
  static_check <- merge(static, static_rederived, by = c("program_id", "trait_scope"))
  for (stem in c("source_loci", "study_loci", "variants", "genes", "studies", "ancestries")) {
    assert(
      all(static_check[[paste0("n_", stem)]] == static_check[[paste0("n_", stem, "_check")]]),
      paste("Static ATAC", stem, "counts do not rederive from hit rows")
    )
  }
  static_global <- fread(file.path(RESULTS, "atac/static_gwas_atac_global_summary.tsv"))
  global_check <- hits[, .(
    n_source_rows = .N,
    n_programs = uniqueN(program_id),
    n_program_gene_incidences = uniqueN(paste(program_id, gene_symbol, sep = "::")),
    n_unique_genes = uniqueN(gene_symbol),
    n_program_variant_incidences = uniqueN(paste(program_id, variant_id, sep = "::")),
    n_unique_variants = uniqueN(variant_id),
    n_unique_studies = uniqueN(study)
  ), by = trait_scope]
  assert(
    isTRUE(all.equal(
      static_global[order(trait_scope)],
      global_check[order(trait_scope)],
      check.attributes = FALSE
    )),
    "Static ATAC global uniqueness summary does not rederive from hit rows"
  )
}

panel4c <- fread(file.path(RESULTS, "proteomics/panel4c_mrna_protein.tsv"))
assert(nrow(panel4c) == 25L, "Protected Panel 5C must contain 25 protein rows")
assert(
  all(panel4c$selection_conditioned == TRUE) &&
    all(panel4c$interpretation == "descriptive_same_cohort_reestimate"),
  "Panel 5C must remain explicitly labeled as selection-conditioned and descriptive"
)
panel4c_contract <- fread(file.path(RESULTS, "proteomics/panel4c_contract_audit.tsv"))
assert(nrow(panel4c_contract) == 1L && panel4c_contract$contract_match, "Panel 5C row contract drifted")
panel4c_sidecar <- fread(file.path(FIG5, "data/composite_mrna_protein_corrected.csv"))
assert(
  nrow(panel4c_sidecar) == 25L &&
    all(panel4c_sidecar$selection_conditioned == TRUE) &&
    all(panel4c_sidecar$interpretation == "descriptive_same_cohort_reestimate"),
  "Rendered Panel 5C sidecar lost its selection-conditioned provenance"
)
module_dia <- fread(file.path(RESULTS, "proteomics/module_protein_results.tsv"))
assert(sum(module_dia$testable) == 17L, "Current frozen-input DIA coverage must contain 17 testable programs")
assert_grouped_bh(module_dia[testable == TRUE], "pvalue", "padj", character(), "Module DIA-MS")
module_dia[, robust_rederived := testable & is.finite(padj) & padj < 0.05 &
  sign(effect) == sign(equal_effect) & sign(effect) == sign(leave_top_effect)]
assert(all(module_dia$robust == module_dia$robust_rederived), "Module DIA-MS robust flags do not rederive")

module_scores <- fread(file.path(RESULTS, "proteomics/module_protein_scores.tsv"))
assert(
  all(c("score", "equal_score", "leave_top_score") %in% names(module_scores)),
  "Module DIA-MS score table lacks weighted/equal/leave-top sample scores"
)
assert(
  nrow(module_scores) == 17L * 58L &&
    uniqueN(module_scores[, paste(program_id, sample_id, sep = "::")]) == nrow(module_scores),
  "Module DIA-MS score table must contain one row per 17 testable programs x 58 samples"
)

module_histology <- fread(file.path(RESULTS, "proteomics/module_protein_histology_masld.tsv"))
severity_features <- c("Steatosis", "Ballooning", "Inflammation", "Fibrosis", "NAS")
assert(
  nrow(module_histology) == 22L * length(severity_features) &&
    uniqueN(module_histology[, paste(program_id, feature, sep = "::")]) == nrow(module_histology),
  "Module protein-histology table must contain 22 programs x five features"
)
assert(
  setequal(module_histology$feature, severity_features) &&
    all(module_histology$cohort == "MASLD_only") &&
    all(module_histology$correction_scope == "17_testable_programs_x_5_histology_features") &&
    all(module_histology$interpretation == "independent_frozen_program_severity_association"),
  "Module protein-histology provenance or feature contract drifted"
)
assert(
  sum(module_histology$testable) == 17L * length(severity_features),
  "Module protein-histology table must contain 85 testable program-feature pairs"
)
assert_grouped_bh(
  module_histology[testable == TRUE], "pvalue", "padj", character(),
  "MASLD-only module protein-histology"
)

severity_meta <- fread(file.path(RESULTS, "proteomics/panel4c_metadata.tsv"))
severity_meta[, `:=`(
  age_z = as.numeric(scale(age)),
  bmi_z = as.numeric(scale(bmi))
)]
severity_meta <- severity_meta[, setdiff(
  names(severity_meta), c("group", "acquisition_batch")
), with = FALSE]
severity_data <- merge(
  module_scores,
  severity_meta,
  by = "sample_id",
  all.x = TRUE,
  sort = FALSE
)[group == "MASLD"]

rank_residualize_severity <- function(y, d, include_batch = TRUE) {
  rhs <- c("age_z", "bmi_z")
  if (include_batch && uniqueN(d$acquisition_batch[!is.na(d$acquisition_batch)]) > 1L) {
    rhs <- c("acquisition_batch", rhs)
  }
  if (uniqueN(d$sex[!is.na(d$sex)]) > 1L) rhs <- c(rhs, "sex")
  nuisance <- model.matrix(reformulate(rhs), data = d)
  ok <- is.finite(y) & apply(nuisance, 1, function(z) all(is.finite(z)))
  out <- rep(NA_real_, length(y))
  if (sum(ok) > ncol(nuisance) + 2L) {
    out[ok] <- residuals(lm.fit(
      nuisance[ok, , drop = FALSE],
      rank(y[ok], ties.method = "average")
    ))
  }
  out
}

rederive_severity_cor <- function(d, score_col, feature, include_batch = TRUE) {
  score_resid <- rank_residualize_severity(d[[score_col]], d, include_batch)
  hist_resid <- rank_residualize_severity(d[[feature]], d, include_batch)
  ok <- is.finite(score_resid) & is.finite(hist_resid)
  if (sum(ok) < 8L) return(c(n = sum(ok), rho = NA_real_, pvalue = NA_real_))
  tst <- suppressWarnings(cor.test(score_resid[ok], hist_resid[ok], method = "pearson"))
  c(n = sum(ok), rho = unname(tst$estimate), pvalue = tst$p.value)
}

severity_rederived <- rbindlist(lapply(unique(module_scores$program_id), function(pid) {
  d <- severity_data[program_id == pid]
  rbindlist(lapply(severity_features, function(feature) {
    primary <- rederive_severity_cor(d, "score", feature)
    equal <- rederive_severity_cor(d, "equal_score", feature)
    leave_top <- rederive_severity_cor(d, "leave_top_score", feature)
    batch_rho <- vapply(c("2019", "2020"), function(batch_id) {
      unname(rederive_severity_cor(
        d[acquisition_batch == batch_id], "score", feature, include_batch = FALSE
      )["rho"])
    }, numeric(1))
    data.table(
      program_id = pid,
      feature = feature,
      n_check = as.integer(primary["n"]),
      rho_check = unname(primary["rho"]),
      pvalue_check = unname(primary["pvalue"]),
      equal_rho_check = unname(equal["rho"]),
      leave_top_rho_check = unname(leave_top["rho"]),
      rho_2019_check = batch_rho[[1]],
      rho_2020_check = batch_rho[[2]]
    )
  }))
}))
severity_check <- merge(
  module_histology[testable == TRUE],
  severity_rederived,
  by = c("program_id", "feature"),
  all.x = TRUE,
  sort = FALSE
)
assert(all(severity_check$n == severity_check$n_check), "Module protein-histology sample counts do not rederive")
for (stem in c("rho", "pvalue", "equal_rho", "leave_top_rho", "rho_2019", "rho_2020")) {
  assert(
    same_numeric(severity_check[[stem]], severity_check[[paste0(stem, "_check")]], tolerance = 1e-10),
    paste("Module protein-histology", stem, "values do not rederive")
  )
}
module_histology[, sensitivity_rederived := testable & is.finite(rho) &
  is.finite(equal_rho) & is.finite(leave_top_rho) &
  sign(rho) == sign(equal_rho) & sign(rho) == sign(leave_top_rho)]
module_histology[, batch_rederived := testable & is.finite(rho) &
  is.finite(rho_2019) & is.finite(rho_2020) &
  sign(rho) == sign(rho_2019) & sign(rho) == sign(rho_2020)]
module_histology[, robust_rederived := testable & is.finite(padj) & padj < 0.05 &
  sensitivity_rederived & batch_rederived]
assert(
  all(module_histology$sensitivity_sign_agree == module_histology$sensitivity_rederived) &&
    all(module_histology$batch_sign_agree == module_histology$batch_rederived) &&
    all(module_histology$robust == module_histology$robust_rederived),
  "Module protein-histology sensitivity, batch, or robust flags do not rederive"
)

burden <- fread(file.path(RESULTS, "proteomics/module_protein_histology_burden_masld.tsv"))
assert(
  nrow(burden) == 22L && uniqueN(burden$program_id) == 22L && sum(burden$testable) == 17L,
  "Histology-burden table must contain 22 programs and 17 testable rows"
)
assert(
  all(burden[testable == TRUE, n_permutations] == 9999L) &&
    all(grepl("rank-residual partial Spearman", burden[testable == TRUE, inference_method], fixed = TRUE)),
  "Histology-burden inference must use 9,999 batch-stratified rank-residual permutations"
)
assert_grouped_bh(burden[testable == TRUE], "pvalue", "padj", character(), "Histology burden")
burden[, sensitivity_rederived := testable & is.finite(rho) &
  is.finite(equal_rho) & is.finite(leave_top_rho) &
  sign(rho) == sign(equal_rho) & sign(rho) == sign(leave_top_rho)]
burden[, batch_rederived := testable & is.finite(rho) &
  is.finite(rho_2019) & is.finite(rho_2020) &
  sign(rho) == sign(rho_2019) & sign(rho) == sign(rho_2020)]
burden[, robust_rederived := testable & is.finite(padj) & padj < 0.05 &
  sensitivity_rederived & batch_rederived]
assert(
  all(burden[testable == TRUE, robust == robust_rederived]) &&
    !any(burden[testable != TRUE, robust %in% TRUE]),
  "Histology-burden robust calls do not rederive from the corrected inference"
)

spatial <- fread(file.path(RESULTS, "spatial/spatial_program_results.tsv"))
spatial_counts <- spatial[, .(
  rows = .N,
  testable = sum(testable),
  robust = sum(robust)
), by = dataset]
assert(nrow(spatial_counts) == 2L, "Expected two spatial datasets")
assert(all(spatial_counts$rows == 22L), "Each spatial dataset must contain 22 rows")
assert(all(spatial_counts$testable == 21L), "Each spatial dataset must contain 21 testable programs")
assert(all(spatial[testable == TRUE, n_null] == 9999L), "Testable spatial rows must use 9,999 primary null draws")
assert_grouped_bh(spatial, "residual_pvalue", "residual_padj", "dataset", "Spatial residual")
assert_grouped_bh(spatial, "raw_pvalue", "raw_padj", "dataset", "Spatial raw")
assert_grouped_bh(spatial, "zonation_pvalue", "zonation_padj", "dataset", "Spatial zonation")
assert(
  all(is.finite(spatial[testable == TRUE, equal_residual_null_mean])) &&
    all(is.finite(spatial[testable == TRUE, leave_top_residual_null_mean])),
  "Spatial sensitivity-specific null means must be finite"
)
spatial[, robust_rederived := testable & residual_padj < 0.05 & sensitivity_sign_agree]
assert(all(spatial$robust == spatial$robust_rederived), "Spatial robust flags do not rederive")
assert(
  all(nchar(spatial[testable == TRUE, matched_set_sha256]) == 64L),
  "Spatial matched-set hashes are missing"
)
matching <- fread(file.path(RESULTS, "spatial/spatial_matching_audit.tsv"))
assert(
  all(c("expression_match_relaxation", "lineage_match_relaxation",
        "mean_absolute_lineage_bin_difference", "max_absolute_lineage_bin_difference",
        "n_unique_control_genes") %in% names(matching)),
  "Spatial matching audit lacks lineage-specific matching fields"
)
assert(
  all(matching$lineage_match_relaxation <= 2L) &&
    all(matching$n_unique_control_genes >= 1L) &&
    all(matching$max_absolute_lineage_bin_difference <= 2 + 1e-12) &&
    all(is.finite(matching$target_lineage_correlation)) &&
    all(is.finite(matching$mean_control_lineage_correlation)),
  "Spatial lineage matching used a control outside the prespecified correlation window"
)
graphs <- fread(file.path(RESULTS, "spatial/spatial_graph_audit.tsv"))
expected_arrays <- c(GSE192741 = 5L, Vu_et_al_2025 = 10L)
observed_arrays <- graphs[, .N, by = dataset]
observed_arrays <- setNames(observed_arrays$N, observed_arrays$dataset)
assert(
  identical(unname(observed_arrays[names(expected_arrays)]), unname(expected_arrays)) &&
    !anyDuplicated(graphs[, paste(dataset, sample_id, sep = "::")]),
  "Spatial graph audit does not contain the expected 5 GSE and 10 Vu arrays"
)
assert(all(graphs$n_tissue_islands >= 1L), "One or more spatial arrays lack a valid tissue-island graph")
assert(
  all(graphs$n_spots_in_graph + graphs$n_spots_excluded_small_islands == graphs$n_spots) &&
    all(graphs$n_edges > 0L),
  "Spatial graph audit does not account for every spot or contains an edgeless array"
)
assert(
  all(graphs$max_observed_edge_distance <= graphs$max_allowed_edge_distance + 1e-8),
  "A spatial graph contains an edge across the tissue-island distance threshold"
)
null_summary <- fread(file.path(RESULTS, "spatial/spatial_null_summary.tsv"))
assert(
  nrow(null_summary) == 3L * spatial[testable == TRUE, .N] &&
    !anyDuplicated(null_summary[, paste(program_id, dataset, statistic, sep = "::")]),
  "Spatial null-summary rows are incomplete or duplicated"
)
assert(
  all(null_summary[statistic %chin% c("residual_moran_i", "zonation_moran_i"),
                   n_null == 9999L]) &&
    all(null_summary[statistic == "raw_moran_i", n_null >= ceiling(0.95 * 9999)]),
  "Spatial primary/sensitivity nulls lack 9,999 finite draws or a raw null lacks 95% finite draws"
)
assert(
  all(null_summary$null_sd > 0) &&
    all(null_summary$q001 <= null_summary$q010) &&
    all(null_summary$q010 <= null_summary$q050) &&
    all(null_summary$q050 <= null_summary$q500) &&
    all(null_summary$q500 <= null_summary$q950) &&
    all(null_summary$q950 <= null_summary$q990) &&
    all(null_summary$q990 <= null_summary$q999),
  "Spatial null-summary variance or quantile ordering is invalid"
)

context <- fread(file.path(RESULTS, "program_context_wide.tsv"))
assert(nrow(context) == 22L && uniqueN(context$program_id) == 22L, "Context matrix must contain 22 unique rows")
context_long <- fread(file.path(RESULTS, "program_context_long.tsv"))
assert(
  all(context_long[status == "no_overlap", value] == 0),
  "Context matrix no-overlap states must represent exact zero counts"
)

firewall <- fread(file.path(FIG5, "panels/data/fig5a_input_firewall.tsv"))
assert(
  nrow(firewall) == 3L &&
    setequal(firewall$branch, c(
      "prioritized_gene_context", "prespecified_program_projection", "fixed_protein_display"
    )),
  "Panel 5A input-firewall sidecar is incomplete"
)

protein_triage <- fread(file.path(FIG5, "panels/data/fig5b_protein_triage_summary.tsv"))
protein_check <- fread(file.path(
  RESULTS, "proteomics/protein_de_adjusted.tsv"
))[, .(gene, protein_logFC = logFC, protein_padj = padj)]
bulk_check <- fread(
  file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"),
  select = c("human_symbol", "bulk_logFC")
)[, .(gene = human_symbol, bulk_logFC)]
universe_check <- trimws(readLines(file.path(
  BASE, "Analysis/Spatial/results/universe_validation/prioritized_universe_FINAL.txt"
)))
universe_check <- universe_check[nzchar(universe_check)]
protein_check <- merge(protein_check, bulk_check, by = "gene", all = FALSE, sort = FALSE)
protein_check <- protein_check[is.finite(bulk_logFC) & is.finite(protein_logFC)]
protein_check[, in_universe := gene %chin% universe_check]
protein_check[, class := fifelse(
  !in_universe, "rest",
  fifelse(
    sign(bulk_logFC) == sign(protein_logFC) & protein_padj < 0.05,
    "confirmed",
    fifelse(
      sign(bulk_logFC) != sign(protein_logFC) & protein_padj < 0.05,
      "discordant", "not_significant"
    )
  )
)]
protein_expected <- protein_check[, .(
  n_full_measured = .N,
  n_prioritized_measured = sum(in_universe),
  rho_full = cor(bulk_logFC, protein_logFC, method = "spearman"),
  direction_full = mean(sign(bulk_logFC) == sign(protein_logFC)),
  rho_prioritized = cor(
    bulk_logFC[in_universe], protein_logFC[in_universe], method = "spearman"
  ),
  direction_prioritized = mean(
    sign(bulk_logFC[in_universe]) == sign(protein_logFC[in_universe])
  ),
  n_significant_concordant = sum(class == "confirmed"),
  n_significant_discordant = sum(class == "discordant"),
  n_not_significant = sum(class == "not_significant")
)]
assert(
  nrow(protein_triage) == 1L &&
    protein_triage$protein_model == "PXD051911_quantile_batch_age_bmi_sex_adjusted" &&
    all(protein_triage[, .(
      n_full_measured, n_prioritized_measured,
      n_significant_concordant, n_significant_discordant, n_not_significant
    )] == protein_expected[, .(
      n_full_measured, n_prioritized_measured,
      n_significant_concordant, n_significant_discordant, n_not_significant
    )]) &&
    all(vapply(
      c("rho_full", "direction_full", "rho_prioritized", "direction_prioritized"),
      function(stem) same_numeric(
        protein_triage[[stem]], protein_expected[[stem]], tolerance = 1e-12
      ),
      logical(1)
    )),
  "Panel 5B sidecar does not rederive from the adjusted protein and bulk inputs"
)

map_selection <- fread(file.path(FIG5, "panels/data/fig5f_spatial_program_maps_selection.tsv"))
assert(
  nrow(map_selection) == 4L &&
    setequal(map_selection$program_id, c("hepatocytes::14", "fibroblasts::6")) &&
    setequal(map_selection$dataset, c("GSE192741", "Vu_et_al_2025")) &&
    all(map_selection[, uniqueN(sample_id), by = dataset]$V1 == 1L),
  "Panel 5F map selection must contain two prespecified programs on one median section per cohort"
)

summary_selection <- fread(file.path(FIG5, "panels/data/fig5e_multimodal_program_summary_selection.tsv"))
summary_expected <- merge(
  context[, .(
    program_id,
    open_promoter_genes = fcoalesce(as.integer(gwas_atac_enzyme_genes), 0L)
  )],
  burden[, .(program_id, histology_burden_robust = robust %in% TRUE)],
  by = "program_id", all.x = TRUE, sort = FALSE
)
summary_expected <- merge(
  summary_expected,
  spatial[, .(
    spatial_shared = all(testable == TRUE & robust == TRUE)
  ), by = program_id],
  by = "program_id", all.x = TRUE, sort = FALSE
)
summary_expected <- summary_expected[
  open_promoter_genes > 0L | histology_burden_robust == TRUE | spatial_shared == TRUE
]
assert(
  setequal(summary_selection$program_id, summary_expected$program_id) &&
    nrow(summary_selection) == nrow(summary_expected) &&
    all(summary_selection$displayed_n == nrow(summary_expected)) &&
    all(summary_selection$frozen_universe_n == 22L) &&
    all(summary_selection$selection_rule == paste0(
      "open_promoter_gene OR robust_histology_burden OR robust_spatial_in_both"
    )),
  "Panel 5E display universe does not rederive from its assay-native inclusion rule"
)

pdfs <- file.path(FIG5, c(
  "panels/fig5a_input_firewall.pdf",
  "panels/fig5b_protein_triage.pdf",
  "panels/fig5c_mrna_protein_composite.pdf",
  "panels/fig5d_snatac_accessibility.pdf",
  "panels/fig5e_multimodal_program_summary.pdf",
  "panels/fig5f_spatial_program_maps.pdf"
))
assert(all(file.exists(pdfs)), "One or more final Figure 5 PDFs are missing")
assert(all(file.info(pdfs)$size > 1000), "One or more final Figure 5 PDFs are empty")
assert(
  length(list.files(file.path(FIG5, "panels"), pattern = "[.]png$", full.names = TRUE, recursive = TRUE)) == 0L,
  "Figure 5 panels directory must remain PDF-only"
)
for (pdf in pdfs) {
  info <- system2("pdfinfo", pdf, stdout = TRUE, stderr = TRUE)
  pages <- grep("^Pages:", info, value = TRUE)
  assert(length(pages) == 1L && grepl("Pages:[[:space:]]+1$", pages), paste("Expected one-page PDF:", basename(pdf)))
}

release <- fread(file.path(FIG5, "CANONICAL_MAIN_PANELS.tsv"))
assert(
  nrow(release) == 6L && identical(release$callout, paste0("5", LETTERS[1:6])) &&
    all(release$source != release$filename),
  "Canonical Figure 5 manifest must map six panels to real generator scripts"
)
release_paths <- file.path(FIG5, release$filename)
release_sha <- vapply(
  release_paths, digest, character(1), algo = "sha256", file = TRUE, serialize = FALSE
)
assert(
  identical(unname(release_sha), release$sha256) && all(file.info(release_paths)$size == release$bytes),
  "Canonical Figure 5 manifest checksums or byte counts do not match the rendered PDFs"
)

cat("[validate] PASS: pinned inputs, corrected proteomics/histology inference, lineage-matched spatial audit, sidecars, and six canonical PDFs\n")
print(spatial_counts)
print(atac_counts)
