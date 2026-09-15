#!/usr/bin/env Rscript
# Additive Gene Catalog evidence tables; does not modify inherited memberships.
source("scripts/analysis/histology_anchored_continuum/translation/lib_translation.R")
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L, !dir.exists(args[[1]]))
out <- args[[1]]
dir.create(out, recursive = TRUE)
set.seed(20260906); setDTthreads(4)
ml <- "RNA-seq/results/histology_anchored_continuum/molecular_layers/hac-molecular-layers-20260818T173348Z"
hac <- "RNA-seq/results/histology_anchored_continuum/candidates/hac-continuum-20260818T024923Z"
ctx <- "Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-candidate-2026-08-11-r1"
pr <- "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot"
sp <- "Analysis/Multimodal_Program_Projection/candidates/spatial-resource-candidate-2026-08-11"
protein_root <- "Analysis/Multimodal_Program_Projection/results/proteomics"
inputs <- character()
read <- function(path, ...) {
  inputs <<- union(inputs, normalizePath(path, mustWork = TRUE))
  fread(path, ...)
}
write <- function(d, name) fwrite(d, file.path(out, name), sep = "\t", na = "NA")
writeLines(c("seed=20260906", "state=candidate_not_adopted", "coloc_threshold=>0.5",
  "direction_mass_threshold=0.95_of_original_total", "marginal_effects_only=true",
  "multi_signal_perturbation_direction=unresolved", "palindromic_direction=unresolved",
  "evaluation_cohorts=GSE162694;GSE213621", "program_BH_family=117_per_axis_per_model",
  "composition=15_orthonormal_logratio_coordinates_of_16_BayesPrism_fractions",
  "composition_zero_floor=1e-6", "primary_model=stage_factor+inferred_sex+axis_z",
  "composition_model=primary+15_logratio_coordinates", "new_nonlinear_tests=none",
  "coherence_excluded=true", "windows=9_frozen_descriptive_only",
  "new_clinical_prediction_or_cutoffs=none"), file.path(out, "analysis_contract.txt"))

registry <- read(file.path(pr, "program_registry_v2.tsv"))
members <- read(file.path(pr, "program_membership_v2.tsv"))
stopifnot(nrow(registry) == 117L, !anyDuplicated(registry$program_uid))
write(registry[, .(program_uid, cell_type, module, module_name, membership_sha256)],
      "frozen_program_registry.tsv")
members <- unique(members[, .(program_uid, gene_name = mapped_symbol,
                              original_l1_weight, membership_sha256)])
members <- members[!is.na(gene_name) & nzchar(gene_name)]
write(members, "gene_program_links.tsv.gz")

# Same participant set for unadjusted and composition-adjusted sensitivity.
# Frozen score units are not recalibrated separately between models.
message("Fitting within-stage composition sensitivity")
scores <- read(file.path(ml, "programs/hotspot/hotspot_participant_scores.tsv.gz"))
axes <- read(file.path(hac, "projection/score_registry.tsv"))
stopifnot(!anyDuplicated(scores[, .(dataset, sample_id, program_uid)]),
          !anyDuplicated(axes[, .(dataset, sample_id)]))
fit_rows <- list(); denominator <- list(); fractions <- list()
for (cohort in c("GSE162694", "GSE213621")) {
  path <- sprintf("Analysis/Deconvolution/results/%s/%s_bayesprism_proportions.tsv", cohort, cohort)
  inputs <- union(inputs, normalizePath(path))
  f <- read.table(path, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
  stopifnot(ncol(f) == 16, all(is.finite(as.matrix(f))), all(as.matrix(f) >= 0),
            max(abs(rowSums(f) - 1)) < 1e-5, !anyDuplicated(rownames(f)))
  fractions[[cohort]] <- data.table(dataset = cohort, lineage = names(f),
    n_participants = nrow(f), median_fraction = apply(f, 2, median),
    n_below_floor = colSums(f < 1e-6))
  h <- contr.helmert(16); h <- sweep(h, 2, sqrt(colSums(h^2)), "/")
  ilr <- as.data.table(log(pmax(as.matrix(f), 1e-6)) %*% h)
  setnames(ilr, sprintf("ilr_%02d", 1:15)); ilr[, sample_id := rownames(f)]
  a <- axes[dataset == cohort]
  d <- merge(scores[dataset == cohort], a[, .(sample_id, fixed_projection_raw, signature_pc1_raw)],
             by = "sample_id", all.x = TRUE)
  d <- merge(d, ilr, by = "sample_id", all.x = TRUE)
  for (axis in c("fixed_projection", "signature_pc1")) {
    for (uid in registry$program_uid) {
      z <- d[program_uid == uid]
      z[, axis_raw := get(paste0(axis, "_raw"))]
      usable <- is.finite(z$outcome_z) & is.finite(z$axis_raw) &
        !is.na(z$fibrosis_stage) & !is.na(z$inferred_sex)
      comp_ok <- complete.cases(z[, grep("^ilr_", names(z), value = TRUE), with = FALSE])
      denominator[[length(denominator) + 1L]] <- data.table(program_uid = uid,
        dataset = cohort, axis_id = axis, n_score_rows = nrow(z),
        n_complete_primary = sum(usable), n_composition_matched = sum(usable & comp_ok))
      z <- z[usable & comp_ok]
      z[, axis_z := as.numeric(scale(axis_raw))]
      for (arm in c("stage_sex", "stage_sex_composition")) {
        r <- tr_fit(z, arm == "stage_sex_composition")
        r[, `:=`(program_uid = uid, dataset = cohort, axis_id = axis, model = arm)]
        fit_rows[[length(fit_rows) + 1L]] <- r
      }
    }
  }
}
fits <- rbindlist(fit_rows)
fits[, `:=`(q_value = tr_bh(p_value), hc3_q_value = tr_bh(hc3_p_value)),
     by = .(dataset, axis_id, model)]
meta <- fits[, tr_meta(.SD), by = .(program_uid, axis_id, model)]
meta[, q_value := tr_bh(p_value), by = .(axis_id, model)]
meta[, supported := is.finite(q_value) & q_value < .05 & direction_consistent]
# Stronger label requires BH in both separate cohorts, in both co-primary axes.
strict <- fits[, .(four_cohort_axis_BH = all(is.finite(q_value) & q_value < .05),
  four_cohort_axis_HC3_BH = all(is.finite(hc3_q_value) & hc3_q_value < .05),
  all_four_directions = all(is.finite(beta)) & uniqueN(sign(beta)) == 1L),
  by = .(program_uid, model)]
summary <- meta[, .(both_axes_supported = all(supported) & .N == 2L,
  fixed_beta = beta[axis_id == "fixed_projection"],
  signature_beta = beta[axis_id == "signature_pc1"]), by = .(program_uid, model)]
summary <- merge(summary, strict, by = c("program_uid", "model"))
summary[, both_axes_supported := both_axes_supported & all_four_directions]
write(fits, "program_cohort_effects.tsv")
write(meta, "program_meta_effects.tsv")
write(summary, "program_composition_summary.tsv")
write(rbindlist(denominator), "program_denominators.tsv")
write(rbindlist(fractions), "composition_observability.tsv")
existing_meta <- read(file.path(ml, "programs/hotspot/hotspot_meta_analysis.tsv"))
check <- merge(meta[model == "stage_sex"], existing_meta[, .(program_uid, axis_id,
  frozen_beta = fixed_beta)], by = c("program_uid", "axis_id"))
check[, absolute_beta_difference := abs(beta - frozen_beta)]
stopifnot(max(check$absolute_beta_difference, na.rm = TRUE) < 1e-8)
write(check, "frozen_effect_reproduction.tsv")

message("Loading complete inherited evidence universe")
aggregate_path <- "GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv"
manifest <- read(file.path(ctx, "genetics/prepared/replay_input_manifest.tsv"))
expected <- manifest[role == "promoted_aggregate", sha256]
actual <- strsplit(system2("sha256sum", aggregate_path, stdout = TRUE), " +")[[1]][[1]]
stopifnot(length(expected) == 1L, actual == expected)
g <- read(aggregate_path, select = c("gwas_name", "gene", "ensembl", "chr",
  "PP.H4.abf", "PP.H4.susie", "ancestry", "method"))
tiers <- read("GWAS/finemapping/config/gwas_trait_tier.tsv")
g <- merge(g, tiers, by.x = "gwas_name", by.y = "study_name", all.x = TRUE)
write(g[is.na(ensembl) | !nzchar(ensembl)], "genetic_rows_without_gene_identifier.tsv")
g <- g[!is.na(ensembl) & nzchar(ensembl)]
stopifnot(!anyNA(g$tier), !anyDuplicated(g[, .(gwas_name, ensembl)]))
g[, `:=`(susie_supported = is.finite(PP.H4.susie) & PP.H4.susie > .5,
          abf_supported = is.finite(PP.H4.abf) & PP.H4.abf > .5,
          trait_group = tr_trait(trait, tier))]
tested <- g[, .(n_gene_trait_tests = .N, n_susie_testable = sum(is.finite(PP.H4.susie)),
  n_abf_testable = sum(is.finite(PP.H4.abf)),
  primary_susie_supported = any(susie_supported & tier <= 2),
  primary_abf_supported = any(abf_supported & tier <= 2),
  any_portfolio_susie_supported = any(susie_supported)), by = .(ensembl)]
write(g[, .(n_gene_trait_tests = .N, n_susie_testable = sum(is.finite(PP.H4.susie)),
  n_abf_testable = sum(is.finite(PP.H4.abf)), n_susie_supported = sum(susie_supported),
  n_abf_supported = sum(abf_supported)), by = .(gwas_name, trait, trait_group, ancestry)],
  "genetic_testing_denominators.tsv")
write(tested, "gene_genetic_testability.tsv")
eligible <- g[susie_supported | abf_supported]
plan <- read(file.path(ctx, "genetics/prepared/replay_plan.tsv"))
stopifnot(nrow(plan) == 816, uniqueN(plan$gene) == 462)
pairs <- list(); variants <- list()
for (i in seq_len(nrow(plan))) {
  p <- plan[i]
  folder <- file.path(ctx, "genetics/replay_execution", paste0("batch_", p$batch_id),
                      "exports", p$gwas_name, paste0("chr", p$chr), p$ensembl)
  s <- read(file.path(folder, "signal_pairs.tsv"))
  s[, `:=`(n_gwas_signals = uniqueN(idx1), n_eqtl_signals = uniqueN(idx2))]
  pairs[[i]] <- s
  v <- read(file.path(folder, "variant_posteriors.tsv.gz"))
  # All exported posterior mass is retained, not only the maximum or credible set.
  variants[[i]] <- v[, .(gwas_name, ensembl, chr = as.integer(chr), signal_pair_index,
     position = as.integer(hg19_position), allele1, allele2, SNP.PP.H4)]
}
pairs <- rbindlist(pairs); variants <- rbindlist(variants)
stopifnot(nrow(pairs) == 5362L, nrow(variants) == 22356697L,
          !anyDuplicated(pairs[, .(gwas_name, ensembl, signal_pair_index)]))
# Read each chromosome and each GWAS only once; retain exact variant/gene keys.
eqtl <- list()
for (chrom in sort(unique(variants$chr))) {
  message("Reading signed eQTL chromosome ", chrom)
  e <- read(sprintf("data/broadaway_eqtl/chr%s_marginal_summary_results.tsv", chrom),
    select = c("ENSG", "CHR", "POS", "EA", "NEA", "Beta", "N"))
  setnames(e, c("ENSG", "CHR", "POS", "EA", "NEA", "Beta", "N"),
    c("ensembl", "chr", "position", "eqtl_ea", "eqtl_nea", "eqtl_beta", "eqtl_n"))
  e[, ensembl := sub("\\..*$", "", ensembl)]
  need <- unique(variants[chr == chrom, .(ensembl, chr, position)])
  e <- merge(e, need, by = c("ensembl", "chr", "position"))
  e[, duplicate_effect := .N > 1L, by = .(ensembl, chr, position)]
  e[duplicate_effect == TRUE, eqtl_beta := NA_real_]
  eqtl[[as.character(chrom)]] <- unique(e, by = c("ensembl", "chr", "position"))
}
eqtl <- rbindlist(eqtl)
variants <- merge(variants, eqtl, by = c("ensembl", "chr", "position"), all.x = TRUE)
rm(eqtl); gc()
gwreg <- read("GWAS/finemapping/config/gwas_registry.tsv")
signed_summaries <- list(); allele_qc <- list(); n_variants_written <- 0L
for (study in unique(variants$gwas_name)) {
  message("Harmonizing ", study)
  path <- file.path("GWAS/finemapping", gwreg[study_name == study, sumstats_path])
  stopifnot(length(path) == 1L)
  b <- read(path, select = c("chromosome", "position", "allele1", "allele2", "beta"))
  setnames(b, c("chromosome", "beta"), c("chr", "gwas_beta"))
  b[, chr := as.integer(chr)]
  v <- variants[gwas_name == study]
  keys <- c("chr", "position", "allele1", "allele2")
  b <- merge(b, unique(v[, ..keys]), by = keys)
  b[, duplicate_gwas := .N > 1L, by = keys]
  b[duplicate_gwas == TRUE, gwas_beta := NA_real_]
  b <- unique(b, by = keys)
  v <- merge(v, b, by = keys, all.x = TRUE)
  h <- tr_harmonize(v$allele1, v$allele2, v$eqtl_ea, v$eqtl_nea, v$gwas_beta, v$eqtl_beta)
  v <- cbind(v, h)
  signed_summaries[[study]] <- v[, c(tr_direction(SNP.PP.H4, marginal_expression_direction),
    list(n_variants = .N, posterior_total = sum(SNP.PP.H4),
         eqtl_n_min = if (all(is.na(eqtl_n))) NA_real_ else min(eqtl_n, na.rm = TRUE),
         eqtl_n_max = if (all(is.na(eqtl_n))) NA_real_ else max(eqtl_n, na.rm = TRUE))),
    by = .(gwas_name, ensembl, signal_pair_index)]
  allele_qc[[study]] <- v[, .(n_variant_instances = .N, posterior_mass_sum = sum(SNP.PP.H4)),
                          by = .(gwas_name, allele_state)]
  # Long variant evidence is partitioned by study to keep memory and I/O bounded.
  dir.create(file.path(out, "signed_variants"), showWarnings = FALSE)
  write(v, paste0("signed_variants/", study, ".tsv.gz"))
  n_variants_written <- n_variants_written + nrow(v)
}
stopifnot(n_variants_written == 22356697L)
rm(variants); gc()
direction <- rbindlist(signed_summaries)
stopifnot(max(abs(direction$posterior_total - 1)) < 1e-5)
pairs <- merge(pairs, direction, by = c("gwas_name", "ensembl", "signal_pair_index"))
setnames(pairs, "PP.H4.abf", "signal_PP_H4")
pairs[, signal_direction := fifelse(n_gwas_signals == 1L & n_eqtl_signals == 1L,
                                     marginal_direction, NA_integer_)]
pairs[, direction_reason := fcase(n_gwas_signals > 1L | n_eqtl_signals > 1L,
  "conditional_effects_unavailable_multiple_signals", is.na(marginal_direction),
  "less_than_95pct_original_mass_in_one_resolved_direction", default =
  "single_exported_signal_pair_marginal_direction_95pct_mass")]
pairs <- merge(pairs, eligible[, .(gwas_name, ensembl, ancestry, trait, tier, trait_group,
  aggregate_PP_H4_susie = PP.H4.susie, aggregate_PP_H4_abf = PP.H4.abf)],
  by = c("gwas_name", "ensembl"), all.x = TRUE)
stopifnot(!anyNA(pairs$trait_group))
# Preserve supported gene-trait records without an exported signal posterior.
absent <- eligible[!pairs, on = .(gwas_name, ensembl)]
absent[, `:=`(signal_pair_index = NA_integer_, signal_PP_H4 = NA_real_,
  signal_direction = NA_integer_, marginal_direction = NA_integer_,
  direction_reason = "signal_posterior_not_exported_for_this_pair")]
signals <- rbindlist(list(pairs, absent), fill = TRUE)
signals[, signal_supported := is.finite(signal_PP_H4) & signal_PP_H4 > .5]
signals <- merge(signals, gwreg[, .(gwas_name = study_name, gwas_n = N_tot,
  gwas_n_cases = N_cases)], by = "gwas_name", all.x = TRUE)
write(signals, "genetic_signal_direction.tsv.gz")
write(rbindlist(allele_qc), "allele_harmonization_qc.tsv")

message("Connecting disease-state and assay-native evidence")
bulk <- read(file.path(ml, "bulk/continuum_membership.tsv.gz"))
stopifnot(nrow(bulk) == 23370L, !anyDuplicated(bulk$gene_id_base))
setnames(bulk, "gene_id_base", "ensembl")
catalog <- merge(bulk, tested, by = "ensembl", all = TRUE)
names_from_genetics <- unique(g[, .(ensembl, genetic_gene_name = gene)], by = "ensembl")
catalog <- merge(catalog, names_from_genetics, by = "ensembl", all.x = TRUE)
catalog[is.na(gene_name) | !nzchar(gene_name), gene_name := genetic_gene_name]
catalog[, genetic_evidence_state := fcase(is.na(n_gene_trait_tests), "not_tested",
  primary_susie_supported, "main_trait_susie_supported", primary_abf_supported,
  "main_trait_abf_sensitivity_supported", n_susie_testable > 0,
  "tested_without_main_trait_support", default = "susie_untestable_abf_tested")]
catalog[, remodeling_state := fcase(is.na(membership_class), "bulk_unobservable",
  axis_component == TRUE, "axis_component_separate_LOO_inference",
  continuum_associated == TRUE, "continuum_associated", default = "tested_without_continuum_support")]
prot <- read(file.path(protein_root, "protein_de_adjusted.tsv"))
stopifnot(!anyDuplicated(prot$gene))
prot[, `:=`(protein_n_samples = 58L, protein_n_control = 12L, protein_n_masld = 46L,
  protein_assay = "PXD051911_unpaired_liver_DIA_MS",
  protein_effect_unit = "log2_abundance_MASLD_minus_control_adjusted_batch_age_BMI_sex")]
setnames(prot, c("gene", "logFC", "padj", "pvalue"),
  c("gene_name", "protein_log2FC", "protein_BH_q", "protein_p"))
catalog <- merge(catalog, prot[, .(gene_name, protein_log2FC, protein_BH_q, protein_p,
  protein_n_samples, protein_n_control, protein_n_masld, protein_assay, protein_effect_unit)],
  by = "gene_name", all.x = TRUE)
catalog[, `:=`(protein_state = fcase(is.na(protein_BH_q), "unobservable_in_source_protein_table",
  protein_BH_q < .05, "tissue_protein_disease_association", default = "measured_without_BH_support"),
  tissue_blood_coupling = "untestable_no_verified_paired_measurements",
  target_readout_identity = "same_gene_tissue_protein_if_measured_otherwise_not_established",
  patient_benefit_prediction = "not_tested")]
write(catalog, "gene_catalog_translation.tsv.gz")
records <- merge(signals, catalog[, .(ensembl, gene_name, remodeling_state,
  beta_fixed_projection, beta_signature_pc1, protein_state, protein_log2FC, protein_BH_q,
  tissue_blood_coupling)], by = "ensembl", all.x = TRUE)
records[, relationship := fcase(!signal_supported, "supported_aggregate_or_sensitivity_without_supported_signal",
  is.na(signal_direction), "direction_unresolved",
  is.na(remodeling_state) | remodeling_state == "bulk_unobservable", "genetics_bulk_unobservable",
  remodeling_state == "axis_component_separate_LOO_inference", "axis_component_use_separate_LOO",
  remodeling_state != "continuum_associated", "genetics_without_supported_bulk_remodeling",
  sign(beta_fixed_projection) == signal_direction, "same_direction",
  default = "opposite_direction")]
records[, perturbation_hypothesis := fcase(relationship == "same_direction",
  "test_reversal_of_expression_association_with_bidirectional_perturbation",
  relationship == "opposite_direction",
  "test_bidirectionally_expression_reversal_alone_may_mislead_no_compensation_claim",
  relationship == "genetics_without_supported_bulk_remodeling",
  "retain_inherited_candidate_test_context_and_direction_without_DEG_filter",
  default = "resolve_variant_transcript_and_direction_before_directional_nomination")]
records[trait_group %in% c("liver_enzyme", "other_trait"), perturbation_hypothesis :=
  paste0("trait_specific_not_MASLD_risk;", perturbation_hypothesis)]
records[, `:=`(candidate_molecular_object = "colocalized_transcript_not_proven_causal_target",
  cell_context = "see_signal_specific_ATAC_and_gene_program_links_no_transfer_of_program_support",
  readout = fifelse(protein_state == "tissue_protein_disease_association",
    "same_gene_tissue_protein_candidate_unpaired_with_bulk", "tissue_RNA_no_established_accessible_readout"),
  unresolved_experiment = "resolve_regulatory_variant_and_transcript_then_test_bidirectional_perturbation_and_independent_readout",
  trial_design_hypothesis = "record_continuous_program_scores_within_histologic_eligibility_groups_no_treatment_response_claim")]
write(records, "therapeutic_hypothesis_records.tsv.gz")

# Native tables preserve their own denominators, effects and testing families.
write(read(file.path(ctx, "genetics/context/genetic_lineage_context_all_pairs.tsv")),
  "signal_ATAC_context.tsv")
write(read(file.path(ctx, "programs/program_atac_results.tsv")), "program_ATAC_native.tsv")
write(read(file.path(sp, "coverage/spatial_program_coverage.tsv")), "program_spatial_coverage.tsv")
write(read(file.path(sp, "effects/spatial_program_effects.tsv")), "program_spatial_effects_native.tsv")
write(read(file.path(ml, "paired/paired_program_results.tsv")), "program_paired_biopsy_native.tsv")
write(read(file.path(ml, "paired/paired_participant_manifest.tsv")), "paired_biopsy_denominators.tsv")
write(read(file.path(ml, "paired/paired_gene_results.tsv")), "gene_paired_biopsy_native.tsv")
write(read(file.path(protein_root, "module_protein_results.tsv")), "legacy_subset_program_protein_native.tsv")
write(prot, "gene_protein_native.tsv")
windows <- read("figures/candidates/histology-continuum-five-cohort-context-20260824T175647Z/source_tables/all117_hotspot_five_cohort_windows.tsv.gz")
stopifnot(uniqueN(windows$program_uid) == 117L, uniqueN(windows$window_id) == 9L,
  all(windows$visualization_only), !anyDuplicated(windows[, .(program_uid, dataset, window_id)]))
windows[, source_role := fifelse(dataset %in% c("GSE162694", "GSE213621"),
  "previously_examined_evaluation_cohort", "source_overlap_descriptive_only")]
write(windows, "frozen_nine_window_profiles.tsv.gz")
write(data.table(path = inputs), "input_paths.tsv")
capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
write(data.table(check = c("frozen_effects_reproduced", "complete_program_family",
  "all_signal_pairs_retained", "all_variant_instances_retained", "genetic_aggregate_hash_matches_replay",
  "complete_bulk_family", "nine_descriptive_windows"),
  value = c(max(check$absolute_beta_difference, na.rm = TRUE), nrow(registry), nrow(pairs),
            n_variants_written, 1, nrow(bulk), uniqueN(windows$window_id))), "build_checks.tsv")
message("Complete evidence build: ", normalizePath(out))
