#!/usr/bin/env Rscript
# KEY MESSAGE: Genetic hypotheses and disease readouts have different evidence,
# and their support must remain visible across all frozen molecular programs.
source("scripts/analysis/histology_anchored_continuum/translation/lib_translation.R")
suppressPackageStartupMessages({library(ggplot2); library(patchwork)})
args <- commandArgs(TRUE); stopifnot(length(args) %in% 1:2)
root <- args[[1]]; out <- file.path(root, if (length(args) == 2L) args[[2]] else "integrated")
stopifnot(!dir.exists(out)); dir.create(out)
read <- function(name) fread(file.path(root, name))
write <- function(d, name) fwrite(d, file.path(out, name), sep = "\t", na = "NA")
registry <- read("frozen_program_registry.tsv")
members <- read("gene_program_links.tsv.gz")
catalog <- read("gene_catalog_translation.tsv.gz")
records <- read("therapeutic_hypothesis_records.tsv.gz")
summary <- read("program_composition_summary.tsv")
windows <- read("frozen_nine_window_profiles.tsv.gz")
atac <- read("program_ATAC_native.tsv")
spatial <- read("program_spatial_effects_native.tsv")
coverage <- read("program_spatial_coverage.tsv")
paired <- read("program_paired_biopsy_native.tsv")
protein <- read("gene_protein_native.tsv")
coupling <- read("paired_protein/paired_tissue_plasma_effects.tsv")
coupling_primary <- coupling[model == "stage_adjusted"]
stopifnot(!anyDuplicated(coupling_primary$gene_name))

# Distinguish main-trait untestability from absence of main-trait support.
# Other-trait tests cannot supply a missing main-trait denominator.
gtest <- fread("GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv",
  select = c("gwas_name", "ensembl", "PP.H4.susie", "PP.H4.abf"))
tiers <- fread("GWAS/finemapping/config/gwas_trait_tier.tsv")
gtest <- merge(gtest, tiers[, .(gwas_name = study_name, tier)], by = "gwas_name")
gtest <- gtest[nzchar(ensembl) & !is.na(ensembl) & tier <= 2,
  .(n_main_gene_trait_tests = .N, n_main_susie_testable = sum(is.finite(PP.H4.susie)),
    n_main_abf_testable = sum(is.finite(PP.H4.abf))), by = ensembl]
catalog <- merge(catalog, gtest, by = "ensembl", all.x = TRUE)
catalog[, genetic_evidence_state := fcase(is.na(n_main_gene_trait_tests), "main_traits_not_tested",
  primary_susie_supported, "main_trait_susie_supported", primary_abf_supported,
  "main_trait_abf_sensitivity_supported", n_main_susie_testable > 0,
  "main_traits_tested_without_support", default = "main_susie_untestable_abf_tested")]
loo <- fread("RNA-seq/results/histology_anchored_continuum/molecular_layers/hac-molecular-layers-20260818T173348Z/bulk/signature_loo_gene_meta.tsv")
write(loo, "axis_component_LOO_native.tsv")
loo_summary <- loo[, .(LOO_both_axes_supported = .N == 2L &&
  all(is.finite(bh_q_value) & bh_q_value < .05 & direction_concordant) && uniqueN(sign(beta)) == 1L,
  LOO_fixed_beta = beta[axis_id == "fixed_projection"],
  LOO_signature_beta = beta[axis_id == "signature_pc1"], LOO_family = unique(bh_family_size)),
  by = .(ensembl = gene_id_base)]
catalog <- merge(catalog, loo_summary, by = "ensembl", all.x = TRUE)

# Expand the Gene Catalog by a sidecar, preserving one stable gene identifier.
catalog[, tissue_blood_coupling := NULL]
catalog <- merge(catalog, coupling_primary[, .(gene_name,
  tissue_blood_coupling = evidence_state, paired_protein_beta = beta,
  paired_protein_q = BH_q, paired_protein_n = n_participants,
  paired_protein_family = family_size)], by = "gene_name", all.x = TRUE)
catalog[is.na(tissue_blood_coupling), tissue_blood_coupling := "unobservable_in_shared_protein_family"]
links <- members[, .(program_uids = paste(sort(unique(program_uid)), collapse = ";")), by = gene_name]
catalog <- merge(catalog, links, by = "gene_name", all.x = TRUE)
catalog[, `:=`(candidate_object = "measured_or_colocalized_transcript_identity_requires_mechanistic_resolution",
  hypothesis = fcase(genetic_evidence_state == "main_trait_susie_supported",
    "retain_signal_specific_inherited_hypothesis_see_therapeutic_records",
    remodeling_state == "continuum_associated", "state_readout_or_downstream_hypothesis_without_supported_main_SuSiE",
    default = "no_supported_translational_nomination"),
  cell_context = "program_source_lineage_and_signal_ATAC_context_are_distinct_evidence",
  measurable_readout = fcase(tissue_blood_coupling == "matched_tissue_plasma_association",
    "same_gene_plasma_protein_candidate_coupled_to_tissue_protein_not_bulk_RNA",
    protein_state == "tissue_protein_disease_association", "same_gene_tissue_protein_candidate",
    default = "RNA_measurement_only_or_unobservable_see_assay_states"),
  unresolved_experiment = "test_transcript_mechanism_and_predefined_readout_in_relevant_cells_or_tissue",
  treatment_response = "not_tested", bulk_RNA_protein_participant_link = "unpaired")]
write(catalog, "gene_catalog_translation_complete.tsv.gz")

records[, tissue_blood_coupling := NULL]
records <- merge(records, catalog[, .(ensembl, tissue_blood_coupling,
  paired_protein_beta, paired_protein_q, paired_protein_n, measurable_readout, program_uids)],
  by = "ensembl", all.x = TRUE)
records[, readout := measurable_readout]
consistency <- records[signal_supported == TRUE & !is.na(signal_direction),
  .(n_studies_with_direction = uniqueN(gwas_name), n_expression_directions = uniqueN(signal_direction),
    study_directions = paste(sort(unique(paste(gwas_name, signal_direction, sep = ":"))), collapse = ";")),
  by = .(ensembl, trait, ancestry)]
consistency[, across_study_direction_status := fcase(n_expression_directions > 1L,
  "conflicting_directions_retained_no_common_perturbation_direction",
  n_studies_with_direction > 1L, "consistent_recorded_directions_source_overlap_not_adjudicated",
  default = "one_study_no_cross_study_replication")]
write(consistency, "genetic_direction_across_studies.tsv")
records <- merge(records, consistency, by = c("ensembl", "trait", "ancestry"), all.x = TRUE)
records[across_study_direction_status == "conflicting_directions_retained_no_common_perturbation_direction",
  perturbation_hypothesis := "resolve_conflicting_study_directions_and_signal_identity_before_common_directional_nomination"]
write(records, "therapeutic_hypothesis_records_complete.tsv.gz")

# All members with observable protein, rather than a selected downstream marker.
# Shared membership is NOT evidence that a target controls a readout.
readouts <- merge(members, protein[, .(gene_name, protein_log2FC, protein_BH_q)],
                  by = "gene_name", all.x = TRUE)
readouts <- merge(readouts, coupling_primary[, .(gene_name, paired_protein_beta = beta,
  paired_protein_q = BH_q, paired_protein_n = n_participants, coupling_state = evidence_state)],
  by = "gene_name", all.x = TRUE)
readouts[, `:=`(target_to_readout_regulation = "not_tested_shared_program_membership_only",
  protein_observability = fifelse(is.finite(protein_BH_q), "measured", "unobservable"))]
readouts[is.na(coupling_state), coupling_state := "unobservable_in_shared_protein_family"]
write(readouts, "program_member_readout_hypotheses.tsv.gz")

# Keep full coverage even where no formal spatial test was undertaken.
assays <- merge(coverage[, .(program_uid, dataset_id, assay_id, coverage_status,
  n_genes_measured, retained_l1_weight, biological_unit, biological_unit_resolution,
  source_dependence)], spatial[, .(program_uid, dataset_id, assay_id, evidence_state,
  estimate, pvalue, qvalue, n_biological, n_technical, estimand, effect_unit,
  multiplicity_family)], by = c("program_uid", "dataset_id", "assay_id"), all.x = TRUE)
assays[is.na(evidence_state), evidence_state := fifelse(coverage_status == "observable",
  "not_tested", "untestable")]
write(assays, "complete_program_spatial_evidence.tsv")

# Native single-cell lineage labels describe the discovery context. They do not
# assign transcript origin after ambient RNA correction or provide replication.
rna <- registry[, .(program_uid, cell_context = cell_type,
  source_RNA_state = "frozen_source_lineage_label_not_corrected_transcript_origin",
  single_cell_evidence_unit = "donor_in_source_tables",
  within_stage_cell_origin = "not_tested",
  ambient_qualification = fifelse(cell_type %in% c("t_cells", "tcells"),
    "native_noninteger_T_cell_transport_requires_source_specific_handling",
    "ambient_origin_not_resolved_by_bulk_composition_adjustment"))]
write(rna, "program_single_cell_context.tsv")

# Count distinct signals for genetic summaries; memberships never turn an
# inherited effect into a program-level genetic test.
supported <- records[tier <= 2 & signal_supported == TRUE]
direction_counts <- supported[, .(n_signals = .N, n_genes = uniqueN(ensembl)),
  by = .(trait_group, relationship)]
write(direction_counts, "genetic_direction_figure_counts.tsv")
gene_relationships <- unique(supported[, .(gene_name, relationship, trait_group)])
program_genetics <- merge(members[, .(program_uid, gene_name)], gene_relationships,
  by = "gene_name", allow.cartesian = TRUE)
program_genetics <- program_genetics[, .(n_genes = uniqueN(gene_name)),
  by = .(program_uid, trait_group, relationship)]
write(program_genetics, "program_member_genetic_relationship_counts.tsv")

feature <- registry[, .(program_uid, cell_type, module, module_name)]
for (arm in c("stage_sex", "stage_sex_composition")) {
  z <- summary[model == arm, .(program_uid, supported = both_axes_supported,
                              separate_cohort_BH = four_cohort_axis_BH & all_four_directions)]
  setnames(z, c("supported", "separate_cohort_BH"), paste0(arm, c("_supported", "_replicated")))
  feature <- merge(feature, z, by = "program_uid", all.x = TRUE)
}
pa <- paired[, .(paired_supported = .N == 2L && all(within_person_supported),
                 paired_n_min = min(n_pairs)), by = .(program_uid = feature_id)]
feature <- merge(feature, pa, by = "program_uid", all.x = TRUE)
as <- atac[, .(ATAC_observed = any(program_score_testable),
              ATAC_testable = any(contrast_testable), ATAC_replicated = any(replicated)), by = program_uid]
feature <- merge(feature, as, by = "program_uid", all.x = TRUE)
ss <- assays[, .(spatial_tested = any(is.finite(pvalue)),
  spatial_supported = any(evidence_state == "supported")), by = program_uid]
feature <- merge(feature, ss, by = "program_uid", all.x = TRUE)
pp <- readouts[, .(n_tissue_proteins_measured = sum(is.finite(protein_BH_q)),
  n_tissue_proteins_BH = sum(protein_BH_q < .05, na.rm = TRUE),
  n_matched_plasma_readouts = sum(coupling_state == "matched_tissue_plasma_association")), by = program_uid]
feature <- merge(feature, pp, by = "program_uid", all.x = TRUE)
# The display uses distinct gene counts (a gene may differ between trait strata).
feature[, n_resolved_genetic_genes := vapply(program_uid, function(uid) {
  symbols <- members[program_uid == uid, gene_name]
  uniqueN(gene_relationships[gene_name %in% symbols & relationship %in%
    c("same_direction", "opposite_direction"), gene_name])
}, integer(1))]
setorder(feature, cell_type, module)
feature[, row_label := paste0(gsub("_", " ", cell_type), " ", module, " | ", module_name)]
feature[, row_label := substr(row_label, 1, 67)]
feature[, row_index := .I]
write(feature, "all117_integrated_program_evidence.tsv")

# Equal-cohort descriptive profiles preserve the nine existing window definitions.
# The individual five-cohort profiles and participant counts remain in source TSV.
w <- windows[, .(mean_score = mean(mean_score, na.rm = TRUE),
  n_cohorts = sum(is.finite(mean_score)), n_participant_instances = sum(n_participants)),
  by = .(program_uid, window_id, window_center)]
w[!is.finite(mean_score), mean_score := NA_real_]
w <- merge(w, feature[, .(program_uid, row_label, row_index)], by = "program_uid")
write(w, "equal_cohort_descriptive_profiles.tsv")
source(path.expand("~/publication_color_themes.R"))
grDevices::pdf.options(useDingbats = FALSE)
base <- theme_classic(base_size = 6, base_family = "Helvetica") +
  theme(text = element_text(size = 6, face = "plain"),
    plot.title = element_text(size = 6, face = "plain"),
    strip.text = element_text(size = 6), legend.title = element_text(size = 6),
    legend.text = element_text(size = 6), axis.text = element_text(size = 6, colour = "black"))
relations <- c("same_direction", "opposite_direction", "genetics_without_supported_bulk_remodeling",
  "axis_component_use_separate_LOO", "genetics_bulk_unobservable", "direction_unresolved")
counts <- merge(CJ(trait_group = c("direct_MASLD", "liver_fat", "liver_enzyme"),
                   relationship = relations), direction_counts,
                by = c("trait_group", "relationship"), all.x = TRUE)
counts[is.na(n_signals), `:=`(n_signals = 0L, n_genes = 0L)]
counts[, relationship := factor(relationship, levels = rev(relations), labels = rev(c(
  "Same direction", "Opposite direction", "No supported bulk remodeling",
  "Axis component: separate LOO", "Bulk unobservable", "Direction unresolved")))]
counts[, trait_group := factor(trait_group, levels = c("direct_MASLD", "liver_fat", "liver_enzyme"),
  labels = c("Direct MASLD", "Liver fat", "Liver enzymes"))]
p1 <- ggplot(counts, aes(n_signals, relationship)) +
  geom_point(size = 1, colour = "#1565C0") + geom_text(aes(label = n_signals),
    hjust = -0.4, size = 6 / ggplot2::.pt) + facet_wrap(~trait_group, scales = "free_x", nrow = 1) +
  scale_x_continuous(expand = expansion(mult = c(.03, .25))) +
  labs(x = "Supported signal pairs (genes may recur)", y = NULL,
       title = "a  Inherited direction and cross-sectional RNA remodeling") + base
p2 <- ggplot(w, aes(window_center, -row_index, fill = mean_score)) +
  geom_tile(width = .1, height = .94) +
  scale_fill_gradient2(low = "#1565C0", mid = "white", high = "#C9265E",
    midpoint = 0, limits = c(-2, 2), oob = scales::squish, na.value = "#D9D9D9",
    name = "Mean RNA score (SD)") +
  scale_y_continuous(breaks = -feature$row_index, labels = feature$row_label,
    expand = expansion(add = .5)) +
  scale_x_continuous(breaks = c(.1, .5, .9), expand = c(0, 0)) +
  labs(x = "Frozen cohort percentile windows", y = NULL,
       title = "b  All 117 programs: frozen lineage/module order") + base +
  theme(axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        legend.position = "bottom", legend.key.width = grid::unit(.2, "in"))
tiles <- melt(feature[, .(row_index, `Stage + sex` = stage_sex_supported,
  `+ composition` = stage_sex_composition_supported,
  `Each cohort BH` = stage_sex_composition_replicated,
  `Paired biopsy` = paired_supported,
  `ATAC observed` = ATAC_observed, `ATAC replicated` = ATAC_replicated,
  `Spatial tested` = spatial_tested, `Spatial support` = spatial_supported)],
  id.vars = "row_index", variable.name = "test", value.name = "supported")
tiles[, value := fifelse(is.na(supported), "Unresolved", fifelse(supported, "Yes", "No"))]
unestimable <- read("program_cohort_effects.tsv")[, .(all_estimable = all(estimable)), by = program_uid]
unestimable <- merge(unestimable, feature[, .(program_uid, row_index)], by = "program_uid")
tiles[row_index %in% unestimable[all_estimable == FALSE, row_index] &
  test %in% c("Stage + sex", "+ composition", "Each cohort BH"), value := "Untestable"]
tiles[row_index %in% feature[ATAC_testable == FALSE, row_index] & test == "ATAC replicated", value := "Untestable"]
tiles[row_index %in% feature[spatial_tested == FALSE, row_index] & test == "Spatial support", value := "Not tested"]
p3 <- ggplot(tiles, aes(test, -row_index, fill = value)) + geom_tile(width = .9, height = .94) +
  scale_fill_manual(values = c(Yes = "#1565C0", No = "#EFEFEF", Unresolved = "#9E9E9E",
                              Untestable = "#9E9E9E", `Not tested` = "#D9D9D9"),
                    drop = FALSE, name = NULL) +
  scale_y_continuous(limits = c(-117.5, -.5), expand = c(0, 0)) +
  labs(x = NULL, y = NULL, title = "c  Distinct tests") + base +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    axis.line.y = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
    legend.position = "bottom") + guides(fill = guide_legend(nrow = 3))
numbers <- melt(feature[, .(row_index, `Genetic direction` = n_resolved_genetic_genes,
  `Tissue protein BH` = n_tissue_proteins_BH,
  `Matched plasma` = n_matched_plasma_readouts)], id.vars = "row_index",
  variable.name = "connection", value.name = "n_genes")
p4 <- ggplot(numbers, aes(connection, -row_index)) +
  geom_text(aes(label = n_genes), size = 6 / ggplot2::.pt) +
  scale_y_continuous(limits = c(-117.5, -.5), expand = c(0, 0)) +
  labs(x = NULL, y = NULL, title = "d  Member genes") + base +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
bottom <- p2 + p3 + p4 + plot_layout(widths = c(1.6, 1, .65))
combined <- p1 / bottom + plot_layout(heights = c(1.5, 12))
ggsave(file.path(out, "integrated_therapeutic_hypotheses.pdf"), combined,
       width = 9.5, height = 14.5, device = grDevices::pdf, useDingbats = FALSE)
ggsave(file.path(out, "panel_a_genetic_direction.pdf"), p1,
       width = 5.5, height = 2, device = grDevices::pdf, useDingbats = FALSE)
ggsave(file.path(out, "panel_b_program_profiles.pdf"), p2,
       width = 5.5, height = 13, device = grDevices::pdf, useDingbats = FALSE)
ggsave(file.path(out, "panel_c_distinct_tests.pdf"), p3,
       width = 2.5, height = 13, device = grDevices::pdf, useDingbats = FALSE)
ggsave(file.path(out, "panel_d_member_readouts.pdf"), p4,
       width = 2, height = 13, device = grDevices::pdf, useDingbats = FALSE)
writeLines(c("# Integrated figure caption", "",
  "Candidate extension; all 117 frozen programs are displayed in source-lineage/module order. Source lineage labels do not resolve ambient RNA or assign transcript origin.",
  "(a) Supported SuSiE signal pairs in the main trait tiers, separated into direct MASLD, liver fat, and liver enzymes. Distinct signals and ancestry contexts are retained; counts are not independent genes. Only a single exported signal pair with at least 95% of original shared posterior mass in one nonambiguous marginal direction is direction-resolved. Conditional effects for multiple signals are unavailable; no drug direction is established.",
  "(b) Equal-cohort averages of the frozen nine overlapping descriptive RNA windows, in SD units; colors clipped to +/-2. The full five-cohort values and participant counts remain in the source table. Percentiles are cohort-relative, not time or shared absolute severity. The three source-overlap cohorts are descriptive.",
  "(c) Stage/sex and composition columns require both co-primary meta-analysis BH q<0.05 and all four cohort/axis directions. Each-cohort BH is the stronger separate-cohort criterion under composition adjustment. Paired biopsy requires both existing within-person support calls. ATAC observability, replicated ATAC association, spatial test availability, and spatial support are different questions. A No in a support column can reflect lack of testing; inspect the adjacent testability column and complete native evidence tables. This is not a modality-voting score.",
  "(d) Counts of member genes with resolved genetic/remodeling relationships, disease-associated tissue protein, or matched tissue/plasma association. Shared membership does not establish target-to-marker regulation, gene-level spatial support, or independent confirmation of an RNA program. Matched plasma is 41 participants in a public source cohort, with BH across all 345 shared eligible proteins; bulk RNA remains unpaired. No treatment-response prediction, source-organ attribution, or clinical cutoff is tested."),
  file.path(out, "FIGURE_CAPTION.md"))

stopifnot(nrow(feature) == 117L, nrow(catalog) == uniqueN(catalog$ensembl),
  nrow(records) == nrow(read("therapeutic_hypothesis_records.tsv.gz")),
  all(catalog$bulk_RNA_protein_participant_link == "unpaired"),
  all(readouts$target_to_readout_regulation == "not_tested_shared_program_membership_only"),
  all(coupling$family_size == 345L), uniqueN(w$window_id) == 9L)
write(data.table(check = c("all_117_programs", "unique_gene_catalog", "signal_records_preserved",
  "no_unpaired_bulk_coupling", "no_target_readout_regulatory_claim", "complete_protein_family"),
  passed = TRUE), "integration_checks.tsv")
capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
cat("Integrated figure and complete evidence tables written to", out, "\n")
