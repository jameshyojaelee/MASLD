#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(TRUE); stopifnot(length(args) == 2L)
root <- args[[1]]; integrated <- file.path(root, args[[2]])
r <- function(x) fread(file.path(root, x))
signals <- r("genetic_signal_direction.tsv.gz")
variant_n <- 0L
for (file in list.files(file.path(root, "signed_variants"), full.names = TRUE)) {
  v <- fread(file)
  stopifnot(!anyDuplicated(v[, .(gwas_name, ensembl, signal_pair_index, chr, position)]))
  direct <- v[allele_state == "direct"]
  swapped <- v[allele_state == "swapped"]
  stopifnot(all(direct$allele1 == direct$eqtl_ea), all(direct$allele2 == direct$eqtl_nea),
    all(direct$aligned_eqtl_beta == direct$eqtl_beta),
    all(swapped$allele1 == swapped$eqtl_nea), all(swapped$allele2 == swapped$eqtl_ea),
    all(swapped$aligned_eqtl_beta == -swapped$eqtl_beta),
    all(is.na(v[!allele_state %in% c("direct", "swapped"), marginal_expression_direction])))
  z <- v[, .(pos = sum(SNP.PP.H4[which(gwas_beta * aligned_eqtl_beta > 0)]),
             neg = sum(SNP.PP.H4[which(gwas_beta * aligned_eqtl_beta < 0)])),
         by = .(gwas_name, ensembl, signal_pair_index)]
  z <- merge(z, signals, by = c("gwas_name", "ensembl", "signal_pair_index"))
  stopifnot(nrow(z) == uniqueN(v[, .(gwas_name, ensembl, signal_pair_index)]),
    max(abs(z$pos - z$positive_mass)) < 1e-10, max(abs(z$neg - z$negative_mass)) < 1e-10)
  variant_n <- variant_n + nrow(v)
}
stopifnot(variant_n == 22356697L,
  all(is.na(signals[n_gwas_signals > 1 | n_eqtl_signals > 1, signal_direction])),
  all(signals[!is.na(signal_direction), pmax(positive_mass, negative_mass)] >= .95))
fits <- r("program_cohort_effects.tsv")
denom <- r("program_denominators.tsv")
stopifnot(uniqueN(fits$program_uid) == 117,
  all(denom$n_complete_primary == denom$n_composition_matched))
families <- fits[, .(n_programs = .N, error = max(abs(q_value -
  p.adjust(p_value, "BH", n = 117)), na.rm = TRUE)), by = .(axis_id, model, dataset)]
stopifnot(all(families$n_programs == 117), all(families$error < 1e-12))
cat <- fread(file.path(integrated, "gene_catalog_translation_complete.tsv.gz"))
stopifnot(!anyDuplicated(cat$ensembl), nrow(cat) >= 23370)
protein <- r("paired_protein/paired_tissue_plasma_effects.tsv")
stopifnot(protein[, all(table(model) == 345)], all(protein$n_participants <= 41))
summary <- r("program_composition_summary.tsv")
stats <- data.table(metric = c("signed_variant_instances", "exported_signal_pairs",
  "supported_main_signal_pairs", "resolved_main_signal_pairs", "unique_catalog_genes",
  "stage_sex_supported_programs", "composition_supported_programs",
  "composition_separate_cohort_BH_programs", "matched_protein_family",
  "matched_protein_stage_adjusted_BH", "main_susie_supported_genes", "main_abf_supported_genes"),
  value = c(variant_n, sum(!is.na(signals$signal_pair_index)),
    signals[tier <= 2 & signal_supported == TRUE, .N],
    signals[tier <= 2 & signal_supported == TRUE & !is.na(signal_direction), .N], nrow(cat),
    summary[model == "stage_sex" & both_axes_supported == TRUE, .N],
    summary[model == "stage_sex_composition" & both_axes_supported == TRUE, .N],
    summary[model == "stage_sex_composition" & four_cohort_axis_BH & all_four_directions, .N],
    345, protein[model == "stage_adjusted" & BH_q < .05, .N],
    sum(cat$primary_susie_supported, na.rm = TRUE), sum(cat$primary_abf_supported, na.rm = TRUE)))
fwrite(stats, file.path(integrated, "verified_numbers.tsv"), sep = "\t")
print(stats)
print(unique(denom[n_complete_primary > 0, .(dataset, n_complete_primary)]))
print(fread(file.path(integrated, "therapeutic_hypothesis_records_complete.tsv.gz"))[
  tier <= 2 & signal_supported == TRUE & relationship %in% c("same_direction", "opposite_direction"),
  .(gene_name, trait_group, gwas_name, ancestry, relationship, signal_direction,
    beta_fixed_projection, positive_mass, negative_mass)])
writeLines("PASS: all signed variants independently checked; families, donor denominators, signal ambiguity, unique catalog and paired protein linkage checked",
  file.path(integrated, "VERIFICATION.txt"))
