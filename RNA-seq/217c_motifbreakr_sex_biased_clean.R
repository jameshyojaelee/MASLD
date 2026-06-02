#!/usr/bin/env Rscript
# 217c_motifbreakr_sex_biased_clean.R -- motifbreakR cross-tab with clean target sets
#
# Agent: B5b (C5 critique fix)
#
# CHANGES vs 217c:
#   - Use clean target sets (217a_clean) implicitly via TF panel; motifbreakR
#     itself is independent of disease_regulons but motif_in_disease_regulon flag
#     IS kept for descriptive context (clearly labeled as descriptive, not
#     a filter).
#   - Adds negative-control TFs (CTCF/MYC/RFX5) to cross-tab.
#   - Fixes 217c bug: max_pip column collision when joining motif with varann
#     (both tables have max_pip). Use a rename in the join.
#
# Inputs (unchanged):
#   - GWAS/finemapping/results/gwas_atac/motif_disruption_scores.csv
#   - GWAS/finemapping/results/gwas_atac/gwas_atac_variant_annotation.csv
#   - RNA-seq/Human/.../results/integration/sex_deg_classification.csv
#   - GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv
#
# Outputs:
#   - RNA-seq/results/stratified_causal/motifbreakr_sex_biased_coloc_variants_clean.csv
#   - RNA-seq/results/stratified_causal/motifbreakr_sex_biased_summary_clean.csv

suppressPackageStartupMessages({
  library(data.table)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
outdir <- file.path(BASE, "RNA-seq/results/stratified_causal")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

HORMONE_TF_PANEL <- c("AR", "ESR1", "ESR2", "FOXA1", "FOXA2",
                      "STAT5A", "STAT5B", "BCL6", "CUX2", "HNF4A")
NEG_CTRL_TFS <- c("CTCF", "MYC", "RFX5")
TF_PANEL <- c(HORMONE_TF_PANEL, NEG_CTRL_TFS)

cat("=== 217c_clean: motifbreakR x sex-biased COLOC loci (with negative controls) ===\n")

# ===========================================================================
# 1. Load motif disruption + variant annotation
# ===========================================================================
motif <- fread(file.path(BASE,
  "GWAS/finemapping/results/gwas_atac/motif_disruption_scores.csv"))
varann <- fread(file.path(BASE,
  "GWAS/finemapping/results/gwas_atac/gwas_atac_variant_annotation.csv"))
cat("  motifbreakR variant-TF rows:", nrow(motif), "\n")
cat("  variant annotation rows:", nrow(varann), "\n")

# ===========================================================================
# 2. Sex DEG classification
# ===========================================================================
# v3 mashr Bayesian preferred; v2 fallback. v3 CSV provides v2-compatible `sex_class` alias.
sex_v3_path <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/sex_v3/sex_deg_classification_v3.csv")
sex_v2_path <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/sex_deg_classification.csv")
sex_file <- if (file.exists(sex_v3_path)) sex_v3_path else sex_v2_path
cat("  Source:", basename(dirname(sex_file)), "/", basename(sex_file), "\n")
sex_degs <- fread(sex_file)
sex_degs[, ensembl_id := sub("\\..*", "", gene)]

atlas <- fread(file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"),
               select = c("ensembl_id", "human_symbol"))
sex_degs <- merge(sex_degs, atlas, by = "ensembl_id", all.x = TRUE)

sex_biased_symbols <- sex_degs[sex_class %in% c("Female_biased", "Male_biased",
                                                  "Divergent",
                                                  "Female_specific", "Male_specific"),
                                human_symbol]
sex_biased_symbols <- unique(sex_biased_symbols[!is.na(sex_biased_symbols)
                                                  & sex_biased_symbols != ""])
cat("  Sex-biased symbols:", length(sex_biased_symbols), "\n")

# ===========================================================================
# 3. COLOC PP4 >= 0.5
# ===========================================================================
coloc <- fread(file.path(BASE,
  "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv"))
if ("coloc_best_susie_pp4" %in% names(coloc)) {
  coloc[, coloc_best_pp4 := fifelse(!is.na(coloc_best_susie_pp4),
                                    coloc_best_susie_pp4, coloc_best_pp4)]
}
coloc <- coloc[gene != "" & !is.na(gene)]
coloc_hits <- unique(coloc[coloc_best_pp4 >= 0.5, gene])
cat("  COLOC PP4>=0.5 genes:", length(coloc_hits), "\n")

sex_biased_coloc_genes <- intersect(sex_biased_symbols, coloc_hits)
cat("  Sex-biased AND COLOC genes:", length(sex_biased_coloc_genes), "\n")

# ===========================================================================
# 4. Variants linked to sex-biased COLOC genes
# ===========================================================================
varann[, linked_any := pmax(
  as.integer(scenic_target_gene %in% sex_biased_coloc_genes),
  as.integer(linked_gene %in% sex_biased_coloc_genes),
  as.integer(nearest_gene %in% sex_biased_coloc_genes),
  na.rm = TRUE)]
sb_coloc_variants <- varann[linked_any == 1, variant_id]
cat("  Variants in CS for sex-biased COLOC genes:", length(sb_coloc_variants), "\n")

motif[, in_sb_coloc := SNP_id %in% sb_coloc_variants]

# ===========================================================================
# 5. Per-TF cross-tab (includes negative controls)
# ===========================================================================
variant_dt <- motif[in_sb_coloc == TRUE & tf_name %in% TF_PANEL]
# FIX: motif already has max_pip; varann also has max_pip -- collision.
# Drop max_pip from varann slice before merging (keep motif's max_pip).
varann_slim <- varann[, .(variant_id, scenic_target_gene, linked_gene,
                          nearest_gene, peak_chr, peak_start, peak_end,
                          cell_type, hep_da_logFC, hep_da_padj)]
variant_dt <- merge(variant_dt, varann_slim,
                    by.x = "SNP_id", by.y = "variant_id", all.x = TRUE)

sex_lookup <- sex_degs[, .(human_symbol, sex_class)]
setnames(sex_lookup, "human_symbol", "linked_gene")
variant_dt <- merge(variant_dt, sex_lookup,
                    by = "linked_gene", all.x = TRUE)

variant_dt[, is_negative_control := tf_name %in% NEG_CTRL_TFS]
variant_dt <- unique(variant_dt[, .(SNP_id, tf_name, is_negative_control,
                                     effect, alleleDiff,
                                     max_pip, scenic_target_gene, linked_gene,
                                     nearest_gene, sex_class, peak_chr,
                                     peak_start, peak_end, cell_type,
                                     hep_da_logFC, hep_da_padj,
                                     motif_in_disease_regulon,
                                     priority_score)])
variant_dt <- variant_dt[order(is_negative_control, tf_name, -max_pip)]
fwrite(variant_dt,
       file.path(outdir, "motifbreakr_sex_biased_coloc_variants_clean.csv"))
cat("  Wrote motifbreakr_sex_biased_coloc_variants_clean.csv (",
    nrow(variant_dt), "rows)\n")

# ===========================================================================
# 6. Per-TF summary
# ===========================================================================
summary_dt <- variant_dt[, .(
  is_negative_control = unique(is_negative_control),
  n_disrupted_variants = uniqueN(SNP_id),
  n_unique_genes = uniqueN(linked_gene[!is.na(linked_gene)]),
  genes_affected = paste(sort(unique(linked_gene[!is.na(linked_gene)])),
                          collapse = ";"),
  female_biased_genes = uniqueN(linked_gene[sex_class %in%
                                              c("Female_biased", "Female_specific")]),
  male_biased_genes = uniqueN(linked_gene[sex_class %in%
                                             c("Male_biased", "Male_specific")]),
  n_in_disease_regulon = sum(motif_in_disease_regulon, na.rm = TRUE),
  max_priority_score = if (all(is.na(priority_score))) NA_real_
                         else max(priority_score, na.rm = TRUE)),
  by = tf_name]

miss <- setdiff(TF_PANEL, summary_dt$tf_name)
if (length(miss) > 0) {
  summary_dt <- rbind(summary_dt,
                      data.table(tf_name = miss,
                                 is_negative_control = miss %in% NEG_CTRL_TFS,
                                 n_disrupted_variants = 0,
                                 n_unique_genes = 0, genes_affected = "",
                                 female_biased_genes = 0, male_biased_genes = 0,
                                 n_in_disease_regulon = 0,
                                 max_priority_score = NA_real_),
                      fill = TRUE)
}
summary_dt <- summary_dt[order(is_negative_control, -n_disrupted_variants)]
fwrite(summary_dt,
       file.path(outdir, "motifbreakr_sex_biased_summary_clean.csv"))
cat("  Wrote motifbreakr_sex_biased_summary_clean.csv (", nrow(summary_dt), "rows)\n")
print(summary_dt)

cat("\n=== 217c_clean complete ===\n")
