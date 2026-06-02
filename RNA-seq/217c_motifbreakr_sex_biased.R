#!/usr/bin/env Rscript
# 217c_motifbreakr_sex_biased.R — motifbreakR cross-tab against sex-biased COLOC loci
#
# Agent: B5 (Sex-Hormone TF x COLOC) -- HEADLINE A7
#
# Goal: For variants in CS for sex-biased COLOC genes, count motif disruptions
# per sex-hormone TF (AR/ESR1/FOXA1/STAT5/...). Bridge GWAS-ATAC motif disruption
# (Script 56 output) to sex-biased causal architecture.
#
# Inputs:
#   - GWAS/finemapping/results/gwas_atac/motif_disruption_scores.csv (Script 56)
#   - GWAS/finemapping/results/gwas_atac/gwas_atac_variant_annotation.csv
#   - RNA-seq/Human/.../results/integration/sex_deg_classification.csv
#   - GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv
#
# Outputs:
#   - RNA-seq/results/stratified_causal/motifbreakr_sex_biased_coloc_variants.csv
#   - RNA-seq/results/stratified_causal/motifbreakr_sex_biased_summary.csv

suppressPackageStartupMessages({
  library(data.table)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
outdir <- file.path(BASE, "RNA-seq/results/stratified_causal")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

TF_PANEL <- c("AR", "ESR1", "ESR2", "FOXA1", "FOXA2",
              "STAT5A", "STAT5B", "BCL6", "CUX2", "HNF4A")

cat("=== 217c: motifbreakR x sex-biased COLOC loci ===\n")

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
sex_file <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/sex_deg_classification.csv")
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
# 3. COLOC PP4 >= 0.5 genes
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

# Sex-biased AND COLOC
sex_biased_coloc_genes <- intersect(sex_biased_symbols, coloc_hits)
cat("  Sex-biased AND COLOC genes:", length(sex_biased_coloc_genes), "\n")
cat("  Examples: ", paste(head(sex_biased_coloc_genes, 10), collapse = ", "), "\n")

# ===========================================================================
# 4. Pick variants in CS for sex-biased COLOC genes
# ===========================================================================
# variant annotation has scenic_target_gene / linked_gene / nearest_gene
# Use the union to identify variants linked to sex-biased COLOC genes
varann[, linked_any := pmax(
  as.integer(scenic_target_gene %in% sex_biased_coloc_genes),
  as.integer(linked_gene %in% sex_biased_coloc_genes),
  as.integer(nearest_gene %in% sex_biased_coloc_genes),
  na.rm = TRUE)]
sb_coloc_variants <- varann[linked_any == 1, variant_id]
cat("  Variants in CS for sex-biased COLOC genes:", length(sb_coloc_variants), "\n")

# Build SNP_id with chr prefix to match motif$SNP_id format (chr:pos:ref:alt)
# motif$SNP_id format examples: "6:32191581:A:T" (no chr prefix)
# variant_id format examples: "1:16505320:A:G" (no chr prefix)
motif[, in_sb_coloc := SNP_id %in% sb_coloc_variants]

# ===========================================================================
# 5. Per-TF cross-tab + the variant-level table
# ===========================================================================
variant_dt <- motif[in_sb_coloc == TRUE & tf_name %in% TF_PANEL]
# Join back linked gene + sex_class
variant_dt <- merge(variant_dt,
                    varann[, .(variant_id, scenic_target_gene, linked_gene,
                                nearest_gene, max_pip, peak_chr,
                                peak_start, peak_end, cell_type,
                                hep_da_logFC, hep_da_padj)],
                    by.x = "SNP_id", by.y = "variant_id", all.x = TRUE)

# Add sex_class label based on which sex-biased gene the variant links to
sex_lookup <- sex_degs[, .(human_symbol, sex_class)]
setnames(sex_lookup, "human_symbol", "linked_gene")
variant_dt <- merge(variant_dt, sex_lookup,
                    by = "linked_gene", all.x = TRUE)

# Persist (deduplicate by SNP+tf)
variant_dt <- unique(variant_dt[, .(SNP_id, tf_name, effect, alleleDiff,
                                     max_pip, scenic_target_gene, linked_gene,
                                     nearest_gene, sex_class, peak_chr,
                                     peak_start, peak_end, cell_type,
                                     hep_da_logFC, hep_da_padj,
                                     motif_in_disease_regulon,
                                     priority_score)])
variant_dt <- variant_dt[order(tf_name, -max_pip)]
fwrite(variant_dt,
       file.path(outdir, "motifbreakr_sex_biased_coloc_variants.csv"))
cat("  Wrote motifbreakr_sex_biased_coloc_variants.csv (",
    nrow(variant_dt), "rows)\n")

# ===========================================================================
# 6. Per-TF summary
# ===========================================================================
summary_dt <- variant_dt[, .(
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
# Add zero rows for TFs in panel with no disruptions
miss <- setdiff(TF_PANEL, summary_dt$tf_name)
if (length(miss) > 0) {
  summary_dt <- rbind(summary_dt,
                      data.table(tf_name = miss, n_disrupted_variants = 0,
                                 n_unique_genes = 0, genes_affected = "",
                                 female_biased_genes = 0, male_biased_genes = 0,
                                 n_in_disease_regulon = 0,
                                 max_priority_score = NA_real_),
                      fill = TRUE)
}
summary_dt <- summary_dt[order(-n_disrupted_variants)]
fwrite(summary_dt,
       file.path(outdir, "motifbreakr_sex_biased_summary.csv"))
cat("  Wrote motifbreakr_sex_biased_summary.csv (", nrow(summary_dt), "rows)\n")
print(summary_dt)

cat("\n=== 217c complete ===\n")
