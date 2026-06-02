#!/usr/bin/env Rscript
# 05_integrate_atlas.R
# Integrates fine-mapping results into the multi-evidence atlas
# Maps fine-mapped variants to genes and adds finemapping columns
# Usage: Rscript 05_integrate_atlas.R

library(data.table)
library(dplyr)
library(readr)

BASE_DIR <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
FM_DIR <- file.path(BASE_DIR, "GWAS/finemapping")
setwd(FM_DIR)

cat("============================================================\n")
cat("Integrating fine-mapping results into multi-evidence atlas\n")
cat("============================================================\n")

# Load fine-mapping results
fm_file <- file.path(FM_DIR, "results/combined_finemapping.csv")
if (!file.exists(fm_file)) {
  stop("Combined finemapping results not found. Run 04_aggregate_results.R first.")
}
fm <- fread(fm_file)
cat("Loaded", nrow(fm), "fine-mapped variant entries\n")

# Load atlas
atlas_file <- file.path(BASE_DIR, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
if (!file.exists(atlas_file)) {
  stop("Multi-evidence atlas not found at ", atlas_file)
}
atlas <- fread(atlas_file)
cat("Loaded atlas:", nrow(atlas), "genes x", ncol(atlas), "columns\n")

# Load gene annotations for variant-to-gene mapping
# Use Broadaway eQTL leads for eQTL-informed mapping
eqtl_file <- file.path(BASE_DIR, "data/broadaway_eqtl/Liver_eQTL_Meta_Leads_ST3_20240530.tsv")
if (file.exists(eqtl_file)) {
  eqtl <- fread(eqtl_file)
  cat("Loaded", nrow(eqtl), "Broadaway eQTL lead variants\n")
} else {
  eqtl <- data.table()
  cat("WARNING: Broadaway eQTL leads not found — using nearest gene only\n")
}

# Strategy: Map variants to genes via two methods:
# 1. eQTL target gene (if variant is an eQTL)
# 2. Nearest gene by position (within 500kb)

# Load gene positions from GENCODE (from atlas or biotype file)
biotype_file <- file.path(BASE_DIR, "streamlit_deg_explorer/data/ensembl_gene_biotypes.tsv.gz")
if (file.exists(biotype_file)) {
  gene_info <- fread(biotype_file)
  cat("Loaded", nrow(gene_info), "gene annotations\n")
} else {
  cat("WARNING: Gene biotype file not found — skipping nearest gene mapping\n")
  gene_info <- data.table()
}

# Method 1: eQTL-based mapping
cat("\n--- Mapping variants to eQTL target genes ---\n")
if (nrow(eqtl) > 0 && "POS" %in% colnames(eqtl)) {
  # Create variant key from eQTL leads
  eqtl_map <- eqtl %>%
    mutate(fm_key = paste(CHR, POS, sep = ":")) %>%
    select(fm_key, eqtl_gene = GeneSymbol, eqtl_ensembl = ENSG) %>%
    distinct(fm_key, .keep_all = TRUE)

  fm <- fm %>%
    mutate(fm_key = paste(chromosome, position, sep = ":")) %>%
    left_join(eqtl_map, by = "fm_key")

  n_eqtl_mapped <- sum(!is.na(fm$eqtl_gene))
  cat("  Mapped", n_eqtl_mapped, "/", nrow(fm), "variants to eQTL target genes\n")
}

# For variants without eQTL mapping, assign to nearest gene
# (simplified: use gene_symbol from atlas if available, otherwise skip)

# Aggregate to gene level
cat("\n--- Aggregating to gene level ---\n")

# For each gene, find finemapping evidence from:
# a) Variants mapped via eQTL
# b) Variants within gene body or ±10kb

# Start with eQTL-mapped variants
if ("eqtl_gene" %in% colnames(fm)) {
  # Use only converged SuSiE PIPs to prevent non-converged PIP=1.0 artifacts
  # from inflating gene-level max PIP. The susie_pip_clean column was set to
  # NA for non-converged loci by 04_aggregate_results.R; use it here instead
  # of raw susie_pip. Also use recommended_pip (which handles EAS preference
  # for CARMA) when available.
  has_clean <- "susie_pip_clean" %in% colnames(fm)
  has_recommended <- "recommended_pip" %in% colnames(fm)

  gene_fm <- fm %>%
    filter(!is.na(eqtl_gene)) %>%
    group_by(gene = eqtl_gene) %>%
    summarise(
      fm_max_pip = if (has_recommended) max(recommended_pip, na.rm = TRUE)
                   else max(max_pip, na.rm = TRUE),
      fm_max_susie_pip = if (has_clean) max(susie_pip_clean, na.rm = TRUE)
                         else max(susie_pip, na.rm = TRUE),
      fm_max_carma_pip = max(carma_pip, na.rm = TRUE),
      fm_concordant_pip = max(concordant_pip, na.rm = TRUE),
      fm_n_cs_variants = sum(either_in_cs, na.rm = TRUE),
      fm_n_concordant_cs = sum(both_in_cs, na.rm = TRUE),
      fm_n_gwas = n_distinct(study),
      fm_n_loci = n_distinct(locus),
      fm_best_study = study[which.max(max_pip)],
      fm_best_locus = locus[which.max(max_pip)],
      fm_best_variant = variant_id[which.max(max_pip)],
      .groups = "drop"
    ) %>%
    mutate(
      fm_concordant = fm_n_concordant_cs > 0
    )
} else {
  gene_fm <- data.table(gene = character(0))
}

cat("Gene-level finemapping evidence for", nrow(gene_fm), "genes\n")

# Remove existing fm_ columns to prevent duplication on re-run
fm_cols <- grep("^fm_", colnames(atlas), value = TRUE)
if (length(fm_cols) > 0) {
  cat("Removing", length(fm_cols), "existing fm_ columns from atlas (re-run detected)\n")
  atlas <- atlas[, !..fm_cols]
}

# Join to atlas
if ("gene_symbol" %in% colnames(atlas)) {
  atlas_updated <- atlas %>%
    left_join(gene_fm, by = c("gene_symbol" = "gene"))
} else if ("gene" %in% colnames(atlas)) {
  atlas_updated <- atlas %>%
    left_join(gene_fm, by = "gene")
} else {
  cat("WARNING: Cannot match gene names between atlas and finemapping\n")
  atlas_updated <- atlas
}

# Summary stats
n_with_fm <- sum(!is.na(atlas_updated$fm_max_pip))
n_high_pip <- sum(atlas_updated$fm_max_pip > 0.5, na.rm = TRUE)
n_concordant <- sum(atlas_updated$fm_concordant, na.rm = TRUE)

cat("\n=== Atlas integration summary ===\n")
cat("Genes with finemapping evidence:", n_with_fm, "/", nrow(atlas_updated), "\n")
cat("Genes with PIP > 0.5:", n_high_pip, "\n")
cat("Genes with concordant credible sets:", n_concordant, "\n")

# Write updated atlas
out_atlas <- file.path(BASE_DIR, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
fwrite(atlas_updated, out_atlas)
cat("\nUpdated atlas written to:", out_atlas, "\n")
cat("New dimensions:", nrow(atlas_updated), "x", ncol(atlas_updated), "\n")

# Also save gene-level finemapping table separately
fwrite(gene_fm, file.path(FM_DIR, "results/gene_level_finemapping.csv"))
cat("Gene-level finemapping table:", file.path(FM_DIR, "results/gene_level_finemapping.csv"), "\n")

cat("\n============================================================\n")
cat("Atlas integration complete\n")
cat("============================================================\n")
