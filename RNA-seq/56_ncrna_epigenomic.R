#!/usr/bin/env Rscript
# =============================================================================
# 56_ncrna_epigenomic.R
# Module 5: Epigenomic Regulation of ncRNA Loci
#
# Integrates ATAC-seq DA peaks, SCENIC+ enhancer-gene links and regulon targets,
# and GWAS variant overlap for ncRNA genes. Compares ncRNA vs protein-coding
# epigenomic accessibility.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr, warn.conflicts = FALSE)
  library(GenomicRanges)
  library(rtracklayer)
})

select <- dplyr::select
filter <- dplyr::filter

# --- Configuration ---
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

ncrna_deg_path     <- file.path(BASE, "RNA-seq/results/ncrna/ncrna_deg_annotated.csv")
atlas_path         <- file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
da_peaks_path      <- file.path(BASE, "Analysis/ATAC/Human_Multiome/results/l8_annotated/scatac_da_gene_annotated.csv")
enhancer_path      <- file.path(BASE, "Analysis/ATAC/Human_Multiome/scenic_plus/enhancer_gene_links.csv")
hep_regulon_path   <- file.path(BASE, "Analysis/ATAC/Human_Multiome/scenic_plus/hepatocyte_regulons.csv")
disease_reg_path   <- file.path(BASE, "Analysis/ATAC/Human_Multiome/scenic_plus/disease_regulons.csv")
human_gtf_path     <- "/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"
out_dir            <- file.path(BASE, "RNA-seq/results/ncrna")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== Module 5: Epigenomic Regulation of ncRNA Loci ===\n")
cat("Start:", format(Sys.time()), "\n\n")

# =============================================================================
# 1. Load Data
# =============================================================================
cat("--- 1. Loading data ---\n")

ncrna <- fread(ncrna_deg_path)
atlas <- fread(atlas_path)
cat(sprintf("  ncRNA genes: %d\n", nrow(ncrna)))
cat(sprintf("  Atlas genes: %d\n", nrow(atlas)))

# Check inputs
input_files <- c(
  "DA peaks" = da_peaks_path,
  "Enhancer links" = enhancer_path,
  "Hepatocyte regulons" = hep_regulon_path,
  "Disease regulons" = disease_reg_path
)

for (nm in names(input_files)) {
  exists <- file.exists(input_files[nm])
  cat(sprintf("  %-25s [%s]\n", nm, ifelse(exists, "OK", "MISSING")))
}

# =============================================================================
# 2. Define ncRNA Promoter Regions
# =============================================================================
cat("\n--- 2. Defining ncRNA promoter regions ---\n")

# Load GTF for coordinates
gtf <- import(human_gtf_path)
gtf_genes <- gtf[gtf$type == "gene"]

# ncRNA biotypes
ncrna_biotypes <- c("lncRNA", "miRNA", "snoRNA", "snRNA", "misc_RNA", "scaRNA")

ncrna_gtf <- gtf_genes[gtf_genes$gene_type %in% ncrna_biotypes]
pc_gtf <- gtf_genes[gtf_genes$gene_type == "protein_coding"]

# Filter to genes in our dataset
ncrna_in_data <- ncrna_gtf[ncrna_gtf$gene_name %in% ncrna$human_symbol]
pc_in_data <- pc_gtf[pc_gtf$gene_name %in% atlas[gene_biotype == "protein_coding"]$human_symbol]

cat(sprintf("  ncRNA genes with coordinates: %d\n", length(ncrna_in_data)))
cat(sprintf("  PC genes with coordinates: %d\n", length(pc_in_data)))

# Define promoter regions (TSS +/- 2kb)
ncrna_promoters <- promoters(ncrna_in_data, upstream = 2000, downstream = 2000)
pc_promoters <- promoters(pc_in_data, upstream = 2000, downstream = 2000)

# =============================================================================
# 3. ATAC-seq DA Peak Overlap
# =============================================================================
cat("\n--- 3. ATAC-seq DA peak analysis ---\n")

if (file.exists(da_peaks_path)) {
  da_peaks <- fread(da_peaks_path)
  cat(sprintf("  DA peaks loaded: %d rows\n", nrow(da_peaks)))
  cat(sprintf("  Columns: %s\n", paste(names(da_peaks), collapse = ", ")))

  # Extract ncRNA-associated DA peaks from atlas columns
  # The atlas already has hepatocyte_da_logFC, hepatocyte_da_padj columns
  ncrna_da <- atlas[gene_biotype %in% ncrna_biotypes & !is.na(hepatocyte_da_logFC),
                     .(human_symbol, gene_biotype, hepatocyte_da_logFC, hepatocyte_da_padj,
                       human_promoter_accessible)]
  pc_da <- atlas[gene_biotype == "protein_coding" & !is.na(hepatocyte_da_logFC),
                  .(human_symbol, gene_biotype, hepatocyte_da_logFC, hepatocyte_da_padj,
                    human_promoter_accessible)]

  cat(sprintf("  ncRNAs with DA data: %d\n", nrow(ncrna_da)))
  cat(sprintf("  PCs with DA data: %d\n", nrow(pc_da)))

  # Compare accessible promoter fractions
  ncrna_accessible <- atlas[gene_biotype %in% ncrna_biotypes & !is.na(human_promoter_accessible)]
  pc_accessible <- atlas[gene_biotype == "protein_coding" & !is.na(human_promoter_accessible)]

  n_ncrna_acc <- sum(ncrna_accessible$human_promoter_accessible == TRUE, na.rm = TRUE)
  n_pc_acc <- sum(pc_accessible$human_promoter_accessible == TRUE, na.rm = TRUE)

  frac_ncrna_acc <- n_ncrna_acc / max(nrow(ncrna_accessible), 1)
  frac_pc_acc <- n_pc_acc / max(nrow(pc_accessible), 1)

  cat(sprintf("  Accessible promoters: ncRNA %.1f%% (%d/%d), PC %.1f%% (%d/%d)\n",
              100 * frac_ncrna_acc, n_ncrna_acc, nrow(ncrna_accessible),
              100 * frac_pc_acc, n_pc_acc, nrow(pc_accessible)))

  # Fisher's test for accessible promoter enrichment
  mat_acc <- matrix(c(
    n_ncrna_acc, nrow(ncrna_accessible) - n_ncrna_acc,
    n_pc_acc, nrow(pc_accessible) - n_pc_acc
  ), nrow = 2, byrow = TRUE)
  fisher_acc <- fisher.test(mat_acc)
  cat(sprintf("  Accessibility enrichment: OR=%.2f, p=%.2e\n",
              fisher_acc$estimate, fisher_acc$p.value))

  # Compare DA logFC distributions
  if (nrow(ncrna_da) > 5 & nrow(pc_da) > 5) {
    da_test <- wilcox.test(abs(ncrna_da$hepatocyte_da_logFC),
                           abs(pc_da$hepatocyte_da_logFC))
    cat(sprintf("  |DA logFC|: ncRNA median=%.3f, PC median=%.3f, p=%.2e\n",
                median(abs(ncrna_da$hepatocyte_da_logFC)),
                median(abs(pc_da$hepatocyte_da_logFC)),
                da_test$p.value))
  }

  # Fraction with significant DA
  ncrna_sig_da <- ncrna_da[hepatocyte_da_padj < 0.05]
  pc_sig_da <- pc_da[hepatocyte_da_padj < 0.05]
  cat(sprintf("  Significant DA: ncRNA %d (%.1f%%), PC %d (%.1f%%)\n",
              nrow(ncrna_sig_da), 100 * nrow(ncrna_sig_da) / max(nrow(ncrna_da), 1),
              nrow(pc_sig_da), 100 * nrow(pc_sig_da) / max(nrow(pc_da), 1)))

  # Build per-ncRNA DA summary
  ncrna_da_summary <- atlas[gene_biotype %in% ncrna_biotypes,
                             .(human_symbol,
                               promoter_accessible = human_promoter_accessible,
                               da_logFC = hepatocyte_da_logFC,
                               da_padj = hepatocyte_da_padj,
                               is_da = !is.na(hepatocyte_da_padj) & hepatocyte_da_padj < 0.05,
                               cell_type_specific_access)]
} else {
  cat("  SKIPPED: DA peaks file not found. Using atlas columns.\n")
  ncrna_da_summary <- atlas[gene_biotype %in% ncrna_biotypes,
                             .(human_symbol,
                               promoter_accessible = human_promoter_accessible,
                               da_logFC = hepatocyte_da_logFC,
                               da_padj = hepatocyte_da_padj,
                               is_da = !is.na(hepatocyte_da_padj) & hepatocyte_da_padj < 0.05,
                               cell_type_specific_access)]
}

# =============================================================================
# 4. SCENIC+ Enhancer-Gene Links
# =============================================================================
cat("\n--- 4. SCENIC+ enhancer-gene links ---\n")

n_enh_links <- 0
if (file.exists(enhancer_path)) {
  enhancers <- fread(enhancer_path)
  cat(sprintf("  Enhancer-gene links loaded: %d\n", nrow(enhancers)))

  # Identify columns
  cat(sprintf("  Columns: %s\n", paste(names(enhancers), collapse = ", ")))

  # Look for target gene column (may be 'target', 'gene', 'target_gene', etc.)
  gene_col <- intersect(names(enhancers), c("target_gene", "gene", "target", "Gene"))
  if (length(gene_col) > 0) {
    gene_col <- gene_col[1]
    # Filter for ncRNA targets
    ncrna_enhancer <- enhancers[get(gene_col) %in% ncrna$human_symbol]
    n_enh_links <- nrow(ncrna_enhancer)
    cat(sprintf("  Enhancer links to ncRNAs: %d\n", n_enh_links))

    # Per-ncRNA enhancer count
    if (n_enh_links > 0) {
      enh_per_gene <- ncrna_enhancer[, .N, by = gene_col]
      setnames(enh_per_gene, c("human_symbol", "n_enhancer_links"))
      ncrna_da_summary <- merge(ncrna_da_summary, enh_per_gene,
                                 by = "human_symbol", all.x = TRUE)
      ncrna_da_summary[is.na(n_enhancer_links), n_enhancer_links := 0L]

      cat(sprintf("  ncRNAs with enhancer links: %d\n",
                  sum(ncrna_da_summary$n_enhancer_links > 0)))
    }
  } else {
    cat("  WARNING: Could not identify target gene column\n")
  }
} else {
  cat("  SKIPPED: Enhancer links file not found\n")
}

if (!"n_enhancer_links" %in% names(ncrna_da_summary)) {
  ncrna_da_summary[, n_enhancer_links := 0L]
}

# =============================================================================
# 5. SCENIC+ Regulon Targets
# =============================================================================
cat("\n--- 5. SCENIC+ regulon targets ---\n")

regulon_tfs <- character(0)
if (file.exists(hep_regulon_path)) {
  hep_reg <- fread(hep_regulon_path)
  cat(sprintf("  Hepatocyte regulons loaded: %d\n", nrow(hep_reg)))
  cat(sprintf("  Columns: %s\n", paste(names(hep_reg), collapse = ", ")))

  # Use atlas regulon columns
  reg_from_atlas <- atlas[gene_biotype %in% ncrna_biotypes & !is.na(scenic_grn_target),
                           .(human_symbol, scenic_grn_target, scenic_regulon_tf,
                             scenic_regulon_activity_diff)]

  # C3 fix: exclude string "FALSE" from count
  n_reg_targets <- sum(!is.na(reg_from_atlas$scenic_grn_target) &
                        reg_from_atlas$scenic_grn_target != "" &
                        reg_from_atlas$scenic_grn_target != "FALSE" &
                        reg_from_atlas$scenic_grn_target != "false")
  cat(sprintf("  ncRNAs as regulon targets (from atlas): %d\n", n_reg_targets))

  # Merge regulon info
  # C3 fix: scenic_grn_target stores string "FALSE" (not logical FALSE),
  # so we must explicitly exclude it to avoid counting all genes as targets
  reg_info <- atlas[gene_biotype %in% ncrna_biotypes,
                     .(human_symbol,
                       regulated_by_tf = scenic_regulon_tf,
                       is_regulon_target = !is.na(scenic_grn_target) & scenic_grn_target != "" &
                                           scenic_grn_target != "FALSE" & scenic_grn_target != "false",
                       regulon_activity_diff = scenic_regulon_activity_diff)]

  ncrna_da_summary <- merge(ncrna_da_summary, reg_info, by = "human_symbol", all.x = TRUE)
} else {
  cat("  SKIPPED: Hepatocyte regulon file not found\n")
  ncrna_da_summary[, regulated_by_tf := NA_character_]
  ncrna_da_summary[, is_regulon_target := FALSE]
  ncrna_da_summary[, regulon_activity_diff := NA_real_]
}

# Disease regulons
if (file.exists(disease_reg_path)) {
  dis_reg <- fread(disease_reg_path)
  cat(sprintf("  Disease regulons loaded: %d\n", nrow(dis_reg)))

  # Check for ncRNA targets in disease regulons
  target_col <- intersect(names(dis_reg), c("target_gene", "gene", "target", "Gene"))
  if (length(target_col) > 0) {
    target_col <- target_col[1]
    ncrna_disease_targets <- dis_reg[get(target_col) %in% ncrna$human_symbol]
    cat(sprintf("  ncRNAs in disease regulons: %d\n", nrow(ncrna_disease_targets)))
  }
}

# =============================================================================
# 6. GWAS Variant Overlap
# =============================================================================
cat("\n--- 6. GWAS variant overlap for ncRNAs ---\n")

# Extract COLOC columns from atlas for ncRNAs
coloc_cols <- grep("coloc_pp4", names(atlas), value = TRUE)
cat(sprintf("  COLOC columns in atlas: %d\n", length(coloc_cols)))

ncrna_gwas <- atlas[gene_biotype %in% ncrna_biotypes,
                     c("human_symbol", coloc_cols), with = FALSE]

# Count ncRNAs with any COLOC signal
ncrna_gwas[, has_any_coloc := apply(.SD, 1, function(x) any(!is.na(x) & x > 0.5)),
            .SDcols = coloc_cols]
ncrna_gwas[, max_coloc_pp4 := apply(.SD, 1, max, na.rm = TRUE),
            .SDcols = coloc_cols]
ncrna_gwas[is.infinite(max_coloc_pp4), max_coloc_pp4 := NA_real_]

n_coloc <- sum(ncrna_gwas$has_any_coloc, na.rm = TRUE)
cat(sprintf("  ncRNAs with any COLOC PP.H4 > 0.5: %d\n", n_coloc))

# Merge GWAS info
ncrna_da_summary <- merge(ncrna_da_summary,
                           ncrna_gwas[, .(human_symbol, has_gwas_coloc = has_any_coloc,
                                          max_coloc_pp4)],
                           by = "human_symbol", all.x = TRUE)
ncrna_da_summary[is.na(has_gwas_coloc), has_gwas_coloc := FALSE]

# Check TWAS signal (MR tally removed 2026-04-22 — MR ditched from paper)
ncrna_twas <- atlas[gene_biotype %in% ncrna_biotypes & !is.na(twas_pval) & twas_pval < 0.05,
                     .(human_symbol)]
cat(sprintf("  ncRNAs with TWAS signal: %d\n", nrow(ncrna_twas)))

# =============================================================================
# 7. Cross-Species Promoter Conservation
# =============================================================================
cat("\n--- 7. Cross-species promoter conservation ---\n")

ncrna_xspecies <- atlas[gene_biotype %in% ncrna_biotypes,
                         .(human_symbol, cross_species_promoter_conserved)]
n_conserved <- sum(ncrna_xspecies$cross_species_promoter_conserved == TRUE, na.rm = TRUE)
cat(sprintf("  ncRNAs with cross-species promoter conservation: %d\n", n_conserved))

ncrna_da_summary <- merge(ncrna_da_summary,
                           ncrna_xspecies,
                           by = "human_symbol", all.x = TRUE)

# =============================================================================
# 8. Compile Epigenomic Summary
# =============================================================================
cat("\n--- 8. Compiling epigenomic summary ---\n")

# Ensure all columns exist
default_cols <- c("promoter_accessible", "is_da", "da_logFC", "da_padj",
                  "n_enhancer_links", "regulated_by_tf", "is_regulon_target",
                  "has_gwas_coloc", "max_coloc_pp4", "cross_species_promoter_conserved")
for (col in default_cols) {
  if (!col %in% names(ncrna_da_summary)) {
    ncrna_da_summary[, (col) := NA]
  }
}

fwrite(ncrna_da_summary, file.path(out_dir, "ncrna_epigenomic_summary.csv"))
cat(sprintf("  Saved ncrna_epigenomic_summary.csv (%d rows)\n", nrow(ncrna_da_summary)))

# Comparison table
atac_comparison <- data.table(
  metric = c("promoter_accessible_frac", "sig_DA_frac", "median_abs_DA_logFC",
             "has_enhancer_link_frac", "is_regulon_target_frac"),
  ncRNA = c(
    mean(ncrna_da_summary$promoter_accessible == TRUE, na.rm = TRUE),
    mean(ncrna_da_summary$is_da == TRUE, na.rm = TRUE),
    median(abs(ncrna_da_summary$da_logFC[!is.na(ncrna_da_summary$da_logFC)])),
    mean(ncrna_da_summary$n_enhancer_links > 0, na.rm = TRUE),
    mean(ncrna_da_summary$is_regulon_target == TRUE, na.rm = TRUE)
  ),
  protein_coding = c(
    mean(atlas[gene_biotype == "protein_coding"]$human_promoter_accessible == TRUE, na.rm = TRUE),
    mean(!is.na(atlas[gene_biotype == "protein_coding"]$hepatocyte_da_padj) &
         atlas[gene_biotype == "protein_coding"]$hepatocyte_da_padj < 0.05, na.rm = TRUE),
    median(abs(atlas[gene_biotype == "protein_coding" & !is.na(hepatocyte_da_logFC)]$hepatocyte_da_logFC)),
    mean(atlas[gene_biotype == "protein_coding"]$scenic_enhancer_link == TRUE, na.rm = TRUE),
    # C3 fix: exclude string "FALSE" from regulon target fraction
    mean(!is.na(atlas[gene_biotype == "protein_coding"]$scenic_grn_target) &
         atlas[gene_biotype == "protein_coding"]$scenic_grn_target != "" &
         atlas[gene_biotype == "protein_coding"]$scenic_grn_target != "FALSE" &
         atlas[gene_biotype == "protein_coding"]$scenic_grn_target != "false", na.rm = TRUE)
  )
)

fwrite(atac_comparison, file.path(out_dir, "ncrna_atac_comparison.csv"))
cat("  Saved ncrna_atac_comparison.csv\n")
print(atac_comparison)

# =============================================================================
# 9. Summary
# =============================================================================
cat("\n--- 9. Summary ---\n")
cat(sprintf("  ncRNA genes profiled: %d\n", nrow(ncrna_da_summary)))
cat(sprintf("  With accessible promoters: %d\n",
            sum(ncrna_da_summary$promoter_accessible == TRUE, na.rm = TRUE)))
cat(sprintf("  With significant DA: %d\n",
            sum(ncrna_da_summary$is_da == TRUE, na.rm = TRUE)))
cat(sprintf("  With enhancer links: %d\n",
            sum(ncrna_da_summary$n_enhancer_links > 0, na.rm = TRUE)))
cat(sprintf("  As regulon targets: %d\n",
            sum(ncrna_da_summary$is_regulon_target == TRUE, na.rm = TRUE)))
cat(sprintf("  With GWAS COLOC: %d\n",
            sum(ncrna_da_summary$has_gwas_coloc == TRUE, na.rm = TRUE)))

cat("\n=== Module 5 Complete ===\n")
cat("End:", format(Sys.time()), "\n")
