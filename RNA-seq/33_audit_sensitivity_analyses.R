#!/usr/bin/env Rscript
# 33_audit_sensitivity_analyses.R
# ---------------------------------------------------------------------------
# Pre-advancement audit sensitivity analyses (Issues 5, 6, 7, 12)
#
# Part A: Sex inference validation on annotated samples (Issue 6)
# Part B: Permissive sex-divergent definition (Issue 6)
# Part C: Liver-specific CGP filter (Issue 7)
# Part D: Conserved non-GWAS enrichment (Issue 12)
#
# Outputs: RNA-seq/results/audit_sensitivity/
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
RDIR <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
OUTDIR <- file.path(BASE, "RNA-seq/results/audit_sensitivity")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

cat("=== Audit Sensitivity Analyses ===\n")
cat("Start time:", format(Sys.time()), "\n\n")

# ================================================================
# Part A: Sex Inference Validation (Issue 6)
# ================================================================
cat("--- Part A: Sex Inference Validation ---\n")

# Load the merged DGEList which contains sex_source and inferred_sex columns
dge_file <- file.path(RDIR, "merged_dge.rds")
if (file.exists(dge_file)) {
  dge <- readRDS(dge_file)
  meta <- as.data.table(dge$samples)
  cat("  Total samples in DGEList:", nrow(meta), "\n")
  cat("  Available columns:", paste(names(meta), collapse = ", "), "\n")

  # The DGEList should have sex_source ("annotated" vs "inferred_kmeans") and
  # inferred_sex and sex_final columns from Script 01
  has_validation <- "sex_source" %in% names(meta) && "inferred_sex" %in% names(meta)

  if (has_validation) {
    annotated <- meta[sex_source == "annotated"]
    cat("  Annotated sex samples:", nrow(annotated), "\n")
    cat("  Inferred sex samples:", sum(meta$sex_source == "inferred_kmeans"), "\n")

    if (nrow(annotated) > 0) {
      # For annotated samples, compare inferred_sex with reported sex
      # reported_sex should be in the metadata; sex_final uses it when available
      if ("reported_sex" %in% names(meta)) {
        both <- annotated[!is.na(reported_sex) & !is.na(inferred_sex)]
        concordant <- sum(both$reported_sex == both$inferred_sex)
        total <- nrow(both)
        accuracy <- concordant / total

        cat(sprintf("  Sex inference concordance: %d/%d (%.1f%%)\n",
                    concordant, total, 100*accuracy))

        # Per-dataset
        per_dataset <- both[, .(
          n = .N,
          concordant = sum(reported_sex == inferred_sex),
          accuracy = round(sum(reported_sex == inferred_sex) / .N, 3)
        ), by = dataset]
        cat("\n  Per-dataset concordance:\n")
        print(per_dataset)
        fwrite(per_dataset, file.path(OUTDIR, "sex_inference_per_dataset.csv"))

        validation_dt <- data.table(
          metric = c("total_samples", "concordant", "accuracy"),
          value = c(total, concordant, round(accuracy, 4))
        )
        fwrite(validation_dt, file.path(OUTDIR, "sex_inference_validation.csv"))
        cat("  Saved: sex_inference_validation.csv\n")
      } else {
        cat("  'reported_sex' column not found; sex concordance already computed in Script 01\n")
        cat("  Check logs/sex_check.pdf for visual validation\n")
      }
    }
  } else {
    cat("  DGEList lacks sex_source/inferred_sex columns.\n")
    cat("  Sex inference concordance was computed during Script 01 QC.\n")
    cat("  Proceeding to Part B.\n")
  }
} else {
  cat("  DGEList not found. Skipping sex validation.\n")
}

# ================================================================
# Part B: Permissive Sex-Divergent Definition (Issue 6)
# ================================================================
cat("\n--- Part B: Permissive Sex-Divergent Definition ---\n")

sex_class_file <- file.path(RDIR, "sex_deg_classification.csv")
if (file.exists(sex_class_file)) {
  sex_dt <- fread(sex_class_file)
  cat("  Total genes in sex classification:", nrow(sex_dt), "\n")

  # Stringent definition (current): BOTH sexes padj<0.1 with opposite signs
  stringent <- sex_dt[sex_class == "Sex_divergent"]
  cat("  Stringent sex-divergent:", nrow(stringent), "\n")

  # Permissive definition: |LFC_female - LFC_male| > 1.0 AND significant in ≥1 sex
  # Need to detect column names
  lfc_m_col <- grep("logFC_M|lfc_M|logFC_male", names(sex_dt), value = TRUE)[1]
  lfc_f_col <- grep("logFC_F|lfc_F|logFC_female", names(sex_dt), value = TRUE)[1]
  padj_m_col <- grep("padj_M|padj_male", names(sex_dt), value = TRUE)[1]
  padj_f_col <- grep("padj_F|padj_female", names(sex_dt), value = TRUE)[1]

  if (!is.na(lfc_m_col) && !is.na(lfc_f_col)) {
    sex_dt[, lfc_diff_abs := abs(get(lfc_f_col) - get(lfc_m_col))]

    sig_either <- !is.na(sex_dt[[padj_m_col]]) & sex_dt[[padj_m_col]] < 0.1 |
                  !is.na(sex_dt[[padj_f_col]]) & sex_dt[[padj_f_col]] < 0.1

    permissive <- sex_dt[sig_either & !is.na(lfc_diff_abs) & lfc_diff_abs > 1.0]
    cat("  Permissive sex-divergent (|LFC diff| > 1.0 + sig in ≥1 sex):", nrow(permissive), "\n")

    # Also try threshold of 0.5
    permissive_05 <- sex_dt[sig_either & !is.na(lfc_diff_abs) & lfc_diff_abs > 0.5]
    cat("  Permissive (|LFC diff| > 0.5):", nrow(permissive_05), "\n")

    # Save permissive results
    if (nrow(permissive) > 0) {
      permissive_out <- permissive[order(-lfc_diff_abs)]
      fwrite(permissive_out, file.path(OUTDIR, "sex_divergent_permissive.csv"))
      cat("  Saved: sex_divergent_permissive.csv\n")
    }

    # Summary comparison
    summary_dt <- data.table(
      definition = c("Stringent (both sig + opposite)",
                      "Permissive (|LFC diff| > 1.0 + sig in 1)",
                      "Permissive (|LFC diff| > 0.5 + sig in 1)"),
      n_genes = c(nrow(stringent), nrow(permissive), nrow(permissive_05))
    )
    fwrite(summary_dt, file.path(OUTDIR, "sex_divergent_definitions_comparison.csv"))
    cat("  Saved: sex_divergent_definitions_comparison.csv\n")
    print(summary_dt)
  } else {
    cat("  Could not find LFC columns:", lfc_m_col, lfc_f_col, "\n")
  }
} else {
  cat("  Sex classification file not found\n")
}

# ================================================================
# Part C: Liver-Specific CGP Filter (Issue 7)
# ================================================================
cat("\n--- Part C: Liver-Specific CGP Filter ---\n")

cgp_file <- file.path(BASE, "RNA-seq/results/drug_repurposing/cgp_reversal_hits.csv")
if (file.exists(cgp_file)) {
  cgp <- fread(cgp_file)
  cat("  Total CGP reversal hits:", nrow(cgp), "\n")

  # Liver/metabolic keywords
  liver_keywords <- "LIVER|HEPAT|STEATO|HCC|NAFLD|NASH|CIRR"
  metabolic_keywords <- "LIPID|CHOLEST|METABOL|INSULIN|ADIPOSE|FAT|OBESE|OBESITY|BMI"
  fibrosis_keywords <- "FIBROSIS|STELLATE|COLLAGEN|TGF"
  all_keywords <- paste(liver_keywords, metabolic_keywords, fibrosis_keywords, sep = "|")

  cgp[, is_liver_metabolic := grepl(all_keywords, toupper(pathway))]
  liver_hits <- cgp[is_liver_metabolic == TRUE]
  cat("  Liver/metabolic-related CGP hits:", nrow(liver_hits), "\n")

  if (nrow(liver_hits) > 0) {
    fwrite(liver_hits, file.path(OUTDIR, "cgp_liver_metabolic_filtered.csv"))
    cat("  Saved: cgp_liver_metabolic_filtered.csv\n")
    cat("  Top liver/metabolic hits:\n")
    print(head(liver_hits[, .(pathway, NES, padj)], 20))
  }

  # Cancer overlap analysis: which cancer CGP sets might relate to MASLD→HCC?
  hcc_related <- cgp[grepl("LIVER|HEPAT|HCC|HEPATOCELLULAR|HEPATOBLAST", toupper(pathway))]
  cat("\n  HCC/hepatocellular-related hits:", nrow(hcc_related), "\n")
  if (nrow(hcc_related) > 0) {
    print(hcc_related[, .(pathway, NES, padj)])
  }

  # Summary by class
  cgp[, cgp_class_refined := fcase(
    grepl(liver_keywords, toupper(pathway)), "Liver/MASLD",
    grepl(metabolic_keywords, toupper(pathway)), "Metabolic",
    grepl(fibrosis_keywords, toupper(pathway)), "Fibrosis",
    grepl("CANCER|CARCINOMA|TUMOR|GLIOMA|LEUKEMIA|LYMPHOMA|MELANOMA|SARCOMA|RCC",
          toupper(pathway)), "Cancer_other",
    grepl("DRUG|TREAT|COMPOUND|DOXO|ETOPO|TAMOX|METFORM", toupper(pathway)), "Drug_treatment",
    default = "Other"
  )]
  class_summary <- cgp[, .N, by = cgp_class_refined][order(-N)]
  cat("\n  CGP reversal hits by refined class:\n")
  print(class_summary)
  fwrite(class_summary, file.path(OUTDIR, "cgp_class_distribution.csv"))
} else {
  cat("  CGP reversal file not found\n")
}

# ================================================================
# Part D: Conserved Non-GWAS Enrichment (Issue 12)
# ================================================================
cat("\n--- Part D: Conserved Non-GWAS Enrichment ---\n")

concordance_file <- file.path(BASE, "Analysis/Cross_Species_Concordance/results/concordance_atlas_unified.csv")
if (file.exists(concordance_file)) {
  concordance <- fread(concordance_file)
  cat("  Total genes in concordance atlas:", nrow(concordance), "\n")

  conserved <- concordance[primary_category == "Conserved"]$human_symbol
  all_genes <- concordance$human_symbol
  cat("  Conserved genes:", length(conserved), "\n")

  # Test 1: Enrichment for pathway membership (using GSEA leading edges)
  gsea_file <- file.path(RDIR, "gsea_results.csv")
  if (file.exists(gsea_file)) {
    gsea <- fread(gsea_file)
    sig_pathways <- gsea[padj < 0.05]

    if (nrow(sig_pathways) > 0 && "leadingEdge" %in% names(sig_pathways)) {
      # Get all leading edge genes
      le_genes <- unique(unlist(strsplit(sig_pathways$leadingEdge, ";")))
      # Map ensembl to symbols
      ann_cache <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/gene_annotation/human_ensg_to_symbol.tsv")
      if (file.exists(ann_cache)) {
        ann <- fread(ann_cache)
        le_symbols <- ann[gene_base %in% sub("\\..*", "", le_genes) & symbol != "", symbol]
      } else {
        le_symbols <- le_genes
      }

      # Fisher's exact: Conserved vs leading edge membership
      in_cc_in_le <- sum(conserved %in% le_symbols)
      in_cc_not_le <- length(conserved) - in_cc_in_le
      not_cc_in_le <- sum(all_genes %in% le_symbols) - in_cc_in_le
      not_cc_not_le <- length(all_genes) - in_cc_in_le - in_cc_not_le - not_cc_in_le

      ft_pathway <- fisher.test(matrix(c(in_cc_in_le, not_cc_in_le,
                                          in_cc_not_le, not_cc_not_le), nrow = 2))
      cat(sprintf("  Pathway enrichment: OR=%.2f, p=%.4f\n",
                  ft_pathway$estimate, ft_pathway$p.value))
    }
  }

  # Test 2: Enrichment for drug targets (DGIdb)
  dgidb_file <- file.path(BASE, "RNA-seq/results/drug_repurposing/dgidb_drug_gene_interactions.csv")
  if (file.exists(dgidb_file)) {
    dgidb <- fread(dgidb_file)
    gene_col <- if ("gene" %in% names(dgidb)) "gene" else if ("symbol" %in% names(dgidb)) "symbol" else NULL
    if (!is.null(gene_col)) {
      drug_targets <- unique(dgidb[[gene_col]])

      in_cc_drug <- sum(conserved %in% drug_targets)
      in_cc_not_drug <- length(conserved) - in_cc_drug
      not_cc_drug <- sum(all_genes %in% drug_targets) - in_cc_drug
      not_cc_not_drug <- length(all_genes) - in_cc_drug - in_cc_not_drug - not_cc_drug

      ft_drug <- fisher.test(matrix(c(in_cc_drug, not_cc_drug,
                                       in_cc_not_drug, not_cc_not_drug), nrow = 2))
      cat(sprintf("  Drug target enrichment: OR=%.2f, p=%.4f\n",
                  ft_drug$estimate, ft_drug$p.value))
    }
  }

  # Test 3: Enrichment for multi-evidence high-scorers
  me_file <- file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_scored_genes.csv")
  if (file.exists(me_file)) {
    me <- fread(me_file)
    top5pct <- me[composite_score >= quantile(composite_score, 0.95)]$human_symbol

    in_cc_top <- sum(conserved %in% top5pct)
    in_cc_not_top <- length(conserved) - in_cc_top
    not_cc_top <- sum(all_genes %in% top5pct) - in_cc_top
    not_cc_not_top <- length(all_genes) - in_cc_top - in_cc_not_top - not_cc_top

    ft_me <- fisher.test(matrix(c(in_cc_top, not_cc_top,
                                    in_cc_not_top, not_cc_not_top), nrow = 2))
    cat(sprintf("  Multi-evidence top 5%% enrichment: OR=%.2f, p=%.4f\n",
                ft_me$estimate, ft_me$p.value))
  }

  # Test 4: Enrichment for Significant DEGs
  consensus_file <- file.path(RDIR, "consensus_degs.csv")
  if (file.exists(consensus_file)) {
    consensus <- fread(consensus_file)
    # Map to symbols
    ortho_file <- file.path(RDIR, "human_mouse_ortholog_comparison.csv")
    if (file.exists(ortho_file)) {
      ortho <- fread(ortho_file)
      ortho[, ensembl_clean := sub("\\..*", "", gene_base)]
      sym_map <- unique(ortho[!is.na(human_symbol) & human_symbol != "",
                               .(ensembl_clean, human_symbol)])
      consensus[, ensembl_clean := sub("\\..*", "", gene)]
      consensus <- merge(consensus, sym_map, by = "ensembl_clean", all.x = TRUE)
    }

    if ("human_symbol" %in% names(consensus)) {
      if (!"dream_padj" %in% names(consensus) && "padj" %in% names(consensus))
        setnames(consensus, "padj", "dream_padj")
      if (!"dream_logFC" %in% names(consensus) && "logFC" %in% names(consensus))
        setnames(consensus, "logFC", "dream_logFC")

      dream_degs <- consensus[!is.na(dream_padj) & dream_padj < 0.1 & abs(dream_logFC) >= 0.58]$human_symbol

      in_cc_dream <- sum(conserved %in% dream_degs)
      in_cc_not_dream <- length(conserved) - in_cc_dream
      not_cc_dream <- sum(all_genes %in% dream_degs) - in_cc_dream
      not_cc_not_dream <- length(all_genes) - in_cc_dream - in_cc_not_dream - not_cc_dream

      ft_dream <- fisher.test(matrix(c(in_cc_dream, not_cc_dream,
                                      in_cc_not_dream, not_cc_not_dream), nrow = 2))
      cat(sprintf("  Significant DEGs enrichment: OR=%.2f, p=%.4f\n",
                  ft_dream$estimate, ft_dream$p.value))
    }
  }

  # Save summary
  enrichment_summary <- data.table(
    test = character(), OR = numeric(), p_value = numeric()
  )
  if (exists("ft_pathway"))
    enrichment_summary <- rbind(enrichment_summary,
      data.table(test = "GSEA_leading_edge", OR = round(ft_pathway$estimate, 3),
                 p_value = ft_pathway$p.value))
  if (exists("ft_drug"))
    enrichment_summary <- rbind(enrichment_summary,
      data.table(test = "DGIdb_drug_targets", OR = round(ft_drug$estimate, 3),
                 p_value = ft_drug$p.value))
  if (exists("ft_me"))
    enrichment_summary <- rbind(enrichment_summary,
      data.table(test = "Multi_evidence_top5pct", OR = round(ft_me$estimate, 3),
                 p_value = ft_me$p.value))
  if (exists("ft_dream"))
    enrichment_summary <- rbind(enrichment_summary,
      data.table(test = "Dream_DEG", OR = round(ft_dream$estimate, 3),
                 p_value = ft_dream$p.value))

  if (nrow(enrichment_summary) > 0) {
    fwrite(enrichment_summary, file.path(OUTDIR, "conserved_enrichment.csv"))
    cat("\n  Enrichment summary:\n")
    print(enrichment_summary)
    cat("  Saved: conserved_enrichment.csv\n")
  }
} else {
  cat("  Concordance atlas not found\n")
}

cat("\n=== Audit Sensitivity Analyses Complete ===\n")
cat("End time:", format(Sys.time()), "\n")
