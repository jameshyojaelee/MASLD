#!/usr/bin/env Rscript
# =============================================================================
# 41_ferroptosis_lipotoxicity.R
# Module B2: Ferroptosis/Lipotoxicity Program Annotation
#
# Annotates DEGs for ferroptosis pathway membership and quantifies ferroptosis
# as a transcriptionally activated MASLD program.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr, warn.conflicts = FALSE)
  library(tidyr)
  library(msigdbr)
  library(fgsea)
})

# Ensure dplyr verbs take priority over data.table
select <- dplyr::select
filter <- dplyr::filter

# --- Configuration ---
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

atlas_path   <- file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
dream_path   <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/canonical_deg_results.csv")
hep_reg_path <- file.path(BASE, "Analysis/ATAC/Human_Multiome/scenic_plus/hepatocyte_regulons.csv")
dis_reg_path <- file.path(BASE, "Analysis/ATAC/Human_Multiome/scenic_plus/disease_regulons.csv")
out_dir      <- file.path(BASE, "RNA-seq/results/pathway_programs")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== Module B2: Ferroptosis/Lipotoxicity Program ===\n")
cat("Start:", format(Sys.time()), "\n\n")

# =============================================================================
# 1. Build Ferroptosis Gene Sets
# =============================================================================
cat("--- Building ferroptosis gene sets ---\n")

# Curated FerrDb v2 gene lists (literature-derived, as direct download is unreliable)
# Sources: Dixon 2012, Stockwell 2017, Jiang 2021, Chen 2021, Liang 2024
ferroptosis_drivers <- c(
  "ACSL4", "LPCAT3", "ALOX5", "ALOX12", "ALOX15", "ALOX15B",
  "TFRC", "SLC11A2", "NCOA4", "IREB2", "HMOX1", "HMOX2",
  "NOX1", "NOX4", "CYBB", "POR", "CBS", "DPP4",
  "VDAC2", "VDAC3", "CARS1", "GLS2", "FDFT1",
  "CHAC1", "ATG5", "ATG7", "BECN1", "SAT1", "RPL8",
  "HSPB1", "FANCD2", "CISD1", "CS", "ACSF2", "ACSL3"
)

ferroptosis_suppressors <- c(
  "GPX4", "SLC7A11", "SLC3A2", "FSP1", "DHODH", "GCH1",
  "NFE2L2", "KEAP1", "NQO1", "FTH1", "FTL", "SLC40A1",
  "GCLC", "GCLM", "GSS", "GSR", "GPX1", "GPX2", "GPX3",
  "TXNRD1", "TXN", "PRDX1", "SOD1", "SOD2", "CAT",
  "HMGCR", "SQLE", "SC5D", "TP53", "CDKN1A", "SLC7A5",
  "HCAR1", "AIFM2", "PROMININ2", "PEBP1"
)

ferroptosis_markers <- c(
  "PTGS2", "CHAC1", "ACSL4", "TFRC", "SLC7A11", "GPX4",
  "NOX1", "HMOX1", "FTH1", "FTL", "SAT1", "GDF15",
  "LAMP2", "NCOA4", "SLC11A2", "LPCAT3", "GCH1"
)

# Lipotoxicity / ER stress genes (curated from MASLD literature)
lipotoxicity_genes <- c(
  # ER stress / UPR
  "DDIT3", "ATF4", "ATF6", "ERN1", "EIF2AK3", "XBP1", "HSPA5",
  "EDEM1", "HERPUD1", "DNAJB9", "DERL1", "CALR", "CANX",
  # Lipotoxic signaling
  "CERS2", "CERS4", "CERS6", "SMPD1", "SMPD3", "ASAH1",
  "DGAT1", "DGAT2", "AGPAT2", "GPAM", "LPIN1",
  # Oxidative stress
  "CYP2E1", "CYP4A11", "SOD2", "CAT", "GPX1",
  # Inflammasome/pyroptosis cross-talk
  "NLRP3", "GSDMD", "CASP1", "IL1B", "IL18",
  # Ceramide/sphingolipid
  "SPTLC1", "SPTLC2", "CERK", "S1PR1", "SGMS1", "UGCG",
  # Mitochondrial dysfunction
  "BNIP3", "PINK1", "PARK2", "MFN1", "MFN2", "DRP1",
  "PPARGC1A", "PPARGC1B", "SIRT1", "SIRT3"
)

# Get ferroptosis pathways from MSigDB (v10+ API: use collection/subcollection)
cat("Fetching MSigDB ferroptosis gene sets...\n")

# WikiPathways Ferroptosis (replaces KEGG which has no ferroptosis set)
wp_sets <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:WIKIPATHWAYS")
wp_ferroptosis <- wp_sets %>%
  filter(gs_name == "WP_FERROPTOSIS") %>%
  pull(gene_symbol) %>% unique()
cat("  WikiPathways Ferroptosis: ", length(wp_ferroptosis), " genes\n")

# GO ferroptosis
go_sets <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP")
go_ferroptosis <- go_sets %>%
  filter(grepl("FERROPTOSIS", gs_name, ignore.case = TRUE)) %>%
  pull(gene_symbol) %>% unique()
cat("  GO Ferroptosis: ", length(go_ferroptosis), " genes\n")

# Reactome ferroptosis (if available)
reactome_sets <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME")
reactome_ferroptosis <- reactome_sets %>%
  filter(grepl("FERROPTOSIS", gs_name, ignore.case = TRUE)) %>%
  pull(gene_symbol) %>% unique()
cat("  Reactome Ferroptosis: ", length(reactome_ferroptosis), " genes\n")

# Hallmark gene sets for pathway context
hallmark_sets <- msigdbr(species = "Homo sapiens", collection = "H")
hallmark_list <- split(hallmark_sets$gene_symbol, hallmark_sets$gs_name)

# Combined ferroptosis universe (union of all sources)
ferroptosis_all <- unique(c(ferroptosis_drivers, ferroptosis_suppressors,
                            wp_ferroptosis, go_ferroptosis, reactome_ferroptosis))
cat("  Combined ferroptosis genes: ", length(ferroptosis_all), "\n")

# =============================================================================
# 2. Load Atlas and Dream Data
# =============================================================================
cat("\n--- Loading data ---\n")
atlas <- fread(atlas_path)
dream <- fread(dream_path)

# Dream results: annotate with gene symbol
# gene column has ENSG IDs - need to map via atlas
dream_with_symbol <- dream %>%
  mutate(ensembl_id = sub("\\..*", "", gene)) %>%
  left_join(atlas %>% select(human_symbol, ensembl_id), by = "ensembl_id")

cat("  Atlas: ", nrow(atlas), " genes\n")
cat("  Dream results: ", nrow(dream), " genes\n")

# DEG definitions
atlas <- atlas %>%
  mutate(
    is_deg = !is.na(bulk_padj) & bulk_padj < 0.1,  # Exploratory annotation threshold; primary DEGs: padj<0.05 + |logFC|>0.5 (Script 05b)
    is_deg_up = is_deg & !is.na(bulk_logFC) & bulk_logFC > 0,
    is_deg_down = is_deg & !is.na(bulk_logFC) & bulk_logFC < 0
  )

n_deg <- sum(atlas$is_deg, na.rm = TRUE)
n_deg_up <- sum(atlas$is_deg_up, na.rm = TRUE)
n_deg_down <- sum(atlas$is_deg_down, na.rm = TRUE)
cat(sprintf("  DEGs: %d total (%d up, %d down)\n", n_deg, n_deg_up, n_deg_down))

# =============================================================================
# 3. Annotate DEGs with Ferroptosis Classification
# =============================================================================
cat("\n--- Annotating DEGs with ferroptosis class ---\n")

atlas <- atlas %>%
  mutate(
    ferroptosis_driver = human_symbol %in% ferroptosis_drivers,
    ferroptosis_suppressor = human_symbol %in% ferroptosis_suppressors,
    ferroptosis_marker = human_symbol %in% ferroptosis_markers,
    ferroptosis_wp = human_symbol %in% wp_ferroptosis,
    ferroptosis_go = human_symbol %in% go_ferroptosis,
    ferroptosis_any = human_symbol %in% ferroptosis_all,
    lipotoxicity = human_symbol %in% lipotoxicity_genes,
    ferroptosis_class = case_when(
      ferroptosis_driver & ferroptosis_suppressor ~ "Dual",
      ferroptosis_driver ~ "Driver",
      ferroptosis_suppressor ~ "Suppressor",
      ferroptosis_wp | ferroptosis_go ~ "Pathway_member",
      ferroptosis_marker & !ferroptosis_driver & !ferroptosis_suppressor ~ "Marker",
      TRUE ~ "None"
    )
  )

# Summary
cat("\nFerroptosis annotation summary:\n")
cat(sprintf("  Drivers in atlas: %d\n", sum(atlas$ferroptosis_driver)))
cat(sprintf("  Suppressors in atlas: %d\n", sum(atlas$ferroptosis_suppressor)))
cat(sprintf("  WikiPathways in atlas: %d\n", sum(atlas$ferroptosis_wp)))
cat(sprintf("  GO ferroptosis in atlas: %d\n", sum(atlas$ferroptosis_go)))
cat(sprintf("  Any ferroptosis in atlas: %d\n", sum(atlas$ferroptosis_any)))
cat(sprintf("  Lipotoxicity in atlas: %d\n", sum(atlas$lipotoxicity)))

# =============================================================================
# 4. Enrichment Testing
# =============================================================================
cat("\n--- Enrichment testing ---\n")

run_fisher <- function(gene_set_name, gene_set, deg_set, universe) {
  # 2x2: gene_set x DEG status
  in_set_deg <- sum(universe %in% gene_set & universe %in% deg_set)
  in_set_nodeg <- sum(universe %in% gene_set & !(universe %in% deg_set))
  noset_deg <- sum(!(universe %in% gene_set) & universe %in% deg_set)
  noset_nodeg <- sum(!(universe %in% gene_set) & !(universe %in% deg_set))

  mat <- matrix(c(in_set_deg, in_set_nodeg, noset_deg, noset_nodeg), nrow = 2)
  ft <- fisher.test(mat, alternative = "greater")

  data.frame(
    gene_set = gene_set_name,
    n_in_set = sum(universe %in% gene_set),
    n_deg_in_set = in_set_deg,
    n_deg_total = length(deg_set),
    n_universe = length(universe),
    odds_ratio = ft$estimate,
    pvalue = ft$p.value,
    ci_lower = ft$conf.int[1],
    ci_upper = ft$conf.int[2],
    stringsAsFactors = FALSE
  )
}

universe <- atlas$human_symbol
deg_all <- atlas$human_symbol[atlas$is_deg]
deg_up <- atlas$human_symbol[atlas$is_deg_up]
deg_down <- atlas$human_symbol[atlas$is_deg_down]

# Hepatocyte-intrinsic DEGs
hep_intrinsic <- atlas$human_symbol[!is.na(atlas$attribution_class) &
  atlas$attribution_class == "Hepatocyte_Intrinsic"]
cat(sprintf("  Hepatocyte-intrinsic DEGs: %d\n", length(hep_intrinsic)))

# Conserved
conserved <- atlas$human_symbol[!is.na(atlas$is_conserved) &
  atlas$is_conserved == TRUE]
cat(sprintf("  Conserved: %d\n", length(conserved)))

# Run Fisher's tests
enrichment_results <- bind_rows(
  # Ferroptosis drivers
  run_fisher("Ferroptosis_Drivers_vs_AllDEG", ferroptosis_drivers, deg_all, universe),
  run_fisher("Ferroptosis_Drivers_vs_UpDEG", ferroptosis_drivers, deg_up, universe),
  run_fisher("Ferroptosis_Drivers_vs_DownDEG", ferroptosis_drivers, deg_down, universe),
  run_fisher("Ferroptosis_Drivers_vs_HepIntrinsic", ferroptosis_drivers, hep_intrinsic, universe),

  # Ferroptosis suppressors
  run_fisher("Ferroptosis_Suppressors_vs_AllDEG", ferroptosis_suppressors, deg_all, universe),
  run_fisher("Ferroptosis_Suppressors_vs_UpDEG", ferroptosis_suppressors, deg_up, universe),
  run_fisher("Ferroptosis_Suppressors_vs_DownDEG", ferroptosis_suppressors, deg_down, universe),

  # Combined ferroptosis
  run_fisher("Ferroptosis_Any_vs_AllDEG", ferroptosis_all, deg_all, universe),
  run_fisher("Ferroptosis_Any_vs_UpDEG", ferroptosis_all, deg_up, universe),
  run_fisher("Ferroptosis_Any_vs_DownDEG", ferroptosis_all, deg_down, universe),
  run_fisher("Ferroptosis_Any_vs_HepIntrinsic", ferroptosis_all, hep_intrinsic, universe),

  # WikiPathways ferroptosis
  run_fisher("WP_Ferroptosis_vs_AllDEG", wp_ferroptosis, deg_all, universe),
  run_fisher("WP_Ferroptosis_vs_UpDEG", wp_ferroptosis, deg_up, universe),

  # Lipotoxicity
  run_fisher("Lipotoxicity_vs_AllDEG", lipotoxicity_genes, deg_all, universe),
  run_fisher("Lipotoxicity_vs_UpDEG", lipotoxicity_genes, deg_up, universe),
  run_fisher("Lipotoxicity_vs_HepIntrinsic", lipotoxicity_genes, hep_intrinsic, universe),

  # Ferroptosis in Conserved
  run_fisher("Ferroptosis_Any_vs_ConservedCore", ferroptosis_all, conserved, universe),
  run_fisher("Ferroptosis_Drivers_vs_ConservedCore", ferroptosis_drivers, conserved, universe)
)

# BH correction
enrichment_results$padj <- p.adjust(enrichment_results$pvalue, method = "BH")

cat("\nEnrichment results:\n")
print(enrichment_results %>%
  select(gene_set, n_in_set, n_deg_in_set, odds_ratio, pvalue, padj) %>%
  arrange(pvalue), row.names = FALSE)

# =============================================================================
# 5. GSEA with Dream t-statistics
# =============================================================================
cat("\n--- GSEA analysis ---\n")

# Build ranked gene list from dream t-statistics
ranked_genes <- dream_with_symbol %>%
  filter(!is.na(human_symbol)) %>%
  group_by(human_symbol) %>%
  slice_max(abs(t), n = 1) %>%
  ungroup() %>%
  arrange(desc(t))

gene_ranks <- setNames(ranked_genes$t, ranked_genes$human_symbol)
gene_ranks <- gene_ranks[!is.na(gene_ranks)]

# Create gene set list for fgsea
ferroptosis_gsets <- list(
  Ferroptosis_Drivers = intersect(ferroptosis_drivers, names(gene_ranks)),
  Ferroptosis_Suppressors = intersect(ferroptosis_suppressors, names(gene_ranks)),
  Ferroptosis_WP = intersect(wp_ferroptosis, names(gene_ranks)),
  Ferroptosis_GO = intersect(go_ferroptosis, names(gene_ranks)),
  Ferroptosis_All = intersect(ferroptosis_all, names(gene_ranks)),
  Lipotoxicity = intersect(lipotoxicity_genes, names(gene_ranks))
)

# Filter out empty sets
ferroptosis_gsets <- ferroptosis_gsets[sapply(ferroptosis_gsets, length) >= 5]

if (length(ferroptosis_gsets) > 0) {
  gsea_results <- fgsea(
    pathways = ferroptosis_gsets,
    stats = gene_ranks,
    minSize = 5,
    maxSize = 500,
    nPermSimple = 10000
  )

  gsea_df <- as.data.frame(gsea_results)
  gsea_df$leadingEdge <- sapply(gsea_df$leadingEdge, function(x) paste(x, collapse = ";"))

  cat("\nGSEA results:\n")
  print(gsea_df %>%
    select(pathway, pval, padj, NES, size) %>%
    arrange(pval), row.names = FALSE)
} else {
  gsea_df <- data.frame()
  cat("  No ferroptosis gene sets with >= 5 genes in ranked list\n")
}

# =============================================================================
# 6. Disease Regulon × Ferroptosis Cross-Reference
# =============================================================================
cat("\n--- Disease regulon × ferroptosis overlap ---\n")

regulons <- fread(dis_reg_path)
hep_regs <- fread(hep_reg_path)

# Get disease regulon TFs
disease_tfs <- regulons$tf_name

# Get target genes per disease regulon
regulon_targets <- hep_regs %>%
  filter(tf_name %in% disease_tfs) %>%
  select(tf_name, target_gene)

# Check overlap with ferroptosis genes
regulon_ferroptosis <- regulon_targets %>%
  mutate(
    is_ferroptosis = target_gene %in% ferroptosis_all,
    ferroptosis_class = case_when(
      target_gene %in% ferroptosis_drivers ~ "Driver",
      target_gene %in% ferroptosis_suppressors ~ "Suppressor",
      target_gene %in% ferroptosis_all ~ "Pathway_member",
      TRUE ~ "None"
    )
  )

n_overlap <- sum(regulon_ferroptosis$is_ferroptosis)
cat(sprintf("  Regulon targets that are ferroptosis genes: %d / %d\n",
            n_overlap, nrow(regulon_ferroptosis)))

if (n_overlap > 0) {
  cat("  Ferroptosis genes in disease regulons:\n")
  print(regulon_ferroptosis %>%
    filter(is_ferroptosis) %>%
    select(tf_name, target_gene, ferroptosis_class), row.names = FALSE)
}

# Also check if any disease TFs ARE ferroptosis genes
tf_ferroptosis <- data.frame(
  tf_name = disease_tfs,
  is_ferroptosis_gene = disease_tfs %in% ferroptosis_all,
  ferroptosis_class = case_when(
    disease_tfs %in% ferroptosis_drivers ~ "Driver",
    disease_tfs %in% ferroptosis_suppressors ~ "Suppressor",
    disease_tfs %in% ferroptosis_all ~ "Pathway_member",
    TRUE ~ "None"
  )
)
cat("\n  Disease TFs that are ferroptosis genes:\n")
print(tf_ferroptosis %>% filter(is_ferroptosis_gene), row.names = FALSE)

# =============================================================================
# 7. Sex Stratification of Ferroptosis Program
# =============================================================================
cat("\n--- Sex stratification ---\n")

ferroptosis_degs <- atlas %>%
  filter(ferroptosis_any & is_deg) %>%
  select(human_symbol, ferroptosis_class, bulk_logFC, bulk_padj,
         sex_class, bulk_logFC_M, bulk_logFC_F, sex_interaction_padj,
         is_conserved, attribution_class)

sex_summary <- ferroptosis_degs %>%
  group_by(sex_class) %>%
  summarise(n = n(), .groups = "drop")

cat("  Ferroptosis DEGs by sex class:\n")
print(sex_summary, row.names = FALSE)

# =============================================================================
# 8. Output: Ferroptosis DEG Table
# =============================================================================
cat("\n--- Writing outputs ---\n")

# Ferroptosis-annotated DEGs
ferroptosis_deg_table <- atlas %>%
  filter(ferroptosis_any | lipotoxicity) %>%
  select(human_symbol, ensembl_id, ferroptosis_class,
         ferroptosis_driver, ferroptosis_suppressor, ferroptosis_wp,
         ferroptosis_go, lipotoxicity,
         bulk_logFC, bulk_padj, bulk_tstat, is_deg, is_deg_up, is_deg_down,
         mouse_meta_logFC, mouse_meta_padj, is_conserved,
         sex_class, attribution_class) %>%
  arrange(bulk_padj)

fwrite(ferroptosis_deg_table,
       file.path(out_dir, "ferroptosis_degs.csv"))
cat("  Written: ferroptosis_degs.csv (", nrow(ferroptosis_deg_table), " genes)\n")

# Lipotoxicity program
lipotox_table <- atlas %>%
  filter(lipotoxicity) %>%
  select(human_symbol, ensembl_id,
         bulk_logFC, bulk_padj, bulk_tstat, is_deg,
         mouse_meta_logFC, is_conserved,
         sex_class, attribution_class) %>%
  arrange(bulk_padj)

fwrite(lipotox_table,
       file.path(out_dir, "lipotoxicity_program.csv"))
cat("  Written: lipotoxicity_program.csv (", nrow(lipotox_table), " genes)\n")

# Enrichment statistics
fwrite(enrichment_results,
       file.path(out_dir, "ferroptosis_enrichment.csv"))
cat("  Written: ferroptosis_enrichment.csv\n")

# GSEA results
if (nrow(gsea_df) > 0) {
  fwrite(gsea_df,
         file.path(out_dir, "ferroptosis_gsea.csv"))
  cat("  Written: ferroptosis_gsea.csv\n")
}

# Regulon × ferroptosis overlap
fwrite(regulon_ferroptosis,
       file.path(out_dir, "regulon_ferroptosis_overlap.csv"))
cat("  Written: regulon_ferroptosis_overlap.csv\n")

# =============================================================================
# 9. Summary Statistics
# =============================================================================
cat("\n--- Final Summary ---\n")

n_ferr_deg <- sum(atlas$ferroptosis_any & atlas$is_deg, na.rm = TRUE)
n_ferr_deg_up <- sum(atlas$ferroptosis_any & atlas$is_deg_up, na.rm = TRUE)
n_ferr_deg_down <- sum(atlas$ferroptosis_any & atlas$is_deg_down, na.rm = TRUE)
n_ferr_hep <- sum(atlas$ferroptosis_any & atlas$human_symbol %in% hep_intrinsic, na.rm = TRUE)
n_ferr_core <- sum(atlas$ferroptosis_any & atlas$human_symbol %in% conserved, na.rm = TRUE)

cat(sprintf("Ferroptosis genes that are DEGs: %d (%d up, %d down)\n",
            n_ferr_deg, n_ferr_deg_up, n_ferr_deg_down))
cat(sprintf("Ferroptosis genes in hepatocyte-intrinsic: %d\n", n_ferr_hep))
cat(sprintf("Ferroptosis genes in Conserved: %d\n", n_ferr_core))

# Check canonical genes
canonical <- c("GPX4", "SLC7A11", "ACSL4", "TFRC", "HMOX1", "NFE2L2", "LPCAT3")
cat("\nCanonical ferroptosis genes:\n")
for (g in canonical) {
  ar <- atlas[atlas$human_symbol == g, ]
  if (nrow(ar) > 0) {
    ar <- ar[1, ]
    deg_str <- if (ar$is_deg) sprintf("DEG (LFC=%.2f, padj=%.2e)", ar$bulk_logFC, ar$bulk_padj) else "Not DEG"
    cat(sprintf("  %s: %s, class=%s\n", g, deg_str,
                ar$ferroptosis_class))
  } else {
    cat(sprintf("  %s: Not in atlas\n", g))
  }
}

cat("\n=== Module B2 Complete ===\n")
cat("End:", format(Sys.time()), "\n")
