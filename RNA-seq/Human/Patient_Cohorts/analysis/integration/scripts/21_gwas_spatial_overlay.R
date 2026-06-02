#!/usr/bin/env Rscript
# 21_gwas_spatial_overlay.R
# ---------------------------------------------------------------------------
# Convergent Evidence: GWAS Risk Loci & Spatial Transcriptomics
# Objective: Test if our consensus DEGs are enriched for human genetic risk
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

cat("=== Phase 9: Genetic Convergence (GWAS Overlay) ===\n\n")

WD  <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq"
RES <- file.path(WD, "Human/Patient_Cohorts/analysis/integration/results")
OUT <- file.path(RES, "gwas_spatial_convergence")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 1. Load Unified Concordance & GWAS Genes
# ============================================================
cat("Loading cross-species concordance atlas...\n")
degs <- fread(file.path(WD, "Analysis/Cross_Species_Concordance/results/concordance_atlas_unified.csv"))

cat("Loading MASLD GWAS target genes (curated)...\n")
# The GWAS file contains raw parsed symbols that might have pipe delimiters (e.g. CLN3|APOBR)
gwas_raw <- fread(file.path(WD, "../GWAS/Closest_genes.csv"), header=FALSE)$V1
gwas_list <- unique(unlist(strsplit(gwas_raw, "\\|")))

cat(sprintf("Parsed %d unique GWAS genes from input.\n", length(gwas_list)))

# ============================================================
# 2. Enrichment Analysis (GWAS)
# ============================================================
cat("\nComputing Fisher Exact Tests for GWAS enrichment...\n")
# Background = all genes in the concordance atlas (the tested gene universe)
n_total <- nrow(degs)
n_gwas <- length(gwas_list)
cat(sprintf("Background universe: %d genes in concordance atlas\n", n_total))
cat(sprintf("GWAS gene list: %d genes\n", n_gwas))

# We want to test enrichment inside the main consensus tiers (or translatability tier)
enrichment_results <- list()

# Function to run Fisher exact test
run_fisher <- function(target_set, set_name, category_name) {
  overlap <- intersect(target_set, gwas_list)
  
  contingency <- matrix(c(
    length(overlap), length(target_set) - length(overlap),
    n_gwas - length(overlap), n_total - n_gwas - length(target_set) + length(overlap)
  ), nrow=2)
  
  test <- fisher.test(contingency, alternative="greater")
  
  data.table(
    Category = category_name,
    SubCategory = set_name,
    Set_Size = length(target_set),
    Overlap = length(overlap),
    Expected = (length(target_set) * n_gwas) / n_total,
    Odds_Ratio = test$estimate,
    P_Value = test$p.value
  )
}

# Test 1: By Translatability Tier
for (t_tier in unique(degs[translatability_tier != "", translatability_tier])) {
  genes <- degs[translatability_tier == t_tier, human_symbol]
  enrichment_results[[paste0("trans_", t_tier)]] <- run_fisher(genes, t_tier, "Translatability Tier")
}

# Test 2: By Primary Category (Conserved vs Discordant, etc)
for (p_cat in unique(degs[primary_category != "", primary_category])) {
  genes <- degs[primary_category == p_cat, human_symbol]
  enrichment_results[[paste0("cat_", p_cat)]] <- run_fisher(genes, p_cat, "Primary Category")
}

# Combine results
enrich_dt <- rbindlist(enrichment_results)
enrich_dt[, FDR := p.adjust(P_Value, method="fdr")]

print(enrich_dt[order(P_Value)])

fwrite(enrich_dt, file=file.path(OUT, "gwas_enrichment_results.csv"))

# ============================================================
# 3. Create Evidence Output Table (The Library Intersect)
# ============================================================
cat("\nGenerating GWAS evidence overlap table...\n")
overlap_dt <- degs[human_symbol %in% gwas_list]

overlapColumns <- c("human_symbol", "translatability_tier", "primary_category", "best_signature")
if("human_consensus_tier" %in% names(degs)) {
    overlapColumns <- c(overlapColumns, "human_consensus_tier")
}

evidence_table <- overlap_dt[, ..overlapColumns, nomatch=NULL]
fwrite(evidence_table, file.path(OUT, "gwas_overlapped_consensus_genes.csv"))

cat(sprintf("Found %d high-confidence consensus DEGs that overlap with GWAS loci.\n", nrow(evidence_table)))

cat("\nGWAS spatial overlay completed successfully.\n")
