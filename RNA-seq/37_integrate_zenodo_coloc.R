#!/usr/bin/env Rscript
# 37_integrate_zenodo_coloc.R
# ---------------------------------------------------------------------------
# Integrate Precomputed Zenodo COLOC Results with MASLD DEG Evidence
#
# Parses cell-type-specific COLOC results from the MASLD sc-eQTL study
# (Nature Genetics 2025) and cross-references with dream DEGs, consensus
# tiers, and multi-evidence atlas. Also extracts FOXO1/siEFHD1 perturbation
# data as functional validation evidence.
#
# Inputs:
#   - Zenodo COLOC files (4 cell types x 25 GWAS traits)
#   - Dream results + consensus DEGs (tier1_high_confidence_degs.csv)
#   - Multi-evidence atlas
#   - Perturbation DEG lists (FOXO1 HepG2 1000nM, FOXO1 organoid, siEFHD1)
#
# Outputs:
#   - results/causal_inference/sceqtl/zenodo_coloc_integrated.csv
#   - results/causal_inference/sceqtl/zenodo_perturbation_validation.csv
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

cat("=== Script 37: Zenodo sc-eQTL COLOC Integration ===\n")
cat("Start time:", format(Sys.time()), "\n\n")

# ==============================================================================
# Configuration
# ==============================================================================
BASE_DIR    <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
ZENODO_DIR  <- file.path(BASE_DIR, "RNA-seq/References/eQTL/MASLD_sc_eQTL_2025/Zenodo")
COLOC_DIR   <- file.path(ZENODO_DIR, "coloc")
PERTURB_DIR <- file.path(ZENODO_DIR, "bulkRNA_seq")
RESULTS_DIR <- file.path(BASE_DIR, "RNA-seq/results/causal_inference/sceqtl")
DREAM_FILE  <- file.path(BASE_DIR,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/dream_results.csv")
GENE_CACHE  <- file.path(BASE_DIR,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/gene_annotation/human_ensg_to_symbol.tsv")
TIER1_FILE  <- file.path(BASE_DIR,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/tier1_high_confidence_degs.csv")
ATLAS_FILE  <- file.path(BASE_DIR,
  "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

# Cell types with precomputed COLOC results
CELL_TYPES <- c("hepatocyte", "endothelial_cell", "cholangiocyte", "stellate_cell")

# NAFLD trait column in COLOC files
NAFLD_COL <- "Nat22_cALT_NAFLD_GWAS"

# Dream significance threshold
PADJ_THR <- 0.1

# Known genes of interest for highlighting
HIGHLIGHT_GENES <- c("HSD17B13", "EFHD1", "ETS2", "PNPLA3", "TM6SF2",
                      "GCKR", "MBOAT7", "FOXO1", "COL1A1", "ACTA2")

# Perturbation experiments to load
PERTURB_FILES <- list(
  foxo1_hepg2_1000nM = list(
    file = file.path(PERTURB_DIR, "hepg2_FOXO1_inhibitor_1000nM_vs_0_DEG_list_Padj_0.05.txt.gz"),
    gene_col = "gene_name",
    label = "FOXO1_inhibitor_HepG2_1000nM"
  ),
  foxo1_organoid = list(
    file = file.path(PERTURB_DIR, "organoid_FOXO1_inhibitor_T_vs_F_DEG_list_Padj_0.05.txt.gz"),
    gene_col = "gene_name",
    label = "FOXO1_inhibitor_organoid"
  ),
  siEFHD1 = list(
    file = file.path(PERTURB_DIR, "organoid_siEFHD1_T_vs_F_DEG_list_Padj_0.1.txt.gz"),
    gene_col = "gene",
    label = "siEFHD1_organoid"
  )
)

# ==============================================================================
# Part 1: Load Gene Annotation (Ensembl -> Symbol)
# ==============================================================================
cat("--- Loading gene annotation ---\n")
gene_map <- fread(GENE_CACHE, header = TRUE, sep = "\t")
# gene_map columns: gene_id, symbol, gene_type, gene_base
# Create lookup: versioned Ensembl ID -> symbol
ensg_to_sym <- setNames(gene_map$symbol, gene_map$gene_id)
cat(sprintf("  Loaded %d gene annotations\n", length(ensg_to_sym)))

# ==============================================================================
# Part 2: Load Dream DEGs
# ==============================================================================
cat("--- Loading dream DEGs ---\n")
dream <- fread(DREAM_FILE, header = TRUE)
cat(sprintf("  Dream results: %d genes\n", nrow(dream)))

# Add symbol column by mapping versioned Ensembl IDs
dream[, symbol := ensg_to_sym[gene]]
dream_with_sym <- dream[!is.na(symbol)]
cat(sprintf("  Mapped to symbol: %d genes\n", nrow(dream_with_sym)))

# Create quick lookup: symbol -> dream row
dream_lookup <- dream_with_sym[, .(dream_logFC = logFC, dream_padj = padj), by = symbol]
# Handle duplicates: keep the one with smallest padj
dream_lookup <- dream_lookup[order(dream_padj)][!duplicated(symbol)]
setkey(dream_lookup, symbol)

# DEG set (padj < 0.1)
deg_symbols <- dream_lookup[dream_padj < PADJ_THR, symbol]
cat(sprintf("  DEGs (padj < %.1f): %d\n", PADJ_THR, length(deg_symbols)))

# ==============================================================================
# Part 3: Load Consensus Tiers
# ==============================================================================
# Primary source: multi-evidence atlas (has all 3 tiers via human_consensus_tier).
# Fallback: tier1_high_confidence_degs.csv (Tier 1 only, Ensembl IDs).
cat("--- Loading consensus tiers ---\n")
tier_lookup <- data.table(symbol = character(), consensus_tier = character())

if (file.exists(ATLAS_FILE)) {
  # Atlas has human_consensus_tier for all genes
  atlas_tiers <- fread(ATLAS_FILE, header = TRUE,
                       select = c("human_symbol", "human_consensus_tier"))
  atlas_tiers <- atlas_tiers[human_consensus_tier != "Not_significant"]
  setnames(atlas_tiers, c("human_symbol", "human_consensus_tier"),
           c("symbol", "consensus_tier"))
  tier_lookup <- atlas_tiers[!duplicated(symbol)]
  cat(sprintf("  Loaded tiers from atlas: %d genes\n", nrow(tier_lookup)))
} else if (file.exists(TIER1_FILE)) {
  tiers <- fread(TIER1_FILE, header = TRUE)
  cat(sprintf("  Loaded %d Tier1 DEGs from fallback file\n", nrow(tiers)))
  tiers[, symbol := ensg_to_sym[gene]]
  tier_lookup <- tiers[!is.na(symbol), .(consensus_tier = tier), by = symbol]
  tier_lookup <- tier_lookup[!duplicated(symbol)]
} else {
  cat("  WARNING: No tier data source found, skipping tier annotation\n")
}

setkey(tier_lookup, symbol)
if (nrow(tier_lookup) > 0) {
  cat("  Tier distribution:\n")
  print(tier_lookup[, .N, by = consensus_tier][order(-N)])
}

# ==============================================================================
# Part 4: Load Multi-Evidence Atlas (layers_active)
# ==============================================================================
cat("\n--- Loading multi-evidence atlas ---\n")
if (file.exists(ATLAS_FILE)) {
  atlas <- fread(ATLAS_FILE, header = TRUE,
                 select = c("human_symbol", "layers_active"))
  cat(sprintf("  Atlas: %d genes\n", nrow(atlas)))
  atlas_lookup <- atlas[!duplicated(human_symbol)]
  setnames(atlas_lookup, "human_symbol", "symbol")
  setkey(atlas_lookup, symbol)
} else {
  cat("  WARNING: Multi-evidence atlas not found, skipping\n")
  atlas_lookup <- data.table(symbol = character(), layers_active = integer())
  setkey(atlas_lookup, symbol)
}

# ==============================================================================
# Part 5: Load and Merge COLOC Files (4 cell types)
# ==============================================================================
cat("\n--- Loading Zenodo COLOC files ---\n")
coloc_list <- list()

for (ct in CELL_TYPES) {
  f <- file.path(COLOC_DIR, sprintf("eGenes_coloc_list_%s.txt.gz", ct))
  if (!file.exists(f)) {
    cat(sprintf("  WARNING: Missing COLOC file for %s\n", ct))
    next
  }
  dt <- fread(f, header = TRUE, sep = "\t")
  dt[, cell_type := ct]
  cat(sprintf("  %s: %d eGenes\n", ct, nrow(dt)))
  coloc_list[[ct]] <- dt
}

coloc_all <- rbindlist(coloc_list, use.names = TRUE, fill = TRUE)
cat(sprintf("  Total eGene-celltype records: %d\n", nrow(coloc_all)))
cat(sprintf("  Unique genes across all cell types: %d\n",
            uniqueN(coloc_all$gene)))

# ==============================================================================
# Part 6: Extract NAFLD-Colocalizing Genes and Build Summary
# ==============================================================================
cat("\n--- Extracting NAFLD-colocalizing genes ---\n")

# Identify the is_liver_eGene column (varies by cell type)
egene_cols <- grep("^is_liver_eGene_", names(coloc_all), value = TRUE)

# Identify GWAS trait columns (all boolean columns except gene, cell_type,
# is_liver_eGene_*)
meta_cols <- c("gene", "cell_type", egene_cols)
trait_cols <- setdiff(names(coloc_all), meta_cols)
cat(sprintf("  GWAS trait columns: %d\n", length(trait_cols)))

# Mark NAFLD colocalization
coloc_all[, nafld_coloc := get(NAFLD_COL) == TRUE]

# Count other colocalizing traits (exclude NAFLD column itself)
other_trait_cols <- setdiff(trait_cols, NAFLD_COL)
coloc_all[, n_traits_coloc := rowSums(.SD == TRUE, na.rm = TRUE),
          .SDcols = other_trait_cols]

# Build trait_list: concatenate names of TRUE trait columns (excluding NAFLD)
coloc_all[, trait_list := apply(.SD, 1, function(row) {
  hits <- other_trait_cols[row == TRUE]
  if (length(hits) == 0) return(NA_character_)
  paste(hits, collapse = ";")
}), .SDcols = other_trait_cols]

# Cross-reference with dream DEGs
coloc_all[, is_deg := gene %in% deg_symbols]

# Merge dream LFC / padj
coloc_merged <- merge(
  coloc_all[, .(gene, cell_type, nafld_coloc, n_traits_coloc, trait_list, is_deg)],
  dream_lookup, by.x = "gene", by.y = "symbol", all.x = TRUE
)

# Merge consensus tier
coloc_merged <- merge(
  coloc_merged, tier_lookup, by.x = "gene", by.y = "symbol", all.x = TRUE
)
coloc_merged[is.na(consensus_tier), consensus_tier := "Not_significant"]

# Merge atlas layers_active
coloc_merged <- merge(
  coloc_merged, atlas_lookup, by.x = "gene", by.y = "symbol", all.x = TRUE
)
coloc_merged[, in_atlas := !is.na(layers_active)]

# Order by NAFLD colocalization, then DEG status, then n_traits
setorder(coloc_merged, -nafld_coloc, -is_deg, -n_traits_coloc)

# Write output
out_coloc <- file.path(RESULTS_DIR, "zenodo_coloc_integrated.csv")
fwrite(coloc_merged, out_coloc)
cat(sprintf("  Wrote: %s (%d rows)\n", out_coloc, nrow(coloc_merged)))

# ==============================================================================
# Part 7: Summary Statistics
# ==============================================================================
cat("\n--- COLOC Summary ---\n")

# NAFLD-colocalizing genes per cell type
nafld_summary <- coloc_merged[nafld_coloc == TRUE, .N, by = cell_type]
cat("  NAFLD-colocalizing eGenes per cell type:\n")
for (i in seq_len(nrow(nafld_summary))) {
  cat(sprintf("    %s: %d\n", nafld_summary$cell_type[i], nafld_summary$N[i]))
}

nafld_genes <- coloc_merged[nafld_coloc == TRUE]
n_nafld_unique <- uniqueN(nafld_genes$gene)
n_nafld_deg <- uniqueN(nafld_genes[is_deg == TRUE, gene])
cat(sprintf("\n  Total unique NAFLD-colocalizing genes: %d\n", n_nafld_unique))
cat(sprintf("  Of which are MASLD DEGs (padj < 0.1): %d (%.1f%%)\n",
            n_nafld_deg, 100 * n_nafld_deg / max(n_nafld_unique, 1)))

# Tier breakdown for NAFLD-colocalizing DEGs
tier_breakdown <- nafld_genes[is_deg == TRUE,
  .(n_genes = uniqueN(gene)), by = consensus_tier]
if (nrow(tier_breakdown) > 0) {
  cat("  Consensus tier breakdown (NAFLD COLOC + DEG):\n")
  for (i in seq_len(nrow(tier_breakdown))) {
    cat(sprintf("    %s: %d\n", tier_breakdown$consensus_tier[i],
                tier_breakdown$n_genes[i]))
  }
}

# Highlight known genes
cat("\n  Known genes of interest in COLOC results:\n")
for (g in HIGHLIGHT_GENES) {
  rows <- coloc_merged[gene == g]
  if (nrow(rows) > 0) {
    cts <- paste(rows$cell_type, collapse = ", ")
    is_nafld <- any(rows$nafld_coloc)
    is_d <- any(rows$is_deg)
    cat(sprintf("    %s: cell_types=[%s], nafld_coloc=%s, is_deg=%s\n",
                g, cts, is_nafld, is_d))
  }
}

# Multi-cell-type colocalization
multi_ct <- coloc_merged[nafld_coloc == TRUE, .(n_cell_types = uniqueN(cell_type)),
                          by = gene][n_cell_types > 1]
if (nrow(multi_ct) > 0) {
  cat(sprintf("\n  Genes colocalizing with NAFLD in >1 cell type: %d\n",
              nrow(multi_ct)))
  setorder(multi_ct, -n_cell_types)
  cat("  Top multi-cell-type genes:\n")
  for (i in seq_len(min(10, nrow(multi_ct)))) {
    cat(sprintf("    %s: %d cell types\n",
                multi_ct$gene[i], multi_ct$n_cell_types[i]))
  }
}

# ==============================================================================
# Part 8: Perturbation Validation
# ==============================================================================
cat("\n--- Loading perturbation data ---\n")

# Total number of testable genes (dream universe size) for Fisher's test
n_universe <- nrow(dream_with_sym)

perturb_results <- list()

for (exp_name in names(PERTURB_FILES)) {
  info <- PERTURB_FILES[[exp_name]]
  f <- info$file
  gene_col <- info$gene_col
  label <- info$label

  if (!file.exists(f)) {
    cat(sprintf("  WARNING: Missing perturbation file: %s\n", f))
    next
  }

  cat(sprintf("  Loading %s...\n", label))
  pdt <- fread(f, header = TRUE, sep = "\t")
  cat(sprintf("    %d DEGs in perturbation experiment\n", nrow(pdt)))

  # Standardize gene column name
  if (gene_col != "gene") {
    setnames(pdt, gene_col, "gene")
  }

  # Find overlap with MASLD DEGs
  pdt[, is_masld_deg := gene %in% deg_symbols]
  overlap_dt <- merge(
    pdt[, .(gene, perturbation_logFC = log2FoldChange,
            perturbation_padj = padj, is_masld_deg)],
    dream_lookup, by.x = "gene", by.y = "symbol", all.x = TRUE
  )

  # Add experiment label
  overlap_dt[, experiment := label]

  # Compute direction concordance: same sign of LFC in perturbation and MASLD
  overlap_dt[, direction_concordant := fifelse(
    is.na(dream_logFC) | is.na(perturbation_logFC), NA,
    sign(perturbation_logFC) == sign(dream_logFC)
  )]

  perturb_results[[exp_name]] <- overlap_dt

  # Summary for this experiment
  n_overlap <- sum(overlap_dt$is_masld_deg, na.rm = TRUE)
  n_perturb <- nrow(pdt)
  n_concordant <- sum(overlap_dt$direction_concordant == TRUE, na.rm = TRUE)
  n_testable <- sum(overlap_dt$is_masld_deg & !is.na(overlap_dt$direction_concordant),
                    na.rm = TRUE)
  cat(sprintf("    Overlap with MASLD DEGs: %d / %d perturbation DEGs\n",
              n_overlap, n_perturb))
  if (n_testable > 0) {
    cat(sprintf("    Direction concordance: %d / %d (%.1f%%)\n",
                n_concordant, n_testable, 100 * n_concordant / n_testable))
  }

  # Fisher's exact test: enrichment of MASLD DEGs among perturbation DEGs
  # 2x2 contingency table:
  #                   MASLD DEG   Not MASLD DEG
  # Perturb DEG         a            b
  # Not Perturb DEG     c            d
  perturb_genes <- unique(pdt$gene)
  # Restrict to genes present in dream universe
  perturb_in_universe <- perturb_genes[perturb_genes %in% dream_lookup$symbol]
  a <- length(intersect(perturb_in_universe, deg_symbols))
  b <- length(perturb_in_universe) - a
  c_val <- length(deg_symbols) - a
  d <- nrow(dream_lookup) - a - b - c_val

  if (a > 0 && d > 0) {
    ft <- fisher.test(matrix(c(a, b, c_val, d), nrow = 2),
                      alternative = "greater")
    cat(sprintf("    Fisher's exact test: OR = %.2f, p = %.2e\n",
                ft$estimate, ft$p.value))
  }
}

# Combine all perturbation results
if (length(perturb_results) > 0) {
  perturb_combined <- rbindlist(perturb_results, use.names = TRUE, fill = TRUE)

  # Reorder columns for output
  out_cols <- c("experiment", "gene", "perturbation_logFC", "perturbation_padj",
                "is_masld_deg", "dream_logFC", "dream_padj", "direction_concordant")
  setcolorder(perturb_combined, intersect(out_cols, names(perturb_combined)))

  out_perturb <- file.path(RESULTS_DIR, "zenodo_perturbation_validation.csv")
  fwrite(perturb_combined, out_perturb)
  cat(sprintf("\n  Wrote: %s (%d rows)\n", out_perturb, nrow(perturb_combined)))

  # Overall perturbation summary
  cat("\n--- Perturbation Validation Summary ---\n")
  summary_dt <- perturb_combined[, .(
    n_perturb_degs  = .N,
    n_masld_overlap = sum(is_masld_deg, na.rm = TRUE),
    n_concordant    = sum(direction_concordant == TRUE, na.rm = TRUE),
    n_discordant    = sum(direction_concordant == FALSE, na.rm = TRUE)
  ), by = experiment]
  print(summary_dt)
} else {
  cat("  No perturbation files loaded\n")
}

# ==============================================================================
# Final summary
# ==============================================================================
cat("\n=== Script 37 Complete ===\n")
cat("End time:", format(Sys.time()), "\n")
cat(sprintf("Outputs:\n  %s\n  %s\n",
  file.path(RESULTS_DIR, "zenodo_coloc_integrated.csv"),
  file.path(RESULTS_DIR, "zenodo_perturbation_validation.csv")
))
