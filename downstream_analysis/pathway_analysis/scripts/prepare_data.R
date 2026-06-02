#!/usr/bin/env Rscript
# prepare_data.R - Prepare gene lists for pathway analysis
# 
# OUTPUTS:
# 1. core_deg_symbols.txt      - DEG symbols for ORA
# 2. background_universe.txt   - All detected genes (for ORA universe)
# 3. gsea_ranked_full.rnk      - FULL transcriptome ranked list for GSEA
# 4. core_degs_with_lfc.csv    - Detailed DEG data for reference

suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
})

# ============================================================================
# Configuration
# ============================================================================
BASE_DIR <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
STREAMLIT_DIR <- file.path(BASE_DIR, "streamlit_deg_explorer")
PATHWAY_DIR <- file.path(BASE_DIR, "downstream_analysis/pathway_analysis")
OUTPUT_DIR <- file.path(PATHWAY_DIR, "data/processed")

# Input files
MASTER_MATRIX <- file.path(STREAMLIT_DIR, "master_ortholog_matrix.csv.gz")
CORE_DEGS <- file.path(BASE_DIR, "final_core_degs.csv")

# Create output directory
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("=== Preparing Data for Pathway Analysis ===\n\n")

# ============================================================================
# Load Data
# ============================================================================
cat("[1/6] Loading master ortholog matrix...\n")
master_df <- fread(MASTER_MATRIX) %>% as_tibble()
cat(sprintf("  Loaded %d genes with expression data\n", nrow(master_df)))

cat("[2/6] Loading core DEGs...\n")
core_degs <- read_csv(CORE_DEGS, show_col_types = FALSE)
cat(sprintf("  Loaded %d core DEGs\n", nrow(core_degs)))

# ============================================================================
# Compute Aggregate Statistics for ALL Genes
# ============================================================================
cat("[3/6] Computing aggregate statistics for full transcriptome...\n")

# Extract LFC columns
lfc_cols <- grep("_lfc$", names(master_df), value = TRUE)
cat(sprintf("  Found %d LFC columns: %s\n", length(lfc_cols), paste(lfc_cols, collapse = ", ")))

# Compute aggregate LFC for ALL genes (not just DEGs)
gene_stats <- master_df %>%
    select(human_id, Symbol, all_of(lfc_cols)) %>%
    rowwise() %>%
    mutate(
        n_datasets = sum(!is.na(c_across(all_of(lfc_cols)))),
        mean_lfc = mean(c_across(all_of(lfc_cols)), na.rm = TRUE),
        max_abs_lfc = max(abs(c_across(all_of(lfc_cols))), na.rm = TRUE),
        sum_lfc = sum(c_across(all_of(lfc_cols)), na.rm = TRUE)
    ) %>%
    ungroup() %>%
    # Filter invalid values but keep ALL genes with valid data
    filter(!is.na(mean_lfc) & !is.infinite(mean_lfc)) %>%
    filter(!is.na(max_abs_lfc) & !is.infinite(max_abs_lfc)) %>%
    # Require at least 1 dataset with data
    filter(n_datasets >= 1)

cat(sprintf("  Computed stats for %d genes (full transcriptome)\n", nrow(gene_stats)))

# ============================================================================
# Create FULL Ranked Gene List for GSEA (Integrative + Individual)
# ============================================================================
cat("[4/6] Creating ranked gene lists for GSEA...\n")

# 1. Integrative Ranked List (Mean LFC across all datasets)
# For GSEA: use ALL genes from the transcriptome, ranked by mean_lfc
gsea_full_ranked <- gene_stats %>%
    filter(!is.na(Symbol) & Symbol != "") %>%
    # Handle duplicates by taking max absolute LFC
    group_by(Symbol) %>%
    slice_max(order_by = abs(mean_lfc), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    # Ranking metric: mean_lfc (preserves direction)
    mutate(gsea_rank = mean_lfc) %>%
    arrange(desc(gsea_rank))

# Save Integrative ranked list
rank_file_full <- file.path(OUTPUT_DIR, "gsea_ranked_full.rnk")
gsea_full_ranked %>%
    select(Symbol, gsea_rank) %>%
    write_tsv(rank_file_full, col_names = FALSE)

cat(sprintf("  Saved INTEGRATIVE ranked list: %s (%d genes)\n", 
            basename(rank_file_full), nrow(gsea_full_ranked)))

# 2. Individual Dataset Ranked Lists
cat("  Generating per-dataset ranked lists...\n")

for (col in lfc_cols) {
    # Extract dataset name from column (e.g., "GSE130970_lfc" -> "GSE130970")
    dataset_name <- sub("_lfc$", "", col)
    
    # Create ranked list for this dataset
    # We use the raw master_df here to get values for this specific column
    dataset_ranked <- master_df %>%
        select(Symbol, all_of(col)) %>%
        rename(lfc = !!col) %>%
        filter(!is.na(Symbol) & Symbol != "") %>%
        filter(!is.na(lfc)) %>%
        group_by(Symbol) %>%
        # If duplicates exist, pick the one with max abs LFC
        slice_max(order_by = abs(lfc), n = 1, with_ties = FALSE) %>%
        ungroup() %>%
        arrange(desc(lfc))
    
    # Save
    out_file <- file.path(OUTPUT_DIR, paste0("gsea_ranked_", dataset_name, ".rnk"))
    dataset_ranked %>%
        select(Symbol, lfc) %>%
        write_tsv(out_file, col_names = FALSE)
        
    cat(sprintf("    -> %s (%d genes)\n", basename(out_file), nrow(dataset_ranked)))
}

# ============================================================================
# Create Background Universe for ORA
# ============================================================================
cat("[5/6] Creating background universe for ORA...\n")

# The universe = all genes that were actually measured/detected
universe_file <- file.path(OUTPUT_DIR, "background_universe.txt")
background_genes <- gsea_full_ranked$Symbol
writeLines(background_genes, universe_file)

cat(sprintf("  Saved background universe: %s (%d genes)\n", 
            basename(universe_file), length(background_genes)))

# ============================================================================
# Extract Core DEGs for ORA
# ============================================================================
cat("[6/6] Preparing core DEG list for ORA...\n")

# Use separate_rows to explode the human_ortholog_ids column (semicolon separated)
# final_core_degs.csv: mouse_gene_id, ..., human_ortholog_ids, human_ortholog_symbols, ...
core_degs_expanded <- core_degs %>%
  filter(!is.na(human_ortholog_ids) & human_ortholog_ids != "") %>%
  separate_rows(human_ortholog_ids, human_ortholog_symbols, sep = ";") %>%
  mutate(
      human_id = str_trim(human_ortholog_ids),
      human_symbol = str_trim(human_ortholog_symbols)
  )

# Join with core DEGs to get overlap_count for the DEG subset
# We join based on the extracted Human ID
merged_df <- core_degs_expanded %>%
    select(gene_id = human_id, gene_symbol = human_symbol, overlap_count = n_analyses, analyses = source_analyses) %>%
    inner_join(
        gene_stats %>% select(human_id, Symbol, mean_lfc, max_abs_lfc, sum_lfc, n_datasets),
        by = c("gene_id" = "human_id")
    )

cat(sprintf("  Matched %d core DEGs with LFC data\n", nrow(merged_df)))

# Prepare DEG list for ORA
ora_genes <- merged_df %>%
    filter(!is.na(gene_symbol) & gene_symbol != "") %>%
    group_by(gene_symbol) %>%
    slice_max(order_by = overlap_count, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
        gsea_rank = mean_lfc,
        gsea_rank_weighted = sum_lfc * log2(overlap_count + 1)
    ) %>%
    arrange(desc(gsea_rank))

# Save processed DEG data (for reference)
output_file <- file.path(OUTPUT_DIR, "core_degs_with_lfc.csv")
write_csv(ora_genes, output_file)
cat(sprintf("  Saved detailed DEG data: %s\n", basename(output_file)))

# Save DEG symbols for ORA
genelist_file <- file.path(OUTPUT_DIR, "core_deg_symbols.txt")
writeLines(ora_genes$gene_symbol, genelist_file)
cat(sprintf("  Saved DEG list for ORA: %s (%d genes)\n", 
            basename(genelist_file), nrow(ora_genes)))

# ============================================================================
# Summary Statistics
# ============================================================================
cat("\n=== Summary ===\n")
cat(sprintf("Full transcriptome (GSEA input): %d genes\n", nrow(gsea_full_ranked)))
cat(sprintf("Background universe (ORA):       %d genes\n", length(background_genes)))
cat(sprintf("Core DEGs (ORA input):           %d genes\n", nrow(ora_genes)))
cat(sprintf("\nGSEA rank range: [%.2f, %.2f]\n", 
            min(gsea_full_ranked$gsea_rank), max(gsea_full_ranked$gsea_rank)))
cat(sprintf("DEGs with positive LFC: %d\n", sum(ora_genes$gsea_rank > 0)))
cat(sprintf("DEGs with negative LFC: %d\n", sum(ora_genes$gsea_rank < 0)))

cat("\n=== Data Preparation Complete ===\n")
cat("\nOutput files:\n")
cat("  - gsea_ranked_full.rnk    : Full transcriptome for GSEA\n")
cat("  - background_universe.txt : Background for ORA enricher()\n")
cat("  - core_deg_symbols.txt    : DEG list for ORA\n")
cat("  - core_degs_with_lfc.csv  : Detailed DEG data\n")
