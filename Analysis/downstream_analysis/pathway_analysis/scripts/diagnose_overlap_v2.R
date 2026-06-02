
suppressPackageStartupMessages({
    library(tidyverse)
    library(msigdbr)
    library(org.Mm.eg.db) # Ensure this is loaded
    library(AnnotationDbi)
})

# Load Hallmark
hallmark_gs <- msigdbr(species = "Homo sapiens", collection = "H")
hallmark_genes <- unique(hallmark_gs$gene_symbol)
cat(sprintf("Hallmark Genes: %d\n", length(hallmark_genes)))

# Function to check file
check_file <- function(name, path) {
    cat(sprintf("\n--- %s ---\n", name))
    if (!file.exists(path)) {
        cat("  File not found.\n")
        return()
    }
    
    counts <- read_csv(path, n_max = 5000, show_col_types = FALSE) # Read first 5000 rows
    cat("  Cols: ", paste(head(colnames(counts), 3), collapse=", "), "...\n")
    
    # Identify ID column (logic from script)
    ids <- counts[[1]] 
    if ("gene_id" %in% colnames(counts)) ids <- counts$gene_id
    
    cat("  IDs Sample: ", paste(head(ids, 3), collapse=", "), "\n")
    
    is_ensembl <- any(grepl("^ENS", ids[1:10]))
    cat(sprintf("  Is Ensembl: %s\n", is_ensembl))
    
    symbols <- ids
    if (is_ensembl) {
        clean_ids <- sub("\\..*", "", ids)
        # Try Mapping
        mapped <- mapIds(org.Mm.eg.db, keys = clean_ids, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
        symbols <- as.character(mapped)
    }
    
    valid_symbols <- symbols[!is.na(symbols)]
    cat(sprintf("  Valid Symbols: %d / %d\n", length(valid_symbols), length(ids)))
    
    # Naive Mapping
    upper_symbols <- toupper(valid_symbols)
    overlap <- intersect(upper_symbols, hallmark_genes)
    
    cat(sprintf("  Overlap (toupper): %d (%.1f%% of valid input)\n", length(overlap), 100 * length(overlap)/length(valid_symbols)))
    cat(sprintf("  Coverage of Hallmark: %.1f%%\n", 100 * length(overlap)/length(hallmark_genes)))
}

base <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq"
check_file("InHouse_MCD", file.path(base, "in-house_MCD_RNAseq/normalized_counts_all_samples.csv"))
check_file("Ext_MCD_1", file.path(base, "other_MCD_RNAseq/GSE156918/analysis_mcd_vs_control/normalized_counts.csv"))
