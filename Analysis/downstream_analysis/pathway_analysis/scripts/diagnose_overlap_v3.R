
# Diagnostic v3 - No suppression
cat("Starting diagnostic...\n")

library(msigdbr)
cat("Loaded msigdbr\n")
library(org.Mm.eg.db)
cat("Loaded org.Mm.eg.db\n")
library(AnnotationDbi)
cat("Loaded AnnotationDbi\n")
library(readr)
cat("Loaded readr\n")

# Load Hallmark
hallmark_gs <- msigdbr(species = "Homo sapiens", collection = "H")
hallmark_genes <- unique(hallmark_gs$gene_symbol)
cat(sprintf("Hallmark Genes: %d\n", length(hallmark_genes)))

check_file <- function(name, path) {
    cat(sprintf("\n--- %s ---\n", name))
    if (!file.exists(path)) {
        cat("  File not found.\n")
        return()
    }
    
    # Use base R for safety
    counts <- read.csv(path, nrows=5000, header=TRUE)
    ids <- counts[,1]
    
    cat("  IDs Sample: ", paste(head(ids, 3), collapse=", "), "\n")
    
    is_ensembl <- any(grepl("^ENS", ids[1:10]))
    cat(sprintf("  Is Ensembl: %s\n", is_ensembl))
    
    symbols <- as.character(ids)
    if (is_ensembl) {
        clean_ids <- sub("\\..*", "", ids)
        mapped <- mapIds(org.Mm.eg.db, keys = clean_ids, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
        symbols <- as.character(mapped)
    }
    
    valid_symbols <- symbols[!is.na(symbols)]
    cat(sprintf("  Valid Symbols: %d / %d\n", length(valid_symbols), length(ids)))
    
    upper_symbols <- toupper(valid_symbols)
    overlap <- intersect(upper_symbols, hallmark_genes)
    
    cat(sprintf("  Overlap (toupper): %d (%.1f%% of valid input)\n", length(overlap), 100 * length(overlap)/length(valid_symbols)))
    cat(sprintf("  Coverage of Hallmark: %.1f%%\n", 100 * length(overlap)/length(hallmark_genes)))
}

base <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq"
check_file("InHouse_MCD", file.path(base, "in-house_MCD_RNAseq/normalized_counts_all_samples.csv"))
