
# Diagnostic Script: Inspect ssGSEA Input Matrices
# Goal: Understand why InHouse_MCD scores are low.
# Checks:
# 1. Number of genes in matrix
# 2. Value distribution (Raw vs Log?)
# 3. Overlap with Hallmark Mouse
# 4. Rank distribution

suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
    library(org.Mm.eg.db)
    library(msigdbr)
    library(AnnotationDbi)
})

# 1. Load Data Function (Simplified from script)
process_counts_smart <- function(file_path, species = "human") {
    if (!file.exists(file_path)) return(NULL)
    counts <- fread(file_path)
    
    # Headers
    if (colnames(counts)[1] == "V1" || colnames(counts)[1] == "") colnames(counts)[1] <- "gene_id"
    if (!"gene_id" %in% colnames(counts)) colnames(counts)[1] <- "gene_id" 
    
    ids <- counts$gene_id
    is_ensembl <- any(grepl("^ENS", ids[1:10]))
    
    if (is_ensembl) {
        counts$gene_id <- sub("\\..*", "", counts$gene_id)
        db <- if (species == "human") org.Hs.eg.db else org.Mm.eg.db
        mapping <- mapIds(db, keys = counts$gene_id, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
        counts$Symbol <- mapping
    } else {
        counts$Symbol <- counts$gene_id
    }
    
    mat_clean <- counts %>%
        as_tibble() %>%
        filter(!is.na(Symbol) & Symbol != "") %>%
        dplyr::select(Symbol, everything(), -gene_id) %>%
        group_by(Symbol) %>%
        summarise(across(where(is.numeric), max), .groups="drop") %>%
        column_to_rownames("Symbol") %>%
        as.matrix()
    
    return(mat_clean)
}

# 2. Load Gene Sets
cat("Loading Hallmark Mouse...\n")
hallmark_mouse <- msigdbr(species = "Mus musculus", collection = "H")
hallmark_genes <- unique(hallmark_mouse$gene_symbol)
cat(sprintf("Hallmark Mouse Genes: %d\n", length(hallmark_genes)))

# 3. Inspect InHouse_MCD
path <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/in-house_MCD_RNAseq/normalized_counts_all_samples.csv"
cat(sprintf("\nProcessing %s...\n", basename(path)))

mat <- process_counts_smart(path, "mouse")
cat(sprintf("Matrix Dimension: %d genes x %d samples\n", nrow(mat), ncol(mat)))

# Stats
cat("Value Distribution (Sample 1):\n")
print(summary(mat[,1]))

zero_pct <- sum(mat == 0) / length(mat) * 100
cat(sprintf("Sparsity (%% Zeros): %.2f%%\n", zero_pct))

# Overlap
genes <- rownames(mat)
overlap <- intersect(genes, hallmark_genes)
cat(sprintf("Overlap with Hallmark: %d genes (%.2f%% of Hallmark)\n", length(overlap), 100 * length(overlap)/length(hallmark_genes)))

# Check for "Title Case" vs "Upper Case"
upper_genes <- toupper(genes)
upper_overlap <- intersect(upper_genes, toupper(hallmark_genes))
cat(sprintf("Overlap (Case-Insensitive): %d genes\n", length(upper_overlap)))

# Sample Ranks (ssGSEA relies on ranks)
# If the matrix is 90% zeros, then 90% of genes have tied ranks at the bottom.
# This pushes variable genes to the top fraction.
ranks <- rank(mat[,1])
cat("Rank Summary:\n")
print(summary(ranks))

# Compare with Ext_MCD_1
path2 <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/other_MCD_RNAseq/GSE156918/analysis_mcd_vs_control/normalized_counts.csv"
cat(sprintf("\nProcessing %s (Comparison)...\n", basename(path2)))
mat2 <- process_counts_smart(path2, "mouse")
cat(sprintf("Matrix Dimension: %d genes x %d samples\n", nrow(mat2), ncol(mat2)))
zero_pct2 <- sum(mat2 == 0) / length(mat2) * 100
cat(sprintf("Sparsity: %.2f%%\n", zero_pct2))
print(summary(mat2[,1]))

overlap2 <- intersect(rownames(mat2), hallmark_genes)
cat(sprintf("Overlap with Hallmark: %d genes (%.2f%%)\n", length(overlap2), 100 * length(overlap2)/length(hallmark_genes)))
