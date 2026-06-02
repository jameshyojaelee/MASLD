
# Diagnostic Script: Check Gene Overlap for ssGSEA
suppressPackageStartupMessages({
    library(tidyverse)
    library(msigdbr)
})

# Load Hallmark Headers
hallmark_gs <- msigdbr(species = "Homo sapiens", collection = "H")
hallmark_genes <- unique(hallmark_gs$gene_symbol)
cat(sprintf("Total Hallmark Genes: %d\n", length(hallmark_genes)))

# Define Process Logic (Minified from run_refactored_ssgsea.R)
process_counts_preview <- function(file_path) {
    if (!file.exists(file_path)) return(NULL)
    counts <- data.table::fread(file_path)
    # Basic cleaning
    if (colnames(counts)[1] == "V1" || colnames(counts)[1] == "") colnames(counts)[1] <- "gene_id"
    if (!"gene_id" %in% colnames(counts)) colnames(counts)[1] <- "gene_id" 
    
    ids <- counts$gene_id
    # Simple check for now
    head(ids, 10)
}

# Check Datasets
base_dir <- "/gpfs/commons/groups/sanjana_lab/Cas13"
rnaseq_dir <- file.path(base_dir, "RNA-seq")
files <- list(
    "InHouse_MCD" = file.path(rnaseq_dir, "in-house_MCD_RNAseq/normalized_counts_all_samples.csv"),
    "Ext_MCD_1" = file.path(rnaseq_dir, "other_MCD_RNAseq/GSE156918/analysis_mcd_vs_control/normalized_counts.csv")
)

for (n in names(files)) {
    cat(sprintf("\n--- Checking %s ---\n", n))
    f <- files[[n]]
    ids <- process_counts_preview(f)
    
    if (is.null(ids)) {
        cat("  File not found.\n")
        next
    }
    
    cat("  First 5 IDs: ", paste(ids[1:5], collapse=", "), "\n")
    
    # Simulate the pipeline's 'toupper' logic
    mapped_ids <- toupper(ids)
    
    overlap <- intersect(mapped_ids, hallmark_genes)
    pct <- length(overlap) / length(hallmark_genes) * 100
    cat(sprintf("  Overlap with Hallmark (via toupper): %d genes (%.1f%%)\n", length(overlap), pct))
    
    # Check if we missed likely orthologs
    # Compare against proper ortholog mapping logic?
    # This diagnostic confirms if 'toupper' is the bottleneck
}
