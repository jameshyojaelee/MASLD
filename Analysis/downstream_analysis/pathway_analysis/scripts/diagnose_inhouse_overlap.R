
suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
    library(org.Mm.eg.db)
    library(msigdbr)
})

# 1. Load Hallmark Mouse
hallmark_mouse <- msigdbr(species = "Mus musculus", collection = "H")
hallmark_genes <- unique(hallmark_mouse$gene_symbol)
cat(sprintf("Mouse Hallmark Genes: %d (e.g., %s)\n", length(hallmark_genes), paste(head(hallmark_genes), collapse=", ")))

# 2. Process InHouse File manually (replicating logic)
f <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/in-house_MCD_RNAseq/normalized_counts_all_samples.csv"
cat(sprintf("\nProcessing %s...\n", basename(f)))

counts <- fread(f)
# Rename first col
colnames(counts)[1] <- "gene_id"
ids <- counts$gene_id
counts$gene_id_clean <- sub("\\..*", "", ids)

cat(sprintf("Input rows: %d\n", nrow(counts)))
cat(sprintf("First 5 IDs: %s\n", paste(head(counts$gene_id_clean), collapse=", ")))

# Map
cat("Mapping Ensembl to SYMBOL...\n")
mapping <- mapIds(org.Mm.eg.db, keys = counts$gene_id_clean, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
counts$Symbol <- mapping

# Check NAs
n_na <- sum(is.na(counts$Symbol))
cat(sprintf("Unmapped IDs: %d (%.1f%%)\n", n_na, 100*n_na/nrow(counts)))

# Check Overlap
mapped_symbols <- unique(na.omit(counts$Symbol))
cat(sprintf("Mapped Unique Symbols: %d(e.g., %s)\n", length(mapped_symbols), paste(head(mapped_symbols), collapse=", ")))

overlap <- intersect(mapped_symbols, hallmark_genes)
cat(sprintf("\nOverlap with Hallmark: %d genes\n", length(overlap)))
cat(sprintf("Hallmark Coverage: %.1f%%\n", 100 * length(overlap) / length(hallmark_genes)))

# Check Case Sensitivity
cat("\nChecking Case Mismatch:\n")
cat(sprintf("  Matrix (head): %s\n", paste(head(mapped_symbols), collapse=", ")))
cat(sprintf("  Hallmark (head): %s\n", paste(head(hallmark_genes), collapse=", ")))

upper_matrix <- toupper(mapped_symbols)
upper_hallmark <- toupper(hallmark_genes)
overlap_upper <- intersect(upper_matrix, upper_hallmark)
cat(sprintf("  Overlap (Case Insensitive): %d genes\n", length(overlap_upper)))

# 3. Check Data Distribution (Log vs Normalized)
# Columns 4 to end are samples
sample_cols <- colnames(counts)[4:ncol(counts)]
vals <- as.matrix(counts[, ..sample_cols])
vals <- as.numeric(vals)
vals <- na.omit(vals)

cat("\nData Distribution Summary:\n")
print(summary(vals))
cat(sprintf("  Max value: %.2f\n", max(vals)))
if (max(vals) > 100) {
    cat("  Likely Linear Scale (needs log1p)\n")
} else {
    cat("  Likely Log Scale\n")
}
