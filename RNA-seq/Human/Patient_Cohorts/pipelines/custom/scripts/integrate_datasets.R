#!/usr/bin/env Rscript

# Integrated MASLD Analysis Script
# Merges counts from GSE167523, GSE126848, and PRJNA512027
# Performs batch correction and unified DE analysis

library(DESeq2)
library(dplyr)
library(sva) # For ComBat-seq
library(ggplot2)
library(pheatmap)

# --- Configuration ---
DIRS <- list(
  GSE167523 = "results/GSE167523/counts/featurecounts/gene_counts.txt",
  GSE126848 = "results/GSE126848/counts/featurecounts/gene_counts.txt",
  PRJNA512027 = "results/PRJNA512027/counts/featurecounts/gene_counts.txt"
)

META_FILES <- list(
  GSE167523 = "pipelines/custom/GSE167523/metadata/samples.tsv",
  GSE126848 = "pipelines/custom/GSE126848/metadata/samples.tsv",
  PRJNA512027 = "pipelines/custom/PRJNA512027/metadata/SraRunTable.csv"
)

OUTPUT_DIR <- "results/integrated_analysis"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# --- Helper Functions ---

load_counts <- function(filepath, dataset_name) {
  if (!file.exists(filepath)) {
    warning(paste("Count file not found:", filepath))
    return(NULL)
  }
  counts <- read.table(filepath, header=TRUE, row.names=1, sep="\t", comment.char="#")
  # Remove annotation cols (Chr, Start, End, Strand, Length) - usually first 5
  if (ncol(counts) > 5) {
    counts <- counts[, 6:ncol(counts)]
  }
  # Clean colnames
  colnames(counts) <- gsub("_Aligned.sortedByCoord.out.bam", "", basename(colnames(counts)))
  colnames(counts) <- gsub(".bam", "", colnames(counts))
  # Add dataset prefix to avoid collisions (though SRR should be unique)
  # colnames(counts) <- paste(dataset_name, colnames(counts), sep="_") 
  return(counts)
}

# --- 1. Load Data ---

cat("Loading count matrices...\n")
counts_list <- list()
genes_list <- list()

for (d in names(DIRS)) {
  mat <- load_counts(DIRS[[d]], d)
  if (!is.null(mat)) {
    counts_list[[d]] <- mat
    genes_list[[d]] <- rownames(mat)
    cat("Loaded", d, ":", ncol(mat), "samples,", nrow(mat), "genes\n")
  }
}

# Intersect Genes
common_genes <- Reduce(intersect, genes_list)
cat("Common Genes:", length(common_genes), "\n")

# Merge Matrices
merged_counts <- do.call(cbind, lapply(counts_list, function(x) x[common_genes, ]))
cat("Merged Matrix:", nrow(merged_counts), "genes x", ncol(merged_counts), "samples\n")

# --- 2. Harmonize Metadata ---

cat("Harmonizing metadata...\n")

# We need a unified metadata frame with: SampleID, Batch, Condition
# Strategies for mapping conditions:
# GSE167523: subtype (NAFL, NASH) -> MASLD
# GSE126848: condition (healthy, obese, NAFLD, NASH) -> healthy=Control, others=MASLD
# PRJNA512027: disease_stage (NORMAL -> Control, others -> MASLD)

meta_rows <- list()

# GSE167523
meta1 <- read.table(META_FILES$GSE167523, header=TRUE, sep="\t")
for (i in 1:nrow(meta1)) {
  sid <- meta1$sample_id[i]
  cond <- meta1$subtype[i]
  # No healthy in this dataset
  unified <- "MASLD" # NAFL or NASH
  meta_rows[[sid]] <- data.frame(sample_id=sid, batch="GSE167523", original=cond, condition=unified, stringsAsFactors=FALSE)
}

# GSE126848
meta2 <- read.table(META_FILES$GSE126848, header=TRUE, sep="\t")
for (i in 1:nrow(meta2)) {
  sid <- meta2$sample_id[i]
  cond <- meta2$condition[i]
  unified <- ifelse(cond %in% c("healthy"), "Control", "MASLD") # Treat obese as MASLD? or exclude? 
  # Let's check user intent. "MASLD vs Control". Usually 'obese' is a confounding control.
  # For now, map 'healthy' to Control, 'NAFLD'/'NASH' to MASLD. 'obese' -> Exclude or separate?
  # Decision: Map 'obese' to 'Obese_Control' class to allow flexible contrast, or exclude.
  # Safest: 'Obese' is NOT 'Healthy Control'.
  if (cond == "obese") unified <- "Obese"
  meta_rows[[sid]] <- data.frame(sample_id=sid, batch="GSE126848", original=cond, condition=unified, stringsAsFactors=FALSE)
}

# PRJNA512027 (CSV)
meta3 <- read.csv(META_FILES$PRJNA512027)
for (i in 1:nrow(meta3)) {
  sid <- meta3$Run[i]
  cond <- meta3$disease_stage[i]
  # disease_stage: NORMAL, Lob Inflam *, Fibrosis *, STEATOSIS *
  unified <- ifelse(cond == "NORMAL", "Control", "MASLD")
  meta_rows[[sid]] <- data.frame(sample_id=sid, batch="PRJNA512027", original=cond, condition=unified, stringsAsFactors=FALSE)
}

metadata <- do.call(rbind, meta_rows)
rownames(metadata) <- metadata$sample_id

# Intersect with counts
valid_samples <- intersect(colnames(merged_counts), rownames(metadata))
metadata <- metadata[valid_samples, ]
merged_counts <- merged_counts[, valid_samples]

# Check balance
print(table(metadata$batch, metadata$condition))

# --- 3. Batch Correction (ComBat-seq) ---
# ComBat-seq is designed for count data (returns "adjusted counts")
# However, DESeq2 prefers creating a design with `~ batch + condition`.
# We will do BOTH: 
# 1. ComBat-seq adjusted counts for PCA/Visualization
# 2. Raw counts + batch design for DE Analysis

cat("Running ComBat-seq for visualization...\n")
adjusted_counts <- ComBat_seq(as.matrix(merged_counts), batch=metadata$batch, group=metadata$condition)

# --- 4. PCA Plot (Adjusted) ---
vsd <- vst(DESeqDataSetFromMatrix(adjusted_counts, metadata, ~condition), blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("batch", "condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

p_pca <- ggplot(pcaData, aes(PC1, PC2, color=condition, shape=batch)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
  theme_bw() +
  ggtitle("PCA (Batch Corrected)")
ggsave(file.path(OUTPUT_DIR, "pca_corrected.pdf"), p_pca)

# --- 5. Differential Expression (DESeq2) ---
cat("Running DESeq2 (Raw Counts + Batch Design)...\n")

metadata$batch <- as.factor(metadata$batch)
metadata$condition <- as.factor(metadata$condition)
# Set reference
metadata$condition <- relevel(metadata$condition, ref="Control")

dds <- DESeqDataSetFromMatrix(countData = merged_counts,
                              colData = metadata,
                              design = ~ batch + condition)

dds <- DESeq(dds)
res <- results(dds, contrast=c("condition", "MASLD", "Control"))

# Save Results
write.csv(as.data.frame(res), file.path(OUTPUT_DIR, "deseq2_results_MASLD_vs_Control.csv"))
saveRDS(dds, file.path(OUTPUT_DIR, "deseq2_dds.rds"))

# Summary Sig
sig <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
cat("Significant Genes (padj<0.05, LFC>1):", nrow(sig), "\n")
write.csv(as.data.frame(sig), file.path(OUTPUT_DIR, "significant_genes.csv"))

cat("Integration Complete!\n")
