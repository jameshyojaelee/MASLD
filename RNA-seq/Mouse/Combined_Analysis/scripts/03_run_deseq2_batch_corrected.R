#!/usr/bin/env Rscript
# 03_run_deseq2_batch_corrected.R
# Run DESeq2 with batch correction for combined MCD vs Control analysis

suppressPackageStartupMessages({
  library(DESeq2)
  library(readr)
  library(dplyr)
  library(tibble)
  library(stringr)
})

root <- normalizePath(".")
metadata_path <- file.path(root, "metadata", "samples.tsv")
counts_path <- file.path(root, "counts", "combined_gene_counts.tsv")
out_dir <- file.path(root, "analysis_mcd_vs_control")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("=== DESeq2 Analysis with Batch Correction ===")

# Read metadata
samples <- read_tsv(metadata_path, col_types = cols()) %>%
  mutate(
    diet = factor(diet, levels = c("Control", "MCD")),
    batch = factor(batch)
  )
samples_df <- as.data.frame(samples)
rownames(samples_df) <- samples_df$sample_id

message("Samples: ", nrow(samples))
message("Diet: ", paste(levels(samples$diet), collapse = " vs "))
message("Batches: ", paste(levels(samples$batch), collapse = ", "))

# Read combined counts
counts_raw <- read_tsv(counts_path, col_types = cols())
count_matrix <- counts_raw %>%
  column_to_rownames("Geneid") %>%
  as.matrix()

# Ensure column order matches metadata
available_samples <- intersect(samples$sample_id, colnames(count_matrix))
samples_df <- samples_df[available_samples, ]
count_matrix <- count_matrix[, available_samples]

message("\nCount matrix: ", nrow(count_matrix), " genes × ", ncol(count_matrix), " samples")

# Filter low-count genes (at least 10 counts in 3+ samples)
keep <- rowSums(count_matrix >= 10) >= 3
count_matrix <- count_matrix[keep, ]
message("After filtering: ", nrow(count_matrix), " genes")

# Create DESeq2 dataset with batch in design
message("\nCreating DESeq2 dataset...")
message("Design: ~ batch + diet")

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = samples_df,
  design = ~ batch + diet
)

# Run DESeq2
message("Running DESeq2...")
dds <- DESeq(dds)

# Extract results for diet effect
res <- results(dds, contrast = c("diet", "MCD", "Control"))
message("\nResults summary:")
summary(res)

# Apply log2FC shrinkage with fallback
message("Applying shrinkage...")
res_shrunk <- tryCatch({
  # Try ashr first (more widely available)
  lfcShrink(dds, coef = "diet_MCD_vs_Control", type = "ashr")
}, error = function(e1) {
  message("ashr failed, trying apeglm...")
  tryCatch({
    lfcShrink(dds, coef = "diet_MCD_vs_Control", type = "apeglm")
  }, error = function(e2) {
    message("Shrinkage not available, using unshrunk results")
    res
  })
})

# Load gene annotations from GTF
gtf_path <- file.path(root, "..", "reference", "raw", "gencode.vM33.annotation.gtf")
gene_annot <- NULL

if (!file.exists(gtf_path)) {
  # Try alternative locations
  alt_paths <- c(
    file.path(root, "..", "InHouse_MCD", "reference", "raw", "gencode.vM33.annotation.gtf"),
  file.path(root, "..", "Public_MCD", "GSE156918", "reference", "raw", "gencode.vM33.annotation.gtf")
  )
  for (p in alt_paths) {
    if (file.exists(p)) {
      gtf_path <- p
      break
    }
  }
}

if (file.exists(gtf_path)) {
  message("Loading gene annotations from: ", gtf_path)
  gtf_data <- read_tsv(gtf_path, comment = "#", col_names = FALSE, col_types = cols())
  gene_annot <- gtf_data %>%
    filter(X3 == "gene") %>%
    mutate(
      gene_id = str_match(X9, 'gene_id "([^"]+)"')[, 2],
      gene_name = str_match(X9, 'gene_name "([^"]+)"')[, 2],
      gene_type = str_match(X9, 'gene_type "([^"]+)"')[, 2]
    ) %>%
    select(gene = gene_id, gene_name, gene_type) %>%
    distinct()
  message("Loaded annotations for ", nrow(gene_annot), " genes")
} else {
  message("GTF not found, skipping gene annotation")
}

# Create results table
res_tbl <- as.data.frame(res_shrunk) %>%
  rownames_to_column("gene") %>%
  arrange(padj)

if (!is.null(gene_annot)) {
  res_tbl <- res_tbl %>%
    left_join(gene_annot, by = "gene") %>%
    select(gene, gene_name, gene_type, everything())
}

# Write results
write_tsv(res_tbl, file.path(out_dir, "deseq2_mcd_vs_control.tsv"))
message("\nWrote DESeq2 results to: ", file.path(out_dir, "deseq2_mcd_vs_control.tsv"))

# Write normalized counts
norm_counts <- counts(dds, normalized = TRUE) %>%
  as.data.frame() %>%
  rownames_to_column("gene")
write_csv(norm_counts, file.path(out_dir, "normalized_counts.csv"))
message("Wrote normalized counts to: ", file.path(out_dir, "normalized_counts.csv"))

# Write VST-transformed counts for visualization
vst_counts <- assay(vst(dds, blind = FALSE)) %>%
  as.data.frame() %>%
  rownames_to_column("gene")
write_csv(vst_counts, file.path(out_dir, "vst_counts.csv"))
message("Wrote VST counts to: ", file.path(out_dir, "vst_counts.csv"))

# Write design info
writeLines(
  c(
    paste0("design: ~ batch + diet"),
    paste0("contrast: diet MCD vs Control"),
    paste0("samples: ", nrow(samples_df)),
    paste0("batches: ", paste(levels(samples_df$batch), collapse = ", ")),
    paste0("genes_tested: ", nrow(res_tbl)),
    paste0("sig_padj_0.1: ", sum(res_tbl$padj < 0.1, na.rm = TRUE)),
    paste0("sig_padj_0.05: ", sum(res_tbl$padj < 0.05, na.rm = TRUE))
  ),
  con = file.path(out_dir, "design.txt")
)

# Summary
message("\n=== Final Summary ===")
message("Significant DEGs (padj < 0.1): ", sum(res_tbl$padj < 0.1, na.rm = TRUE))
message("Significant DEGs (padj < 0.05): ", sum(res_tbl$padj < 0.05, na.rm = TRUE))
message("Up in MCD (padj < 0.1, LFC > 0): ", sum(res_tbl$padj < 0.1 & res_tbl$log2FoldChange > 0, na.rm = TRUE))
message("Down in MCD (padj < 0.1, LFC < 0): ", sum(res_tbl$padj < 0.1 & res_tbl$log2FoldChange < 0, na.rm = TRUE))
