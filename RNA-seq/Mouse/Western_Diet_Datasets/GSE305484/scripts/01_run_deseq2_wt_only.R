#!/usr/bin/env Rscript
# GSE305484: DESeq2 analysis — WT mice only (WD vs Chow)
# Skytthe et al. (PMID 41234701): CD163KO mouse liver, Western diet, 8 weeks
# Design: 15 WT WD-fed vs 5 WT chow-fed (all male, C57BL/6NTac)
#
# Raw counts already aligned to GRCm39 by depositors (STAR 2.7.8a + featureCounts 2.0).
# Gene IDs are Ensembl (GENCODE vM38-era). We map to gene symbols via GENCODE vM38 GTF.

library(DESeq2)
library(dplyr)
library(readr)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)

# ------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------
BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Mouse/Western_Diet_Datasets/GSE305484"
COUNTS_FILE  <- file.path(BASE, "counts", "GSE305484_GEO_CD163KO_rawcounts.txt.gz")
META_FILE    <- file.path(BASE, "metadata", "sample_metadata.csv")
RESULTS_DIR  <- file.path(BASE, "results")
## Use vM37 GTF (also GRCm39-based); version-stripped Ensembl IDs match depositor's vM38-era IDs
GENCODE_GTF  <- "/gpfs/commons/home/jameslee/reference_genome/refdata-gex-GRCm39-2024-A/annotation/gencode.vM37.chr_patch_hapl_scaff.annotation.gtf.gz"

dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------
# 1. Load counts and metadata
# ------------------------------------------------------------------
cat("Loading counts...\n")
counts_raw <- read.delim(gzfile(COUNTS_FILE), row.names = 1, check.names = FALSE)
cat(sprintf("  Count matrix: %d genes x %d samples\n", nrow(counts_raw), ncol(counts_raw)))

meta <- read_csv(META_FILE, show_col_types = FALSE) %>%
  mutate(sample_id = as.character(sample_id))
cat(sprintf("  Metadata: %d samples\n", nrow(meta)))

# ------------------------------------------------------------------
# 2. Subset to WT only
# ------------------------------------------------------------------
wt_meta <- meta %>% filter(genotype == "CD163WT")
cat(sprintf("  WT samples: %d (Chow=%d, WD=%d)\n",
            nrow(wt_meta),
            sum(wt_meta$condition == "Chow"),
            sum(wt_meta$condition == "WD")))

# Match column names (library names) to WT samples
wt_libs <- wt_meta$library_name
stopifnot(all(wt_libs %in% colnames(counts_raw)))
counts_wt <- counts_raw[, wt_libs]

# ------------------------------------------------------------------
# 3. Build gene annotation from GENCODE vM38
# ------------------------------------------------------------------
cat("Building gene annotation from GENCODE vM38...\n")
if (file.exists(GENCODE_GTF)) {
  gtf_lines <- readLines(gzfile(GENCODE_GTF))
  gene_lines <- gtf_lines[grepl("^\tchr.*\tgene\t", paste0("\t", gtf_lines))]
  # Actually parse properly
  gene_lines <- gtf_lines[grepl("\tgene\t", gtf_lines) & !grepl("^#", gtf_lines)]

  parse_attr <- function(line, key) {
    m <- regmatches(line, regexpr(paste0(key, ' "([^"]+)"'), line))
    if (length(m) == 0 || m == "") return(NA_character_)
    gsub(paste0(key, ' "'), "", gsub('"$', "", m))
  }

  gene_ann <- data.frame(
    ensembl_gene_id = sapply(gene_lines, parse_attr, key = "gene_id", USE.NAMES = FALSE),
    gene_symbol = sapply(gene_lines, parse_attr, key = "gene_name", USE.NAMES = FALSE),
    gene_biotype = sapply(gene_lines, parse_attr, key = "gene_type", USE.NAMES = FALSE),
    stringsAsFactors = FALSE
  )
  # Strip version suffix from Ensembl IDs
  gene_ann$ensembl_gene_id_noversion <- gsub("\\.\\d+$", "", gene_ann$ensembl_gene_id)
  gene_ann <- gene_ann[!duplicated(gene_ann$ensembl_gene_id_noversion), ]
  rownames(gene_ann) <- gene_ann$ensembl_gene_id_noversion
  cat(sprintf("  Annotation: %d genes\n", nrow(gene_ann)))
} else {
  cat("  WARNING: GENCODE vM38 GTF not found. Using Ensembl IDs only.\n")
  gene_ann <- data.frame(
    ensembl_gene_id_noversion = rownames(counts_wt),
    gene_symbol = rownames(counts_wt),
    gene_biotype = NA_character_,
    stringsAsFactors = FALSE
  )
  rownames(gene_ann) <- gene_ann$ensembl_gene_id_noversion
}

# ------------------------------------------------------------------
# 4. Pre-filter: remove genes with < 10 counts total
# ------------------------------------------------------------------
keep <- rowSums(counts_wt) >= 10
counts_filt <- counts_wt[keep, ]
cat(sprintf("  After filtering (rowSum>=10): %d genes\n", nrow(counts_filt)))

# ------------------------------------------------------------------
# 5. DESeq2
# ------------------------------------------------------------------
cat("Running DESeq2...\n")
col_data <- data.frame(
  row.names = wt_meta$library_name,
  condition = factor(wt_meta$condition, levels = c("Chow", "WD"))
)

dds <- DESeqDataSetFromMatrix(countData = counts_filt,
                               colData = col_data,
                               design = ~ condition)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "WD", "Chow"), alpha = 0.05)
res_df <- as.data.frame(res) %>%
  tibble::rownames_to_column("ensembl_gene_id") %>%
  arrange(padj)

# Add gene symbols
res_df <- res_df %>%
  left_join(gene_ann %>% select(ensembl_gene_id_noversion, gene_symbol, gene_biotype),
            by = c("ensembl_gene_id" = "ensembl_gene_id_noversion"))

cat(sprintf("  DESeq2 results: %d genes tested\n", nrow(res_df)))
cat(sprintf("  padj < 0.05: %d\n", sum(res_df$padj < 0.05, na.rm = TRUE)))
cat(sprintf("  padj < 0.05 & |LFC| > 0.5: %d\n",
            sum(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 0.5, na.rm = TRUE)))
cat(sprintf("  padj < 0.05 & |LFC| > 1.0: %d\n",
            sum(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1.0, na.rm = TRUE)))

# Breakdown by direction
sig <- res_df %>% filter(padj < 0.05)
cat(sprintf("  Up in WD: %d (padj<0.05)\n", sum(sig$log2FoldChange > 0, na.rm = TRUE)))
cat(sprintf("  Down in WD: %d (padj<0.05)\n", sum(sig$log2FoldChange < 0, na.rm = TRUE)))

sig_lfc <- res_df %>% filter(padj < 0.05, abs(log2FoldChange) > 0.5)
cat(sprintf("  Up in WD (|LFC|>0.5): %d\n", sum(sig_lfc$log2FoldChange > 0, na.rm = TRUE)))
cat(sprintf("  Down in WD (|LFC|>0.5): %d\n", sum(sig_lfc$log2FoldChange < 0, na.rm = TRUE)))

# ------------------------------------------------------------------
# 6. Write results
# ------------------------------------------------------------------
write_csv(res_df, file.path(RESULTS_DIR, "deseq2_wt_wd_vs_chow.csv"))
cat(sprintf("  Results written to %s\n", file.path(RESULTS_DIR, "deseq2_wt_wd_vs_chow.csv")))

# Also save significant genes only
write_csv(sig, file.path(RESULTS_DIR, "deseq2_wt_wd_vs_chow_sig005.csv"))

# ------------------------------------------------------------------
# 7. Normalized counts (for downstream use)
# ------------------------------------------------------------------
norm_counts <- counts(dds, normalized = TRUE)
write.csv(norm_counts, file.path(RESULTS_DIR, "normalized_counts_wt.csv"))

# ------------------------------------------------------------------
# 8. QC Plots
# ------------------------------------------------------------------
cat("Generating QC plots...\n")

# 8a. PCA
vsd <- vst(dds, blind = TRUE)
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pca_var <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", pca_var[1], "% variance")) +
  ylab(paste0("PC2: ", pca_var[2], "% variance")) +
  ggtitle("GSE305484 — WT mice PCA (WD vs Chow)") +
  scale_color_manual(values = c("Chow" = "#9E9E9E", "WD" = "#D32F2F")) +
  theme_bw(base_size = 14)
ggsave(file.path(RESULTS_DIR, "pca_wt.pdf"), p_pca, width = 7, height = 5)

# 8b. Volcano
res_vol <- res_df %>% filter(!is.na(padj))
res_vol$significance <- case_when(
  res_vol$padj < 0.05 & abs(res_vol$log2FoldChange) > 1 ~ "Sig (|LFC|>1)",
  res_vol$padj < 0.05 & abs(res_vol$log2FoldChange) > 0.5 ~ "Sig (|LFC|>0.5)",
  res_vol$padj < 0.05 ~ "Sig (padj<0.05)",
  TRUE ~ "NS"
)
# Top genes to label
top_up <- res_vol %>% filter(padj < 0.05, log2FoldChange > 0) %>% slice_min(padj, n = 10)
top_dn <- res_vol %>% filter(padj < 0.05, log2FoldChange < 0) %>% slice_min(padj, n = 10)
top_genes <- bind_rows(top_up, top_dn)

p_vol <- ggplot(res_vol, aes(log2FoldChange, -log10(padj), color = significance)) +
  geom_point(alpha = 0.5, size = 0.8) +
  geom_text_repel(data = top_genes, aes(label = gene_symbol), size = 3, max.overlaps = 20, color = "black") +
  scale_color_manual(values = c("NS" = "grey70",
                                "Sig (padj<0.05)" = "steelblue",
                                "Sig (|LFC|>0.5)" = "darkorange",
                                "Sig (|LFC|>1)" = "firebrick")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey40") +
  ggtitle("GSE305484 — WT WD vs Chow Volcano") +
  xlab("log2 Fold Change (WD / Chow)") +
  ylab("-log10(padj)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")
ggsave(file.path(RESULTS_DIR, "volcano_wt.pdf"), p_vol, width = 8, height = 6)

# 8c. Sample distance heatmap
sample_dists <- dist(t(assay(vsd)))
dist_mat <- as.matrix(sample_dists)
rownames(dist_mat) <- paste0(col_data$condition, "_", rownames(col_data))
colnames(dist_mat) <- rownames(dist_mat)
ann_col <- data.frame(Condition = col_data$condition, row.names = rownames(dist_mat))
pdf(file.path(RESULTS_DIR, "sample_distance_heatmap_wt.pdf"), width = 8, height = 7)
pheatmap(dist_mat,
         clustering_distance_rows = sample_dists,
         clustering_distance_cols = sample_dists,
         annotation_col = ann_col,
         color = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
         main = "GSE305484 — WT Sample Distance")
dev.off()

# ------------------------------------------------------------------
# 9. Summary stats file
# ------------------------------------------------------------------
summary_lines <- c(
  "GSE305484 DESeq2 Summary — WT mice only (WD vs Chow)",
  paste0("Date: ", Sys.time()),
  "",
  paste0("Total samples: ", ncol(counts_wt)),
  paste0("  Chow (control): ", sum(col_data$condition == "Chow")),
  paste0("  WD (disease): ", sum(col_data$condition == "WD")),
  "",
  paste0("Genes in count matrix: ", nrow(counts_raw)),
  paste0("Genes after filtering (rowSum>=10): ", nrow(counts_filt)),
  paste0("Genes tested by DESeq2: ", nrow(res_df)),
  "",
  paste0("DEGs (padj < 0.05): ", sum(res_df$padj < 0.05, na.rm = TRUE)),
  paste0("  Up in WD: ", sum(sig$log2FoldChange > 0, na.rm = TRUE)),
  paste0("  Down in WD: ", sum(sig$log2FoldChange < 0, na.rm = TRUE)),
  "",
  paste0("DEGs (padj < 0.05, |LFC| > 0.5): ", nrow(sig_lfc)),
  paste0("  Up: ", sum(sig_lfc$log2FoldChange > 0, na.rm = TRUE)),
  paste0("  Down: ", sum(sig_lfc$log2FoldChange < 0, na.rm = TRUE)),
  "",
  paste0("DEGs (padj < 0.05, |LFC| > 1.0): ",
         sum(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1.0, na.rm = TRUE)),
  "",
  "Top 20 DEGs by padj:",
  paste0("  ", head(res_df$gene_symbol, 20), " (LFC=", round(head(res_df$log2FoldChange, 20), 2),
         ", padj=", formatC(head(res_df$padj, 20), format = "e", digits = 2), ")")
)
writeLines(summary_lines, file.path(RESULTS_DIR, "deseq2_summary.txt"))

cat("\nDone. All outputs in:", RESULTS_DIR, "\n")
