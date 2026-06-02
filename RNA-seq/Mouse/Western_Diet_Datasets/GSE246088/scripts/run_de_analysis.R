#!/usr/bin/env Rscript
# GSE246088 — Western diet / HFD Plvap KO study (Hansen et al. Cell Metabolism 2025)
# Whole-liver bulk RNA-seq from C57BL/6JBomTac mice
#
# Design: 2 x 3 factorial (genotype: Plvap_Control vs Plvap_Knockout;
#         diet: Chow vs WD vs HFD), all male, 12h fasted before sacrifice.
#
# Key contrasts for MASLD library design:
#   C1: WD_Control vs Chow_Control — Western diet effect in WT liver
#   C2: HFD_Control vs Chow_Control — HFD effect in WT liver
#   C3: WD_KO vs WD_Control — Plvap KO effect under WD stress
#   C4: HFD_KO vs HFD_Control — Plvap KO effect under HFD stress
#
# CRITICAL NOTE: tissue = whole liver, NOT sorted stellate cells.
#   The "HSC Plvap knockout" in sample titles = Lrat-Cre-driven HSC-specific
#   Plvap deletion in vivo, but RNA was extracted from whole liver.
#   This makes the diet-effect contrasts (C1, C2) fully comparable to our
#   existing mouse diet datasets. The KO contrasts (C3, C4) capture
#   cell-autonomous HSC effects + downstream non-cell-autonomous changes.

suppressPackageStartupMessages({
  library(DESeq2)
  library(data.table)
  library(ggplot2)
})

cat("=== GSE246088 DE Analysis ===\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# ---- Paths ----
base_dir <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Mouse/Western_Diet_Datasets/GSE246088"
wd_hfd_counts_file <- file.path(base_dir, "counts", "GSE246086_6w_WD_HFD_12h_fasted_Plvap_RNA-seq_noadj.txt.gz")
chow_counts_file   <- file.path(base_dir, "counts", "GSE246086_12h_fasted_Plvap_RNA-seq_noadj.txt.gz")
meta_file           <- file.path(base_dir, "metadata", "sample_metadata.csv")
out_dir             <- file.path(base_dir, "results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Load metadata ----
meta <- fread(meta_file)
cat("Metadata loaded:", nrow(meta), "samples\n")
cat("Condition breakdown:\n")
print(table(meta$condition))

# ---- Load and merge count matrices ----
cat("\nLoading count matrices...\n")

# WD/HFD file (clean column names)
wd_hfd <- fread(wd_hfd_counts_file)
gene_info <- wd_hfd[, .(Geneid, Chr, Start, End, Strand, Length)]
wd_hfd_counts <- as.matrix(wd_hfd[, 7:ncol(wd_hfd), with = FALSE])
rownames(wd_hfd_counts) <- wd_hfd$Geneid
cat("WD/HFD matrix:", nrow(wd_hfd_counts), "genes x", ncol(wd_hfd_counts), "samples\n")

# Chow file (extract sample IDs from BAM paths)
chow <- fread(chow_counts_file)
chow_counts <- as.matrix(chow[, 7:ncol(chow), with = FALSE])
rownames(chow_counts) <- chow$Geneid
# Clean column names from BAM paths -> animal IDs
colnames(chow_counts) <- gsub("^\\./align/(A\\d+)_S\\d+_R.*", "\\1", colnames(chow_counts))
cat("Chow matrix:", nrow(chow_counts), "genes x", ncol(chow_counts), "samples\n")

# Verify gene order matches
stopifnot(identical(rownames(wd_hfd_counts), rownames(chow_counts)))

# Merge matrices
all_counts <- cbind(chow_counts, wd_hfd_counts)
cat("Merged matrix:", nrow(all_counts), "genes x", ncol(all_counts), "samples\n")

# Keep only samples in metadata
keep_samples <- intersect(colnames(all_counts), meta$sample_id)
all_counts <- all_counts[, keep_samples]
meta <- meta[match(keep_samples, sample_id)]
cat("After metadata match:", ncol(all_counts), "samples\n\n")

# ---- Pre-filtering ----
# Keep genes with >= 10 counts in >= 3 samples
keep_genes <- rowSums(all_counts >= 10) >= 3
all_counts <- all_counts[keep_genes, ]
gene_info  <- gene_info[keep_genes]
cat("After pre-filtering:", nrow(all_counts), "genes retained\n")

# ---- DESeq2 full model ----
meta$group <- factor(meta$condition,
                     levels = c("Chow_Control", "Chow_KO",
                                "WD_Control", "WD_KO",
                                "HFD_Control", "HFD_KO"))
dds <- DESeqDataSetFromMatrix(
  countData = all_counts,
  colData   = meta,
  design    = ~ group
)
cat("\nRunning DESeq2...\n")
dds <- DESeq(dds, parallel = FALSE)
cat("DESeq2 done.\n\n")

# ---- Extract contrasts ----
contrasts <- list(
  WD_vs_Chow     = c("group", "WD_Control",  "Chow_Control"),
  HFD_vs_Chow    = c("group", "HFD_Control", "Chow_Control"),
  WD_KO_vs_Ctrl  = c("group", "WD_KO",       "WD_Control"),
  HFD_KO_vs_Ctrl = c("group", "HFD_KO",      "HFD_Control")
)

results_list <- list()
for (nm in names(contrasts)) {
  cat("Contrast:", nm, "\n")
  res <- results(dds, contrast = contrasts[[nm]])
  res_dt <- as.data.table(as.data.frame(res), keep.rownames = "gene_id")
  res_dt <- merge(res_dt, gene_info[, .(Geneid, Length)], by.x = "gene_id", by.y = "Geneid", all.x = TRUE)
  res_dt <- res_dt[order(padj)]
  results_list[[nm]] <- res_dt

  n_up   <- sum(res_dt$padj < 0.05 & res_dt$log2FoldChange > 0.5, na.rm = TRUE)
  n_down <- sum(res_dt$padj < 0.05 & res_dt$log2FoldChange < -0.5, na.rm = TRUE)
  n_sig  <- sum(res_dt$padj < 0.05, na.rm = TRUE)
  cat(sprintf("  padj<0.05: %d total (%d up LFC>0.5, %d down LFC<-0.5)\n", n_sig, n_up, n_down))

  out_file <- file.path(out_dir, paste0("de_results_", nm, ".csv"))
  fwrite(res_dt, out_file)
  cat("  Saved:", out_file, "\n\n")
}

# ---- Summary table ----
cat("=== SUMMARY ===\n")
summary_rows <- lapply(names(results_list), function(nm) {
  dt <- results_list[[nm]]
  data.table(
    contrast       = nm,
    n_tested       = sum(!is.na(dt$padj)),
    n_sig_005      = sum(dt$padj < 0.05, na.rm = TRUE),
    n_up_lfc05     = sum(dt$padj < 0.05 & dt$log2FoldChange > 0.5, na.rm = TRUE),
    n_down_lfc05   = sum(dt$padj < 0.05 & dt$log2FoldChange < -0.5, na.rm = TRUE),
    n_sig_01       = sum(dt$padj < 0.1, na.rm = TRUE),
    median_baseMean = median(dt$baseMean, na.rm = TRUE)
  )
})
summary_dt <- rbindlist(summary_rows)
print(summary_dt)
fwrite(summary_dt, file.path(out_dir, "de_summary.csv"))

# ---- Volcano plots ----
cat("\nGenerating volcano plots...\n")
for (nm in names(results_list)) {
  dt <- results_list[[nm]]
  dt$sig <- ifelse(is.na(dt$padj), "NS",
                   ifelse(dt$padj < 0.05 & abs(dt$log2FoldChange) > 0.5, "Sig", "NS"))

  p <- ggplot(dt, aes(x = log2FoldChange, y = -log10(pvalue), color = sig)) +
    geom_point(size = 0.5, alpha = 0.5) +
    scale_color_manual(values = c(NS = "grey70", Sig = "firebrick3")) +
    labs(title = paste("GSE246088:", gsub("_", " ", nm)),
         subtitle = sprintf("padj<0.05 & |LFC|>0.5: %d up, %d down",
                            sum(dt$sig == "Sig" & dt$log2FoldChange > 0, na.rm = TRUE),
                            sum(dt$sig == "Sig" & dt$log2FoldChange < 0, na.rm = TRUE)),
         x = "log2 Fold Change", y = "-log10(p-value)") +
    theme_classic(base_size = 12) +
    theme(legend.position = "none")

  ggsave(file.path(out_dir, paste0("volcano_", nm, ".pdf")), p, width = 6, height = 5)
}

# ---- Normalized counts for downstream use ----
norm_counts <- counts(dds, normalized = TRUE)
fwrite(as.data.table(norm_counts, keep.rownames = "gene_id"),
       file.path(out_dir, "normalized_counts.csv.gz"))
cat("Normalized counts saved.\n")

# ---- PCA ----
vsd <- vst(dds, blind = TRUE)
pca_data <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
pca_var  <- round(100 * attr(pca_data, "percentVar"), 1)

p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = group, label = name)) +
  geom_point(size = 3) +
  labs(title = "GSE246088 — PCA (VST-transformed)",
       x = paste0("PC1 (", pca_var[1], "%)"),
       y = paste0("PC2 (", pca_var[2], "%)")) +
  theme_classic(base_size = 12) +
  scale_color_brewer(palette = "Set2")
ggsave(file.path(out_dir, "pca_all_samples.pdf"), p_pca, width = 8, height = 6)

cat("\n=== DONE ===\n")
cat("All outputs in:", out_dir, "\n")
sessionInfo()
