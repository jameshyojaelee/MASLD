
# 01_per_study_plots.R
# Generates core visualizations (PCA, Volcano, Heatmap) for all Human Patient Cohorts
# Updated: Feb 2026 (Sanjana Lab Theme)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(DESeq2)
  library(ComplexHeatmap)
  library(circlize)
  library(EnhancedVolcano)
  library(grid) # For gpar
})

# Paths
BASE_DIR <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts"
INT_DIR <- file.path(BASE_DIR, "analysis/integration")
RESULTS_DIR <- file.path(BASE_DIR, "analysis/visualization/plots/per_study")
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

# Source plotting functions & theme
FUNC_DIR <- file.path(BASE_DIR, "analysis/visualization/functions")
source(file.path(FUNC_DIR, "theme_publication.R"))
source(file.path(FUNC_DIR, "plot_volcano.R"))
source(file.path(FUNC_DIR, "plot_pca.R"))
source(file.path(FUNC_DIR, "plot_heatmap.R"))

# Load Integrated Data
message("Loading integrated counts and metadata...")
counts_path <- file.path(INT_DIR, "results/integration/merged_counts_raw.rds")
meta_path <- file.path(INT_DIR, "results/integration/meta_matched.rds")

if (!file.exists(counts_path) || !file.exists(meta_path)) {
  stop("Integrated counts/metadata not found. Run integration pipeline first.")
}

counts_all <- readRDS(counts_path)
meta_all <- readRDS(meta_path)

# Datasets to process (only those present in metadata)
DATASETS <- unique(meta_all$dataset)
message(paste("Found datasets:", paste(DATASETS, collapse=", ")))

for (ds in DATASETS) {
  message(paste0("\nProcessing ", ds, "..."))
  
  # create output dir
  out_dir <- file.path(RESULTS_DIR, ds)
  dir.create(out_dir, showWarnings = FALSE)

  # 1. Load DE Results
  de_file <- file.path(INT_DIR, "results/per_study", paste0(ds, "_de_results.csv"))
  if (!file.exists(de_file)) {
    message(paste0("  DE results not found for ", ds, ". Skipping."))
    next
  }
  
  res <- fread(de_file)
  setDF(res) # Ensure data.frame
  
  # Load mapping if not already loaded
  if (!exists("id_map")) {
      map_file <- file.path(INT_DIR, "results/integration/human_mouse_ortholog_comparison.csv")
      if (file.exists(map_file)) {
          map_df <- fread(map_file)
          id_map <- map_df[, c("gene_base", "human_symbol")]
          id_map <- id_map[!duplicated(id_map$gene_base), ]
      } else {
          id_map <- NULL
      }
  }
  
  # Map IDs to Symbols
  # Current 'gene' column is Ensembl ID (e.g., ENSG00000225880.7). 
  # Need to strip version for mapping if map has no version.
  # Map file has ENSG00000... (no version? Check first row: ENSG00000198830).
  # Yes, likely no version.
  
  if ("gene" %in% names(res) && ! "symbol" %in% names(res)) {
      res$gene_base <- gsub("\\..*", "", res$gene)
      if (!is.null(id_map)) {
          res <- merge(res, id_map, by="gene_base", all.x=TRUE)
          # If symbol is NA, fallback to gene ID
          res$symbol <- ifelse(is.na(res$human_symbol) | res$human_symbol == "", res$gene, res$human_symbol)
      } else {
          res$symbol <- res$gene # Fallback if map missing
      }
  }
  
  # 2. Subset Counts & Metadata
  meta_sub <- meta_all[meta_all$dataset == ds, ]
  # Align counts
  samples <- intersect(colnames(counts_all), meta_sub$sample_id)
  if (length(samples) < 3) {
    message("  Too few samples. Skipping.")
    next
  }
  
  counts_sub <- counts_all[, samples]
  meta_sub <- as.data.frame(meta_sub[match(samples, meta_sub$sample_id), ])
  
  # 3. Normalize (VST) for PCA/Heatmap
  # Create DESeqDataSet just for VST (design doesn't matter much for blind VST)
  if (length(unique(meta_sub$condition)) < 2) {
    dds <- DESeqDataSetFromMatrix(countData = counts_sub, colData = meta_sub, design = ~ 1)
  } else {
    dds <- DESeqDataSetFromMatrix(countData = counts_sub, colData = meta_sub, design = ~ condition)
  }
  
  # Filter low count genes
  keep <- rowSums(counts(dds) >= 10) >= 3
  dds <- dds[keep, ]
  
  message("  Running VST...")
  vst_obj <- vst(dds, blind = TRUE)
  vst_mat <- assay(vst_obj)
  
  # 4. Generate Plots using helper functions (which now use theme_publication)
  
  # Volcano
  message("  Plotting Volcano...")
  # Map columns
  res_plot <- as.data.frame(res)
  if ("adj.P.Val" %in% names(res_plot)) names(res_plot)[names(res_plot) == "adj.P.Val"] <- "padj"
  if ("logFC" %in% names(res_plot)) names(res_plot)[names(res_plot) == "logFC"] <- "log2FoldChange"
  if (!"symbol" %in% names(res_plot) && "gene" %in% names(res_plot)) res_plot$symbol <- res_plot$gene
  
  p_vol <- plot_volcano(res_plot, 
                        title=paste0(ds, ": Disease vs Control"), 
                        subtitle = paste0("FDR < 0.1, |LFC| > 1"),
                        fc_cutoff=1, p_cutoff=0.1) 
  save_pdf(p_vol, file.path(out_dir, "volcano.pdf"), width=5, height=5) # 5x5 inches fits Guidelines

  # PCA
  message("  Plotting PCA...")
  p_pca <- plot_pca(vst_mat, meta_sub, color_by="condition", title=paste0(ds, " PCA"))
  save_pdf(p_pca, file.path(out_dir, "pca.pdf"), width=5, height=4) # 4:3 ratio mostly
  
  # Heatmap (Top 50 DEGs)
  message("  Plotting Heatmap...")
  if ("padj" %in% names(res_plot)) {
    # Sort and take top 50 symbols
    sorted <- res_plot[order(res_plot$padj), ]
    top_genes <- sorted$symbol[1:50] 
    top_genes <- top_genes[!is.na(top_genes)]
    
    rows_to_plot <- intersect(top_genes, rownames(vst_mat))
    
    if (length(rows_to_plot) > 5) {
      heatmap_mat <- vst_mat[rows_to_plot, ]
      
      # For ComplexHeatmap, we create pdf device directly
      # theme_publication.R defines save_pdf for ggplot, but Heatmap is grid-based.
      # We can use pdf() with correct params manually or verify if save_pdf handles general grid objects.
      # ggsave can handle grid objects? Not reliably for ComplexHeatmap.
      # Use manual pdf() with useDingbats = FALSE
      
      pdf(file.path(out_dir, "heatmap_top50.pdf"), width=6, height=8, useDingbats = FALSE)
      draw(plot_heatmap(heatmap_mat, meta_sub, top_annotation_col="condition", title=paste0(ds, " Top 50 DEGs")))
      dev.off()
    }
  }
}

message("Done.")
