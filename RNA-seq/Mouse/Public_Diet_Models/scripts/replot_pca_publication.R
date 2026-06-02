#!/usr/bin/env Rscript
# Publication-Quality PCA Plots for Diet Mouse RNA-seq
# Re-runs VST and plots with Helvetica fonts and Sanjana Lab color scheme

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(grid)
})

# ==============================================================================
# PUBLICATION THEME CONFIGURATION
# ==============================================================================

# Font Configuration for Illustrator (Type 42)
pdf.options(useDingbats = FALSE)

# Source Color Themes
source("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/publication_color_themes.R")

# Defined in publication_color_themes.R:
# palette1, palette2, pink_gradient, etc.
# Global Lab Colors:
COLOR_CASE <- "#e14b9d"   # Magenta (Primary Highlight)
COLOR_CTRL <- "#4baeef"   # Blue (Non-essential/Background)
COLOR_CASE_ALT <- "#e35070" # Pink
COLOR_CTRL_ALT <- "#253148" # Dark Blue/Grey

# Publication Theme
theme_publication <- function() {
  theme_minimal(base_family = "Helvetica", base_size = 7) +
    theme(
      plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 6),
      legend.text = element_text(size = 6),
      legend.title = element_text(size = 7),
      legend.position = "right", 
      panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
      plot.margin = margin(10, 10, 10, 10)
    )
}

# ==============================================================================
# CONFIGURATION
# ==============================================================================

ROOT_DIR <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/other_diet_mouse_RNAseq"
COUNTS_FILE <- file.path(ROOT_DIR, "counts", "featurecounts", "gene_counts.txt")
METADATA_FILE <- file.path(ROOT_DIR, "metadata", "samples.tsv")
OUTPUT_DIR <- file.path(ROOT_DIR, "analysis", "deseq2")

# Contrast Definitions (Same as run_deseq2_varied_datasets.R)
CONTRASTS <- list(
  list(dataset = "GSE162876", id = "20wks_FPC_vs_LFD", case = "20wks_FPC", control = "20wks_LFD_Ctrl", design = "~ diet"),
  list(dataset = "GSE162876", id = "7wks_CDAHFD_vs_LFD", case = "7wks_CDAHFD", control = "7wks_LFD_Ctrl", design = "~ diet"),
  list(dataset = "GSE224069", id = "HFD_vs_Chow", case = "High_Fat_Diet_Vehicle", control = "Chow_Diet_Vehicle", design = "~ diet"),
  list(dataset = "GSE274914", id = "52w_HFD_vs_LFD", case = "52w_HFD", control = "52w_LFD", design = "~ diet", subset_sex = "female"),
  list(dataset = "GSE274914", id = "7w_HFD_vs_LFD", case = "7w_HFD", control = "7w_LFD", design = "~ sex + diet"),
  list(dataset = "GSE159911", id = "LIDPAD_vs_Control", case = "LIDPAD", control = "Control", design = "~ diet"),
  list(dataset = "GSE225616", id = "GAN_siRNA_vs_Vehicle", case = "GAN_siRNA", control = "GAN_Vehicle", design = "~ condition"),
  list(dataset = "GSE263273", id = "OCA_vs_Vehicle", case = "OCA_0.03%", control = "Vehicle", design = "~ condition"),
  list(dataset = "GSE263273", id = "INT787_vs_Vehicle", case = "INT-787_1.5mg/kg", control = "Vehicle", design = "~ condition")
)

# Nice Titles Map
DATASET_TITLES <- list(
  "GSE162876_20wks_FPC_vs_LFD" = "GSE162876: 20wk FPC vs LFD",
  "GSE162876_7wks_CDAHFD_vs_LFD" = "GSE162876: 7wk CDAHFD vs LFD",
  "GSE224069_HFD_vs_Chow" = "GSE224069: HFD vs Chow",
  "GSE274914_52w_HFD_vs_LFD" = "GSE274914: 52w HFD vs LFD (Female)",
  "GSE274914_7w_HFD_vs_LFD" = "GSE274914: 7w HFD vs LFD",
  "GSE159911_LIDPAD_vs_Control" = "GSE159911: LIDPAD vs Control",
  "GSE225616_GAN_siRNA_vs_Vehicle" = "GSE225616: GAN siRNA vs Vehicle",
  "GSE263273_OCA_vs_Vehicle" = "GSE263273: OCA vs Vehicle",
  "GSE263273_INT787_vs_Vehicle" = "GSE263273: INT787 vs Vehicle"
)

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

load_data <- function() {
  colData <- read.delim(METADATA_FILE, header = TRUE, stringsAsFactors = FALSE)
  rownames(colData) <- colData$sample_id
  
  if (!file.exists(COUNTS_FILE)) stop("Counts file not found")
  
  counts <- read.table(COUNTS_FILE, header = TRUE, row.names = 1, comment.char = "#", check.names = FALSE)
  colnames(counts) <- gsub(".*/(SRR[0-9]+)/.*", "\\1", colnames(counts))
  
  common_samples <- intersect(rownames(colData), colnames(counts))
  if (length(common_samples) == 0) stop("No common samples")
  
  counts <- counts[, common_samples]
  colData <- colData[common_samples, ]
  
  return(list(counts = counts, colData = colData))
}

run_contrast_pca <- function(config, all_counts, all_colData) {
  # Subset Logic
  keep_samples <- all_colData$dataset == config$dataset
  keep_samples <- keep_samples & (all_colData$condition %in% c(config$case, config$control))
  
  if (!is.null(config$subset_sex)) {
    keep_samples <- keep_samples & (all_colData$sex == config$subset_sex)
  }
  
  sub_colData <- all_colData[keep_samples, ]
  sub_counts <- all_counts[, keep_samples]
  
  if (nrow(sub_colData) < 2) return(NULL)
  
  # Factorize
  sub_colData$condition <- factor(sub_colData$condition, levels = c(config$control, config$case))
  # Map condition to 'Group' for simpler plotting
  sub_colData$Group <- ifelse(sub_colData$condition == config$case, "Case", "Control")
  sub_colData$Group <- factor(sub_colData$Group, levels = c("Control", "Case"))
  
  # DESeq Object & VST
  # We use the design just for the object creation, but VST is blind=FALSE usually
  formula_str <- config$design
  # Fix: If design uses 'diet' but we have 'condition', we need to ensure colData matches
  # In original script, they renamed 'condition' to 'diet'. Let's do that if needed.
  if (grepl("diet", formula_str) && !"diet" %in% colnames(sub_colData)) {
     sub_colData$diet <- sub_colData$condition
  }
  
  dds <- DESeqDataSetFromMatrix(countData = sub_counts, colData = sub_colData, design = as.formula(formula_str))
  keep <- rowSums(counts(dds) >= 10) >= 3
  dds <- dds[keep, ]
  
  vsd <- vst(dds, blind = FALSE)
  
  # PCA Data
  # Determine groups for coloring
  # Use 'condition' (the biological group) for color
  # Use 'sex' for shape if multiple sexes exist
  
  pcaData <- plotPCA(vsd, intgroup = c("condition", "sex"), returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  
  # Improve labels in pcaData
  pcaData$Condition_Label <- ifelse(pcaData$condition == config$case, config$case, config$control)
  # Shorten labels if too long?
  # Maybe just use Case/Control mapping for color consistency
  pcaData$Type <- ifelse(pcaData$condition == config$case, "Case", "Control")
  
  # Title
  plot_title <- DATASET_TITLES[[paste0(config$dataset, "_", config$id)]]
  if (is.null(plot_title)) plot_title <- paste(config$dataset, config$id)
  
  # Plot
  p <- ggplot(pcaData, aes(PC1, PC2, color = condition, shape = sex)) +
    geom_point(size = 3, alpha = 0.8) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
    coord_fixed() +
    theme_publication() +
    labs(title = plot_title, color = "Condition", shape = "Sex") +
    scale_color_manual(values = c(COLOR_CTRL, COLOR_CASE, COLOR_CTRL_ALT, COLOR_CASE_ALT)) # Fallback palette
    
    # Try to be smart about colors:
    # If we have exactly 2 conditions, map them to Control/Case colors
    curr_levels <- levels(pcaData$condition)
    # Identifying which is which
    # We set levels earlier: Control, Case
    # So index 1 = Control, index 2 = Case
    
    p <- p + scale_color_manual(values = c(COLOR_CTRL, COLOR_CASE))
  
  
  # Output
  out_file <- file.path(OUTPUT_DIR, paste0(config$dataset, "_", config$id, "_pca.pdf"))
  
  cairo_pdf(out_file, width = 5, height = 4, family = "Helvetica")
  print(p)
  dev.off()
  
  cat("  Created PCA:", out_file, "\n")
}

# ==============================================================================
# MAIN
# ==============================================================================

cat("========================================\n")
cat("Publication PCA Plot Generator\n")
cat("========================================\n\n")

data <- load_data()

for (contrast in CONTRASTS) {
  tryCatch({
    cat("Processing:", contrast$id, "\n")
    run_contrast_pca(contrast, data$counts, data$colData)
  }, error = function(e) {
    cat("  Error:", e$message, "\n")
  })
}

cat("\nDone!\n")
