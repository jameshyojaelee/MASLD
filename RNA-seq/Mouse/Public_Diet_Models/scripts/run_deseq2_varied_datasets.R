
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(stringr)
library(pheatmap)

# Configuration
ROOT_DIR <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/other_diet_mouse_RNAseq"
COUNTS_FILE <- file.path(ROOT_DIR, "counts", "featurecounts", "gene_counts.txt") # Actual location after pipeline
METADATA_FILE <- file.path(ROOT_DIR, "metadata", "samples.tsv")
OUTPUT_DIR <- file.path(ROOT_DIR, "analysis", "deseq2")

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# --- Contrast Definitions ---
# List of comparisons to run. 
# Logic: Dataset -> list(Case, Control, Design)
CONTRASTS <- list(
  # GSE162876: Timepoint matched
  list(
    dataset = "GSE162876",
    id = "20wks_FPC_vs_LFD",
    case = "20wks_FPC",
    control = "20wks_LFD_Ctrl",
    design = "~ diet"
  ),
  list(
    dataset = "GSE162876",
    id = "7wks_CDAHFD_vs_LFD",
    case = "7wks_CDAHFD",
    control = "7wks_LFD_Ctrl",
    design = "~ diet"
  ),
  
  # GSE224069: Standard HFD vs Chow (Vehicle groups)
  list(
    dataset = "GSE224069",
    id = "HFD_vs_Chow",
    case = "High_Fat_Diet_Vehicle",
    control = "Chow_Diet_Vehicle",
    design = "~ diet"
  ),
  
  # GSE274914: 52w (Females Only)
  list(
    dataset = "GSE274914",
    id = "52w_HFD_vs_LFD",
    case = "52w_HFD",
    control = "52w_LFD",
    design = "~ diet",
    subset_sex = "female" # Critical: Control has no males
  ),
  
  # GSE274914: 7w (Mixed Sex)
  list(
    dataset = "GSE274914",
    id = "7w_HFD_vs_LFD",
    case = "7w_HFD",
    control = "7w_LFD",
    design = "~ sex + diet"
  ),

  # GSE159911: LIDPAD (N=152)
  list(
    dataset = "GSE159911",
    id = "LIDPAD_vs_Control",
    case = "LIDPAD",
    control = "Control",
    design = "~ diet"
  )
  # Removed: GSE225616 (siRNA vs Vehicle within GAN diet - not diet vs control)
  # Removed: GSE263273 (OCA, INT787 - pharmacological interventions)
)

# --- Helper Functions ---

load_data <- function() {
  # 1. Load Metadata
  colData <- read.delim(METADATA_FILE, header = TRUE, stringsAsFactors = FALSE)
  rownames(colData) <- colData$sample_id
  
  # 2. Load Counts
  # Use check because file might not exist yet if pipeline running
  if (!file.exists(COUNTS_FILE)) {
    stop(paste("Counts file not found at:", COUNTS_FILE))
  }
  
  # featureCounts format: Geneid, Chr, Start, End, Strand, Length, [Samples...]
  counts <- read.table(COUNTS_FILE, header = TRUE, row.names = 1, comment.char = "#", check.names = FALSE)
  
  # Clean sample names in counts (remove path/extensions)
  # featureCounts uses full BAM paths like: /path/SRR12883362/SRR12883362.Aligned.sortedByCoord.out.bam
  # Extract SRR IDs using regex
  colnames(counts) <- gsub(".*/(SRR[0-9]+)/.*", "\\1", colnames(counts))
  # Keep annotation columns as-is (they don't match the pattern)
  

  # Filter metadata to intersect with counts
  common_samples <- intersect(rownames(colData), colnames(counts))
  
  if (length(common_samples) == 0) {
    stop("No common samples between counts and metadata!")
  }
  
  counts <- counts[, common_samples]
  colData <- colData[common_samples, ]
  
  # Ensure order matches
  counts <- counts[, rownames(colData)]
  
  return(list(counts = counts, colData = colData))
}

run_contrast <- function(config, all_counts, all_colData) {
  message(paste("Processing:", config$id))
  
  # 1. Subset Data
  # Filter by Dataset
  keep_samples <- all_colData$dataset == config$dataset
  
  # Filter by Condition (Case + Control only? Or keep all for variance?)
  # Best practice: Keep only relevant groups to avoid dispersion inflation from unrelated groups
  keep_samples <- keep_samples & (all_colData$condition %in% c(config$case, config$control))
  
  # Filter by Sex (if requested)
  if (!is.null(config$subset_sex)) {
    keep_samples <- keep_samples & (all_colData$sex == config$subset_sex)
  }
  
  sub_colData <- all_colData[keep_samples, ]
  sub_counts <- all_counts[, keep_samples]
  
  # 2. Factor Levels
  # Rename 'condition' to 'diet' for formula consistency standard
  # Or just use the 'condition' column but rename the variable in formula
  sub_colData$diet <- factor(sub_colData$condition, levels = c(config$control, config$case))
  sub_colData$sex <- factor(sub_colData$sex)
  
  # Check if we have enough samples
  if (nrow(sub_colData) < 2) {
    warning(paste("Skipping", config$id, "- Not enough samples."))
    return(NULL)
  }
  
  # 3. Create DESeq Object
  formula <- as.formula(config$design)
  dds <- DESeqDataSetFromMatrix(countData = sub_counts,
                                colData = sub_colData,
                                design = formula)
  
  # 4. Filter Low Counts
  keep <- rowSums(counts(dds) >= 10) >= 3
  dds <- dds[keep, ]
  
  # 5. Run DESeq
  dds <- DESeq(dds)
  
  # 6. Results
  res <- results(dds, contrast = c("diet", config$case, config$control))
  
  # Shrinkage (try apeglm, fallback to normal if not available)
  # Need to check coef name first
  coef_name <- paste0("diet_", config$case, "_vs_", config$control)
  
  if (!coef_name %in% resultsNames(dds)) {
    warning(paste("Expected coefficient", coef_name, "not found. Available:", paste(resultsNames(dds), collapse=", ")))
    # Fallback: use unshrunken results
    resLFC <- res
  } else {
    # Try apeglm, fallback to normal if unavailable
    resLFC <- tryCatch({
      lfcShrink(dds, coef = coef_name, type = "apeglm")
    }, error = function(e) {
      message(paste("apeglm not available, using 'normal' shrinkage:", e$message))
      lfcShrink(dds, coef = coef_name, type = "normal")
    })
  }
  
  # 7. Export
  out_prefix <- file.path(OUTPUT_DIR, paste0(config$dataset, "_", config$id))
  
  # CSV
  res_df <- as.data.frame(resLFC)
  res_df$gene_id <- rownames(res_df)
  write.csv(res_df, paste0(out_prefix, "_results.csv"), row.names = FALSE)
  
  # Volcano Plot
  pdf(paste0(out_prefix, "_volcano.pdf"), width = 6, height = 5)
  
  # Create simple volcano dataframe
  volc_data <- res_df %>% 
    mutate(sig = padj < 0.05 & abs(log2FoldChange) > 1) %>%
    mutate(label = ifelse(sig & rank(padj) < 20, gene_id, NA))
  
  p <- ggplot(volc_data, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
    geom_point(alpha = 0.5) +
    scale_color_manual(values = c("grey", "red")) +
    geom_text_repel(aes(label = label), max.overlaps = 20) +
    theme_minimal() +
    labs(title = paste(config$dataset, config$id),
         subtitle = paste("Case:", config$case, "| Control:", config$control))
  print(p)
  dev.off()
  
  # PCA Plot
  vsd <- vst(dds, blind = FALSE)
  pcaData <- plotPCA(vsd, intgroup = c("diet", "sex"), returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  
  pdf(paste0(out_prefix, "_pca.pdf"), width = 6, height = 5)
  p_pca <- ggplot(pcaData, aes(PC1, PC2, color = diet, shape = sex)) +
    geom_point(size = 3) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
    coord_fixed() +
    theme_bw() +
    labs(title = paste(config$dataset, config$id))
  print(p_pca)
  dev.off()
  
  message(paste("Finished:", config$id))
}

# --- Main Execution ---

main <- function() {
  data <- load_data()
  
  for (contrast in CONTRASTS) {
    tryCatch({
      run_contrast(contrast, data$counts, data$colData)
    }, error = function(e) {
      message(paste("Error processing", contrast$id, ":", e$message))
    })
  }
}

main()
