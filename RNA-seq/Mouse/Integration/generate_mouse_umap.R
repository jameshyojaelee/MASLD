
# generate_mouse_umap.R
# Goal: Integrate In-house and Other Diet Mouse RNA-seq datasets and generate UMAP plots.

# Load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(DESeq2)
    library(limma)
    # library(umap) # Package not available
    library(ggplot2)
    library(RColorBrewer)
    library(dplyr)
    library(tibble)
})

# Define Paths
BASE_DIR <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq"
INHOUSE_COUNTS <- file.path(BASE_DIR, "in-house_MCD_RNAseq/counts/featurecounts/gene_counts.txt")
OTHER_COUNTS <- file.path(BASE_DIR, "other_diet_mouse_RNAseq/counts/featurecounts/gene_counts.txt")
INHOUSE_META <- file.path(BASE_DIR, "in-house_MCD_RNAseq/metadata/samples.tsv")
OTHER_META <- file.path(BASE_DIR, "other_diet_mouse_RNAseq/metadata/samples.tsv")
OUT_DIR <- file.path(BASE_DIR, "mouse_integration")

# Ensure output directory exists
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# ==============================================================================
# 1. Load and Clean Counts
# ==============================================================================

message("Loading count matrices...")

# Function to load featureCounts
load_counts <- function(path) {
    # Skip first line if it's a command info line (starts with #)
    # featureCounts often has a command line on line 1, header on line 2
    # But fread is smart. Let's inspect if needed.
    # We'll rely on "Geneid" being the first column name.
    
    # Read header manually to check safely
    header_line <- system(paste("grep -v '^#'", path, "| head -n 1"), intern=TRUE)
    col_names <- unlist(strsplit(header_line, "\t"))
    
    # Read data
    counts <- fread(path, skip = "Geneid") # Skips until it finds Geneid
    
    # Filter for standard chromosomes if needed, but for now just keep all
    # Select only Geneid and sample columns (columns 7 onwards)
    # Actually, keep Geneid, drop Chr, Start, End, Strand, Length
    counts_clean <- counts[, .SD, .SDcols = c("Geneid", colnames(counts)[7:ncol(counts)])]
    
    return(counts_clean)
}

cnts_inhouse <- load_counts(INHOUSE_COUNTS)
cnts_other <- load_counts(OTHER_COUNTS)

# Clean Clean Column Names
# In-house: remove path, keep basename (e.g., 1CF)
clean_cols_inhouse <- colnames(cnts_inhouse)
clean_cols_inhouse[-1] <- sapply(clean_cols_inhouse[-1], function(x) {
    # Assuming path ends with .../SampleID.Aligned...bam or SampleID/...
    # Based on previous `ls`, it was .../1CF.Aligned.sortedByCoord.out.bam
    # We want "1CF"
    # Using gsub to extract the sample ID.
    # Pattern: 1CF.Aligned... -> take ^[^.]+
    basename_x <- basename(x)
    return(sub("\\.Aligned.*", "", basename_x))
})
colnames(cnts_inhouse) <- clean_cols_inhouse

# Other: remove path, keep basename (e.g., SRR...)
clean_cols_other <- colnames(cnts_other)
clean_cols_other[-1] <- sapply(clean_cols_other[-1], function(x) {
    # .../SRR12883362.Aligned...bam
    basename_x <- basename(x)
    return(sub("\\.Aligned.*", "", basename_x))
})
colnames(cnts_other) <- clean_cols_other

message("Counts loaded. In-house samples: ", ncol(cnts_inhouse)-1, " Other samples: ", ncol(cnts_other)-1)

# Merge Counts
counts_merged <- merge(cnts_inhouse, cnts_other, by = "Geneid")
message("Merged dimensions: ", nrow(counts_merged), " genes x ", ncol(counts_merged)-1, " samples")

# Convert to Matrix
mat <- as.matrix(counts_merged[,-1])
rownames(mat) <- counts_merged$Geneid

# ==============================================================================
# 2. Load and Prepare Metadata
# ==============================================================================

message("Loading or constructing metadata...")

# Load Metadata Files
meta_inhouse <- fread(INHOUSE_META)
meta_other <- fread(OTHER_META)

# Prepare In-house Metadata
# Need: sample_id, Diet, Dataset
meta_inhouse_clean <- meta_inhouse %>%
    transmute(
        sample_id = as.character(sample_id),
        Condition = diet,
        Dataset = "In-house"
    )

# Prepare Other Metadata
# Need: sample_id, Condition, Dataset
meta_other_clean <- meta_other %>%
    transmute(
        sample_id = as.character(sample_id),
        Condition = condition,
        Dataset = dataset
    )

# Combine Metadata
metadata_combined <- bind_rows(meta_inhouse_clean, meta_other_clean)

# Align Metadata with Count Matrix
common_samples <- intersect(colnames(mat), metadata_combined$sample_id)
message("Number of matching samples: ", length(common_samples))

if (length(common_samples) < ncol(mat)) {
    missing <- setdiff(colnames(mat), metadata_combined$sample_id)
    warning("Missing metadata for samples: ", paste(head(missing), collapse=", "))
    # Subset matrix to common samples
    mat <- mat[, common_samples]
}

# Ensure metadata order matches matrix columns
metadata_combined <- metadata_combined %>%
    filter(sample_id %in% common_samples) %>%
    arrange(match(sample_id, colnames(mat)))

if (!all(metadata_combined$sample_id == colnames(mat))) {
    stop("Metadata and Matrix sample order mismatch!")
}

# ==============================================================================
# 3. Normalization and Batch Correction
# ==============================================================================

message("Normalizing and correcting batch effects...")

# Create DESeq2 dataset for VST
# Round counts to integer
dds <- DESeqDataSetFromMatrix(
    countData = round(mat),
    colData = metadata_combined,
    design = ~ 1
)

# VST Normalization
vst_data <- vst(dds, blind = TRUE)
vst_mat <- assay(vst_data)

# Remove Batch Effects using Limma
# We want to remove the 'Dataset' effect but preserve 'Condition' if possible?
# Actually, if we just remove Dataset effect, we might remove biological signal if Condition is perfectly confounded with Dataset.
# But here, we have Controls in most datasets.
# So we should model Condition (if possible) to preserve it, or just regress out Dataset.
# Given the heterogeneity, removing Dataset effect is key.
# We will use design matrix for Condition to preserve it.
# Check if Condition is evaluable (e.g. not just 1 level per batch).
# In-house: Control, MCD.
# Other: Control, LIDPAD, etc.
# So we have Controls across datasets. This is good.

# Create design matrix for Condition (to preserve it)
# Handle potential problematic chars in Condition
metadata_combined$Condition <- make.names(metadata_combined$Condition)
design <- model.matrix(~ Condition, data = metadata_combined)

# Remove batch effect
corrected_mat <- removeBatchEffect(vst_mat, batch = metadata_combined$Dataset, design = design)

# ==============================================================================
# 4. PCA and Plotting
# ==============================================================================

message("Running PCA...")

# Run PCA
# Use prcomp on transposed matrix (samples as rows)
pca_result <- prcomp(t(corrected_mat))
pca_coords <- as.data.frame(pca_result$x)
# Calculate variance explained
var_explained <- round(pca_result$sdev^2 / sum(pca_result$sdev^2) * 100, 1)

# Bind with metadata
plot_df <- cbind(pca_coords, metadata_combined)

# Plotting
message("Generating Plot...")

p <- ggplot(plot_df, aes(x = PC1, y = PC2, color = Condition, shape = Dataset)) +
    geom_point(alpha = 0.8, size = 3) +
    theme_bw() +
    labs(
        title = "Integrated Mouse RNA-seq PCA",
        subtitle = "Batch Corrected (Dataset Effect Removed)",
        x = paste0("PC1 (", var_explained[1], "%)"),
        y = paste0("PC2 (", var_explained[2], "%)")
    ) +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "right"
    ) +
    scale_shape_manual(values = 1:length(unique(plot_df$Dataset)))

# Save Plot
out_file <- file.path(OUT_DIR, "mouse_integrated_pca_diet.pdf")
ggsave(out_file, plot = p, width = 10, height = 8, device = "pdf")

message("Done! Saved plot to: ", out_file)
