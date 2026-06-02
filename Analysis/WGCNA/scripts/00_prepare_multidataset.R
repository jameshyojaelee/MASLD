# WGCNA Pipeline - Step 00: Multi-Dataset Preparation
# This script ingests "patient", "other_MCD", and "in-house_MCD" datasets,
# standardizes them (VST or Log2TPM), and saves them for WGCNA.

library(DESeq2)
library(tidyverse)
library(limma) # for normalizeQuantiles if needed, or just standard use

options(stringsAsFactors = FALSE)

# ==============================================================================
# CONFIGURATION
# ==============================================================================
ROOT_DIR <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq"
OUTPUT_DIR <- file.path(ROOT_DIR, "WGCNA_pipeline/results/01_preprocessing")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# 1. Patient Data Helpers
PATIENT_HELPER <- file.path(ROOT_DIR, "patient_RNAseq/analysis/downstream/_shared_cross_dataset_helpers.R")

# 2. Other MCD (FeatureCounts output)
OTHER_MCD_COUNTS <- file.path(ROOT_DIR, "other_MCD_RNAseq/counts/featurecounts/gene_counts.txt")
OTHER_MCD_META <- file.path(ROOT_DIR, "other_MCD_RNAseq/metadata/samples.tsv")

# 3. In-House MCD (Normalized TPMs)
INHOUSE_MCD_FILE <- file.path(ROOT_DIR, "in-house_MCD_RNAseq/normalized_counts_all_samples.csv")
INHOUSE_MCD_META <- file.path(ROOT_DIR, "in-house_MCD_RNAseq/metadata/samples.tsv") # Assuming metadata exists here

# ==============================================================================
# DATASET 1: PATIENT (GSE130970 + GSE135251) - VST
# ==============================================================================
message("--- Processing Patient Data ---")
source(PATIENT_HELPER)
# root argument for combine_counts_and_metadata needs to point such that default paths work
# Helper expects to find "results" and "data" relative to project root.
# We pass the absolute paths directly to be safe.
patient_counts_gse130970 <- file.path(ROOT_DIR, "patient_RNAseq/results/GSE130970/counts/gene_counts_matrix.txt")
patient_counts_gse135251 <- file.path(ROOT_DIR, "patient_RNAseq/results/GSE135251/counts/gene_counts_matrix.txt")
patient_meta_gse130970 <- file.path(ROOT_DIR, "patient_RNAseq/data/metadata/GSE130970_SraRunTable.csv")
patient_meta_gse135251 <- file.path(ROOT_DIR, "patient_RNAseq/data/metadata/GSE135251_SraRunTable.csv")

patient_combined <- combine_counts_and_metadata(
  counts_paths = c(GSE130970 = patient_counts_gse130970, GSE135251 = patient_counts_gse135251),
  metadata_paths = c(GSE130970 = patient_meta_gse130970, GSE135251 = patient_meta_gse135251),
  dataset_ids = c("GSE130970", "GSE135251")
)

cnts <- patient_combined$count_matrix
meta <- patient_combined$metadata
rownames(cnts) <- cnts$GeneID
cnts$GeneID <- NULL
cnts <- as.matrix(cnts)

# Prepare DESeq2 object for VST
# Simplified design: ~ dataset + condition
# Using 'lenient_control' to define condition
meta$condition <- ifelse(meta$lenient_control, "Control", "MASLD")
dds <- DESeqDataSetFromMatrix(countData = cnts, colData = meta, design = ~ dataset + condition)

# Filter low counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# VST
message("Running VST on Patient data...")
vst_patient <- vst(dds, blind = FALSE)
datExpr_patient <- t(assay(vst_patient))
metadata_patient <- colData(vst_patient) %>% as.data.frame()

# Save
save(datExpr_patient, metadata_patient, file = file.path(OUTPUT_DIR, "Patient.RData"))


# ==============================================================================
# DATASET 2: OTHER MCD (Mouse) - VST
# ==============================================================================
message("--- Processing Other MCD Data ---")
# FeatureCounts format: Geneid, Chr, Start, End, Strand, Length, [Samples...]
other_counts <- read.table(OTHER_MCD_COUNTS, header = TRUE, comment.char = "#", check.names = FALSE)
rownames(other_counts) <- other_counts$Geneid
# Remove annotation columns
sample_cols <- colnames(other_counts)[7:ncol(other_counts)]
other_mat <- as.matrix(other_counts[, sample_cols])

# Load Metadata
# Inspect metadata format. Assuming tab separated
other_meta_raw <- read.table(OTHER_MCD_META, header = TRUE, sep = "\t")
# Verify sample matching
common_samples <- intersect(colnames(other_mat), other_meta_raw$sample) # Adjust 'sample' column name if needed

if (length(common_samples) == 0) {
  # Fallback: Extract GSM IDs from long paths
  # Column names look like: .../GSM4748152/GSM4748152.Aligned.sortedByCoord.out.bam
  # We want "GSM4748152"
  clean_cols <- gsub(".*(GSM[0-9]+).*", "\\1", colnames(other_mat))
  colnames(other_mat) <- clean_cols
  
  common_samples <- intersect(colnames(other_mat), other_meta_raw$sample_id)
}

if (length(common_samples) == 0) {
  # Debugging info
  cat("Column names (head): ", head(colnames(other_mat)), "\n")
  cat("Metadata sample IDs (head): ", head(other_meta_raw$sample_id), "\n")
  stop("No matching samples for Other MCD data after cleanup")
}

other_mat <- other_mat[, common_samples]
other_meta <- other_meta_raw[match(common_samples, other_meta_raw$sample_id), ]
rownames(other_meta) <- other_meta$sample_id

# DESeq2 for VST
dds_other <- DESeqDataSetFromMatrix(countData = other_mat, colData = other_meta, design = ~ 1)
keep <- rowSums(counts(dds_other)) >= 10
dds_other <- dds_other[keep, ]

message("Running VST on Other MCD data...")
vst_other <- vst(dds_other, blind = TRUE)
datExpr_other <- t(assay(vst_other))
metadata_other <- colData(vst_other) %>% as.data.frame()

save(datExpr_other, metadata_other, file = file.path(OUTPUT_DIR, "Other_MCD.RData"))


# ==============================================================================
# DATASET 3: IN-HOUSE MCD - Log2(TPM+1)
# ==============================================================================
message("--- Processing In-House MCD Data ---")
# CSV: gene, gene_name, gene_type, [Samples...]
inhouse_raw <- read.csv(INHOUSE_MCD_FILE, check.names = FALSE)
# Filter for unique genes if duplicates exist (though header check didn't show obvious ones)
# Just in case, aggregate or unique
inhouse_raw <- inhouse_raw[!duplicated(inhouse_raw$gene), ]
rownames(inhouse_raw) <- inhouse_raw$gene

# Extract numeric data
sample_cols_inhouse <- grep("^[1-3][CM][FM]", colnames(inhouse_raw), value = TRUE)
inhouse_tpm <- as.matrix(inhouse_raw[, sample_cols_inhouse])

# Log Transformation
datExpr_inhouse <- t(log2(inhouse_tpm + 1))

# Construct metadata from sample names
# Format: [Replicate][Condition][Sex] e.g., 1CF = Rep1, Control, Female
metadata_inhouse <- data.frame(sample_id = colnames(datExpr_inhouse))
metadata_inhouse$replicate <- substr(metadata_inhouse$sample_id, 1, 1)
metadata_inhouse$condition_code <- substr(metadata_inhouse$sample_id, 2, 2) # C or M
metadata_inhouse$sex <- substr(metadata_inhouse$sample_id, 3, 3) # F or M
metadata_inhouse$condition <- ifelse(metadata_inhouse$condition_code == "C", "Control", "MCD")
rownames(metadata_inhouse) <- metadata_inhouse$sample_id

save(datExpr_inhouse, metadata_inhouse, file = file.path(OUTPUT_DIR, "InHouse_MCD.RData"))

message("All datasets processed and saved to ", OUTPUT_DIR)
