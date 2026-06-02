#!/usr/bin/env Rscript
# M01_mouse_sample_qc.R
# ---------------------------------------------------------------------------
# Merge all mouse featureCounts matrices (vM38 GTF) and perform QC:
#   - PCA outlier detection per dataset (> 3 SD from centroid)
#   - Library size filter (< 5M counts)
# Input:  Gene count files from recount (InHouse, GSE156918, GSE205974)
#         + Diet Models featureCounts, unified_mouse_metadata.csv
# Output: merged_counts_raw.rds, meta_matched.rds, sample_qc_report.csv
# ---------------------------------------------------------------------------

# ---- Seed pinning (T2.4, 2026-04-22) -----
set.seed(42)
Sys.setenv(R_PARALLEL_SEED = "42")

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(ggplot2)
})

MOUSE  <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Mouse"
INT    <- file.path(MOUSE, "Unified_Integration")
RDIR   <- file.path(INT, "results")
QCDIR  <- file.path(INT, "qc")
dir.create(RDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(QCDIR, recursive = TRUE, showWarnings = FALSE)

cat("=== M01: Mouse Sample QC ===\n\n")

# --- Load metadata ---
meta <- fread(file.path(INT, "metadata/unified_mouse_metadata.csv"))
cat("Loaded metadata:", nrow(meta), "samples\n")

# --- tximport code path ---
# If a pre-built tximport RDS exists (Kallisto pipeline), load counts + lengths
# from it instead of featureCounts TSVs. The length matrix is saved separately
# so M02/M03 can apply transcript-length offsets to voom/dream.
tximport_rds_path <- file.path(INT, "counts/tximport/txi.rds")
USE_TXIMPORT <- file.exists(tximport_rds_path)

if (USE_TXIMPORT) {
  cat("\n*** tximport RDS detected:", tximport_rds_path, "***\n")
  cat("Loading counts + gene lengths from tximport (Kallisto pipeline)\n")
  txi <- readRDS(tximport_rds_path)

  # Validate required slots
  stopifnot(
    "tximport object must contain $counts"  = !is.null(txi$counts),
    "tximport object must contain $length"   = !is.null(txi$length),
    "counts and length dims must match"      = identical(dim(txi$counts), dim(txi$length))
  )

  # Verify countsFromAbundance was "no" — offset-based analysis requires raw
  # estimated counts; "lengthScaledTPM" bakes lengths in, double-correcting.
  if (!is.null(txi$countsFromAbundance) && txi$countsFromAbundance != "no") {
    warning("tximport countsFromAbundance = '", txi$countsFromAbundance,
            "'; expected 'no' for offset-based DE. Proceeding anyway.")
  }

  merged_txi <- round(txi$counts)  # tximport counts are non-integer estimates
  lengths_txi <- txi$length

  # Strip Ensembl version suffix for cross-tool mergeability
  rownames(merged_txi)  <- sub("\\.[0-9]+$", "", rownames(merged_txi))
  rownames(lengths_txi) <- sub("\\.[0-9]+$", "", rownames(lengths_txi))

  cat("  Genes:", nrow(merged_txi), "  Samples:", ncol(merged_txi), "\n")

  # Save the gene-length matrix for M02/M03 offset correction
  saveRDS(lengths_txi, file.path(RDIR, "merged_gene_lengths.rds"))
  cat("  Saved gene-length matrix:", file.path(RDIR, "merged_gene_lengths.rds"), "\n\n")
}

# --- Helper: load featureCounts, extract sample IDs from BAM paths ---
load_counts <- function(path, label) {
  cat("  Loading:", label, "from", basename(path), "\n")
  
  # Auto-detect whether first line is a comment (starts with #)
  first_line <- readLines(path, n = 1)
  skip_n <- ifelse(grepl("^#", first_line), 1, 0)
  
  raw <- fread(path, skip = skip_n)
  gene_cols <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
  sample_cols <- setdiff(names(raw), gene_cols)
  
  # Extract sample ID from BAM path
  # Handles both: .../star/SAMPLE/SAMPLE.Aligned... and .../star/SAMPLE.Aligned...
  new_names <- basename(sample_cols)
  new_names <- gsub("\\.Aligned\\.sortedByCoord\\.out\\.bam$", "", new_names)
  setnames(raw, sample_cols, new_names)
  
  counts <- as.matrix(raw[, ..new_names])
  # Strip Ensembl version suffix (e.g. ENSMUSG00000102693.2 -> ENSMUSG00000102693)
  # so that datasets quantified with different GTF versions or tools (featureCounts
  # vs Kallisto) produce mergeable gene IDs.  Idempotent: no-op if already bare.
  rownames(counts) <- sub("\\.[0-9]+$", "", raw$Geneid)
  cat("    Genes:", nrow(counts), "  Samples:", ncol(counts), "\n")
  return(counts)
}

# --- Load all count matrices ---
if (USE_TXIMPORT) {
  # tximport path: merged matrix already built above
  cat("\nUsing tximport counts (skipping featureCounts loading)\n")
  merged <- merged_txi
  cat("Matrix:", nrow(merged), "genes x", ncol(merged), "samples\n")
} else {
  cat("\nLoading featureCounts matrices...\n")
  count_list <- list()

  # Re-counted with vM38 GTF
  count_list[["InHouse_MCD"]] <- load_counts(
    file.path(INT, "counts/featurecounts/gene_counts_inhouse.txt"), "InHouse_MCD")
  count_list[["GSE156918"]] <- load_counts(
    file.path(INT, "counts/featurecounts/gene_counts_gse156918.txt"), "GSE156918")
  count_list[["GSE205974"]] <- load_counts(
    file.path(INT, "counts/featurecounts/gene_counts_gse205974.txt"), "GSE205974")

  # Diet Models — already vM38
  count_list[["DietModels"]] <- load_counts(
    file.path(MOUSE, "Public_Diet_Models/counts/featurecounts/gene_counts.txt"), "DietModels")

  # --- Find common genes ---
  common_genes <- Reduce(intersect, lapply(count_list, rownames))
  cat("\nCommon genes across datasets:", length(common_genes), "\n")

  # --- Merge into single matrix ---
  merged <- do.call(cbind, lapply(count_list, function(x) x[common_genes, ]))
  cat("Merged matrix:", nrow(merged), "genes x", ncol(merged), "samples\n")
}

# --- Match metadata ---
shared_samples <- intersect(colnames(merged), meta$sample_id)
cat("Samples in counts AND metadata:", length(shared_samples), "\n")

merged <- merged[, shared_samples]
meta_matched <- meta[match(shared_samples, sample_id)]

# ===== PCA OUTLIER DETECTION =====
cat("\n===== PCA OUTLIER DETECTION =====\n")

# TMM normalize + logCPM for PCA
dge <- DGEList(counts = merged)
keep <- filterByExpr(dge, group = meta_matched$group_binary)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)
logcpm <- cpm(dge, log = TRUE, prior.count = 1)

outlier_flags <- rep(FALSE, ncol(merged))
names(outlier_flags) <- colnames(merged)

for (ds in unique(meta_matched$dataset)) {
  idx <- which(meta_matched$dataset == ds)
  if (length(idx) < 5) {
    cat("  ", ds, ": too few samples (", length(idx), "), skipping PCA\n")
    next
  }
  
  pca <- prcomp(t(logcpm[, idx]), center = TRUE, scale. = FALSE)
  pc12 <- pca$x[, 1:2]
  centroid <- colMeans(pc12)
  dists <- sqrt(rowSums((t(t(pc12) - centroid))^2))
  threshold <- mean(dists) + 3 * sd(dists)
  is_outlier <- dists > threshold
  outlier_flags[meta_matched$sample_id[idx][is_outlier]] <- TRUE
  
  n_out <- sum(is_outlier)
  if (n_out > 0) {
    cat("  ", ds, ":", n_out, "outlier(s) flagged\n")
  } else {
    cat("  ", ds, ": no outliers\n")
  }
}

# PCA plot across all datasets
pca_all <- prcomp(t(logcpm), center = TRUE, scale. = FALSE)
pca_df <- data.frame(
  PC1 = pca_all$x[, 1], PC2 = pca_all$x[, 2],
  dataset = meta_matched$dataset,
  group = meta_matched$group_binary,
  outlier = outlier_flags[meta_matched$sample_id]
)
var_expl <- summary(pca_all)$importance[2, 1:2] * 100

pdf(file.path(QCDIR, "pca_all_datasets.pdf"), width = 10, height = 7)
p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = dataset, shape = group)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_point(data = pca_df[pca_df$outlier, ], size = 4, shape = 1, color = "red", stroke = 1.5) +
  labs(
    title = "All Mouse Samples — PCA Before Batch Correction",
    subtitle = sprintf("Outliers circled in red (>3 SD from per-dataset centroid)"),
    x = sprintf("PC1 (%.1f%%)", var_expl[1]),
    y = sprintf("PC2 (%.1f%%)", var_expl[2])
  ) +
  theme_minimal(base_size = 12) +
  theme(text = element_text(family = "sans"))
print(p)
dev.off()
cat("PCA plot saved\n")

# ===== LIBRARY SIZE QC =====
cat("\n===== LIBRARY SIZE QC =====\n")
lib_sizes <- colSums(merged)
low_lib <- lib_sizes < 1e6
cat("Samples with < 1M total counts:", sum(low_lib), "\n")

# ===== QC SUMMARY =====
qc_report <- data.table(
  sample_id    = meta_matched$sample_id,
  dataset      = meta_matched$dataset,
  diet_model   = meta_matched$diet_model,
  group_binary = meta_matched$group_binary,
  lib_size     = lib_sizes,
  pca_outlier  = outlier_flags[meta_matched$sample_id],
  low_lib      = low_lib,
  pass_qc      = !outlier_flags[meta_matched$sample_id] & !low_lib
)

cat("\n===== QC SUMMARY =====\n")
cat("Total samples:", nrow(qc_report), "\n")
cat("Pass PCA:", sum(!qc_report$pca_outlier), "\n")
cat("Pass library size:", sum(!qc_report$low_lib), "\n")
cat("Pass ALL:", sum(qc_report$pass_qc), "\n")
cat("Failed:", sum(!qc_report$pass_qc), "\n")

fwrite(qc_report, file.path(QCDIR, "sample_qc_report.csv"))
cat("\nQC report written:", file.path(QCDIR, "sample_qc_report.csv"), "\n")

# ===== SAVE RDS =====
saveRDS(merged, file.path(RDIR, "merged_counts_raw.rds"))
saveRDS(meta_matched, file.path(RDIR, "meta_matched.rds"))
cat("Saved merged counts and matched metadata as RDS\n")
