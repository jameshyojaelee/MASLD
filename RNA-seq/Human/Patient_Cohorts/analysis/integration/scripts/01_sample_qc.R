#!/usr/bin/env Rscript
# 01_sample_qc.R
# ---------------------------------------------------------------------------
# Sample-level QC: Sex check, PCA outlier detection, alignment QC.
# Input:  unified_metadata.csv + count matrices
# Output: qc/sample_qc_report.csv, qc/sex_check.pdf, qc/pca_outliers.pdf
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(ggplot2)
  library(yaml)
})

# Resolve project root: env var override or default absolute path
PROJECT_ROOT <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
BASE <- file.path(PROJECT_ROOT, "RNA-seq/Human/Patient_Cohorts")
SOURCE_INT <- file.path(BASE, "analysis/integration")
RUN_ROOT <- Sys.getenv("MASLD_RUN_ROOT", "")
REMEDIATION_MODE <- nzchar(RUN_ROOT)
if (REMEDIATION_MODE) {
  RUN_ROOT <- normalizePath(RUN_ROOT, mustWork = TRUE)
  if (!file.exists(file.path(RUN_ROOT, ".bg001_candidate_root"))) {
    stop("MASLD_RUN_ROOT lacks .bg001_candidate_root sentinel: ", RUN_ROOT)
  }
  forbidden <- c(
    normalizePath(SOURCE_INT, mustWork = TRUE),
    normalizePath(file.path(BASE, "results"), mustWork = FALSE)
  )
  if (RUN_ROOT %in% forbidden) stop("Refusing canonical remediation output root: ", RUN_ROOT)
  INT <- RUN_ROOT
} else {
  INT <- SOURCE_INT
}
QC <- file.path(INT, "qc")
dir.create(QC, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(INT, "results/integration"), recursive = TRUE, showWarnings = FALSE)

# --- Load dataset config ---
# Count file paths are auto-derived from config so new datasets are picked up
# automatically when added to config/human_datasets.yaml.
cfg_path <- Sys.getenv("MASLD_CONFIG_PATH", file.path(PROJECT_ROOT, "config/human_datasets.yaml"))
if (!file.exists(cfg_path)) stop("Config not found: ", cfg_path)
cfg <- yaml.load_file(cfg_path)$datasets

# --- Load metadata ---
metadata_path <- Sys.getenv("MASLD_METADATA_PATH", file.path(SOURCE_INT, "metadata/unified_metadata.csv"))
meta <- fread(metadata_path, na.strings = c("", "NA"))
cat("Loaded metadata:", nrow(meta), "samples\n")

# --- Load count matrices (auto-derived from config) ---
count_files <- lapply(cfg, function(ds) file.path(BASE, ds$counts_path))
names(count_files) <- names(cfg)

count_overrides_path <- Sys.getenv("MASLD_COUNT_OVERRIDES", "")
if (nzchar(count_overrides_path)) {
  overrides <- fread(count_overrides_path)
  required_override_cols <- c("dataset", "count_path")
  if (!all(required_override_cols %in% names(overrides))) {
    stop("MASLD_COUNT_OVERRIDES requires columns: ", paste(required_override_cols, collapse = ", "))
  }
  if (anyDuplicated(overrides$dataset) || any(!overrides$dataset %in% names(cfg))) {
    stop("Count overrides contain duplicate or unknown datasets")
  }
  if (any(!grepl("^/", overrides$count_path)) || any(!file.exists(overrides$count_path))) {
    stop("Every count override must be an existing absolute path")
  }
  for (i in seq_len(nrow(overrides))) count_files[[overrides$dataset[[i]]]] <- overrides$count_path[[i]]
}
if (REMEDIATION_MODE) {
  if (!exists("overrides") || !setequal(overrides$dataset, names(cfg))) {
    stop("Remediation mode requires frozen overrides for every active dataset")
  }
}

# Keep only datasets present in metadata (skip if counts not yet generated)
datasets_in_meta <- unique(meta$dataset)
count_files <- count_files[names(count_files) %in% datasets_in_meta]
if (REMEDIATION_MODE && any(!file.exists(unlist(count_files, use.names = FALSE)))) {
  stop("Remediation mode requires every effective count matrix to exist")
}

cat("Datasets to load (in metadata):", paste(names(count_files), collapse = ", "), "\n")

load_counts <- function(path) {
  # featureCounts output: first 6 cols are annotation, rest are sample columns
  ct <- fread(path, skip = 1)  # skip the command line header
  # Gene ID is column 1 ("Geneid")
  genes <- ct$Geneid
  # Sample columns: everything after the first 6 annotation columns
  counts_mat <- as.matrix(ct[, 7:ncol(ct), with = FALSE])
  rownames(counts_mat) <- genes
  # Clean column names — featureCounts uses full BAM paths as column names
  colnames(counts_mat) <- gsub(".*/", "", gsub("\\.Aligned\\.sortedByCoord\\.out\\.bam$", "", colnames(counts_mat)))
  if (anyDuplicated(genes) || anyDuplicated(colnames(counts_mat))) {
    stop("Duplicate gene or cleaned sample ID in count matrix: ", path)
  }
  if (REMEDIATION_MODE && nrow(counts_mat) != 86369L) {
    stop("Remediation count matrix does not contain exactly 86,369 genes: ", path)
  }
  counts_mat
}

cat("\nLoading count matrices...\n")
counts_list <- lapply(count_files, function(f) {
  if (file.exists(f)) {
    cat("  Loading:", basename(dirname(dirname(dirname(f)))), "\n")
    load_counts(f)
  } else {
    cat("  MISSING:", f, "\n")
    NULL
  }
})
counts_list <- counts_list[!sapply(counts_list, is.null)]

# --- Merge into single matrix ---
# Find common genes
common_genes <- Reduce(intersect, lapply(counts_list, rownames))
cat("Common genes across datasets:", length(common_genes), "\n")

merged <- do.call(cbind, lapply(counts_list, function(m) m[common_genes, , drop = FALSE]))
if (anyDuplicated(colnames(merged))) stop("Merged count matrices contain duplicate sample IDs")

# Match metadata to count columns
sample_ids_in_counts <- colnames(merged)
meta_matched <- meta[sample_id %in% sample_ids_in_counts]
cat("Samples in counts AND metadata:", nrow(meta_matched), "\n")

# Align
merged <- merged[, meta_matched$sample_id]

if (REMEDIATION_MODE) {
  provenance_dir <- file.path(RUN_ROOT, "provenance")
  dir.create(provenance_dir, recursive = TRUE, showWarnings = FALSE)
  effective_sources <- data.table(
    dataset = names(count_files),
    count_path = unlist(count_files, use.names = FALSE),
    overridden = names(count_files) %in% if (exists("overrides")) overrides$dataset else character()
  )
  effective_sources[, size_bytes := file.info(count_path)$size]
  effective_sources[, mtime_utc := format(file.info(count_path)$mtime, tz = "UTC", usetz = TRUE)]
  hash_file <- function(path) {
    output <- system2("sha256sum", path, stdout = TRUE, stderr = TRUE)
    status <- attr(output, "status")
    if (!is.null(status) && status != 0L) stop("sha256sum failed for ", path)
    strsplit(output[[1]], "[[:space:]]+")[[1]][[1]]
  }
  effective_sources[, sha256 := vapply(count_path, hash_file, character(1))]
  fwrite(effective_sources, file.path(provenance_dir, "effective_count_sources.tsv"), sep = "\t")
}

qc_mode <- Sys.getenv("MASLD_QC_MODE", "recompute")
if (!qc_mode %in% c("locked", "recompute")) stop("MASLD_QC_MODE must be locked or recompute")
if (qc_mode == "locked") {
  if (!REMEDIATION_MODE) stop("MASLD_QC_MODE=locked is remediation-only")
  locked_qc_path <- Sys.getenv("MASLD_LOCKED_QC", file.path(SOURCE_INT, "qc/sample_qc_report.csv"))
  locked_meta_path <- Sys.getenv("MASLD_LOCKED_META", file.path(SOURCE_INT, "results/integration/meta_matched.rds"))
  locked_qc <- fread(locked_qc_path)
  locked_meta <- as.data.table(readRDS(locked_meta_path))
  if (!setequal(colnames(merged), locked_meta$sample_id) || !setequal(colnames(merged), locked_qc$sample_id)) {
    stop("Locked QC/meta sample set does not match effective count matrices")
  }
  locked_meta <- locked_meta[match(colnames(merged), sample_id)]
  locked_qc <- locked_qc[match(colnames(merged), sample_id)]
  if (anyNA(locked_meta$sample_id) || anyNA(locked_qc$sample_id)) stop("Locked sample alignment failed")
  locked_qc[, canonical_total_counts := total_counts]
  locked_qc[, candidate_total_counts := as.numeric(colSums(merged))[match(sample_id, colnames(merged))]]
  dir.create(file.path(INT, "results/integration"), recursive = TRUE, showWarnings = FALSE)
  fwrite(locked_qc, file.path(QC, "sample_qc_report.csv"))
  saveRDS(merged, file.path(INT, "results/integration/merged_counts_raw.rds"))
  saveRDS(locked_meta, file.path(INT, "results/integration/meta_matched.rds"))
  writeLines(capture.output(sessionInfo()), file.path(QC, "session_info.txt"))
  cat("Locked QC/meta preserved; candidate count totals recorded.\n")
  quit(save = "no", status = 0L)
}

# ===== 1. SEX CHECK =====
cat("\n===== SEX CHECK =====\n")
sex_genes <- c("ENSG00000229807.12", "ENSG00000229807")  # XIST
y_genes   <- c("ENSG00000012817.16", "ENSG00000012817")  # DDX3Y

# Try to find the gene IDs (may have version numbers)
find_gene <- function(pattern, genes) {
  hits <- grep(paste0("^", pattern), genes, value = TRUE)
  if (length(hits) > 0) return(hits[1])
  return(NA)
}

xist_id <- find_gene("ENSG00000229807", common_genes)
ddx3y_id <- find_gene("ENSG00000012817", common_genes)
rps4y1_id <- find_gene("ENSG00000129824", common_genes)

sex_check_df <- data.table(
  sample_id = meta_matched$sample_id,
  dataset   = meta_matched$dataset,
  reported_sex = meta_matched$sex
)

if (!is.na(xist_id)) {
  sex_check_df[, xist_cpm := cpm(merged)[xist_id, ]]
  cat("XIST gene found:", xist_id, "\n")
} else {
  cat("WARNING: XIST gene not found in count matrix\n")
  sex_check_df[, xist_cpm := NA_real_]
}

y_check_id <- if (!is.na(ddx3y_id)) ddx3y_id else rps4y1_id
if (!is.na(y_check_id)) {
  sex_check_df[, y_gene_cpm := cpm(merged)[y_check_id, ]]
  cat("Y-chromosome gene found:", y_check_id, "\n")
} else {
  cat("WARNING: No Y-chromosome gene found\n")
  sex_check_df[, y_gene_cpm := NA_real_]
}

# Infer sex using k-means clustering on XIST + Y-gene expression
# (data-driven — no hardcoded thresholds)
if (!is.na(xist_id) && !is.na(y_check_id)) {
  # k-means on log1p-transformed CPM for two clusters (M vs F)
  sex_mat <- cbind(log1p(sex_check_df$xist_cpm), log1p(sex_check_df$y_gene_cpm))
  set.seed(42)
  km <- kmeans(sex_mat, centers = 2, nstart = 25)

  # Identify which cluster is Female (higher XIST) vs Male (higher Y-gene)
  cluster_xist_means <- tapply(sex_check_df$xist_cpm, km$cluster, mean)
  female_cluster <- which.max(cluster_xist_means)
  sex_check_df[, inferred_sex := fifelse(km$cluster == female_cluster, "F", "M")]

  # Track sex source: annotated (has reported_sex) vs inferred_kmeans (no annotation)
  sex_check_df[, sex_source := fifelse(!is.na(reported_sex), "annotated", "inferred_kmeans")]

  # sex_final: use reported_sex when available, inferred_sex otherwise
  sex_check_df[, sex_final := fifelse(!is.na(reported_sex), reported_sex, inferred_sex)]

  # Only flag mismatch if sex is reported AND conflicts with cluster
  sex_check_df[, pass_sex := is.na(reported_sex) | reported_sex == inferred_sex]

  # Print per-dataset sex source summary
  cat("\n--- Sex source per dataset ---\n")
  sex_src_summary <- sex_check_df[, .(
    n_total      = .N,
    n_annotated  = sum(sex_source == "annotated"),
    n_inferred   = sum(sex_source == "inferred_kmeans"),
    n_F          = sum(sex_final == "F"),
    n_M          = sum(sex_final == "M"),
    n_mismatch   = sum(!is.na(reported_sex) & reported_sex != inferred_sex)
  ), by = dataset]
  sex_src_summary[, pct_inferred := round(100 * n_inferred / n_total, 1)]
  print(sex_src_summary)
  cat("NOTE: Datasets with pct_inferred=100 are using expression-inferred sex only.\n")
  cat("      Concordance validation is only possible for annotated datasets.\n\n")

  # Concordance check: for annotated datasets, how well does k-means match?
  annotated_samples <- sex_check_df[sex_source == "annotated"]
  if (nrow(annotated_samples) > 0) {
    concordance <- sum(annotated_samples$reported_sex == annotated_samples$inferred_sex) /
                   nrow(annotated_samples)
    cat(sprintf("k-means concordance with annotated sex: %.1f%% (%d/%d samples)\n",
                100 * concordance, round(concordance * nrow(annotated_samples)),
                nrow(annotated_samples)))
  }

  # Plot
  pdf(file.path(QC, "sex_check.pdf"), width = 12, height = 6)
  p <- ggplot(sex_check_df, aes(x = log1p(xist_cpm), y = log1p(y_gene_cpm),
                                 color = inferred_sex, shape = sex_source)) +
    geom_point(size = 2, alpha = 0.7) +
    scale_color_manual(values = c("F" = "#E91E63", "M" = "#2196F3")) +
    scale_shape_manual(values = c("annotated" = 16, "inferred_kmeans" = 17),
                       labels = c("annotated" = "Reported sex", "inferred_kmeans" = "Inferred (no annotation)")) +
    facet_wrap(~ dataset, nrow = 1) +
    labs(title = "Sex Check: k-means clustering on XIST vs Y-gene (log1p CPM)",
         subtitle = sprintf("Annotated: %d samples | Inferred (XIST/DDX3Y): %d samples",
                            sum(sex_check_df$sex_source == "annotated"),
                            sum(sex_check_df$sex_source == "inferred_kmeans")),
         x = "log1p(XIST CPM)", y = "log1p(Y-gene CPM)",
         color = "Inferred sex", shape = "Sex source") +
    theme_minimal(base_size = 12)
  # Mark annotation mismatches with X
  mismatches <- sex_check_df[!is.na(reported_sex) & reported_sex != inferred_sex]
  if (nrow(mismatches) > 0) {
    p <- p + geom_point(data = mismatches,
                        aes(x = log1p(xist_cpm), y = log1p(y_gene_cpm)),
                        color = "red", size = 4, shape = 4, stroke = 2)
  }
  print(p)
  dev.off()
  cat("Sex check plot saved\n")

  n_mismatch <- sex_check_df[pass_sex == FALSE, .N]
  cat("Sex annotation mismatches:", n_mismatch, "/", nrow(sex_check_df[!is.na(reported_sex)]), "\n")
  n_imputed <- sum(sex_check_df$sex_source == "inferred_kmeans")
  cat("Sex inferred from XIST/DDX3Y (no annotation):", n_imputed, "samples\n")
} else {
  sex_check_df[, inferred_sex := NA_character_]
  sex_check_df[, sex_source  := NA_character_]
  sex_check_df[, sex_final   := NA_character_]
  sex_check_df[, pass_sex := TRUE]
}

# ===== 2. PER-DATASET PCA OUTLIER DETECTION =====
cat("\n===== PCA OUTLIER DETECTION =====\n")
dge <- DGEList(counts = merged)
dge <- calcNormFactors(dge, method = "TMM")
logcpm <- cpm(dge, log = TRUE, prior.count = 1)

pca_results <- list()
outlier_flags <- rep(FALSE, ncol(merged))
names(outlier_flags) <- colnames(merged)

for (ds in unique(meta_matched$dataset)) {
  idx <- which(meta_matched$dataset == ds)
  if (length(idx) < 5) next

  # Top 1000 most variable genes within this dataset
  vars <- apply(logcpm[, idx], 1, var)
  top_genes <- names(sort(vars, decreasing = TRUE))[1:min(1000, length(vars))]

  pca <- prcomp(t(logcpm[top_genes, idx]), scale. = TRUE, center = TRUE)
  pc_scores <- pca$x[, 1:2]

  # Flag outliers: >3 SD from centroid
  centroid <- colMeans(pc_scores)
  dist_from_centroid <- sqrt(rowSums((sweep(pc_scores, 2, centroid))^2))
  threshold <- mean(dist_from_centroid) + 3 * sd(dist_from_centroid)
  is_outlier <- dist_from_centroid > threshold

  if (any(is_outlier)) {
    cat("  ", ds, ": ", sum(is_outlier), " outlier(s) flagged\n", sep = "")
    outlier_flags[colnames(merged)[idx][is_outlier]] <- TRUE
  } else {
    cat("  ", ds, ": no outliers\n", sep = "")
  }

  pca_results[[ds]] <- data.table(
    sample_id = colnames(merged)[idx],
    PC1 = pc_scores[, 1],
    PC2 = pc_scores[, 2],
    dist = dist_from_centroid,
    outlier = is_outlier,
    dataset = ds
  )
}

pca_all <- rbindlist(pca_results)

# Plot PCA per dataset
pdf(file.path(QC, "pca_outliers.pdf"), width = 12, height = 8)
p <- ggplot(pca_all, aes(x = PC1, y = PC2, color = outlier)) +
  geom_point(size = 1.5, alpha = 0.7) +
  facet_wrap(~ dataset, scales = "free") +
  scale_color_manual(values = c("FALSE" = "grey50", "TRUE" = "red")) +
  labs(title = "Per-Dataset PCA — Outlier Detection (>3 SD)") +
  theme_minimal(base_size = 11)
print(p)
dev.off()
cat("PCA outlier plot saved\n")

# ===== 3. ALIGNMENT QC =====
# (Alignment QC from MultiQC is complementary; for now flag based on total counts)
cat("\n===== LIBRARY SIZE QC =====\n")
lib_sizes <- colSums(merged)
low_lib <- lib_sizes < 5e6  # < 5M total counts
cat("Samples with < 5M total counts:", sum(low_lib), "\n")

# ===== COMPILE QC REPORT =====
qc_report <- data.table(
  sample_id = meta_matched$sample_id,
  dataset   = meta_matched$dataset,
  pass_sex  = sex_check_df$pass_sex,
  pass_pca  = !outlier_flags[meta_matched$sample_id],
  pass_libsize = !low_lib[meta_matched$sample_id],
  total_counts = lib_sizes[meta_matched$sample_id]
)
qc_report[, pass_technical := pass_pca & pass_libsize]
qc_report[, pass_all := pass_sex & pass_pca & pass_libsize]

cat("\n===== QC SUMMARY =====\n")
cat("Total samples:", nrow(qc_report), "\n")
cat("Pass sex check:", sum(qc_report$pass_sex), "\n")
cat("Pass PCA:", sum(qc_report$pass_pca), "\n")
cat("Pass library size:", sum(qc_report$pass_libsize), "\n")
cat("Pass technical (PCA + libsize):", sum(qc_report$pass_technical), "\n")
cat("Pass ALL (incl. sex):", sum(qc_report$pass_all), "\n")
cat("Failed technical:", sum(!qc_report$pass_technical), "\n")
cat("Failed sex only:", sum(!qc_report$pass_sex & qc_report$pass_technical), "\n")

fwrite(qc_report, file.path(QC, "sample_qc_report.csv"))
cat("\nQC report written:", file.path(QC, "sample_qc_report.csv"), "\n")

# Merge sex columns into meta_matched
# sex_source: "annotated" (has reported sex) | "inferred_kmeans" (expression-based)
# sex_final:  use reported_sex when available, inferred_sex otherwise
meta_matched <- merge(
  meta_matched,
  sex_check_df[, .(sample_id, inferred_sex, sex_source, sex_final)],
  by = "sample_id", all.x = TRUE
)
cat("Merged inferred_sex / sex_source / sex_final into metadata\n")
cat("  sex_source='annotated'      :", sum(meta_matched$sex_source == "annotated", na.rm = TRUE), "samples\n")
cat("  sex_source='inferred_kmeans':", sum(meta_matched$sex_source == "inferred_kmeans", na.rm = TRUE), "samples\n")

# Also save the merged count matrix and matched metadata for downstream use
saveRDS(merged, file.path(INT, "results/integration/merged_counts_raw.rds"))
saveRDS(meta_matched, file.path(INT, "results/integration/meta_matched.rds"))
if (REMEDIATION_MODE) writeLines(capture.output(sessionInfo()), file.path(QC, "session_info.txt"))
cat("Saved merged counts and matched metadata as RDS\n")
