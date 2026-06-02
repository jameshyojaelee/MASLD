#!/usr/bin/env Rscript
# 118b_bayesprism_hires.R
# ---------------------------------------------------------------------------
# BayesPrism high-resolution deconvolution:
# Infer cell-type-specific gene expression for each of 1,444 bulk samples.
#
# BayesPrism jointly estimates:
#   1. Cell-type proportions (theta) — fraction of each cell type per sample
#   2. Cell-type-specific expression (Z) — per-sample, per-cell-type expression
#
# This gives us the ability to ask: "What is the hepatocyte-specific expression
# of gene X in sample Y?" — enabling cell-type-resolved progression analysis.
#
# Input:
#   - results/progression/scrna_reference/scrna_counts.csv.gz
#   - results/progression/scrna_reference/scrna_cell_labels.csv
#   - results/integration/merged_dge.rds (1,444 bulk samples)
#   - results/staging_classifier/modeling_metadata.csv
#
# Output (to results/progression/cibersortx_celltype_expression/):
#   - bayesprism_proportions.csv     (1444 x n_celltypes)
#   - bayesprism_<celltype>.csv.gz   (per-cell-type expression, 1444 x genes)
#   - bayesprism_summary.csv         (QC metrics per cell type)
#
# SLURM: cpu, 16 CPUs, 210G RAM, 48h
# Env:   micromamba activate rnaseq
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(BayesPrism)
  library(data.table)
  library(parallel)
})

cat("=== 118b: BayesPrism High-Resolution Deconvolution ===\n")
cat(sprintf("Started: %s\n", Sys.time()))
cat(sprintf("Cores: %d\n\n", parallel::detectCores()))

# ── Paths ──────────────────────────────────────────────────────────────────
BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
INTEG <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
REF_DIR <- file.path(INTEG, "results/progression/scrna_reference")
OUTDIR <- file.path(INTEG, "results/progression/cibersortx_celltype_expression")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

n_cores <- min(parallel::detectCores(), 16)

# ── 1. Load scRNA reference ───────────────────────────────────────────────
cat("Loading scRNA reference...\n")
sc_counts <- fread(file.path(REF_DIR, "scrna_counts.csv.gz"), header = TRUE)
sc_labels <- fread(file.path(REF_DIR, "scrna_cell_labels.csv"))

# First column is cell_id (row names)
cell_ids_sc <- sc_counts[[1]]
sc_counts <- as.matrix(sc_counts[, -1, with = FALSE])
rownames(sc_counts) <- cell_ids_sc

cell_types <- sc_labels$cell_type
names(cell_types) <- sc_labels$cell_id

# Ensure alignment
stopifnot(all(rownames(sc_counts) %in% names(cell_types)))
cell_types <- cell_types[rownames(sc_counts)]

cat(sprintf("  scRNA reference: %d cells x %d genes\n", nrow(sc_counts), ncol(sc_counts)))
cat(sprintf("  Cell types: %s\n", paste(sort(unique(cell_types)), collapse = ", ")))

# ── 2. Load bulk expression ───────────────────────────────────────────────
cat("\nLoading bulk expression...\n")
dge <- readRDS(file.path(INTEG, "results/integration/merged_dge.rds"))

# Extract raw counts (BayesPrism needs counts, not logCPM)
if ("counts" %in% names(dge)) {
  bulk_counts <- dge$counts
} else {
  # DGEList object
  bulk_counts <- dge$counts
}

cat(sprintf("  Bulk: %d genes x %d samples\n", nrow(bulk_counts), ncol(bulk_counts)))

# Load metadata
meta <- fread(file.path(INTEG, "results/staging_classifier/modeling_metadata.csv"))
cat(sprintf("  Metadata: %d samples\n", nrow(meta)))

# ── 2b. Map bulk Ensembl IDs to gene symbols ──────────────────────────────
cat("\nMapping bulk Ensembl IDs to gene symbols...\n")

# Use the gene annotation from the DGE object if available, otherwise use GTF
if (!is.null(dge$genes) && "external_gene_name" %in% colnames(dge$genes)) {
  id_to_symbol <- dge$genes$external_gene_name
  names(id_to_symbol) <- rownames(dge$genes)
} else if (!is.null(dge$genes) && "gene_name" %in% colnames(dge$genes)) {
  id_to_symbol <- dge$genes$gene_name
  names(id_to_symbol) <- rownames(dge$genes)
} else {
  # Build mapping from GENCODE GTF annotation
  gtf_file <- "/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"
  if (file.exists(gtf_file)) {
    cat("  Building Ensembl→symbol map from GENCODE v49 GTF...\n")
    gtf_lines <- fread(cmd = paste("zcat", gtf_file, "| grep -P '^chr.*\\tgene\\t' | head -70000"),
                       header = FALSE, sep = "\t")
    # Parse gene_id and gene_name from attributes (col 9)
    gene_ids <- gsub('.*gene_id "([^"]+)".*', "\\1", gtf_lines$V9)
    gene_names <- gsub('.*gene_name "([^"]+)".*', "\\1", gtf_lines$V9)
    id_to_symbol <- gene_names
    names(id_to_symbol) <- gene_ids
    cat(sprintf("  GTF mapping: %d genes\n", length(id_to_symbol)))
  } else {
    stop("Cannot build gene name mapping: no GTF file found")
  }
}

# Map bulk rownames (Ensembl with version) to symbols
bulk_symbols <- id_to_symbol[rownames(bulk_counts)]
# For unmapped, try without version
unmapped <- is.na(bulk_symbols)
if (any(unmapped)) {
  base_ids <- sub("\\.[0-9]+$", "", rownames(bulk_counts)[unmapped])
  bulk_symbols[unmapped] <- id_to_symbol[base_ids]
}

n_mapped <- sum(!is.na(bulk_symbols))
cat(sprintf("  Mapped: %d / %d bulk genes to symbols\n", n_mapped, length(bulk_symbols)))

# Replace rownames with symbols (drop unmapped and duplicates)
valid <- !is.na(bulk_symbols) & !duplicated(bulk_symbols)
bulk_counts <- bulk_counts[valid, ]
rownames(bulk_counts) <- bulk_symbols[valid]
cat(sprintf("  Bulk after symbol mapping: %d genes x %d samples\n",
            nrow(bulk_counts), ncol(bulk_counts)))

# ── 3. Align genes between bulk and scRNA ──────────────────────────────────
cat("\nAligning genes...\n")
common_genes <- intersect(colnames(sc_counts), rownames(bulk_counts))
cat(sprintf("  Common genes: %d\n", length(common_genes)))

sc_counts <- sc_counts[, common_genes]
bulk_counts <- bulk_counts[common_genes, ]

# Transpose bulk to samples x genes (BayesPrism convention)
bulk_matrix <- t(as.matrix(bulk_counts))
cat(sprintf("  Bulk matrix: %d samples x %d genes\n", nrow(bulk_matrix), ncol(bulk_matrix)))
cat(sprintf("  scRNA matrix: %d cells x %d genes\n", nrow(sc_counts), ncol(sc_counts)))

# ── 4. Build BayesPrism reference ─────────────────────────────────────────
cat("\nBuilding BayesPrism reference...\n")

# Select marker genes (BayesPrism's built-in selection)
# select.marker identifies genes that distinguish cell types
sc_dat <- new.prism(
  reference = sc_counts,
  mixture = bulk_matrix,
  input.type = "count.matrix",
  cell.type.labels = cell_types,
  cell.state.labels = cell_types,  # Same as cell.type for simplicity
  key = NULL,
  outlier.cut = 0.01,
  outlier.fraction = 0.1
)

cat(sprintf("  Prism object created\n"))

# ── 5. Run BayesPrism ─────────────────────────────────────────────────────
cat("\nRunning BayesPrism (this will take a while)...\n")
cat(sprintf("  Using %d cores\n", n_cores))
cat(sprintf("  Time: %s\n", Sys.time()))

bp_result <- run.prism(
  prism = sc_dat,
  n.cores = n_cores,
  update.gibbs = TRUE
)

cat(sprintf("  BayesPrism complete at %s\n", Sys.time()))

# Save BayesPrism object immediately (prevents loss if extraction crashes)
saveRDS(bp_result, file.path(OUTDIR, "bayesprism_result.rds"))
cat("  Saved bayesprism_result.rds (for recovery)\n")

# ── 6. Extract results ────────────────────────────────────────────────────
cat("\n=== Extracting Results ===\n")

# Cell-type proportions (theta)
proportions <- get.fraction(
  bp = bp_result,
  which.theta = "final",
  state.or.type = "type"
)
cat(sprintf("  Proportions: %d samples x %d types\n",
            nrow(proportions), ncol(proportions)))

# Save proportions
prop_df <- as.data.frame(proportions)
prop_df$sample_id <- rownames(prop_df)
prop_df <- prop_df[, c("sample_id", setdiff(names(prop_df), "sample_id"))]
fwrite(prop_df, file.path(OUTDIR, "bayesprism_proportions.csv"))
cat("  Saved bayesprism_proportions.csv\n")

# Cell-type-specific expression (Z matrix)
cell_type_names <- colnames(proportions)

summary_rows <- list()
for (ct in cell_type_names) {
  cat(sprintf("  Extracting %s expression...\n", ct))

  # get.exp returns samples x genes matrix for a specific cell type
  ct_expr <- get.exp(
    bp = bp_result,
    state.or.type = "type",
    cell.name = ct
  )

  # Save
  ct_safe <- gsub("[^A-Za-z0-9_]", "_", ct)
  ct_df <- as.data.frame(ct_expr)
  ct_df$sample_id <- rownames(ct_df)
  ct_df <- ct_df[, c("sample_id", setdiff(names(ct_df), "sample_id"))]

  outfile <- file.path(OUTDIR, sprintf("bayesprism_%s.csv.gz", ct_safe))
  fwrite(ct_df, outfile)
  cat(sprintf("    Saved %s (%d x %d)\n", basename(outfile), nrow(ct_df), ncol(ct_df) - 1))

  # QC: correlation with bulk (align genes first)
  ct_genes <- colnames(ct_expr)
  common_qc_genes <- intersect(colnames(bulk_matrix), ct_genes)
  if (length(common_qc_genes) > 100) {
    bulk_mean <- colMeans(bulk_matrix[, common_qc_genes, drop = FALSE])
    ct_mean <- colMeans(ct_expr[, common_qc_genes, drop = FALSE])
    cor_val <- cor(bulk_mean, ct_mean, method = "spearman")
  } else {
    cor_val <- NA
  }

  summary_rows[[ct]] <- data.frame(
    cell_type = ct,
    n_genes = ncol(ct_expr),
    mean_proportion = mean(proportions[, ct]),
    median_proportion = median(proportions[, ct]),
    cor_with_bulk = cor_val,
    stringsAsFactors = FALSE
  )
}

# ── 7. Save summary ───────────────────────────────────────────────────────
summary_df <- do.call(rbind, summary_rows)
fwrite(summary_df, file.path(OUTDIR, "bayesprism_summary.csv"))
cat("\n  Saved bayesprism_summary.csv\n")

# ── 8. Compare with existing MuSiC proportions ────────────────────────────
existing_prop <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/cell_type_proportions.csv")
if (file.exists(existing_prop)) {
  cat("\nComparing with existing MuSiC proportions...\n")
  music_prop <- fread(existing_prop)

  # Align samples
  common_samples <- intersect(rownames(proportions), music_prop$sample)
  if (length(common_samples) > 10) {
    cat(sprintf("  Common samples: %d\n", length(common_samples)))
    # Compare hepatocyte proportions
    bp_hep <- proportions[common_samples, "Hepatocyte"]
    music_hep_col <- grep("Hepatocyte|hepatocyte", names(music_prop), value = TRUE)[1]
    if (!is.na(music_hep_col)) {
      music_aligned <- music_prop[match(common_samples, music_prop$sample), ]
      music_hep <- music_aligned[[music_hep_col]]
      rho <- cor(bp_hep, music_hep, method = "spearman", use = "complete.obs")
      cat(sprintf("  Hepatocyte proportion correlation (BayesPrism vs MuSiC): rho=%.3f\n", rho))
    }
  }
}

# ── 9. BayesPrism object already saved after Gibbs step ───────────────────
cat(sprintf("\n=== 118b: COMPLETE (%s) ===\n", Sys.time()))
