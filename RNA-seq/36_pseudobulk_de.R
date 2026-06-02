#!/usr/bin/env Rscript
# 36_pseudobulk_de.R
# ---------------------------------------------------------------------------
# Pseudobulk Differential Expression from scVI-Integrated Single-Cell Atlas
#
# Validates hepatocyte-intrinsic DEGs from bulk deconvolution (Script 25)
# by comparing with pseudobulk DE at single-cell resolution.
#
# Method: muscat framework (Crowell et al., 2020, Nature Communications)
#   1. Load scVI-integrated h5ad via zellkonverter
#   2. Create SingleCellExperiment with CellTypist labels
#   3. Aggregate to pseudobulk per sample x cell type
#   4. Run edgeR-based DE (muscat::pbDS) for MASLD vs Control
#   5. Compare pseudobulk DEGs with bulk dream DEGs
#
# Inputs:
#   - Analysis/SingleCell/results/integrated_atlas.h5ad
#   - RNA-seq/Human/.../results/integration/dream_results.csv
#   - RNA-seq/results/causal_inference/deconv_attribution_scores.csv
#
# Outputs (results/pseudobulk/):
#   - pseudobulk_de_hepatocytes.csv
#   - pseudobulk_de_all_celltypes.csv
#   - pseudobulk_vs_bulk_concordance.csv
#   - figures/pseudobulk_volcano_{celltype}.pdf
#   - figures/pseudobulk_bulk_concordance.pdf
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(muscat)
  library(scater)
  library(edgeR)
  library(data.table)
  library(ggplot2)
  library(zellkonverter)
})

cat("=== Script 36: Pseudobulk DE from scVI Atlas ===\n")
cat("Start time:", format(Sys.time()), "\n\n")

# ==============================================================================
# Configuration
# ==============================================================================
BASE_DIR    <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SC_DIR      <- file.path(BASE_DIR, "Analysis/SingleCell/results")
RESULTS_DIR <- file.path(BASE_DIR, "RNA-seq/results/pseudobulk")
FIG_DIR     <- file.path(BASE_DIR, "figures")

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

H5AD_FILE <- file.path(SC_DIR, "integrated_atlas.h5ad")
DREAM_FILE <- file.path(BASE_DIR,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/dream_results.csv")
ATTRIB_FILE <- file.path(BASE_DIR,
  "RNA-seq/results/causal_inference/deconv_attribution_scores.csv")
GENE_CACHE <- file.path(BASE_DIR,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/gene_annotation/human_ensg_to_symbol.tsv")

# Minimum cells per pseudobulk sample
MIN_CELLS <- 10
# Minimum samples per condition for DE
MIN_SAMPLES <- 3

# ==============================================================================
# 1. Load scVI-integrated atlas
# ==============================================================================
cat("--- Step 1: Loading scVI-integrated atlas ---\n")

if (!file.exists(H5AD_FILE)) {
  cat("  ERROR: Integrated atlas not found:", H5AD_FILE, "\n")
  cat("  scVI integration must complete first. Exiting.\n")
  quit(save = "no", status = 1)
}

cat("  Reading h5ad file...\n")
sce <- readH5AD(H5AD_FILE)
cat("  Cells:", ncol(sce), "\n")
cat("  Genes:", nrow(sce), "\n")
cat("  Metadata columns:", paste(names(colData(sce)), collapse = ", "), "\n")

# ==============================================================================
# 2. Prepare SCE for muscat
# ==============================================================================
cat("\n--- Step 2: Preparing SCE for muscat ---\n")

# Identify cell type and condition columns
# CellTypist labels are typically in 'cell_type' or 'predicted_labels'
ct_col <- if ("cell_type" %in% names(colData(sce))) "cell_type" else
          if ("predicted_labels" %in% names(colData(sce))) "predicted_labels" else
          stop("Cannot find cell type column in SCE metadata")
cat("  Cell type column:", ct_col, "\n")

# Condition column: 'condition' or 'disease'
cond_col <- if ("condition" %in% names(colData(sce))) "condition" else
            if ("disease" %in% names(colData(sce))) "disease" else
            if ("disease_status" %in% names(colData(sce))) "disease_status" else NULL

if (is.null(cond_col)) {
  cat("  WARNING: No condition column found. Attempting to infer from dataset metadata.\n")
  cat("  Available columns:", paste(names(colData(sce)), collapse = ", "), "\n")
  quit(save = "no", status = 1)
}
cat("  Condition column:", cond_col, "\n")

# Sample ID column
sample_col <- if ("sample_id" %in% names(colData(sce))) "sample_id" else
              if ("sample" %in% names(colData(sce))) "sample" else
              if ("donor_id" %in% names(colData(sce))) "donor_id" else NULL

if (is.null(sample_col)) {
  cat("  WARNING: No sample_id column found.\n")
  cat("  Available columns:", paste(names(colData(sce)), collapse = ", "), "\n")
  quit(save = "no", status = 1)
}
cat("  Sample column:", sample_col, "\n")

# Set required muscat fields
sce$cluster_id <- as.factor(colData(sce)[[ct_col]])
sce$sample_id  <- as.factor(colData(sce)[[sample_col]])
sce$group_id   <- as.factor(colData(sce)[[cond_col]])

# Report cell type distribution
ct_counts <- table(sce$cluster_id)
cat("\n  Cell type distribution:\n")
for (ct in sort(names(ct_counts))) {
  cat("    ", ct, ":", ct_counts[ct], "\n")
}

# Report condition distribution
cond_counts <- table(sce$group_id)
cat("\n  Condition distribution:\n")
for (cond in names(cond_counts)) {
  cat("    ", cond, ":", cond_counts[cond], "\n")
}

# ==============================================================================
# 3. Aggregate to pseudobulk
# ==============================================================================
cat("\n--- Step 3: Aggregating to pseudobulk ---\n")

# Use raw counts (X layer from h5ad)
# muscat expects integer counts in the 'counts' assay
if (!"counts" %in% assayNames(sce)) {
  if ("X" %in% assayNames(sce)) {
    assay(sce, "counts") <- assay(sce, "X")
  } else {
    cat("  WARNING: No counts assay found. Using first assay.\n")
    assay(sce, "counts") <- assay(sce, 1)
  }
}

# Prepare for aggregation
sce <- prepSCE(sce,
  kid = "cluster_id",
  sid = "sample_id",
  gid = "group_id",
  drop = FALSE
)

# Aggregate
pb <- aggregateData(sce,
  assay = "counts",
  fun = "sum",
  by = c("cluster_id", "sample_id")
)

cat("  Pseudobulk samples:", ncol(pb), "\n")
cat("  Genes:", nrow(pb), "\n")

# Check sample counts per cell type x condition
ei <- metadata(pb)$experiment_info
cat("\n  Samples per cell type x condition:\n")
for (ct in levels(sce$cluster_id)) {
  n_ctrl <- sum(ei$group_id[ei$cluster_id == ct] %in% c("Control", "Healthy", "Normal"))
  n_disease <- sum(!ei$group_id[ei$cluster_id == ct] %in% c("Control", "Healthy", "Normal"))
  cat("    ", ct, ": Control =", n_ctrl, ", Disease =", n_disease, "\n")
}

# ==============================================================================
# 4. Run DE with muscat::pbDS
# ==============================================================================
cat("\n--- Step 4: Running pseudobulk DE ---\n")

# Determine reference level (Control/Healthy)
ctrl_levels <- c("Control", "Healthy", "Normal")
ref_level <- intersect(ctrl_levels, levels(sce$group_id))
if (length(ref_level) == 0) {
  cat("  WARNING: No control level found in:", paste(levels(sce$group_id), collapse = ", "), "\n")
  cat("  Using first level as reference.\n")
  ref_level <- levels(sce$group_id)[1]
} else {
  ref_level <- ref_level[1]
}
cat("  Reference level:", ref_level, "\n")

# Run DE
res <- pbDS(pb,
  method = "edgeR",
  min_cells = MIN_CELLS,
  filter = "both",
  verbose = TRUE
)

# Extract results
cat("\n  Extracting results...\n")
all_results <- list()
cell_types_tested <- names(res$table)

for (ct in cell_types_tested) {
  contrasts <- names(res$table[[ct]])
  for (contrast in contrasts) {
    de_table <- res$table[[ct]][[contrast]]
    if (!is.null(de_table) && nrow(de_table) > 0) {
      de_dt <- as.data.table(de_table)
      de_dt[, cell_type := ct]
      de_dt[, contrast := contrast]
      all_results[[paste(ct, contrast)]] <- de_dt
    }
  }
}

if (length(all_results) == 0) {
  cat("  WARNING: No DE results produced.\n")
  quit(save = "no", status = 1)
}

all_de <- rbindlist(all_results, fill = TRUE)
cat("  Total DE results:", nrow(all_de), "\n")

# Save all cell type results
fwrite(all_de, file.path(RESULTS_DIR, "pseudobulk_de_all_celltypes.csv"))
cat("  Saved: pseudobulk_de_all_celltypes.csv\n")

# ==============================================================================
# 5. Extract hepatocyte results
# ==============================================================================
cat("\n--- Step 5: Hepatocyte-specific results ---\n")

# Find hepatocyte cell type (name may vary)
hep_patterns <- c("Hepatocyte", "hepatocyte", "Hep")
hep_ct <- cell_types_tested[grepl(paste(hep_patterns, collapse = "|"), cell_types_tested)]

if (length(hep_ct) == 0) {
  cat("  WARNING: No hepatocyte cell type found in:", paste(cell_types_tested, collapse = ", "), "\n")
  cat("  Using all results for concordance analysis.\n")
  hep_de <- all_de
} else {
  hep_ct <- hep_ct[1]
  cat("  Hepatocyte cell type:", hep_ct, "\n")
  hep_de <- all_de[cell_type == hep_ct]
  cat("  Hepatocyte DE genes:", nrow(hep_de), "\n")
  cat("  Significant (p_adj < 0.05):", nrow(hep_de[p_adj.loc < 0.05]), "\n")
  cat("  Significant (p_adj < 0.1):", nrow(hep_de[p_adj.loc < 0.1]), "\n")
}

fwrite(hep_de, file.path(RESULTS_DIR, "pseudobulk_de_hepatocytes.csv"))
cat("  Saved: pseudobulk_de_hepatocytes.csv\n")

# ==============================================================================
# 6. Compare with bulk dream DEGs
# ==============================================================================
cat("\n--- Step 6: Concordance with bulk dream DEGs ---\n")

if (file.exists(DREAM_FILE)) {
  dream <- fread(DREAM_FILE)
  dream[, ensembl_id := sub("\\.\\d+$", "", gene)]

  # Map gene IDs — pseudobulk may use symbols or Ensembl
  # Try matching on gene symbol first
  if ("gene" %in% names(hep_de)) {
    # Map Ensembl -> symbol for dream
    if (file.exists(GENE_CACHE)) {
      ann <- fread(GENE_CACHE)
      ann_map <- ann[symbol != "" & !is.na(symbol) & !duplicated(gene_base),
                     .(ensembl_id = gene_base, symbol)]
      dream <- merge(dream, ann_map, by = "ensembl_id", all.x = TRUE)
    }

    # Merge pseudobulk with bulk dream on gene name
    concordance <- merge(
      hep_de[, .(gene, pb_logFC = logFC, pb_padj = p_adj.loc)],
      dream[!is.na(symbol), .(gene = symbol, bulk_logFC = logFC, bulk_padj = padj)],
      by = "gene"
    )

    if (nrow(concordance) == 0) {
      # Try on Ensembl IDs
      concordance <- merge(
        hep_de[, .(gene, pb_logFC = logFC, pb_padj = p_adj.loc)],
        dream[, .(gene = ensembl_id, bulk_logFC = logFC, bulk_padj = padj)],
        by = "gene"
      )
    }
  } else {
    concordance <- data.table()
  }

  if (nrow(concordance) > 0) {
    cat("  Shared genes:", nrow(concordance), "\n")

    # LFC correlation
    r_lfc <- cor(concordance$pb_logFC, concordance$bulk_logFC, use = "complete.obs")
    rho_lfc <- cor(concordance$pb_logFC, concordance$bulk_logFC,
                   method = "spearman", use = "complete.obs")
    cat("  LFC Pearson r:", round(r_lfc, 3), "\n")
    cat("  LFC Spearman rho:", round(rho_lfc, 3), "\n")

    # Direction concordance
    both_sig <- concordance[pb_padj < 0.1 & bulk_padj < 0.1]
    if (nrow(both_sig) > 0) {
      dir_concordance <- mean(sign(both_sig$pb_logFC) == sign(both_sig$bulk_logFC))
      cat("  Direction concordance (both sig):", round(dir_concordance * 100, 1), "%\n")
    }

    # Enrichment of hepatocyte-intrinsic DEGs
    if (file.exists(ATTRIB_FILE)) {
      attrib <- fread(ATTRIB_FILE)
      hep_intrinsic <- attrib[attribution_class == "Hepatocyte_intrinsic", gene]

      pb_sig <- concordance[pb_padj < 0.1, gene]
      in_hep_intrinsic <- sum(pb_sig %in% hep_intrinsic)
      total_tested <- nrow(concordance)
      hep_total <- sum(concordance$gene %in% hep_intrinsic)

      # Fisher's test
      a <- in_hep_intrinsic
      b <- length(pb_sig) - a
      c <- hep_total - a
      d <- total_tested - a - b - c
      fisher <- fisher.test(matrix(c(a, b, c, d), nrow = 2))

      cat("\n  Hepatocyte-intrinsic enrichment in pseudobulk DEGs:\n")
      cat("    Pseudobulk DEGs in hepatocyte-intrinsic:", a, "/", length(pb_sig), "\n")
      cat("    OR:", round(fisher$estimate, 2), "\n")
      cat("    p-value:", format(fisher$p.value, digits = 3), "\n")
    }

    fwrite(concordance, file.path(RESULTS_DIR, "pseudobulk_vs_bulk_concordance.csv"))
    cat("  Saved: pseudobulk_vs_bulk_concordance.csv\n")

    # Concordance scatter plot
    pdf(file.path(FIG_DIR, "pseudobulk_bulk_concordance.pdf"), width = 7, height = 6)
    p <- ggplot(concordance, aes(x = bulk_logFC, y = pb_logFC)) +
      geom_point(alpha = 0.3, size = 0.8, color = "grey40") +
      geom_point(data = concordance[pb_padj < 0.1 & bulk_padj < 0.1],
                 aes(color = "Both significant"), size = 1.2) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
      scale_color_manual(values = c("Both significant" = "#D7191C")) +
      labs(title = "Pseudobulk vs Bulk DE Concordance (Hepatocytes)",
           subtitle = paste0("r = ", round(r_lfc, 3), ", rho = ", round(rho_lfc, 3),
                            ", n = ", nrow(concordance)),
           x = "Bulk Dream log2FC",
           y = "Pseudobulk log2FC",
           color = NULL) +
      theme_bw(base_size = 11) +
      theme(plot.title = element_text(face = "bold", size = 13),
            legend.position = "bottom")
    print(p)
    dev.off()
    cat("  Saved: figures/pseudobulk_bulk_concordance.pdf\n")
  } else {
    cat("  WARNING: No shared genes between pseudobulk and bulk results.\n")
  }
} else {
  cat("  Dream results not found. Skipping concordance.\n")
}

# ==============================================================================
# 7. Per-cell-type volcano plots
# ==============================================================================
cat("\n--- Step 7: Volcano plots ---\n")

for (ct in unique(all_de$cell_type)) {
  ct_data <- all_de[cell_type == ct]
  ct_clean <- gsub("[^[:alnum:]]", "_", ct)

  if (nrow(ct_data) > 0 && "p_adj.loc" %in% names(ct_data)) {
    ct_data[, sig := ifelse(p_adj.loc < 0.05 & abs(logFC) > 0.5, "Significant", "NS")]

    pdf(file.path(FIG_DIR, paste0("pseudobulk_volcano_", ct_clean, ".pdf")),
        width = 7, height = 6)
    p <- ggplot(ct_data, aes(x = logFC, y = -log10(p_val))) +
      geom_point(aes(color = sig), alpha = 0.5, size = 0.8) +
      scale_color_manual(values = c("Significant" = "#D7191C", "NS" = "grey60")) +
      labs(title = paste0("Pseudobulk DE: ", ct),
           x = "log2 Fold Change", y = "-log10(p-value)") +
      theme_bw(base_size = 11) +
      theme(legend.position = "none")
    print(p)
    dev.off()
  }
}
cat("  Volcano plots saved.\n")

# ==============================================================================
# Summary
# ==============================================================================
cat("\n=== Script 36: Pseudobulk DE Summary ===\n")
cat("  Cell types tested:", length(cell_types_tested), "\n")
cat("  Total DE genes (all cell types):", nrow(all_de), "\n")
if (exists("hep_de") && nrow(hep_de) > 0) {
  cat("  Hepatocyte DEGs (p_adj < 0.1):", nrow(hep_de[p_adj.loc < 0.1]), "\n")
}
if (exists("r_lfc")) {
  cat("  Bulk concordance r:", round(r_lfc, 3), "\n")
}
cat("End time:", format(Sys.time()), "\n")
