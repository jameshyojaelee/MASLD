#!/usr/bin/env Rscript
# 16_disease_signatures_consensus.R
# ---------------------------------------------------------------------------
# Consensus disease signature integration:
#   - Combines NAFL-vs-NASH, fibrosis progression, and NAS components
#   - Classifies genes into stage-specific categories
#   - Overlaps with existing library
#   - Sanity checks on known markers
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

RDIR <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/results/disease_signatures"

# ============================================================
#  Load all results
# ============================================================
cat("=== LOADING RESULTS ===\n")

# NAFL vs NASH
nn_dream  <- fread(file.path(RDIR, "nafl_vs_nash_dream.csv"))
nn_meta   <- fread(file.path(RDIR, "nafl_vs_nash_meta.csv"))
nn_cons   <- fread(file.path(RDIR, "nafl_vs_nash_consensus.csv"))

# Fibrosis progression
fib_dream <- fread(file.path(RDIR, "fibrosis_dream.csv"))
fib_slopes <- fread(file.path(RDIR, "fibrosis_slopes_meta.csv"))
fib_early_late <- tryCatch(fread(file.path(RDIR, "fibrosis_early_late.csv")), error = function(e) NULL)

# NAS components
nas_all   <- fread(file.path(RDIR, "nas_components_all.csv"))

# Existing disease-vs-control consensus
disease_cons <- tryCatch(
  fread("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/consensus_degs.csv"),
  error = function(e) NULL
)

cat("Loaded:\n")
cat("  NAFL-vs-NASH dream:", nrow(nn_dream), "genes\n")
cat("  NAFL-vs-NASH meta:", nrow(nn_meta), "genes\n")
cat("  Fibrosis dream:", nrow(fib_dream), "genes\n")
cat("  Fibrosis slopes meta:", nrow(fib_slopes), "genes\n")
cat("  NAS components:", nrow(nas_all), "genes across", length(unique(nas_all$contrast)), "contrasts\n")

# ============================================================
#  Build unified signature table
# ============================================================
cat("\n=== BUILDING UNIFIED SIGNATURE TABLE ===\n")

# Start with all unique genes
all_genes <- unique(c(nn_dream$gene, fib_dream$gene, nas_all$gene))
cat("Total unique genes:", length(all_genes), "\n")

# NAFL vs NASH columns
nn <- nn_dream[, .(gene, nn_lfc = logFC, nn_padj = adj.P.Val)]

# Fibrosis dream columns  
fib <- fib_dream[, .(gene, fib_slope = logFC, fib_padj = adj.P.Val)]

# NAS component columns (pivot wider)
nas_wide <- dcast(nas_all, gene ~ contrast,
  value.var = c("logFC", "adj.P.Val"))

# Merge everything
unified <- merge(data.table(gene = all_genes), nn, by = "gene", all.x = TRUE)
unified <- merge(unified, fib, by = "gene", all.x = TRUE)
unified <- merge(unified, nas_wide, by = "gene", all.x = TRUE)

# Add disease-vs-control info if available
if (!is.null(disease_cons)) {
  dvc_cols <- intersect(names(disease_cons), c("gene", "dream_logFC", "dream_padj", "dream_sig"))
  if (length(dvc_cols) >= 2) {
    dvc <- disease_cons[, ..dvc_cols]
    setnames(dvc, setdiff(dvc_cols, "gene"), paste0("dvc_", setdiff(dvc_cols, "gene")))
    unified <- merge(unified, dvc, by = "gene", all.x = TRUE)
    cat("  Merged disease-vs-control consensus:", nrow(disease_cons), "genes\n")
  }
}

# ============================================================
#  Stage classification
# ============================================================
cat("\n=== STAGE CLASSIFICATION ===\n")

# Significance flags (padj < 0.05)
unified[, `:=`(
  nn_sig = !is.na(nn_padj) & nn_padj < 0.05,
  fib_sig = !is.na(fib_padj) & fib_padj < 0.05
)]

# Classification logic
unified[, stage_category := "Unclassified"]

# 1. Progression: significant in both NAFL-vs-NASH AND fibrosis
unified[nn_sig == TRUE & fib_sig == TRUE, stage_category := "Progression"]

# 2. NASH-specific: significant in NAFL→NASH, NOT significant fibrosis slope
unified[nn_sig == TRUE & fib_sig == FALSE, stage_category := "NASH_specific"]

# 3. Fibrosis-specific: significant fibrosis slope, NOT NAFL→NASH
unified[fib_sig == TRUE & nn_sig == FALSE, stage_category := "Fibrosis_specific"]

# 4. For genes with NAS component data, sub-classify
has_steatosis <- !is.na(unified[["adj.P.Val_Steatosis_vs_Normal"]])
has_inflam <- !is.na(unified[["adj.P.Val_Inflammation_vs_Steatosis"]])

steatosis_sig <- has_steatosis & unified[["adj.P.Val_Steatosis_vs_Normal"]] < 0.05
inflam_sig <- has_inflam & unified[["adj.P.Val_Inflammation_vs_Steatosis"]] < 0.05

unified[steatosis_sig & !inflam_sig & stage_category == "Unclassified",
  stage_category := "Steatosis_only"]
unified[inflam_sig & stage_category == "Unclassified",
  stage_category := "Inflammation_driven"]

cat("Stage classification:\n")
print(unified[, .N, by = stage_category][order(-N)])

# ============================================================
#  Sanity checks
# ============================================================
cat("\n=== SANITY CHECKS ===\n")

# Known NASH markers — search by symbol column if available, otherwise skip
# (gene column contains ENSG IDs; grepping symbols against ENSG IDs always fails)
nash_markers <- c("CYP7A1", "CYP2E1", "PNPLA3", "FASN", "SCD", "SREBF1",
                   "CCL2", "TNF", "IL6", "IL1B", "TGFB1")
if ("symbol" %in% names(nn_dream)) {
  cat("\nKnown NASH markers in NAFL-vs-NASH results (by symbol):\n")
  for (marker in nash_markers) {
    hits <- nn_dream[symbol == marker]
    if (nrow(hits) > 0) {
      cat(sprintf("  %s: LFC=%.3f, padj=%.2e\n", marker, hits$logFC[1], hits$adj.P.Val[1]))
    } else {
      cat(sprintf("  %s: not detected\n", marker))
    }
  }
} else {
  cat("\nNOTE: 'symbol' column not yet present in nn_dream (run script 17 to annotate).\n")
  cat("      Skipping marker sanity check — markers cannot be found by ENSG ID.\n")
}

# Known fibrosis markers
fib_markers <- c("COL1A1", "COL3A1", "ACTA2", "TGFB1", "TIMP1", "LOX",
                  "COL1A2", "MMP2", "SERPINE1")
if ("symbol" %in% names(fib_dream)) {
  cat("\nKnown fibrosis markers in fibrosis dream (by symbol):\n")
  for (marker in fib_markers) {
    hits <- fib_dream[symbol == marker]
    if (nrow(hits) > 0) {
      cat(sprintf("  %s: slope=%.3f, padj=%.2e\n", marker, hits$logFC[1], hits$adj.P.Val[1]))
    } else {
      cat(sprintf("  %s: not detected\n", marker))
    }
  }
} else {
  cat("\nNOTE: 'symbol' column not yet present in fib_dream (run script 17 to annotate).\n")
}

# ============================================================
#  Overlap with existing library
# ============================================================
cat("\n=== LIBRARY OVERLAP ===\n")

lib_file <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/library_intersection.csv"
if (file.exists(lib_file)) {
  lib <- fread(lib_file)

  # The library contains mouse ENSMUSG IDs; nn_dream and fib_dream contain human ENSG IDs.
  # Direct overlap (lib$gene %in% human_ENSG) will always be 0.
  # Use the mouse_gene_id column (added by script 17) to map library → human genes.
  if ("mouse_gene_id" %in% names(unified)) {
    lib_mouse_base <- gsub("\\..*", "", lib[[1]])  # first column = gene ID
    lib_in_unified <- unified[gsub("\\..*", "", mouse_gene_id) %in% lib_mouse_base]
    cat(sprintf("Library genes mapped to unified disease signatures: %d / %d\n",
                nrow(lib_in_unified), nrow(lib)))
    cat("\nLibrary genes by stage category:\n")
    print(lib_in_unified[, .N, by = stage_category][order(-N)])
  } else {
    cat("NOTE: mouse_gene_id column not present in unified table.\n")
    cat("      Run script 17 to annotate, then re-run script 16 for library overlap.\n")
  }
}

# ============================================================
#  Save
# ============================================================
fwrite(unified, file.path(RDIR, "unified_disease_signatures.csv"))
cat("\nSaved: unified_disease_signatures.csv\n")

# Top genes per category
for (cat_name in unique(unified$stage_category)) {
  if (cat_name == "Unclassified") next
  n <- unified[stage_category == cat_name, .N]
  cat(sprintf("  %s: %d genes\n", cat_name, n))
}

# ============================================================
#  Summary visualization
# ============================================================
cat("\nGenerating summary plots...\n")

pdf(file.path(RDIR, "disease_signatures_summary.pdf"), width = 14, height = 10)

# 1. LFC correlation: NAFL-vs-NASH vs Fibrosis slope
common <- merge(nn[!is.na(nn_lfc)], fib[!is.na(fib_slope)], by = "gene")
rho <- cor(common$nn_lfc, common$fib_slope, method = "spearman", use = "complete.obs")

p1 <- ggplot(common, aes(x = nn_lfc, y = fib_slope)) +
  geom_point(size = 0.3, alpha = 0.2, color = "grey40") +
  geom_point(data = common[nn_padj < 0.05 & fib_padj < 0.05],
             size = 0.5, alpha = 0.6, color = "#E91E63") +
  geom_smooth(method = "lm", se = FALSE, color = "#2196F3", linewidth = 0.5) +
  labs(title = sprintf("NAFL→NASH vs Fibrosis Progression (ρ=%.3f)", rho),
       subtitle = "Pink = significant in both",
       x = "NAFL→NASH logFC (Dream)", y = "Fibrosis slope (Dream)") +
  theme_minimal(base_size = 11)
print(p1)

# 2. Stage category bar chart
cat_counts <- unified[stage_category != "Unclassified", .N, by = stage_category]
cat_counts[, stage_category := factor(stage_category,
  levels = c("Progression", "NASH_specific", "Fibrosis_specific",
             "Steatosis_only", "Inflammation_driven"))]
cat_counts <- cat_counts[!is.na(stage_category)]

if (nrow(cat_counts) > 0) {
  p2 <- ggplot(cat_counts, aes(x = reorder(stage_category, -N), y = N, fill = stage_category)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = c(
      "Progression" = "#E91E63",
      "NASH_specific" = "#9C27B0",
      "Fibrosis_specific" = "#FF5722",
      "Steatosis_only" = "#FF9800",
      "Inflammation_driven" = "#3F51B5"
    ), guide = "none") +
    labs(title = "Gene Stage Classification",
         x = NULL, y = "Number of genes") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  print(p2)
}

# 3. Venn-style overlap counts
nn_sig_set <- nn_dream[adj.P.Val < 0.05, gene]
fib_sig_set <- fib_dream[adj.P.Val < 0.05, gene]
both <- intersect(nn_sig_set, fib_sig_set)
nn_only <- setdiff(nn_sig_set, fib_sig_set)
fib_only <- setdiff(fib_sig_set, nn_sig_set)

cat(sprintf("\n  NAFL-vs-NASH only: %d genes\n", length(nn_only)))
cat(sprintf("  Fibrosis only: %d genes\n", length(fib_only)))
cat(sprintf("  Both (progression): %d genes\n", length(both)))

dev.off()
cat("Saved: disease_signatures_summary.pdf\n")

cat("\n=== Script 16 complete ===\n")
