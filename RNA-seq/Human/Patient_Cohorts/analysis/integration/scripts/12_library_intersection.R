#!/usr/bin/env Rscript
# 12_library_intersection.R
# ---------------------------------------------------------------------------
# Intersect the new unified consensus DEGs with the existing
# final_core_degs.csv library. Identify overlap, new candidates, and
# genes no longer supported.
# Input:  consensus_degs.csv, human_mouse_ortholog_comparison.csv,
#         results/library/final_core_degs.csv
# Output: library_intersection.csv, library_overlap_summary.pdf
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
INT  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR <- file.path(INT, "results/integration")

# --- Load data ---
consensus <- fread(file.path(RDIR, "consensus_degs.csv"))

# NOTE: human_mouse_ortholog_comparison.csv was produced by 11_ortholog_mapping.R,
# which was removed. The cross-species concordance atlas now handles ortholog mapping.
# If the file is absent, use the annotated consensus (with mouse_gene_id from script 17).
ortho_path <- file.path(RDIR, "human_mouse_ortholog_comparison.csv")
if (!file.exists(ortho_path)) {
  cat("WARNING: human_mouse_ortholog_comparison.csv not found.\n")
  cat("  Attempting to use annotated consensus_degs.csv (requires script 17 to have run).\n")
  # Fallback: read mouse_gene_id/mouse_symbol from the annotated consensus CSV
  if ("mouse_gene_id" %in% names(consensus)) {
    ortho <- consensus[!is.na(mouse_gene_id) & mouse_gene_id != "",
                       .(gene_base = gsub("\\..*", "", gene), mouse_gene_id, human_gene_id = gsub("\\..*", "", gene))]
    cat("  Using", nrow(ortho), "orthologs from annotated consensus.\n")
  } else {
    stop("human_mouse_ortholog_comparison.csv missing and consensus has no mouse_gene_id column. ",
         "Run script 17 first, or regenerate the ortholog comparison file.")
  }
} else {
  ortho <- fread(ortho_path)
}
library   <- fread(file.path(BASE, "results/library/final_core_degs.csv"))

cat("Consensus DEGs:", nrow(consensus), "\n")
cat("Ortholog comparison:", nrow(ortho), "\n")
cat("Current library:", nrow(library), "\n")

# --- Identify human genes with mouse orthologs in consensus ---
# Consensus uses human ENSG IDs (with version); library uses mouse ENSMUSG IDs
consensus[, gene_base := gsub("\\..*", "", gene)]

# Get the mouse_gene_id mapping from ortho table
ortho_map <- ortho[!is.na(mouse_gene_id) & mouse_gene_id != "",
                    .(gene_base, mouse_gene_id, human_symbol, mouse_symbol)]
ortho_map <- unique(ortho_map)

# --- Categories ---
# 1. Library genes that are in the consensus (validated by new analysis)
library_mouse_ids <- unique(library$mouse_gene_id)
consensus_mouse_ids <- unique(ortho_map$mouse_gene_id)

# Overlap
overlap_ids <- intersect(library_mouse_ids, consensus_mouse_ids)

# Significant DEGs
sig_mouse <- ortho_map[gene_base %in% consensus[bulk_sig == TRUE, gene_base], mouse_gene_id]

overlap_sig <- intersect(library_mouse_ids, sig_mouse)

# New candidates (in consensus but not in library)
new_candidates <- setdiff(consensus_mouse_ids, library_mouse_ids)

# Library genes not supported (in library but no consensus support)
unsupported <- setdiff(library_mouse_ids, consensus_mouse_ids)

cat("\n===== LIBRARY INTERSECTION SUMMARY =====\n")
cat("Current library size:", length(library_mouse_ids), "\n")
cat("Consensus genes with mouse orthologs:", length(consensus_mouse_ids), "\n")
cat("\nOverlap (library ∩ consensus):", length(overlap_ids), "\n")
cat("  Validated by significant DEGs:", length(overlap_sig), "\n")
cat("\nNew candidates (consensus only):", length(new_candidates), "\n")
cat("Unsupported (library only, no consensus):", length(unsupported), "\n")

# Coverage percentage
coverage <- length(overlap_ids) / length(library_mouse_ids) * 100
cat("\nLibrary coverage by new analysis:", round(coverage, 1), "%\n")

# --- Build detailed intersection table ---
# Start from library
lib_dt <- copy(library)
lib_dt[, in_library := TRUE]

# Merge with consensus info via mouse_gene_id
consensus_info <- ortho_map[, .(mouse_gene_id, gene_base, human_symbol)]
consensus_info <- merge(consensus_info,
                         consensus[, .(gene_base, bulk_sig, bulk_logFC, bulk_padj)],
                         by = "gene_base", all.x = TRUE)
consensus_info <- unique(consensus_info)

intersection <- merge(
  lib_dt[, .(mouse_gene_id, mouse_gene_symbol, mouse_biotype, 
             source_analyses, n_analyses, has_mouse_data, has_human_data, in_library)],
  consensus_info[, .(mouse_gene_id, bulk_sig, bulk_logFC, bulk_padj)],
  by = "mouse_gene_id", all = TRUE
)

# Fill NAs
intersection[is.na(in_library), in_library := FALSE]
intersection[is.na(bulk_sig), bulk_sig := FALSE]

# Classify
intersection[, status := fcase(
  in_library & bulk_sig == TRUE, "Validated",
  in_library & bulk_sig == FALSE, "Unsupported",
  !in_library & bulk_sig == TRUE, "New_candidate",
  default = "Other"
)]

cat("\nDetailed status breakdown:\n")
print(intersection[, .N, by = status][order(-N)])

# --- Summary bar plot ---
status_summary <- intersection[, .N, by = status]
status_summary[, status := factor(status, levels = c("Validated", 
                                                       "New_candidate", "Unsupported", "Other"))]

pdf(file.path(RDIR, "library_overlap_summary.pdf"), width = 8, height = 5)
p <- ggplot(status_summary[!is.na(status)], aes(x = status, y = N, fill = status)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = N), vjust = -0.3, size = 4) +
  scale_fill_manual(values = c(
    "Validated" = "#4CAF50",
    "New_candidate" = "#2196F3",
    "Unsupported" = "#E91E63",
    "Other" = "grey70"
  )) +
  labs(
    title = "Library Intersection: Existing vs New Consensus",
    subtitle = sprintf("Library: %d genes | Consensus (mapped): %d genes | Overlap: %d",
                       length(library_mouse_ids), length(consensus_mouse_ids), length(overlap_ids)),
    x = "", y = "Number of Genes"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    text = element_text(family = "sans"),
    legend.position = "none",
    axis.text.x = element_text(angle = 20, hjust = 1)
  )
print(p)
dev.off()
cat("\nSaved: library_overlap_summary.pdf\n")

# --- Save ---
fwrite(intersection[order(status, bulk_padj)], file.path(RDIR, "library_intersection.csv"))
cat("Saved: library_intersection.csv\n")

# Also save the new candidates separately
new_cands <- intersection[status == "New_candidate"][order(bulk_padj)]
if (nrow(new_cands) > 0) {
  fwrite(new_cands, file.path(RDIR, "new_candidate_genes.csv"))
  cat("Saved: new_candidate_genes.csv (", nrow(new_cands), "genes)\n")
}

cat("\nLibrary intersection analysis complete.\n")
