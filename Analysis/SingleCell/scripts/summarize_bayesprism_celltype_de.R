#!/usr/bin/env Rscript
# summarize_bayesprism_celltype_de.R
# ---------------------------------------------------------------------------
# Combine the 13 cell types x 4 F-transition consecutive-contrast DE outputs
# from run_bayesprism_celltype_de.R into a long-form table:
#   cell_type, contrast, gene, logFC, padj
#
# Output:
#   Analysis/SingleCell/results_gpu_v2/pseudobulk_de/bayesprism_celltype_de_summary.csv
#
# Also prints per-(cell_type, contrast) and per-cell_type DEG counts
# (padj < 0.05, |logFC| > 0.5) and writes a counts CSV alongside the summary.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

BASE <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)

OUT_DIR <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/pseudobulk_de")

CELL_TYPES <- c("Hepatocyte", "Macrophage", "Stellate", "Endothelial",
                "Cholangiocyte", "T_cell", "B_cell", "Plasma_cell",
                "NK_cell", "Monocyte", "DC", "Neutrophil", "Other_immune")
CONTRASTS <- c("F1_vs_F0", "F2_vs_F1", "F3_vs_F2", "F4_vs_F3")

cat("=== summarize_bayesprism_celltype_de ===\n")
cat(sprintf("Started: %s\n", Sys.time()))

rows <- list()
missing <- character()
for (ct in CELL_TYPES) {
  for (contrast in CONTRASTS) {
    fp <- file.path(OUT_DIR, sprintf("%s_%s_bayesprism_de.csv", ct, contrast))
    if (!file.exists(fp)) {
      missing <- c(missing, fp)
      next
    }
    dt <- fread(fp, select = c("gene", "logFC", "padj"))
    dt[, cell_type := ct]
    dt[, contrast := contrast]
    rows[[length(rows) + 1L]] <- dt[, .(cell_type, contrast, gene, logFC, padj)]
  }
}

if (length(missing) > 0) {
  cat(sprintf("WARN: %d expected DE files missing (e.g. low-power contrasts):\n",
              length(missing)))
  for (m in missing) cat("  ", m, "\n")
}

if (length(rows) == 0) {
  stop("No DE result files found in: ", OUT_DIR)
}

big <- rbindlist(rows)
out_path <- file.path(OUT_DIR, "bayesprism_celltype_de_summary.csv")
fwrite(big, out_path)
cat(sprintf("Wrote: %s (rows=%d)\n", out_path, nrow(big)))

# Per-(cell_type, contrast) and per-cell_type significant counts
counts_pair <- big[, .(
  n_genes = .N,
  n_sig   = sum(padj < 0.05 & abs(logFC) > 0.5, na.rm = TRUE),
  n_up    = sum(padj < 0.05 & logFC >  0.5, na.rm = TRUE),
  n_down  = sum(padj < 0.05 & logFC < -0.5, na.rm = TRUE)
), by = .(cell_type, contrast)][order(cell_type, contrast)]

counts_ct <- big[, .(
  n_sig_total = sum(padj < 0.05 & abs(logFC) > 0.5, na.rm = TRUE),
  n_unique_sig_genes = uniqueN(gene[padj < 0.05 & abs(logFC) > 0.5])
), by = cell_type][order(-n_sig_total)]

fwrite(counts_pair,
       file.path(OUT_DIR, "bayesprism_celltype_de_counts_per_contrast.csv"))
fwrite(counts_ct,
       file.path(OUT_DIR, "bayesprism_celltype_de_counts_per_celltype.csv"))

cat("\nPer-(cell_type, contrast) significant gene counts (padj<0.05, |logFC|>0.5):\n")
print(counts_pair)

cat("\nPer-cell_type total DEG sums (padj<0.05, |logFC|>0.5):\n")
print(counts_ct)

cat(sprintf("\nDone: %s\n", Sys.time()))
