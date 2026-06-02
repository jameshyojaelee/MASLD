#!/usr/bin/env Rscript
# Combine all coarse-stage pseudobulk DE outputs into a long-form summary.
#
# Inputs:  results_gpu_v2/pseudobulk_de/{CellType}_{Contrast}_de.csv
# Output:  results_gpu_v2/pseudobulk_de/coarse_stage_de_summary.csv
# Schema:  cell_type, contrast, gene, logFC, padj

suppressPackageStartupMessages({
  library(data.table)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
de_dir <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/pseudobulk_de")
out_file <- file.path(de_dir, "coarse_stage_de_summary.csv")

CONTRASTS <- c("Steatosis_vs_Healthy",
               "Steatohepatitis_vs_Steatosis",
               "Cirrhosis_vs_Steatohepatitis")

CELL_TYPES <- c(
  "Hepatocytes", "Macrophages", "Fibroblasts", "Cholangiocytes",
  "Endothelial_cells", "T_cells", "B_cells", "Mono+mono_derived_cells",
  "Plasma_cells", "Resident_NK", "Circulating_NK_NKT"
)

all_rows <- list()
missing <- c()
counts <- list()

for (ct in CELL_TYPES) {
  for (cn in CONTRASTS) {
    f <- file.path(de_dir, paste0(ct, "_", cn, "_de.csv"))
    if (!file.exists(f)) {
      missing <- c(missing, basename(f))
      next
    }
    dt <- fread(f)
    keep_cols <- c("cell_type", "contrast", "gene", "logFC", "padj")
    if (!all(keep_cols %in% colnames(dt))) {
      message("Schema mismatch: ", basename(f),
              " (cols: ", paste(colnames(dt), collapse = ","), ")")
      next
    }
    dt <- dt[, ..keep_cols]
    all_rows[[length(all_rows) + 1L]] <- dt
    n_sig <- sum(dt$padj < 0.05 & abs(dt$logFC) > 0.5, na.rm = TRUE)
    counts[[length(counts) + 1L]] <- data.table(
      cell_type = ct, contrast = cn,
      n_genes = nrow(dt), n_sig_padj05_lfc05 = n_sig
    )
  }
}

if (length(all_rows) == 0) {
  stop("No DE files found in ", de_dir)
}

summary_dt <- rbindlist(all_rows, use.names = TRUE, fill = TRUE)
fwrite(summary_dt, out_file)
message(sprintf("Wrote %s (%d rows from %d files)",
                out_file, nrow(summary_dt), length(all_rows)))

if (length(missing) > 0) {
  message("Missing files:")
  for (m in missing) message("  ", m)
}

# Print per-contrast significance counts
sig_dt <- rbindlist(counts)
message("\nPer-(cell_type, contrast) significance counts (padj<0.05, |logFC|>0.5):")
print(sig_dt)

# Also save the counts table for downstream use
counts_file <- file.path(de_dir, "coarse_stage_de_sig_counts.csv")
fwrite(sig_dt, counts_file)
message("Wrote ", counts_file)
