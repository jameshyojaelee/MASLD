#!/usr/bin/env Rscript
# ==========================================================================
# Cell-type-specific pathway enrichment via fgsea (Hallmark)
# Reads pseudobulk DE t-stats per cell type, runs fgsea, exports CSV.
# Output: results_gpu_v2/fig2_data/celltype_pathway_enrichment.csv
# ==========================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
  library(msigdbr)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
SC_DIR   <- file.path(BASE, "Analysis/SingleCell")
PB_DIR   <- file.path(SC_DIR, "results_gpu_v2/pseudobulk_de")
OUT_DIR  <- file.path(SC_DIR, "results_gpu_v2/fig2_data")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# --------------------------------------------------------------------------
# Load MSigDB Hallmark gene sets (human)
# --------------------------------------------------------------------------
message("Loading MSigDB Hallmark gene sets...")
msig <- msigdbr(species = "Homo sapiens", collection = "H")
pathways <- split(msig$gene_symbol, msig$gs_name)
message("  ", length(pathways), " Hallmark pathways loaded")

# --------------------------------------------------------------------------
# Gene symbol mapping from the MASLD Gene Catalog
# --------------------------------------------------------------------------
atlas_path <- file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
gene_map <- NULL
if (file.exists(atlas_path)) {
  atlas <- fread(atlas_path, select = c("ensembl_id", "human_symbol"))
  atlas[, ensembl_clean := sub("\\..*", "", ensembl_id)]
  gene_map <- unique(atlas[!is.na(human_symbol) & human_symbol != "",
                           .(ensembl_clean, symbol = human_symbol)],
                     by = "ensembl_clean")
  message("  Gene map: ", nrow(gene_map), " Ensembl -> symbol mappings")
}

# --------------------------------------------------------------------------
# Run fgsea per cell type
# --------------------------------------------------------------------------
de_files <- list.files(PB_DIR, pattern = "_de\\.csv$", full.names = TRUE)
message("Found ", length(de_files), " pseudobulk DE files")

all_results <- list()

for (f in de_files) {
  ct_name <- gsub("_de\\.csv$", "", basename(f))
  ct_clean <- gsub("_", " ", gsub("\\+", "+", ct_name))
  message("Processing: ", ct_clean)

  dt <- fread(f)

  # Need a gene identifier and t-statistic
  # Columns: logFC, AveExpr, t_stat, pvalue, padj, B, gene, cell_type, contrast
  if (!"t_stat" %in% names(dt) && "t" %in% names(dt))
    setnames(dt, "t", "t_stat")

  if (!"t_stat" %in% names(dt)) {
    message("  WARNING: No t-stat column, skipping ", ct_clean)
    next
  }

  # Map Ensembl to symbols
  dt[, ensembl_clean := sub("\\..*", "", gene)]
  if (!is.null(gene_map)) {
    dt <- merge(dt, gene_map, by = "ensembl_clean", all.x = TRUE)
    dt[is.na(symbol), symbol := ensembl_clean]
  } else {
    dt[, symbol := ensembl_clean]
  }

  # Remove duplicates (keep highest abs t-stat per symbol)
  dt <- dt[order(-abs(t_stat))]
  dt <- dt[!duplicated(symbol)]

  # Build named ranking vector
  ranks <- setNames(dt$t_stat, dt$symbol)
  ranks <- ranks[!is.na(ranks)]
  ranks <- sort(ranks, decreasing = TRUE)

  if (length(ranks) < 100) {
    message("  WARNING: Only ", length(ranks), " genes ranked, skipping")
    next
  }

  # Run fgsea
  res <- fgsea(pathways = pathways, stats = ranks, minSize = 15, maxSize = 500)
  res <- as.data.table(res)
  res[, cell_type := ct_clean]
  res[, pathway_short := gsub("^HALLMARK_", "", pathway)]
  all_results[[ct_name]] <- res

  n_sig <- sum(res$padj < 0.05, na.rm = TRUE)
  message("  ", nrow(res), " pathways tested, ", n_sig, " significant (padj < 0.05)")
}

# --------------------------------------------------------------------------
# Combine and save
# --------------------------------------------------------------------------
if (length(all_results) > 0) {
  combined <- rbindlist(all_results, fill = TRUE)

  # Drop leadingEdge column (list column, hard to write to CSV)
  if ("leadingEdge" %in% names(combined)) {
    combined[, leadingEdge_genes := sapply(leadingEdge, function(x) paste(head(x, 10), collapse = ";"))]
    combined[, leadingEdge := NULL]
  }

  out_file <- file.path(OUT_DIR, "celltype_pathway_enrichment.csv")
  fwrite(combined, out_file)
  message("\nSaved: ", out_file)
  message("  ", nrow(combined), " rows, ", combined[, uniqueN(cell_type)], " cell types, ",
          combined[, uniqueN(pathway)], " pathways")
  message("  Significant (padj < 0.05): ", sum(combined$padj < 0.05, na.rm = TRUE))
} else {
  message("ERROR: No results to save")
  quit(status = 1)
}

message("\nDone.")
