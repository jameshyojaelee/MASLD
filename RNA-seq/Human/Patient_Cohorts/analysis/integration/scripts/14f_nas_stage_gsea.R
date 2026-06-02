#!/usr/bin/env Rscript
# 14f_nas_stage_gsea.R
# Per-NAS-score pathway enrichment (fgsea) for Figure 2 panel e
# Mirrors 14c_fibrosis_stage_gsea.R but for NAS scores
# Requires: nas_score_dream.csv from 14b
#
# Usage: Rscript 14f_nas_stage_gsea.R
# SLURM: cpu, 4 CPU, 32GB RAM, ~30min

suppressPackageStartupMessages({
  library(fgsea)
  library(msigdbr)
  library(data.table)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
SIGS <- file.path(INT, "results/disease_signatures")

cat("=== 14f: Per-NAS-Score GSEA ===\n")
cat("Started:", as.character(Sys.time()), "\n\n")

# --- Load NAS score dream results ---
nas_file <- file.path(SIGS, "nas_score_dream.csv")
if (!file.exists(nas_file)) {
  stop("nas_score_dream.csv not found — run 14b first")
}
nas <- fread(nas_file)
cat("Loaded NAS score dream:", nrow(nas), "rows\n")
cat("Contrasts:", paste(unique(nas$contrast), collapse = ", "), "\n\n")

# --- Load Hallmark gene sets (Ensembl IDs) ---
cat("Loading MSigDB Hallmark gene sets...\n")
hallmark_df <- msigdbr(species = "Homo sapiens", collection = "H")
hallmark <- split(hallmark_df$ensembl_gene, hallmark_df$gs_name)
hallmark <- lapply(hallmark, function(x) unique(x[!is.na(x) & x != ""]))
cat("Hallmark pathways:", length(hallmark), "\n\n")

# --- Run fgsea per NAS score contrast ---
# Strip Ensembl version for matching
nas[, gene_base := gsub("\\..*", "", gene)]

all_gsea <- list()
for (ctr in unique(nas$contrast)) {
  cat("Running fgsea for", ctr, "...\n")
  sub <- nas[contrast == ctr]

  # Rank by t-statistic
  sub <- sub[order(-t)]
  ranks <- sub$t
  names(ranks) <- sub$gene_base

  # Remove duplicates (keep first = highest |t|)
  ranks <- ranks[!duplicated(names(ranks))]

  res <- fgsea(
    pathways = hallmark,
    stats    = ranks,
    minSize  = 15,
    maxSize  = 500,
    nPermSimple = 10000
  )
  res$contrast <- ctr

  # Convert leadingEdge list to semicolon-separated string
  res[, leadingEdge := sapply(leadingEdge, paste, collapse = ";")]

  sig <- res[padj < 0.05]
  cat("  Significant (padj < 0.05):", nrow(sig), "of", nrow(res), "\n")

  all_gsea[[ctr]] <- res
}

gsea_combined <- rbindlist(all_gsea, fill = TRUE)
fwrite(gsea_combined, file.path(SIGS, "nas_score_gsea.csv"))
cat("\nGSEA results saved:", nrow(gsea_combined), "rows\n")
cat("=== 14f completed:", as.character(Sys.time()), "===\n")
