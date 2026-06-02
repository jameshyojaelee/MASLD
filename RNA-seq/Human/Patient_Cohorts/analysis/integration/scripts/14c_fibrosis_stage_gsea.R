#!/usr/bin/env Rscript
# 14c_fibrosis_stage_gsea.R
# Per-fibrosis-stage pathway enrichment (fgsea) for Figure 2 panel e
# Requires: fibrosis_stage_dream.csv from 14b
#
# Usage: Rscript 14c_fibrosis_stage_gsea.R
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

cat("=== 14c: Per-Fibrosis-Stage GSEA ===\n")
cat("Started:", as.character(Sys.time()), "\n\n")

# --- Load fibrosis stage dream results ---
fib_file <- file.path(SIGS, "fibrosis_stage_dream.csv")
if (!file.exists(fib_file)) {
  stop("fibrosis_stage_dream.csv not found — run 14b first")
}
fib <- fread(fib_file)
cat("Loaded fibrosis stage dream:", nrow(fib), "rows\n")
cat("Contrasts:", paste(unique(fib$contrast), collapse = ", "), "\n\n")

# --- Load Hallmark gene sets (Ensembl IDs) ---
cat("Loading MSigDB Hallmark gene sets...\n")
hallmark_df <- msigdbr(species = "Homo sapiens", collection = "H")
hallmark <- split(hallmark_df$ensembl_gene, hallmark_df$gs_name)
hallmark <- lapply(hallmark, function(x) unique(x[!is.na(x) & x != ""]))
cat("Hallmark pathways:", length(hallmark), "\n\n")

# --- Run fgsea per fibrosis stage contrast ---
# Strip Ensembl version for matching
fib[, gene_base := gsub("\\..*", "", gene)]

all_gsea <- list()
for (ctr in unique(fib$contrast)) {
  cat("Running fgsea for", ctr, "...\n")
  sub <- fib[contrast == ctr]

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
fwrite(gsea_combined, file.path(SIGS, "fibrosis_stage_gsea.csv"))
cat("\nGSEA results saved:", nrow(gsea_combined), "rows\n")
cat("=== 14c completed:", as.character(Sys.time()), "===\n")
