#!/usr/bin/env Rscript
# 14g_grouped_nas_gsea.R
# Pathway enrichment (fgsea) for the grouped-NAS-activity dose-response.
# Mirrors 14c_fibrosis_stage_gsea.R / 14f_nas_stage_gsea.R exactly, but consumes
# the 3 grouped contrasts (NAS1-2 / NAS3-4 / NAS5-8 vs NAS0) from 14e so that a
# pathway's NES across the three bins reads as a NAS-activity dose-response.
# Requires: nas_grouped_vs_nas0_lvqw.csv from 14e_grouped_nas_stage_lvqw.R
#
# Usage: Rscript 14g_grouped_nas_gsea.R
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

cat("=== 14g: Grouped-NAS-Activity GSEA ===\n")
cat("Started:", as.character(Sys.time()), "\n\n")

# --- Load grouped-NAS LVQW results ---
nas_file <- file.path(SIGS, "nas_grouped_vs_nas0_lvqw.csv")
if (!file.exists(nas_file)) {
  stop("nas_grouped_vs_nas0_lvqw.csv not found — run 14e_grouped_nas_stage_lvqw.R first")
}
nas <- fread(nas_file)
cat("Loaded grouped-NAS LVQW:", nrow(nas), "rows\n")
cat("Contrasts:", paste(unique(nas$contrast), collapse = ", "), "\n\n")

# --- Load Hallmark gene sets (Ensembl IDs) ---
cat("Loading MSigDB Hallmark gene sets...\n")
hallmark_df <- msigdbr(species = "Homo sapiens", collection = "H")
hallmark <- split(hallmark_df$ensembl_gene, hallmark_df$gs_name)
hallmark <- lapply(hallmark, function(x) unique(x[!is.na(x) & x != ""]))
cat("Hallmark pathways:", length(hallmark), "\n\n")

# --- Run fgsea per grouped-NAS contrast ---
# Strip Ensembl version for matching
nas[, gene_base := gsub("\\..*", "", gene)]

# Order the bins so the dose-response reads low -> high NAS activity
ctr_levels <- c("NAS1-2_vs_NAS0", "NAS3-4_vs_NAS0", "NAS5-8_vs_NAS0")
contrasts  <- intersect(ctr_levels, unique(nas$contrast))

all_gsea <- list()
for (ctr in contrasts) {
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
gsea_combined[, contrast := factor(contrast, levels = ctr_levels)]
fwrite(gsea_combined, file.path(SIGS, "nas_grouped_gsea.csv"))
cat("\nGSEA results saved:", nrow(gsea_combined), "rows\n")
cat("=== 14g completed:", as.character(Sys.time()), "===\n")
