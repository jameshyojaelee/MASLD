#!/usr/bin/env Rscript
# M07_mouse_pathway_analysis.R
# ---------------------------------------------------------------------------
# Gene set enrichment analysis (fgsea) on dream mega-analysis results.
# Uses Mus musculus gene sets from MSigDB.
# Input:  results/meta_analysis/dream_pooled_results.csv
# Output: results/gsea_results.csv, results/gsea_hallmark_barplot.pdf
# ---------------------------------------------------------------------------

# ---- Seed pinning (T2.4, 2026-04-22) -----
set.seed(42)
Sys.setenv(R_PARALLEL_SEED = "42")

suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
  library(msigdbr)
  library(ggplot2)
})

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Mouse/Unified_Integration"
RDIR <- file.path(BASE, "results")

# --- Load dream results ---
cat("Loading dream mega-analysis results...\n")
dream <- fread(file.path(RDIR, "meta_analysis/dream_pooled_results.csv"))
cat("  Genes:", nrow(dream), "\n")

# --- Prepare ranked gene list ---
# Dream output has ENSEMBL gene IDs — strip version suffix for mapping
dream[, gene_base := gsub("\\..*", "", gene)]

# Rank: -log10(pvalue) * sign(logFC)
dream[, rank_metric := -log10(pmax(P.Value, 1e-300)) * sign(logFC)]
dream <- dream[order(-rank_metric)]

# Named vector for fgsea
ranks <- dream$rank_metric
names(ranks) <- dream$gene_base

# --- Get MSigDB gene sets (Mus musculus) ---
cat("Loading MSigDB gene sets for Mus musculus...\n")
hallmark_df <- msigdbr(species = "Mus musculus", collection = "H")
hallmark <- split(hallmark_df$ensembl_gene, hallmark_df$gs_name)

kegg_df <- msigdbr(species = "Mus musculus", collection = "C2", subcollection = "CP:KEGG_MEDICUS")
kegg <- split(kegg_df$ensembl_gene, kegg_df$gs_name)

reactome_df <- msigdbr(species = "Mus musculus", collection = "C2", subcollection = "CP:REACTOME")
reactome <- split(reactome_df$ensembl_gene, reactome_df$gs_name)

cat("  Hallmark:", length(hallmark), "| KEGG:", length(kegg),
    "| Reactome:", length(reactome), "gene sets\n")

# --- Run fgsea ---
run_gsea <- function(pathways, name) {
  cat("\nRunning fgsea on", name, ":", length(pathways), "gene sets...\n")
  res <- fgsea(pathways = pathways, stats = ranks, minSize = 15, maxSize = 500)
  res$collection <- name
  res <- res[order(pval)]
  sig <- res[padj < 0.05]
  cat("  Significant (padj < 0.05):", nrow(sig), "\n")
  if (nrow(sig) > 0) {
    cat("  Top 5:\n")
    for (i in 1:min(5, nrow(sig))) {
      cat("    ", sig$pathway[i], ": NES=", round(sig$NES[i], 2), "\n")
    }
  }
  return(res)
}

gsea_hallmark <- run_gsea(hallmark, "Hallmark")
gsea_kegg     <- run_gsea(kegg, "KEGG")
gsea_reactome <- run_gsea(reactome, "Reactome")

# Combine and save
all_gsea <- rbindlist(list(gsea_hallmark, gsea_kegg, gsea_reactome), fill = TRUE)
# Convert leadingEdge list column to string for CSV
all_gsea[, leadingEdge := sapply(leadingEdge, paste, collapse = ";")]
fwrite(all_gsea, file.path(RDIR, "gsea_results.csv"))
cat("\nSaved:", file.path(RDIR, "gsea_results.csv"), "\n")

# --- Plot top Hallmark pathways ---
top_hallmark <- gsea_hallmark[padj < 0.05][order(NES)]
if (nrow(top_hallmark) > 0) {
  top_hallmark <- head(top_hallmark, 20)
  top_hallmark[, pathway_short := gsub("HALLMARK_", "", pathway)]
  top_hallmark[, pathway_short := gsub("_", " ", pathway_short)]

  pdf(file.path(RDIR, "gsea_hallmark_barplot.pdf"), width = 10, height = 7)
  p <- ggplot(top_hallmark, aes(x = reorder(pathway_short, NES), y = NES, fill = NES > 0)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = c("TRUE" = "#00BCD4", "FALSE" = "#1B5E20"),
                      labels = c("Downregulated", "Upregulated"),
                      name = "Direction") +
    labs(title = "GSEA: Hallmark Pathways — Mouse MASLD (Disease vs Control)",
         subtitle = "Dream mega-analysis (443 samples, 5 diet models)",
         x = "", y = "Normalized Enrichment Score") +
    theme_minimal(base_size = 12)
  print(p)
  dev.off()
  cat("Saved: gsea_hallmark_barplot.pdf\n")
} else {
  cat("No significant Hallmark pathways found.\n")
}

cat("\nPathway analysis complete.\n")
