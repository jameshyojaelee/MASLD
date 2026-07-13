#!/usr/bin/env Rscript
# 08_pathway_analysis.R
# ---------------------------------------------------------------------------
# Gene set enrichment analysis (fgsea) and overrepresentation analysis
# (clusterProfiler) on canonical bulk DEGs and Tier 1 consensus DEGs.
# Output: results/integration/gsea_results.csv, enrichment plots
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
  library(msigdbr)
  library(org.Hs.eg.db)
  library(ggplot2)
})

# Optional dependency for ORA
has_clusterProfiler <- requireNamespace("clusterProfiler", quietly = TRUE)
if (has_clusterProfiler) {
  library(clusterProfiler)
} else {
  warning("clusterProfiler not found. Overrepresentation analysis (ORA) will be skipped.")
}

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts"
INT  <- file.path(BASE, "analysis/integration")
RDIR <- file.path(INT, "results/integration")

# --- Load canonical bulk DEGs (limma-voom-qw C2) ---
dream <- fread(file.path(RDIR, "canonical_deg_results.csv"))

# --- Prepare ranked gene list (strip version from ENSEMBL IDs) ---
dream[, gene_base := gsub("\\..*", "", gene)]
# Rank: t-statistic (more robust than p-value for GSEA)
dream[, rank_metric := t]
dream <- dream[order(-rank_metric)]

# Named vector for fgsea
ranks <- dream$rank_metric
names(ranks) <- dream$gene_base

# --- Get MSigDB gene sets ---
cat("Loading MSigDB gene sets...\n")
# Use ensembl_gene for matching (same namespace as our ranked list).
# Validate that ensembl_gene is populated; fall back to gene_symbol with a warning.
check_ensembl <- function(df, name) {
  n_mapped <- sum(!is.na(df$ensembl_gene) & df$ensembl_gene != "")
  n_total  <- nrow(df)
  cat(sprintf("  %s: %d/%d members have Ensembl IDs (%.0f%%)\n",
              name, n_mapped, n_total, 100 * n_mapped / max(n_total, 1)))
  if (n_mapped / max(n_total, 1) < 0.5) {
    warning(name, ": fewer than 50% of gene set members have Ensembl IDs. ",
            "GSEA results may be unreliable. Check msigdbr version.")
  }
}

hallmark_df <- msigdbr(species = "Homo sapiens", collection = "H")
check_ensembl(hallmark_df, "Hallmark")
hallmark <- split(hallmark_df$ensembl_gene, hallmark_df$gs_name)
hallmark <- lapply(hallmark, function(x) x[!is.na(x) & x != ""])

# CP:KEGG_MEDICUS introduced in MSigDB v2023.1; fall back to CP:KEGG if absent
kegg_df <- tryCatch(
  msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_MEDICUS"),
  error = function(e) {
    cat("  CP:KEGG_MEDICUS not found, falling back to CP:KEGG\n")
    msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG")
  }
)
if (is.null(kegg_df) || nrow(kegg_df) == 0) {
  cat("WARNING: No KEGG gene sets found. Skipping KEGG GSEA.\n")
  kegg <- list()
} else {
  check_ensembl(kegg_df, "KEGG")
  kegg <- split(kegg_df$ensembl_gene, kegg_df$gs_name)
  kegg <- lapply(kegg, function(x) x[!is.na(x) & x != ""])
}

reactome_df <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME")
check_ensembl(reactome_df, "Reactome")
reactome <- split(reactome_df$ensembl_gene, reactome_df$gs_name)
reactome <- lapply(reactome, function(x) x[!is.na(x) & x != ""])

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
top_hallmark <- gsea_hallmark[padj < 0.05][order(-abs(NES))]
if (nrow(top_hallmark) > 0) {
  top_hallmark <- head(top_hallmark, 20)  # Top 20 by absolute NES (strongest effects, bidirectional)
  top_hallmark[, pathway_short := gsub("HALLMARK_", "", pathway)]
  top_hallmark[, pathway_short := gsub("_", " ", pathway_short)]

  pdf(file.path(RDIR, "gsea_hallmark_barplot.pdf"), width = 10, height = 7)
  p <- ggplot(top_hallmark, aes(x = reorder(pathway_short, NES), y = NES, fill = NES > 0)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = c("TRUE" = "#E91E63", "FALSE" = "#2196F3"),
                      labels = c("Downregulated", "Upregulated"),
                      name = "Direction") +
    labs(title = "GSEA: Hallmark Pathways (Disease vs Control)",
         x = "", y = "Normalized Enrichment Score") +
    theme_minimal(base_size = 12)
  print(p)
  dev.off()
  cat("Saved: gsea_hallmark_barplot.pdf\n")
}

# --- Overrepresentation on Significant DEGs ---
if (has_clusterProfiler) {
  tier1 <- fread(file.path(RDIR, "dream_significant_degs.csv"))
  if (nrow(tier1) > 0) {
    cat("\n===== ORA on Significant DEGs =====\n")
    tier1[, gene_base := gsub("\\..*", "", gene)]

    # Map ENSEMBL to ENTREZ — both query genes AND background universe
    # Background = all tested genes from canonical bulk DEGs (not the full genome)
    dream_bg <- fread(file.path(RDIR, "canonical_deg_results.csv"))
    dream_bg[, gene_base := gsub("\\..*", "", gene)]
    bg_mapping <- bitr(dream_bg$gene_base, fromType = "ENSEMBL", toType = "ENTREZID",
                       OrgDb = org.Hs.eg.db)

    mapping <- bitr(tier1$gene_base, fromType = "ENSEMBL", toType = "ENTREZID",
                    OrgDb = org.Hs.eg.db)

    if (nrow(mapping) > 0) {
      ego <- enrichGO(gene = mapping$ENTREZID,
                      universe = bg_mapping$ENTREZID,  # tested genes, not full genome
                      OrgDb = org.Hs.eg.db,
                      ont = "BP",
                      pAdjustMethod = "BH",
                      qvalueCutoff = 0.05,
                      readable = TRUE)

      if (!is.null(ego) && nrow(ego@result[ego@result$p.adjust < 0.05, ]) > 0) {
        cat("Significant GO BP terms:", nrow(ego@result[ego@result$p.adjust < 0.05, ]), "\n")
        pdf(file.path(RDIR, "ora_go_bp_dotplot.pdf"), width = 10, height = 8)
        print(dotplot(ego, showCategory = 20, title = "GO Biological Process (Tier 1 DEGs)"))
        dev.off()
        cat("Saved: ora_go_bp_dotplot.pdf\n")
      } else {
        cat("No significant GO terms found.\n")
      }
    }
  }
} else {
  cat("\nSkipping ORA (clusterProfiler not available).\n")
}

cat("\nPathway analysis complete.\n")
