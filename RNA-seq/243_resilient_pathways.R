#!/usr/bin/env Rscript
# 243_resilient_pathways.R
#
# Healthy-control characterization pipeline, step 4 of 6.
#
# Pathway enrichment for the 3-way resilience DE results from Script 242.
# fgsea on MSigDB Hallmark + Reactome collections, ranked by t-stat.
#
# Inputs:
#   resilient_de_pairwise_GSE126848.csv  (within-cohort, primary)
#   resilient_de_pairwise_crosscohort.csv (cross-cohort, sensitivity)
# Output:
#   resilient_gsea.csv (NES, padj, leading_edge per pathway × contrast × analysis)
#
# Spec: docs/superpowers/specs/2026-04-27-healthy-control-audit-design.md
# Env: rnaseq

suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
  library(msigdbr)
})

BASE  <- Sys.getenv("MASLD_PROJECT_ROOT",
                    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
HCDIR <- file.path(BASE, "RNA-seq/results/audit_sensitivity/healthy_control_audit")
DEDIR <- file.path(HCDIR, "resilience_de")

# Gene-symbol map from atlas (ensembl -> human_symbol)
ATLAS <- file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
atlas <- fread(ATLAS, select = c("ensembl_id", "human_symbol"))
atlas[, ensembl_clean := sub("\\..*", "", ensembl_id)]
gene_map <- unique(atlas[, .(ensembl_clean, symbol = human_symbol)])

# MSigDB collections (Hallmark + selected Reactome)
cat("Loading MSigDB collections...\n")
hallmark <- as.data.table(msigdbr(species = "Homo sapiens", collection = "H"))
reactome <- as.data.table(msigdbr(species = "Homo sapiens",
                                   collection = "C2", subcollection = "CP:REACTOME"))
hallmark[, gs_collection := "H_HALLMARK"]
reactome[, gs_collection := "C2_REACTOME"]
gs_long <- rbindlist(list(hallmark, reactome), use.names = TRUE, fill = TRUE)

# Build pathway list keyed by gene_symbol
pathways <- split(gs_long$gene_symbol, gs_long$gs_name)
pathways <- lapply(pathways, unique)
cat(sprintf("  Pathways: Hallmark=%d, Reactome=%d, total=%d\n",
            length(unique(hallmark$gs_name)),
            length(unique(reactome$gs_name)),
            length(pathways)))

# Run fgsea on a DE contrast results table
run_fgsea <- function(de_dt, contrast_name, label, min_size = 10, max_size = 1000) {
  d <- de_dt[contrast == contrast_name]
  if (nrow(d) == 0) {
    cat(sprintf("  [skip] No rows for contrast '%s' in %s\n", contrast_name, label))
    return(NULL)
  }
  # Map ensembl gene_id -> symbol via atlas
  d[, ensembl_clean := sub("\\..*", "", gene_id)]
  d <- merge(d, gene_map, by = "ensembl_clean", all.x = TRUE)
  d <- d[!is.na(symbol) & symbol != ""]
  d <- d[order(-abs(t))][!duplicated(symbol)]
  ranks <- d$t
  names(ranks) <- d$symbol

  cat(sprintf("  [%s | %s] %d ranked genes; running fgsea...\n",
              label, contrast_name, length(ranks)))
  set.seed(42)
  fg <- fgsea(pathways = pathways, stats = ranks,
              minSize = min_size, maxSize = max_size,
              eps = 0)
  fg <- as.data.table(fg)
  fg[, contrast := contrast_name]
  fg[, analysis := label]
  # Compress leading_edge (originally a list column) to comma-separated
  fg[, leading_edge_str := vapply(leadingEdge,
                                    function(x) paste(head(x, 50), collapse = ";"),
                                    character(1))]
  fg[, leadingEdge := NULL]
  setnames(fg, "leading_edge_str", "leading_edge")
  fg
}

# === Run for each available DE result ===
all_gsea <- list()

# Analysis A: within-GSE126848
de_A <- file.path(DEDIR, "resilient_de_pairwise_GSE126848.csv")
if (file.exists(de_A)) {
  cat(sprintf("\n=== Analysis A: within-GSE126848 ===\n"))
  d_A <- fread(de_A)
  for (cn in unique(d_A$contrast)) {
    fg <- run_fgsea(d_A, cn, "GSE126848_within")
    if (!is.null(fg)) all_gsea[[length(all_gsea) + 1]] <- fg
  }
} else {
  cat("  [warn] Script 242 GSE126848 output not yet present. Run sbatch first.\n")
}

# Analysis B: cross-cohort
de_B <- file.path(DEDIR, "resilient_de_pairwise_crosscohort.csv")
if (file.exists(de_B)) {
  cat(sprintf("\n=== Analysis B: cross-cohort ===\n"))
  d_B <- fread(de_B)
  for (cn in unique(d_B$contrast)) {
    fg <- run_fgsea(d_B, cn, "cross_cohort")
    if (!is.null(fg)) all_gsea[[length(all_gsea) + 1]] <- fg
  }
} else {
  cat("  [warn] Script 242 cross-cohort output not yet present. Run sbatch first.\n")
}

if (length(all_gsea) == 0) {
  cat("\nNo DE files found; exiting (no GSEA computed).\n")
  quit(save = "no", status = 1)
}

gsea_all <- rbindlist(all_gsea, use.names = TRUE, fill = TRUE)
fwrite(gsea_all, file.path(DEDIR, "resilient_gsea.csv"))
cat(sprintf("\nWrote: resilient_gsea.csv (%d rows)\n", nrow(gsea_all)))

# Summary: top pathways per contrast (FDR<0.05)
cat("\nTop pathways by contrast (FDR<0.05, sorted by NES):\n")
for (a in unique(gsea_all$analysis)) {
  for (cn in unique(gsea_all$contrast[gsea_all$analysis == a])) {
    sig <- gsea_all[analysis == a & contrast == cn & padj < 0.05][order(-abs(NES))][1:10]
    cat(sprintf("\n[%s | %s] top 10 by |NES|:\n", a, cn))
    if (nrow(sig) > 0 && !is.na(sig$pathway[1])) {
      print(sig[, .(pathway, NES, padj, size)])
    } else {
      cat("  (none significant)\n")
    }
  }
}

cat("\nDone.\n")
