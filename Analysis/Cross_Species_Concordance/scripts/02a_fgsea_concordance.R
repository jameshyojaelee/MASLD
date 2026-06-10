#!/usr/bin/env Rscript
# 02a_fgsea_concordance.R
# ---------------------------------------------------------------------------
# Pathway-level concordance via FGSEA: rank genes by signed significance,
# test enrichment in MSigDB Hallmark + KEGG + Reactome, compare NES across species
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
  library(msigdbr)
  library(ggplot2)
})

pdf.options(useDingbats = FALSE)

cat("=== Phase 2a: FGSEA Pathway Concordance ===\n\n")

BASE    <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq"
H_INT   <- file.path(BASE, "Human/Patient_Cohorts/analysis/integration")
ANNOT   <- file.path(H_INT, "results/gene_annotation")
DS_DIR  <- file.path(H_INT, "results/disease_signatures")
INT_DIR <- file.path(H_INT, "results/integration")
MOUSE_PD <- file.path(BASE, "Mouse/Unified_Integration/results/per_diet")
WD      <- file.path(BASE, "Analysis/Cross_Species_Concordance")
RES     <- file.path(WD, "results")
PLOTS   <- file.path(WD, "plots")

DIETS <- c("MCD", "HFD", "CDAHFD", "FPC")  # 2026-05-29: LIDPAD archived/dropped
ortho <- fread(file.path(ANNOT, "ortholog_mapping.tsv"))

# ============================================================
#  Gene sets: Hallmark + KEGG + Reactome (human symbols)
# ============================================================
cat("Loading MSigDB gene sets...\n")
hallmark <- as.data.table(msigdbr(species = "Homo sapiens", collection = "H"))
kegg     <- as.data.table(msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY"))
reactome <- as.data.table(msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME"))

all_sets <- rbindlist(list(
  hallmark[, .(gs_name, gene_symbol, collection = "Hallmark")],
  kegg[, .(gs_name, gene_symbol, collection = "KEGG")],
  reactome[, .(gs_name, gene_symbol, collection = "Reactome")]
))

pathways <- split(all_sets$gene_symbol, all_sets$gs_name)
cat(sprintf("  %d gene sets (%d Hallmark, %d KEGG, %d Reactome)\n",
  length(pathways),
  length(unique(hallmark$gs_name)),
  length(unique(kegg$gs_name)),
  length(unique(reactome$gs_name))))

# Mouse gene sets (using ortholog-mapped symbols)
mouse_hallmark  <- as.data.table(msigdbr(species = "Mus musculus", collection = "H"))
mouse_kegg      <- as.data.table(msigdbr(species = "Mus musculus", collection = "C2", subcollection = "CP:KEGG_LEGACY"))
mouse_reactome  <- as.data.table(msigdbr(species = "Mus musculus", collection = "C2", subcollection = "CP:REACTOME"))

mouse_sets <- rbindlist(list(
  mouse_hallmark[, .(gs_name, gene_symbol, collection = "Hallmark")],
  mouse_kegg[, .(gs_name, gene_symbol, collection = "KEGG")],
  mouse_reactome[, .(gs_name, gene_symbol, collection = "Reactome")]
))
mouse_pathways <- split(mouse_sets$gene_symbol, mouse_sets$gs_name)

# ============================================================
#  Helper: create ranked gene list for FGSEA
# ============================================================
make_ranks <- function(dt, lfc_col = "logFC", pval_col, symbol_col = "symbol") {
  dt <- copy(dt)
  # Resolve the p-value column robustly: C2/LVQW outputs use `padj`, dream/limma
  # outputs use `adj.P.Val`. Fall back across common names if the requested one is
  # absent (prevents the chain aborting when a signature file's columns have drifted).
  if (!(pval_col %in% names(dt))) {
    alt <- intersect(c("padj", "adj.P.Val", "FDR", "P.Value", "pvalue", "pval", "P.value"), names(dt))
    if (length(alt) == 0)
      stop(sprintf("make_ranks: no usable p-value column; have: %s", paste(names(dt), collapse = ", ")))
    pval_col <- alt[1]
  }
  dt <- dt[!is.na(get(symbol_col)) & get(symbol_col) != ""]
  dt[, rank_score := sign(get(lfc_col)) * -log10(pmax(get(pval_col), 1e-300))]
  ranks <- dt[, .(score = mean(rank_score)), by = symbol_col]
  r <- setNames(ranks$score, ranks[[symbol_col]])
  sort(r, decreasing = TRUE)
}

# ============================================================
#  Run FGSEA on human signatures
# ============================================================
cat("\n=== Running FGSEA on Human Signatures ===\n")

human_files <- list(
  disease_vs_ctrl = list(path = file.path(INT_DIR, "canonical_deg_results.csv"), padj = "padj"),
  nafl_specific   = list(path = file.path(RES, "nafl_vs_ctrl_dream.csv"), padj = "adj.P.Val"),
  nafl_vs_nash    = list(path = file.path(DS_DIR, "nafl_vs_nash_dream.csv"), padj = "adj.P.Val"),
  fibrosis        = list(path = file.path(DS_DIR, "fibrosis_dream.csv"), padj = "adj.P.Val")
)

human_fgsea <- list()
for (sig_name in names(human_files)) {
  info <- human_files[[sig_name]]
  dt <- fread(info$path)
  if (!"symbol" %in% names(dt)) {
    annot_h <- fread(file.path(ANNOT, "human_ensg_to_symbol.tsv"))
    dt[, gene_base := gsub("\\..*", "", gene)]
    dt <- merge(dt, annot_h[, .(gene_base, symbol)], by = "gene_base", all.x = TRUE)
  }
  ranks <- make_ranks(dt, "logFC", info$padj, "symbol")
  cat(sprintf("  %s: %d ranked genes\n", sig_name, length(ranks)))

  res <- fgsea(pathways = pathways, stats = ranks, minSize = 15, maxSize = 500, nPermSimple = 10000)
  res <- as.data.table(res)
  res[, source := sig_name]
  human_fgsea[[sig_name]] <- res
}

human_fgsea_all <- rbindlist(human_fgsea, fill = TRUE)

# ============================================================
#  Run FGSEA on mouse diets
# ============================================================
cat("\n=== Running FGSEA on Mouse Diets ===\n")

mouse_fgsea <- list()
for (diet in DIETS) {
  dt <- fread(file.path(MOUSE_PD, paste0(diet, "_de_results.csv")))
  # Get mouse symbols
  if (!"symbol" %in% names(dt)) {
    m_annot <- fread(file.path(ANNOT, "mouse_ensmusg_to_symbol.tsv"))
    dt[, gene_base := gsub("\\..*", "", gene)]
    dt <- merge(dt, m_annot[, .(gene_base, symbol)], by = "gene_base", all.x = TRUE)
  }
  ranks <- make_ranks(dt, "logFC", "adj.P.Val", "symbol")
  cat(sprintf("  %s: %d ranked genes\n", diet, length(ranks)))

  res <- fgsea(pathways = mouse_pathways, stats = ranks, minSize = 15, maxSize = 500, nPermSimple = 10000)
  res <- as.data.table(res)
  res[, source := diet]
  mouse_fgsea[[diet]] <- res
}

mouse_fgsea_all <- rbindlist(mouse_fgsea, fill = TRUE)

# ============================================================
#  Concordance: compare NES across species
# ============================================================
cat("\n=== Computing Pathway-Level Concordance ===\n")

# Find pathways present in both
common_pathways <- intersect(
  human_fgsea_all[, unique(pathway)],
  mouse_fgsea_all[, unique(pathway)]
)
cat("Common pathways:", length(common_pathways), "\n")

pathway_concordance <- data.table()
for (sig_name in names(human_files)) {
  h_nes <- human_fgsea_all[source == sig_name & pathway %in% common_pathways,
    .(pathway, h_NES = NES, h_padj = padj)]

  for (diet in DIETS) {
    m_nes <- mouse_fgsea_all[source == diet & pathway %in% common_pathways,
      .(pathway, m_NES = NES, m_padj = padj)]

    paired <- merge(h_nes, m_nes, by = "pathway")
    if (nrow(paired) < 10) next

    rho <- cor(paired$h_NES, paired$m_NES, method = "spearman", use = "complete.obs")
    both_sig <- paired[h_padj < 0.05 & m_padj < 0.05]
    conc <- if (nrow(both_sig) > 0) sum(sign(both_sig$h_NES) == sign(both_sig$m_NES)) else 0
    disc <- if (nrow(both_sig) > 0) sum(sign(both_sig$h_NES) != sign(both_sig$m_NES)) else 0

    row <- data.table(
      human_signature = sig_name,
      diet = diet,
      n_pathways = nrow(paired),
      rho_NES = round(rho, 4),
      n_both_sig = nrow(both_sig),
      n_concordant = conc,
      n_discordant = disc,
      concordance_pct = round(conc / max(conc + disc, 1) * 100, 1)
    )
    pathway_concordance <- rbindlist(list(pathway_concordance, row))

    cat(sprintf("  %s × %s: ρ_NES=%.3f, pathway_conc=%d/%d (%.1f%%)\n",
      sig_name, diet, rho, conc, conc + disc,
      conc / max(conc + disc, 1) * 100))
  }
}

fwrite(pathway_concordance, file.path(RES, "fgsea_pathway_concordance.csv"))

# Save full FGSEA results
fwrite(human_fgsea_all[, .(pathway, padj, NES, size, source)],
  file.path(RES, "fgsea_human_results.csv"))
fwrite(mouse_fgsea_all[, .(pathway, padj, NES, size, source)],
  file.path(RES, "fgsea_mouse_results.csv"))

cat("\nSaved: fgsea_pathway_concordance.csv, fgsea_human/mouse_results.csv\n")
cat("=== Phase 2a complete ===\n")
