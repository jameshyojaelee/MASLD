#!/usr/bin/env Rscript
# 02c_ora_concordance.R
# ---------------------------------------------------------------------------
# ORA concordance: run enrichGO + enrichKEGG on DEG lists, compute Jaccard
# overlap of significant terms between species
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  library(ggplot2)
})

pdf.options(useDingbats = FALSE)

cat("=== Phase 2c: ORA Concordance ===\n\n")

BASE    <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq"
H_INT   <- file.path(BASE, "Human/Patient_Cohorts/analysis/integration")
ANNOT   <- file.path(H_INT, "results/gene_annotation")
DS_DIR  <- file.path(H_INT, "results/disease_signatures")
INT_DIR <- file.path(H_INT, "results/integration")
MOUSE_PD <- file.path(BASE, "Mouse/Unified_Integration/results/per_diet")
WD      <- file.path(BASE, "Analysis/Cross_Species_Concordance")
RES     <- file.path(WD, "results")

DIETS <- c("MCD", "HFD", "CDAHFD", "FPC")  # 2026-05-29: LIDPAD archived/dropped
h_annot <- fread(file.path(ANNOT, "human_ensg_to_symbol.tsv"))
m_annot <- fread(file.path(ANNOT, "mouse_ensmusg_to_symbol.tsv"))

# ============================================================
#  Helper: run ORA for a DE result
# ============================================================
run_ora <- function(symbols, universe_symbols, org_db, label) {
  # Convert symbols to Entrez IDs
  ids <- tryCatch({
    bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org_db)
  }, error = function(e) data.frame(SYMBOL = character(0), ENTREZID = character(0)))

  univ_ids <- tryCatch({
    bitr(universe_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org_db)
  }, error = function(e) data.frame(SYMBOL = character(0), ENTREZID = character(0)))

  if (nrow(ids) < 10) {
    cat(sprintf("  %s: too few genes mapped (%d), skipping\n", label, nrow(ids)))
    return(NULL)
  }

  # enrichGO
  ego <- tryCatch({
    enrichGO(gene = ids$ENTREZID, universe = univ_ids$ENTREZID,
             OrgDb = org_db, ont = "BP", pAdjustMethod = "BH",
             pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
  }, error = function(e) NULL)

  # enrichKEGG
  org_code <- if (identical(org_db, org.Hs.eg.db)) "hsa" else "mmu"
  ekegg <- tryCatch({
    enrichKEGG(gene = ids$ENTREZID, universe = univ_ids$ENTREZID,
               organism = org_code, pAdjustMethod = "BH",
               pvalueCutoff = 0.05, qvalueCutoff = 0.2)
  }, error = function(e) NULL)

  results <- data.table()
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    go_res <- as.data.table(as.data.frame(ego))
    go_res[, db := "GO_BP"]
    setnames(go_res, "p.adjust", "padj", skip_absent = TRUE)
    results <- rbindlist(list(results, go_res[, .(ID, Description, padj, Count, db)]))
  }
  if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
    kegg_res <- as.data.table(as.data.frame(ekegg))
    kegg_res[, db := "KEGG"]
    setnames(kegg_res, "p.adjust", "padj", skip_absent = TRUE)
    results <- rbindlist(list(results, kegg_res[, .(ID, Description, padj, Count, db)]))
  }

  if (nrow(results) > 0) {
    results[, source := label]
    cat(sprintf("  %s: %d GO_BP + %d KEGG terms\n", label,
      sum(results$db == "GO_BP"), sum(results$db == "KEGG")))
  } else {
    cat(sprintf("  %s: no significant terms\n", label))
  }
  return(results)
}

# ============================================================
#  Run ORA on human signatures
# ============================================================
cat("=== Running ORA on Human Signatures ===\n")

human_ora <- list()
human_files <- list(
  nafl_vs_nash = list(path = file.path(DS_DIR, "nafl_vs_nash_dream.csv"), padj = "adj.P.Val"),
  fibrosis     = list(path = file.path(DS_DIR, "fibrosis_dream.csv"), padj = "adj.P.Val")
)

for (sig_name in names(human_files)) {
  info <- human_files[[sig_name]]
  dt <- fread(info$path)
  if (!"symbol" %in% names(dt)) {
    dt[, gene_base := gsub("\\..*", "", gene)]
    dt <- merge(dt, h_annot[, .(gene_base, symbol)], by = "gene_base", all.x = TRUE)
  }
  sig_symbols <- dt[get(info$padj) < 0.05 & !is.na(symbol), unique(symbol)]
  all_symbols <- dt[!is.na(symbol), unique(symbol)]
  res <- run_ora(sig_symbols, all_symbols, org.Hs.eg.db, sig_name)
  if (!is.null(res)) human_ora[[sig_name]] <- res
}

# ============================================================
#  Run ORA on mouse diets
# ============================================================
cat("\n=== Running ORA on Mouse Diets ===\n")

mouse_ora <- list()
for (diet in DIETS) {
  dt <- fread(file.path(MOUSE_PD, paste0(diet, "_de_results.csv")))
  if (!"symbol" %in% names(dt)) {
    dt[, gene_base := gsub("\\..*", "", gene)]
    dt <- merge(dt, m_annot[, .(gene_base, symbol)], by = "gene_base", all.x = TRUE)
  }
  sig_symbols <- dt[adj.P.Val < 0.05 & !is.na(symbol), unique(symbol)]
  all_symbols <- dt[!is.na(symbol), unique(symbol)]
  res <- run_ora(sig_symbols, all_symbols, org.Mm.eg.db, diet)
  if (!is.null(res)) mouse_ora[[diet]] <- res
}

# ============================================================
#  Concordance: Jaccard overlap of significant terms
# ============================================================
cat("\n=== Computing Term Overlap Concordance ===\n")

all_human_ora <- rbindlist(human_ora, fill = TRUE)
all_mouse_ora <- rbindlist(mouse_ora, fill = TRUE)

fwrite(all_human_ora, file.path(RES, "ora_human_results.csv"))
fwrite(all_mouse_ora, file.path(RES, "ora_mouse_results.csv"))

# Use IDs for cross-species comparison (more reliable than description strings)
# GO IDs (GO:0006955) are species-independent. KEGG IDs have organism prefixes
# (hsa04060 vs mmu04060) that must be stripped for cross-species matching.
strip_kegg_prefix <- function(ids) gsub("^[a-z]+", "", ids)

term_conc <- data.table()
for (sig_name in names(human_ora)) {
  h_raw <- human_ora[[sig_name]]
  h_dt <- if (nrow(h_raw) > 0 && "padj" %in% names(h_raw)) h_raw[padj < 0.05] else data.table(ID=character(), Description=character(), padj=numeric(), Count=integer(), db=character())
  h_ids <- if (nrow(h_dt) > 0) c(
    h_dt[db == "GO_BP", unique(ID)],
    strip_kegg_prefix(h_dt[db == "KEGG", unique(ID)])
  ) else character()
  for (diet in names(mouse_ora)) {
    m_raw <- mouse_ora[[diet]]
    m_dt <- if (nrow(m_raw) > 0 && "padj" %in% names(m_raw)) m_raw[padj < 0.05] else data.table(ID=character(), Description=character(), padj=numeric(), Count=integer(), db=character())
    m_ids <- if (nrow(m_dt) > 0) c(
      m_dt[db == "GO_BP", unique(ID)],
      strip_kegg_prefix(m_dt[db == "KEGG", unique(ID)])
    ) else character()
    overlap <- length(intersect(h_ids, m_ids))
    union_n <- length(union(h_ids, m_ids))
    jacc <- if (union_n > 0) overlap / union_n else 0

    row <- data.table(
      human_signature = sig_name,
      diet = diet,
      n_human_terms = length(h_ids),
      n_mouse_terms = length(m_ids),
      n_overlap = overlap,
      jaccard_terms = round(jacc, 4)
    )
    term_conc <- rbindlist(list(term_conc, row))

    cat(sprintf("  %s × %s: %d/%d human terms, %d overlap (J=%.3f)\n",
      sig_name, diet, length(h_ids), length(m_ids), overlap, jacc))
  }
}

fwrite(term_conc, file.path(RES, "ora_term_concordance.csv"))

# Identify universally enriched terms
if (nrow(all_mouse_ora) > 0) {
  if ("padj" %in% names(all_mouse_ora) && nrow(all_mouse_ora) > 0) {
    all_mouse_ora[, padj_num := as.numeric(padj)]
  } else {
    all_mouse_ora[, padj_num := NA_real_]
  }
  mouse_terms_by_diet <- tryCatch({
    all_mouse_ora[padj_num < 0.05, .(n_diets = uniqueN(source)),
      by = .(Description_lower = tolower(Description))]
  }, error = function(e) {
    cat("  Warning: universally enriched terms calculation failed:", e$message, "\n")
    data.table()
  })
  if (nrow(mouse_terms_by_diet) > 0) {
    universal <- mouse_terms_by_diet[n_diets >= 3]
    cat(sprintf("\nUniversally enriched terms (≥3 diets): %d\n", nrow(universal)))
    if (nrow(universal) > 0) {
      fwrite(universal[order(-n_diets)], file.path(RES, "universally_enriched_terms.csv"))
    }
  }
}

cat("\nSaved: ora_term_concordance.csv\n")
cat("=== Phase 2c complete ===\n")
