#!/usr/bin/env Rscript
# 133_progression_pathway_enrichment.R
# ---------------------------------------------------------------------------
# Pathway enrichment (fgsea) across ALL progression contrasts.
#
# Runs fgsea with t-statistic ranking for three MSigDB collections:
#   - Hallmark (collection = "H")
#   - KEGG (collection = "C2", subcollection = "CP:KEGG_LEGACY")
#   - Reactome (collection = "C2", subcollection = "CP:REACTOME")
#
# Contrasts loaded:
#   C1: MASLD vs Control      (dream_results.csv from Script 05)
#   C2: NAFL vs NASH           (nafl_vs_nash_dream.csv from Script 13)
#   C3: F3-F4 vs F0-F2         (c3_adv_vs_early_fib_dream.csv from Script 130)
#   C4: NAFL vs Control        (c4_nafl_vs_ctrl_dream.csv from Script 130)
#   C5: NAS >= 5 vs NAS < 5    (c5_nas_ge5_vs_lt5_dream.csv from Script 130)
#   C6: Extreme endpoints      (c6_extreme_endpoints_dream.csv from Script 130)
#   C8: Cirrhosis binary       (c8_cirrhosis_dream.csv from Script 130)
#   C9: F2 inflection point    (c9_f2_inflection_dream.csv from Script 130)
#   C7a: Steatosis ordinal     (c7a_steatosis_ordinal_dream.csv from Script 131)
#   C7b: Inflammation ordinal  (c7b_inflammation_ordinal_dream.csv from Script 131)
#   C7c: Ballooning ordinal    (c7c_ballooning_ordinal_dream.csv from Script 131)
#   C11: NASH vs Control       (c11_nash_vs_ctrl_dream.csv from Script 130)
#   C12: Within-NASH fib prog  (c12_early_vs_late_nash_dream.csv from Script 130)
#   C13: NASH vs NAFL fib-adj  (c13_nash_vs_nafl_fib_adj_dream.csv from Script 130)
#   C17: Fibrosis ordinal      (c17_fibrosis_ordinal_dream.csv from Script 131)
#
# Output (to results/progression/):
#   progression_gsea_all.csv         — all fgsea results stacked
#   progression_gsea_nes_matrix.csv  — NES heatmap matrix (pathways x contrasts)
#   progression_gsea_summary.csv     — significant pathway counts per contrast x collection
#
# Usage: Rscript 133_progression_pathway_enrichment.R
# SLURM: cpu, 4 CPU, 32GB RAM, ~45min
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(fgsea)
  library(msigdbr)
  library(data.table)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
ODIR <- file.path(INT, "results/progression")
dir.create(ODIR, recursive = TRUE, showWarnings = FALSE)

cat("=== 133: Progression Pathway Enrichment (fgsea) ===\n")
cat("Started:", as.character(Sys.time()), "\n\n")

# ============================================================
# 1. Load gene set collections from MSigDB
# ============================================================
cat("Loading MSigDB gene sets...\n")

# Hallmark
hallmark_df <- msigdbr(species = "Homo sapiens", collection = "H")
hallmark <- split(hallmark_df$ensembl_gene, hallmark_df$gs_name)
hallmark <- lapply(hallmark, function(x) unique(x[!is.na(x) & x != ""]))
cat(sprintf("  Hallmark pathways: %d\n", length(hallmark)))

# KEGG
kegg_df <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY")
kegg <- split(kegg_df$ensembl_gene, kegg_df$gs_name)
kegg <- lapply(kegg, function(x) unique(x[!is.na(x) & x != ""]))
cat(sprintf("  KEGG pathways: %d\n", length(kegg)))

# Reactome
reactome_df <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME")
reactome <- split(reactome_df$ensembl_gene, reactome_df$gs_name)
reactome <- lapply(reactome, function(x) unique(x[!is.na(x) & x != ""]))
cat(sprintf("  Reactome pathways: %d\n", length(reactome)))

collections <- list(
  Hallmark = hallmark,
  KEGG     = kegg,
  Reactome = reactome
)

# ============================================================
# 2. Define contrast file registry
# ============================================================
# Each entry: contrast_id, file path, t-statistic column name, gene column name
# Note: column names differ across scripts:
#   - Script 05 (dream_results.csv): padj, gene col = "gene"
#   - Script 13 (nafl_vs_nash_dream.csv): adj.P.Val, gene col = "gene"
#   - Script 130 outputs: padj (renamed in script), gene col = "gene"
#   - Script 131 outputs: padj (renamed in script), gene col = "gene"
#   - Transition dream results: adj.P.Val, gene col = "gene"

PROG <- file.path(INT, "results/progression")
SIGS <- file.path(INT, "results/disease_signatures")
INTG <- file.path(INT, "results/integration")

contrast_registry <- list(
  list(id = "C1_MASLD_vs_Ctrl",
       file = file.path(INTG, "dream_results.csv"),
       description = "MASLD vs Control (Script 05)"),
  list(id = "C2_NASH_vs_NAFL",
       file = file.path(SIGS, "nafl_vs_nash_dream.csv"),
       description = "NASH vs NAFL (Script 13; coef=nafl_nashNASH, positive=UP in NASH)"),
  list(id = "C3_Adv_vs_Early_Fib",
       file = file.path(PROG, "c3_adv_vs_early_fib_dream.csv"),
       description = "F3-F4 vs F0-F2 (Script 130)"),
  list(id = "C4_NAFL_vs_Ctrl",
       file = file.path(PROG, "c4_nafl_vs_ctrl_dream.csv"),
       description = "NAFL vs Control (Script 130)"),
  list(id = "C5_NAS_ge5_vs_lt5",
       file = file.path(PROG, "c5_nas_ge5_vs_lt5_dream.csv"),
       description = "NAS >= 5 vs NAS < 5 (Script 130)"),
  list(id = "C6_Extreme_Endpoints",
       file = file.path(PROG, "c6_extreme_endpoints_dream.csv"),
       description = "NASH+F3-F4 vs NAFL+F0-F1 (Script 130)"),
  list(id = "C8_Cirrhosis",
       file = file.path(PROG, "c8_cirrhosis_dream.csv"),
       description = "F4 vs F0-F3 (Script 130)"),
  list(id = "C9_F2_Inflection",
       file = file.path(PROG, "c9_f2_inflection_dream.csv"),
       description = "F2-F4 vs F0-F1 (Script 130)"),
  list(id = "C7a_Steatosis",
       file = file.path(PROG, "c7a_steatosis_ordinal_dream.csv"),
       description = "Steatosis grade ordinal (Script 131)"),
  list(id = "C7b_Inflammation",
       file = file.path(PROG, "c7b_inflammation_ordinal_dream.csv"),
       description = "Lobular inflammation ordinal (Script 131)"),
  list(id = "C7c_Ballooning",
       file = file.path(PROG, "c7c_ballooning_ordinal_dream.csv"),
       description = "Ballooning grade ordinal (Script 131)"),
  list(id = "C11_NASH_vs_Ctrl",
       file = file.path(PROG, "c11_nash_vs_ctrl_dream.csv"),
       description = "NASH vs Control (Script 130)"),
  list(id = "C12_Early_vs_Late_NASH",
       file = file.path(PROG, "c12_early_vs_late_nash_dream.csv"),
       description = "Within-NASH fibrosis progression (Script 130)"),
  list(id = "C13_NASH_vs_NAFL_FibAdj",
       file = file.path(PROG, "c13_nash_vs_nafl_fib_adj_dream.csv"),
       description = "NASH vs NAFL fibrosis-adjusted (Script 130)"),
  list(id = "C17_Fibrosis_Ordinal",
       file = file.path(PROG, "c17_fibrosis_ordinal_dream.csv"),
       description = "Fibrosis dose-response ordinal (Script 131)")
)

# ============================================================
# 3. Load contrast results and run fgsea
# ============================================================

# Helper: load dream results, normalize column names, return data.table with
# gene_base (Ensembl ID without version) and t-statistic
load_dream <- function(entry) {
  if (!file.exists(entry$file)) {
    cat(sprintf("  WARNING: %s not found — SKIPPED\n", basename(entry$file)))
    return(NULL)
  }
  dt <- fread(entry$file)

  # Normalize column names: ensure we have "t" and "gene"
  if (!"t" %in% names(dt)) {
    cat(sprintf("  WARNING: no 't' column in %s — SKIPPED\n", basename(entry$file)))
    return(NULL)
  }

  # Strip Ensembl version numbers for MSigDB matching
  dt[, gene_base := gsub("\\..*", "", gene)]
  dt
}

# Run fgsea for one contrast across all collections
run_fgsea_all_collections <- function(dt, contrast_id, collections) {
  # Rank by t-statistic
  dt <- dt[order(-t)]
  ranks <- dt$t
  names(ranks) <- dt$gene_base

  # Remove duplicates (keep first = highest |t|)
  ranks <- ranks[!duplicated(names(ranks))]

  # Remove NA/NaN/Inf
  ranks <- ranks[is.finite(ranks)]

  results <- list()
  for (coll_name in names(collections)) {
    pathways <- collections[[coll_name]]

    res <- fgsea(
      pathways    = pathways,
      stats       = ranks,
      minSize     = 15,
      maxSize     = 500,
      nPermSimple = 10000
    )
    res$contrast_id <- contrast_id
    res$collection  <- coll_name

    # Convert leadingEdge list to semicolon-separated string
    res[, leadingEdge := sapply(leadingEdge, paste, collapse = ";")]

    sig <- res[padj < 0.05]
    cat(sprintf("    %s: %d sig / %d total (padj < 0.05)\n",
      coll_name, nrow(sig), nrow(res)))

    results[[coll_name]] <- res
  }

  rbindlist(results, fill = TRUE)
}

# --- Main loop ---
all_gsea <- list()
n_loaded <- 0
n_skipped <- 0

for (entry in contrast_registry) {
  cat(sprintf("\n--- %s: %s ---\n", entry$id, entry$description))

  dt <- load_dream(entry)
  if (is.null(dt)) {
    n_skipped <- n_skipped + 1
    next
  }

  n_genes <- uniqueN(dt$gene_base)
  cat(sprintf("  Loaded: %d rows, %d unique genes\n", nrow(dt), n_genes))

  res <- run_fgsea_all_collections(dt, entry$id, collections)
  all_gsea[[entry$id]] <- res
  n_loaded <- n_loaded + 1
}

cat(sprintf("\n\nContrasts loaded: %d, skipped: %d\n", n_loaded, n_skipped))

if (n_loaded == 0) {
  cat("ERROR: No contrasts loaded — no upstream dream results found.\n")
  cat("Run Scripts 05, 13, 130, and 131 first.\n")
  quit(status = 1)
}

# ============================================================
# 4. Combine and save all results
# ============================================================
gsea_all <- rbindlist(all_gsea, fill = TRUE)

# Reorder columns for clarity
col_order <- c("pathway", "pval", "padj", "NES", "size", "leadingEdge",
               "contrast_id", "collection", "ES", "log2err")
col_order <- intersect(col_order, names(gsea_all))
extra_cols <- setdiff(names(gsea_all), col_order)
setcolorder(gsea_all, c(col_order, extra_cols))

fwrite(gsea_all, file.path(ODIR, "progression_gsea_all.csv"))
cat(sprintf("\nSaved: progression_gsea_all.csv (%d rows)\n", nrow(gsea_all)))

# ============================================================
# 5. NES heatmap matrix: pathways (rows) x contrasts (columns)
# ============================================================
# Use all pathway-contrast combinations (fill missing with NA)
nes_wide <- dcast(gsea_all, pathway + collection ~ contrast_id, value.var = "NES")
fwrite(nes_wide, file.path(ODIR, "progression_gsea_nes_matrix.csv"))
cat(sprintf("Saved: progression_gsea_nes_matrix.csv (%d pathways x %d contrasts)\n",
  nrow(nes_wide), n_loaded))

# ============================================================
# 6. Summary: significant pathway counts per contrast x collection
# ============================================================
summary_dt <- gsea_all[, .(
  n_sig_005    = sum(padj < 0.05, na.rm = TRUE),
  n_sig_01     = sum(padj < 0.10, na.rm = TRUE),
  n_total      = .N,
  n_up         = sum(padj < 0.05 & NES > 0, na.rm = TRUE),
  n_down       = sum(padj < 0.05 & NES < 0, na.rm = TRUE),
  top_pathway  = pathway[which.min(padj)],
  top_NES      = NES[which.min(padj)],
  top_padj     = min(padj, na.rm = TRUE)
), by = .(contrast_id, collection)]

setorder(summary_dt, contrast_id, collection)
fwrite(summary_dt, file.path(ODIR, "progression_gsea_summary.csv"))
cat(sprintf("Saved: progression_gsea_summary.csv (%d rows)\n", nrow(summary_dt)))

# --- Print summary table ---
cat("\n")
cat(strrep("=", 80), "\n")
cat("  PROGRESSION GSEA SUMMARY\n")
cat(strrep("=", 80), "\n\n")

for (cid in unique(summary_dt$contrast_id)) {
  sub <- summary_dt[contrast_id == cid]
  total_sig <- sum(sub$n_sig_005)
  cat(sprintf("  %s: %d sig pathways (padj < 0.05)\n", cid, total_sig))
  for (i in seq_len(nrow(sub))) {
    cat(sprintf("    %s: %d up, %d down (top: %s, NES=%.2f, padj=%.2e)\n",
      sub$collection[i], sub$n_up[i], sub$n_down[i],
      sub$top_pathway[i], sub$top_NES[i], sub$top_padj[i]))
  }
}

cat(sprintf("\n=== Script 133 completed: %s ===\n", as.character(Sys.time())))
