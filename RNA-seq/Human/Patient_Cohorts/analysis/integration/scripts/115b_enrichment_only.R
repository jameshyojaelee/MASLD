#!/usr/bin/env Rscript
# 115b_enrichment_only.R — Recovery script
# Picks up from 115's saved intermediate results (transition_fib/nas_dream_results.csv,
# transition_tau_index.csv) and does the merge + pathway enrichment step.
# Avoids re-running the 3-hour dream contrasts.

suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
  library(msigdbr)
})

cat("=== 115b: Recovery — Merge + Pathway Enrichment ===\n")

INTEG <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration"
OUTDIR <- file.path(INTEG, "results/progression")

# ── Load saved intermediate results ───────────────────────────────────────
fib_dt <- fread(file.path(OUTDIR, "transition_fib_dream_results.csv"))
nas_dt <- fread(file.path(OUTDIR, "transition_nas_dream_results.csv"))
all_tau <- fread(file.path(OUTDIR, "transition_tau_index.csv"))

cat(sprintf("  Fib dream: %d rows\n", nrow(fib_dt)))
cat(sprintf("  NAS dream: %d rows\n", nrow(nas_dt)))
cat(sprintf("  Tau index: %d rows\n", nrow(all_tau)))

# ── Load Ensembl → gene symbol mapping (from Script 17 cache) ────────────
ANNOT_FILE <- file.path(INTEG, "results/gene_annotation/human_ensg_to_symbol.tsv")
stopifnot(file.exists(ANNOT_FILE))
annot <- fread(ANNOT_FILE, select = c("gene_id", "symbol"))
cat(sprintf("  Gene annotation: %d Ensembl IDs loaded\n", nrow(annot)))

# Map Ensembl IDs to gene symbols in dream results
fib_dt[annot, gene_symbol := i.symbol, on = .(gene = gene_id)]
nas_dt[annot, gene_symbol := i.symbol, on = .(gene = gene_id)]
all_tau[annot, gene_symbol := i.symbol, on = .(gene = gene_id)]

fib_mapped <- sum(!is.na(fib_dt$gene_symbol))
nas_mapped <- sum(!is.na(nas_dt$gene_symbol))
cat(sprintf("  Fib mapped: %d / %d (%.1f%%)\n", fib_mapped, nrow(fib_dt), 100 * fib_mapped / nrow(fib_dt)))
cat(sprintf("  NAS mapped: %d / %d (%.1f%%)\n", nas_mapped, nrow(nas_dt), 100 * nas_mapped / nrow(nas_dt)))

# ── Merge into transition programs ────────────────────────────────────────
all_dt <- rbind(fib_dt, nas_dt, fill = TRUE)
tau_cols <- intersect(c("gene", "gene_symbol", "tau", "peak_transition", "max_abs_t", "peak_logFC", "peak_padj"),
                      names(all_tau))

# Ensure tau has unique genes (take first if duplicates)
tau_unique <- all_tau[!duplicated(gene), ..tau_cols]
# Drop gene_symbol from tau_unique to avoid duplication on merge (all_dt already has it)
tau_merge_cols <- setdiff(tau_cols, "gene_symbol")
tp <- merge(all_dt, tau_unique[, ..tau_merge_cols], by = "gene", all.x = TRUE, allow.cartesian = FALSE)
fwrite(tp, file.path(OUTDIR, "transition_programs.csv"))
cat(sprintf("  Saved transition_programs.csv (%d rows, gene_symbol mapped: %d)\n",
            nrow(tp), sum(!is.na(tp$gene_symbol))))

# ── Pathway enrichment per transition ─────────────────────────────────────
cat("\n=== Pathway Enrichment ===\n")

hallmark <- msigdbr(species = "Homo sapiens", collection = "H")
hallmark_list <- split(hallmark$gene_symbol, hallmark$gs_name)

kegg <- tryCatch(
  msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_MEDICUS"),
  error = function(e) {
    tryCatch(
      msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY"),
      error = function(e2) {
        cat("  WARNING: KEGG not available, using Hallmark only\n")
        NULL
      }
    )
  }
)

if (!is.null(kegg)) {
  kegg_list <- split(kegg$gene_symbol, kegg$gs_name)
  pathway_db <- c(hallmark_list, kegg_list)
} else {
  pathway_db <- hallmark_list
}
cat(sprintf("  Pathway sets: %d\n", length(pathway_db)))

# Run fgsea per transition
transitions <- unique(tp$transition)
all_enrich <- list()

for (trans in transitions) {
  sub <- copy(tp[transition == trans])

  # Build ranked gene list (t-statistic)
  t_col <- grep("^t$|^tstat|^t\\..*", names(sub), value = TRUE)
  if (length(t_col) == 0) {
    # Use logFC / SE approximation
    if ("logFC" %in% names(sub) && "P.Value" %in% names(sub)) {
      sub[, .rank_stat := sign(logFC) * -log10(pmax(P.Value, 1e-300))]
    } else if ("logFC" %in% names(sub)) {
      sub[, .rank_stat := logFC]
    } else {
      cat(sprintf("  %s: no ranking stat available, skipping\n", trans))
      next
    }
  } else {
    sub[, .rank_stat := get(t_col[1])]
  }

  # Use gene symbols for pathway matching; drop unmapped genes and NAs
  sub_sym <- sub[!is.na(gene_symbol) & gene_symbol != "" & !is.na(.rank_stat)]

  # Deduplicate symbols: keep the entry with the highest |stat| per symbol
  # (multiple Ensembl IDs can map to the same symbol, e.g. Y_RNA, U6)
  sub_sym[, .abs_stat := abs(.rank_stat)]
  sub_dedup <- sub_sym[sub_sym[, .I[which.max(.abs_stat)], by = gene_symbol]$V1]

  stats <- setNames(sub_dedup$.rank_stat, sub_dedup$gene_symbol)
  stats <- sort(stats, decreasing = TRUE)

  if (length(stats) < 100) {
    cat(sprintf("  %s: too few genes (%d), skipping\n", trans, length(stats)))
    next
  }

  res <- fgsea(pathways = pathway_db, stats = stats, minSize = 15, maxSize = 500)
  res$transition <- trans
  n_sig <- sum(res$padj < 0.05)
  cat(sprintf("  %s: %d pathways tested, %d significant (padj<0.05), %d genes in ranked list\n",
              trans, nrow(res), n_sig, length(stats)))

  all_enrich[[trans]] <- res
}

if (length(all_enrich) > 0) {
  enrich_dt <- rbindlist(all_enrich, fill = TRUE)
  # Remove leadingEdge column (list column) for CSV
  enrich_dt[, leadingEdge := sapply(leadingEdge, function(x) paste(head(x, 10), collapse = ";"))]
  fwrite(enrich_dt, file.path(OUTDIR, "transition_pathway_enrichment.csv"))
  cat(sprintf("\n  Saved transition_pathway_enrichment.csv (%d rows)\n", nrow(enrich_dt)))
} else {
  cat("\n  No enrichment results to save\n")
}

cat(sprintf("\n=== 115b: COMPLETE (%s) ===\n", Sys.time()))
