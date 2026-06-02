#!/usr/bin/env Rscript
# Compare old library (Govaere/Hoang + MCD mouse) vs Dream DEGs at various LFC cutoffs
# Library uses mouse Ensembl IDs with human ortholog mappings
# Dream uses versioned human Ensembl IDs

library(data.table)

cat("=== Library vs Dream DEG Overlap Analysis ===\n\n")

# --- 1. Read library ---
lib_core <- fread("results/library/final_core_degs.csv")
lib_lnc  <- fread("results/library/final_lncrna_degs.csv")

cat(sprintf("Library: %d protein-coding + %d lncRNA = %d total targets\n",
            nrow(lib_core), nrow(lib_lnc), nrow(lib_core) + nrow(lib_lnc)))

# Extract all unique human Ensembl IDs from library (split multi-orthologs)
lib_human_ids <- unique(unlist(strsplit(lib_core$human_ortholog_ids, ";")))
lib_human_ids <- lib_human_ids[lib_human_ids != "" & !is.na(lib_human_ids)]

# Also extract human gene symbols for display
lib_symbol_map <- rbindlist(lapply(seq_len(nrow(lib_core)), function(i) {
  ids <- unlist(strsplit(lib_core$human_ortholog_ids[i], ";"))
  syms <- unlist(strsplit(lib_core$human_ortholog_symbols[i], ";"))
  if (length(ids) == 0 || all(ids == "")) return(NULL)
  data.table(ensembl_id = ids, symbol = syms)
}))
lib_symbol_map <- unique(lib_symbol_map)

cat(sprintf("Library maps to %d unique human Ensembl IDs (%d unique gene symbols)\n",
            length(lib_human_ids), length(unique(lib_symbol_map$symbol))))

# How many library entries have NO human ortholog?
no_ortho <- sum(lib_core$n_human_orthologs == 0)
cat(sprintf("Library entries with no human ortholog (cannot compare): %d\n", no_ortho))
cat(sprintf("lncRNA library: all %d have 0 human orthologs (mouse-only, not comparable)\n\n",
            nrow(lib_lnc)))

# --- 2. Read dream results ---
dream <- fread("RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/dream_results.csv")
cat(sprintf("Dream results: %d genes total\n", nrow(dream)))

# Strip version from Ensembl IDs
dream[, ensembl_id := sub("\\.[0-9]+$", "", gene)]

# Quick stats
cat(sprintf("Dream DEGs (padj < 0.1): %d\n", sum(dream$padj < 0.1, na.rm = TRUE)))
cat(sprintf("Dream DEGs (padj < 0.05): %d\n\n", sum(dream$padj < 0.05, na.rm = TRUE)))

# --- 3. LFC distribution in dream for library genes ---
dream_lib <- dream[ensembl_id %in% lib_human_ids]
cat(sprintf("Library human orthologs found in dream results: %d / %d (%.1f%%)\n",
            nrow(dream_lib), length(lib_human_ids),
            100 * nrow(dream_lib) / length(lib_human_ids)))

# Add symbols
dream_lib <- merge(dream_lib, lib_symbol_map, by = "ensembl_id", all.x = TRUE)

cat(sprintf("\nLFC distribution of library genes in dream:\n"))
cat(sprintf("  Min: %.3f, Q1: %.3f, Median: %.3f, Q3: %.3f, Max: %.3f\n",
            min(dream_lib$logFC), quantile(dream_lib$logFC, 0.25),
            median(dream_lib$logFC), quantile(dream_lib$logFC, 0.75),
            max(dream_lib$logFC)))
cat(sprintf("  Mean |logFC|: %.3f\n", mean(abs(dream_lib$logFC))))
cat(sprintf("  Fraction with padj < 0.1: %.1f%%\n",
            100 * mean(dream_lib$padj < 0.1, na.rm = TRUE)))
cat(sprintf("  Fraction with padj < 0.05: %.1f%%\n\n",
            100 * mean(dream_lib$padj < 0.05, na.rm = TRUE)))

# --- 4. Overlap at various LFC cutoffs ---
cutoffs <- c(0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0)

cat("=== Overlap at Different |logFC| Cutoffs (padj < 0.1) ===\n\n")
cat(sprintf("%-10s %10s %10s %10s %10s %10s\n",
            "|LFC| >=", "Dream_DEGs", "Overlap", "Pct_of_Lib", "Pct_of_Dream", "Jaccard"))

results <- list()
for (lfc_cut in cutoffs) {
  # Dream DEGs at this cutoff
  dream_degs <- dream[padj < 0.1 & abs(logFC) >= lfc_cut, ensembl_id]

  # Overlap
  overlap <- intersect(lib_human_ids, dream_degs)
  lib_only <- setdiff(lib_human_ids, dream_degs)
  dream_only <- setdiff(dream_degs, lib_human_ids)

  n_overlap <- length(overlap)
  n_dream <- length(dream_degs)
  pct_lib <- 100 * n_overlap / length(lib_human_ids)
  pct_dream <- ifelse(n_dream > 0, 100 * n_overlap / n_dream, 0)
  jaccard <- n_overlap / (length(lib_human_ids) + n_dream - n_overlap)

  cat(sprintf("%-10s %10d %10d %9.1f%% %11.1f%% %9.3f\n",
              sprintf("%.2f", lfc_cut), n_dream, n_overlap, pct_lib, pct_dream, jaccard))

  results[[as.character(lfc_cut)]] <- data.table(
    lfc_cutoff = lfc_cut,
    dream_degs = n_dream,
    overlap = n_overlap,
    lib_only = length(lib_only),
    dream_only = length(dream_only),
    pct_of_library = round(pct_lib, 1),
    pct_of_dream = round(pct_dream, 1),
    jaccard = round(jaccard, 3)
  )
}

# --- 5. Also do UP-only comparison (library was upregulated targets) ---
cat("\n=== Overlap at Different logFC Cutoffs (UP only, padj < 0.1) ===\n")
cat("(Library targeted upregulated genes)\n\n")
cat(sprintf("%-10s %10s %10s %10s %10s %10s\n",
            "LFC >=", "Dream_UP", "Overlap", "Pct_of_Lib", "Pct_of_Dream", "Jaccard"))

for (lfc_cut in cutoffs) {
  dream_up <- dream[padj < 0.1 & logFC >= lfc_cut, ensembl_id]
  overlap <- intersect(lib_human_ids, dream_up)
  n_overlap <- length(overlap)
  n_dream <- length(dream_up)
  pct_lib <- 100 * n_overlap / length(lib_human_ids)
  pct_dream <- ifelse(n_dream > 0, 100 * n_overlap / n_dream, 0)
  jaccard <- n_overlap / (length(lib_human_ids) + n_dream - n_overlap)

  cat(sprintf("%-10s %10d %10d %9.1f%% %11.1f%% %9.3f\n",
              sprintf("%.2f", lfc_cut), n_dream, n_overlap, pct_lib, pct_dream, jaccard))
}

# --- 6. Direction check ---
cat("\n=== Direction Concordance ===\n")
# Library was all upregulated targets — check dream direction
dream_lib_sig <- dream_lib[padj < 0.1]
n_up <- sum(dream_lib_sig$logFC > 0)
n_down <- sum(dream_lib_sig$logFC < 0)
cat(sprintf("Library genes significant in dream (padj<0.1): %d\n", nrow(dream_lib_sig)))
cat(sprintf("  Upregulated in dream (concordant): %d (%.1f%%)\n",
            n_up, 100 * n_up / nrow(dream_lib_sig)))
cat(sprintf("  Downregulated in dream (discordant): %d (%.1f%%)\n\n",
            n_down, 100 * n_down / nrow(dream_lib_sig)))

# --- 7. Top library genes by dream effect size ---
cat("=== Top 30 Library Genes by Dream |logFC| (padj < 0.1) ===\n\n")
top_genes <- dream_lib_sig[order(-abs(logFC))][1:min(30, nrow(dream_lib_sig))]
cat(sprintf("%-15s %-18s %8s %12s\n", "Symbol", "Ensembl", "logFC", "padj"))
for (i in seq_len(nrow(top_genes))) {
  cat(sprintf("%-15s %-18s %8.3f %12.2e\n",
              top_genes$symbol[i], top_genes$ensembl_id[i],
              top_genes$logFC[i], top_genes$padj[i]))
}

# --- 8. Library genes NOT in dream (padj < 0.1) ---
dream_sig_ids <- dream[padj < 0.1, ensembl_id]
lib_not_in_dream <- setdiff(lib_human_ids, dream_sig_ids)
cat(sprintf("\n=== Library genes NOT significant in dream (padj >= 0.1): %d / %d ===\n",
            length(lib_not_in_dream), length(lib_human_ids)))

# Get their dream stats
not_sig <- dream_lib[!ensembl_id %in% dream_sig_ids]
if (nrow(not_sig) > 0) {
  cat(sprintf("  Their LFC distribution: median=%.3f, mean |LFC|=%.3f\n",
              median(not_sig$logFC), mean(abs(not_sig$logFC))))
  cat(sprintf("  Their padj distribution: median=%.3f, mean=%.3f\n",
              median(not_sig$padj, na.rm = TRUE), mean(not_sig$padj, na.rm = TRUE)))
}

# Library genes not even found in dream
lib_not_found <- setdiff(lib_human_ids, dream$ensembl_id)
cat(sprintf("  Not found in dream results at all: %d\n", length(lib_not_found)))
if (length(lib_not_found) > 0 && length(lib_not_found) <= 20) {
  found_syms <- lib_symbol_map[ensembl_id %in% lib_not_found, symbol]
  cat(sprintf("  Missing IDs: %s\n", paste(found_syms, collapse = ", ")))
}

# --- 9. Save results ---
results_dt <- rbindlist(results)
fwrite(results_dt, "results/library/library_vs_dream_overlap.csv")
cat("\nResults saved to results/library/library_vs_dream_overlap.csv\n")

cat("\n=== Summary ===\n")
cat(sprintf("Old library: %d unique human gene targets (from %d mouse genes + orthologs)\n",
            length(lib_human_ids), nrow(lib_core)))
cat(sprintf("Dream (10-cohort, padj<0.1): %d DEGs\n", sum(dream$padj < 0.1, na.rm = TRUE)))
cat(sprintf("Overlap at LFC>=0: %s / %d library genes (%.1f%%)\n",
            results[["0"]]$overlap, length(lib_human_ids),
            results[["0"]]$pct_of_library))
