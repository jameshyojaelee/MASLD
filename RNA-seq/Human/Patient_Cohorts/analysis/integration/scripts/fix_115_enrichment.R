#!/usr/bin/env Rscript
# Quick fix: re-run only the fgsea enrichment + summary from Script 115
# Uses already-saved intermediate results from the dream run
library(data.table)
library(fgsea)
library(msigdbr)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
OUTDIR <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/progression")

cat("=== Re-running Script 115 enrichment (post-dream) ===\n")

# Load saved intermediate results
fib_dt <- fread(file.path(OUTDIR, "transition_fib_dream_results.csv"))
nas_dt <- fread(file.path(OUTDIR, "transition_nas_dream_results.csv"))
all_tau <- fread(file.path(OUTDIR, "transition_tau_index.csv"))
all_dt <- rbind(fib_dt, nas_dt, fill = TRUE)

cat("Loaded: fib=", nrow(fib_dt), ", nas=", nrow(nas_dt), ", tau=", nrow(all_tau), "\n")

# Fix: merge with allow.cartesian (multiple rows per gene across transitions/stage types)
tau_cols <- intersect(c("gene", "tau", "peak_transition", "max_transition", "max_abs_t", "peak_logFC", "peak_padj"),
                      names(all_tau))
tp <- merge(all_dt, all_tau[, ..tau_cols], by = "gene", all.x = TRUE, allow.cartesian = TRUE)
fwrite(tp, file.path(OUTDIR, "transition_programs_full.csv"))
cat("Saved transition_programs_full.csv (", nrow(tp), "rows)\n")

# Pathway enrichment per transition (fgsea)
cat("\n=== Pathway Enrichment per Transition ===\n")
hallmark <- msigdbr(species = "Homo sapiens", collection = "H")
kegg <- tryCatch(
  msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_MEDICUS"),
  error = function(e) msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY")
)
pathways <- rbind(hallmark[, c("gs_name", "gene_symbol")],
                  kegg[, c("gs_name", "gene_symbol")])
pathway_list <- split(pathways$gene_symbol, pathways$gs_name)

run_fgsea_transition <- function(dt, transition_name) {
  sub <- dt[transition == transition_name]
  ranks <- setNames(sub$t, sub$gene)
  ranks <- ranks[!is.na(ranks)]
  ranks <- sort(ranks, decreasing = TRUE)
  res <- fgsea(pathways = pathway_list, stats = ranks,
               minSize = 15, maxSize = 500, nPermSimple = 10000)
  res$transition <- transition_name
  return(as.data.table(res))
}

fib_pathways <- rbindlist(lapply(unique(fib_dt$transition), function(tr) {
  cat("  fgsea:", tr, "\n")
  run_fgsea_transition(fib_dt, tr)
}), fill = TRUE)

nas_pathways <- rbindlist(lapply(unique(nas_dt$transition), function(tr) {
  cat("  fgsea:", tr, "\n")
  run_fgsea_transition(nas_dt, tr)
}), fill = TRUE)

all_pathways <- rbind(fib_pathways, nas_pathways)
all_pathways$leadingEdge <- sapply(all_pathways$leadingEdge, paste, collapse = ";")
cat("Pathway enrichment: total", nrow(all_pathways), "pathway-transition tests\n")
cat("Significant (padj < 0.05):", sum(all_pathways$padj < 0.05, na.rm = TRUE), "\n")
fwrite(all_pathways, file.path(OUTDIR, "transition_pathway_enrichment.csv"))
cat("Saved transition_pathway_enrichment.csv\n")

# Rebuild combined tau with updated flags (padj<0.05 instead of 0.1)
fib_tau <- all_tau[grepl("^F", peak_transition)]
nas_tau <- all_tau[grepl("^NAS", peak_transition)]
fib_tau$stage_type <- "fibrosis"
nas_tau$stage_type <- "NAS"
combined_tau <- rbind(fib_tau, nas_tau, fill = TRUE)
combined_tau$is_transition_unique <- combined_tau$tau > 0.8 & combined_tau$peak_padj < 0.05
combined_tau$is_moderately_specific <- combined_tau$tau > 0.6 & combined_tau$peak_padj < 0.05
fwrite(combined_tau, file.path(OUTDIR, "transition_programs.csv"))
cat("Saved transition_programs.csv (", nrow(combined_tau), "genes)\n")

# Summary table
summary_dt <- rbind(
  fib_dt[, .(n_deg_05 = sum(adj.P.Val < 0.05, na.rm=T),
             n_deg_05_lfc02 = sum(adj.P.Val < 0.05 & abs(logFC) > 0.2, na.rm=T),
             n_deg_05_lfc03 = sum(adj.P.Val < 0.05 & abs(logFC) > 0.3, na.rm=T),
             n_deg_05_lfc05 = sum(adj.P.Val < 0.05 & abs(logFC) > 0.5, na.rm=T),  # New canonical (LOO-CV stability)
             n_deg_01 = sum(adj.P.Val < 0.1, na.rm=T),
             n_up = sum(adj.P.Val < 0.05 & logFC > 0, na.rm=T),
             n_down = sum(adj.P.Val < 0.05 & logFC < 0, na.rm=T),
             n_genes = .N,
             type = "fibrosis"), by = transition],
  nas_dt[, .(n_deg_05 = sum(adj.P.Val < 0.05, na.rm=T),
             n_deg_05_lfc02 = sum(adj.P.Val < 0.05 & abs(logFC) > 0.2, na.rm=T),
             n_deg_05_lfc03 = sum(adj.P.Val < 0.05 & abs(logFC) > 0.3, na.rm=T),
             n_deg_05_lfc05 = sum(adj.P.Val < 0.05 & abs(logFC) > 0.5, na.rm=T),  # New canonical (LOO-CV stability)
             n_deg_01 = sum(adj.P.Val < 0.1, na.rm=T),
             n_up = sum(adj.P.Val < 0.05 & logFC > 0, na.rm=T),
             n_down = sum(adj.P.Val < 0.05 & logFC < 0, na.rm=T),
             n_genes = .N,
             type = "NAS"), by = transition]
)
fwrite(summary_dt, file.path(OUTDIR, "transition_summary.csv"))
cat("\nTransition summary (padj<0.05 primary):\n")
print(summary_dt[, .(transition, n_deg_05, n_deg_05_lfc05, n_up, n_down)])

cat("\n=== Script 115 enrichment fix complete ===\n")
