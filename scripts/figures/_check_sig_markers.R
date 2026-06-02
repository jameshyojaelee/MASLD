library(data.table)
BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

nas_dt <- load_nas_score_progression()
fib_dt <- load_fibrosis_stage_dream()

targets <- c("PNPLA3", "TM6SF2", "HSD17B13", "MBOAT7", "GCKR",
             "CXCL10", "CCL2", "TREM2", "SPP1", "IL1B",
             "COL1A1", "TIMP1", "TGFB1", "ACTA2", "LUM",
             "FASN", "SCD", "PPARA", "CYP7A1", "AKR1B10")

cat("\n=== NAS DEGs (padj<0.1, |LFC|>0.5) ===\n")
nas_sig <- nas_dt[padj < 0.1 & abs(logFC) > 0.5]
for (g in targets) {
  hits <- nas_sig[symbol == g]
  if (nrow(hits) > 0) {
    cat(sprintf("  %s: %d stages sig, max|LFC|=%.2f, min padj=%.1e\n",
                g, nrow(hits), max(abs(hits$logFC)), min(hits$padj)))
  } else {
    cat(sprintf("  %s: NOT significant\n", g))
  }
}

cat("\n=== Fibrosis DEGs (padj<0.1, |LFC|>0.5) ===\n")
fib_sig <- fib_dt[padj < 0.1 & abs(logFC) > 0.5]
for (g in targets) {
  hits <- fib_sig[symbol == g]
  if (nrow(hits) > 0) {
    cat(sprintf("  %s: %d stages sig, max|LFC|=%.2f, min padj=%.1e\n",
                g, nrow(hits), max(abs(hits$logFC)), min(hits$padj)))
  } else {
    cat(sprintf("  %s: NOT significant\n", g))
  }
}

cat("\n=== Combined: significant in either ===\n")
all_sig_sym <- union(nas_sig[symbol %in% targets, unique(symbol)],
                     fib_sig[symbol %in% targets, unique(symbol)])
cat("Significant:", paste(sort(all_sig_sym), collapse=", "), "\n")
cat("Not significant:", paste(sort(setdiff(targets, all_sig_sym)), collapse=", "), "\n")
