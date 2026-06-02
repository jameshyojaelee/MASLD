#!/usr/bin/env Rscript
# 217a_hormone_tf_targetsets_clean.R -- NON-CIRCULAR TF target sets (B5b sensitivity)
#
# Agent: B5b (C5 critique fix)
#
# CRITICAL CHANGE vs 217a: REMOVE disease_regulons.csv ingestion entirely.
# disease_regulons.csv is by-construction filtered for MASLD-vs-Normal differential
# activity (Mann-Whitney in 04b_scenic_grn_from_activity.py L580-620). Testing
# "are these TFs' targets enriched in COLOC genes" using disease-filtered targets
# is tautological -- the targets were pre-selected to "move in disease".
#
# Clean target-set sources:
#   1. SCENIC+ hepatocyte_regulons.csv (MASLD-agnostic, hepatocyte regulon discovery
#      based on activity in hepatocytes irrespective of disease state)
#   2. TRRUSTv2 (literature-curated, disease-agnostic)
#
# Also adds negative-control TFs (CTCF, MYC, RFX5) -- sex-agnostic generic TFs.
# If these show comparable enrichment, the signal is generic, not hormone-specific.
#
# Inputs:
#   - Analysis/ATAC/Human_Multiome/scenic_plus/hepatocyte_regulons.csv
#   - data/trrust_rawdata.human.tsv
#
# Outputs:
#   - RNA-seq/results/stratified_causal/hormone_tf_targetsets_long_clean.csv
#   - RNA-seq/results/stratified_causal/hormone_tf_targetsets_summary_clean.csv

suppressPackageStartupMessages({
  library(data.table)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
outdir <- file.path(BASE, "RNA-seq/results/stratified_causal")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Primary TF panel (sex-dimorphic hepatic TFs)
HORMONE_TF_PANEL <- c("AR", "ESR1", "ESR2", "FOXA1", "FOXA2",
                      "STAT5A", "STAT5B", "BCL6", "CUX2", "HNF4A")
# Negative-control TFs (sex-agnostic, broadly expressed)
NEG_CTRL_TFS <- c("CTCF", "MYC", "RFX5")
TF_PANEL <- c(HORMONE_TF_PANEL, NEG_CTRL_TFS)

cat("=== 217a_clean: NON-CIRCULAR TF target sets ===\n")
cat("  Hormone TFs: ", paste(HORMONE_TF_PANEL, collapse = ", "), "\n")
cat("  Negative controls: ", paste(NEG_CTRL_TFS, collapse = ", "), "\n")
cat("  *** disease_regulons.csv is INTENTIONALLY NOT LOADED ***\n")

# ===========================================================================
# 1. SCENIC+ hepatocyte regulons (MASLD-agnostic, hepatocyte-defined)
# ===========================================================================
cat("\n--- SCENIC+ hepatocyte regulons (clean) ---\n")
scenic_hep_path <- file.path(BASE,
  "Analysis/ATAC/Human_Multiome/scenic_plus/hepatocyte_regulons.csv")
scenic_long <- data.table()
if (file.exists(scenic_hep_path)) {
  sh <- fread(scenic_hep_path)
  scenic_long <- unique(sh[tf_name %in% TF_PANEL,
                            .(tf = tf_name, target_gene, source = "SCENIC_hep")])
  cat("  Hepatocyte regulons:", nrow(scenic_long), "edges\n")
  cat("  TFs with hep regulon: ",
      paste(sort(unique(scenic_long$tf)), collapse = ", "), "\n")
} else {
  cat("  WARNING: hepatocyte_regulons.csv not found at", scenic_hep_path, "\n")
}

# ===========================================================================
# 2. TRRUSTv2 (literature-curated)
# ===========================================================================
cat("\n--- TRRUSTv2 (literature) ---\n")
trrust_path <- file.path(BASE, "data/trrust_rawdata.human.tsv")
trrust_dt <- data.table()
if (file.exists(trrust_path)) {
  trrust <- fread(trrust_path, header = FALSE,
                  col.names = c("tf", "target_gene", "mode", "pmid"))
  trrust <- trrust[tf %in% TF_PANEL]
  trrust_dt <- unique(trrust[, .(tf, target_gene, source = "TRRUSTv2")])
  cat("  TRRUSTv2:", nrow(trrust_dt), "edges for panel TFs\n")
  cat("  TFs with TRRUST: ",
      paste(sort(unique(trrust_dt$tf)), collapse = ", "), "\n")
}

# ===========================================================================
# 3. Combine + self-loops
# ===========================================================================
cat("\n--- Combining ---\n")
combined <- rbindlist(list(scenic_long, trrust_dt), use.names = TRUE, fill = TRUE)
combined <- combined[!is.na(target_gene) & target_gene != "" & target_gene != "."]
self_loops <- data.table(tf = TF_PANEL, target_gene = TF_PANEL,
                          source = "self_loop")
combined <- unique(rbind(combined, self_loops, fill = TRUE))

fwrite(combined,
       file.path(outdir, "hormone_tf_targetsets_long_clean.csv"))
cat("  Wrote hormone_tf_targetsets_long_clean.csv (", nrow(combined), "rows)\n")

# Summary
summary_dt <- combined[, .(n_unique_targets = uniqueN(target_gene),
                            sources = paste(sort(unique(source)), collapse = ";")),
                       by = tf]
missing_tfs <- setdiff(TF_PANEL, summary_dt$tf)
if (length(missing_tfs) > 0) {
  summary_dt <- rbind(summary_dt,
                      data.table(tf = missing_tfs, n_unique_targets = 0,
                                 sources = "none"))
}
summary_dt[, is_negative_control := tf %in% NEG_CTRL_TFS]
summary_dt <- summary_dt[order(is_negative_control, -n_unique_targets)]
print(summary_dt)
fwrite(summary_dt,
       file.path(outdir, "hormone_tf_targetsets_summary_clean.csv"))
cat("\n=== 217a_clean complete ===\n")
