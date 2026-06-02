#!/usr/bin/env Rscript
# 217a_hormone_tf_targetsets.R — Build sex-hormone TF target sets
#
# Agent: B5 (Sex-Hormone TF x COLOC) -- HEADLINE A7
#
# For each sex-hormone TF (AR, ESR1, ESR2, FOXA1, FOXA2, STAT5A, STAT5B, BCL6,
# CUX2, HNF4A), build a target set from SCENIC+ (preferred) or TRRUSTv2 fallback.
# Output: long-format CSV (tf, target_gene, source) for downstream enrichment.
#
# Inputs:
#   - Analysis/ATAC/Human_Multiome/scenic_plus/hepatocyte_regulons.csv
#   - Analysis/ATAC/Human_Multiome/scenic_plus/disease_regulons.csv
#   - data/trrust_rawdata.human.tsv (TRRUSTv2 cached 2026-05-11)
#
# Outputs:
#   - RNA-seq/results/stratified_causal/hormone_tf_targetsets_long.csv
#   - RNA-seq/results/stratified_causal/hormone_tf_targetsets_summary.csv

suppressPackageStartupMessages({
  library(data.table)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
outdir <- file.path(BASE, "RNA-seq/results/stratified_causal")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Primary TF panel (sex-dimorphic hepatic TFs)
TF_PANEL <- c("AR", "ESR1", "ESR2", "FOXA1", "FOXA2",
              "STAT5A", "STAT5B", "BCL6", "CUX2", "HNF4A")

cat("=== 217a: Building sex-hormone TF target sets ===\n")
cat("TF panel:", paste(TF_PANEL, collapse = ", "), "\n")

# ===========================================================================
# 1. SCENIC+ hepatocyte regulons (primary source)
# ===========================================================================
cat("\n--- SCENIC+ hepatocyte regulons ---\n")
scenic_hep_path <- file.path(BASE,
  "Analysis/ATAC/Human_Multiome/scenic_plus/hepatocyte_regulons.csv")
scenic_disease_path <- file.path(BASE,
  "Analysis/ATAC/Human_Multiome/scenic_plus/disease_regulons.csv")

scenic_long <- list()
if (file.exists(scenic_hep_path)) {
  sh <- fread(scenic_hep_path)
  # tf_name, target_gene columns are long-format already
  scenic_long[["hep"]] <- unique(sh[tf_name %in% TF_PANEL,
                                    .(tf = tf_name, target_gene, source = "SCENIC_hep")])
  cat("  Hepatocyte regulons: ", nrow(scenic_long[["hep"]]), " edges\n")
}
if (file.exists(scenic_disease_path)) {
  sd <- fread(scenic_disease_path)
  # disease_regulons is wide-format: target_genes is semicolon-separated
  sd <- sd[tf_name %in% TF_PANEL]
  if (nrow(sd) > 0) {
    sd_long <- sd[, .(target_gene = unlist(strsplit(target_genes, ";"))),
                  by = .(tf = tf_name)]
    sd_long[, source := "SCENIC_disease"]
    scenic_long[["dis"]] <- unique(sd_long)
    cat("  Disease regulons: ", nrow(scenic_long[["dis"]]), " edges\n")
  }
}
scenic_dt <- rbindlist(scenic_long, use.names = TRUE, fill = TRUE)
if (nrow(scenic_dt) > 0) {
  cat("  TFs with SCENIC+ targets: ",
      paste(sort(unique(scenic_dt$tf)), collapse = ", "), "\n")
} else {
  cat("  No SCENIC+ targets found.\n")
}

# ===========================================================================
# 2. TRRUSTv2 fallback
# ===========================================================================
cat("\n--- TRRUSTv2 fallback ---\n")
trrust_path <- file.path(BASE, "data/trrust_rawdata.human.tsv")
if (!file.exists(trrust_path)) {
  cat("  TRRUSTv2 not cached -- attempting fetch\n")
  tryCatch({
    download.file("https://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv",
                  trrust_path, quiet = TRUE, mode = "wb")
  }, error = function(e) cat("  TRRUSTv2 download failed:", conditionMessage(e), "\n"))
}
trrust_dt <- data.table()
if (file.exists(trrust_path)) {
  trrust <- fread(trrust_path, header = FALSE,
                  col.names = c("tf", "target_gene", "mode", "pmid"))
  trrust <- trrust[tf %in% TF_PANEL]
  trrust_dt <- unique(trrust[, .(tf, target_gene, source = "TRRUSTv2")])
  cat("  TRRUSTv2: ", nrow(trrust_dt), " edges for panel TFs\n")
  cat("  TFs with TRRUST targets: ",
      paste(sort(unique(trrust_dt$tf)), collapse = ", "), "\n")
}

# ===========================================================================
# 3. Combine and write long format
# ===========================================================================
cat("\n--- Combining ---\n")
combined <- rbindlist(list(scenic_dt, trrust_dt), use.names = TRUE, fill = TRUE)
combined <- combined[!is.na(target_gene) & target_gene != "" & target_gene != "."]
# Add self-loop (TF -> TF) for completeness if not present
self_loops <- data.table(tf = TF_PANEL, target_gene = TF_PANEL,
                          source = "self_loop")
combined <- unique(rbind(combined, self_loops, fill = TRUE))

fwrite(combined,
       file.path(outdir, "hormone_tf_targetsets_long.csv"))
cat("  Wrote hormone_tf_targetsets_long.csv (", nrow(combined), "rows)\n")

# Summary by TF
summary_dt <- combined[, .(n_unique_targets = uniqueN(target_gene),
                            sources = paste(sort(unique(source)), collapse = ";")),
                       by = tf]
# Make sure every TF in panel has a row even if 0 targets
missing_tfs <- setdiff(TF_PANEL, summary_dt$tf)
if (length(missing_tfs) > 0) {
  summary_dt <- rbind(summary_dt,
                      data.table(tf = missing_tfs, n_unique_targets = 0,
                                 sources = "none"))
}
summary_dt <- summary_dt[order(-n_unique_targets)]
print(summary_dt)
fwrite(summary_dt,
       file.path(outdir, "hormone_tf_targetsets_summary.csv"))
cat("\n=== 217a complete ===\n")
