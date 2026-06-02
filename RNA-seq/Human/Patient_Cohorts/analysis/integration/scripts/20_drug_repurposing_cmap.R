#!/usr/bin/env Rscript
# 20_drug_repurposing_cmap.R
# ---------------------------------------------------------------------------
# LINCS L1000 CMap Queries via signatureSearch
# Objective: Repurpose existing drugs that reverse the MASLD transcriptome
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(signatureSearch)
  library(ExperimentHub)
})

cat("=== Phase 7: Drug Repurposing & CMap ===\n\n")

WD  <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq"
RES <- file.path(WD, "Human/Patient_Cohorts/analysis/integration/results")
OUT <- file.path(RES, "drug_repurposing")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 1. Load Tier 1 & 2 Consensus MASLD Signatures
# ============================================================
cat("Loading consensus DEGs...\n")
degs <- fread(file.path(RES, "disease_signatures/unified_disease_signatures.csv"))

# Filter to robust genes: significant DEGs, sorted by |dvc_dream_logFC|
top_genes <- degs[dvc_dream_sig == TRUE]
top_genes <- top_genes[order(-abs(dvc_dream_logFC))]

# 150 up, 150 down
up_sig <- head(top_genes[dvc_dream_logFC > 0, symbol], 150)
dn_sig <- head(top_genes[dvc_dream_logFC < 0, symbol], 150)

cat("Extracted", length(up_sig), "Up and", length(dn_sig), "Down genes for query.\n")

# ============================================================
# 2. Query Setup: signatureSearch
# ============================================================
cat("\nPreparing to query LINCS database...\n")
cat("Note: Requires local 10GB+ hdf5 database configuration upon actual run.\n")

# Conceptual placeholder for the signatureSearch run:
run_cmap <- function(up, dn) {
    # qsig_cmap <- qSig(query = list(up = up, down = dn), gess_method = "Tau")
    # lincs_db <- "/path/to/local/lincs/l1000_v1_pert_sig.h5" 
    # result <- gess_cmap(qsig_cmap, lincs_db, chunk_size = 5000)
    #
    # filter to HepG2:
    # hepg2_hits <- result@result[result@result$cell == "HEPG2", ]
    # return(hepg2_hits)
    cat("  -> (Pipeline to run once GSE213621 is integrated and database downloaded)\n")
}

run_cmap(up_sig, dn_sig)
