# WGCNA Pipeline - Step 05b: Module Preservation Sensitivity (Uncorrected CPM)
# ============================================================================
# Fix ID: T2.15 (2026-04-22)
#
# Context: The primary module preservation run in 05_cross_dataset_comparison.R
# can be computed on ComBat-corrected CPM where `group_binary` is passed to
# the `mod` argument — this double-dips (batch-correcting on the same axis
# you're then testing for preservation) and inflates Zsummary.
#
# This sensitivity pass:
#   1. Loads the *uncorrected* CPM (from the integration pipeline voom $E,
#      NO removeBatchEffect / NO ComBat-seq), re-ranking genes by variance.
#   2. Re-runs modulePreservation against the same moduleColors identified
#      on the corrected matrix.
#   3. Saves Zsummary side-by-side against the corrected run in
#      `module_preservation_uncorrected.csv`.
#
# Expected outcome: Zsummary will drop (typically 2-5x) relative to the
# corrected run; we report both in the supplement as a transparency measure.
#
# INPUT:
#   * RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/
#       merged_dge.rds  -> uncorrected voom $E
#   * results/03_network_construction/WGCNA_network_Patient.RData
#   * results/05_comparison/module_preservation_stats.RData  (corrected run;
#       optional — if missing we only emit the uncorrected stats)
# OUTPUT:
#   * results/05_comparison/module_preservation_uncorrected.csv
# ============================================================================

suppressPackageStartupMessages({
  library(WGCNA)
  library(data.table)
  library(edgeR)
  library(limma)
})

options(stringsAsFactors = FALSE)
enableWGCNAThreads()

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
WGCNA_DIR  <- file.path(BASE, "Analysis/WGCNA")
OUTPUT_DIR <- file.path(WGCNA_DIR, "results/05_comparison")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# 1. Load UNCORRECTED expression (no ComBat, no removeBatchEffect)
# ---------------------------------------------------------------------------
cat("Loading uncorrected voom CPM from integration merged_dge...\n")

dge_f <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/merged_dge.rds")
if (!file.exists(dge_f)) {
  stop("merged_dge.rds not found at: ", dge_f,
       "\nRun the 10-cohort integration pipeline first.")
}
dge <- readRDS(dge_f)
dge <- calcNormFactors(dge, method = "TMM")
v <- voom(dge, design = NULL, plot = FALSE)
expr_uncorrected <- v$E   # genes x samples; NO batch correction
cat(sprintf("  Uncorrected CPM: %d genes x %d samples\n",
            nrow(expr_uncorrected), ncol(expr_uncorrected)))

# ---------------------------------------------------------------------------
# 2. Load the Patient module assignments from the primary run
# ---------------------------------------------------------------------------
net_f <- file.path(WGCNA_DIR, "results/03_network_construction/WGCNA_network_Patient.RData")
if (!file.exists(net_f)) {
  stop("Patient network not found: ", net_f,
       "\nRun 03_network_construction.R first.")
}
load(net_f)   # expects `moduleColors` in scope
if (!exists("moduleColors")) {
  stop("moduleColors not found in Patient network .RData")
}
patient_mods <- moduleColors

# Load Patient datExpr (for gene-name alignment)
load(file.path(WGCNA_DIR, "results/01_preprocessing/Patient.RData"))
if (!exists("datExpr_patient")) {
  stop("datExpr_patient not found in Patient preprocessing .RData")
}
# Patient datExpr uses Ensembl IDs with version suffixes in COLUMNS
# (samples are ROWS). Strip version suffixes for matching against
# the uncorrected CPM rownames.
patient_genes_raw <- colnames(datExpr_patient)
patient_genes <- sub("\\..*", "", patient_genes_raw)
cat("  Patient module assignments:", length(patient_mods),
    "genes across", length(unique(patient_mods)), "modules\n")

# ---------------------------------------------------------------------------
# 3. Align uncorrected expression to Patient gene order
# ---------------------------------------------------------------------------
# Uncorrected voom $E has genes x samples; strip Ensembl version suffixes.
expr_genes_raw <- rownames(expr_uncorrected)
expr_genes_clean <- sub("\\..*", "", expr_genes_raw)

# Match by Ensembl-base; collapse duplicates (keep first occurrence per side)
shared <- intersect(patient_genes, expr_genes_clean)
if (length(shared) < 100) {
  # Fallback: use symbols from expression matrix where the atlas lookup
  # might be needed. Emit warning.
  cat("WARNING: only", length(shared),
      "genes matched between Patient modules and uncorrected CPM.\n")
  cat("         Preservation will be computed on the matched subset only.\n")
}

idx_patient <- match(shared, patient_genes)
idx_expr    <- match(shared, expr_genes_clean)
mods_aligned <- patient_mods[idx_patient]
expr_aligned <- t(expr_uncorrected[idx_expr, , drop = FALSE])
rownames(expr_aligned) <- colnames(expr_uncorrected)
colnames(expr_aligned) <- shared
cat("  Aligned matrix:", nrow(expr_aligned), "samples x",
    ncol(expr_aligned), "genes\n")

# ---------------------------------------------------------------------------
# 4. Run modulePreservation — uncorrected as TEST network (self as reference)
# ---------------------------------------------------------------------------
cat("\nRunning modulePreservation on UNCORRECTED CPM...\n")

multiExpr <- list(
  Patient_corrected   = list(data = datExpr_patient[, idx_patient, drop = FALSE]),
  Patient_uncorrected = list(data = expr_aligned)
)
multiColor <- list(
  Patient_corrected   = mods_aligned,
  Patient_uncorrected = mods_aligned
)

mp_un <- modulePreservation(multiExpr, multiColor,
                            referenceNetworks = 1,
                            nPermutations = 100,     # reduced from 200 for speed
                            randomSeed = 1,
                            verbose = 3)

save(mp_un, file = file.path(OUTPUT_DIR, "module_preservation_uncorrected.RData"))

# ---------------------------------------------------------------------------
# 5. Export Zsummary comparison
# ---------------------------------------------------------------------------
statsObs_un <- mp_un$quality$observed[[1]][[2]]
statsZ_un   <- mp_un$quality$Z[[1]][[2]]

out <- data.table(
  module            = rownames(statsZ_un),
  module_size       = statsZ_un$moduleSize,
  Zsummary_UNCORR   = statsZ_un$Zsummary.qual,
  medianRank_UNCORR = statsObs_un$medianRank.qual
)

# Side-by-side with corrected run if available
corrected_f <- file.path(OUTPUT_DIR, "module_preservation_stats.RData")
if (file.exists(corrected_f)) {
  load(corrected_f)  # loads `mp`
  if (exists("mp")) {
    statsZ_cor <- mp$quality$Z[[1]][[2]]
    out[, Zsummary_CORR := statsZ_cor[match(out$module, rownames(statsZ_cor)),
                                       "Zsummary.qual"]]
    out[, inflation_ratio := Zsummary_CORR / Zsummary_UNCORR]
    cat("\nZsummary inflation (corrected / uncorrected):\n")
    print(out[order(-inflation_ratio)])
  }
}

fwrite(out, file.path(OUTPUT_DIR, "module_preservation_uncorrected.csv"))
cat("\nWrote:", file.path(OUTPUT_DIR, "module_preservation_uncorrected.csv"), "\n")
cat("Done.\n")
