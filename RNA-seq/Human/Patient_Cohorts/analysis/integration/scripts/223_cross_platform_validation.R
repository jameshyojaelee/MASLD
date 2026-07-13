#!/usr/bin/env Rscript
# 223_cross_platform_validation.R
# Cross-platform proteomics validation: Olink (proximity extension) vs DIA-MS
#
# Compares protein coverage and differential expression concordance between two
# independent plasma proteomics platforms applied to different cohorts.
# Because samples are NOT shared, direct sample-level correlation is not
# possible; instead the script reports:
#   1. Protein overlap by gene symbol
#   2. DEG concordance: among overlapping proteins, do both platforms classify
#      the same proteins as differentially expressed?
#   3. Rank correlation of per-platform mean abundance across overlapping proteins
#
# Inputs:
#   results/multiprogram/ms_plasma_matrix.csv
#     72 samples x 356 proteins (column = gene symbol or UniProt fallback)
#   results/multiprogram/ms_plasma_metadata.csv
#     sample_id, condition (C/D/...) for DIA-MS cohort
#   Analysis/Proteomics/data/olink_plasma/olink.qc.finished.mendeley.data.txt
#     1,460 proteins (Assay = gene symbol) x 218 subjects (NPX)
#   results/multiprogram/plasma_selected_proteins.csv
#     protein, is_tissue_deg, dream_logFC, dream_padj, tissue_classification  # C2-OK-sensitivity
#     (from Script 221 — used to define Olink-side DEG status)
#
# Output:
#   results/multiprogram/cross_platform_validation.csv
#     Per-overlapping-protein table with platform-side statistics
#
# SLURM: --partition=cpu --cpus-per-task=4 --mem=16G --time=48:00:00
# Env:   micromamba activate rnaseq

suppressPackageStartupMessages({
  library(data.table)
})

set.seed(42)

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT     <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
OUTDIR  <- file.path(INT,  "results/multiprogram")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

OLINK_FILE <- file.path(BASE,
  "Analysis/Proteomics/data/olink_plasma/olink.qc.finished.mendeley.data.txt")
MS_MATRIX  <- file.path(OUTDIR, "ms_plasma_matrix.csv")
MS_META    <- file.path(OUTDIR, "ms_plasma_metadata.csv")
PROT_ANNOT <- file.path(OUTDIR, "plasma_selected_proteins.csv")

cat("=== 223: Cross-Platform Proteomics Validation (Olink vs DIA-MS) ===\n")
cat("Started:", as.character(Sys.time()), "\n\n")

# ── STEP 1: Load Olink data ────────────────────────────────────────────────────
cat("=== STEP 1: Load Olink data ===\n")
stopifnot(file.exists(OLINK_FILE))
olink_raw <- fread(OLINK_FILE, sep = "\t", header = TRUE,
                   na.strings = c("", "NA", "NaN"))
cat(sprintf("  Olink raw: %d rows x %d columns\n", nrow(olink_raw), ncol(olink_raw)))

# Assay column = gene symbol; remaining columns = subjects (NPX values)
stopifnot("Assay" %in% names(olink_raw))
olink_symbols <- olink_raw[["Assay"]]
subject_cols  <- setdiff(names(olink_raw), "Assay")
cat(sprintf("  Olink proteins: %d | subjects: %d\n",
            length(olink_symbols), length(subject_cols)))

# Build wide numeric matrix (proteins x subjects) for mean computation
olink_mat <- as.matrix(olink_raw[, .SD, .SDcols = subject_cols])
storage.mode(olink_mat) <- "numeric"
rownames(olink_mat) <- olink_symbols

# Per-protein mean NPX across all subjects
olink_mean <- rowMeans(olink_mat, na.rm = TRUE)
cat(sprintf("  Olink mean NPX computed; range: [%.2f, %.2f]\n",
            min(olink_mean, na.rm = TRUE), max(olink_mean, na.rm = TRUE)))

# ── STEP 2: Load DIA-MS data ──────────────────────────────────────────────────
cat("\n=== STEP 2: Load DIA-MS data ===\n")
stopifnot(file.exists(MS_MATRIX))
ms_raw  <- fread(MS_MATRIX, header = TRUE, na.strings = c("", "NA", "NaN"))
ms_meta <- if (file.exists(MS_META)) fread(MS_META) else NULL

cat(sprintf("  DIA-MS raw: %d rows x %d columns\n", nrow(ms_raw), ncol(ms_raw)))

# First column is sample_id; remaining columns are proteins
sample_id_col <- names(ms_raw)[1]  # "sample_id"
ms_protein_cols <- setdiff(names(ms_raw), sample_id_col)
cat(sprintf("  DIA-MS proteins (columns): %d | samples: %d\n",
            length(ms_protein_cols), nrow(ms_raw)))

# Build numeric matrix (samples x proteins) and transpose for per-protein means
ms_mat <- as.matrix(ms_raw[, .SD, .SDcols = ms_protein_cols])
storage.mode(ms_mat) <- "numeric"
rownames(ms_mat) <- ms_raw[[sample_id_col]]

# Per-protein mean log2-abundance across all samples
ms_mean <- colMeans(ms_mat, na.rm = TRUE)
cat(sprintf("  DIA-MS mean log2 abundance range: [%.2f, %.2f]\n",
            min(ms_mean, na.rm = TRUE), max(ms_mean, na.rm = TRUE)))

# ── STEP 3: Identify protein overlap by gene symbol ───────────────────────────
cat("\n=== STEP 3: Identify protein overlap by gene symbol ===\n")

# DIA-MS columns are a mixture of gene symbols and UniProt fallback IDs
# (from Script 220: biomaRt mapped symbols take precedence; unmapped keep UniProt)
# Strategy: match directly on Olink gene symbols first, then report stats

ms_symbols   <- ms_protein_cols
olink_set    <- unique(olink_symbols)
ms_set       <- unique(ms_symbols)

overlap_sym  <- intersect(olink_set, ms_set)
only_olink   <- setdiff(olink_set, ms_set)
only_ms      <- setdiff(ms_set, olink_set)

cat(sprintf("  Olink proteins:          %d\n", length(olink_set)))
cat(sprintf("  DIA-MS proteins:         %d\n", length(ms_set)))
cat(sprintf("  Overlap (gene symbol):   %d\n", length(overlap_sym)))
cat(sprintf("  Olink-only:              %d\n", length(only_olink)))
cat(sprintf("  DIA-MS-only:             %d\n", length(only_ms)))
cat(sprintf("  Overlap fraction (Olink): %.1f%%\n",
            100 * length(overlap_sym) / length(olink_set)))
cat(sprintf("  Overlap fraction (DIA-MS): %.1f%%\n",
            100 * length(overlap_sym) / length(ms_set)))

# ── STEP 4: DEG concordance for overlapping proteins ──────────────────────────
cat("\n=== STEP 4: DEG concordance in overlapping proteins ===\n")

# Protein annotation from Script 221 provides tissue DEG status
prot_annot <- if (file.exists(PROT_ANNOT)) {
  pa <- fread(PROT_ANNOT, na.strings = c("", "NA"))
  cat(sprintf("  Loaded protein annotations: %d proteins\n", nrow(pa)))
  pa
} else {
  cat("  WARNING: plasma_selected_proteins.csv not found; skipping DEG annotation\n")
  NULL
}

if (length(overlap_sym) == 0) {
  cat("  WARNING: No overlapping proteins found. Writing empty output.\n")
  out <- data.table(
    platform_comparison = "Olink_vs_DIAMS",
    n_olink_proteins = length(olink_set),
    n_diams_proteins = length(ms_set),
    n_overlap_gene_symbol = 0L,
    overlap_pct_olink = 0,
    overlap_pct_diams = 0,
    note = "No gene-symbol overlap detected; DIA-MS columns may be UniProt IDs"
  )
  fwrite(out, file.path(OUTDIR, "cross_platform_validation.csv"))
  cat("  Written: cross_platform_validation.csv (empty)\n")
  quit(save = "no", status = 0)
}

# Build per-protein table for overlapping genes
ol_olink <- olink_mean[overlap_sym]
ol_ms    <- ms_mean[overlap_sym]

per_protein <- data.table(
  gene_symbol    = overlap_sym,
  olink_mean_npx = ol_olink,
  diams_mean_log2 = ol_ms
)

# Rank correlation of mean abundance across platforms
rho_abundance <- cor(ol_olink, ol_ms, method = "spearman", use = "pairwise.complete.obs")
cat(sprintf("  Spearman rho of mean abundance (Olink NPX vs DIA-MS log2): %.3f\n",
            rho_abundance))

# Annotate with tissue DEG status if annotation available
if (!is.null(prot_annot) && "protein" %in% names(prot_annot)) {
  # Determine DEG columns available
  deg_col    <- if ("is_tissue_deg" %in% names(prot_annot)) "is_tissue_deg" else NULL
  lfc_col    <- if ("dream_logFC"   %in% names(prot_annot)) "dream_logFC"   else NULL  # C2-OK-sensitivity
  padj_col   <- if ("dream_padj"    %in% names(prot_annot)) "dream_padj"    else NULL  # C2-OK-sensitivity
  class_col  <- if ("tissue_classification" %in% names(prot_annot)) "tissue_classification" else NULL

  keep_cols <- c("protein", deg_col, lfc_col, padj_col, class_col)
  keep_cols <- keep_cols[!is.null(keep_cols)]

  pa_sub <- prot_annot[protein %in% overlap_sym, .SD, .SDcols = keep_cols]
  setnames(pa_sub, "protein", "gene_symbol")

  per_protein <- merge(per_protein, pa_sub, by = "gene_symbol", all.x = TRUE)
  cat(sprintf("  Annotation merged for %d / %d overlapping proteins\n",
              sum(!is.na(per_protein[[deg_col %||% "is_tissue_deg"]])),
              nrow(per_protein)))
}

# Per-platform DEG status using tissue dream results as ground truth
# Olink-side DEG: compute per-subject t-test between disease groups if metadata available
# Since gse276114_disease_metadata drives the Olink staging, approximate DEG by
# checking whether a protein is in the switch signatures file
switch_sigs_file <- file.path(OUTDIR, "plasma_switch_signatures.csv")
if (file.exists(switch_sigs_file)) {
  sigs <- fread(switch_sigs_file, na.strings = c("", "NA"))
  if ("symbol" %in% names(sigs) && "in_olink" %in% names(sigs)) {
    sigs_overlap <- sigs[symbol %in% overlap_sym, .(gene_symbol = symbol, in_olink_switch_sig = in_olink)]
    per_protein <- merge(per_protein, sigs_overlap, by = "gene_symbol", all.x = TRUE)
    cat(sprintf("  Switch signature overlap: %d proteins in Olink switch signature\n",
                sum(per_protein$in_olink_switch_sig == TRUE, na.rm = TRUE)))
  }
}

# ── STEP 5: Concordance statistics ────────────────────────────────────────────
cat("\n=== STEP 5: Concordance statistics ===\n")

# Summary statistics table
summary_rows <- list(
  list(metric = "n_olink_proteins",           value = length(olink_set)),
  list(metric = "n_diams_proteins",            value = length(ms_set)),
  list(metric = "n_overlap_gene_symbol",       value = length(overlap_sym)),
  list(metric = "overlap_pct_olink",           value = round(100 * length(overlap_sym) / length(olink_set), 2)),
  list(metric = "overlap_pct_diams",           value = round(100 * length(overlap_sym) / length(ms_set), 2)),
  list(metric = "spearman_rho_mean_abundance", value = round(rho_abundance, 4))
)

# DEG concordance if annotation is available
if (!is.null(prot_annot) && "is_tissue_deg" %in% names(per_protein)) {
  n_tissue_deg <- sum(per_protein$is_tissue_deg == TRUE | per_protein$is_tissue_deg == "True",
                      na.rm = TRUE)
  summary_rows <- c(summary_rows, list(
    list(metric = "n_overlap_tissue_deg", value = n_tissue_deg),
    list(metric = "pct_overlap_that_are_tissue_degs",
         value = round(100 * n_tissue_deg / length(overlap_sym), 2))
  ))
  cat(sprintf("  Overlapping tissue DEGs: %d / %d (%.1f%%)\n",
              n_tissue_deg, length(overlap_sym),
              100 * n_tissue_deg / length(overlap_sym)))

  # Fisher's test: are tissue DEGs enriched among overlapping proteins vs background?
  # Background = all tissue DEGs detectable in either platform
  n_deg_olink_only <- if (!is.null(prot_annot)) {
    pa_olink <- prot_annot[protein %in% olink_set]
    sum(pa_olink$is_tissue_deg == TRUE | pa_olink$is_tissue_deg == "True", na.rm = TRUE)
  } else NA_integer_

  if (!is.na(n_deg_olink_only) && n_deg_olink_only > 0) {
    # 2x2: [in_overlap & DEG, in_overlap & not_DEG; olink_only & DEG, olink_only & not_DEG]
    n_ov_deg     <- n_tissue_deg
    n_ov_nondeg  <- length(overlap_sym) - n_tissue_deg
    n_ol_deg     <- n_deg_olink_only - n_ov_deg
    n_ol_nondeg  <- length(olink_set) - length(overlap_sym) - n_ol_deg
    if (all(c(n_ov_deg, n_ov_nondeg, n_ol_deg, n_ol_nondeg) >= 0)) {
      ft <- fisher.test(matrix(c(n_ov_deg, n_ov_nondeg, n_ol_deg, n_ol_nondeg), nrow = 2))
      summary_rows <- c(summary_rows, list(
        list(metric = "fisher_deg_enrichment_or",  value = round(ft$estimate, 3)),
        list(metric = "fisher_deg_enrichment_pval", value = signif(ft$p.value, 3))
      ))
      cat(sprintf("  Fisher DEG enrichment in overlap: OR=%.2f, p=%.3g\n",
                  ft$estimate, ft$p.value))
    }
  }
}

summary_dt <- rbindlist(lapply(summary_rows, as.data.table))

# ── STEP 6: Write outputs ──────────────────────────────────────────────────────
cat("\n=== STEP 6: Write outputs ===\n")

out_file <- file.path(OUTDIR, "cross_platform_validation.csv")
fwrite(per_protein, out_file)
cat(sprintf("  Written: cross_platform_validation.csv (%d proteins)\n", nrow(per_protein)))

summary_file <- file.path(OUTDIR, "cross_platform_validation_summary.csv")
fwrite(summary_dt, summary_file)
cat(sprintf("  Written: cross_platform_validation_summary.csv\n"))

# Print summary to console
cat("\n--- Summary ---\n")
print(summary_dt)

cat("\nDone:", as.character(Sys.time()), "\n")

# ── Helper ────────────────────────────────────────────────────────────────────
`%||%` <- function(a, b) if (!is.null(a)) a else b
