#!/usr/bin/env Rscript
# 138_progression_drug_repurposing.R
# ---------------------------------------------------------------------------
# Progression-specific drug repurposing via fgsea with MSigDB C2:CGP.
#
# For each progression contrast, ranks genes by t-statistic and runs fgsea
# against Chemical and Genetic Perturbation (CGP) gene sets. Identifies
# reversal hits (NES < 0, padj < 0.05) — perturbation signatures that
# anti-correlate with the disease signature, suggesting therapeutic reversal.
#
# Contrasts analysed:
#   C1: MASLD vs Control             (onset — from Script 05)
#   C2: NASH vs NAFL                 (steatohepatitis transition — from Script 13)
#   C3: F3-F4 vs F0-F2              (advanced fibrosis — from Script 130/15b)
#   C5: NAS >= 5 vs NAS < 5         (clinical NASH threshold — from Script 130)
#   C6: NASH+F3-F4 vs NAFL+F0-F1   (extreme endpoints — from Script 130)
#   C11: NASH vs Control             (onset — from Script 130)
#   C12: Early vs Late NASH          (progression — from Script 130)
#   C13: NASH vs NAFL (fib-adj)      (progression — from Script 130)
#
# Comparison logic:
#   - Which drugs reverse onset (C1) but NOT progression (C3/C5/C6)?
#   - Which drugs specifically reverse the progression signature?
#   - Overlap counts between onset and each progression contrast
#
# Inputs:
#   - Dream results for each contrast (Ensembl-ID-keyed, with t-statistic)
#   - MSigDB C2:CGP gene sets via msigdbr
#
# Outputs (results/progression/):
#   - progression_drug_reversal_all.csv      — all CGP fgsea results
#   - progression_drug_reversal_hits.csv     — significant reversal hits
#   - progression_drug_comparison.csv        — onset vs progression overlap
#
# Follows patterns from RNA-seq/20_drug_repurposing_cmap.R (Script 20).
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
  library(msigdbr)
})

set.seed(42)

cat("=== Script 138: Progression-Specific Drug Repurposing ===\n")
cat("Start time:", format(Sys.time()), "\n\n")

# ==============================================================================
# Configuration
# ==============================================================================
BASE_DIR <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
INT_DIR  <- file.path(BASE_DIR,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results")
ODIR     <- file.path(INT_DIR, "progression")

FGSEA_PADJ   <- 0.05
FGSEA_MIN    <- 15
FGSEA_MAX    <- 500
FGSEA_NPERM  <- 10000

dir.create(ODIR, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# Define contrast file mapping
# ==============================================================================
# For each contrast, list candidate file paths in priority order.
# Script 130 outputs to results/progression/; earlier scripts output to
# results/disease_signatures/ or results/integration/.
contrast_sources <- list(
  C1 = list(
    label       = "MASLD vs Control (onset)",
    short_label = "C1_onset",
    files = c(
      file.path(INT_DIR, "integration/dream_results.csv")
    )
  ),
  C2 = list(
    label       = "NASH vs NAFL",
    short_label = "C2_NASH_vs_NAFL",
    files = c(
      file.path(INT_DIR, "progression/c2_nafl_vs_nash_dream.csv"),
      file.path(INT_DIR, "disease_signatures/nafl_vs_nash_dream.csv")
    )
  ),
  C3 = list(
    label       = "Advanced vs Early Fibrosis (F3-F4 vs F0-F2)",
    short_label = "C3_adv_fibrosis",
    files = c(
      file.path(INT_DIR, "progression/c3_adv_vs_early_fib_dream.csv"),
      file.path(INT_DIR, "disease_signatures/adv_vs_early_fibrosis_dream.csv")
    )
  ),
  C5 = list(
    label       = "NAS >= 5 vs NAS < 5",
    short_label = "C5_nas_ge5",
    files = c(
      file.path(INT_DIR, "progression/c5_nas_ge5_vs_lt5_dream.csv")
    )
  ),
  C6 = list(
    label       = "Extreme Endpoints (NASH+F3-F4 vs NAFL+F0-F1)",
    short_label = "C6_extreme",
    files = c(
      file.path(INT_DIR, "progression/c6_extreme_endpoints_dream.csv")
    )
  ),
  C11 = list(
    label       = "NASH vs Control (onset)",
    short_label = "C11_nash_vs_ctrl",
    files = c(
      file.path(INT_DIR, "progression/c11_nash_vs_ctrl_dream.csv")
    )
  ),
  C12 = list(
    label       = "Early vs Late NASH",
    short_label = "C12_early_late_nash",
    files = c(
      file.path(INT_DIR, "progression/c12_early_vs_late_nash_dream.csv")
    )
  ),
  C13 = list(
    label       = "NASH vs NAFL (fibrosis-adjusted)",
    short_label = "C13_nash_nafl_fibadj",
    files = c(
      file.path(INT_DIR, "progression/c13_nash_vs_nafl_fib_adj_dream.csv")
    )
  )
)

# ==============================================================================
# Helper: load dream results and build ranked t-statistic vector
# ==============================================================================
load_ranked_stats <- function(spec) {
  # Try each candidate file in order
 found_file <- NULL
  for (f in spec$files) {
    if (file.exists(f)) {
      found_file <- f
      break
    }
  }
  if (is.null(found_file)) {
    cat(sprintf("  WARNING: No file found for %s. Tried:\n", spec$label))
    for (f in spec$files) cat(sprintf("    - %s\n", f))
    return(NULL)
  }

  cat(sprintf("  Loading: %s\n", basename(found_file)))
  dt <- fread(found_file)

  # Identify the gene ID column — may be 'gene' as header or first column
  if ("gene" %in% names(dt)) {
    gene_col <- "gene"
  } else {
    # Some files have gene as the first column without header name
    gene_col <- names(dt)[1]
  }

  # Identify t-statistic column
  if ("t" %in% names(dt)) {
    t_col <- "t"
  } else {
    cat(sprintf("  ERROR: No 't' column in %s. Columns: %s\n",
      basename(found_file), paste(names(dt), collapse = ", ")))
    return(NULL)
  }

  # Handle multi-contrast files (e.g., nas_score_dream.csv with 'contrast' column)
  # These are not expected here, but guard against it
  if ("contrast" %in% names(dt)) {
    cat(sprintf("  NOTE: Multi-contrast file detected. Using first contrast only.\n"))
    first_contrast <- dt[, unique(contrast)][1]
    dt <- dt[contrast == first_contrast]
    cat(sprintf("  Using contrast: %s (%d genes)\n", first_contrast, nrow(dt)))
  }

  # Strip Ensembl version suffix for consistency
  dt[, ensembl_id := gsub("\\..*", "", get(gene_col))]

  # Deduplicate: keep most significant per Ensembl ID
  p_col <- intersect(c("P.Value", "P.value", "pvalue"), names(dt))
  if (length(p_col) > 0) {
    setorderv(dt, p_col[1])
  }
  dt <- dt[!duplicated(ensembl_id)]

  # Build named vector
  stats <- dt[[t_col]]
  names(stats) <- dt$ensembl_id
  stats <- stats[!is.na(stats)]
  stats <- sort(stats, decreasing = TRUE)

  cat(sprintf("  Ranked list: %d genes, t-stat range [%.2f, %.2f]\n",
    length(stats), min(stats), max(stats)))

  return(list(stats = stats, file = found_file, n_genes = length(stats)))
}

# ==============================================================================
# Load MSigDB C2:CGP gene sets
# ==============================================================================
cat("--- Loading MSigDB C2:CGP gene sets ---\n")

cgp_df <- as.data.table(msigdbr(species = "Homo sapiens",
                                 collection = "C2",
                                 subcollection = "CGP"))
cat(sprintf("  C2:CGP gene sets: %d\n", uniqueN(cgp_df$gs_name)))

# Build pathway list using Ensembl IDs
cgp_df <- cgp_df[ensembl_gene != ""]
pathways <- split(cgp_df$ensembl_gene, cgp_df$gs_name)
cat(sprintf("  Total gene-set-to-gene mappings: %d\n", nrow(cgp_df)))

# ==============================================================================
# CGP class annotation (from Script 20)
# ==============================================================================
classify_cgp <- function(pathway_name) {
  pn <- toupper(pathway_name)
  if (grepl("LIVER|HEPAT|HCC|NAFLD|NASH|STEATO|CIRR", pn)) return("Liver/MASLD")
  if (grepl("DRUG|TREAT|COMPOUND|DOXO|ETOPO|TAMOX|METFORM", pn)) return("Drug_Treatment")
  if (grepl("OBESE|OBESITY|BMI|ADIPOSE|FAT|LIPID|CHOLEST", pn)) return("Metabolic")
  if (grepl("INFLAM|TNF|IL6|NFKB|IMMUNE|MACROPHAGE", pn)) return("Inflammatory")
  if (grepl("FIBROSIS|STELLATE|COLLAGEN|TGF", pn)) return("Fibrosis")
  return("Other_CGP")
}

# ==============================================================================
# Run fgsea for each contrast
# ==============================================================================
cat("\n--- Running fgsea for each progression contrast ---\n")

all_results   <- list()
loaded_contrasts <- character()

for (cid in names(contrast_sources)) {
  spec <- contrast_sources[[cid]]
  cat(sprintf("\n[%s] %s\n", cid, spec$label))

  loaded <- load_ranked_stats(spec)
  if (is.null(loaded)) {
    cat(sprintf("  SKIPPED: %s (no input file available)\n", cid))
    next
  }

  # Count pathway overlap with this contrast's gene list
  n_usable <- sum(vapply(pathways, function(g) {
    sum(g %in% names(loaded$stats)) >= FGSEA_MIN
  }, logical(1)))
  cat(sprintf("  Pathways with >= %d overlapping genes: %d\n", FGSEA_MIN, n_usable))

  # Run fgsea
  cat(sprintf("  Running fgsea (nPermSimple=%d)...\n", FGSEA_NPERM))
  set.seed(42)
  fres <- fgsea(pathways    = pathways,
                stats       = loaded$stats,
                minSize     = FGSEA_MIN,
                maxSize     = FGSEA_MAX,
                nPermSimple = FGSEA_NPERM)
  fres_dt <- as.data.table(fres)

  # Annotate
  fres_dt[, contrast_id    := cid]
  fres_dt[, contrast_label := spec$label]
  fres_dt[, short_label    := spec$short_label]
  fres_dt[, cgp_class      := vapply(pathway, classify_cgp, character(1))]

  n_sig       <- sum(fres_dt$padj < FGSEA_PADJ, na.rm = TRUE)
  n_reversal  <- sum(fres_dt$padj < FGSEA_PADJ & fres_dt$NES < 0, na.rm = TRUE)
  n_concordant <- sum(fres_dt$padj < FGSEA_PADJ & fres_dt$NES > 0, na.rm = TRUE)
  cat(sprintf("  Results: %d tested, %d significant (padj < %.2f)\n",
    nrow(fres_dt), n_sig, FGSEA_PADJ))
  cat(sprintf("  Reversal (NES<0): %d | Concordant (NES>0): %d\n",
    n_reversal, n_concordant))

  if (n_reversal > 0) {
    liver_rev <- sum(fres_dt$padj < FGSEA_PADJ & fres_dt$NES < 0 &
      fres_dt$cgp_class == "Liver/MASLD", na.rm = TRUE)
    cat(sprintf("  Liver/MASLD reversal hits: %d\n", liver_rev))
  }

  all_results[[cid]] <- fres_dt
  loaded_contrasts   <- c(loaded_contrasts, cid)

  gc()
}

if (length(all_results) == 0) {
  cat("\nERROR: No contrasts could be loaded. Exiting.\n")
  quit(status = 1)
}

# ==============================================================================
# Combine results
# ==============================================================================
cat("\n--- Combining results across contrasts ---\n")

# Flatten leadingEdge to semicolon-separated string for CSV output
flatten_le <- function(dt) {
  out <- copy(dt)
  if ("leadingEdge" %in% names(out)) {
    out[, leadingEdge := vapply(leadingEdge, function(x) {
      paste(x, collapse = ";")
    }, character(1))]
  }
  out
}

combined <- rbindlist(all_results, fill = TRUE)
cat(sprintf("  Total rows: %d across %d contrasts\n",
  nrow(combined), length(all_results)))

# Save all results
fwrite(flatten_le(combined[order(contrast_id, padj)]),
  file.path(ODIR, "progression_drug_reversal_all.csv"))
cat("  Saved: progression_drug_reversal_all.csv\n")

# Save reversal hits only
hits <- combined[!is.na(NES) & NES < 0 & !is.na(padj) & padj < FGSEA_PADJ]
setorder(hits, contrast_id, NES)
fwrite(flatten_le(hits),
  file.path(ODIR, "progression_drug_reversal_hits.csv"))
cat(sprintf("  Saved: progression_drug_reversal_hits.csv (%d rows)\n", nrow(hits)))

# ==============================================================================
# Onset vs Progression comparison
# ==============================================================================
cat("\n--- Comparing onset vs progression drug hits ---\n")

if ("C1" %in% loaded_contrasts) {
  onset_hits <- hits[contrast_id == "C1", unique(pathway)]
  cat(sprintf("  C1 (onset) reversal hits: %d\n", length(onset_hits)))

  # Compare with each progression contrast
  comparison_rows <- list()

  onset_ids <- c("C1", "C11")
  progression_ids <- setdiff(loaded_contrasts, onset_ids)
  for (pid in progression_ids) {
    prog_spec <- contrast_sources[[pid]]
    prog_hits <- hits[contrast_id == pid, unique(pathway)]
    cat(sprintf("  %s reversal hits: %d\n", pid, length(prog_hits)))

    shared       <- intersect(onset_hits, prog_hits)
    onset_only   <- setdiff(onset_hits, prog_hits)
    prog_only    <- setdiff(prog_hits, onset_hits)

    cat(sprintf("    Shared with C1: %d | Onset-only: %d | Progression-only: %d\n",
      length(shared), length(onset_only), length(prog_only)))

    comparison_rows[[pid]] <- data.table(
      comparison         = paste0("C1_vs_", pid),
      onset_label        = "C1_onset",
      progression_label  = prog_spec$short_label,
      n_onset_hits       = length(onset_hits),
      n_progression_hits = length(prog_hits),
      n_shared           = length(shared),
      n_onset_only       = length(onset_only),
      n_progression_only = length(prog_only),
      jaccard            = round(length(shared) /
        length(union(onset_hits, prog_hits)), 3),
      shared_pathways    = paste(head(shared, 20), collapse = ";"),
      progression_specific_top10 = paste(head(prog_only, 10), collapse = ";")
    )
  }

  comparison_dt <- rbindlist(comparison_rows)
  fwrite(comparison_dt, file.path(ODIR, "progression_drug_comparison.csv"))
  cat(sprintf("  Saved: progression_drug_comparison.csv (%d rows)\n",
    nrow(comparison_dt)))

  # Print summary table
  cat("\n  === Onset vs Progression Drug Hit Comparison ===\n")
  cat(sprintf("  %-20s %6s %6s %6s %7s\n",
    "Comparison", "Onset", "Prog", "Shared", "Jaccard"))
  cat(sprintf("  %s\n", strrep("-", 55)))
  for (i in seq_len(nrow(comparison_dt))) {
    r <- comparison_dt[i]
    cat(sprintf("  %-20s %6d %6d %6d %7.3f\n",
      r$comparison, r$n_onset_hits, r$n_progression_hits,
      r$n_shared, r$jaccard))
  }
} else {
  cat("  WARNING: C1 (onset) not loaded — skipping onset vs progression comparison.\n")
  fwrite(data.table(note = "C1 onset not available; comparison skipped"),
    file.path(ODIR, "progression_drug_comparison.csv"))
  cat("  Saved: progression_drug_comparison.csv (placeholder)\n")
}

# ==============================================================================
# Per-contrast summary
# ==============================================================================
cat("\n--- Per-contrast reversal summary ---\n")
cat(sprintf("  %-20s %8s %8s %8s %8s\n",
  "Contrast", "Tested", "Sig", "Reversal", "Liver"))
cat(sprintf("  %s\n", strrep("-", 60)))

for (cid in loaded_contrasts) {
  cdt <- combined[contrast_id == cid]
  n_tested   <- nrow(cdt)
  n_sig      <- sum(cdt$padj < FGSEA_PADJ, na.rm = TRUE)
  n_rev      <- sum(cdt$padj < FGSEA_PADJ & cdt$NES < 0, na.rm = TRUE)
  n_liver    <- sum(cdt$padj < FGSEA_PADJ & cdt$NES < 0 &
    cdt$cgp_class == "Liver/MASLD", na.rm = TRUE)
  cat(sprintf("  %-20s %8d %8d %8d %8d\n", cid, n_tested, n_sig, n_rev, n_liver))
}

# ==============================================================================
# Progression-specific pathway highlight
# ==============================================================================
cat("\n--- Top 10 reversal hits per progression contrast ---\n")

for (cid in loaded_contrasts) {
  cdt_hits <- hits[contrast_id == cid]
  if (nrow(cdt_hits) == 0) {
    cat(sprintf("\n  [%s] No reversal hits.\n", cid))
    next
  }
  top10 <- head(cdt_hits, 10)
  cat(sprintf("\n  [%s] %s — Top 10 reversal hits:\n",
    cid, contrast_sources[[cid]]$label))
  for (i in seq_len(nrow(top10))) {
    cat(sprintf("    %2d. %-55s NES=%6.2f  padj=%.2e  [%s]\n",
      i, substr(top10$pathway[i], 1, 55),
      top10$NES[i], top10$padj[i], top10$cgp_class[i]))
  }
}

# ==============================================================================
# Summary
# ==============================================================================
cat("\n=== Script 138: Progression Drug Repurposing Complete ===\n")
cat(sprintf("  Contrasts loaded:       %d / %d\n",
  length(loaded_contrasts), length(contrast_sources)))
cat(sprintf("  Loaded:                 %s\n", paste(loaded_contrasts, collapse = ", ")))
skipped <- setdiff(names(contrast_sources), loaded_contrasts)
if (length(skipped) > 0) {
  cat(sprintf("  Skipped (no file):      %s\n", paste(skipped, collapse = ", ")))
}
cat(sprintf("  Total reversal hits:    %d\n", nrow(hits)))
cat(sprintf("  Output dir:             %s\n", ODIR))
cat("End time:", format(Sys.time()), "\n")
