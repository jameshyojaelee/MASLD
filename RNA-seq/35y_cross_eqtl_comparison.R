#!/usr/bin/env Rscript
# 35y_cross_eqtl_comparison.R
# ---------------------------------------------------------------------------
# Cross-eQTL Source Comparison: Broadaway (N=1,183) vs GTEx v8 Liver (N=208)
#
# Compares COLOC results for UKBB ALT across eQTL sources:
#   1. Broadaway 2024 bulk liver eQTL (N=1,183, largest liver eQTL)
#   2. GTEx v8 Liver (N=208, widely used reference)
#   3. sc-eQTL multi-cell-type (N=312, cell-type-resolved, Ghodsian GWAS)
#
# Analyses:
#   - PP.H4 correlation between eQTL sources (Spearman)
#   - Concordance at PP.H4 > 0.5
#   - Cross-eQTL validated genes (PP.H4 thresholds adjusted for sample size)
#   - Source-unique discoveries (Broadaway-only, GTEx-only)
#   - Three-way comparison when sc-eQTL available
#
# Inputs:
#   - RNA-seq/results/causal_inference/broadaway_ukbb/coloc_results.csv
#   - RNA-seq/results/causal_inference/gtex_v8_ukbb_alt/coloc_results.csv
#   - RNA-seq/results/causal_inference/sceqtl/coloc_results.csv (optional)
#   - RNA-seq/results/causal_inference/sceqtl_ukbb/coloc_results.csv (optional)
#
# Outputs (RNA-seq/results/causal_inference/cross_eqtl_comparison/):
#   - cross_eqtl_comparison.csv  — per-gene comparison (all sources)
#   - cross_eqtl_summary.csv     — summary statistics
#   - cross_eqtl_validated_genes.csv — genes validated across eQTL sources
#
# No figures — see figS_cross_eqtl.R for visualization.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

cat("=== Script 35y: Cross-eQTL Source Comparison ===\n")
cat("Start time:", format(Sys.time()), "\n\n")

# ==============================================================================
# Configuration
# ==============================================================================
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
OUTDIR <- file.path(BASE, "RNA-seq/results/causal_inference/cross_eqtl_comparison")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# Source gene-symbol mapping utility
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

# Input paths
BROADAWAY_FILE <- file.path(BASE, "RNA-seq/results/causal_inference/broadaway_ukbb/coloc_results.csv")
GTEX_FILE      <- file.path(BASE, "RNA-seq/results/causal_inference/gtex_v8_ukbb_alt/coloc_results.csv")
SCEQTL_FILE    <- file.path(BASE, "RNA-seq/results/causal_inference/sceqtl/coloc_results.csv")
SCEQTL_UKBB_FILE <- file.path(BASE, "RNA-seq/results/causal_inference/sceqtl_ukbb/coloc_results.csv")

# Thresholds
# Broadaway PP.H4 > 0.5 (standard, N=1,183 gives strong power)
# GTEx PP.H4 > 0.3 (relaxed, N=208 is 6x smaller → less power)
# sc-eQTL PP.H4 > 0.3 (relaxed, N=312 is 4x smaller)
BROADAWAY_THRESH <- 0.5
GTEX_THRESH      <- 0.3
SCEQTL_THRESH    <- 0.3
STANDARD_THRESH  <- 0.5   # for concordance counting

# ==============================================================================
# Helper: load COLOC results from a file
# ==============================================================================
load_coloc_source <- function(path, source_name, gene_col = "gene") {
  if (!file.exists(path)) {
    cat("  ", source_name, ": FILE NOT FOUND (", path, ")\n")
    return(NULL)
  }
  dt <- fread(path)
  if (nrow(dt) == 0) {
    cat("  ", source_name, ": empty file\n")
    return(NULL)
  }
  cat("  ", source_name, ":", nrow(dt), "rows loaded\n")
  return(dt)
}

# ==============================================================================
# 1. Load Broadaway x UKBB ALT
# ==============================================================================
cat("--- Step 1: Loading eQTL COLOC results ---\n")

brd <- load_coloc_source(BROADAWAY_FILE, "Broadaway x UKBB_ALT")
if (is.null(brd)) {
  stop("Broadaway COLOC results are required. Run 35b_broadaway_coloc_ukbb.R first.")
}

# Broadaway gene col contains gene symbols already; ensembl col has Ensembl IDs
# Standardize: strip version from ensembl, create clean ID for matching
brd[, ensembl_clean := sub("\\..*", "", ensembl)]
brd_slim <- brd[, .(
  ensembl_clean,
  gene_symbol_brd  = gene,
  PP.H4_broadaway  = PP.H4,
  PP.H3_broadaway  = PP.H3,
  n_snps_broadaway = n_snps,
  top_snp_broadaway = top_snp,
  eqtl_N_broadaway = eqtl_N
)]
cat("    Broadaway genes:", nrow(brd_slim), "\n")
cat("    Broadaway PP.H4 > 0.5:", sum(brd_slim$PP.H4_broadaway > STANDARD_THRESH), "\n\n")

# ==============================================================================
# 2. Load GTEx v8 x UKBB ALT
# ==============================================================================
gtex <- load_coloc_source(GTEX_FILE, "GTEx v8 x UKBB_ALT")
has_gtex <- !is.null(gtex)

if (has_gtex) {
  # GTEx gene col may have versioned Ensembl IDs; strip version for matching
  gtex[, ensembl_clean := sub("\\..*", "", gene)]
  # If ensembl column exists, prefer that; otherwise use gene
  if ("ensembl" %in% names(gtex)) {
    gtex[, ensembl_clean := sub("\\..*", "", ensembl)]
  }

  eqtl_n_col <- if ("eqtl_N" %in% names(gtex)) "eqtl_N" else NA
  gtex_slim <- gtex[, .(
    ensembl_clean,
    PP.H4_gtex  = PP.H4,
    PP.H3_gtex  = PP.H3,
    n_snps_gtex = n_snps
  )]
  if (!is.na(eqtl_n_col)) {
    gtex_slim[, eqtl_N_gtex := gtex[[eqtl_n_col]]]
  } else {
    gtex_slim[, eqtl_N_gtex := 208L]
  }

  # Deduplicate (take best PP.H4 per gene if multiple entries)
  gtex_slim <- gtex_slim[order(-PP.H4_gtex)]
  gtex_slim <- gtex_slim[!duplicated(ensembl_clean)]
  cat("    GTEx genes:", nrow(gtex_slim), "\n")
  cat("    GTEx PP.H4 > 0.5:", sum(gtex_slim$PP.H4_gtex > STANDARD_THRESH), "\n")
  cat("    GTEx PP.H4 > 0.3:", sum(gtex_slim$PP.H4_gtex > GTEX_THRESH), "\n\n")
} else {
  cat("    GTEx results not available — pairwise Broadaway-GTEx comparison skipped\n\n")
}

# ==============================================================================
# 3. Load sc-eQTL COLOC results (if available)
# ==============================================================================
# Try UKBB ALT version first, fall back to Ghodsian version
sceqtl_raw <- load_coloc_source(SCEQTL_UKBB_FILE, "sc-eQTL x UKBB_ALT")
sceqtl_gwas_label <- "UKBB_ALT"
if (is.null(sceqtl_raw)) {
  sceqtl_raw <- load_coloc_source(SCEQTL_FILE, "sc-eQTL x Ghodsian")
  sceqtl_gwas_label <- "Ghodsian"
}
has_sceqtl <- !is.null(sceqtl_raw)

if (has_sceqtl) {
  # sc-eQTL has per-cell-type results; take best PP.H4 across cell types per gene
  # Determine which PP.H4 column to use (depends on GWAS source)
  pp4_cols <- grep("^PP\\.H4", names(sceqtl_raw), value = TRUE)
  cat("    sc-eQTL PP.H4 columns:", paste(pp4_cols, collapse = ", "), "\n")

  if ("PP.H4" %in% pp4_cols) {
    # Simple format (e.g., sceqtl_ukbb)
    sceqtl_raw[, best_PP.H4 := PP.H4]
  } else if ("PP.H4_all" %in% pp4_cols) {
    # Multi-condition format — use the "all" (all-cells) condition
    sceqtl_raw[, best_PP.H4 := PP.H4_all]
  } else {
    # Fall back to first PP.H4 column
    sceqtl_raw[, best_PP.H4 := get(pp4_cols[1])]
  }

  # Best PP.H4 per gene across cell types + record the best cell type
  sceqtl_raw[, .best := best_PP.H4 == max(best_PP.H4), by = gene]
  sceqtl_best <- sceqtl_raw[sceqtl_raw[, .I[which.max(best_PP.H4)], by = gene]$V1]

  # sc-eQTL gene col has gene symbols (not Ensembl IDs)
  # Use add_symbols in reverse: need symbol -> ensembl mapping
  # Build reverse map from load_gene_map()
  gm <- load_gene_map()  # ensembl_clean, symbol
  sceqtl_best <- merge(sceqtl_best, gm, by.x = "gene", by.y = "symbol", all.x = TRUE)

  sceqtl_slim <- sceqtl_best[, .(
    ensembl_clean,
    gene_symbol_sceqtl = gene,
    PP.H4_sceqtl       = best_PP.H4,
    best_cell_type      = cell_type,
    sceqtl_gwas         = sceqtl_gwas_label
  )]
  # Drop genes that couldn't map to Ensembl
  n_before <- nrow(sceqtl_slim)
  sceqtl_slim <- sceqtl_slim[!is.na(ensembl_clean) & ensembl_clean != ""]
  if (n_before > nrow(sceqtl_slim)) {
    cat("    sc-eQTL: dropped", n_before - nrow(sceqtl_slim),
        "genes without Ensembl mapping\n")
  }
  sceqtl_slim <- sceqtl_slim[!duplicated(ensembl_clean)]
  cat("    sc-eQTL best-per-gene:", nrow(sceqtl_slim), "genes\n")
  cat("    sc-eQTL PP.H4 > 0.5:", sum(sceqtl_slim$PP.H4_sceqtl > STANDARD_THRESH), "\n")
  cat("    sc-eQTL PP.H4 > 0.3:", sum(sceqtl_slim$PP.H4_sceqtl > SCEQTL_THRESH), "\n")
  cat("    sc-eQTL GWAS source:", sceqtl_gwas_label, "\n\n")
} else {
  cat("    sc-eQTL results not available — three-way comparison skipped\n\n")
}

# ==============================================================================
# 4. Merge: Broadaway + GTEx (+ sc-eQTL if available)
# ==============================================================================
cat("--- Step 2: Merging eQTL sources ---\n")

# Start with Broadaway (always present)
merged <- copy(brd_slim)

# Merge GTEx
if (has_gtex) {
  merged <- merge(merged, gtex_slim, by = "ensembl_clean", all = TRUE)
  cat("  After Broadaway + GTEx merge:", nrow(merged), "genes\n")
  cat("    In both:", sum(!is.na(merged$PP.H4_broadaway) & !is.na(merged$PP.H4_gtex)), "\n")
  cat("    Broadaway-only:", sum(!is.na(merged$PP.H4_broadaway) & is.na(merged$PP.H4_gtex)), "\n")
  cat("    GTEx-only:", sum(is.na(merged$PP.H4_broadaway) & !is.na(merged$PP.H4_gtex)), "\n")
}

# Merge sc-eQTL
if (has_sceqtl) {
  merged <- merge(merged, sceqtl_slim, by = "ensembl_clean", all = TRUE)
  cat("  After sc-eQTL merge:", nrow(merged), "genes\n")
}

# Add gene symbols via add_symbols (fills in missing symbols from atlas)
if (!"gene" %in% names(merged)) {
  merged[, gene := ensembl_clean]
}
# add_symbols() drops its key column (`ensembl_clean := NULL`) internally; the
# merged table's join key IS `ensembl_clean` and is referenced by every
# downstream step (out_cols/comp_cols), so preserve it across the call.
merged[, .ensembl_clean_keep := ensembl_clean]
merged <- add_symbols(merged, "ensembl_clean")
merged[, ensembl_clean := .ensembl_clean_keep]
merged[, .ensembl_clean_keep := NULL]
# Prefer existing symbol columns over the atlas-mapped one
if ("gene_symbol_brd" %in% names(merged)) {
  merged[!is.na(gene_symbol_brd) & gene_symbol_brd != "", symbol := gene_symbol_brd]
}
if ("gene_symbol_sceqtl" %in% names(merged)) {
  merged[!is.na(gene_symbol_sceqtl) & gene_symbol_sceqtl != "" & (is.na(symbol) | symbol == ensembl_clean),
         symbol := gene_symbol_sceqtl]
}

cat("  Final merged table:", nrow(merged), "genes\n\n")

# ==============================================================================
# 5. Pairwise Analyses: Broadaway vs GTEx
# ==============================================================================
cat("--- Step 3: Pairwise analysis — Broadaway vs GTEx ---\n")

summary_rows <- list()

if (has_gtex) {
  shared_bg <- merged[!is.na(PP.H4_broadaway) & !is.na(PP.H4_gtex)]
  n_shared <- nrow(shared_bg)
  cat("  Shared genes:", n_shared, "\n")

  if (n_shared >= 10) {
    # Spearman correlation of PP.H4
    rho_test <- cor.test(shared_bg$PP.H4_broadaway, shared_bg$PP.H4_gtex,
                         method = "spearman", exact = FALSE)
    rho <- rho_test$estimate
    rho_p <- rho_test$p.value
    cat("  Spearman rho (PP.H4):", round(rho, 4), " p =", format(rho_p, digits = 3), "\n")

    # Concordance at standard threshold
    both_sig <- sum(shared_bg$PP.H4_broadaway > STANDARD_THRESH &
                    shared_bg$PP.H4_gtex > STANDARD_THRESH)
    brd_only_sig <- sum(shared_bg$PP.H4_broadaway > STANDARD_THRESH &
                        shared_bg$PP.H4_gtex <= STANDARD_THRESH)
    gtex_only_sig <- sum(shared_bg$PP.H4_broadaway <= STANDARD_THRESH &
                         shared_bg$PP.H4_gtex > STANDARD_THRESH)
    neither_sig <- sum(shared_bg$PP.H4_broadaway <= STANDARD_THRESH &
                       shared_bg$PP.H4_gtex <= STANDARD_THRESH)

    concordance_rate <- (both_sig + neither_sig) / n_shared
    cat("  Concordance at PP.H4 > 0.5: ", round(concordance_rate * 100, 1), "%\n")
    cat("    Both > 0.5:", both_sig, "\n")
    cat("    Broadaway-only:", brd_only_sig, "\n")
    cat("    GTEx-only:", gtex_only_sig, "\n")
    cat("    Neither:", neither_sig, "\n")

    # Cross-validated: Broadaway > 0.5 AND GTEx > 0.3 (relaxed for smaller N)
    cross_validated <- sum(shared_bg$PP.H4_broadaway > BROADAWAY_THRESH &
                           shared_bg$PP.H4_gtex > GTEX_THRESH)
    cat("  Cross-validated (Broadaway > 0.5 & GTEx > 0.3):", cross_validated, "\n")

    # Fisher's exact test: is PP.H4 > 0.5 in Broadaway associated with
    # PP.H4 > 0.3 in GTEx?
    fisher_mat <- matrix(c(
      sum(shared_bg$PP.H4_broadaway > BROADAWAY_THRESH & shared_bg$PP.H4_gtex > GTEX_THRESH),
      sum(shared_bg$PP.H4_broadaway > BROADAWAY_THRESH & shared_bg$PP.H4_gtex <= GTEX_THRESH),
      sum(shared_bg$PP.H4_broadaway <= BROADAWAY_THRESH & shared_bg$PP.H4_gtex > GTEX_THRESH),
      sum(shared_bg$PP.H4_broadaway <= BROADAWAY_THRESH & shared_bg$PP.H4_gtex <= GTEX_THRESH)
    ), nrow = 2, byrow = TRUE)
    fisher_res <- fisher.test(fisher_mat)
    cat("  Fisher's exact OR:", round(fisher_res$estimate, 2),
        " p =", format(fisher_res$p.value, digits = 3), "\n\n")

    summary_rows[["Broadaway_vs_GTEx"]] <- data.table(
      comparison       = "Broadaway_vs_GTEx",
      gwas             = "UKBB_ALT",
      n_shared_genes   = n_shared,
      spearman_rho     = round(rho, 4),
      spearman_p       = rho_p,
      concordance_rate = round(concordance_rate, 4),
      both_sig         = both_sig,
      brd_only_sig     = brd_only_sig,
      gtex_only_sig    = gtex_only_sig,
      neither_sig      = neither_sig,
      cross_validated  = cross_validated,
      fisher_OR        = round(fisher_res$estimate, 2),
      fisher_p         = fisher_res$p.value,
      brd_thresh       = BROADAWAY_THRESH,
      gtex_thresh      = GTEX_THRESH
    )
  } else {
    cat("  Too few shared genes for correlation analysis\n\n")
    summary_rows[["Broadaway_vs_GTEx"]] <- data.table(
      comparison = "Broadaway_vs_GTEx", gwas = "UKBB_ALT",
      n_shared_genes = n_shared, spearman_rho = NA_real_,
      spearman_p = NA_real_, concordance_rate = NA_real_,
      both_sig = NA_integer_, brd_only_sig = NA_integer_,
      gtex_only_sig = NA_integer_, neither_sig = NA_integer_,
      cross_validated = NA_integer_, fisher_OR = NA_real_,
      fisher_p = NA_real_, brd_thresh = BROADAWAY_THRESH,
      gtex_thresh = GTEX_THRESH
    )
  }
} else {
  cat("  Skipped (GTEx results not available)\n\n")
}

# ==============================================================================
# 6. Pairwise Analyses: Broadaway vs sc-eQTL
# ==============================================================================
cat("--- Step 4: Pairwise analysis — Broadaway vs sc-eQTL ---\n")

if (has_sceqtl) {
  shared_bs <- merged[!is.na(PP.H4_broadaway) & !is.na(PP.H4_sceqtl)]
  n_shared_bs <- nrow(shared_bs)
  cat("  Shared genes:", n_shared_bs, "\n")
  cat("  NOTE: sc-eQTL uses", sceqtl_gwas_label, "GWAS (different from UKBB ALT)\n")

  if (n_shared_bs >= 10) {
    rho_bs <- cor.test(shared_bs$PP.H4_broadaway, shared_bs$PP.H4_sceqtl,
                       method = "spearman", exact = FALSE)
    cat("  Spearman rho (PP.H4):", round(rho_bs$estimate, 4),
        " p =", format(rho_bs$p.value, digits = 3), "\n")

    both_bs <- sum(shared_bs$PP.H4_broadaway > BROADAWAY_THRESH &
                   shared_bs$PP.H4_sceqtl > SCEQTL_THRESH)
    cat("  Cross-validated (Broadaway > 0.5 & sc-eQTL > 0.3):", both_bs, "\n\n")

    summary_rows[["Broadaway_vs_sceqtl"]] <- data.table(
      comparison       = "Broadaway_vs_sceqtl",
      gwas             = paste0("UKBB_ALT vs ", sceqtl_gwas_label),
      n_shared_genes   = n_shared_bs,
      spearman_rho     = round(rho_bs$estimate, 4),
      spearman_p       = rho_bs$p.value,
      concordance_rate = round(
        (sum(shared_bs$PP.H4_broadaway > STANDARD_THRESH &
             shared_bs$PP.H4_sceqtl > STANDARD_THRESH) +
         sum(shared_bs$PP.H4_broadaway <= STANDARD_THRESH &
             shared_bs$PP.H4_sceqtl <= STANDARD_THRESH)) / n_shared_bs, 4),
      both_sig         = sum(shared_bs$PP.H4_broadaway > STANDARD_THRESH &
                             shared_bs$PP.H4_sceqtl > STANDARD_THRESH),
      brd_only_sig     = sum(shared_bs$PP.H4_broadaway > STANDARD_THRESH &
                             shared_bs$PP.H4_sceqtl <= STANDARD_THRESH),
      gtex_only_sig    = sum(shared_bs$PP.H4_broadaway <= STANDARD_THRESH &
                             shared_bs$PP.H4_sceqtl > STANDARD_THRESH),
      neither_sig      = sum(shared_bs$PP.H4_broadaway <= STANDARD_THRESH &
                             shared_bs$PP.H4_sceqtl <= STANDARD_THRESH),
      cross_validated  = both_bs,
      fisher_OR        = NA_real_,
      fisher_p         = NA_real_,
      brd_thresh       = BROADAWAY_THRESH,
      gtex_thresh      = SCEQTL_THRESH
    )
  } else {
    cat("  Too few shared genes for correlation analysis\n\n")
  }
} else {
  cat("  Skipped (sc-eQTL results not available)\n\n")
}

# ==============================================================================
# 7. Pairwise Analyses: GTEx vs sc-eQTL
# ==============================================================================
cat("--- Step 5: Pairwise analysis — GTEx vs sc-eQTL ---\n")

if (has_gtex && has_sceqtl) {
  shared_gs <- merged[!is.na(PP.H4_gtex) & !is.na(PP.H4_sceqtl)]
  n_shared_gs <- nrow(shared_gs)
  cat("  Shared genes:", n_shared_gs, "\n")

  if (n_shared_gs >= 10) {
    rho_gs <- cor.test(shared_gs$PP.H4_gtex, shared_gs$PP.H4_sceqtl,
                       method = "spearman", exact = FALSE)
    cat("  Spearman rho (PP.H4):", round(rho_gs$estimate, 4),
        " p =", format(rho_gs$p.value, digits = 3), "\n")

    both_gs <- sum(shared_gs$PP.H4_gtex > GTEX_THRESH &
                   shared_gs$PP.H4_sceqtl > SCEQTL_THRESH)
    cat("  Cross-validated (GTEx > 0.3 & sc-eQTL > 0.3):", both_gs, "\n\n")

    summary_rows[["GTEx_vs_sceqtl"]] <- data.table(
      comparison       = "GTEx_vs_sceqtl",
      gwas             = paste0("UKBB_ALT vs ", sceqtl_gwas_label),
      n_shared_genes   = n_shared_gs,
      spearman_rho     = round(rho_gs$estimate, 4),
      spearman_p       = rho_gs$p.value,
      concordance_rate = round(
        (sum(shared_gs$PP.H4_gtex > STANDARD_THRESH &
             shared_gs$PP.H4_sceqtl > STANDARD_THRESH) +
         sum(shared_gs$PP.H4_gtex <= STANDARD_THRESH &
             shared_gs$PP.H4_sceqtl <= STANDARD_THRESH)) / n_shared_gs, 4),
      both_sig         = sum(shared_gs$PP.H4_gtex > STANDARD_THRESH &
                             shared_gs$PP.H4_sceqtl > STANDARD_THRESH),
      brd_only_sig     = sum(shared_gs$PP.H4_gtex > STANDARD_THRESH &
                             shared_gs$PP.H4_sceqtl <= STANDARD_THRESH),
      gtex_only_sig    = sum(shared_gs$PP.H4_gtex <= STANDARD_THRESH &
                             shared_gs$PP.H4_sceqtl > STANDARD_THRESH),
      neither_sig      = sum(shared_gs$PP.H4_gtex <= STANDARD_THRESH &
                             shared_gs$PP.H4_sceqtl <= STANDARD_THRESH),
      cross_validated  = both_gs,
      fisher_OR        = NA_real_,
      fisher_p         = NA_real_,
      brd_thresh       = GTEX_THRESH,
      gtex_thresh      = SCEQTL_THRESH
    )
  } else {
    cat("  Too few shared genes for correlation analysis\n\n")
  }
} else {
  cat("  Skipped (one or both sources not available)\n\n")
}

# ==============================================================================
# 8. Classify genes by eQTL source support
# ==============================================================================
cat("--- Step 6: Classifying genes by eQTL source support ---\n")

merged[, brd_sig  := !is.na(PP.H4_broadaway) & PP.H4_broadaway > BROADAWAY_THRESH]
if (has_gtex) {
  merged[, gtex_sig := !is.na(PP.H4_gtex) & PP.H4_gtex > GTEX_THRESH]
} else {
  merged[, gtex_sig := FALSE]
}
if (has_sceqtl) {
  merged[, sceqtl_sig := !is.na(PP.H4_sceqtl) & PP.H4_sceqtl > SCEQTL_THRESH]
} else {
  merged[, sceqtl_sig := FALSE]
}

# Count supporting eQTL sources
merged[, n_eqtl_sources := as.integer(brd_sig) + as.integer(gtex_sig) + as.integer(sceqtl_sig)]

# Classification
merged[, eqtl_class := fcase(
  n_eqtl_sources >= 3, "triple_validated",
  n_eqtl_sources == 2, "double_validated",
  n_eqtl_sources == 1 & brd_sig,  "broadaway_only",
  n_eqtl_sources == 1 & gtex_sig, "gtex_only",
  n_eqtl_sources == 1 & sceqtl_sig, "sceqtl_only",
  default = "none"
)]

cat("  Gene classification:\n")
class_counts <- merged[, .N, by = eqtl_class][order(-N)]
for (i in seq_len(nrow(class_counts))) {
  cat("    ", class_counts$eqtl_class[i], ":", class_counts$N[i], "\n")
}
cat("\n")

# ==============================================================================
# 9. Build validated-genes table
# ==============================================================================
cat("--- Step 7: Building validated-genes table ---\n")

validated <- merged[n_eqtl_sources >= 1][order(-n_eqtl_sources, -PP.H4_broadaway)]

# Select output columns
out_cols <- c("ensembl_clean", "symbol", "eqtl_class", "n_eqtl_sources",
              "PP.H4_broadaway", "PP.H3_broadaway", "n_snps_broadaway")
if (has_gtex) {
  out_cols <- c(out_cols, "PP.H4_gtex", "PP.H3_gtex", "n_snps_gtex")
}
if (has_sceqtl) {
  out_cols <- c(out_cols, "PP.H4_sceqtl", "best_cell_type", "sceqtl_gwas")
}
# Only keep columns that actually exist
out_cols <- intersect(out_cols, names(validated))
validated_out <- validated[, ..out_cols]

cat("  Validated genes (n_eqtl_sources >= 1):", nrow(validated_out), "\n")
cat("    Double-validated:", sum(validated_out$n_eqtl_sources >= 2), "\n")
cat("    Triple-validated:", sum(validated_out$n_eqtl_sources >= 3), "\n\n")

# ==============================================================================
# 10. Build full comparison table
# ==============================================================================
cat("--- Step 8: Building full comparison table ---\n")

comp_cols <- c("ensembl_clean", "symbol", "eqtl_class", "n_eqtl_sources",
               "PP.H4_broadaway", "PP.H3_broadaway", "n_snps_broadaway",
               "eqtl_N_broadaway")
if (has_gtex) {
  comp_cols <- c(comp_cols, "PP.H4_gtex", "PP.H3_gtex", "n_snps_gtex", "eqtl_N_gtex")
}
if (has_sceqtl) {
  comp_cols <- c(comp_cols, "PP.H4_sceqtl", "best_cell_type", "sceqtl_gwas")
}
comp_cols <- c(comp_cols, "brd_sig", "gtex_sig", "sceqtl_sig")
comp_cols <- intersect(comp_cols, names(merged))

comparison_out <- merged[, ..comp_cols][order(-n_eqtl_sources, -PP.H4_broadaway)]

# Alias to the column names the figure (figS_cross_eqtl.R) expects:
#   broadaway_pp4 / gtex_pp4 / sceqtl_pp4 (documented in that script's header).
# The PP.H4_* names are retained above for backward compatibility.
if ("PP.H4_broadaway" %in% names(comparison_out))
  comparison_out[, broadaway_pp4 := PP.H4_broadaway]
if ("PP.H4_gtex" %in% names(comparison_out))
  comparison_out[, gtex_pp4 := PP.H4_gtex]
if ("PP.H4_sceqtl" %in% names(comparison_out))
  comparison_out[, sceqtl_pp4 := PP.H4_sceqtl]

cat("  Full comparison table:", nrow(comparison_out), "genes\n\n")

# ==============================================================================
# 11. Compile summary statistics
# ==============================================================================
cat("--- Step 9: Compiling summary statistics ---\n")

summary_dt <- rbindlist(summary_rows, fill = TRUE)

# Add overall gene counts
overall <- data.table(
  comparison     = "overall_counts",
  gwas           = "mixed",
  n_shared_genes = nrow(merged),
  spearman_rho   = NA_real_,
  spearman_p     = NA_real_,
  concordance_rate = NA_real_,
  both_sig       = sum(merged$n_eqtl_sources >= 2, na.rm = TRUE),
  brd_only_sig   = sum(merged$eqtl_class == "broadaway_only", na.rm = TRUE),
  gtex_only_sig  = sum(merged$eqtl_class == "gtex_only", na.rm = TRUE),
  neither_sig    = sum(merged$eqtl_class == "none", na.rm = TRUE),
  cross_validated = sum(merged$n_eqtl_sources >= 2, na.rm = TRUE),
  fisher_OR      = NA_real_,
  fisher_p       = NA_real_,
  brd_thresh     = BROADAWAY_THRESH,
  gtex_thresh    = GTEX_THRESH
)

# Add source-level counts
source_counts <- data.table(
  comparison = "source_counts",
  gwas = "mixed",
  n_shared_genes = nrow(merged),
  spearman_rho = NA_real_,
  spearman_p = NA_real_,
  concordance_rate = NA_real_,
  both_sig = sum(merged$brd_sig, na.rm = TRUE),       # n Broadaway sig
  brd_only_sig = if (has_gtex) sum(merged$gtex_sig, na.rm = TRUE) else 0L,  # n GTEx sig
  gtex_only_sig = if (has_sceqtl) sum(merged$sceqtl_sig, na.rm = TRUE) else 0L,  # n sc-eQTL sig
  neither_sig = sum(merged$n_eqtl_sources >= 2, na.rm = TRUE),   # n multi-source
  cross_validated = sum(merged$n_eqtl_sources >= 3, na.rm = TRUE), # n triple
  fisher_OR = NA_real_,
  fisher_p = NA_real_,
  brd_thresh = BROADAWAY_THRESH,
  gtex_thresh = GTEX_THRESH
)

summary_dt <- rbindlist(list(summary_dt, overall, source_counts), fill = TRUE)
cat("  Summary rows:", nrow(summary_dt), "\n\n")

# ==============================================================================
# 12. Highlight known MASLD anchor genes
# ==============================================================================
cat("--- Step 10: Known MASLD anchor genes ---\n")

anchors <- c("PNPLA3", "TM6SF2", "HSD17B13", "MBOAT7", "MARC1", "GCKR",
             "THRB", "NR1H4", "PPARA", "DGAT2", "SLC39A8")
anchor_dt <- merged[symbol %in% anchors]
if (nrow(anchor_dt) > 0) {
  anchor_show <- anchor_dt[, .(symbol, eqtl_class, n_eqtl_sources)]
  if ("PP.H4_broadaway" %in% names(anchor_dt))
    anchor_show[, PP.H4_brd := round(anchor_dt$PP.H4_broadaway, 3)]
  if (has_gtex && "PP.H4_gtex" %in% names(anchor_dt))
    anchor_show[, PP.H4_gtex := round(anchor_dt$PP.H4_gtex, 3)]
  if (has_sceqtl && "PP.H4_sceqtl" %in% names(anchor_dt))
    anchor_show[, PP.H4_sceqtl := round(anchor_dt$PP.H4_sceqtl, 3)]

  cat("  Found", nrow(anchor_show), "of", length(anchors), "anchor genes:\n")
  print(anchor_show, topn = nrow(anchor_show))
} else {
  cat("  No anchor genes found in COLOC results\n")
}
cat("\n")

# ==============================================================================
# 13. Save outputs
# ==============================================================================
cat("--- Step 11: Saving outputs ---\n")

fwrite(comparison_out,
       file.path(OUTDIR, "cross_eqtl_comparison.csv"))
cat("  cross_eqtl_comparison.csv:", nrow(comparison_out), "genes\n")

fwrite(summary_dt,
       file.path(OUTDIR, "cross_eqtl_summary.csv"))
cat("  cross_eqtl_summary.csv:", nrow(summary_dt), "rows\n")

fwrite(validated_out,
       file.path(OUTDIR, "cross_eqtl_validated_genes.csv"))
cat("  cross_eqtl_validated_genes.csv:", nrow(validated_out), "genes\n")

cat("\n=== Script 35y complete ===\n")
cat("End time:", format(Sys.time()), "\n")
cat("Output directory:", OUTDIR, "\n")
