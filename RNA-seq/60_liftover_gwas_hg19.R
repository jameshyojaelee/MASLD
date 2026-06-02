#!/usr/bin/env Rscript
# 60_liftover_gwas_hg19.R
# ---------------------------------------------------------------------------
# Phase 0D: Liftover GWAS summary statistics from hg38 to hg19
#
# Required for: OTTERS (Broadaway TWAS), FOCUS, MAGMA (scDRS), LD-based
# clumping with 1000G EUR (hg19).
#
# Inputs:  GWAS/MR_Data/*.tsv.gz (hg38 coordinates)
# Outputs: GWAS/MR_Data/hg19/*.tsv.gz (hg19 coordinates)
#
# Downloads hg38ToHg19.over.chain from UCSC if not present.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(rtracklayer)
  library(GenomicRanges)
})

# ============================================================================
# Configuration
# ============================================================================
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

GWAS_DIR <- file.path(BASE, "GWAS/MR_Data")
OUT_DIR  <- file.path(GWAS_DIR, "hg19")
CHAIN_DIR <- file.path(BASE, "data/broadaway_eqtl")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# Get chain file (hg38 → hg19)
# ============================================================================
chain_file <- file.path(CHAIN_DIR, "hg38ToHg19.over.chain")
if (!file.exists(chain_file)) {
  chain_gz <- paste0(chain_file, ".gz")
  if (!file.exists(chain_gz)) {
    cat("Downloading hg38ToHg19.over.chain.gz from UCSC...\n")
    download.file(
      "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz",
      chain_gz, mode = "wb", quiet = TRUE
    )
  }
  cat("Decompressing chain file...\n")
  system(paste("gunzip -k", chain_gz))
}
stopifnot(file.exists(chain_file))
chain <- import.chain(chain_file)
cat("Chain file loaded\n")

# ============================================================================
# GWAS file definitions (hg38 → hg19)
# ============================================================================
gwas_defs <- list(
  list(
    name = "ghodsian",
    file = "Ghodsian_2021_NAFLD_harmonised.tsv.gz",
    chr_col = "hm_chrom", pos_col = "hm_pos",
    cols_keep = c("hm_variant_id", "hm_rsid", "hm_chrom", "hm_pos",
                  "hm_other_allele", "hm_effect_allele", "hm_beta",
                  "hm_effect_allele_frequency", "standard_error", "p_value")
  ),
  list(
    name = "ukbb_alt",
    file = "GCST90019492_UKBB_ALT_harmonised.tsv.gz",
    chr_col = "chromosome", pos_col = "base_pair_location",
    cols_keep = c("chromosome", "base_pair_location", "effect_allele",
                  "other_allele", "beta", "standard_error",
                  "effect_allele_frequency", "p_value", "variant_id", "rsid")
  ),
  list(
    name = "ukbb_ast",
    file = "GCST90019497_UKBB_AST_harmonised.tsv.gz",
    chr_col = "chromosome", pos_col = "base_pair_location",
    cols_keep = c("chromosome", "base_pair_location", "effect_allele",
                  "other_allele", "beta", "standard_error",
                  "effect_allele_frequency", "p_value", "variant_id", "rsid")
  ),
  list(
    name = "ukbb_ggt",
    file = "GCST90019507_UKBB_GGT_harmonised.tsv.gz",
    chr_col = "chromosome", pos_col = "base_pair_location",
    cols_keep = c("chromosome", "base_pair_location", "effect_allele",
                  "other_allele", "beta", "standard_error",
                  "effect_allele_frequency", "p_value", "variant_id", "rsid")
  ),
  list(
    name = "pdff",
    file = "GCST90267352_PDFF_Pazoki2022.tsv.gz",
    chr_col = "chromosome", pos_col = "base_pair_location",
    cols_keep = c("chromosome", "base_pair_location", "effect_allele",
                  "other_allele", "beta", "standard_error",
                  "effect_allele_frequency", "p_value", "variant_id")
  )
)

# ============================================================================
# Liftover function
# ============================================================================
liftover_gwas <- function(def) {
  in_file  <- file.path(GWAS_DIR, def$file)
  out_file <- file.path(OUT_DIR, gsub("\\.tsv\\.gz$", "_hg19.tsv.gz", def$file))

  if (file.exists(out_file)) {
    cat("  Already exists:", basename(out_file), "\n")
    return(invisible(NULL))
  }

  if (!file.exists(in_file)) {
    cat("  SKIP (file not found):", def$file, "\n")
    return(invisible(NULL))
  }

  cat("  Reading", def$name, "...\n")
  dt <- fread(in_file)
  n_start <- nrow(dt)
  cat("    Rows:", format(n_start, big.mark = ","), "\n")

  # Extract chr/pos — filter NAs before GRanges creation
  chr_vals <- dt[[def$chr_col]]
  pos_vals <- as.integer(dt[[def$pos_col]])

  # Remove rows with NA chr or pos
  valid_idx <- which(!is.na(chr_vals) & !is.na(pos_vals) & chr_vals != "" & pos_vals > 0)
  if (length(valid_idx) < nrow(dt)) {
    cat("    Filtered", nrow(dt) - length(valid_idx), "rows with NA chr/pos\n")
    dt <- dt[valid_idx]
    chr_vals <- chr_vals[valid_idx]
    pos_vals <- pos_vals[valid_idx]
  }

  # Add chr prefix if missing
  chr_str <- ifelse(grepl("^chr", chr_vals), chr_vals, paste0("chr", chr_vals))

  # Create GRanges (hg38)
  gr <- GRanges(seqnames = chr_str,
                ranges   = IRanges(start = pos_vals, width = 1))

  # Liftover to hg19
  cat("    Lifting over...\n")
  lifted <- liftOver(gr, chain)

  # Get successfully mapped indices
  n_mapped <- lengths(lifted)
  good_idx <- which(n_mapped == 1)  # only keep unambiguous 1:1 mappings

  cat("    Mapped:", format(length(good_idx), big.mark = ","), "/",
      format(n_start, big.mark = ","),
      "(", round(100 * length(good_idx) / n_start, 1), "%)\n")

  # Extract new coordinates
  new_gr  <- unlist(lifted[good_idx])
  new_chr <- as.integer(gsub("chr", "", as.character(seqnames(new_gr))))
  new_pos <- start(new_gr)

  # Subset original data and add hg19 coordinates
  dt_out <- dt[good_idx]
  dt_out[, hg19_chr := new_chr]
  dt_out[, hg19_pos := new_pos]

  # Write output
  cat("    Writing", basename(out_file), "...\n")
  fwrite(dt_out, out_file, sep = "\t", compress = "gzip")
  cat("    Done:", format(nrow(dt_out), big.mark = ","), "variants\n")
}

# ============================================================================
# Run all liftoversa
# ============================================================================
cat("=== Liftover GWAS hg38 → hg19 ===\n")
cat("Output directory:", OUT_DIR, "\n\n")

for (def in gwas_defs) {
  cat("Processing:", def$name, "\n")
  tryCatch(
    liftover_gwas(def),
    error = function(e) cat("  ERROR:", conditionMessage(e), "\n")
  )
  cat("\n")
}

cat("=== All liftover complete ===\n")
