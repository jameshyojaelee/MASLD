#!/usr/bin/env Rscript
# Harmonize Sveinbjornsson 2022 NAFLD GWAS summary stats.
# 2026-04-22: Previously produced the input format for the (now archived)
# Script 19 MR driver. MR was ditched from the paper; output is now consumed
# by the TWAS / COLOC pipeline.
#
# Output columns:
#   rsid, p_value, lnOR, standard_error, effect_allele_frequency,
#   effect_allele, other_allele, hm_chrom, hm_pos
#
# Usage: Rscript harmonize_sveinbjornsson_gwas.R

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

INFILE  <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/MR_Data/GCST90294579-buildGRCh38.tsv.gz"
OUTFILE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/MR_Data/sveinbjornsson_2022_nafld_gwas_harmonized.tsv.gz"

cat("Reading GWAS file...\n")
gwas <- fread(INFILE)
cat("Columns:", paste(names(gwas), collapse = ", "), "\n")
cat("Rows:", nrow(gwas), "\n")

# --- Detect and standardize column names ---
# GWAS Catalog harmonized files typically have:
#   hm_rsid, hm_chrom, hm_pos, hm_effect_allele, hm_other_allele,
#   hm_beta (log OR for binary; beta for continuous), hm_odds_ratio,
#   standard_error, p_value, hm_effect_allele_frequency
#
# Map to script 19 expected names:
col_map <- list(
  rsid                   = c("hm_rsid", "rsid", "SNP", "snp", "variant_id"),
  p_value                = c("p_value", "p-value", "pvalue", "P", "p"),
  lnOR                   = c("hm_beta", "beta", "lnOR", "log_odds_ratio"),
  standard_error         = c("standard_error", "se", "SE"),
  effect_allele_frequency = c("hm_effect_allele_frequency", "effect_allele_frequency", "eaf", "EAF"),
  effect_allele          = c("hm_effect_allele", "effect_allele", "A1"),
  other_allele           = c("hm_other_allele", "other_allele", "A2"),
  hm_chrom               = c("hm_chrom", "chromosome", "chr", "CHR"),
  hm_pos                 = c("hm_pos", "base_pair_location", "position", "POS", "BP")
)

find_col <- function(df, candidates) {
  for (col in candidates) if (col %in% names(df)) return(col)
  return(NA_character_)
}

rename_vec <- sapply(names(col_map), function(target) {
  src <- find_col(gwas, col_map[[target]])
  if (is.na(src)) {
    cat("WARNING: Could not find column for", target, ". Candidates:", paste(col_map[[target]], collapse=","), "\n")
  }
  src
})

cat("\nColumn mapping:\n")
print(rename_vec)

# Keep only columns we can map
valid <- !is.na(rename_vec)
gwas_out <- gwas[, rename_vec[valid], with = FALSE]
setnames(gwas_out, rename_vec[valid], names(rename_vec)[valid])

# If hm_odds_ratio present and lnOR is missing or 1.0, compute log(OR)
if ("hm_odds_ratio" %in% names(gwas) && (!"lnOR" %in% names(gwas_out) || all(is.na(gwas_out$lnOR)))) {
  cat("Computing lnOR from hm_odds_ratio\n")
  gwas_out$lnOR <- log(gwas[["hm_odds_ratio"]])
}

# Basic filters
gwas_out <- gwas_out[!is.na(gwas_out$rsid) & gwas_out$rsid != "", ]
gwas_out <- gwas_out[!is.na(gwas_out$p_value), ]

cat("\nOutput preview:\n")
print(head(gwas_out))
cat("Output rows:", nrow(gwas_out), "\n")

# Write output
fwrite(gwas_out, OUTFILE)
cat("Written to:", OUTFILE, "\n")

cat("\n--- Next steps ---\n")
cat("Feed formatted GWAS to the TWAS / COLOC pipeline.\n")
cat("  ANSTEE_GWAS_PATH  <- '", OUTFILE, "'\n", sep="")
cat("  ANSTEE_N_TOTAL    <- <check paper for total N>\n")
cat("  ANSTEE_N_CASES    <- <check paper for case N>\n")
# 2026-04-22: MR ditched from paper; the Script 19 constants this block used
# to populate are no longer consumed. Output is now used by TWAS / COLOC.
