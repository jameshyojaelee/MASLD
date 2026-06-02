#!/usr/bin/env Rscript
# format_glgc_for_spredixcan.R
# ---------------------------------------------------------------------------
# Format GLGC 2021 trans-ancestry Triglycerides for S-PrediXcan
#
# Input:  with_BF_meta-analysis_AFR_EAS_EUR_HIS_SAS_logTG_INV_ALL_with_N_1.gz
# Format: rsID, CHROM, POS_b37, REF, ALT, N, N_studies, POOLED_ALT_AF,
#         pvalue_neg_log10, pvalue, lnBF, ..., METAL_Effect, METAL_StdErr, METAL_Pvalue
# Output: GLGC_TG_for_spredixcan.csv.gz
# ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(data.table))

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
GLGC_DIR <- file.path(BASE, "GWAS/MR_Data/GLGC")

infile  <- file.path(GLGC_DIR, "with_BF_meta-analysis_AFR_EAS_EUR_HIS_SAS_logTG_INV_ALL_with_N_1.gz")
outfile <- file.path(GLGC_DIR, "GLGC_TG_for_spredixcan.csv.gz")

cat("=== GLGC TG Formatting for S-PrediXcan ===\n")

dt <- fread(infile)
cat("  Raw rows:", format(nrow(dt), big.mark = ","), "\n")
cat("  Columns:", paste(names(dt), collapse = ", "), "\n")

# Use METAL_Effect/METAL_StdErr (fixed-effects meta-analysis)
dt <- dt[!is.na(METAL_Effect) & !is.na(METAL_StdErr) & METAL_StdErr > 0]
dt <- dt[!is.na(rsID) & rsID != "" & grepl("^rs", rsID)]

# S-PrediXcan format
out <- dt[, .(
  snp                 = rsID,
  chromosome          = CHROM,
  position            = POS_b37,
  non_effect_allele   = REF,
  effect_allele       = ALT,
  frequency           = POOLED_ALT_AF,
  zscore              = METAL_Effect / METAL_StdErr,
  pvalue              = METAL_Pvalue,
  beta                = METAL_Effect,
  standard_error      = METAL_StdErr,
  sample_size         = N
)]

# Remove duplicates (keep first by rsID)
out <- out[!duplicated(snp)]

cat("  Output rows:", format(nrow(out), big.mark = ","), "\n")
cat("  Median N:", format(median(out$sample_size, na.rm = TRUE), big.mark = ","), "\n")
cat("  GW-sig (p<5e-8):", sum(out$pvalue < 5e-8, na.rm = TRUE), "\n")

fwrite(out, outfile)
cat("  Saved:", basename(outfile), "\n")
cat("Done.\n")
