#!/usr/bin/env Rscript
# format_magic_for_spredixcan.R
# ---------------------------------------------------------------------------
# Format MAGIC HOMA-IR for S-PrediXcan
#
# Input:  MAGIC_Manning_et_al_HOMAIR_MainEffect.txt.gz
# Format: Snp, maf, effect_allele, other_allele, MainEffects, MainSE, MainP,
#         BMIadjMainEffects, BMIadjMainSE, BMIadjMainP
# Output: MAGIC_HOMAIR_for_spredixcan.csv.gz
# ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(data.table))

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
MAGIC_DIR <- file.path(BASE, "GWAS/MR_Data/MAGIC")

infile  <- file.path(MAGIC_DIR, "MAGIC_Manning_et_al_HOMAIR_MainEffect.txt.gz")
outfile <- file.path(MAGIC_DIR, "MAGIC_HOMAIR_for_spredixcan.csv.gz")

cat("=== MAGIC HOMA-IR Formatting for S-PrediXcan ===\n")

dt <- fread(infile)
cat("  Raw rows:", format(nrow(dt), big.mark = ","), "\n")
cat("  Columns:", paste(names(dt), collapse = ", "), "\n")

# Use MainEffects/MainSE (unadjusted for BMI — captures full IR effect)
dt <- dt[!is.na(MainEffects) & !is.na(MainSE) & MainSE > 0]
dt <- dt[!is.na(Snp) & Snp != "" & grepl("^rs", Snp)]

# S-PrediXcan format (note: MAGIC Manning has no chr/pos — rsID-based matching)
out <- dt[, .(
  snp                 = Snp,
  non_effect_allele   = other_allele,
  effect_allele       = toupper(effect_allele),
  frequency           = maf,
  zscore              = MainEffects / MainSE,
  pvalue              = MainP,
  beta                = MainEffects,
  standard_error      = MainSE,
  sample_size         = 46186L  # Manning et al. 2012 / Dupuis et al. 2010
)]

# Fix allele case (MAGIC uses lowercase)
out[, non_effect_allele := toupper(non_effect_allele)]

out <- out[!duplicated(snp)]

cat("  Output rows:", format(nrow(out), big.mark = ","), "\n")
cat("  GW-sig (p<5e-8):", sum(out$pvalue < 5e-8, na.rm = TRUE), "\n")

fwrite(out, outfile)
cat("  Saved:", basename(outfile), "\n")
cat("Done.\n")
