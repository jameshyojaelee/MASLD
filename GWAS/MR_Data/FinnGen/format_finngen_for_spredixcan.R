#!/usr/bin/env Rscript
# format_finngen_for_spredixcan.R
# ---------------------------------------------------------------------------
# Format FinnGen R12 summary statistics for S-PrediXcan
#
# S-PrediXcan expects a CSV with columns:
#   panel_variant_id, chromosome, position, effect_allele,
#   non_effect_allele, zscore, pvalue, sample_size
#
# FinnGen R12 is GRCh38, which matches PredictDB hg38 models.
# The panel_variant_id uses chr_pos_ref_alt_b38 format since FinnGen
# rsIDs are semicolon-separated (multiple per position).
#
# Input:  GWAS/MR_Data/FinnGen/finngen_R12_{PHENOTYPE}.gz
# Output: GWAS/MR_Data/FinnGen/FinnGen_{PHENOTYPE}_for_spredixcan.csv.gz
#
# Usage: PHENOTYPE=NAFLD Rscript format_finngen_for_spredixcan.R
#        (or run without PHENOTYPE to process NAFLD and NASH)
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

BASE_DIR <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
FG_DIR   <- file.path(BASE_DIR, "GWAS/MR_Data/FinnGen")

PHENOTYPE <- Sys.getenv("PHENOTYPE", "ALL")

cat("=== FinnGen → S-PrediXcan Formatting ===\n")
cat("Start time:", format(Sys.time()), "\n\n")

# Define phenotypes and their FinnGen file names + sample sizes
phenotypes <- list(
  NAFLD = list(
    file     = "finngen_R12_NAFLD.gz",
    N_cases  = 4614L,
    N_controls = 434243L
  ),
  NASH = list(
    file     = "finngen_R12_CHIRHEP_NAS.gz",
    N_cases  = 1823L,
    N_controls = 434243L
  )
  # HCC excluded from MR/TWAS — different biology from MASLD
)

if (PHENOTYPE != "ALL") {
  if (!PHENOTYPE %in% names(phenotypes)) {
    stop("Unknown PHENOTYPE: ", PHENOTYPE, ". Use: ", paste(names(phenotypes), collapse = ", "))
  }
  phenotypes <- phenotypes[PHENOTYPE]
}

for (pheno_name in names(phenotypes)) {
  cfg <- phenotypes[[pheno_name]]
  infile <- file.path(FG_DIR, cfg$file)
  outfile <- file.path(FG_DIR, paste0("FinnGen_", pheno_name, "_for_spredixcan.csv.gz"))

  if (!file.exists(infile)) {
    cat("  SKIP:", pheno_name, "- file not found:", infile, "\n")
    next
  }

  cat("--- Processing FinnGen", pheno_name, "---\n")

  dt <- fread(infile)
  cat("  Raw rows:", nrow(dt), "\n")

  # Detect chrom column name
  chrom_col <- grep("chrom", names(dt), value = TRUE)[1]

  # Extract columns
  dt[, chr := as.integer(sub("^chr", "", get(chrom_col)))]
  dt <- dt[!is.na(chr) & chr %in% 1:22]
  dt <- dt[!is.na(beta) & !is.na(sebeta) & sebeta > 0]

  # Compute z-score
  dt[, zscore := beta / sebeta]

  # Total sample size
  N_total <- cfg$N_cases + cfg$N_controls

  # Build S-PrediXcan format
  out_dt <- dt[, .(
    panel_variant_id = paste(chr, pos, ref, alt, "b38", sep = "_"),
    chromosome       = chr,
    position         = as.integer(pos),
    effect_allele    = alt,
    non_effect_allele = ref,
    zscore           = zscore,
    pvalue           = pval,
    sample_size      = N_total
  )]

  # Deduplicate by position (keep lowest p-value)
  out_dt <- out_dt[order(chromosome, position, pvalue)]
  out_dt <- out_dt[!duplicated(paste(chromosome, position))]

  fwrite(out_dt, outfile)
  cat("  Output:", outfile, "\n")
  cat("  Variants:", nrow(out_dt), "\n")
  cat("  Sample size:", format(N_total, big.mark = ","),
      "(cases:", format(cfg$N_cases, big.mark = ","),
      "controls:", format(cfg$N_controls, big.mark = ","), ")\n\n")
}

cat("=== Formatting complete ===\n")
cat("End time:", format(Sys.time()), "\n")
cat("\nNext step: TWAS / COLOC pipeline (MR ditched 2026-04-22; Script 19 archived).\n")
