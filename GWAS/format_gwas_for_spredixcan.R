#!/usr/bin/env Rscript
# format_gwas_for_spredixcan.R
#
# Format GWAS summary statistics for S-PrediXcan (rsID-based matching).
# Output: per-GWAS CSV.gz files in GWAS/MR_Data/ ready for S-PrediXcan.
#
# Usage: Rscript format_gwas_for_spredixcan.R <gwas_name>
#   gwas_name: "ghodsian", "chen", or "namjou"
#
# Actual column names observed 2026-02-27:
#
#   Ghodsian 2021 (harmonised GWAS-SSF):
#     hm_variant_id, hm_rsid, hm_chrom, hm_pos, hm_other_allele, hm_effect_allele,
#     hm_beta, hm_odds_ratio, hm_ci_lower, hm_ci_upper, hm_effect_allele_frequency,
#     hm_code, variant_id, chromosome, base_pair_location, effect_allele, other_allele,
#     effect_allele_frequency, beta, standard_error, p_value, odds_ratio, ci_lower, ci_upper
#     => Use: hm_rsid, hm_effect_allele, hm_other_allele, hm_effect_allele_frequency,
#             hm_beta (harmonised), standard_error, p_value
#
#   Chen 2023 GOLDPlus (harmonised, but beta/se=1 — Z-score meta, TWAS only):
#     chromosome, base_pair_location, effect_allele, other_allele, beta, standard_error,
#     effect_allele_frequency, p_value, rsid, hm_coordinate_conversion, hm_code, variant_id
#     => Use: rsid, effect_allele, other_allele, effect_allele_frequency, beta, standard_error, p_value
#
#   Namjou 2019 eMERGE (plain GWAS, OR not beta, no SE):
#     SNP, A1, A2, CHR, BP, N, OR, P
#     => beta = log(OR); no SE available => SKIP_MR=TRUE for this GWAS

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript format_gwas_for_spredixcan.R <gwas_name>")
gwas_name <- tolower(args[1])
GWAS_DIR  <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/MR_Data"

cat("=== Formatting", gwas_name, "GWAS for S-PrediXcan ===\n")

if (gwas_name == "ghodsian") {
  raw <- fread(file.path(GWAS_DIR, "Ghodsian_2021_NAFLD_harmonised.tsv.gz"))
  cat("  Loaded", nrow(raw), "variants\n")
  cat("  Columns:", paste(names(raw)[1:12], collapse = ", "), "...\n")
  # Use harmonised columns (hm_*) for allele-flipped consistency.
  # hm_beta is the beta on the hm_effect_allele strand.
  # standard_error and p_value come from the original study columns (same values post-harmonisation).
  formatted <- raw[
    !is.na(hm_rsid) & grepl("^rs", hm_rsid) &
    !is.na(hm_beta) & !is.na(standard_error) & !is.na(p_value),
    .(
      snp              = hm_rsid,
      effect_allele    = hm_effect_allele,
      non_effect_allele = hm_other_allele,
      frequency        = hm_effect_allele_frequency,
      beta             = hm_beta,
      standard_error   = standard_error,
      pvalue           = p_value,
      sample_size      = 778614L
    )
  ]
  formatted <- unique(formatted, by = "snp")
  out_file  <- file.path(GWAS_DIR, "Ghodsian_2021_for_spredixcan.csv.gz")

} else if (gwas_name == "chen") {
  raw <- fread(file.path(GWAS_DIR, "Chen_2023_GOLDPlus_harmonised.tsv.gz"))
  cat("  Loaded", nrow(raw), "variants\n")
  cat("  Columns:", paste(names(raw), collapse = ", "), "\n")
  cat("  WARNING: Chen 2023 betas/SEs are Z-score placeholders (all = 1). TWAS only — MR will be skipped.\n")
  # Verify Chen 2023 uses Z-score placeholders (beta=SE=1)
  n_real_betas <- sum(raw$beta != 1, na.rm = TRUE)
  if (n_real_betas > 0) warning("Chen 2023: ", n_real_betas, " rows have beta != 1 — file may have been updated with real effect sizes. Review SKIP_MR flag.")
  # rsid column exists directly (no hm_ prefix in this file)
  formatted <- raw[
    !is.na(rsid) & grepl("^rs", rsid) & !is.na(p_value),
    .(
      snp              = rsid,
      effect_allele    = effect_allele,
      non_effect_allele = other_allele,
      frequency        = effect_allele_frequency,
      beta             = beta,           # placeholder = 1; required by S-PrediXcan but not real
      standard_error   = standard_error, # placeholder = 1
      pvalue           = p_value,
      sample_size      = 691479L
    )
  ]
  formatted <- unique(formatted, by = "snp")
  out_file  <- file.path(GWAS_DIR, "Chen_2023_for_spredixcan.csv.gz")

} else if (gwas_name == "namjou") {
  raw <- fread(file.path(GWAS_DIR, "Namjou_2019_NAFLD_GWAS.txt"))
  cat("  Loaded", nrow(raw), "variants\n")
  cat("  Columns:", paste(names(raw), collapse = ", "), "\n")
  cat("  NOTE: Namjou 2019 provides OR, not beta; no SE column present.\n")
  cat("  beta = log(OR); SE will be NA -> coloc will be skipped for this GWAS.\n")
  if (!"N" %in% names(raw)) stop("Namjou: expected column 'N' not found. Available columns: ", paste(names(raw), collapse=", "))
  # OR format: no SE. Included anyway for TWAS (MR ditched 2026-04-22).
  # S-PrediXcan requires beta+se; without SE we approximate se=NA and accept that
  # TWAS will fail gracefully. Primary use: catalogue entry and completeness.
  formatted <- raw[
    !is.na(SNP) & grepl("^rs", SNP) & !is.na(OR) & OR > 0 & !is.na(P),
    .(
      snp              = SNP,
      effect_allele    = A1,
      non_effect_allele = A2,
      frequency        = NA_real_,   # not available
      beta             = log(OR),
      standard_error   = NA_real_,   # not available; TWAS/MR will skip
      pvalue           = P,
      sample_size      = N           # per-row N (varies)
    )
  ]
  formatted <- unique(formatted, by = "snp")
  out_file  <- file.path(GWAS_DIR, "Namjou_2019_for_spredixcan.csv.gz")

} else {
  stop("Unknown GWAS: '", gwas_name, "'. Use: ghodsian, chen, or namjou")
}

cat("  Formatted", nrow(formatted), "SNPs with rsIDs\n")
cat("  AF available:", sum(!is.na(formatted$frequency)), "/", nrow(formatted), "\n")
cat("  SE available:", sum(!is.na(formatted$standard_error)), "/", nrow(formatted), "\n")
if (nrow(formatted) == 0) stop("GWAS '", gwas_name, "': formatted output has 0 rows after filtering. Check rsID column and input file.")
fwrite(formatted, out_file, sep = ",", compress = "gzip")
cat("  Saved:", out_file, "\n")
