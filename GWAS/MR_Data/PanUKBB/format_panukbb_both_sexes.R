#!/usr/bin/env Rscript
# format_panukbb_both_sexes.R
# ---------------------------------------------------------------------------
# Format Neale Lab Round 2 BOTH-SEXES UKBB GWAS (ALT/AST/GGT, EUR) as the
# matched baseline for the B1 sex-stratified COLOC chain.
#
# Rationale: F-only and M-only Neale R2 are used in B1 sex-strat COLOC, but
# the canonical `UKBB_ALT/AST/GGT` (GWAS Catalog harmonised, Pazoki) is a
# DIFFERENT release. To attribute differences to sex (not release), need a
# matched-release both-sexes baseline.
#
# Input format (identical to F/M Neale R2):
#   variant (chr:pos:ref:alt hg19) | minor_allele | minor_AF | low_confidence
#   | n_complete_samples | AC | ytx | beta | se | tstat | pval
#
# Output: GWAS/finemapping/data/sumstats/PanUKBB_BS_{ALT,AST,GGT}_reformatted_hg19.tsv
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({ library(data.table) })

BASE_DIR <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
IN_DIR  <- file.path(BASE_DIR, "GWAS/MR_Data/PanUKBB/sex_stratified")
OUT_DIR <- file.path(BASE_DIR, "GWAS/finemapping/data/sumstats")
LEAD_DIR <- file.path(BASE_DIR, "GWAS/finemapping/data/lead_snps")
REG_FILE <- file.path(BASE_DIR, "GWAS/finemapping/config/gwas_registry.tsv")  # FIX (review A01#5): was data/ (non-existent) -> crashed; canonical registry is config/

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LEAD_DIR, recursive = TRUE, showWarnings = FALSE)

cat("=== Neale R2 BOTH_SEXES UKBB → finemapping reformatting ===\n")
cat("Start:", format(Sys.time()), "\n\n")

phenocodes <- list(ALT = "30620", AST = "30650", GGT = "30730")
summary_rows <- list()

for (trait in names(phenocodes)) {
  pc <- phenocodes[[trait]]
  infile <- file.path(IN_DIR,
                      sprintf("%s_irnt.gwas.imputed_v3.both_sexes.varorder.tsv.bgz", pc))
  out_tag <- sprintf("PanUKBB_BS_%s", trait)
  outfile <- file.path(OUT_DIR, paste0(out_tag, "_reformatted_hg19.tsv"))

  cat("--- ", out_tag, "---\n")
  if (!file.exists(infile)) { cat("  SKIP — input not found:", basename(infile), "\n\n"); next }

  dt <- fread(infile, showProgress = FALSE)
  cat("  Raw rows:", format(nrow(dt), big.mark = ","), "\n")

  dt[, c("chr_str","pos_str","ref","alt") := tstrsplit(variant, ":", fixed = TRUE)]
  dt[, chromosome := suppressWarnings(as.integer(chr_str))]
  dt[, position   := suppressWarnings(as.integer(pos_str))]

  dt <- dt[!is.na(chromosome) & chromosome %in% 1:22]
  dt <- dt[!is.na(position) & position > 0]
  dt[, beta := suppressWarnings(as.numeric(beta))]
  dt[, se   := suppressWarnings(as.numeric(se))]
  dt[, pval := suppressWarnings(as.numeric(pval))]
  dt <- dt[!is.na(beta) & !is.na(se) & se > 0 & !is.na(pval)]

  if ("low_confidence_variant" %in% names(dt)) {
    dt[, low_confidence_variant := as.character(low_confidence_variant)]
    n_before <- nrow(dt)
    dt <- dt[is.na(low_confidence_variant) | low_confidence_variant != "true"]
    cat("  Low-confidence variants removed:", n_before - nrow(dt), "\n")
  }

  if ("minor_AF" %in% names(dt)) {
    dt[, minor_AF := suppressWarnings(as.numeric(minor_AF))]
    n_before <- nrow(dt)
    dt <- dt[is.na(minor_AF) | minor_AF >= 0.01]
    cat("  MAF<0.01 removed:", n_before - nrow(dt), "\n")
  }

  cat("  After QC:", format(nrow(dt), big.mark = ","), "\n")

  out_dt <- dt[, .(
    chromosome = chromosome, position = position,
    allele1 = alt, allele2 = ref,
    beta = beta, se = se, pval = pval
  )]
  out_dt <- out_dt[order(chromosome, position, pval)]
  out_dt <- out_dt[!duplicated(paste(chromosome, position))]

  cat("  After dedup:", format(nrow(out_dt), big.mark = ","), "\n")
  cat("  Genome-wide sig (p<5e-8):", sum(out_dt$pval < 5e-8), "\n")

  n_eff <- if ("n_complete_samples" %in% names(dt)) {
    round(median(suppressWarnings(as.numeric(dt$n_complete_samples)), na.rm = TRUE))
  } else NA_integer_

  fwrite(out_dt, outfile, sep = "\t")
  cat("  Wrote:", outfile, "\n")
  cat("  N (median samples):", format(n_eff, big.mark = ","), "\n\n")

  # Build lead SNP file (p<5e-8) for finemapping pipeline
  lead_dt <- out_dt[pval < 5e-8, .(chromosome, position, allele1, allele2, beta, se, pval)]
  if (nrow(lead_dt) > 0) {
    lead_file <- file.path(LEAD_DIR, paste0(out_tag, "_leadSNPs.tsv"))
    fwrite(lead_dt, lead_file, sep = "\t")
    cat("  Lead SNPs written:", lead_file, "(", nrow(lead_dt), "rows)\n\n")
  }

  summary_rows[[length(summary_rows) + 1L]] <- data.table(
    stratum = out_tag, trait = trait, sex = "BS",
    n_variants = nrow(out_dt), n_sig = sum(out_dt$pval < 5e-8),
    n_samples_median = n_eff
  )
  rm(dt, out_dt); gc(verbose = FALSE)
}

# Append registry rows (skip if already present)
if (length(summary_rows) > 0) {
  summary_dt <- rbindlist(summary_rows)
  fwrite(summary_dt, file.path(IN_DIR, "format_both_sexes_summary.tsv"), sep = "\t")
  print(summary_dt)

  reg <- fread(REG_FILE)
  new_rows <- list()
  for (i in seq_len(nrow(summary_dt))) {
    tag <- summary_dt$stratum[i]
    if (!tag %in% reg$study_name) {  # FIX (review A01#5): registry key is study_name, not gwas_name
      n <- summary_dt$n_samples_median[i]
      # Column names MUST match the registry schema exactly, else rbindlist(fill=TRUE)
      # corrupts the registry with a union of mismatched columns.
      new_rows[[length(new_rows) + 1L]] <- data.table(
        study_name = tag,
        sumstats_path = sprintf("data/sumstats/%s_reformatted_hg19.tsv", tag),
        leadsnps_path = sprintf("data/lead_snps/%s_leadSNPs.tsv", tag),
        ancestry = "EUR", trait_type = "quantitative",
        N_tot = n, N_cases = 0,
        ld_panel = "ukbb_eur", window_mb = 0.5
      )
      cat("  Appending registry:", tag, "(N=", n, ")\n")
    } else {
      cat("  Registry already has:", tag, "\n")
    }
  }
  if (length(new_rows) > 0) {
    reg <- rbindlist(list(reg, rbindlist(new_rows)), use.names = TRUE, fill = TRUE)
    fwrite(reg, REG_FILE, sep = "\t")
    cat("\n  Registry updated:", REG_FILE, "\n")
  }
}

cat("\n=== Done:", format(Sys.time()), "===\n")
