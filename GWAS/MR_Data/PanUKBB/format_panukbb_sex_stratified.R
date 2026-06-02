#!/usr/bin/env Rscript
# format_panukbb_sex_stratified.R
# ---------------------------------------------------------------------------
# Format Neale Lab Round 2 sex-stratified UKBB GWAS (ALT/AST/GGT, EUR) for
# both the existing COLOC pipeline (hg19 input) and the registry/atlas.
#
# Why Neale Lab Round 2 not Pan-UKBB v0.4?
#   Pan-UKBB v0.4 manifest contains only `both_sexes` rows for biomarkers
#   (verified 2026-05-11). Neale Round 2 publishes sex-stratified files
#   (female/male) for the same phenocodes — EUR-only, hg19. Matches PolyFun
#   EUR LD ref used by canonical SuSiE-COLOC.
#
# Input file format (Neale Lab Round 2 `varorder.tsv.bgz`):
#   variant (chr:pos:ref:alt hg19) | minor_allele | minor_AF | low_confidence
#   | n_complete_samples | AC | ytx | beta | se | tstat | pval
#
# Output (one per stratum, 6 total):
#   GWAS/finemapping/data/sumstats/PanUKBB_{F,M}_{ALT,AST,GGT}_reformatted_hg19.tsv
#   Columns: chromosome | position | allele1 | allele2 | beta | se | pval
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

BASE_DIR <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
IN_DIR  <- file.path(BASE_DIR, "GWAS/MR_Data/PanUKBB/sex_stratified")
OUT_DIR <- file.path(BASE_DIR, "GWAS/finemapping/data/sumstats")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("=== Neale Lab sex-stratified UKBB → finemapping reformatting ===\n")
cat("Start:", format(Sys.time()), "\n\n")

phenocodes <- list(ALT = "30620", AST = "30650", GGT = "30730")
sex_map    <- list(F = "female", M = "male")

summary_rows <- list()

for (trait in names(phenocodes)) {
  pc <- phenocodes[[trait]]
  for (sex_short in names(sex_map)) {
    sex_long <- sex_map[[sex_short]]
    infile <- file.path(IN_DIR,
                        sprintf("%s_irnt.gwas.imputed_v3.%s.varorder.tsv.bgz",
                                pc, sex_long))
    out_tag <- sprintf("PanUKBB_%s_%s", sex_short, trait)
    outfile <- file.path(OUT_DIR, paste0(out_tag, "_reformatted_hg19.tsv"))

    cat("--- ", out_tag, "---\n")
    if (!file.exists(infile)) {
      cat("  SKIP — input not found:", basename(infile), "\n\n")
      next
    }

    dt <- fread(infile, showProgress = FALSE)
    cat("  Raw rows:", format(nrow(dt), big.mark = ","), "\n")

    # Parse variant column: chr:pos:ref:alt (hg19)
    dt[, c("chr_str", "pos_str", "ref", "alt") := tstrsplit(variant, ":", fixed = TRUE)]
    dt[, chromosome := suppressWarnings(as.integer(chr_str))]
    dt[, position   := suppressWarnings(as.integer(pos_str))]

    # QC filters
    dt <- dt[!is.na(chromosome) & chromosome %in% 1:22]
    dt <- dt[!is.na(position)   & position > 0]
    dt[, beta := suppressWarnings(as.numeric(beta))]
    dt[, se   := suppressWarnings(as.numeric(se))]
    dt[, pval := suppressWarnings(as.numeric(pval))]
    dt <- dt[!is.na(beta) & !is.na(se) & se > 0 & !is.na(pval)]

    # low_confidence_variant filter (Neale lab QC flag)
    if ("low_confidence_variant" %in% names(dt)) {
      dt[, low_confidence_variant := as.character(low_confidence_variant)]
      n_before <- nrow(dt)
      dt <- dt[is.na(low_confidence_variant) | low_confidence_variant != "true"]
      cat("  Low-confidence variants removed:", n_before - nrow(dt), "\n")
    }

    # MAF >= 0.01 (COLOC robustness)
    if ("minor_AF" %in% names(dt)) {
      dt[, minor_AF := suppressWarnings(as.numeric(minor_AF))]
      n_before <- nrow(dt)
      dt <- dt[is.na(minor_AF) | minor_AF >= 0.01]
      cat("  MAF<0.01 removed:", n_before - nrow(dt), "\n")
    }

    cat("  After QC:", format(nrow(dt), big.mark = ","), "\n")

    # In Neale Round 2 `varorder` files, `beta` is on the alt allele.
    # The COLOC pipeline expects allele1 = effect_allele, allele2 = other_allele.
    # We set allele1 = alt (effect) and allele2 = ref (other).
    out_dt <- dt[, .(
      chromosome = chromosome,
      position   = position,
      allele1    = alt,
      allele2    = ref,
      beta       = beta,
      se         = se,
      pval       = pval
    )]

    # Deduplicate by position (keep best p-value)
    out_dt <- out_dt[order(chromosome, position, pval)]
    out_dt <- out_dt[!duplicated(paste(chromosome, position))]

    cat("  After dedup:", format(nrow(out_dt), big.mark = ","), "\n")
    cat("  Genome-wide sig (p<5e-8):", sum(out_dt$pval < 5e-8), "\n")

    # Sample size = mean of n_complete_samples (round-2 keeps it ~constant)
    n_eff <- if ("n_complete_samples" %in% names(dt)) {
      round(median(suppressWarnings(as.numeric(dt$n_complete_samples)), na.rm = TRUE))
    } else NA_integer_

    fwrite(out_dt, outfile, sep = "\t")
    cat("  Wrote:", outfile, "\n")
    cat("  N (median samples):", format(n_eff, big.mark = ","), "\n\n")

    summary_rows[[length(summary_rows) + 1L]] <- data.table(
      stratum  = out_tag,
      trait    = trait,
      sex      = sex_short,
      n_variants = nrow(out_dt),
      n_sig    = sum(out_dt$pval < 5e-8),
      n_samples_median = n_eff
    )

    rm(dt, out_dt); gc(verbose = FALSE)
  }
}

if (length(summary_rows) > 0) {
  summary_dt <- rbindlist(summary_rows)
  summary_file <- file.path(IN_DIR, "format_summary.tsv")
  fwrite(summary_dt, summary_file, sep = "\t")
  cat("\nSummary written:", summary_file, "\n")
  print(summary_dt)
}

cat("\n=== Reformatting complete ===\n")
cat("End:", format(Sys.time()), "\n")
