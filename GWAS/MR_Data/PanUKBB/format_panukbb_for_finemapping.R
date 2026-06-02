#!/usr/bin/env Rscript
# format_panukbb_for_finemapping.R
# ---------------------------------------------------------------------------
# Reformat Pan-UKBB AFR & CSA liver enzyme GWAS for the finemapping pipeline.
#
# The harmonised hg38 files (produced by format_panukbb_for_coloc.R) already
# carry the original hg19 coordinates in the `pos_hg19` column, so no
# liftOver is needed -- we simply select and rename columns.
#
# Input:  GWAS/MR_Data/PanUKBB/PanUKBB_{POP}_{TRAIT}_harmonised_hg38.tsv.gz
# Output: GWAS/finemapping/data/sumstats/PanUKBB_{POP}_{TRAIT}_reformatted_hg19.tsv
#
# Expected finemapping columns (tab-separated, uncompressed):
#   chromosome  position  allele1  allele2  beta  se  pval
#
# Usage: Rscript format_panukbb_for_finemapping.R
#        Or: POP=AFR TRAIT=ALT Rscript format_panukbb_for_finemapping.R
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

BASE_DIR <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
PANUKBB_DIR <- file.path(BASE_DIR, "GWAS/MR_Data/PanUKBB")
OUT_DIR     <- file.path(BASE_DIR, "GWAS/finemapping/data/sumstats")

# Optionally restrict to a single pop/trait
SEL_POP   <- Sys.getenv("POP",   "ALL")
SEL_TRAIT <- Sys.getenv("TRAIT", "ALL")

cat("=== Pan-UKBB GWAS → Finemapping Reformatting ===\n")
cat("Start time:", format(Sys.time()), "\n\n")

# ============================================================================
# Define populations and traits
# ============================================================================
populations <- c("AFR", "CSA")
traits      <- c("ALT", "AST", "GGT")

if (SEL_POP   != "ALL") populations <- SEL_POP
if (SEL_TRAIT != "ALL") traits      <- SEL_TRAIT

# Sample sizes (Pan-UKBB biomarker phenotypes, both sexes, IRNT)
pop_N <- list(
  AFR = list(ALT = 6636L,  AST = 6636L,  GGT = 6636L),
  CSA = list(ALT = 8876L,  AST = 8876L,  GGT = 8876L)
)

# ============================================================================
# Process each population x trait
# ============================================================================
for (pop in populations) {
  for (tr in traits) {
    infile  <- file.path(PANUKBB_DIR,
                         paste0("PanUKBB_", pop, "_", tr, "_harmonised_hg38.tsv.gz"))
    outfile <- file.path(OUT_DIR,
                         paste0("PanUKBB_", pop, "_", tr, "_reformatted_hg19.tsv"))

    if (!file.exists(infile)) {
      cat("  SKIP:", pop, tr, "- input not found:", basename(infile), "\n")
      next
    }

    n <- pop_N[[pop]][[tr]]
    cat("--- Processing Pan-UKBB", pop, tr, "(N ~", format(n, big.mark = ","), ") ---\n")

    dt <- fread(infile)
    cat("  Loaded:", format(nrow(dt), big.mark = ","), "variants\n")

    # Harmonised columns:
    #   chromosome, base_pair_location (hg38), pos_hg19,
    #   effect_allele, other_allele, beta, standard_error, p_value, EAF

    # Drop rows where pos_hg19 is missing (should not happen, but be safe)
    dt <- dt[!is.na(pos_hg19) & pos_hg19 > 0]

    # Basic QC: require valid beta/SE/p, autosomes only
    dt <- dt[!is.na(beta) & !is.na(standard_error) & standard_error > 0 & !is.na(p_value)]
    dt <- dt[!is.na(chromosome) & chromosome %in% 1:22]

    # Remove variants with zero effect (finemapping won't use them)
    # Keep beta == 0 though — some may be real nulls at tagged loci

    cat("  After QC:", format(nrow(dt), big.mark = ","), "variants\n")

    # Deduplicate by hg19 position (keep lowest p-value)
    dt <- dt[order(chromosome, pos_hg19, p_value)]
    dt <- dt[!duplicated(paste(chromosome, pos_hg19))]
    cat("  After dedup:", format(nrow(dt), big.mark = ","), "variants\n")

    # Format to finemapping standard: chromosome, position, allele1, allele2, beta, se, pval
    out_dt <- dt[, .(
      chromosome = chromosome,
      position   = pos_hg19,
      allele1    = effect_allele,
      allele2    = other_allele,
      beta       = beta,
      se         = standard_error,
      pval       = p_value
    )]

    # Sort by genomic position
    out_dt <- out_dt[order(chromosome, position)]

    n_sig <- sum(out_dt$pval < 5e-8)
    cat("  Genome-wide significant (p<5e-8):", n_sig, "\n")

    fwrite(out_dt, outfile, sep = "\t")
    cat("  Written:", outfile, "\n")
    cat("  Final:", format(nrow(out_dt), big.mark = ","), "variants\n\n")

    rm(dt, out_dt)
    gc(verbose = FALSE)
  }
}

cat("=== Reformatting complete ===\n")
cat("End time:", format(Sys.time()), "\n\n")
cat("Next steps:\n")
cat("  1. Identify lead SNPs:  sbatch src/01c_panukbb_lead_snps.sh\n")
cat("     (requires AFR/SAS PLINK LD panels for clumping)\n")
cat("  2. Add registry entries to config/gwas_registry.tsv\n")
cat("  3. Run finemapping:     sbatch src/02_run_susie.sh\n")
