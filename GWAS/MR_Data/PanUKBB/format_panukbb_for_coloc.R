#!/usr/bin/env Rscript
# format_panukbb_for_coloc.R
# ---------------------------------------------------------------------------
# Format Pan-UKBB multi-ancestry GWAS for COLOC
#
# Pan-UKBB files contain all ancestries in one file per phenotype.
# Columns: chr, pos (hg19), ref, alt, beta_{POP}, se_{POP}, neglog10_pval_{POP}
# We extract AFR and CSA ancestry-specific columns and liftover to hg38.
#
# Input:  GWAS/MR_Data/PanUKBB/biomarkers-{phenocode}-{trait}-both_sexes-irnt.tsv.bgz
# Output: GWAS/MR_Data/PanUKBB/PanUKBB_{POP}_{TRAIT}_harmonised_hg38.tsv.gz
#
# Usage: Rscript format_panukbb_for_coloc.R
#        Or: POP=AFR TRAIT=ALT Rscript format_panukbb_for_coloc.R
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(rtracklayer)
  library(GenomicRanges)
})

BASE_DIR   <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
PANUKBB_DIR <- file.path(BASE_DIR, "GWAS/MR_Data/PanUKBB")
CHAIN_FILE <- file.path(BASE_DIR, "data/broadaway_eqtl/hg19ToHg38.over.chain")
CHAIN_GZ   <- paste0(CHAIN_FILE, ".gz")

# Pop/trait selection (default: process all)
SEL_POP   <- Sys.getenv("POP", "ALL")
SEL_TRAIT <- Sys.getenv("TRAIT", "ALL")

cat("=== Pan-UKBB GWAS Formatting & Liftover ===\n")
cat("Start time:", format(Sys.time()), "\n\n")

# ============================================================================
# 1. Load chain file
# ============================================================================
cat("--- Loading liftover chain ---\n")
if (!file.exists(CHAIN_FILE)) {
  if (!file.exists(CHAIN_GZ)) {
    download.file(
      "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz",
      CHAIN_GZ, mode = "wb", quiet = FALSE
    )
  }
  system2("gunzip", args = c("-k", CHAIN_GZ))
}
chain <- import.chain(CHAIN_FILE)
cat("  Chain loaded:", length(chain), "chains\n\n")

# ============================================================================
# 2. Define traits and sample sizes
# ============================================================================
# Pan-UKBB sample sizes by ancestry (from Pan-UKBB docs)
# AFR: ~6,636; CSA: ~8,876 (for biomarker phenotypes)
pop_N <- list(
  AFR = list(ALT = 6636L, AST = 6636L, GGT = 6636L),
  CSA = list(ALT = 8876L, AST = 8876L, GGT = 8876L)
)

traits <- data.table(
  trait     = c("ALT", "AST", "GGT"),
  phenocode = c("30620", "30650", "30730")
)

populations <- c("AFR", "CSA")

if (SEL_POP != "ALL") populations <- SEL_POP
if (SEL_TRAIT != "ALL") traits <- traits[trait == SEL_TRAIT]

# ============================================================================
# 3. Process each population x trait
# ============================================================================
for (pop in populations) {
  for (i in seq_len(nrow(traits))) {
    tr <- traits$trait[i]
    pc <- traits$phenocode[i]
    n  <- pop_N[[pop]][[tr]]

    infile <- file.path(PANUKBB_DIR,
                        paste0("biomarkers-", pc, "-", tr, "-both_sexes-irnt.tsv.bgz"))
    outfile <- file.path(PANUKBB_DIR,
                         paste0("PanUKBB_", pop, "_", tr, "_harmonised_hg38.tsv.gz"))

    if (!file.exists(infile)) {
      cat("  SKIP:", pop, tr, "- input not found:", basename(infile), "\n")
      next
    }

    cat("--- Processing Pan-UKBB", pop, tr, "(N ~", format(n, big.mark = ","), ") ---\n")

    # Read full file (all ancestries)
    dt <- fread(infile)
    cat("  Raw rows:", format(nrow(dt), big.mark = ","), "\n")

    # Extract population-specific columns
    beta_col <- paste0("beta_", pop)
    se_col   <- paste0("se_", pop)
    pval_col <- paste0("neglog10_pval_", pop)
    af_col   <- paste0("af_", pop)
    lc_col   <- paste0("low_confidence_", pop)

    if (!beta_col %in% names(dt)) {
      cat("  SKIP:", pop, "columns not found in file\n")
      next
    }

    # Select and rename columns
    sub_dt <- dt[, .(
      chr      = as.integer(chr),
      pos_hg19 = as.integer(pos),
      ref, alt,
      beta     = get(beta_col),
      se       = get(se_col),
      neglog10_pval = get(pval_col),
      af       = get(af_col),
      low_confidence = get(lc_col)
    )]

    # Ensure numeric types (Pan-UKBB may have string "NA" for missing ancestries)
    sub_dt[, beta := as.numeric(beta)]
    sub_dt[, se := as.numeric(se)]
    sub_dt[, neglog10_pval := as.numeric(neglog10_pval)]
    sub_dt[, af := as.numeric(af)]

    # Filter: require valid beta/SE, exclude low-confidence, autosomes only
    sub_dt <- sub_dt[!is.na(beta) & !is.na(se) & se > 0 & !is.na(neglog10_pval)]
    sub_dt <- sub_dt[!is.na(chr) & chr %in% 1:22]
    sub_dt[, low_confidence := as.character(low_confidence)]
    sub_dt <- sub_dt[is.na(low_confidence) | low_confidence == "false"]

    # Convert neglog10 p-value to p-value
    sub_dt[, p_value := 10^(-neglog10_pval)]

    # MAF filter (>=0.01 for COLOC robustness)
    if (sum(!is.na(sub_dt$af)) > 0) {
      sub_dt[, maf := pmin(af, 1 - af)]
      n_before <- nrow(sub_dt)
      sub_dt <- sub_dt[is.na(maf) | maf >= 0.01]
      cat("  MAF filter (>=0.01):", n_before - nrow(sub_dt), "removed\n")
    }

    cat("  After filtering:", format(nrow(sub_dt), big.mark = ","), "variants\n")

    # Liftover hg19 → hg38
    cat("  Lifting over to GRCh38...\n")
    gr <- GRanges(
      seqnames = paste0("chr", sub_dt$chr),
      ranges = IRanges(start = sub_dt$pos_hg19, width = 1)
    )
    lifted <- liftOver(gr, chain)
    n_mapped <- lengths(lifted)

    sub_dt[, pos_hg38 := NA_integer_]
    unique_idx <- which(n_mapped == 1L)
    if (length(unique_idx) > 0) {
      sub_dt$pos_hg38[unique_idx] <- start(unlist(lifted[unique_idx]))
    }
    n_lifted <- sum(!is.na(sub_dt$pos_hg38))
    cat("  Liftover:", n_lifted, "/", nrow(sub_dt),
        "(", round(100 * n_lifted / nrow(sub_dt), 1), "%)\n")
    sub_dt <- sub_dt[!is.na(pos_hg38)]

    # Deduplicate
    sub_dt <- sub_dt[order(chr, pos_hg38, p_value)]
    sub_dt <- sub_dt[!duplicated(paste(chr, pos_hg38))]

    cat("  After dedup:", format(nrow(sub_dt), big.mark = ","), "\n")
    cat("  Genome-wide significant (p<5e-8):", sum(sub_dt$p_value < 5e-8), "\n")

    # Write harmonised output
    out_dt <- sub_dt[, .(
      chromosome = chr,
      base_pair_location = pos_hg38,
      pos_hg19,
      effect_allele = alt,
      other_allele = ref,
      beta,
      standard_error = se,
      p_value,
      EAF = af
    )]

    fwrite(out_dt, outfile, sep = "\t")
    cat("  Saved:", basename(outfile), "\n")

    # Free memory before next iteration
    rm(sub_dt, out_dt, gr, lifted)
    gc(verbose = FALSE)
    cat("\n")
  }
  # Free the large input data.table between phenocode files
  rm(dt)
  gc(verbose = FALSE)
}

cat("=== Pan-UKBB formatting complete ===\n")
cat("End time:", format(Sys.time()), "\n")
cat("\nNext steps:\n")
cat("  Run COLOC: GWAS_NAME=PANUKBB_AFR_ALT Rscript 50_panukbb_coloc.R\n")
