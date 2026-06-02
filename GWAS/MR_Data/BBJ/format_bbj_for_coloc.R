#!/usr/bin/env Rscript
# format_bbj_for_coloc.R
# ---------------------------------------------------------------------------
# Format Biobank Japan (BBJ) liver enzyme GWAS for COLOC and S-PrediXcan
#
# BBJ GWAS (Kanai et al. 2018) are in GRCh37 (hg19) and need:
#   1. Liftover to GRCh38 (to match Broadaway eQTLs post-liftover)
#   2. Column harmonization to project standard format
#   3. Optional rsID annotation (for S-PrediXcan)
#
# Input:  GWAS/MR_Data/BBJ/raw/BBJ_{ALT,AST,GGT}_raw.tsv.gz
# Output: GWAS/MR_Data/BBJ/BBJ_{ALT,AST,GGT}_harmonised_hg38.tsv.gz
#         GWAS/MR_Data/BBJ/BBJ_{ALT,AST,GGT}_for_spredixcan.csv.gz
#
# Usage: Rscript format_bbj_for_coloc.R
#        Or: TRAIT=ALT Rscript format_bbj_for_coloc.R
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(rtracklayer)
  library(GenomicRanges)
})

BASE_DIR   <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
BBJ_DIR    <- file.path(BASE_DIR, "GWAS/MR_Data/BBJ")
RAW_DIR    <- file.path(BBJ_DIR, "raw")
CHAIN_FILE <- file.path(BASE_DIR, "data/broadaway_eqtl/hg19ToHg38.over.chain")
CHAIN_GZ   <- paste0(CHAIN_FILE, ".gz")

# Optionally process a single trait
TRAIT <- Sys.getenv("TRAIT", "ALL")

cat("=== BBJ GWAS Formatting & Liftover ===\n")
cat("Start time:", format(Sys.time()), "\n\n")

# ============================================================================
# 1. Load chain file
# ============================================================================
cat("--- Loading liftover chain ---\n")
if (!file.exists(CHAIN_FILE)) {
  if (!file.exists(CHAIN_GZ)) {
    cat("  Downloading hg19ToHg38.over.chain.gz...\n")
    download.file(
      "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz",
      CHAIN_GZ, mode = "wb", quiet = FALSE
    )
  }
  cat("  Decompressing chain file...\n")
  system2("gunzip", args = c("-k", CHAIN_GZ))
}
chain <- import.chain(CHAIN_FILE)
cat("  Chain loaded:", length(chain), "chains\n\n")

# ============================================================================
# 2. Define BBJ traits and sample sizes
# ============================================================================
traits <- data.table(
  trait = c("ALT", "AST", "GGT"),
  N     = c(261406L, 261270L, 164006L)
)

if (TRAIT != "ALL") {
  traits <- traits[trait == TRAIT]
  if (nrow(traits) == 0) stop("Unknown TRAIT: ", TRAIT)
}

# ============================================================================
# 3. Process each trait
# ============================================================================
for (i in seq_len(nrow(traits))) {
  tr <- traits$trait[i]
  n  <- traits$N[i]

  zip_file <- file.path(RAW_DIR, paste0("BBJ_", tr, "_raw.tsv.gz"))
  out_coloc <- file.path(BBJ_DIR, paste0("BBJ_", tr, "_harmonised_hg38.tsv.gz"))
  out_spredixcan <- file.path(BBJ_DIR, paste0("BBJ_", tr, "_for_spredixcan.csv.gz"))

  if (!file.exists(zip_file)) {
    cat("  SKIP:", tr, "- raw file not found:", zip_file, "\n")
    cat("    Run download_bbj_liver_enzymes.sh first\n\n")
    next
  }

  cat("--- Processing BBJ", tr, "(N=", format(n, big.mark = ","), ") ---\n")

  # BBJ files are ZIP archives containing gzipped BOLT-LMM output
  # Extract to raw directory if not already done
  inner_dir <- file.path(RAW_DIR, paste0("hum0197.v3.BBJ.", tr, ".v1"))
  inner_auto <- file.path(inner_dir, paste0("GWASsummary_", tr, "_Japanese_SakaueKanai2020.auto.txt.gz"))

  if (!file.exists(inner_auto)) {
    cat("  Extracting ZIP archive...\n")
    unzip(zip_file, exdir = RAW_DIR)
  }

  if (!file.exists(inner_auto)) {
    cat("  ERROR: Expected inner file not found:", inner_auto, "\n")
    # Try to find it
    extracted <- list.files(RAW_DIR, pattern = paste0(".*", tr, ".*auto.*\\.gz$"), recursive = TRUE, full.names = TRUE)
    if (length(extracted) > 0) {
      inner_auto <- extracted[1]
      cat("  Found at:", inner_auto, "\n")
    } else {
      cat("  Skipping", tr, "\n\n")
      next
    }
  }

  # Load autosomal data (skip chrX — COLOC is autosomal only)
  dt <- fread(inner_auto)
  cat("  Raw rows:", nrow(dt), "\n")
  cat("  Columns:", paste(names(dt), collapse = ", "), "\n")

  # ========================================================================
  # Auto-detect column format
  # PheWeb format: chrom, pos, ref, alt, pval, beta, sebeta, ...
  # GWAS Catalog format: chromosome, base_pair_location, ...
  # JENGER format may differ
  # ========================================================================

  # Standardize column names
  # BBJ BOLT-LMM format has: SNP, CHR, BP, ALLELE1 (effect), ALLELE0 (other),
  # A1FREQ, BETA, SE, P_BOLT_LMM_INF (preferred p-value)
  col_map <- list(
    chr = c("CHR", "#chrom", "chrom", "chromosome", "chr"),
    pos = c("BP", "pos", "base_pair_location", "position", "POS"),
    ref = c("ALLELE0", "ref", "other_allele", "REF", "non_effect_allele", "A2"),
    alt = c("ALLELE1", "alt", "effect_allele", "ALT", "EA", "A1"),
    pval = c("P_BOLT_LMM_INF", "P_LINREG", "pval", "p_value", "P", "p", "PVAL"),
    beta = c("BETA", "beta", "b", "effect"),
    se = c("SE", "sebeta", "standard_error", "se", "stderr"),
    maf = c("A1FREQ", "af_alt", "maf", "MAF", "eaf", "EAF", "freq"),
    rsid = c("SNP", "rsids", "rsid", "variant_id", "MarkerName")
  )

  detected <- list()
  for (field in names(col_map)) {
    match <- intersect(col_map[[field]], names(dt))
    if (length(match) > 0) {
      detected[[field]] <- match[1]
    } else {
      detected[[field]] <- NA
    }
  }

  cat("  Detected columns:\n")
  for (field in names(detected)) {
    cat("    ", field, "->", ifelse(is.na(detected[[field]]), "NOT FOUND", detected[[field]]), "\n")
  }

  # Validate required columns exist
  required <- c("chr", "pos", "ref", "alt", "pval", "beta", "se")
  missing <- required[sapply(required, function(f) is.na(detected[[f]]))]
  if (length(missing) > 0) {
    cat("  ERROR: Missing required columns:", paste(missing, collapse = ", "), "\n")
    cat("    Available columns:", paste(names(dt), collapse = ", "), "\n")
    cat("    Skipping", tr, "\n\n")
    next
  }

  # Rename to standard
  setnames(dt, detected$chr, "chr_raw")
  setnames(dt, detected$pos, "pos_hg19")
  setnames(dt, detected$ref, "other_allele")
  setnames(dt, detected$alt, "effect_allele")
  setnames(dt, detected$pval, "p_value")
  setnames(dt, detected$beta, "beta")
  setnames(dt, detected$se, "standard_error")
  if (!is.na(detected$maf)) setnames(dt, detected$maf, "EAF")
  if (!is.na(detected$rsid)) setnames(dt, detected$rsid, "rsid_raw")

  # Clean chromosome
  dt[, chr := as.integer(sub("^chr", "", chr_raw))]
  dt[, pos_hg19 := as.integer(pos_hg19)]
  dt <- dt[!is.na(chr) & chr %in% 1:22]
  dt <- dt[!is.na(beta) & !is.na(standard_error) & standard_error > 0 & !is.na(p_value)]

  cat("  After filtering:", nrow(dt), "variants\n")

  # ========================================================================
  # Liftover hg19 → hg38
  # ========================================================================
  cat("  Lifting over to GRCh38...\n")

  gr <- GRanges(
    seqnames = paste0("chr", dt$chr),
    ranges = IRanges(start = dt$pos_hg19, width = 1)
  )
  lifted <- liftOver(gr, chain)
  n_mapped <- lengths(lifted)

  dt[, pos_hg38 := NA_integer_]
  unique_idx <- which(n_mapped == 1L)
  if (length(unique_idx) > 0) {
    dt$pos_hg38[unique_idx] <- start(unlist(lifted[unique_idx]))
  }

  n_lifted <- sum(!is.na(dt$pos_hg38))
  cat("  Liftover success:", n_lifted, "/", nrow(dt),
      "(", round(100 * n_lifted / nrow(dt), 1), "%)\n")

  dt <- dt[!is.na(pos_hg38)]

  # Deduplicate by hg38 position (keep lowest p-value)
  dt <- dt[order(chr, pos_hg38, p_value)]
  dt <- dt[!duplicated(paste(chr, pos_hg38))]

  cat("  After dedup:", nrow(dt), "variants\n")
  cat("  Genome-wide significant (p<5e-8):", sum(dt$p_value < 5e-8), "\n")

  # ========================================================================
  # Write harmonised output (COLOC-ready)
  # ========================================================================
  out_dt <- dt[, .(
    chromosome = chr,
    base_pair_location = pos_hg38,
    pos_hg19 = pos_hg19,
    effect_allele,
    other_allele,
    beta,
    standard_error,
    p_value,
    EAF = if ("EAF" %in% names(dt)) EAF else NA_real_
  )]

  fwrite(out_dt, out_coloc, sep = "\t")
  cat("  Saved COLOC format:", out_coloc, "\n")

  # ========================================================================
  # Write S-PrediXcan format (CSV with rsID)
  # ========================================================================
  # S-PrediXcan expects: panel_variant_id (rsID or chr_pos_ref_alt_build),
  # effect_allele, non_effect_allele, zscore
  # Since BBJ may not have rsIDs, use chr_pos_ref_alt_b38 as variant ID
  spredixcan_dt <- dt[, .(
    panel_variant_id = paste(chr, pos_hg38, other_allele, effect_allele, "b38", sep = "_"),
    chromosome = chr,
    position = pos_hg38,
    effect_allele,
    non_effect_allele = other_allele,
    zscore = beta / standard_error,
    pvalue = p_value,
    sample_size = n
  )]

  fwrite(spredixcan_dt, out_spredixcan)
  cat("  Saved S-PrediXcan format:", out_spredixcan, "\n")

  cat("  Done:", tr, "\n\n")
}

cat("=== BBJ formatting complete ===\n")
cat("End time:", format(Sys.time()), "\n")
cat("\nNext steps:\n")
cat("  1. Run COLOC: GWAS_NAME=BBJ_ALT Rscript 47_bbj_coloc.R\n")
# 2026-04-22: TWAS-via-Script 19 step retired. MR ditched from paper; feed the
# formatted BBJ sumstats directly to the TWAS / COLOC pipeline.

