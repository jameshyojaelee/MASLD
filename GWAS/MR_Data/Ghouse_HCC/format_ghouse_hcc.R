#!/usr/bin/env Rscript
# format_ghouse_hcc.R
# ---------------------------------------------------------------------------
# Format Ghouse 2025 HCC GWAS for COLOC (GRCh37 → GRCh38 liftover)
#
# Ghouse et al. 2025, Nat Genet — "Genome-wide meta-analysis identifies nine
# loci associated with higher risk of hepatocellular carcinoma development"
# GCST90809296: European-only (3,748 cases / 1,861,536 controls)
#
# Input:  GCST90809296_EUR_HCC_raw.tsv (GRCh37)
# Output: GCST90809296_EUR_HCC_harmonised_hg38.tsv.gz
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(rtracklayer)
  library(GenomicRanges)
})

BASE_DIR   <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
HCC_DIR    <- file.path(BASE_DIR, "GWAS/MR_Data/Ghouse_HCC")
CHAIN_FILE <- file.path(BASE_DIR, "data/broadaway_eqtl/hg19ToHg38.over.chain")
CHAIN_GZ   <- paste0(CHAIN_FILE, ".gz")

infile  <- file.path(HCC_DIR, "GCST90809296_EUR_HCC_raw.tsv")
outfile <- file.path(HCC_DIR, "GCST90809296_EUR_HCC_harmonised_hg38.tsv.gz")

cat("=== Ghouse HCC GWAS Formatting & Liftover ===\n")
cat("Start time:", format(Sys.time()), "\n")

# Chain file
if (!file.exists(CHAIN_FILE)) {
  if (!file.exists(CHAIN_GZ)) {
    download.file("https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz",
                  CHAIN_GZ, mode = "wb", quiet = FALSE)
  }
  system2("gunzip", args = c("-k", CHAIN_GZ))
}
chain <- import.chain(CHAIN_FILE)

# Load
dt <- fread(infile)
cat("  Raw rows:", format(nrow(dt), big.mark = ","), "\n")
cat("  Columns:", paste(names(dt), collapse = ", "), "\n")

# Auto-detect column mapping (GWAS Catalog SSF or custom)
# Expected columns might be: chromosome/chr, base_pair_location/pos,
# effect_allele, other_allele, beta, standard_error, p_value
col_candidates <- list(
  chr  = c("chromosome", "chr", "#chrom", "CHROM"),
  pos  = c("base_pair_location", "pos", "position", "BP", "POS"),
  ea   = c("effect_allele", "alt", "ALT", "A1", "EA"),
  oa   = c("other_allele", "ref", "REF", "A2", "NEA"),
  beta = c("beta", "BETA", "effect"),
  se   = c("standard_error", "se", "SE", "sebeta"),
  pval = c("p_value", "pval", "P", "p", "PVAL")
)

detected <- list()
for (field in names(col_candidates)) {
  match <- intersect(col_candidates[[field]], names(dt))
  detected[[field]] <- if (length(match) > 0) match[1] else NA
  cat("  ", field, "->", ifelse(is.na(detected[[field]]), "NOT FOUND", detected[[field]]), "\n")
}

# Validate
missing <- names(detected)[sapply(detected, is.na)]
if (length(missing) > 0) stop("Missing columns: ", paste(missing, collapse = ", "))

# Standardize
setnames(dt, detected$chr, "chr_raw")
setnames(dt, detected$pos, "pos_hg19")
setnames(dt, detected$ea, "effect_allele")
setnames(dt, detected$oa, "other_allele")
setnames(dt, detected$beta, "beta")
setnames(dt, detected$se, "standard_error")
setnames(dt, detected$pval, "p_value")

dt[, chr := as.integer(sub("^chr", "", chr_raw))]
dt[, pos_hg19 := as.integer(pos_hg19)]
dt <- dt[!is.na(chr) & chr %in% 1:22]
dt <- dt[!is.na(beta) & !is.na(standard_error) & standard_error > 0 & !is.na(p_value)]
cat("  After QC:", format(nrow(dt), big.mark = ","), "\n")

# Liftover
cat("  Lifting over to GRCh38...\n")
gr <- GRanges(seqnames = paste0("chr", dt$chr), ranges = IRanges(start = dt$pos_hg19, width = 1))
lifted <- liftOver(gr, chain)
n_mapped <- lengths(lifted)
dt[, pos_hg38 := NA_integer_]
idx <- which(n_mapped == 1L)
if (length(idx) > 0) dt$pos_hg38[idx] <- start(unlist(lifted[idx]))
n_lifted <- sum(!is.na(dt$pos_hg38))
cat("  Liftover:", n_lifted, "/", nrow(dt), "(", round(100*n_lifted/nrow(dt), 1), "%)\n")
dt <- dt[!is.na(pos_hg38)]

# Dedup
dt <- dt[order(chr, pos_hg38, p_value)]
dt <- dt[!duplicated(paste(chr, pos_hg38))]
cat("  Final:", format(nrow(dt), big.mark = ","), "\n")
cat("  GW-sig:", sum(dt$p_value < 5e-8), "\n")

# Write
eaf_col <- intersect(c("effect_allele_frequency", "EAF", "eaf"), names(dt))
out_dt <- dt[, .(chromosome = chr, base_pair_location = pos_hg38, pos_hg19,
                  effect_allele, other_allele, beta, standard_error, p_value,
                  EAF = if (length(eaf_col) > 0) get(eaf_col[1]) else NA_real_)]
fwrite(out_dt, outfile, sep = "\t")
cat("  Saved:", basename(outfile), "\n")
cat("Done. End time:", format(Sys.time()), "\n")
