#!/usr/bin/env Rscript
# build_topld_blocks.R
# Convert TOP-LD pairwise R²/D' files into per-block PLINK-compatible LD matrices,
# matching the layout of GWAS/finemapping/data/ld_ref/{ukbb_eur,1kg_*}/chr{N}/{start}.{stop}/.
#
# Input format (per chr × ancestry):
#   {POP}_chr{N}_Dprime_extended_0.2_1000000.csv.gz
#     columns: SNP1, SNP2, R2, Dprime, +/-corr   (SNP1/SNP2 are POSITIONS)
#   {POP}_chr{N}_Dprime_extended_0.2_1000000_info.csv.gz
#     columns: Position, rsID, MAF, REF, ALT
#
# Output (per block):
#   {ld_dir}/chr{N}/{start}.{stop}/{start}.{stop}.bim   # PLINK .bim
#   {ld_dir}/chr{N}/{start}.{stop}/{start}.{stop}.ld    # PLINK --r square output (signed)
#
# Notes:
# - TOP-LD only stores pairs with R² ≥ 0.2 (sparse). Missing pairs assumed R = 0.
# - The +/-corr column gives the sign of correlation; use it with sqrt(R²) to get r.
# - Diagonal = 1.0.
# - LD blocks come from the corresponding 1kg_*/approx_LD_blocks.txt for the ancestry.
#
# Usage:
#   Rscript build_topld_blocks.R <ancestry> <chr_num>
# Example:
#   Rscript build_topld_blocks.R EUR 22

suppressPackageStartupMessages({
  library(data.table)
  library(R.utils)  # for fast .gz reading via fread
  library(rtracklayer)   # liftOver (TOP-LD is hg38; our 1kg panels are hg19)
  library(GenomicRanges)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript build_topld_blocks.R <ancestry> <chr>")
POP   <- toupper(args[1])
CHR   <- as.integer(args[2])
stopifnot(POP %in% c("EUR", "AFR", "EAS", "SAS"))
stopifnot(CHR >= 1, CHR <= 22)

PROJ <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
LD_REF <- file.path(PROJ, "GWAS/finemapping/data/ld_ref")
RAW_DIR <- file.path(LD_REF, "topld_raw", POP)
OUT_DIR <- file.path(LD_REF, paste0("topld_", tolower(POP)))
BLOCKS_FILE <- file.path(LD_REF, paste0("1kg_", tolower(POP)), "approx_LD_blocks.txt")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Copy block file to topld dir for self-containment (one-time)
target_blocks <- file.path(OUT_DIR, "approx_LD_blocks.txt")
if (!file.exists(target_blocks)) file.copy(BLOCKS_FILE, target_blocks)

cat("============================================================\n")
cat("Building TOP-LD ", POP, " chr", CHR, " — ", as.character(Sys.time()), "\n", sep = "")
cat("============================================================\n")
cat("RAW_DIR : ", RAW_DIR, "\n")
cat("OUT_DIR : ", OUT_DIR, "\n")
cat("BLOCKS  : ", BLOCKS_FILE, "\n\n")

pair_path <- file.path(RAW_DIR, sprintf("%s_chr%d_Dprime_extended_0.2_1000000.csv.gz",  POP, CHR))
info_path <- file.path(RAW_DIR, sprintf("%s_chr%d_Dprime_extended_0.2_1000000_info.csv.gz", POP, CHR))
for (f in c(pair_path, info_path, BLOCKS_FILE)) {
  if (!file.exists(f)) stop("Missing input: ", f)
}

cat("Loading info file...\n")
info <- fread(cmd = paste("zcat", shQuote(info_path)))
setnames(info, c("Position", "rsID", "MAF", "REF", "ALT"))
info[, Position := as.integer(Position)]
# Keep a copy of the hg38 position as the key into the pairs file.
info[, Position_hg38 := Position]
setkey(info, Position_hg38)
cat("  ", nrow(info), " variants loaded (raw, hg38)\n", sep = "")

# LiftOver hg38 → hg19 (our GWAS sumstats + 1kg_eur panel are hg19).
# Chain file vendored locally at data/broadaway_eqtl/ (hg38ToHg19.over.chain.gz).
cat("LiftOver hg38 → hg19...\n")
chain_path <- file.path(PROJ, "data/broadaway_eqtl/hg38ToHg19.over.chain")
if (!file.exists(chain_path)) {
  # uncompress if only .gz is present
  gz_path <- paste0(chain_path, ".gz")
  if (file.exists(gz_path)) {
    system(paste("gunzip -kc", shQuote(gz_path), ">", shQuote(chain_path)))
  } else {
    stop("LiftOver chain missing: ", chain_path)
  }
}
chain <- import.chain(chain_path)
gr38 <- GRanges(seqnames = paste0("chr", CHR),
                ranges = IRanges(start = info$Position_hg38, width = 1))
lifted <- liftOver(gr38, chain)
lifted_ok <- lengths(lifted) == 1L
info_lifted <- info[lifted_ok]
info_lifted[, Position := as.integer(start(unlist(lifted[lifted_ok])))]
cat("  ", sum(lifted_ok), " / ", length(lifted_ok), " positions lifted hg38→hg19 (",
    round(100 * sum(lifted_ok) / length(lifted_ok), 1), "%)\n", sep = "")

# Subset to positions that also exist in the corresponding 1kg_<pop> PLINK panel (hg19).
# Why: TOP-LD has 2M+ variants/chr at MAF≥0.01 — block matrices would be 10s of GB.
# Intersecting with 1kg keeps the panel size manageable (~400k/chr) AND makes the
# subsequent LD-panel comparison (Phase 3) apples-to-apples (same variant set).
kg_bim_path <- file.path(LD_REF, paste0("1kg_", tolower(POP)),
                         sprintf("chr%d_%s.bim", CHR, tolower(POP)))
if (file.exists(kg_bim_path)) {
  kg_bim <- fread(kg_bim_path, header = FALSE,
                  col.names = c("chr", "snp", "cm", "bp", "a1", "a2"))
  kept <- info_lifted$Position %in% kg_bim$bp
  cat("  intersect with 1kg_", tolower(POP), " bim (hg19): ",
      sum(kept), " / ", length(kept), " variants kept\n", sep = "")
  info <- info_lifted[kept]
  setkey(info, Position_hg38)
} else {
  cat("  WARNING: 1kg bim not found at ", kg_bim_path, " — keeping all lifted TOP-LD variants\n", sep = "")
  info <- info_lifted
  setkey(info, Position_hg38)
}
cat("  ", nrow(info), " variants after liftover + 1kg intersection\n", sep = "")

cat("Loading pairwise file (this may take a few minutes)...\n")
# Just read the columns we need (SNP1, SNP2, R2, +/-corr)
pairs <- fread(cmd = paste("zcat", shQuote(pair_path)),
               select = c("SNP1", "SNP2", "R2", "+/-corr"))
setnames(pairs, c("p1", "p2", "r2", "sign_str"))
pairs[, p1 := as.integer(p1)]
pairs[, p2 := as.integer(p2)]
pairs[, r2 := as.numeric(r2)]
pairs[, sgn := ifelse(sign_str == "+", 1, -1)]
pairs[, r := sgn * sqrt(r2)]
pairs[, c("sign_str", "sgn") := NULL]
cat("  ", nrow(pairs), " pairs loaded (raw)\n", sep = "")

# Prune pairs to those where both endpoints survived liftover + 1kg intersection.
# Pairs are keyed by hg38 positions (TOP-LD's native coord system); convert to hg19
# using the info lookup so all downstream block-based extraction uses hg19 positions.
kept_hg38 <- info$Position_hg38
pairs <- pairs[p1 %in% kept_hg38 & p2 %in% kept_hg38]
cat("  ", nrow(pairs), " pairs after liftover + 1kg-intersect prune\n", sep = "")
hg38_to_hg19 <- setNames(info$Position, as.character(info$Position_hg38))
pairs[, p1 := hg38_to_hg19[as.character(p1)]]
pairs[, p2 := hg38_to_hg19[as.character(p2)]]
pairs[, p1 := as.integer(p1)]
pairs[, p2 := as.integer(p2)]
setkey(pairs, p1, p2)
invisible(gc(verbose = FALSE))
# From here on, info$Position is hg19 and pairs p1/p2 are hg19 — the script can
# treat everything as a standard hg19 panel.

cat("Loading LD blocks for chr", CHR, "...\n", sep = "")
blocks <- fread(BLOCKS_FILE)
setnames(blocks, tolower(names(blocks)))
blocks <- blocks[chr == CHR]
cat("  ", nrow(blocks), " blocks for chr", CHR, "\n", sep = "")

n_ok <- 0; n_skip <- 0; n_empty <- 0

for (b in seq_len(nrow(blocks))) {
  bs <- blocks$start[b]
  be <- blocks$stop[b]
  block_dir <- file.path(OUT_DIR, paste0("chr", CHR), paste0(bs, ".", be))
  block_prefix <- file.path(block_dir, paste0(bs, ".", be))
  ld_path <- paste0(block_prefix, ".ld")
  bim_path <- paste0(block_prefix, ".bim")

  if (file.exists(ld_path) && file.exists(bim_path)) { n_skip <- n_skip + 1; next }

  # Variants in this block
  binfo <- info[Position >= bs & Position < be]
  if (nrow(binfo) == 0) { n_empty <- n_empty + 1; next }

  dir.create(block_dir, recursive = TRUE, showWarnings = FALSE)

  # Build position-to-index lookup
  pos_to_idx <- setNames(seq_len(nrow(binfo)), as.character(binfo$Position))
  n_var <- nrow(binfo)

  # Subset pairs to those within the block (both endpoints in binfo).
  # Use data.table join on p1 first (indexed), then filter by p2.
  bpairs <- pairs[p1 >= bs & p1 < be & p2 >= bs & p2 < be]

  # Build dense symmetric matrix; default 0, diagonal 1
  M <- matrix(0, nrow = n_var, ncol = n_var)
  diag(M) <- 1
  if (nrow(bpairs) > 0) {
    i_idx <- pos_to_idx[as.character(bpairs$p1)]
    j_idx <- pos_to_idx[as.character(bpairs$p2)]
    # Drop any pair where either endpoint isn't in the block's info list
    valid <- !is.na(i_idx) & !is.na(j_idx)
    i_idx <- i_idx[valid]; j_idx <- j_idx[valid]
    r_vals <- bpairs$r[valid]
    M[cbind(i_idx, j_idx)] <- r_vals
    M[cbind(j_idx, i_idx)] <- r_vals
  }

  # Write .ld (PLINK --r square format: space-separated, no header)
  fwrite(as.data.table(M), ld_path, sep = " ", col.names = FALSE)

  # Write .bim (PLINK format: chr, SNP_ID, cm, bp, allele1, allele2)
  bim <- data.table(
    chr   = CHR,
    snp   = ifelse(binfo$rsID == "" | is.na(binfo$rsID),
                   paste0(CHR, ":", binfo$Position, ":", binfo$REF, ":", binfo$ALT),
                   binfo$rsID),
    cm    = 0,
    bp    = binfo$Position,
    a1    = binfo$ALT,
    a2    = binfo$REF
  )
  fwrite(bim, bim_path, sep = "\t", col.names = FALSE)

  n_ok <- n_ok + 1
  if (n_ok %% 25 == 0) {
    cat("  block ", n_ok, "/", nrow(blocks), " (", bs, "-", be, ", ", n_var, " vars)\n", sep = "")
  }
  rm(M, binfo, bpairs, pos_to_idx); gc(verbose = FALSE)
}

cat("\nDone chr", CHR, ": ", n_ok, " built, ", n_skip, " skipped, ", n_empty, " empty.\n", sep = "")
cat("Total output dir size: ",
    system(paste("du -sh", shQuote(OUT_DIR), "| cut -f1"), intern = TRUE), "\n")
