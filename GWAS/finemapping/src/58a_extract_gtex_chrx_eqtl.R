#!/usr/bin/env Rscript
# 58a_extract_gtex_chrx_eqtl.R — Extract GTEx v8 Liver chrX eQTLs and
# liftover hg38 -> hg19 to match GWAS coordinate system (Broadaway-like schema).
#
# B6 — chrX × MASLD COLOC pipeline (Team B).
# Broadaway eQTL has no chrX coverage; GTEx v8 Liver (N=208) is the substitute.
# See outputs/team_C/C2_critique_sex_inference.md Issue 3/4 for the rationale.
#
# Inputs:
#   GWAS/MR_Data/GTEx_v8_Liver_eQTL.tsv.gz (hg38; chromosome col is numeric/X)
#   data/broadaway_eqtl/hg38ToHg19.over.chain
# Outputs:
#   data/gtex_v8_liver_eqtl/chrX_marginal_summary_results.tsv  (Broadaway schema)
#
# Schema match (Broadaway):
#   Entrez Variant CHR POS NEA EA EAF Beta SE PVAL N Studies GeneSymbol ENSG Gene_Biotype

suppressPackageStartupMessages({
  library(data.table)
  library(rtracklayer)
  library(GenomicRanges)
})

BASE_DIR <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

GTEX_FILE  <- file.path(BASE_DIR, "GWAS/MR_Data/GTEx_v8_Liver_eQTL.tsv.gz")
CHAIN_FILE <- file.path(BASE_DIR, "data/broadaway_eqtl/hg38ToHg19.over.chain")
OUT_DIR    <- file.path(BASE_DIR, "data/gtex_v8_liver_eqtl")
OUT_FILE   <- file.path(OUT_DIR, "chrX_marginal_summary_results.tsv")
SAMPLE_N   <- 208L  # GTEx v8 Liver eQTL discovery N

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("============================================================\n")
cat("58a: extract GTEx v8 Liver chrX eQTLs + liftover hg38 -> hg19\n")
cat("============================================================\n")
cat("Input:  ", GTEX_FILE, "\n")
cat("Chain:  ", CHAIN_FILE, "\n")
cat("Output: ", OUT_FILE, "\n\n")

# ---- 1. Read GTEx Liver eQTL, filter to chrX -------------------------------
cat("[1/4] Reading GTEx v8 Liver eQTL (may be slow; ~3.3GB compressed)...\n")
t0 <- Sys.time()
# The 'chromosome' column is numeric (1-22) or 'X' — read all then filter
gtex <- fread(cmd = paste0("zcat ", shQuote(GTEX_FILE)),
              showProgress = FALSE, data.table = TRUE)
cat("  Loaded", nrow(gtex), "rows in", round(as.numeric(Sys.time() - t0, units = "secs"), 1), "s\n")
cat("  Columns:", paste(names(gtex), collapse = ", "), "\n")

# Filter to chrX
chrx <- gtex[chromosome == "X" | chromosome == "chrX" | chromosome == "23"]
rm(gtex); invisible(gc())
cat("  chrX rows:", nrow(chrx), "\n")

if (nrow(chrx) == 0) {
  stop("ERROR: No chrX rows found in GTEx v8 Liver eQTL. Check column convention.")
}

# Sanity: report unique chromosomes encountered
cat("  Unique chromosome values in filtered table:",
    paste(unique(chrx$chromosome), collapse = ", "), "\n")

# ---- 2. Liftover hg38 -> hg19 ----------------------------------------------
cat("\n[2/4] Liftover hg38 -> hg19...\n")
chain <- import.chain(CHAIN_FILE)

# Variant id format: chrX_pos_ref_alt_b38 (hg38). Use the explicit position col.
chrx[, hg38_pos := as.integer(position)]
gr_hg38 <- GRanges(
  seqnames = "chrX",
  ranges   = IRanges(start = chrx$hg38_pos, end = chrx$hg38_pos),
  hg38_pos = chrx$hg38_pos
)

cat("  hg38 records:", length(gr_hg38), "\n")
gr_hg19_list <- liftOver(gr_hg38, chain)
n_mapped <- sum(lengths(gr_hg19_list) == 1L)
cat("  Mapped uniquely to hg19:", n_mapped, "/", length(gr_hg38),
    sprintf("(%.1f%%)", 100 * n_mapped / length(gr_hg38)), "\n")

# Keep only uniquely-mapped variants on chrX in hg19
ok_idx <- which(lengths(gr_hg19_list) == 1L)
gr_ok  <- unlist(gr_hg19_list[ok_idx])
keep_seq <- as.character(seqnames(gr_ok)) == "chrX"
ok_idx   <- ok_idx[keep_seq]
gr_ok    <- gr_ok[keep_seq]

cat("  Retained chrX-mapped variants:", length(gr_ok), "\n")

chrx_lifted <- chrx[ok_idx]
chrx_lifted[, POS := start(gr_ok)]
chrx_lifted[, CHR := "X"]

rm(chain, gr_hg38, gr_hg19_list, gr_ok); invisible(gc())

# ---- 3. Reformat to Broadaway schema --------------------------------------
cat("\n[3/4] Reformatting to Broadaway-compatible schema...\n")

# GTEx columns: variant r2 pvalue molecular_trait_object_id molecular_trait_id
#               maf gene_id median_tpm beta se an ac chromosome position ref alt
#               type rsid
# Broadaway target: Entrez Variant CHR POS NEA EA EAF Beta SE PVAL N Studies
#                   GeneSymbol ENSG Gene_Biotype
# ENSG conversion: strip version suffix if present (e.g., ENSG0000xxxx.1 -> ENSG0000xxxx)
chrx_lifted[, ENSG_BASE := sub("\\..*$", "", gene_id)]

out <- chrx_lifted[, .(
  Entrez       = NA_character_,
  Variant      = paste0(CHR, "_", POS, "_", ref, "_", alt),
  CHR          = CHR,
  POS          = POS,
  NEA          = ref,    # GTEx: ref allele = non-effect
  EA           = alt,    # GTEx: alt allele = effect allele (beta is per alt)
  EAF          = maf,    # MAF is reported; not strictly EAF but downstream COLOC
                         # uses only beta/se (sdY=1), EAF is informational
  Beta         = beta,
  SE           = se,
  PVAL         = pvalue,
  N            = SAMPLE_N,
  Studies      = "GTEx_v8",
  GeneSymbol   = ENSG_BASE,  # fall back to ENSG until we add a symbol map
  ENSG         = ENSG_BASE,
  Gene_Biotype = "unknown"
)]

# Sanity filters
out <- out[!is.na(Beta) & !is.na(SE) & SE > 0]
cat("  Final variant count:", nrow(out), "\n")
cat("  Unique eGenes:", length(unique(out$ENSG)), "\n")

# Per-gene SNP count summary
gene_counts <- out[, .N, by = ENSG]
cat("  Genes with >=100 SNPs:", sum(gene_counts$N >= 100), "\n")
cat("  Genes with >=50 SNPs:",  sum(gene_counts$N >= 50),  "\n")

# ---- 4. Write output -------------------------------------------------------
cat("\n[4/4] Writing", OUT_FILE, "\n")
fwrite(out, OUT_FILE, sep = "\t")
cat("  Done. Size:", round(file.info(OUT_FILE)$size / 1e6, 1), "MB\n")

cat("\n============================================================\n")
cat("58a: complete\n")
cat("============================================================\n")
