#!/usr/bin/env Rscript
# 35s_susie_coloc_broadaway.R
# ---------------------------------------------------------------------------
# SuSiE-COLOC: Broadaway Liver eQTL (N=1,183) x EUR GWAS with LD from 1000G
#
# Extends Script 35b/35g (ABF-only COLOC) by adding per-gene LD extraction
# from 1000 Genomes EUR PLINK2 files and running coloc.susie() alongside
# coloc.abf() for direct comparison. SuSiE fine-mapping resolves multiple
# causal signals per locus, improving colocalization at complex loci.
#
# Supports 6 EUR GWAS via GWAS_NAME env var:
#   UKBB_ALT, UKBB_AST, UKBB_GGT (quantitative, N=343,850)
#   FINNGEN_NAFLD, FINNGEN_NASH, FINNGEN_HCC (case-control)
#
# Key additions over 35b/35g:
#   - Per-gene LD extraction via PLINK2 from 1000G EUR (hg19, 379 samples)
#   - Triple-intersect: SNPs in eQTL AND GWAS AND LD matrix
#   - susie_rss() fine-mapping on both GWAS and eQTL z-scores
#   - coloc.susie() for multi-signal colocalization
#   - ABF fallback when SuSiE fails or < 50 triple-intersect SNPs
#
# Requires: PLINK2 in PATH (module load PLINK/2.0a5.13 in sbatch wrapper)
#           PLINK 1.9 at /nfs/sw/easybuild/software/PLINK/1.9b_6.21-x86_64/plink
#           (used for --r square text output; PLINK2 a5.13 outputs binary .vcor1)
#
# Usage:
#   GWAS_NAME=UKBB_ALT Rscript 35s_susie_coloc_broadaway.R
#   GWAS_NAME=FINNGEN_NAFLD Rscript 35s_susie_coloc_broadaway.R
#
# Inputs:
#   - data/broadaway_eqtl/Liver_eQTL_Meta_Leads_ST3_20240530.tsv
#   - data/broadaway_eqtl/chr{1-22}_marginal_summary_results.tsv
#   - GWAS/MR_Data/{GWAS_FILE}
#   - data/broadaway_eqtl/hg19ToHg38.over.chain
#   - data/1kg_eur/chr{1-22}_eur.{pgen,psam,pvar}
#
# Outputs:
#   - RNA-seq/results/causal_inference/susie_broadaway_{gwas_name}/
#     - susie_coloc_results.csv
#     - susie_summary.csv
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(coloc)
  library(susieR)
  library(data.table)
  library(rtracklayer)
  library(GenomicRanges)
  library(ggplot2)
})

# ==============================================================================
# Configuration — parameterized by GWAS_NAME env var
# ==============================================================================
GWAS_NAME  <- Sys.getenv("GWAS_NAME", "UKBB_ALT")
CHR_FILTER <- as.integer(Sys.getenv("CHR_FILTER", "0"))
if (is.na(CHR_FILTER) || CHR_FILTER == 0L) CHR_FILTER <- NULL

gwas_config <- list(
  UKBB_ALT = list(
    file      = "GCST90019492_UKBB_ALT_harmonised.tsv.gz",
    N         = 343850L,
    type      = "quant",
    label     = "UKBB ALT (Sinnott-Armstrong 2021)",
    is_finngen = FALSE
  ),
  UKBB_AST = list(
    file      = "GCST90019497_UKBB_AST_harmonised.tsv.gz",
    N         = 343850L,
    type      = "quant",
    label     = "UKBB AST (Sinnott-Armstrong 2021)",
    is_finngen = FALSE
  ),
  UKBB_GGT = list(
    file      = "GCST90019507_UKBB_GGT_harmonised.tsv.gz",
    N         = 343850L,
    type      = "quant",
    label     = "UKBB GGT (Sinnott-Armstrong 2021)",
    is_finngen = FALSE
  ),
  FINNGEN_NAFLD = list(
    file      = "FinnGen/finngen_R12_NAFLD.gz",
    N_cases   = 4614L,
    N_controls = 434243L,
    type      = "cc",
    label     = "FinnGen R12 NAFLD (K11_NAFLD)",
    is_finngen = TRUE
  ),
  FINNGEN_NASH = list(
    file      = "FinnGen/finngen_R12_CHIRHEP_NAS.gz",
    N_cases   = 1823L,
    N_controls = 434243L,
    type      = "cc",
    label     = "FinnGen R12 NASH/Cirrhosis-Hepatitis NAS",
    is_finngen = TRUE
  ),
  FINNGEN_HCC = list(
    file      = "FinnGen/finngen_R12_C3_HEPATOCELLU_CARC_EXALLC.gz",
    N_cases   = 1156L,
    N_controls = 291387L,
    type      = "cc",
    label     = "FinnGen R12 HCC",
    is_finngen = TRUE
  )
)

if (!GWAS_NAME %in% names(gwas_config)) {
  stop("Unknown GWAS_NAME: ", GWAS_NAME,
       ". Valid options: ", paste(names(gwas_config), collapse = ", "))
}

cfg <- gwas_config[[GWAS_NAME]]

BASE_DIR    <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
EQTL_DIR    <- file.path(BASE_DIR, "data/broadaway_eqtl")
GWAS_FILE   <- file.path(BASE_DIR, "GWAS/MR_Data", cfg$file)
LD_DIR      <- file.path(BASE_DIR, "data/1kg_eur")

gwas_short  <- tolower(GWAS_NAME)
RESULTS_DIR <- file.path(BASE_DIR, "RNA-seq/results/causal_inference",
                          paste0("susie_broadaway_", gwas_short))

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

# Per-chromosome output suffix when running in array mode
CHR_SUFFIX <- if (!is.null(CHR_FILTER)) paste0("_chr", CHR_FILTER) else ""

# Compute effective N for case-control or quantitative
if (cfg$type == "cc") {
  GWAS_N_CASES    <- cfg$N_cases
  GWAS_N_CONTROLS <- cfg$N_controls
  GWAS_N          <- GWAS_N_CASES + GWAS_N_CONTROLS
  GWAS_S          <- GWAS_N_CASES / GWAS_N
} else {
  GWAS_N    <- cfg$N
  GWAS_S    <- NULL
}
GWAS_TYPE <- cfg$type

# Broadaway eQTL parameters
EQTL_N    <- 1183L

# COLOC priors
COLOC_P1  <- 1e-4
COLOC_P2  <- 1e-4
COLOC_P12 <- 5e-6

# SuSiE parameters
SUSIE_L          <- 10L    # max number of causal signals
MIN_TRIPLE_SNPS  <- 50L    # minimum triple-intersect SNPs for SuSiE
LD_REGULARIZE    <- 1e-3   # ridge regularization for LD matrix

# Liftover chain file
CHAIN_FILE <- file.path(EQTL_DIR, "hg19ToHg38.over.chain")
CHAIN_GZ   <- paste0(CHAIN_FILE, ".gz")

cat("=== Script 35s: SuSiE-COLOC — Broadaway eQTL x", GWAS_NAME, "===\n")
if (!is.null(CHR_FILTER)) cat("CHR_FILTER:", CHR_FILTER, "(processing single chromosome)\n")
cat("GWAS:", cfg$label, "\n")
cat("GWAS file:", cfg$file, "\n")
cat("GWAS N:", format(GWAS_N, big.mark = ","), "\n")
if (cfg$type == "cc") {
  cat("GWAS type: case-control (s =", round(GWAS_S, 4), ")\n")
} else {
  cat("GWAS type: quantitative\n")
}
cat("LD reference:", LD_DIR, "\n")
cat("SuSiE L:", SUSIE_L, " | Min triple SNPs:", MIN_TRIPLE_SNPS, "\n")
cat("Results dir:", RESULTS_DIR, "\n")
cat("Start time:", format(Sys.time()), "\n\n")

# ==============================================================================
# 1. Download chain file if needed
# ==============================================================================
cat("--- Step 1: Preparing liftover chain ---\n")

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
cat("  Chain loaded:", length(chain), "chains\n")

# ==============================================================================
# 2. Load GWAS
# ==============================================================================
cat("\n--- Step 2: Loading", GWAS_NAME, "GWAS ---\n")

if (!file.exists(GWAS_FILE)) {
  stop("GWAS file not found: ", GWAS_FILE)
}

gwas_raw <- fread(GWAS_FILE)
cat("  Raw GWAS rows:", nrow(gwas_raw), "\n")

if (cfg$is_finngen) {
  # FinnGen format: #chrom/chrom, pos, ref, alt, beta, sebeta, pval
  chrom_col <- grep("chrom", names(gwas_raw), value = TRUE)[1]
  setnames(gwas_raw, chrom_col, "chr_raw", skip_absent = TRUE)

  gwas_raw[, chr := as.integer(sub("^chr", "", chr_raw))]
  gwas_raw[, pos_hg38 := as.integer(pos)]
  gwas_raw <- gwas_raw[!is.na(chr) & chr %in% 1:22]
  gwas_raw <- gwas_raw[!is.na(beta) & !is.na(sebeta) & sebeta > 0 & !is.na(pval)]

  # Standardize allele column names
  setnames(gwas_raw, c("alt", "ref"), c("effect_allele", "other_allele"))

  # Standardize SE column name
  gwas_raw[, standard_error := sebeta]

  # Deduplicate
  gwas_raw <- gwas_raw[order(chr, pos_hg38, pval)]
  gwas_raw <- gwas_raw[!duplicated(paste(chr, pos_hg38))]
  gwas_raw[, merge_key := paste0(chr, ":", pos_hg38)]
  gwas_raw[, p_value := pval]

} else {
  # UKBB format: chromosome, base_pair_location, effect_allele, other_allele, beta, standard_error, p_value
  gwas_raw <- gwas_raw[!is.na(chromosome) & !is.na(base_pair_location) & !is.na(p_value)]
  gwas_raw[, chr := as.integer(chromosome)]
  gwas_raw[, pos_hg38 := as.integer(base_pair_location)]
  gwas_raw <- gwas_raw[!is.na(chr) & chr %in% 1:22]
  gwas_raw <- gwas_raw[!is.na(beta) & !is.na(standard_error) & standard_error > 0]

  # Deduplicate
  gwas_raw <- gwas_raw[order(chr, pos_hg38, p_value)]
  gwas_raw <- gwas_raw[!duplicated(paste(chr, pos_hg38))]
  gwas_raw[, merge_key := paste0(chr, ":", pos_hg38)]
}

gwas <- gwas_raw
rm(gwas_raw)

cat("  GWAS after dedup:", nrow(gwas), "variants\n")
cat("  Genome-wide significant (p<5e-8):", sum(gwas$p_value < 5e-8), "\n")

# ==============================================================================
# 3. Load Broadaway leads — ALL eGenes (no DEG filter)
# ==============================================================================
cat("\n--- Step 3: Loading Broadaway eGene leads ---\n")

leads <- fread(file.path(EQTL_DIR, "Liver_eQTL_Meta_Leads_ST3_20240530.tsv"))
cat("  Total eQTL signals:", nrow(leads), "\n")

all_eGenes <- unique(leads$Gene)
cat("  Unique eGenes (ALL, no DEG filter):", length(all_eGenes), "\n")

# Build gene -> chr mapping from leads Variant field (format: CHR_POS_REF_ALT)
leads[, eqtl_chr := as.integer(sub("_.*", "", Variant))]
gene_chr_map <- unique(leads[!is.na(eqtl_chr), .(Gene, chr = eqtl_chr)])
gene_chr_map <- gene_chr_map[, .N, by = .(Gene, chr)][order(-N)][!duplicated(Gene)]
cat("  Gene-to-chromosome mapping:", nrow(gene_chr_map), "genes\n")

# Filter to single chromosome when running in array mode
if (!is.null(CHR_FILTER)) {
  chr_genes <- gene_chr_map[chr == CHR_FILTER]$Gene
  all_eGenes <- intersect(all_eGenes, chr_genes)
  cat("  CHR_FILTER=", CHR_FILTER, ": ", length(all_eGenes), " eGenes on this chromosome\n")
  if (length(all_eGenes) == 0) {
    cat("  No eGenes on chr", CHR_FILTER, ". Exiting.\n")
    quit(save = "no", status = 0)
  }
}

# Build gene symbol <-> Ensembl lookup from leads
leads_lookup <- unique(leads[, .(Gene, Ensembl)])

# Verify LD reference exists
ld_test <- file.path(LD_DIR, "chr1_eur.pgen")
if (!file.exists(ld_test)) {
  stop("1000G EUR LD reference not found: ", ld_test,
       "\nExpected PLINK2 pgen/psam/pvar files at ", LD_DIR)
}
cat("  LD reference verified:", LD_DIR, "\n")

# Verify PLINK available (PLINK/2.0a5.13 module installs as 'plink' not 'plink2')
plink_check <- system("which plink", intern = TRUE, ignore.stderr = TRUE)
if (length(plink_check) == 0 || !nzchar(plink_check)) {
  stop("plink not found in PATH. Load it via: module load PLINK/2.0a5.13")
}
cat("  PLINK found:", plink_check, "\n")

# ==============================================================================
# 4. Helper functions
# ==============================================================================

# Liftover hg19 -> hg38 (vectorized)
liftover_positions <- function(chr_num, positions) {
  gr <- GRanges(
    seqnames = paste0("chr", chr_num),
    ranges = IRanges(start = positions, width = 1)
  )
  lifted <- liftOver(gr, chain)
  n_mapped <- lengths(lifted)
  hg38_pos <- rep(NA_integer_, length(positions))
  unique_idx <- which(n_mapped == 1L)
  if (length(unique_idx) > 0) {
    hg38_pos[unique_idx] <- start(unlist(lifted[unique_idx]))
  }
  hg38_pos
}

# Per-chromosome eQTL cache (stores BOTH hg19 POS and hg38 pos_hg38)
chr_eqtl_cache <- list()

load_chr_eqtl <- function(chr_num) {
  key <- as.character(chr_num)
  if (!is.null(chr_eqtl_cache[[key]])) return(chr_eqtl_cache[[key]])

  fname <- file.path(EQTL_DIR, paste0("chr", chr_num, "_marginal_summary_results.tsv"))
  if (!file.exists(fname)) return(NULL)

  cat("    Loading chr", chr_num, "eQTL data...\n")
  dt <- fread(fname)

  # Liftover all positions for this chromosome at once
  cat("    Lifting over", nrow(dt), "positions from hg19 to hg38...\n")
  dt[, pos_hg38 := liftover_positions(chr_num, POS)]
  dt <- dt[!is.na(pos_hg38)]
  cat("    Successful liftover:", nrow(dt), "variants\n")

  # Create merge key (hg38 for matching GWAS)
  dt[, merge_key := paste0(CHR, ":", pos_hg38)]

  chr_eqtl_cache[[key]] <<- dt
  return(dt)
}

# Allele harmonization (Broadaway EA/NEA vs GWAS effect_allele/other_allele)
harmonize_alleles <- function(dt) {
  dt <- copy(dt)
  comp <- c("A" = "T", "T" = "A", "C" = "G", "G" = "C")

  dt[, match_type := "none"]

  # Direct: EA == effect_allele AND NEA == other_allele
  dt[EA == effect_allele & NEA == other_allele, match_type := "direct"]
  # Flipped: EA == other_allele AND NEA == effect_allele
  dt[EA == other_allele & NEA == effect_allele, match_type := "flipped"]

  # Remove strand-ambiguous variants that failed direct/flipped match
  dt[, is_ambiguous := (EA %in% c("A","T") & NEA %in% c("A","T")) |
                        (EA %in% c("C","G") & NEA %in% c("C","G"))]
  n_ambig_unmatched <- sum(dt$is_ambiguous & dt$match_type == "none", na.rm = TRUE)
  cat("    Strand-ambiguous unmatched:", n_ambig_unmatched, "\n")
  dt <- dt[!(is_ambiguous & match_type == "none")]

  # Complement matching for remaining non-ambiguous strand SNPs
  dt[match_type == "none" & nchar(EA) == 1 & nchar(NEA) == 1,
     match_type := fifelse(
       comp[EA] == effect_allele & comp[NEA] == other_allele, "direct",
       fifelse(comp[EA] == other_allele & comp[NEA] == effect_allele, "flipped",
               "none")
     )]

  n_direct  <- sum(dt$match_type == "direct")
  n_flipped <- sum(dt$match_type == "flipped")
  n_none    <- sum(dt$match_type == "none")
  cat("    Allele harmonization: direct=", n_direct, " flipped=", n_flipped,
      " dropped=", n_none, "\n")

  dt <- dt[match_type != "none"]

  # Flip eQTL beta when alleles are swapped
  if (any(dt$match_type == "flipped")) {
    dt[match_type == "flipped", Beta := -Beta]
  }

  dt[, is_ambiguous := NULL]
  return(dt)
}

# Extract LD matrix for a genomic region using PLINK2
# Returns list(R = correlation matrix, vars = data.table with chr, pos columns)
# or NULL on failure
extract_ld_matrix <- function(chr_num, min_pos, max_pos) {
  pfile_prefix <- file.path(LD_DIR, paste0("chr", chr_num, "_eur"))

  # Check pgen file exists
  if (!file.exists(paste0(pfile_prefix, ".pgen"))) {
    cat("    LD: pgen not found for chr", chr_num, "\n")
    return(NULL)
  }

  # Create temp prefix for PLINK2 output
  tmp_prefix <- file.path(tempdir(), paste0("ld_chr", chr_num, "_",
                                             min_pos, "_", max_pos))

  # Two-step LD computation:
  #   Step 1: PLINK2 extracts region from pgen → bed/bim/fam
  #   Step 2: PLINK 1.9 computes --r square (text output)
  # Reason: PLINK2 a5.13 --r square outputs binary .vcor1 which fread can't parse.
  # PLINK 1.9 --r square outputs tab-delimited text .ld files.

  # Step 1: Extract region to bed format using PLINK2
  bed_prefix <- paste0(tmp_prefix, "_bed")
  plink2_bin <- Sys.getenv("PLINK2_BIN", unset = "/nfs/sw/easybuild/software/PLINK/2.0a5.13/bin/plink")  # review hardcoded-path: env-overridable, default unchanged
  cmd_extract <- paste0(plink2_bin, " --pfile ", pfile_prefix,
                        " --chr ", chr_num,
                        " --from-bp ", min_pos,
                        " --to-bp ", max_pos,
                        " --make-bed",
                        " --out ", bed_prefix,
                        " --threads 1 --memory 2000")

  ret1 <- system(cmd_extract, intern = FALSE, ignore.stdout = TRUE, ignore.stderr = TRUE)

  if (ret1 != 0 || !file.exists(paste0(bed_prefix, ".bed"))) {
    unlink(list.files(tempdir(), pattern = basename(tmp_prefix), full.names = TRUE))
    return(NULL)
  }

  # Step 2: Compute LD matrix using PLINK 1.9 (produces text .ld file)
  plink19 <- Sys.getenv("PLINK19_BIN", unset = "/nfs/sw/easybuild/software/PLINK/1.9b_6.21-x86_64/plink")  # review hardcoded-path: env-overridable, default unchanged
  cmd_ld <- paste0(plink19, " --bfile ", bed_prefix,
                   " --r square",
                   " --out ", tmp_prefix,
                   " --threads 1 --memory 2000")

  ret2 <- system(cmd_ld, intern = FALSE, ignore.stdout = TRUE, ignore.stderr = TRUE)

  # Check output files exist (PLINK 1.9 outputs .ld and uses .bim for variant IDs)
  ld_file <- paste0(tmp_prefix, ".ld")
  bim_file <- paste0(bed_prefix, ".bim")

  if (ret2 != 0 || !file.exists(ld_file) || !file.exists(bim_file)) {
    # Clean up partial files
    unlink(list.files(tempdir(), pattern = basename(tmp_prefix), full.names = TRUE))
    return(NULL)
  }

  # Read LD matrix (PLINK 1.9 --r square: tab-delimited text, no header)
  ld_mat <- tryCatch(
    as.matrix(fread(ld_file, header = FALSE)),
    error = function(e) NULL
  )

  # Read variant info from .bim file (chr, var_id, cM, pos, A1, A2)
  bim_dt <- tryCatch(
    fread(bim_file, header = FALSE, col.names = c("chr", "var_id", "cM", "pos", "A1", "A2")),
    error = function(e) NULL
  )

  # Clean up temp files
  unlink(list.files(tempdir(), pattern = basename(tmp_prefix), full.names = TRUE))

  if (is.null(ld_mat) || is.null(bim_dt) || nrow(bim_dt) == 0) {
    return(NULL)
  }

  ld_vars_raw <- bim_dt$var_id

  # Verify dimensions match
  if (nrow(ld_mat) != nrow(bim_dt) || ncol(ld_mat) != nrow(bim_dt)) {
    cat("    LD: dimension mismatch (matrix:", nrow(ld_mat), "x", ncol(ld_mat),
        " vars:", nrow(bim_dt), ")\n")
    return(NULL)
  }

  # Parse variant IDs to extract positions
  # .bim has positions directly in column 4, but var_id format is chr:pos:ref:alt
  # Carry .bim A1/A2 (review A04#1): PLINK --r square signs correlations to A1, so the
  # LD must be re-signed to the z-score (effect-allele) orientation before susie_rss.
  ld_vars_dt <- data.table(var_id = ld_vars_raw, ld_pos = bim_dt$pos,
                           ld_a1 = toupper(bim_dt$A1), ld_a2 = toupper(bim_dt$A2))

  # Remove variants with failed position parsing
  valid_idx <- which(!is.na(ld_vars_dt$ld_pos))
  if (length(valid_idx) < 10) {
    return(NULL)
  }

  # Subset to valid variants
  ld_mat <- ld_mat[valid_idx, valid_idx, drop = FALSE]
  ld_vars_dt <- ld_vars_dt[valid_idx]

  # Deduplicate by position (keep first occurrence)
  dup_idx <- duplicated(ld_vars_dt$ld_pos)
  if (any(dup_idx)) {
    keep_idx <- which(!dup_idx)
    ld_mat <- ld_mat[keep_idx, keep_idx, drop = FALSE]
    ld_vars_dt <- ld_vars_dt[keep_idx]
  }

  return(list(R = ld_mat, vars = ld_vars_dt))
}

# ==============================================================================
# 5. Run SuSiE-COLOC + ABF per eGene
# ==============================================================================
cat("\n--- Step 5: Running SuSiE-COLOC + ABF ---\n")

candidate_chrs <- unique(gene_chr_map[Gene %in% all_eGenes]$chr)
cat("  Chromosomes to process:", paste(sort(candidate_chrs), collapse = ", "), "\n")
cat("  Total eGenes to test:", length(all_eGenes), "\n\n")

results_list <- list()
n_tested       <- 0L
n_skipped      <- 0L
n_low_snps     <- 0L
n_no_harm      <- 0L
n_susie_ok     <- 0L
n_abf_fallback <- 0L
n_ld_fail      <- 0L

# Checkpointing: resume from previous partial run if checkpoint exists
checkpoint_file <- file.path(RESULTS_DIR, paste0("susie_coloc_checkpoint", CHR_SUFFIX, ".csv"))
processed_genes <- character(0)
if (file.exists(checkpoint_file)) {
  prev_results <- fread(checkpoint_file)
  if (nrow(prev_results) > 0 && "gene" %in% names(prev_results)) {
    processed_genes <- unique(prev_results$gene)
    results_list <- split(prev_results, seq_len(nrow(prev_results)))
    n_tested <- length(processed_genes)
    cat("  Resuming from checkpoint:", n_tested, "genes already processed\n")
  }
}
checkpoint_interval <- 50L
last_checkpoint <- n_tested

for (gene_name in all_eGenes) {
  # Skip if already processed in checkpoint
  if (gene_name %in% processed_genes) next
  # Get chromosome for this gene
  chr_info <- gene_chr_map[Gene == gene_name]
  if (nrow(chr_info) == 0) {
    n_skipped <- n_skipped + 1L
    next
  }
  gene_chr <- chr_info$chr[1]

  # Load chromosome eQTL data (cached)
  eqtl_data <- load_chr_eqtl(gene_chr)
  if (is.null(eqtl_data)) {
    n_skipped <- n_skipped + 1L
    next
  }

  # Get eQTL data for this gene
  gene_eqtl <- eqtl_data[GeneSymbol == gene_name]
  if (nrow(gene_eqtl) == 0) {
    # Try matching on Ensembl ID
    gene_ensg <- leads_lookup$Ensembl[leads_lookup$Gene == gene_name]
    if (length(gene_ensg) > 0) {
      gene_eqtl <- eqtl_data[ENSG %in% gene_ensg]
    }
  }

  if (nrow(gene_eqtl) < 10) {
    n_low_snps <- n_low_snps + 1L
    next
  }

  # ---- (a) Get hg19 position range from eQTL native POS column ----
  eqtl_pos_hg19 <- gene_eqtl$POS
  pos_min_hg19  <- min(eqtl_pos_hg19)
  pos_max_hg19  <- max(eqtl_pos_hg19)

  # ---- (e) Liftover already done in load_chr_eqtl; create hg38 merge key ----
  # gene_eqtl already has pos_hg38 and merge_key from cache

  # ---- (f) Merge eQTL (hg38) with GWAS (hg38) on chr:pos ----
  merged <- merge(gene_eqtl, gwas, by = "merge_key", suffixes = c(".eqtl", ".gwas"))

  if (nrow(merged) < 10) {
    n_low_snps <- n_low_snps + 1L
    next
  }

  # Filter to rows with valid beta/SE in both datasets
  # Note: FinnGen sebeta already copied to standard_error in GWAS loading step
  merged <- merged[!is.na(Beta) & !is.na(SE) & SE > 0 &
                   !is.na(beta) & !is.na(standard_error) & standard_error > 0]

  if (nrow(merged) < 10) {
    n_low_snps <- n_low_snps + 1L
    next
  }

  # ---- (g) Allele harmonization ----
  merged <- harmonize_alleles(merged)

  if (nrow(merged) < 10) {
    n_no_harm <- n_no_harm + 1L
    next
  }

  n_snps_merged <- nrow(merged)

  # ---- (c) Extract LD from 1000G EUR for hg19 range ----
  ld_result <- extract_ld_matrix(gene_chr, pos_min_hg19, pos_max_hg19)

  # ---- (h) Triple-intersect: eQTL + GWAS + LD matrix ----
  # Match LD matrix SNP positions (hg19) to eQTL native POS (hg19)
  use_susie   <- FALSE
  n_triple    <- 0L
  triple_merged <- NULL

  if (!is.null(ld_result) && nrow(ld_result$vars) >= 10) {
    # merged has both POS (hg19) and pos_hg38 columns from the eQTL side
    # LD vars have ld_pos (hg19 positions)
    # Match on hg19 positions
    merged[, pos_hg19 := POS]  # eQTL native hg19 position
    ld_positions <- ld_result$vars$ld_pos

    # Find which merged rows have matching LD positions
    triple_idx <- which(merged$pos_hg19 %in% ld_positions)
    n_triple <- length(triple_idx)

    if (n_triple >= MIN_TRIPLE_SNPS) {
      triple_merged <- merged[triple_idx]

      # Subset LD matrix to matching positions (preserve order)
      ld_pos_match <- match(triple_merged$pos_hg19, ld_positions)

      # Remove any NAs from match
      valid_match <- !is.na(ld_pos_match)
      if (sum(valid_match) >= MIN_TRIPLE_SNPS) {
        triple_merged <- triple_merged[valid_match]
        ld_pos_match  <- ld_pos_match[valid_match]
        n_triple      <- nrow(triple_merged)

        # Subset and reorder LD matrix
        R <- ld_result$R[ld_pos_match, ld_pos_match, drop = FALSE]

        # FIX (review A04#1/B5#1/B2#2): re-sign LD to the harmonized GWAS effect-allele
        # orientation of the z-scores. PLINK 1.9 --r square signs r to .bim A1; where the
        # 1000G A1 != effect allele the off-diagonal signs are inverted relative to z,
        # corrupting susie_rss/coloc.susie. Flip discordant SNPs (diag(s) R diag(s)) and
        # drop allele-mismatched variants.
        ld_a1 <- ld_result$vars$ld_a1[ld_pos_match]
        ld_a2 <- ld_result$vars$ld_a2[ld_pos_match]
        ea_z  <- toupper(triple_merged$effect_allele)
        oa_z  <- toupper(triple_merged$other_allele)
        concord <- (ea_z == ld_a1 & oa_z == ld_a2) | (ea_z == ld_a2 & oa_z == ld_a1)
        sgn     <- ifelse(ea_z == ld_a1, 1, ifelse(ea_z == ld_a2, -1, NA_real_))
        concord[is.na(concord)] <- FALSE
        if (sum(concord) >= MIN_TRIPLE_SNPS) {
          triple_merged <- triple_merged[concord]
          R   <- R[concord, concord, drop = FALSE]
          sgn <- sgn[concord]
          R   <- sweep(sweep(R, 1, sgn, `*`), 2, sgn, `*`)  # diag(s) %*% R %*% diag(s)
          n_triple  <- nrow(triple_merged)
          use_susie <- TRUE
        } else {
          use_susie <- FALSE
        }
      }
    }
  }

  if (!use_susie && is.null(ld_result)) {
    n_ld_fail <- n_ld_fail + 1L
  }

  # ---- Initialize result placeholders ----
  pp_abf    <- rep(NA_real_, 5)
  pp_susie  <- rep(NA_real_, 5)
  n_cs_pairs <- NA_integer_
  method    <- "abf_only"
  top_snp   <- NA_character_
  top_snp_PP <- NA_real_

  # ---- (o) Always run ABF for comparison ----
  dataset1_abf <- list(
    snp     = merged$merge_key,
    beta    = merged$beta,
    varbeta = (merged$standard_error)^2,
    type    = GWAS_TYPE,
    N       = GWAS_N
  )
  if (GWAS_TYPE == "quant") {
    dataset1_abf$sdY <- 1
  } else {
    dataset1_abf$s <- GWAS_S
  }

  dataset2_abf <- list(
    snp     = merged$merge_key,
    beta    = merged$Beta,
    varbeta = (merged$SE)^2,
    type    = "quant",
    sdY     = 1,
    N       = EQTL_N
  )

  abf_res <- tryCatch(
    coloc.abf(dataset1_abf, dataset2_abf,
              p1 = COLOC_P1, p2 = COLOC_P2, p12 = COLOC_P12),
    error = function(e) {
      cat("  coloc.abf error:", gene_name, "-", conditionMessage(e), "\n")
      NULL
    }
  )

  if (is.null(abf_res)) next

  pp_abf <- c(abf_res$summary["PP.H0.abf"],
               abf_res$summary["PP.H1.abf"],
               abf_res$summary["PP.H2.abf"],
               abf_res$summary["PP.H3.abf"],
               abf_res$summary["PP.H4.abf"])

  # Extract top SNP from ABF
  top_snp_row <- abf_res$results[which.max(abf_res$results$SNP.PP.H4), ]
  top_snp     <- top_snp_row$snp
  top_snp_PP  <- top_snp_row$SNP.PP.H4

  # ---- (j-n) SuSiE fine-mapping + coloc.susie ----
  if (use_susie) {
    susie_success <- FALSE

    tryCatch({
      # Compute z-scores for the triple-intersect subset
      gwas_z  <- triple_merged$beta / triple_merged$standard_error
      eqtl_z  <- triple_merged$Beta / triple_merged$SE

      # SNP names are REQUIRED for coloc.susie — without them,
      # susie_rss lbf_variable has no colnames → coloc.bf_bf fails
      snp_ids <- paste0(gene_chr, ":", triple_merged$pos_hg19)
      names(gwas_z)  <- snp_ids
      names(eqtl_z)  <- snp_ids

      # LD regularization
      R_reg <- R + LD_REGULARIZE * diag(nrow(R))
      colnames(R_reg) <- snp_ids
      rownames(R_reg) <- snp_ids

      # SuSiE fine-mapping on GWAS
      # estimate_residual_variance=FALSE + residual_variance=1 required for
      # out-of-sample LD (1000G reference, not in-sample; see coloc vignette)
      s_gwas <- susie_rss(
        z       = gwas_z,
        R       = R_reg,
        n       = GWAS_N,
        L       = SUSIE_L,
        estimate_residual_variance = FALSE,
        check_R = FALSE
      )

      # SuSiE fine-mapping on eQTL
      s_eqtl <- susie_rss(
        z       = eqtl_z,
        R       = R_reg,
        n       = EQTL_N,
        L       = SUSIE_L,
        estimate_residual_variance = FALSE,
        check_R = FALSE
      )

      # coloc.susie
      susie_res <- coloc.susie(s_gwas, s_eqtl)

      if (!is.null(susie_res) && !is.null(susie_res$summary)) {
        # Extract best credible set pair PP.H4
        if (nrow(susie_res$summary) > 0) {
          best_row <- susie_res$summary[which.max(susie_res$summary$PP.H4.abf), ]
          pp_susie <- c(best_row$PP.H0.abf, best_row$PP.H1.abf,
                        best_row$PP.H2.abf, best_row$PP.H3.abf,
                        best_row$PP.H4.abf)
          n_cs_pairs <- nrow(susie_res$summary)
          method <- "susie"
          susie_success <- TRUE
        }
      }
    }, error = function(e) {
      cat("  SuSiE error:", gene_name, "-", conditionMessage(e), "\n")
    })

    # ---- (p) On SuSiE failure: flag as ABF fallback ----
    if (!susie_success) {
      method <- "abf_fallback"
      n_abf_fallback <- n_abf_fallback + 1L
    } else {
      n_susie_ok <- n_susie_ok + 1L
    }
  }

  # ---- Store results ----
  n_tested <- n_tested + 1L

  results_list[[gene_name]] <- data.table(
    gene          = gene_name,
    ensembl       = paste(unique(gene_eqtl$ENSG[gene_eqtl$ENSG != ""]), collapse = ";"),
    chr           = gene_chr,
    PP.H0.abf     = pp_abf[1],
    PP.H1.abf     = pp_abf[2],
    PP.H2.abf     = pp_abf[3],
    PP.H3.abf     = pp_abf[4],
    PP.H4.abf     = pp_abf[5],
    PP.H0.susie   = pp_susie[1],
    PP.H1.susie   = pp_susie[2],
    PP.H2.susie   = pp_susie[3],
    PP.H3.susie   = pp_susie[4],
    PP.H4.susie   = pp_susie[5],
    n_cs_pairs    = n_cs_pairs,
    n_snps_merged = n_snps_merged,
    n_snps_triple = as.integer(n_triple),
    method        = method,
    top_snp       = top_snp,
    top_snp_PP    = top_snp_PP
  )

  if (n_tested %% 50 == 0) {
    cat("  [", n_tested, "tested |", n_susie_ok, "susie |",
        n_abf_fallback, "abf_fallback |", n_ld_fail, "ld_fail |",
        n_low_snps, "low_snps |", n_skipped, "skipped ]\n")
  }

  # Checkpoint: write partial results every 50 genes
  if (n_tested - last_checkpoint >= checkpoint_interval) {
    cat("  Writing checkpoint at", n_tested, "genes...\n")
    tryCatch({
      checkpoint_dt <- rbindlist(results_list, fill = TRUE)
      fwrite(checkpoint_dt, checkpoint_file)
      last_checkpoint <- n_tested
    }, error = function(e) {
      cat("  WARNING: checkpoint write failed:", conditionMessage(e), "\n")
    })
  }
}

# ==============================================================================
# 6. Compile and save results
# ==============================================================================
cat("\n--- Step 6: Compiling results ---\n")

cat("  Genes tested:", n_tested, "\n")
cat("  SuSiE succeeded:", n_susie_ok, "\n")
cat("  ABF fallback (SuSiE error):", n_abf_fallback, "\n")
cat("  ABF only (LD fail/low triple):", n_tested - n_susie_ok - n_abf_fallback, "\n")
cat("  Genes skipped (no chr mapping):", n_skipped, "\n")
cat("  Genes skipped (< 10 shared SNPs):", n_low_snps, "\n")
cat("  Genes skipped (allele harmonization):", n_no_harm, "\n")
cat("  LD extraction failures:", n_ld_fail, "\n")

if (length(results_list) == 0) {
  cat("  No COLOC results.\n")
  fwrite(data.table(gene = character(), PP.H4.abf = numeric(), PP.H4.susie = numeric()),
         file.path(RESULTS_DIR, paste0("susie_coloc_results", CHR_SUFFIX, ".csv")))
  quit(save = "no", status = 0)
}

coloc_res <- rbindlist(results_list, fill = TRUE)
setorder(coloc_res, -PP.H4.abf)

fwrite(coloc_res, file.path(RESULTS_DIR, paste0("susie_coloc_results", CHR_SUFFIX, ".csv")))

cat("\n  COLOC Summary (Broadaway N=1,183 x", GWAS_NAME, "N=",
    format(GWAS_N, big.mark = ","), "):\n")
cat("    Genes tested:", n_tested, "\n")
cat("    ABF PP.H4 > 0.8:", nrow(coloc_res[PP.H4.abf > 0.8]), "\n")
cat("    ABF PP.H4 > 0.5:", nrow(coloc_res[PP.H4.abf > 0.5]), "\n")
cat("    ABF PP.H4 > 0.3:", nrow(coloc_res[PP.H4.abf > 0.3]), "\n")

susie_genes <- coloc_res[!is.na(PP.H4.susie)]
if (nrow(susie_genes) > 0) {
  cat("    SuSiE PP.H4 > 0.8:", nrow(susie_genes[PP.H4.susie > 0.8]), "\n")
  cat("    SuSiE PP.H4 > 0.5:", nrow(susie_genes[PP.H4.susie > 0.5]), "\n")
  cat("    SuSiE PP.H4 > 0.3:", nrow(susie_genes[PP.H4.susie > 0.3]), "\n")
}

# Top ABF hits
top_abf <- coloc_res[PP.H4.abf > 0.3]
if (nrow(top_abf) > 0) {
  cat("\n  Top ABF hits (PP.H4 > 0.3):\n")
  for (i in seq_len(min(nrow(top_abf), 20))) {
    susie_pp <- if (!is.na(top_abf$PP.H4.susie[i])) {
      sprintf(" | SuSiE=%.3f", top_abf$PP.H4.susie[i])
    } else ""
    cat("    ", top_abf$gene[i], ": ABF=", round(top_abf$PP.H4.abf[i], 3),
        susie_pp, " (", top_abf$n_snps_merged[i], "/",
        top_abf$n_snps_triple[i], " SNPs, ", top_abf$method[i], ")\n")
  }
}

# ABF vs SuSiE comparison
if (nrow(susie_genes) > 10) {
  cat("\n  ABF vs SuSiE comparison (", nrow(susie_genes), "genes with both):\n")
  rho <- cor(susie_genes$PP.H4.abf, susie_genes$PP.H4.susie,
             use = "complete.obs", method = "spearman")
  cat("    Spearman rho(PP.H4):", round(rho, 3), "\n")

  # Concordance at PP.H4 > 0.5 threshold
  abf_sig    <- susie_genes$PP.H4.abf > 0.5
  susie_sig  <- susie_genes$PP.H4.susie > 0.5
  both_sig   <- sum(abf_sig & susie_sig)
  abf_only   <- sum(abf_sig & !susie_sig)
  susie_only <- sum(!abf_sig & susie_sig)
  cat("    PP.H4 > 0.5: both=", both_sig, " ABF-only=", abf_only,
      " SuSiE-only=", susie_only, "\n")

  # Genes reclassified by SuSiE
  reclass_up   <- susie_genes[PP.H4.abf <= 0.5 & PP.H4.susie > 0.5]
  reclass_down <- susie_genes[PP.H4.abf > 0.5 & PP.H4.susie <= 0.5]
  if (nrow(reclass_up) > 0) {
    cat("    Reclassified UP by SuSiE (ABF<=0.5, SuSiE>0.5):\n")
    for (i in seq_len(min(nrow(reclass_up), 10))) {
      cat("      ", reclass_up$gene[i], ": ABF=", round(reclass_up$PP.H4.abf[i], 3),
          " -> SuSiE=", round(reclass_up$PP.H4.susie[i], 3), "\n")
    }
  }
  if (nrow(reclass_down) > 0) {
    cat("    Reclassified DOWN by SuSiE (ABF>0.5, SuSiE<=0.5):\n")
    for (i in seq_len(min(nrow(reclass_down), 10))) {
      cat("      ", reclass_down$gene[i], ": ABF=", round(reclass_down$PP.H4.abf[i], 3),
          " -> SuSiE=", round(reclass_down$PP.H4.susie[i], 3), "\n")
    }
  }
}

# Save summary
summary_dt <- data.table(
  metric = c("genes_tested", "susie_succeeded", "abf_fallback", "abf_only",
             "genes_skipped_no_chr", "genes_skipped_low_snps",
             "genes_skipped_harmonization", "ld_extraction_failures",
             "ABF_PP.H4_gt_0.8", "ABF_PP.H4_gt_0.5", "ABF_PP.H4_gt_0.3",
             "SuSiE_PP.H4_gt_0.8", "SuSiE_PP.H4_gt_0.5", "SuSiE_PP.H4_gt_0.3",
             "eqtl_source", "eqtl_N", "gwas_source", "gwas_N", "gwas_type",
             "ld_reference", "susie_L", "min_triple_snps"),
  value = c(n_tested, n_susie_ok, n_abf_fallback,
            n_tested - n_susie_ok - n_abf_fallback,
            n_skipped, n_low_snps, n_no_harm, n_ld_fail,
            nrow(coloc_res[PP.H4.abf > 0.8]),
            nrow(coloc_res[PP.H4.abf > 0.5]),
            nrow(coloc_res[PP.H4.abf > 0.3]),
            nrow(susie_genes[PP.H4.susie > 0.8]),
            nrow(susie_genes[PP.H4.susie > 0.5]),
            nrow(susie_genes[PP.H4.susie > 0.3]),
            "Broadaway2024", EQTL_N,
            GWAS_NAME, GWAS_N, GWAS_TYPE,
            "1000G_EUR_379", SUSIE_L, MIN_TRIPLE_SNPS)
)
fwrite(summary_dt, file.path(RESULTS_DIR, paste0("susie_summary", CHR_SUFFIX, ".csv")))

# ==============================================================================
# 7. Comparison with existing ABF-only results
# ==============================================================================
cat("\n--- Step 7: Cross-validation with existing ABF results ---\n")

# Try to load the matching ABF-only results (from 35b/35g/46)
if (grepl("^UKBB", GWAS_NAME)) {
  abf_dir <- file.path(BASE_DIR, "RNA-seq/results/causal_inference",
                        paste0("broadaway_", gwas_short))
} else if (grepl("^FINNGEN", GWAS_NAME)) {
  fg_suffix <- tolower(sub("FINNGEN_", "", GWAS_NAME))
  abf_dir <- file.path(BASE_DIR, "RNA-seq/results/causal_inference",
                        paste0("finngen_", fg_suffix))
} else {
  abf_dir <- NULL
}

if (!is.null(abf_dir)) {
  abf_file <- file.path(abf_dir, "coloc_results.csv")
  if (file.exists(abf_file)) {
    prev_abf <- fread(abf_file)
    if (nrow(prev_abf) > 0 && "PP.H4" %in% names(prev_abf)) {
      cat("  Previous ABF results found:", abf_file, "\n")
      cat("  Previous PP.H4 > 0.5:", nrow(prev_abf[PP.H4 > 0.5]), "\n")
      cat("  Current ABF PP.H4 > 0.5:", nrow(coloc_res[PP.H4.abf > 0.5]), "\n")

      # Correlation check
      shared <- merge(
        prev_abf[, .(gene, PP.H4_prev = PP.H4)],
        coloc_res[, .(gene, PP.H4_curr = PP.H4.abf)],
        by = "gene"
      )
      if (nrow(shared) > 10) {
        r <- cor(shared$PP.H4_prev, shared$PP.H4_curr, use = "complete.obs")
        cat("  Pearson r(PP.H4) vs previous:", round(r, 4), "\n")
        if (r < 0.99) {
          cat("  WARNING: ABF results differ from previous run (r < 0.99).\n")
          cat("  This may indicate GWAS/eQTL data changes.\n")
        }
      }
    }
  } else {
    cat("  No previous ABF results found at:", abf_file, "\n")
  }
}

# ==============================================================================
# 8. Figures
# ==============================================================================
cat("\n--- Step 8: Generating figures ---\n")

FIG_DIR <- file.path(RESULTS_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

if (nrow(coloc_res) > 0) {
  # (a) ABF PP.H4 distribution
  pdf(file.path(FIG_DIR, "abf_pp4_distribution.pdf"), width = 8, height = 5)
  p1 <- ggplot(coloc_res, aes(x = PP.H4.abf)) +
    geom_histogram(bins = 50, fill = "#2C7BB6", alpha = 0.8) +
    geom_vline(xintercept = c(0.5, 0.8), linetype = "dashed", color = c("orange", "red")) +
    labs(title = paste0("ABF COLOC PP.H4: Broadaway eQTL x ", GWAS_NAME),
         subtitle = paste0(n_tested, " genes tested; ",
                          nrow(coloc_res[PP.H4.abf > 0.5]), " with PP.H4 > 0.5"),
         x = "PP.H4 (ABF)", y = "Count") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 13))
  print(p1)
  dev.off()
  cat("  Saved: abf_pp4_distribution.pdf\n")

  # (b) ABF vs SuSiE scatter
  if (nrow(susie_genes) > 10) {
    pdf(file.path(FIG_DIR, "abf_vs_susie_scatter.pdf"), width = 7, height = 7)
    p2 <- ggplot(susie_genes, aes(x = PP.H4.abf, y = PP.H4.susie)) +
      geom_point(alpha = 0.4, size = 1.2, color = "#2C7BB6") +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
      geom_hline(yintercept = 0.5, linetype = "dotted", color = "orange", linewidth = 0.4) +
      geom_vline(xintercept = 0.5, linetype = "dotted", color = "orange", linewidth = 0.4) +
      labs(title = paste0("ABF vs SuSiE PP.H4: ", GWAS_NAME),
           subtitle = paste0(nrow(susie_genes), " genes with both methods"),
           x = "PP.H4 (ABF)", y = "PP.H4 (SuSiE, best CS pair)") +
      coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
      theme_bw(base_size = 11) +
      theme(plot.title = element_text(face = "bold", size = 13))
    print(p2)
    dev.off()
    cat("  Saved: abf_vs_susie_scatter.pdf\n")
  }

  # (c) Method breakdown
  method_counts <- coloc_res[, .N, by = method]
  pdf(file.path(FIG_DIR, "method_breakdown.pdf"), width = 6, height = 4)
  p3 <- ggplot(method_counts, aes(x = reorder(method, -N), y = N, fill = method)) +
    geom_col(width = 0.6) +
    scale_fill_manual(values = c("susie" = "#2C7BB6", "abf_fallback" = "#FDAE61",
                                  "abf_only" = "#D7191C")) +
    labs(title = paste0("SuSiE-COLOC Method Breakdown: ", GWAS_NAME),
         x = NULL, y = "Number of genes") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 13),
          legend.position = "none")
  print(p3)
  dev.off()
  cat("  Saved: method_breakdown.pdf\n")
}

cat("\n=== Script 35s: Complete ===\n")
cat("GWAS:", GWAS_NAME, "\n")
cat("Results:", RESULTS_DIR, "\n")
cat("End time:", format(Sys.time()), "\n")
