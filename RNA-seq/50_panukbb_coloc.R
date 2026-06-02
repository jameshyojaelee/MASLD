#!/usr/bin/env Rscript
# 50_panukbb_coloc.R
# ---------------------------------------------------------------------------
# Pan-UKBB Multi-Ancestry COLOC with Broadaway Liver eQTLs
#
# Cross-ancestry colocalization: European eQTLs (Broadaway N=1,183) x
# Pan-UKBB GWAS for AFR (~6.6K) and CSA (~8.9K) ancestry groups.
#
# ANCESTRY MISMATCH CAVEAT (same as Script 47 / BBJ):
#   European eQTLs with non-European GWAS → conservative (attenuated, not inflated).
#   ABF-only (no LD reference panels for AFR/CSA). Small N means fewer hits —
#   the attempt demonstrates thoroughness; null results are informative.
#
# Adapted from Script 47 (BBJ cross-ancestry COLOC).
#
# Usage:
#   GWAS_NAME=PANUKBB_AFR_ALT Rscript 50_panukbb_coloc.R
#   GWAS_NAME=PANUKBB_CSA_GGT Rscript 50_panukbb_coloc.R
#
# Inputs:
#   - data/broadaway_eqtl/{leads, chr* marginals}
#   - GWAS/MR_Data/PanUKBB/PanUKBB_{POP}_{TRAIT}_harmonised_hg38.tsv.gz
#
# Outputs:
#   - RNA-seq/results/causal_inference/panukbb_{pop}_{trait}/coloc_results.csv
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(coloc)
  library(data.table)
  library(rtracklayer)
  library(GenomicRanges)
  library(ggplot2)
})

# ==============================================================================
# Configuration
# ==============================================================================
GWAS_NAME <- Sys.getenv("GWAS_NAME", "PANUKBB_AFR_ALT")

gwas_config <- list(
  PANUKBB_AFR_ALT = list(file = "PanUKBB_AFR_ALT_harmonised_hg38.tsv.gz", N = 6636L,
                          type = "quant", ancestry = "African", pop = "AFR", trait = "ALT"),
  PANUKBB_AFR_AST = list(file = "PanUKBB_AFR_AST_harmonised_hg38.tsv.gz", N = 6636L,
                          type = "quant", ancestry = "African", pop = "AFR", trait = "AST"),
  PANUKBB_AFR_GGT = list(file = "PanUKBB_AFR_GGT_harmonised_hg38.tsv.gz", N = 6636L,
                          type = "quant", ancestry = "African", pop = "AFR", trait = "GGT"),
  PANUKBB_CSA_ALT = list(file = "PanUKBB_CSA_ALT_harmonised_hg38.tsv.gz", N = 8876L,
                          type = "quant", ancestry = "Central_South_Asian", pop = "CSA", trait = "ALT"),
  PANUKBB_CSA_AST = list(file = "PanUKBB_CSA_AST_harmonised_hg38.tsv.gz", N = 8876L,
                          type = "quant", ancestry = "Central_South_Asian", pop = "CSA", trait = "AST"),
  PANUKBB_CSA_GGT = list(file = "PanUKBB_CSA_GGT_harmonised_hg38.tsv.gz", N = 8876L,
                          type = "quant", ancestry = "Central_South_Asian", pop = "CSA", trait = "GGT")
)

if (!GWAS_NAME %in% names(gwas_config)) {
  stop("Unknown GWAS_NAME: ", GWAS_NAME,
       ". Valid: ", paste(names(gwas_config), collapse = ", "))
}

cfg <- gwas_config[[GWAS_NAME]]

BASE_DIR    <- Sys.getenv("MASLD_PROJECT_ROOT",
                          "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
EQTL_DIR   <- file.path(BASE_DIR, "data/broadaway_eqtl")
GWAS_FILE   <- file.path(BASE_DIR, "GWAS/MR_Data/PanUKBB", cfg$file)

out_suffix <- tolower(paste0("panukbb_", cfg$pop, "_", cfg$trait))
RESULTS_DIR <- file.path(BASE_DIR, "RNA-seq/results/causal_inference", out_suffix)
SENS_DIR    <- file.path(RESULTS_DIR, "coloc_sensitivity")
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SENS_DIR, recursive = TRUE, showWarnings = FALSE)

GWAS_N    <- cfg$N
GWAS_TYPE <- cfg$type
EQTL_N    <- 1183L

COLOC_P1  <- 1e-4
COLOC_P2  <- 1e-4
COLOC_P12 <- 5e-6
MIN_MAF   <- 0.01

# Minimum overlapping SNPs for COLOC (raised from 10; see threshold audit)
MIN_SNPS  <- 100L

CHAIN_FILE <- file.path(EQTL_DIR, "hg19ToHg38.over.chain")
CHAIN_GZ   <- paste0(CHAIN_FILE, ".gz")

cat("=== Script 50: Pan-UKBB Multi-Ancestry COLOC x Broadaway ===\n")
cat("GWAS:", GWAS_NAME, "(", cfg$ancestry, cfg$trait, ")\n")
cat("GWAS N:", format(GWAS_N, big.mark = ","), "\n")
cat("ANCESTRY MISMATCH: YES — EUR eQTLs x", cfg$ancestry, "GWAS\n")
cat("COLOC method: ABF only (no", cfg$ancestry, "LD panel)\n")
cat("NOTE: Small N (", format(GWAS_N, big.mark = ","), ") — expect fewer hits.\n")
cat("      Null results still contribute to cross-ancestry comparison.\n")
cat("Results dir:", RESULTS_DIR, "\n")
cat("Start time:", format(Sys.time()), "\n\n")

# ==============================================================================
# 1. Chain file
# ==============================================================================
cat("--- Step 1: Preparing liftover chain ---\n")
if (!file.exists(CHAIN_FILE)) {
  if (!file.exists(CHAIN_GZ)) {
    download.file(
      "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz",
      CHAIN_GZ, mode = "wb", quiet = FALSE)
  }
  system2("gunzip", args = c("-k", CHAIN_GZ))
}
chain <- import.chain(CHAIN_FILE)
cat("  Chain loaded:", length(chain), "chains\n")

# ==============================================================================
# 2. Load GWAS (already hg38 from format_panukbb_for_coloc.R)
# ==============================================================================
cat("\n--- Step 2: Loading Pan-UKBB GWAS ---\n")
if (!file.exists(GWAS_FILE)) {
  stop("GWAS file not found: ", GWAS_FILE,
       "\nRun GWAS/MR_Data/PanUKBB/format_panukbb_for_coloc.R first.")
}

gwas <- fread(GWAS_FILE)
cat("  Raw rows:", format(nrow(gwas), big.mark = ","), "\n")

# Validate required GWAS columns (B4: catch format changes early)
required_cols <- c("chromosome", "base_pair_location", "beta",
                   "standard_error", "p_value", "effect_allele", "other_allele")
missing <- setdiff(required_cols, names(gwas))
if (length(missing) > 0) {
  stop("GWAS missing columns: ", paste(missing, collapse = ", "),
       "\nFound: ", paste(names(gwas), collapse = ", "))
}

gwas[, chr := as.integer(chromosome)]
gwas[, pos_hg38 := as.integer(base_pair_location)]
gwas <- gwas[!is.na(chr) & chr %in% 1:22]
gwas <- gwas[!is.na(beta) & !is.na(standard_error) & standard_error > 0]

# MAF filter
if ("EAF" %in% names(gwas) && sum(!is.na(gwas$EAF)) > 0) {
  gwas[, maf := pmin(EAF, 1 - EAF)]
  n_before <- nrow(gwas)
  gwas <- gwas[is.na(maf) | maf >= MIN_MAF]
  cat("  MAF filter:", n_before - nrow(gwas), "removed\n")
}

gwas <- gwas[order(chr, pos_hg38, p_value)]
gwas <- gwas[!duplicated(paste(chr, pos_hg38))]
gwas[, merge_key := paste0(chr, ":", pos_hg38)]

cat("  GWAS after QC:", format(nrow(gwas), big.mark = ","), "\n")
cat("  GW-significant (p<5e-8):", sum(gwas$p_value < 5e-8, na.rm = TRUE), "\n")

# ==============================================================================
# 3. Load Broadaway leads
# ==============================================================================
cat("\n--- Step 3: Loading Broadaway eGene leads ---\n")
leads <- fread(file.path(EQTL_DIR, "Liver_eQTL_Meta_Leads_ST3_20240530.tsv"))
all_eGenes <- unique(leads$Gene)
cat("  Unique eGenes:", length(all_eGenes), "\n")

leads[, eqtl_chr := as.integer(sub("_.*", "", Variant))]
gene_chr_map <- unique(leads[!is.na(eqtl_chr), .(Gene, chr = eqtl_chr)])
gene_chr_map <- gene_chr_map[, .N, by = .(Gene, chr)][order(-N)][!duplicated(Gene)]
leads_lookup <- unique(leads[, .(Gene, Ensembl)])

# ==============================================================================
# 4. Helper functions (identical to Scripts 46/47/49)
# ==============================================================================
liftover_positions <- function(chr_num, positions) {
  gr <- GRanges(seqnames = paste0("chr", chr_num), ranges = IRanges(start = positions, width = 1))
  lifted <- liftOver(gr, chain)
  n_mapped <- lengths(lifted)
  hg38_pos <- rep(NA_integer_, length(positions))
  idx <- which(n_mapped == 1L)
  if (length(idx) > 0) hg38_pos[idx] <- start(unlist(lifted[idx]))
  hg38_pos
}

chr_eqtl_cache <- list()
load_chr_eqtl <- function(chr_num) {
  key <- as.character(chr_num)
  if (!is.null(chr_eqtl_cache[[key]])) return(chr_eqtl_cache[[key]])
  fname <- file.path(EQTL_DIR, paste0("chr", chr_num, "_marginal_summary_results.tsv"))
  if (!file.exists(fname)) return(NULL)
  cat("    Loading chr", chr_num, "eQTL data...\n")
  dt <- fread(fname)
  cat("    Lifting over", nrow(dt), "positions...\n")
  dt[, pos_hg38 := liftover_positions(chr_num, POS)]
  n_before_lift <- nrow(dt)
  dt <- dt[!is.na(pos_hg38)]
  cat("    Liftover:", nrow(dt), "/", n_before_lift,
      "(", round(100 * nrow(dt) / n_before_lift, 1), "% )\n")
  dt[, merge_key := paste0(CHR, ":", pos_hg38)]
  chr_eqtl_cache[[key]] <<- dt
  return(dt)
}

harmonize_alleles <- function(dt) {
  dt <- copy(dt)
  comp <- c("A" = "T", "T" = "A", "C" = "G", "G" = "C")
  dt[, match_type := "none"]
  dt[EA == effect_allele & NEA == other_allele, match_type := "direct"]
  dt[EA == other_allele & NEA == effect_allele, match_type := "flipped"]
  # Remove strand-ambiguous variants that failed direct/flipped match
  # (complement matching unreliable for A/T and C/G pairs)
  dt[, is_ambiguous := (EA %in% c("A","T") & NEA %in% c("A","T")) |
                        (EA %in% c("C","G") & NEA %in% c("C","G"))]
  n_ambig_unmatched <- sum(dt$is_ambiguous & dt$match_type == "none", na.rm = TRUE)
  cat("    Strand-ambiguous unmatched:", n_ambig_unmatched, "\n")
  dt <- dt[!(is_ambiguous & match_type == "none")]
  dt[match_type == "none" & nchar(EA) == 1 & nchar(NEA) == 1,
     match_type := fifelse(
       comp[EA] == effect_allele & comp[NEA] == other_allele, "direct",
       fifelse(comp[EA] == other_allele & comp[NEA] == effect_allele, "flipped", "none"))]
  cat("    Harmonization: direct=", sum(dt$match_type == "direct"),
      " flipped=", sum(dt$match_type == "flipped"),
      " dropped=", sum(dt$match_type == "none"), "\n")
  dt <- dt[match_type != "none"]
  if (any(dt$match_type == "flipped")) dt[match_type == "flipped", Beta := -Beta]
  dt[, is_ambiguous := NULL]
  return(dt)
}

# ==============================================================================
# 5. Run COLOC per eGene
# ==============================================================================
cat("\n--- Step 5: Running COLOC ---\n")
cat("  Total eGenes:", length(all_eGenes), "\n\n")

results_list <- list()
n_tested <- 0L; n_skipped <- 0L; n_low_snps <- 0L; n_no_harm <- 0L

for (gene_name in all_eGenes) {
  chr_info <- gene_chr_map[Gene == gene_name]
  if (nrow(chr_info) == 0) { n_skipped <- n_skipped + 1L; next }
  gene_chr <- chr_info$chr[1]

  eqtl_data <- load_chr_eqtl(gene_chr)
  if (is.null(eqtl_data)) { n_skipped <- n_skipped + 1L; next }

  gene_eqtl <- eqtl_data[GeneSymbol == gene_name]
  if (nrow(gene_eqtl) == 0) {
    gene_ensg <- leads_lookup$Ensembl[leads_lookup$Gene == gene_name]
    if (length(gene_ensg) > 0) gene_eqtl <- eqtl_data[ENSG %in% gene_ensg]
  }
  if (nrow(gene_eqtl) < MIN_SNPS) { n_low_snps <- n_low_snps + 1L; next }

  merged <- merge(gene_eqtl, gwas, by = "merge_key", suffixes = c(".eqtl", ".gwas"))
  if (nrow(merged) < MIN_SNPS) { n_low_snps <- n_low_snps + 1L; next }

  merged <- merged[!is.na(Beta) & !is.na(SE) & SE > 0 &
                   !is.na(beta) & !is.na(standard_error) & standard_error > 0]
  if (nrow(merged) < MIN_SNPS) { n_low_snps <- n_low_snps + 1L; next }

  merged <- harmonize_alleles(merged)
  if (nrow(merged) < MIN_SNPS) { n_no_harm <- n_no_harm + 1L; next }

  # Cross-population EAF concordance diagnostic (M4)
  if (n_tested == 0L) {
    # Report once for the first successfully tested gene and save to file
    if ("EAF.eqtl" %in% names(merged) && "EAF.gwas" %in% names(merged)) {
      eaf_comp <- merged[!is.na(EAF.eqtl) & !is.na(EAF.gwas)]
      if (nrow(eaf_comp) > 10) {
        eaf_r <- round(cor(eaf_comp$EAF.eqtl, eaf_comp$EAF.gwas), 3)
        eaf_n_large <- sum(abs(eaf_comp$EAF.eqtl - eaf_comp$EAF.gwas) > 0.3)
        cat("    EAF concordance (EUR vs", cfg$ancestry, "):",
            "r =", eaf_r, ", |diff|>0.3:", eaf_n_large, "\n")
        # Save to file for manuscript reporting
        eaf_out <- data.table(
          ancestry = cfg$ancestry, gene = gene_name,
          n_variants = nrow(eaf_comp), pearson_r = eaf_r,
          n_large_diff = eaf_n_large)
        eaf_file <- file.path(RESULTS_DIR, "eaf_concordance_first_gene.csv")
        fwrite(eaf_out, eaf_file)
      }
    }
  }

  dataset1 <- list(snp = merged$merge_key, beta = merged$beta,
                   varbeta = (merged$standard_error)^2,
                   type = "quant", sdY = 1, N = GWAS_N)
  dataset2 <- list(snp = merged$merge_key, beta = merged$Beta,
                   varbeta = (merged$SE)^2,
                   type = "quant", sdY = 1, N = EQTL_N)

  my.res <- tryCatch(
    coloc.abf(dataset1, dataset2, p1 = COLOC_P1, p2 = COLOC_P2, p12 = COLOC_P12),
    error = function(e) { cat("  Error:", gene_name, conditionMessage(e), "\n"); NULL })

  if (!is.null(my.res)) {
    n_tested <- n_tested + 1L
    top_snp_row <- my.res$results[which.max(my.res$results$SNP.PP.H4), ]
    results_list[[gene_name]] <- data.table(
      gene = gene_name,
      ensembl = paste(unique(gene_eqtl$ENSG[gene_eqtl$ENSG != ""]), collapse = ";"),
      chr = gene_chr,
      PP.H0 = my.res$summary["PP.H0.abf"], PP.H1 = my.res$summary["PP.H1.abf"],
      PP.H2 = my.res$summary["PP.H2.abf"], PP.H3 = my.res$summary["PP.H3.abf"],
      PP.H4 = my.res$summary["PP.H4.abf"],
      n_snps = nrow(merged), n_eqtl_snps = nrow(gene_eqtl),
      top_snp = top_snp_row$snp, top_snp_PP = top_snp_row$SNP.PP.H4,
      eqtl_source = "Broadaway2024", eqtl_N = EQTL_N, eqtl_ancestry = "European",
      gwas_source = GWAS_NAME, gwas_ancestry = cfg$ancestry,
      gwas_N = GWAS_N, gwas_type = "quant",
      ancestry_match = FALSE, method = "abf_beta_varbeta"
    )

    if (my.res$summary["PP.H4.abf"] > 0.3) {
      tryCatch({
        pdf(file.path(SENS_DIR, paste0("sensitivity_", gene_name, ".pdf")), width = 10, height = 8)
        sensitivity(my.res, rule = "H4 > 0.5"); dev.off()
      }, error = function(e) tryCatch(dev.off(), error = function(e2) NULL))
    }

    if (n_tested %% 200 == 0) cat("  [", n_tested, "tested]\n")
  }
}

# ==============================================================================
# 6. Results
# ==============================================================================
cat("\n--- Step 6: Results ---\n")
cat("  Tested:", n_tested, " Skipped:", n_skipped,
    " Low-SNP:", n_low_snps, " No-harmonize:", n_no_harm, "\n")

if (length(results_list) == 0) {
  cat("  No results (expected for small N).\n")
  fwrite(data.table(gene = character(), PP.H4 = numeric()),
         file.path(RESULTS_DIR, "coloc_results.csv"))
  quit(save = "no", status = 0)
}

coloc_res <- rbindlist(results_list, fill = TRUE)
setorder(coloc_res, -PP.H4)
fwrite(coloc_res, file.path(RESULTS_DIR, "coloc_results.csv"))

cat("  PP.H4 > 0.8:", nrow(coloc_res[PP.H4 > 0.8]), "\n")
cat("  PP.H4 > 0.5:", nrow(coloc_res[PP.H4 > 0.5]), "\n")
cat("  PP.H4 > 0.3:", nrow(coloc_res[PP.H4 > 0.3]), "\n")
cat("  NOTE: Fewer hits expected due to small N and ancestry mismatch.\n")

# Summary
summary_dt <- data.table(
  metric = c("genes_tested", "PP.H4_gt_0.5", "PP.H4_gt_0.3",
             "gwas_source", "gwas_ancestry", "gwas_N",
             "eqtl_ancestry", "ancestry_match", "maf_filter"),
  value = c(n_tested, nrow(coloc_res[PP.H4 > 0.5]), nrow(coloc_res[PP.H4 > 0.3]),
            GWAS_NAME, cfg$ancestry, GWAS_N,
            "European", "FALSE", MIN_MAF)
)
fwrite(summary_dt, file.path(RESULTS_DIR, "coloc_summary.csv"))

cat("\n=== Script 50: Complete ===\n")
cat("End time:", format(Sys.time()), "\n")
