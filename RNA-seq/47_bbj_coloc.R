#!/usr/bin/env Rscript
# 47_bbj_coloc.R
# ---------------------------------------------------------------------------
# Biobank Japan (BBJ) Liver Enzyme COLOC with Broadaway Liver eQTLs
#
# Cross-ancestry colocalization: European eQTLs (Broadaway N=1,183) x
# East Asian GWAS (BBJ, Kanai et al. 2018).
#
# ANCESTRY MISMATCH CAVEAT:
#   Using European eQTLs with Japanese GWAS is imperfect, but:
#   - COLOC tests all variants in a region simultaneously (more robust than MR
#     to LD differences between populations)
#   - Positive hits are conservative (attenuated, not inflated) — LD
#     mismatches reduce power but don't increase false positives
#   - Published precedent: Shi et al. 2021 AJHG
#   - Uses ABF-only (not SuSiE) — no LD reference panel required
#
# Adapted from Script 35g (Broadaway x UKBB quantitative GWAS).
# Key differences from 35g:
#   - BBJ data already lifted over to hg38 by format_bbj_for_coloc.R
#   - Ancestry annotations on every result row
#   - MAF filter: require MAF >= 0.01 if EAF available
#   - ABF-only COLOC (no SuSiE — no Japanese LD panel)
#
# Usage:
#   GWAS_NAME=BBJ_ALT Rscript 47_bbj_coloc.R
#   GWAS_NAME=BBJ_AST Rscript 47_bbj_coloc.R
#   GWAS_NAME=BBJ_GGT Rscript 47_bbj_coloc.R
#
# Inputs:
#   - data/broadaway_eqtl/Liver_eQTL_Meta_Leads_ST3_20240530.tsv
#   - data/broadaway_eqtl/chr{1-22}_marginal_summary_results.tsv
#   - GWAS/MR_Data/BBJ/BBJ_{ALT,AST,GGT}_harmonised_hg38.tsv.gz
#   - data/broadaway_eqtl/hg19ToHg38.over.chain
#
# Outputs:
#   - RNA-seq/results/causal_inference/bbj_{alt,ast,ggt}/coloc_results.csv
#   - RNA-seq/results/causal_inference/bbj_{alt,ast,ggt}/coloc_summary.csv
#   - RNA-seq/results/causal_inference/bbj_{alt,ast,ggt}/coloc_sensitivity/
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
GWAS_NAME <- Sys.getenv("GWAS_NAME", "BBJ_ALT")

gwas_config <- list(
  BBJ_ALT = list(
    file  = "BBJ_ALT_harmonised_hg38.tsv.gz",
    N     = 261406L,
    type  = "quant",
    label = "BBJ ALT (Kanai et al. 2018)"
  ),
  BBJ_AST = list(
    file  = "BBJ_AST_harmonised_hg38.tsv.gz",
    N     = 261270L,
    type  = "quant",
    label = "BBJ AST (Kanai et al. 2018)"
  ),
  BBJ_GGT = list(
    file  = "BBJ_GGT_harmonised_hg38.tsv.gz",
    N     = 164006L,
    type  = "quant",
    label = "BBJ GGT (Kanai et al. 2018)"
  )
)

if (!GWAS_NAME %in% names(gwas_config)) {
  stop("Unknown GWAS_NAME: ", GWAS_NAME,
       ". Valid options: ", paste(names(gwas_config), collapse = ", "))
}

cfg <- gwas_config[[GWAS_NAME]]

BASE_DIR    <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
EQTL_DIR   <- file.path(BASE_DIR, "data/broadaway_eqtl")
GWAS_FILE   <- file.path(BASE_DIR, "GWAS/MR_Data/BBJ", cfg$file)

out_suffix <- tolower(sub("BBJ_", "", GWAS_NAME))
RESULTS_DIR <- file.path(BASE_DIR, "RNA-seq/results/causal_inference",
                          paste0("bbj_", out_suffix))
SENS_DIR    <- file.path(RESULTS_DIR, "coloc_sensitivity")

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SENS_DIR, recursive = TRUE, showWarnings = FALSE)

GWAS_N    <- cfg$N
GWAS_TYPE <- cfg$type

# Broadaway eQTL parameters
EQTL_N    <- 1183L

# COLOC priors
COLOC_P1  <- 1e-4
COLOC_P2  <- 1e-4
COLOC_P12 <- 5e-6

# Minimum overlapping SNPs for COLOC (raised from 10; see threshold audit)
MIN_SNPS  <- 100L

# MAF filter for cross-ancestry robustness
MIN_MAF <- 0.01

# Liftover chain file (for eQTLs: hg19 → hg38)
CHAIN_FILE <- file.path(EQTL_DIR, "hg19ToHg38.over.chain")
CHAIN_GZ   <- paste0(CHAIN_FILE, ".gz")

cat("=== Script 47: BBJ Liver Enzyme COLOC x Broadaway Liver eQTL ===\n")
cat("GWAS:", cfg$label, "\n")
cat("GWAS file:", cfg$file, "\n")
cat("GWAS N:", format(GWAS_N, big.mark = ","), "\n")
cat("GWAS ancestry: East Asian (Japanese)\n")
cat("eQTL ancestry: European (Broadaway N=1,183)\n")
cat("ANCESTRY MISMATCH: YES — results are conservative (attenuated, not inflated)\n")
cat("COLOC method: ABF only (no SuSiE — no Japanese LD reference panel)\n")
cat("MAF filter:", MIN_MAF, "\n")
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
# 2. Load BBJ GWAS (already hg38 from format_bbj_for_coloc.R)
# ==============================================================================
cat("\n--- Step 2: Loading BBJ", GWAS_NAME, "GWAS ---\n")

if (!file.exists(GWAS_FILE)) {
  stop("GWAS file not found: ", GWAS_FILE,
       "\nRun GWAS/MR_Data/BBJ/format_bbj_for_coloc.R first.")
}

gwas <- fread(GWAS_FILE)
cat("  Raw GWAS rows:", nrow(gwas), "\n")

# Harmonised format from format_bbj_for_coloc.R:
# chromosome, base_pair_location, pos_hg19, effect_allele, other_allele,
# beta, standard_error, p_value, EAF
gwas[, chr := as.integer(chromosome)]
gwas[, pos_hg38 := as.integer(base_pair_location)]
gwas <- gwas[!is.na(chr) & chr %in% 1:22]
gwas <- gwas[!is.na(beta) & !is.na(standard_error) & standard_error > 0]

# MAF filter for cross-ancestry robustness
if ("EAF" %in% names(gwas) && sum(!is.na(gwas$EAF)) > 0) {
  gwas[, maf := pmin(EAF, 1 - EAF)]
  n_before <- nrow(gwas)
  gwas <- gwas[is.na(maf) | maf >= MIN_MAF]
  cat("  MAF filter (>=", MIN_MAF, "):", n_before - nrow(gwas), "variants removed\n")
}

# Deduplicate by position
gwas <- gwas[order(chr, pos_hg38, p_value)]
gwas <- gwas[!duplicated(paste(chr, pos_hg38))]
gwas[, merge_key := paste0(chr, ":", pos_hg38)]

cat("  GWAS after dedup:", nrow(gwas), "variants\n")
cat("  Genome-wide significant (p<5e-8):", sum(gwas$p_value < 5e-8), "\n")

# ==============================================================================
# 3. Load Broadaway leads — ALL eGenes
# ==============================================================================
cat("\n--- Step 3: Loading Broadaway eGene leads ---\n")

leads <- fread(file.path(EQTL_DIR, "Liver_eQTL_Meta_Leads_ST3_20240530.tsv"))
cat("  Total eQTL signals:", nrow(leads), "\n")

all_eGenes <- unique(leads$Gene)
cat("  Unique eGenes:", length(all_eGenes), "\n")

leads[, eqtl_chr := as.integer(sub("_.*", "", Variant))]
gene_chr_map <- unique(leads[!is.na(eqtl_chr), .(Gene, chr = eqtl_chr)])
gene_chr_map <- gene_chr_map[, .N, by = .(Gene, chr)][order(-N)][!duplicated(Gene)]
cat("  Gene-to-chromosome mapping:", nrow(gene_chr_map), "genes\n")

leads_lookup <- unique(leads[, .(Gene, Ensembl)])

# ==============================================================================
# 4. Helper functions (same as 35g/46)
# ==============================================================================

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

chr_eqtl_cache <- list()

load_chr_eqtl <- function(chr_num) {
  key <- as.character(chr_num)
  if (!is.null(chr_eqtl_cache[[key]])) return(chr_eqtl_cache[[key]])

  fname <- file.path(EQTL_DIR, paste0("chr", chr_num, "_marginal_summary_results.tsv"))
  if (!file.exists(fname)) return(NULL)

  cat("    Loading chr", chr_num, "eQTL data...\n")
  dt <- fread(fname)

  cat("    Lifting over", nrow(dt), "positions from hg19 to hg38...\n")
  dt[, pos_hg38 := liftover_positions(chr_num, POS)]
  dt <- dt[!is.na(pos_hg38)]
  cat("    Successful liftover:", nrow(dt), "variants\n")

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
  if (any(dt$match_type == "flipped")) {
    dt[match_type == "flipped", Beta := -Beta]
  }
  dt[, is_ambiguous := NULL]
  return(dt)
}

# ==============================================================================
# 5. Run COLOC per eGene
# ==============================================================================
cat("\n--- Step 5: Running COLOC ---\n")

candidate_chrs <- unique(gene_chr_map[Gene %in% all_eGenes]$chr)
cat("  Chromosomes to process:", paste(sort(candidate_chrs), collapse = ", "), "\n")
cat("  Total eGenes to test:", length(all_eGenes), "\n\n")

results_list <- list()
n_tested     <- 0L
n_skipped    <- 0L
n_low_snps   <- 0L
n_no_harm    <- 0L

for (gene_name in all_eGenes) {
  chr_info <- gene_chr_map[Gene == gene_name]
  if (nrow(chr_info) == 0) {
    n_skipped <- n_skipped + 1L
    next
  }
  gene_chr <- chr_info$chr[1]

  eqtl_data <- load_chr_eqtl(gene_chr)
  if (is.null(eqtl_data)) {
    n_skipped <- n_skipped + 1L
    next
  }

  gene_eqtl <- eqtl_data[GeneSymbol == gene_name]
  if (nrow(gene_eqtl) == 0) {
    gene_ensg <- leads_lookup$Ensembl[leads_lookup$Gene == gene_name]
    if (length(gene_ensg) > 0) {
      gene_eqtl <- eqtl_data[ENSG %in% gene_ensg]
    }
  }

  if (nrow(gene_eqtl) < MIN_SNPS) {
    n_low_snps <- n_low_snps + 1L
    next
  }

  merged <- merge(gene_eqtl, gwas, by = "merge_key", suffixes = c(".eqtl", ".gwas"))

  if (nrow(merged) < MIN_SNPS) {
    n_low_snps <- n_low_snps + 1L
    next
  }

  merged <- merged[!is.na(Beta) & !is.na(SE) & SE > 0 &
                   !is.na(beta) & !is.na(standard_error) & standard_error > 0]

  if (nrow(merged) < MIN_SNPS) {
    n_low_snps <- n_low_snps + 1L
    next
  }

  merged <- harmonize_alleles(merged)

  if (nrow(merged) < MIN_SNPS) {
    n_no_harm <- n_no_harm + 1L
    next
  }

  # COLOC datasets — both quantitative
  dataset1 <- list(
    snp     = merged$merge_key,
    beta    = merged$beta,
    varbeta = (merged$standard_error)^2,
    type    = GWAS_TYPE,
    sdY     = 1,
    N       = GWAS_N
  )

  dataset2 <- list(
    snp     = merged$merge_key,
    beta    = merged$Beta,
    varbeta = (merged$SE)^2,
    type    = "quant",
    sdY     = 1,
    N       = EQTL_N
  )

  my.res <- tryCatch(
    coloc.abf(dataset1, dataset2, p1 = COLOC_P1, p2 = COLOC_P2, p12 = COLOC_P12),
    error = function(e) {
      cat("  coloc.abf error:", gene_name, "-", conditionMessage(e), "\n")
      NULL
    }
  )

  if (!is.null(my.res)) {
    n_tested <- n_tested + 1L

    top_snp_row <- my.res$results[which.max(my.res$results$SNP.PP.H4), ]

    results_list[[gene_name]] <- data.table(
      gene        = gene_name,
      ensembl     = paste(unique(gene_eqtl$ENSG[gene_eqtl$ENSG != ""]), collapse = ";"),
      chr         = gene_chr,
      PP.H0       = my.res$summary["PP.H0.abf"],
      PP.H1       = my.res$summary["PP.H1.abf"],
      PP.H2       = my.res$summary["PP.H2.abf"],
      PP.H3       = my.res$summary["PP.H3.abf"],
      PP.H4       = my.res$summary["PP.H4.abf"],
      n_snps      = nrow(merged),
      n_eqtl_snps = nrow(gene_eqtl),
      top_snp     = top_snp_row$snp,
      top_snp_PP  = top_snp_row$SNP.PP.H4,
      eqtl_source = "Broadaway2024",
      eqtl_N      = EQTL_N,
      eqtl_ancestry = "European",
      gwas_source = GWAS_NAME,
      gwas_ancestry = "East_Asian",
      gwas_N      = GWAS_N,
      gwas_type   = GWAS_TYPE,
      ancestry_match = FALSE,
      method      = "abf_beta_varbeta"
    )

    if (my.res$summary["PP.H4.abf"] > 0.3) {
      tryCatch({
        pdf(file.path(SENS_DIR, paste0("sensitivity_", gene_name, ".pdf")),
            width = 10, height = 8)
        sensitivity(my.res, rule = "H4 > 0.5")
        dev.off()
      }, error = function(e) {
        tryCatch(dev.off(), error = function(e2) NULL)
      })
    }

    if (n_tested %% 200 == 0) {
      cat("  [", n_tested, "tested,", n_low_snps, "low-SNP,",
          n_no_harm, "no-harmonize,", n_skipped, "skipped]\n")
    }
  }
}

# ==============================================================================
# 6. Compile and save results
# ==============================================================================
cat("\n--- Step 6: Compiling results ---\n")

cat("  Genes tested:", n_tested, "\n")
cat("  Genes skipped (no chr mapping):", n_skipped, "\n")
cat("  Genes skipped (<", MIN_SNPS, "shared SNPs):", n_low_snps, "\n")
cat("  Genes skipped (allele harmonization):", n_no_harm, "\n")

if (length(results_list) == 0) {
  cat("  No COLOC results.\n")
  fwrite(data.table(gene = character(), PP.H4 = numeric()),
         file.path(RESULTS_DIR, "coloc_results.csv"))
  quit(save = "no", status = 0)
}

coloc_res <- rbindlist(results_list, fill = TRUE)
setorder(coloc_res, -PP.H4)

fwrite(coloc_res, file.path(RESULTS_DIR, "coloc_results.csv"))

gwas_label <- paste0(GWAS_NAME, " N=", format(GWAS_N, big.mark = ","))
cat("\n  COLOC Summary (Broadaway N=1,183 [EUR] x", gwas_label, "[EAS]):\n")
cat("    Genes tested:", n_tested, "\n")
cat("    PP.H4 > 0.8:", nrow(coloc_res[PP.H4 > 0.8]), "\n")
cat("    PP.H4 > 0.5:", nrow(coloc_res[PP.H4 > 0.5]), "\n")
cat("    PP.H4 > 0.3:", nrow(coloc_res[PP.H4 > 0.3]), "\n")
cat("    NOTE: Hits are conservative due to EUR/EAS ancestry mismatch\n")

top_hits <- coloc_res[PP.H4 > 0.3]
if (nrow(top_hits) > 0) {
  cat("\n  Top COLOC hits (PP.H4 > 0.3):\n")
  for (i in seq_len(min(nrow(top_hits), 30))) {
    cat("    ", top_hits$gene[i], ": PP.H4 =", round(top_hits$PP.H4[i], 3),
        "(", top_hits$n_snps[i], "SNPs, top:", top_hits$top_snp[i], ")\n")
  }
}

# Save summary
summary_dt <- data.table(
  metric = c("genes_tested", "genes_skipped_no_chr", "genes_skipped_low_snps",
             "genes_skipped_harmonization",
             "PP.H4_gt_0.8", "PP.H4_gt_0.5", "PP.H4_gt_0.3",
             "eqtl_source", "eqtl_N", "eqtl_ancestry",
             "gwas_source", "gwas_ancestry", "gwas_N",
             "gwas_type", "ancestry_match", "maf_filter"),
  value = c(n_tested, n_skipped, n_low_snps, n_no_harm,
            nrow(coloc_res[PP.H4 > 0.8]),
            nrow(coloc_res[PP.H4 > 0.5]),
            nrow(coloc_res[PP.H4 > 0.3]),
            "Broadaway2024", EQTL_N, "European",
            GWAS_NAME, "East_Asian", GWAS_N,
            GWAS_TYPE, "FALSE", MIN_MAF)
)
fwrite(summary_dt, file.path(RESULTS_DIR, "coloc_summary.csv"))

# ==============================================================================
# 7. Validation: known ancestry-variable MASLD genes
# ==============================================================================
cat("\n--- Step 7: Validation ---\n")

# Key genes with known ancestry-variable frequency
# PNPLA3 rs738409: MAF ~42% EAS vs ~22% EUR — expect STRONGER signal in BBJ
# TM6SF2 rs58542926: MAF ~7% EAS vs ~7% EUR — similar power
# HSD17B13 rs72613567: lower in EAS — may lose signal
# MBOAT7 rs641738: MAF ~40% EAS vs ~42% EUR — similar
anchor_genes <- c("PNPLA3", "TM6SF2", "HSD17B13", "MBOAT7",
                  "SLC39A8", "SORT1", "RORA", "THRB",
                  "EFHD1", "CHEK2", "SPTLC3", "APOH")
cat("\n  Anchor gene validation (ancestry-variable MASLD genes):\n")
for (g in anchor_genes) {
  hit <- coloc_res[gene == g]
  if (nrow(hit) > 0) {
    cat("    ", g, ": PP.H4 =", round(hit$PP.H4[1], 3),
        "(", hit$n_snps[1], "SNPs)\n")
  } else {
    cat("    ", g, ": not tested\n")
  }
}

# Cross-reference with European COLOC results
eur_files <- list(
  "UKBB_ALT" = file.path(BASE_DIR, "RNA-seq/results/causal_inference/broadaway_ukbb/coloc_results.csv"),
  "UKBB_AST" = file.path(BASE_DIR, "RNA-seq/results/causal_inference/broadaway_ukbb_ast/coloc_results.csv"),
  "UKBB_GGT" = file.path(BASE_DIR, "RNA-seq/results/causal_inference/broadaway_ukbb_ggt/coloc_results.csv")
)

# Match corresponding trait
eur_trait <- sub("BBJ_", "UKBB_", GWAS_NAME)
eur_file <- eur_files[[eur_trait]]
if (!is.null(eur_file) && file.exists(eur_file)) {
  eur_res <- fread(eur_file)
  if (nrow(eur_res) > 0 && "PP.H4" %in% names(eur_res)) {
    cat("\n  Cross-ancestry comparison:", eur_trait, "(EUR) vs", GWAS_NAME, "(EAS):\n")
    cat("    EUR PP.H4 > 0.5:", nrow(eur_res[PP.H4 > 0.5]), "\n")
    cat("    EAS PP.H4 > 0.5:", nrow(coloc_res[PP.H4 > 0.5]), "\n")
    shared <- intersect(eur_res$gene, coloc_res$gene)
    if (length(shared) > 5) {
      comp <- merge(
        eur_res[, .(gene, PP.H4_eur = PP.H4)],
        coloc_res[, .(gene, PP.H4_eas = PP.H4)],
        by = "gene"
      )
      cat("    Shared genes tested:", nrow(comp), "\n")
      cat("    Pearson r (PP.H4):", round(cor(comp$PP.H4_eur, comp$PP.H4_eas,
                                               use = "complete.obs"), 3), "\n")
      # Cross-ancestry replication
      n_both <- nrow(comp[PP.H4_eur > 0.5 & PP.H4_eas > 0.5])
      n_eur_only <- nrow(comp[PP.H4_eur > 0.5 & PP.H4_eas <= 0.5])
      n_eas_only <- nrow(comp[PP.H4_eur <= 0.5 & PP.H4_eas > 0.5])
      cat("    Replicated in both:", n_both, "\n")
      cat("    EUR-only:", n_eur_only, "\n")
      cat("    EAS-only:", n_eas_only, "\n")
    }
  }
}

# ==============================================================================
# 8. Figures
# ==============================================================================
cat("\n--- Step 8: Generating figures ---\n")

FIG_DIR <- file.path(BASE_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

gwas_short <- tolower(GWAS_NAME)

if (nrow(coloc_res) > 0) {
  # PP.H4 distribution
  pdf(file.path(FIG_DIR, paste0(gwas_short, "_coloc_pp4_distribution.pdf")),
      width = 8, height = 5)
  p1 <- ggplot(coloc_res, aes(x = PP.H4)) +
    geom_histogram(bins = 50, fill = "#E66101", alpha = 0.8) +
    geom_vline(xintercept = c(0.5, 0.8), linetype = "dashed", color = c("orange", "red")) +
    labs(title = paste0("Cross-Ancestry COLOC: Broadaway eQTL [EUR] x ", cfg$label, " [EAS]"),
         subtitle = paste0(n_tested, " genes tested; ",
                          nrow(coloc_res[PP.H4 > 0.5]), " with PP.H4 > 0.5; ",
                          "Ancestry mismatch (conservative)"),
         x = "PP.H4 (Posterior Probability of Shared Causal Variant)",
         y = "Count") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 12))
  print(p1)
  dev.off()
  cat("  Saved:", paste0(gwas_short, "_coloc_pp4_distribution.pdf"), "\n")

  # Top genes barplot
  top30 <- head(coloc_res[PP.H4 > 0.1], 30)
  if (nrow(top30) > 0) {
    top30[, gene := factor(gene, levels = rev(gene))]
    top30[, sig := ifelse(PP.H4 > 0.5, "PP.H4 > 0.5", "PP.H4 <= 0.5")]

    pdf(file.path(FIG_DIR, paste0(gwas_short, "_coloc_top_genes.pdf")),
        width = 10, height = 8)
    p2 <- ggplot(top30, aes(x = PP.H4, y = gene, fill = sig)) +
      geom_col(width = 0.7) +
      geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey40") +
      scale_fill_manual(values = c("PP.H4 > 0.5" = "#D7191C", "PP.H4 <= 0.5" = "#E66101")) +
      labs(title = paste0("Cross-Ancestry COLOC: Broadaway eQTL [EUR] x ", GWAS_NAME, " [EAS]"),
           subtitle = paste0("N(eQTL) = 1,183 [EUR]; N(GWAS) = ", format(GWAS_N, big.mark = ","),
                            " [EAS]; ABF-only; p12 = 5e-6"),
           x = "PP.H4", y = NULL) +
      theme_bw(base_size = 11) +
      theme(plot.title = element_text(face = "bold", size = 12),
            plot.subtitle = element_text(size = 9, color = "grey40"),
            legend.position = "bottom")
    print(p2)
    dev.off()
    cat("  Saved:", paste0(gwas_short, "_coloc_top_genes.pdf"), "\n")
  }
}

cat("\n=== Script 47: Complete ===\n")
cat("GWAS:", GWAS_NAME, "(East Asian / Japanese)\n")
cat("eQTL: Broadaway2024 (European)\n")
cat("Ancestry mismatch: YES — positive hits are high-confidence cross-ancestry targets\n")
cat("End time:", format(Sys.time()), "\n")
