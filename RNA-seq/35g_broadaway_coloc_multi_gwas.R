#!/usr/bin/env Rscript
# 35g_broadaway_coloc_multi_gwas.R
# ---------------------------------------------------------------------------
# Parameterized Broadaway Liver eQTL COLOC with UKBB quantitative GWAS
#
# Generalization of Script 35b (Broadaway x UKBB ALT) to support AST and GGT.
# All three GWAS are from the same study (Sinnott-Armstrong 2021, PMID 33462484),
# use identical column formats, and are hg38 harmonised.
#
# Usage: GWAS_NAME=UKBB_AST Rscript 35g_broadaway_coloc_multi_gwas.R
#        GWAS_NAME=UKBB_GGT Rscript 35g_broadaway_coloc_multi_gwas.R
#        GWAS_NAME=UKBB_ALT Rscript 35g_broadaway_coloc_multi_gwas.R  (equivalent to 35b)
#
# Inputs:
#   - data/broadaway_eqtl/Liver_eQTL_Meta_Leads_ST3_20240530.tsv
#   - data/broadaway_eqtl/chr{1-22}_marginal_summary_results.tsv
#   - GWAS/MR_Data/{GWAS_FILE} (from gwas_config)
#   - data/broadaway_eqtl/hg19ToHg38.over.chain
#
# Outputs:
#   - RNA-seq/results/causal_inference/broadaway_{gwas_name}/coloc_results.csv
#   - RNA-seq/results/causal_inference/broadaway_{gwas_name}/coloc_summary.csv
#   - RNA-seq/results/causal_inference/broadaway_{gwas_name}/coloc_sensitivity/
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(coloc)
  library(data.table)
  library(rtracklayer)
  library(GenomicRanges)
  library(ggplot2)
})

# ==============================================================================
# Configuration — parameterized by GWAS_NAME env var
# ==============================================================================
GWAS_NAME <- Sys.getenv("GWAS_NAME", "UKBB_AST")

gwas_config <- list(
  UKBB_ALT = list(
    file = "GCST90019492_UKBB_ALT_harmonised.tsv.gz",
    N    = 343850L,
    type = "quant",
    label = "UKBB ALT (Sinnott-Armstrong 2021)"
  ),
  UKBB_AST = list(
    file = "GCST90019497_UKBB_AST_harmonised.tsv.gz",
    N    = 343850L,
    type = "quant",
    label = "UKBB AST (Sinnott-Armstrong 2021)"
  ),
  UKBB_GGT = list(
    file = "GCST90019507_UKBB_GGT_harmonised.tsv.gz",
    N    = 343850L,
    type = "quant",
    label = "UKBB GGT (Sinnott-Armstrong 2021)"
  )
)

if (!GWAS_NAME %in% names(gwas_config)) {
  stop("Unknown GWAS_NAME: ", GWAS_NAME,
       ". Valid options: ", paste(names(gwas_config), collapse = ", "))
}

cfg <- gwas_config[[GWAS_NAME]]

BASE_DIR    <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
EQTL_DIR   <- file.path(BASE_DIR, "data/broadaway_eqtl")
GWAS_FILE   <- file.path(BASE_DIR, "GWAS/MR_Data", cfg$file)
RESULTS_DIR <- file.path(BASE_DIR, "RNA-seq/results/causal_inference",
                          paste0("broadaway_", tolower(GWAS_NAME)))
SENS_DIR    <- file.path(RESULTS_DIR, "coloc_sensitivity")

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SENS_DIR, recursive = TRUE, showWarnings = FALSE)

GWAS_N    <- cfg$N
GWAS_TYPE <- cfg$type

# Broadaway eQTL parameters
EQTL_N    <- 1183L

# Minimum overlapping SNPs for COLOC (raised from 10; see threshold audit)
MIN_SNPS <- 100L

# COLOC priors
COLOC_P1  <- 1e-4
COLOC_P2  <- 1e-4
COLOC_P12 <- 5e-6

# Liftover chain file
CHAIN_FILE <- file.path(EQTL_DIR, "hg19ToHg38.over.chain")
CHAIN_GZ   <- paste0(CHAIN_FILE, ".gz")

cat("=== Script 35g: Broadaway eQTL COLOC x", GWAS_NAME, "GWAS ===\n")
cat("GWAS:", cfg$label, "\n")
cat("GWAS file:", cfg$file, "\n")
cat("GWAS N:", GWAS_N, "\n")
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
  stop("GWAS file not found: ", GWAS_FILE,
       "\nRun GWAS/download_ukbb_liver_enzymes.sh first.")
}

gwas <- fread(GWAS_FILE)
cat("  Raw GWAS rows:", nrow(gwas), "\n")

# Filter to valid autosomal rows
gwas <- gwas[!is.na(chromosome) & !is.na(base_pair_location) & !is.na(p_value)]
gwas[, chr := as.integer(chromosome)]
gwas[, pos_hg38 := as.integer(base_pair_location)]
gwas <- gwas[!is.na(chr) & chr %in% 1:22]

# Need beta + SE for beta/varbeta COLOC
gwas <- gwas[!is.na(beta) & !is.na(standard_error) & standard_error > 0]

# Deduplicate by position (keep lowest p-value)
gwas <- gwas[order(chr, pos_hg38, p_value)]
gwas <- gwas[!duplicated(paste(chr, pos_hg38))]
gwas[, merge_key := paste0(chr, ":", pos_hg38)]

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

# Build gene → chr mapping from leads Variant field (format: CHR_POS_REF_ALT)
leads[, eqtl_chr := as.integer(sub("_.*", "", Variant))]
gene_chr_map <- unique(leads[!is.na(eqtl_chr), .(Gene, chr = eqtl_chr)])
gene_chr_map <- gene_chr_map[, .N, by = .(Gene, chr)][order(-N)][!duplicated(Gene)]
cat("  Gene-to-chromosome mapping:", nrow(gene_chr_map), "genes\n")

# Build gene symbol <-> Ensembl lookup from leads
leads_lookup <- unique(leads[, .(Gene, Ensembl)])

# ==============================================================================
# 4. Helper functions
# ==============================================================================

# Liftover hg19 → hg38 (vectorized)
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

# Per-chromosome eQTL cache
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

  # Create merge key
  dt[, merge_key := paste0(CHR, ":", pos_hg38)]

  chr_eqtl_cache[[key]] <<- dt
  return(dt)
}

# Allele harmonization
harmonize_alleles <- function(dt) {
  dt <- copy(dt)
  comp <- c("A" = "T", "T" = "A", "C" = "G", "G" = "C")

  dt[, match_type := "none"]

  # Broadaway: EA = effect allele, NEA = non-effect allele
  # UKBB: effect_allele, other_allele
  # Direct: EA == effect_allele AND NEA == other_allele
  dt[EA == effect_allele & NEA == other_allele, match_type := "direct"]
  # Flipped: EA == other_allele AND NEA == effect_allele
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

  # Flip eQTL beta when alleles are swapped
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

  if (nrow(gene_eqtl) < MIN_SNPS) {
    n_low_snps <- n_low_snps + 1L
    next
  }

  # Merge eQTL with GWAS on chr:pos (hg38)
  merged <- merge(gene_eqtl, gwas, by = "merge_key", suffixes = c(".eqtl", ".gwas"))

  if (nrow(merged) < MIN_SNPS) {
    n_low_snps <- n_low_snps + 1L
    next
  }

  # Filter to rows with valid beta/SE in both datasets
  merged <- merged[!is.na(Beta) & !is.na(SE) & SE > 0 &
                   !is.na(beta) & !is.na(standard_error) & standard_error > 0]

  if (nrow(merged) < MIN_SNPS) {
    n_low_snps <- n_low_snps + 1L
    next
  }

  # Allele harmonization
  merged <- harmonize_alleles(merged)

  if (nrow(merged) < MIN_SNPS) {
    n_no_harm <- n_no_harm + 1L
    next
  }

  # Prepare COLOC datasets with beta/varbeta
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

    # Extract top SNP
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
      gwas_source = GWAS_NAME,
      gwas_N      = GWAS_N,
      method      = "abf_beta_varbeta"
    )

    # Sensitivity analysis for promising hits
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
cat("\n  COLOC Summary (Broadaway N=1,183 x", gwas_label, "):\n")
cat("    Genes tested:", n_tested, "\n")
cat("    PP.H4 > 0.8:", nrow(coloc_res[PP.H4 > 0.8]), "\n")
cat("    PP.H4 > 0.5:", nrow(coloc_res[PP.H4 > 0.5]), "\n")
cat("    PP.H4 > 0.3:", nrow(coloc_res[PP.H4 > 0.3]), "\n")

# Top hits
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
             "eqtl_source", "eqtl_N", "gwas_source", "gwas_N"),
  value = c(n_tested, n_skipped, n_low_snps, n_no_harm,
            nrow(coloc_res[PP.H4 > 0.8]),
            nrow(coloc_res[PP.H4 > 0.5]),
            nrow(coloc_res[PP.H4 > 0.3]),
            "Broadaway2024", EQTL_N,
            GWAS_NAME, GWAS_N)
)
fwrite(summary_dt, file.path(RESULTS_DIR, "coloc_summary.csv"))

# ==============================================================================
# 7. Validation checks
# ==============================================================================
cat("\n--- Step 7: Validation ---\n")

# Known MASLD/liver COLOC genes (from Zenodo sc-eQTL and UKBB ALT run)
known_genes <- c("HSD17B13", "EFHD1", "CHEK2", "MAP1LC3A", "SPTLC3",
                 "APOH", "MFHAS1", "ETS2", "LCOR", "USP40", "ALAD",
                 "FABP1", "RORA", "HKDC1")
cat("\n  Known MASLD gene validation:\n")
for (g in known_genes) {
  hit <- coloc_res[gene == g]
  if (nrow(hit) > 0) {
    cat("    ", g, ": PP.H4 =", round(hit$PP.H4[1], 3),
        "(", hit$n_snps[1], "SNPs)\n")
  } else {
    cat("    ", g, ": not tested (no eQTL or insufficient overlap)\n")
  }
}

# Compare with UKBB ALT results (Script 35b)
alt_file <- file.path(BASE_DIR, "RNA-seq/results/causal_inference/broadaway_ukbb/coloc_results.csv")
if (file.exists(alt_file) && GWAS_NAME != "UKBB_ALT") {
  alt_res <- fread(alt_file)
  if (nrow(alt_res) > 0 && "PP.H4" %in% names(alt_res)) {
    cat("\n  Comparison with Broadaway x UKBB ALT:\n")
    cat("    ALT PP.H4 > 0.5:", nrow(alt_res[PP.H4 > 0.5]), "\n")
    cat("    ", GWAS_NAME, "PP.H4 > 0.5:", nrow(coloc_res[PP.H4 > 0.5]), "\n")
    shared <- intersect(alt_res$gene, coloc_res$gene)
    if (length(shared) > 5) {
      comp <- merge(
        alt_res[, .(gene, PP.H4_alt = PP.H4)],
        coloc_res[, .(gene, PP.H4_new = PP.H4)],
        by = "gene"
      )
      cat("    Shared genes tested:", nrow(comp), "\n")
      cat("    Pearson r (PP.H4):", round(cor(comp$PP.H4_alt, comp$PP.H4_new,
                                               use = "complete.obs"), 3), "\n")
      # Unique to this GWAS
      new_unique <- setdiff(coloc_res[PP.H4 > 0.5]$gene, alt_res[PP.H4 > 0.5]$gene)
      cat("    Unique to", GWAS_NAME, "(PP.H4>0.5):", length(new_unique), "\n")
      if (length(new_unique) > 0) {
        cat("      ", paste(head(new_unique, 20), collapse = ", "), "\n")
      }
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
  pdf(file.path(FIG_DIR, paste0("broadaway_", gwas_short, "_coloc_pp4_distribution.pdf")),
      width = 8, height = 5)
  p1 <- ggplot(coloc_res, aes(x = PP.H4)) +
    geom_histogram(bins = 50, fill = "#2C7BB6", alpha = 0.8) +
    geom_vline(xintercept = c(0.5, 0.8), linetype = "dashed", color = c("orange", "red")) +
    labs(title = paste0("COLOC PP.H4: Broadaway Liver eQTL (N=1,183) x ", cfg$label),
         subtitle = paste0(n_tested, " genes tested; ",
                          nrow(coloc_res[PP.H4 > 0.5]), " with PP.H4 > 0.5"),
         x = "PP.H4 (Posterior Probability of Shared Causal Variant)",
         y = "Count") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 13))
  print(p1)
  dev.off()
  cat("  Saved:", paste0("broadaway_", gwas_short, "_coloc_pp4_distribution.pdf"), "\n")

  # Top genes barplot
  top30 <- head(coloc_res[PP.H4 > 0.1], 30)
  if (nrow(top30) > 0) {
    top30[, gene := factor(gene, levels = rev(gene))]
    top30[, sig := ifelse(PP.H4 > 0.5, "PP.H4 > 0.5", "PP.H4 <= 0.5")]

    pdf(file.path(FIG_DIR, paste0("broadaway_", gwas_short, "_coloc_top_genes.pdf")),
        width = 10, height = 8)
    p2 <- ggplot(top30, aes(x = PP.H4, y = gene, fill = sig)) +
      geom_col(width = 0.7) +
      geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey40") +
      scale_fill_manual(values = c("PP.H4 > 0.5" = "#D7191C", "PP.H4 <= 0.5" = "#2C7BB6")) +
      labs(title = paste0("Top COLOC Hits: Broadaway Liver eQTL x ", GWAS_NAME, " GWAS"),
           subtitle = paste0("N(eQTL) = 1,183; N(GWAS) = ", format(GWAS_N, big.mark = ","),
                            "; priors: p12 = 5e-6; beta/varbeta"),
           x = "PP.H4", y = NULL) +
      theme_bw(base_size = 11) +
      theme(plot.title = element_text(face = "bold", size = 13),
            plot.subtitle = element_text(size = 9, color = "grey40"),
            legend.position = "bottom")
    print(p2)
    dev.off()
    cat("  Saved:", paste0("broadaway_", gwas_short, "_coloc_top_genes.pdf"), "\n")
  }
}

cat("\n=== Script 35g: Complete ===\n")
cat("GWAS:", GWAS_NAME, "\n")
cat("End time:", format(Sys.time()), "\n")
