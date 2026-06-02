#!/usr/bin/env Rscript
# 35_broadaway_coloc.R
# ---------------------------------------------------------------------------
# Colocalization with Broadaway et al. 2024 Liver eQTL (N=1,183)
#
# Replaces underpowered GTEx liver eQTL (N=208) COLOC from Script 19.
# Broadaway eQTL data: 6,564 eGenes, 9,013 signals, 22 per-chr marginal files.
# Genome build: eQTL = hg19, GWAS (Ghodsian) = hg38 → liftover required.
#
# Steps:
#   1. Load GWAS summary statistics (hg38)
#   2. Load dream DEGs to select candidate genes
#   3. Load Broadaway leads file to identify testable eGenes
#   4. For each candidate gene:
#      a. Load per-chr marginal eQTL stats
#      b. Liftover eQTL positions hg19 → hg38
#      c. Merge with GWAS on chr:pos
#      d. Run coloc.abf()
#   5. Save results and sensitivity plots
#
# Inputs:
#   - data/broadaway_eqtl/Liver_eQTL_Meta_Leads_ST3_20240530.tsv
#   - data/broadaway_eqtl/chr{1-22}_marginal_summary_results.tsv
#   - GWAS/MR_Data/Ghodsian_2021_NAFLD_harmonised.tsv.gz
#   - RNA-seq/Human/.../dream_results.csv
#   - data/broadaway_eqtl/hg19ToHg38.over.chain (downloaded if missing)
#
# Outputs:
#   - results/causal_inference/broadaway/coloc_results.csv
#   - results/causal_inference/broadaway/coloc_sensitivity/ (PDFs)
#   - results/causal_inference/broadaway/coloc_summary.csv
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(coloc)
  library(data.table)
  library(rtracklayer)
  library(GenomicRanges)
  library(ggplot2)
})

cat("=== Script 35: Broadaway eQTL Colocalization ===\n")
cat("Start time:", format(Sys.time()), "\n\n")

# ==============================================================================
# Configuration
# ==============================================================================
BASE_DIR    <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
EQTL_DIR   <- file.path(BASE_DIR, "data/broadaway_eqtl")
GWAS_DIR    <- file.path(BASE_DIR, "GWAS/MR_Data")
RESULTS_DIR <- file.path(BASE_DIR, "RNA-seq/results/causal_inference/broadaway")
SENS_DIR    <- file.path(RESULTS_DIR, "coloc_sensitivity")
DREAM_FILE  <- file.path(BASE_DIR,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/dream_results.csv")
GENE_CACHE  <- file.path(BASE_DIR,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/gene_annotation/human_ensg_to_symbol.tsv")

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SENS_DIR, recursive = TRUE, showWarnings = FALSE)

# GWAS parameters (Ghodsian 2021: EHR-based NAFLD case-control, N=778,614)
GWAS_FILE <- file.path(GWAS_DIR, "Ghodsian_2021_NAFLD_harmonised.tsv.gz")
GWAS_N    <- 778614L
GWAS_TYPE <- "cc"     # Ghodsian is case-control (EHR NAFLD diagnosis)
GWAS_S    <- 0.0108   # Case fraction = 8434/778614 (review B2#4; was 0.08 ~7x too high — Script 36 is pseudobulk DE, not a case-fraction source)

# Broadaway eQTL parameters
EQTL_N    <- 1183L  # meta-analysis across 4 studies

# COLOC priors (same as Script 19)
COLOC_P1  <- 1e-4
COLOC_P2  <- 1e-4
COLOC_P12 <- 5e-6

# Minimum overlapping SNPs for COLOC (raised from 10; see threshold audit)
MIN_SNPS <- 100L

# Candidate selection threshold
DREAM_PADJ_THR <- 0.1

# Liftover chain file
CHAIN_FILE <- file.path(EQTL_DIR, "hg19ToHg38.over.chain")
CHAIN_GZ   <- paste0(CHAIN_FILE, ".gz")

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
# 2. Load GWAS summary statistics
# ==============================================================================
cat("\n--- Step 2: Loading GWAS (Ghodsian 2021) ---\n")

gwas <- fread(GWAS_FILE)
cat("  GWAS rows:", nrow(gwas), "\n")

# Standardize column names
# hm_chrom, hm_pos (hg38), hm_beta, p_value, standard_error
gwas <- gwas[!is.na(hm_chrom) & !is.na(hm_pos) & !is.na(p_value)]
gwas[, chr := as.integer(hm_chrom)]
gwas[, pos_hg38 := as.integer(hm_pos)]

# Deduplicate by chr:pos (keep lowest p-value)
gwas <- gwas[order(chr, pos_hg38, p_value)]
gwas <- gwas[!duplicated(paste(chr, pos_hg38))]
cat("  GWAS after dedup:", nrow(gwas), "variants\n")

# Create chr:pos key for merging
gwas[, merge_key := paste0(chr, ":", pos_hg38)]

# ==============================================================================
# 3. Load dream DEGs for candidate filtering
# ==============================================================================
cat("\n--- Step 3: Loading dream DEGs ---\n")

dream <- fread(DREAM_FILE)
dream[, ensembl_id := sub("\\.\\d+$", "", gene)]

# Map Ensembl -> symbol
if (file.exists(GENE_CACHE)) {
  ann <- fread(GENE_CACHE)
  ann_map <- ann[symbol != "" & !is.na(symbol) & !duplicated(gene_base),
                 .(ensembl_id = gene_base, symbol)]
  dream <- merge(dream, ann_map, by = "ensembl_id", all.x = TRUE)
} else {
  cat("  WARNING: Gene cache not found. Using Ensembl IDs.\n")
  dream[, symbol := ensembl_id]
}

deg_genes <- dream[padj < DREAM_PADJ_THR & !is.na(symbol), symbol]
deg_ensembl <- dream[padj < DREAM_PADJ_THR & !is.na(ensembl_id), ensembl_id]
cat("  Dream DEGs (padj <", DREAM_PADJ_THR, "):", length(deg_genes), "symbols,",
    length(deg_ensembl), "Ensembl IDs\n")

# ==============================================================================
# 4. Load Broadaway leads to identify testable eGenes
# ==============================================================================
cat("\n--- Step 4: Loading Broadaway eGene leads ---\n")

leads <- fread(file.path(EQTL_DIR, "Liver_eQTL_Meta_Leads_ST3_20240530.tsv"))
cat("  Total eQTL signals:", nrow(leads), "\n")
cat("  Unique eGenes:", length(unique(leads$Gene)), "\n")

# Intersect eGenes with DEGs
# Match on both symbol and Ensembl ID
eGenes_symbol  <- unique(leads$Gene)
eGenes_ensembl <- unique(leads$Ensembl)

deg_eGenes_sym  <- intersect(eGenes_symbol, deg_genes)
deg_eGenes_ensg <- intersect(eGenes_ensembl, deg_ensembl)

# Combine — some genes match on symbol, some on Ensembl
# Build a lookup from Ensembl → symbol from leads
leads_lookup <- unique(leads[, .(Gene, Ensembl)])
candidates_all <- unique(c(deg_eGenes_sym,
                           leads_lookup$Gene[leads_lookup$Ensembl %in% deg_eGenes_ensg]))
cat("  DEG ∩ eGene candidates:", length(candidates_all), "\n")

if (length(candidates_all) == 0) {
  cat("  No candidates to test. Exiting.\n")
  fwrite(data.table(gene = character(), PP.H4 = numeric()),
         file.path(RESULTS_DIR, "coloc_results.csv"))
  quit(save = "no", status = 0)
}

# Build gene → chr mapping from leads
gene_chr <- unique(leads[, .(Gene, chr = as.integer(sub("_.*", "",
  sub("^([0-9]+)_.*", "\\1", Variant))))])
# Actually, leads have no explicit CHR column. Extract from Variant format: "CHR_POS_REF_ALT"
# But Variant format varies. Use the Ensembl ID to get chr from marginal files instead.
# For efficiency, build gene-to-chromosome map from leads Variant field
leads[, eqtl_chr := as.integer(sub("_.*", "", Variant))]
gene_chr_map <- unique(leads[!is.na(eqtl_chr), .(Gene, chr = eqtl_chr)])
# Some genes span multiple chromosomes (trans-eQTL); take the most common
gene_chr_map <- gene_chr_map[, .N, by = .(Gene, chr)][order(-N)][!duplicated(Gene)]

cat("  Gene-to-chromosome mapping:", nrow(gene_chr_map), "genes\n")

# ==============================================================================
# 5. Run COLOC per candidate gene
# ==============================================================================
cat("\n--- Step 5: Running COLOC ---\n")

# Liftover function: convert hg19 positions to hg38 (vectorized)
liftover_positions <- function(chr_num, positions) {
  # Create GRanges in hg19
  gr <- GRanges(
    seqnames = paste0("chr", chr_num),
    ranges = IRanges(start = positions, width = 1)
  )
  # Liftover
  lifted <- liftOver(gr, chain)

  # Vectorized extraction (replaces slow for-loop)
  n_mapped <- lengths(lifted)
  hg38_pos <- rep(NA_integer_, length(positions))
  unique_idx <- which(n_mapped == 1L)
  if (length(unique_idx) > 0) {
    hg38_pos[unique_idx] <- start(unlist(lifted[unique_idx]))
  }
  hg38_pos
}

# Cache for per-chromosome eQTL data (avoid reloading)
chr_eqtl_cache <- list()

load_chr_eqtl <- function(chr_num) {
  key <- as.character(chr_num)
  if (!is.null(chr_eqtl_cache[[key]])) return(chr_eqtl_cache[[key]])

  fname <- file.path(EQTL_DIR, paste0("chr", chr_num, "_marginal_summary_results.tsv"))
  if (!file.exists(fname)) return(NULL)

  cat("    Loading chr", chr_num, "eQTL data...\n")
  dt <- fread(fname)
  # Columns: Entrez | Variant | CHR | POS(hg19) | NEA | EA | EAF | Beta | SE | PVAL | N | Studies | GeneSymbol | ENSG | Gene_Biotype

  # Liftover all positions for this chromosome at once (much more efficient)
  cat("    Lifting over", nrow(dt), "positions from hg19 to hg38...\n")
  dt[, pos_hg38 := liftover_positions(chr_num, POS)]
  dt <- dt[!is.na(pos_hg38)]
  cat("    Successful liftover:", nrow(dt), "variants\n")

  # Create merge key
  dt[, merge_key := paste0(CHR, ":", pos_hg38)]

  chr_eqtl_cache[[key]] <<- dt
  return(dt)
}

# Only process chromosomes that have candidate genes
candidate_chrs <- unique(gene_chr_map[Gene %in% candidates_all]$chr)
cat("  Chromosomes to process:", paste(sort(candidate_chrs), collapse = ", "), "\n")
cat("  (Skipping chromosomes without candidate genes)\n\n")

# Process candidates grouped by chromosome
results_list <- list()
n_tested  <- 0L
n_skipped <- 0L
n_low_snps <- 0L

for (gene_name in candidates_all) {
  # Get chromosome for this gene
  chr_info <- gene_chr_map[Gene == gene_name]
  if (nrow(chr_info) == 0) {
    n_skipped <- n_skipped + 1L
    next
  }
  gene_chr <- chr_info$chr[1]

  # Load chromosome eQTL data
  eqtl_data <- load_chr_eqtl(gene_chr)
  if (is.null(eqtl_data)) {
    n_skipped <- n_skipped + 1L
    next
  }

  # Get eQTL data for this gene
  gene_eqtl <- eqtl_data[GeneSymbol == gene_name]
  if (nrow(gene_eqtl) == 0) {
    # Try matching on ENSG
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

  # MAF: use eQTL EAF, clamp to (0.001, 0.499)
  use_maf <- pmin(pmax(merged$EAF, 0.001), 0.499)

  # Prepare coloc datasets
  dataset1 <- list(
    snp     = merged$merge_key,
    pvalues = merged$p_value,
    type    = GWAS_TYPE,
    N       = GWAS_N,
    MAF     = use_maf
  )
  # s (case fraction) required for case-control GWAS
  if (GWAS_TYPE == "cc") dataset1$s <- GWAS_S

  dataset2 <- list(
    snp     = merged$merge_key,
    pvalues = merged$PVAL,
    type    = "quant",
    N       = EQTL_N,
    MAF     = use_maf
  )

  # Run coloc.abf
  my.res <- tryCatch(
    coloc.abf(dataset1, dataset2, p1 = COLOC_P1, p2 = COLOC_P2, p12 = COLOC_P12),
    error = function(e) {
      cat("  coloc.abf error:", gene_name, "-", conditionMessage(e), "\n")
      NULL
    }
  )

  if (!is.null(my.res)) {
    n_tested <- n_tested + 1L

    results_list[[gene_name]] <- data.table(
      gene       = gene_name,
      ensembl    = paste(unique(gene_eqtl$ENSG[gene_eqtl$ENSG != ""]), collapse = ";"),
      chr        = gene_chr,
      PP.H0      = my.res$summary["PP.H0.abf"],
      PP.H1      = my.res$summary["PP.H1.abf"],
      PP.H2      = my.res$summary["PP.H2.abf"],
      PP.H3      = my.res$summary["PP.H3.abf"],
      PP.H4      = my.res$summary["PP.H4.abf"],
      n_snps     = nrow(merged),
      n_eqtl_snps = nrow(gene_eqtl),
      eqtl_source = "Broadaway2024",
      eqtl_N     = EQTL_N,
      method     = "abf"
    )

    # Sensitivity analysis for promising hits
    if (my.res$summary["PP.H4.abf"] > 0.3) {
      tryCatch({
        pdf(file.path(SENS_DIR, paste0("sensitivity_", gene_name, ".pdf")),
            width = 10, height = 8)
        sensitivity(my.res, rule = "H4 > 0.5")
        dev.off()
        cat("  Sensitivity plot:", gene_name,
            "(PP.H4 =", round(my.res$summary["PP.H4.abf"], 3), ")\n")
      }, error = function(e) NULL)
    }

    if (n_tested %% 50 == 0) {
      cat("  [", n_tested, "/", length(candidates_all), "tested,",
          n_low_snps, "low-SNP,", n_skipped, "skipped]\n")
    }
  }
}

# ==============================================================================
# 6. Compile and save results
# ==============================================================================
cat("\n--- Step 6: Compiling results ---\n")

if (length(results_list) == 0) {
  cat("  No COLOC results. Tested:", n_tested, "Skipped:", n_skipped,
      "Low SNPs:", n_low_snps, "\n")
  fwrite(data.table(gene = character(), PP.H4 = numeric()),
         file.path(RESULTS_DIR, "coloc_results.csv"))
  quit(save = "no", status = 0)
}

coloc_res <- rbindlist(results_list, fill = TRUE)
setorder(coloc_res, -PP.H4)

# Save full results
fwrite(coloc_res, file.path(RESULTS_DIR, "coloc_results.csv"))

# Summary statistics
cat("\n  COLOC Summary (Broadaway N=1,183 vs Ghodsian N=778,614):\n")
cat("    Genes tested:", n_tested, "\n")
cat("    Genes skipped (no chr mapping):", n_skipped, "\n")
cat("    Genes skipped (<", MIN_SNPS, "shared SNPs):", n_low_snps, "\n")
cat("    PP.H4 > 0.8:", nrow(coloc_res[PP.H4 > 0.8]), "\n")
cat("    PP.H4 > 0.5:", nrow(coloc_res[PP.H4 > 0.5]), "\n")
cat("    PP.H4 > 0.3:", nrow(coloc_res[PP.H4 > 0.3]), "\n")

# Top hits
if (nrow(coloc_res[PP.H4 > 0.3]) > 0) {
  cat("\n  Top COLOC hits (PP.H4 > 0.3):\n")
  top_hits <- coloc_res[PP.H4 > 0.3]
  for (i in seq_len(min(nrow(top_hits), 20))) {
    cat("    ", top_hits$gene[i], ": PP.H4 =", round(top_hits$PP.H4[i], 3),
        "(", top_hits$n_snps[i], "SNPs)\n")
  }
}

# Save summary
summary_dt <- data.table(
  metric = c("genes_tested", "genes_skipped", "genes_low_snps",
             "PP.H4_gt_0.8", "PP.H4_gt_0.5", "PP.H4_gt_0.3",
             "eqtl_source", "eqtl_N", "gwas_source", "gwas_N"),
  value = c(n_tested, n_skipped, n_low_snps,
            nrow(coloc_res[PP.H4 > 0.8]),
            nrow(coloc_res[PP.H4 > 0.5]),
            nrow(coloc_res[PP.H4 > 0.3]),
            "Broadaway2024", EQTL_N,
            "Ghodsian2021", GWAS_N)
)
fwrite(summary_dt, file.path(RESULTS_DIR, "coloc_summary.csv"))

# Compare with GTEx COLOC if available
gtex_coloc_file <- file.path(BASE_DIR, "RNA-seq/results/causal_inference/ghodsian/coloc_results.csv")
if (file.exists(gtex_coloc_file)) {
  gtex_coloc <- fread(gtex_coloc_file)
  if (nrow(gtex_coloc) > 0 && "PP.H4" %in% names(gtex_coloc)) {
    cat("\n  Comparison with GTEx (N=208) COLOC:\n")
    cat("    GTEx PP.H4 > 0.5:", nrow(gtex_coloc[PP.H4 > 0.5]), "\n")
    cat("    Broadaway PP.H4 > 0.5:", nrow(coloc_res[PP.H4 > 0.5]), "\n")
    shared_genes <- intersect(gtex_coloc$gene, coloc_res$gene)
    if (length(shared_genes) > 0) {
      comp <- merge(
        gtex_coloc[, .(gene, PP.H4_gtex = PP.H4)],
        coloc_res[, .(gene, PP.H4_broadaway = PP.H4)],
        by = "gene"
      )
      cat("    Shared genes:", nrow(comp), "\n")
      cat("    Pearson r (PP.H4):", round(cor(comp$PP.H4_gtex, comp$PP.H4_broadaway,
                                               use = "complete.obs"), 3), "\n")
    }
  }
}

# ==============================================================================
# 7. Visualization
# ==============================================================================
cat("\n--- Step 7: Generating figures ---\n")

FIG_DIR <- file.path(BASE_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

if (nrow(coloc_res) > 0) {
  # PP.H4 distribution
  pdf(file.path(FIG_DIR, "broadaway_coloc_pp4_distribution.pdf"), width = 8, height = 5)
  p1 <- ggplot(coloc_res, aes(x = PP.H4)) +
    geom_histogram(bins = 50, fill = "#2C7BB6", alpha = 0.8) +
    geom_vline(xintercept = c(0.5, 0.8), linetype = "dashed", color = c("orange", "red")) +
    labs(title = "COLOC PP.H4 Distribution (Broadaway Liver eQTL, N=1,183)",
         subtitle = paste0(n_tested, " genes tested; ",
                          nrow(coloc_res[PP.H4 > 0.5]), " with PP.H4 > 0.5"),
         x = "PP.H4 (Posterior Probability of Shared Causal Variant)",
         y = "Count") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 13))
  print(p1)
  dev.off()

  # Top genes barplot
  top30 <- head(coloc_res[PP.H4 > 0.1], 30)
  if (nrow(top30) > 0) {
    top30[, gene := factor(gene, levels = rev(gene))]
    top30[, sig := ifelse(PP.H4 > 0.5, "PP.H4 > 0.5", "PP.H4 <= 0.5")]

    pdf(file.path(FIG_DIR, "broadaway_coloc_top_genes.pdf"), width = 10, height = 8)
    p2 <- ggplot(top30, aes(x = PP.H4, y = gene, fill = sig)) +
      geom_col(width = 0.7) +
      geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey40") +
      scale_fill_manual(values = c("PP.H4 > 0.5" = "#D7191C", "PP.H4 <= 0.5" = "#2C7BB6")) +
      labs(title = "Top COLOC Hits: Broadaway Liver eQTL x Ghodsian MASLD GWAS",
           subtitle = paste0("N(eQTL) = 1,183; N(GWAS) = 778,614; priors: p12 = 5e-6"),
           x = "PP.H4", y = NULL) +
      theme_bw(base_size = 11) +
      theme(plot.title = element_text(face = "bold", size = 13),
            plot.subtitle = element_text(size = 9, color = "grey40"),
            legend.position = "bottom")
    print(p2)
    dev.off()
    cat("  Saved: broadaway_coloc_top_genes.pdf\n")
  }
}

cat("\n=== Script 35: Broadaway COLOC Complete ===\n")
cat("End time:", format(Sys.time()), "\n")
