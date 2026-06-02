#!/usr/bin/env Rscript
# 50b_strand_ambiguity_sensitivity.R
# ---------------------------------------------------------------------------
# Strand-Ambiguity Sensitivity Analysis for COLOC
#
# Quantifies the impact of filtering A/T and C/G SNPs before complement
# matching. Runs COLOC with and without the strand-ambiguous filter for
# a completed Pan-UKBB analysis to measure impact on PP.H4.
#
# Expects panukbb_afr_alt to be complete (first completed Pan-UKBB run).
#
# Output: RNA-seq/results/causal_inference/strand_sensitivity/
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(coloc)
  library(data.table)
  library(rtracklayer)
  library(GenomicRanges)
  library(ggplot2)
})

BASE_DIR    <- Sys.getenv("MASLD_PROJECT_ROOT",
                          "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
EQTL_DIR   <- file.path(BASE_DIR, "data/broadaway_eqtl")
GWAS_FILE   <- file.path(BASE_DIR, "GWAS/MR_Data/PanUKBB/PanUKBB_AFR_ALT_harmonised_hg38.tsv.gz")
RESULTS_DIR <- file.path(BASE_DIR, "RNA-seq/results/causal_inference/strand_sensitivity")
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

GWAS_N    <- 6636L
EQTL_N    <- 1183L
COLOC_P1  <- 1e-4
COLOC_P2  <- 1e-4
COLOC_P12 <- 5e-6
MIN_MAF   <- 0.01

cat("=== Script 50b: Strand-Ambiguity Sensitivity ===\n")
cat("Test GWAS: PanUKBB AFR ALT (N =", GWAS_N, ")\n")
cat("Start time:", format(Sys.time()), "\n\n")

# ==============================================================================
# 1. Chain file + GWAS
# ==============================================================================
CHAIN_FILE <- file.path(EQTL_DIR, "hg19ToHg38.over.chain")
CHAIN_GZ   <- paste0(CHAIN_FILE, ".gz")
if (!file.exists(CHAIN_FILE)) {
  if (!file.exists(CHAIN_GZ)) {
    download.file(
      "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz",
      CHAIN_GZ, mode = "wb", quiet = FALSE)
  }
  system2("gunzip", args = c("-k", CHAIN_GZ))
}
chain <- import.chain(CHAIN_FILE)

gwas <- fread(GWAS_FILE)
gwas[, chr := as.integer(chromosome)]
gwas[, pos_hg38 := as.integer(base_pair_location)]
gwas <- gwas[!is.na(chr) & chr %in% 1:22]
gwas <- gwas[!is.na(beta) & !is.na(standard_error) & standard_error > 0]
if ("EAF" %in% names(gwas)) {
  gwas[, maf := pmin(EAF, 1 - EAF)]
  gwas <- gwas[is.na(maf) | maf >= MIN_MAF]
}
gwas <- gwas[order(chr, pos_hg38, p_value)]
gwas <- gwas[!duplicated(paste(chr, pos_hg38))]
gwas[, merge_key := paste0(chr, ":", pos_hg38)]
cat("  GWAS variants:", format(nrow(gwas), big.mark = ","), "\n")

# ==============================================================================
# 2. Leads + helpers
# ==============================================================================
leads <- fread(file.path(EQTL_DIR, "Liver_eQTL_Meta_Leads_ST3_20240530.tsv"))
all_eGenes <- unique(leads$Gene)
leads[, eqtl_chr := as.integer(sub("_.*", "", Variant))]
gene_chr_map <- unique(leads[!is.na(eqtl_chr), .(Gene, chr = eqtl_chr)])
gene_chr_map <- gene_chr_map[, .N, by = .(Gene, chr)][order(-N)][!duplicated(Gene)]
leads_lookup <- unique(leads[, .(Gene, Ensembl)])

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
  dt <- fread(fname)
  dt[, pos_hg38 := liftover_positions(chr_num, POS)]
  dt <- dt[!is.na(pos_hg38)]
  dt[, merge_key := paste0(CHR, ":", pos_hg38)]
  chr_eqtl_cache[[key]] <<- dt
  return(dt)
}

# ==============================================================================
# 3. Two harmonization modes
# ==============================================================================
harmonize_with_filter <- function(dt) {
  dt <- copy(dt)
  comp <- c("A" = "T", "T" = "A", "C" = "G", "G" = "C")
  dt[, match_type := "none"]
  dt[EA == effect_allele & NEA == other_allele, match_type := "direct"]
  dt[EA == other_allele & NEA == effect_allele, match_type := "flipped"]
  # Filter strand-ambiguous before complement
  dt[, is_ambiguous := (EA %in% c("A","T") & NEA %in% c("A","T")) |
                        (EA %in% c("C","G") & NEA %in% c("C","G"))]
  dt <- dt[!(is_ambiguous & match_type == "none")]
  dt[match_type == "none" & nchar(EA) == 1 & nchar(NEA) == 1,
     match_type := fifelse(
       comp[EA] == effect_allele & comp[NEA] == other_allele, "direct",
       fifelse(comp[EA] == other_allele & comp[NEA] == effect_allele, "flipped", "none"))]
  dt <- dt[match_type != "none"]
  if (any(dt$match_type == "flipped")) dt[match_type == "flipped", Beta := -Beta]
  dt[, is_ambiguous := NULL]
  return(dt)
}

harmonize_without_filter <- function(dt) {
  dt <- copy(dt)
  comp <- c("A" = "T", "T" = "A", "C" = "G", "G" = "C")
  dt[, match_type := "none"]
  dt[EA == effect_allele & NEA == other_allele, match_type := "direct"]
  dt[EA == other_allele & NEA == effect_allele, match_type := "flipped"]
  # No filter — complement all unmatched
  dt[match_type == "none" & nchar(EA) == 1 & nchar(NEA) == 1,
     match_type := fifelse(
       comp[EA] == effect_allele & comp[NEA] == other_allele, "direct",
       fifelse(comp[EA] == other_allele & comp[NEA] == effect_allele, "flipped", "none"))]
  dt <- dt[match_type != "none"]
  if (any(dt$match_type == "flipped")) dt[match_type == "flipped", Beta := -Beta]
  return(dt)
}

run_coloc <- function(merged) {
  dataset1 <- list(snp = merged$merge_key, beta = merged$beta,
                   varbeta = (merged$standard_error)^2,
                   type = "quant", sdY = 1, N = GWAS_N)
  dataset2 <- list(snp = merged$merge_key, beta = merged$Beta,
                   varbeta = (merged$SE)^2,
                   type = "quant", sdY = 1, N = EQTL_N)
  tryCatch(
    coloc.abf(dataset1, dataset2, p1 = COLOC_P1, p2 = COLOC_P2, p12 = COLOC_P12),
    error = function(e) NULL
  )
}

# ==============================================================================
# 4. Pre-filter: only test genes with non-trivial PP.H4 in existing results
# ==============================================================================
EXISTING_RESULTS <- file.path(BASE_DIR, "RNA-seq/results/causal_inference/panukbb_afr_alt/coloc_results.csv")
PP4_THRESHOLD <- 0.05  # only test genes where filtering could plausibly change status

if (file.exists(EXISTING_RESULTS)) {
  existing <- fread(EXISTING_RESULTS)
  candidate_genes <- existing[PP.H4 >= PP4_THRESHOLD, unique(gene)]
  cat("  Pre-filter: ", length(candidate_genes), "genes with PP.H4 >=", PP4_THRESHOLD,
      "(from", nrow(existing), "total)\n")
  # Also include all eGenes with PP.H4 > 0.3 for thorough analysis
  genes_to_test <- intersect(all_eGenes, candidate_genes)
} else {
  cat("  No existing results found — testing all eGenes (slow)\n")
  genes_to_test <- all_eGenes
}
cat("  Genes to test:", length(genes_to_test), "\n\n")

# ==============================================================================
# 4b. Run COLOC both ways (pre-filtered)
# ==============================================================================
cat("--- Running paired COLOC (with/without strand filter) ---\n")

results <- list()
n_tested <- 0L
total_ambig_removed <- 0L

for (gene_name in genes_to_test) {
  chr_info <- gene_chr_map[Gene == gene_name]
  if (nrow(chr_info) == 0) next
  gene_chr <- chr_info$chr[1]

  eqtl_data <- load_chr_eqtl(gene_chr)
  if (is.null(eqtl_data)) next

  gene_eqtl <- eqtl_data[GeneSymbol == gene_name]
  if (nrow(gene_eqtl) == 0) {
    gene_ensg <- leads_lookup$Ensembl[leads_lookup$Gene == gene_name]
    if (length(gene_ensg) > 0) gene_eqtl <- eqtl_data[ENSG %in% gene_ensg]
  }
  if (nrow(gene_eqtl) < 10) next

  merged_raw <- merge(gene_eqtl, gwas, by = "merge_key", suffixes = c(".eqtl", ".gwas"))
  merged_raw <- merged_raw[!is.na(Beta) & !is.na(SE) & SE > 0 &
                            !is.na(beta) & !is.na(standard_error) & standard_error > 0]
  if (nrow(merged_raw) < 10) next

  # With filter
  m_filt <- harmonize_with_filter(merged_raw)
  # Without filter (original behavior)
  m_nofilt <- harmonize_without_filter(merged_raw)

  if (nrow(m_filt) < 10 && nrow(m_nofilt) < 10) next

  n_removed <- nrow(m_nofilt) - nrow(m_filt)
  total_ambig_removed <- total_ambig_removed + max(0, n_removed)

  res_filt <- if (nrow(m_filt) >= 10) run_coloc(m_filt) else NULL
  res_nofilt <- if (nrow(m_nofilt) >= 10) run_coloc(m_nofilt) else NULL

  if (!is.null(res_filt) || !is.null(res_nofilt)) {
    n_tested <- n_tested + 1L
    results[[gene_name]] <- data.table(
      gene = gene_name,
      pp4_filtered = if (!is.null(res_filt)) res_filt$summary["PP.H4.abf"] else NA_real_,
      pp4_unfiltered = if (!is.null(res_nofilt)) res_nofilt$summary["PP.H4.abf"] else NA_real_,
      n_snps_filtered = nrow(m_filt),
      n_snps_unfiltered = nrow(m_nofilt),
      n_ambig_removed = max(0, n_removed)
    )
  }

  if (n_tested %% 200 == 0) cat("  [", n_tested, "tested]\n")
}

# ==============================================================================
# 5. Analyze impact
# ==============================================================================
cat("\n--- Results ---\n")
cat("  Genes tested (paired):", n_tested, "\n")
cat("  Total ambiguous SNPs removed:", total_ambig_removed, "\n")

if (length(results) == 0) {
  cat("  No results to compare.\n")
  quit(save = "no", status = 0)
}

comp <- rbindlist(results, fill = TRUE)
comp <- comp[!is.na(pp4_filtered) & !is.na(pp4_unfiltered)]

cat("  Genes with both results:", nrow(comp), "\n")

# Correlation
if (nrow(comp) > 5) {
  rho <- cor(comp$pp4_filtered, comp$pp4_unfiltered, method = "spearman")
  r   <- cor(comp$pp4_filtered, comp$pp4_unfiltered, method = "pearson")
  cat("  Spearman rho:", round(rho, 4), "\n")
  cat("  Pearson r:", round(r, 4), "\n")
}

# Status changes
sig_filt   <- comp$pp4_filtered > 0.5
sig_nofilt <- comp$pp4_unfiltered > 0.5

gained  <- sum(!sig_nofilt & sig_filt)   # new sig after filtering
lost    <- sum(sig_nofilt & !sig_filt)    # lost sig after filtering
stable  <- sum(sig_filt == sig_nofilt)

cat("  PP.H4 > 0.5 (filtered):", sum(sig_filt), "\n")
cat("  PP.H4 > 0.5 (unfiltered):", sum(sig_nofilt), "\n")
cat("  Status changed:", gained + lost, "/", nrow(comp),
    "(", round(100 * (gained + lost) / nrow(comp), 1), "%)\n")
cat("    Gained significance:", gained, "\n")
cat("    Lost significance:", lost, "\n")
cat("    Stable:", stable, "\n")

# Genes with largest PP.H4 shift
comp[, pp4_diff := pp4_filtered - pp4_unfiltered]
comp[, abs_diff := abs(pp4_diff)]
setorder(comp, -abs_diff)

cat("\n  Top 10 genes by |PP.H4 change|:\n")
for (i in seq_len(min(10, nrow(comp)))) {
  cat("    ", comp$gene[i], ": filtered=", round(comp$pp4_filtered[i], 3),
      " unfiltered=", round(comp$pp4_unfiltered[i], 3),
      " (diff=", round(comp$pp4_diff[i], 3),
      ", removed=", comp$n_ambig_removed[i], ")\n")
}

# ==============================================================================
# 6. Save + figure
# ==============================================================================
fwrite(comp, file.path(RESULTS_DIR, "strand_sensitivity_comparison.csv"))
cat("\n  Saved: strand_sensitivity_comparison.csv\n")

# Summary
summary_dt <- data.table(
  metric = c("genes_paired", "spearman_rho", "pearson_r",
             "n_sig_filtered", "n_sig_unfiltered",
             "status_changed", "pct_changed",
             "total_ambig_removed"),
  value = c(nrow(comp),
            round(cor(comp$pp4_filtered, comp$pp4_unfiltered, method = "spearman"), 4),
            round(cor(comp$pp4_filtered, comp$pp4_unfiltered, method = "pearson"), 4),
            sum(sig_filt), sum(sig_nofilt),
            gained + lost,
            round(100 * (gained + lost) / nrow(comp), 1),
            total_ambig_removed)
)
fwrite(summary_dt, file.path(RESULTS_DIR, "strand_sensitivity_summary.csv"))

# Figure: scatter
if (nrow(comp) > 5) {
  pdf(file.path(RESULTS_DIR, "strand_sensitivity_scatter.pdf"), width = 7, height = 7)
  p <- ggplot(comp, aes(x = pp4_unfiltered, y = pp4_filtered)) +
    geom_point(alpha = 0.3, size = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = 0.5, linetype = "dotted", color = "red", alpha = 0.5) +
    geom_vline(xintercept = 0.5, linetype = "dotted", color = "red", alpha = 0.5) +
    labs(title = "Strand-Ambiguity Sensitivity: COLOC PP.H4",
         subtitle = paste0("PanUKBB AFR ALT; rho=",
                          round(cor(comp$pp4_filtered, comp$pp4_unfiltered, method = "spearman"), 3),
                          "; status changed: ", gained + lost, "/", nrow(comp)),
         x = "PP.H4 (without strand filter)", y = "PP.H4 (with strand filter)") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  print(p)
  dev.off()
  cat("  Saved: strand_sensitivity_scatter.pdf\n")
}

cat("\n=== Script 50b: Complete ===\n")
cat("VERDICT:", ifelse((gained + lost) / nrow(comp) < 0.05,
                       "NEGLIGIBLE (<5% status change) — no re-run needed",
                       "SIGNIFICANT (>=5%) — consider re-running affected COLOC"), "\n")
cat("End time:", format(Sys.time()), "\n")
