#!/usr/bin/env Rscript
# 06b_coloc_window_sensitivity.R
# Window size sensitivity analysis for SuSiE-COLOC ABF
#
# Compares PP.H4 at three window sizes: ±250kb, ±500kb, ±1000kb (current)
# Window is centered on the top_snp (highest SNP.PP.H4) from the original run.
# This follows field-standard practice (OpenTargets, eQTL Catalogue, GTEx COLOC papers).
#
# Input:
#   GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv
#   GWAS/finemapping/config/gwas_registry.tsv
#   data/broadaway_eqtl/chr{N}_marginal_summary_results.tsv
#   GWAS/finemapping/data/sumstats/{study}_preprocessed.tsv
#
# Output:
#   GWAS/finemapping/results/coloc_threshold_audit/window_sensitivity.csv
#
# Run via: sbatch run_window_sensitivity.sbatch
# Runtime: ~2-4h on cpu partition (8 CPUs)

suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
})

setDTthreads(8)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE     <- Sys.getenv("MASLD_PROJECT_ROOT",
                       "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
FM_DIR   <- file.path(BASE, "GWAS/finemapping")
EQTL_DIR <- file.path(BASE, "data/broadaway_eqtl")

COLOC_FILE <- file.path(FM_DIR, "results/susie_coloc/susie_coloc_all_gwas.csv")
REGISTRY   <- file.path(FM_DIR, "config/gwas_registry.tsv")
OUT_DIR    <- file.path(FM_DIR, "results/coloc_threshold_audit")
OUT_FILE   <- file.path(OUT_DIR, "window_sensitivity.csv")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
# Windows to test (half-window in bp, centered on top_snp)
WINDOWS_BP <- c(250e3, 500e3, 1000e3)   # 250kb, 500kb, 1Mb
EQTL_N     <- 1183L
MIN_SNPS   <- 50L    # relaxed for narrower windows; standard filter is 100 for 1Mb
COLOC_P1   <- 1e-4
COLOC_P2   <- 1e-4
COLOC_P12  <- 5e-6   # match 06_susie_coloc.R exactly

# ---------------------------------------------------------------------------
# Load inputs
# ---------------------------------------------------------------------------
cat("Loading COLOC results...\n")
coloc_all <- fread(COLOC_FILE)
cat("  Total pairs:", nrow(coloc_all), "| GWAS:", uniqueN(coloc_all$gwas_name), "\n")

# Focus on pairs with PP.H4 > 0.3 at the current 1Mb window
# (includes hits + near-misses to characterise the full landscape)
hits <- coloc_all[PP.H4.abf > 0.3 & !is.na(top_snp) & top_snp != ""]
cat("  Pairs to retest (PP.H4 > 0.3):", nrow(hits), "\n")
cat("  Unique genes:", uniqueN(hits$gene),
    "| Unique GWAS:", uniqueN(hits$gwas_name), "\n")

# Parse top_snp "chr:pos" -> integer position
hits[, c("top_chr", "top_pos") := tstrsplit(top_snp, ":", type.convert = TRUE)]

# Load GWAS registry — comment.char="#" skips archived study lines
registry <- fread(REGISTRY, sep = "\t", fill = TRUE, comment.char = "#")
registry <- registry[!is.na(ancestry)]

# Keep only EUR GWAS (EAS excluded: different LD, not comparable sensitivity)
eur_gwas <- registry[ancestry == "EUR", study_name]
hits <- hits[gwas_name %in% eur_gwas]
cat("  EUR GWAS pairs after filtering:", nrow(hits), "\n")

# ---------------------------------------------------------------------------
# Main loop: chromosome by chromosome
# ---------------------------------------------------------------------------
all_results <- vector("list", 22L)

for (chr_num in 1:22) {
  hits_chr <- hits[top_chr == chr_num]
  if (nrow(hits_chr) == 0) next

  cat(sprintf("\n=== Chr %d: %d pairs across %d GWAS ===\n",
              chr_num, nrow(hits_chr), uniqueN(hits_chr$gwas_name)))

  # Load eQTL data for this chromosome once
  eqtl_file <- file.path(EQTL_DIR,
                          paste0("chr", chr_num, "_marginal_summary_results.tsv"))
  if (!file.exists(eqtl_file)) {
    cat("  No eQTL file, skipping\n"); next
  }
  eqtl_chr <- fread(eqtl_file,
                    select = c("ENSG", "CHR", "POS", "NEA", "EA", "Beta", "SE"))
  setnames(eqtl_chr, c("Beta", "SE"), c("eqtl_beta", "eqtl_se"))
  cat("  eQTL rows:", nrow(eqtl_chr), "\n")

  # Process each GWAS separately (different sumstat files)
  chr_results <- vector("list", length(unique(hits_chr$gwas_name)))
  g_idx <- 0L

  for (gwas in unique(hits_chr$gwas_name)) {
    g_idx <- g_idx + 1L
    hits_gw <- hits_chr[gwas_name == gwas]
    reg <- registry[study_name == gwas]
    if (nrow(reg) == 0) next

    sumstats_path <- file.path(FM_DIR, reg$sumstats_path)
    if (!file.exists(sumstats_path)) {
      cat("  MISSING sumstats:", sumstats_path, "\n"); next
    }

    gwas_type  <- ifelse(reg$trait_type == "binary", "cc", "quant")
    gwas_n     <- reg$N_tot
    gwas_s     <- if (gwas_type == "cc") reg$N_cases / reg$N_tot else NULL

    # Load GWAS chromosome slice
    gwas_dat <- fread(sumstats_path,
                      select = c("chromosome", "position",
                                 "allele1", "allele2", "beta", "se", "pval"))
    gwas_dat <- gwas_dat[chromosome == chr_num]
    gwas_dat[, merge_key := paste(chromosome, position, sep = ":")]
    cat(sprintf("  %s: %d GWAS variants, %d gene-pairs\n",
                gwas, nrow(gwas_dat), nrow(hits_gw)))

    # For each gene × window size
    pair_results <- vector("list", nrow(hits_gw) * length(WINDOWS_BP))
    p_idx <- 0L

    for (i in seq_len(nrow(hits_gw))) {
      row       <- hits_gw[i]
      ensg      <- row$ensembl
      gene_sym  <- row$gene
      center    <- row$top_pos   # SNP with highest SNP.PP.H4 in original run
      pp4_orig  <- row$PP.H4.abf

      eqtl_gene <- eqtl_chr[ENSG == ensg]
      if (nrow(eqtl_gene) == 0) next

      for (half_win in WINDOWS_BP) {
        p_idx <- p_idx + 1L
        lo <- center - half_win
        hi <- center + half_win

        # Subset both datasets to window
        eq_w  <- eqtl_gene[POS >= lo & POS <= hi]
        gw_w  <- gwas_dat[position >= lo & position <= hi]

        if (nrow(eq_w) == 0 || nrow(gw_w) == 0) {
          pair_results[[p_idx]] <- data.table(
            gene = gene_sym, ensembl = ensg, gwas_name = gwas,
            half_window_kb = half_win / 1e3, n_snps = 0L,
            PP.H4.abf = NA_real_, pp4_1mb = pp4_orig)
          next
        }

        # Merge on position
        eq_w[,  merge_key := paste(CHR, POS, sep = ":")]
        merged <- merge(
          eq_w[,  .(merge_key, eqtl_ea = EA, eqtl_nea = NEA,
                    eqtl_beta, eqtl_se)],
          gw_w[,  .(merge_key, gwas_a1 = allele1, gwas_a2 = allele2,
                    gwas_beta = beta, gwas_se = se, gwas_pval = pval)],
          by = "merge_key"
        )
        # Keep best p-value per position, resolve alleles
        merged <- merged[order(gwas_pval)][!duplicated(merge_key)]
        merged[, allele_match := (eqtl_ea == gwas_a1 & eqtl_nea == gwas_a2)]
        merged[, allele_flip  := (eqtl_ea == gwas_a2 & eqtl_nea == gwas_a1)]
        merged <- merged[allele_match | allele_flip]
        merged[allele_flip == TRUE, eqtl_beta := -eqtl_beta]

        # Remove palindromic SNPs (A/T and C/G pairs)
        merged <- merged[!((gwas_a1 %in% c("A","T") & gwas_a2 %in% c("A","T")) |
                           (gwas_a1 %in% c("C","G") & gwas_a2 %in% c("C","G")))]

        if (nrow(merged) < MIN_SNPS) {
          pair_results[[p_idx]] <- data.table(
            gene = gene_sym, ensembl = ensg, gwas_name = gwas,
            half_window_kb = half_win / 1e3, n_snps = nrow(merged),
            PP.H4.abf = NA_real_, pp4_1mb = pp4_orig)
          next
        }

        d1 <- list(beta    = merged$gwas_beta,
                   varbeta = merged$gwas_se^2,
                   N       = gwas_n,
                   type    = gwas_type,
                   snp     = merged$merge_key)
        if (gwas_type == "cc" && gwas_s > 0) d1$s <- gwas_s
        if (gwas_type == "quant") d1$sdY <- 1   # match 06_susie_coloc.R

        d2 <- list(beta    = merged$eqtl_beta,
                   varbeta = merged$eqtl_se^2,
                   N       = EQTL_N,
                   type    = "quant",
                   sdY     = 1,
                   snp     = merged$merge_key)

        res  <- tryCatch(
          suppressMessages(suppressWarnings(
            coloc.abf(d1, d2, p1 = COLOC_P1, p2 = COLOC_P2, p12 = COLOC_P12)
          )),
          error = function(e) NULL
        )
        pp4  <- if (!is.null(res)) unname(res$summary["PP.H4.abf"]) else NA_real_

        pair_results[[p_idx]] <- data.table(
          gene      = gene_sym,
          ensembl   = ensg,
          gwas_name = gwas,
          half_window_kb = half_win / 1e3,
          n_snps    = nrow(merged),
          PP.H4.abf = pp4,
          pp4_1mb   = pp4_orig)
      }
    }

    chr_results[[g_idx]] <- rbindlist(pair_results[!sapply(pair_results, is.null)],
                                      fill = TRUE)
  }

  all_results[[chr_num]] <- rbindlist(chr_results[!sapply(chr_results, is.null)],
                                      fill = TRUE)
  cat(sprintf("  Chr %d done: %d results\n", chr_num,
              nrow(all_results[[chr_num]])))

  # Release chromosome eQTL from memory
  rm(eqtl_chr); gc()
}

# ---------------------------------------------------------------------------
# Save and summarise
# ---------------------------------------------------------------------------
final <- rbindlist(all_results[!sapply(all_results, is.null)], fill = TRUE)
fwrite(final, OUT_FILE)
cat("\nSaved:", nrow(final), "rows to", OUT_FILE, "\n")

cat("\n=== Window Sensitivity Summary (EUR GWAS only) ===\n")
cat(sprintf("%-12s  %6s  %6s  %6s  %8s\n",
            "Window", ">0.5", ">0.8", ">0.9", "n_pairs"))
cat(strrep("-", 50), "\n")
for (w in WINDOWS_BP / 1e3) {
  sub <- final[half_window_kb == w & !is.na(PP.H4.abf)]
  # Count unique genes (gene-level, not pair-level)
  gl  <- sub[, .(max_pp4 = max(PP.H4.abf, na.rm = TRUE)), by = gene]
  cat(sprintf("±%-10gkb  %6d  %6d  %6d  %8d\n",
              w, gl[max_pp4 > 0.5, .N], gl[max_pp4 > 0.8, .N],
              gl[max_pp4 > 0.9, .N], nrow(sub)))
}

# Retention of original PP.H4>0.5 hits at narrower windows
orig_hits <- unique(final[pp4_1mb > 0.5, gene])
cat(sprintf("\nOriginal PP.H4>0.5 genes (1Mb): %d\n", length(orig_hits)))
for (w in WINDOWS_BP / 1e3) {
  sub <- final[half_window_kb == w & gene %in% orig_hits & !is.na(PP.H4.abf)]
  gl  <- sub[, .(max_pp4 = max(PP.H4.abf, na.rm = TRUE)), by = gene]
  retained <- gl[max_pp4 > 0.5, .N]
  cat(sprintf("  Retained at ±%gkb: %d / %d (%.1f%%)\n",
              w, retained, length(orig_hits),
              100 * retained / length(orig_hits)))
}

cat("\nDone.\n")
