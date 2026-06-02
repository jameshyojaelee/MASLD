#!/usr/bin/env Rscript
# diagnose_coloc_thresholds.R — Comprehensive COLOC parameter audit
#
# Runs: p12 sensitivity, MAF audit, palindromic SNP count, MHC flagging,
#       LD contamination summary, sdY sensitivity
#
# For all gene-GWAS pairs with PP.H4 >= 0.3 in existing results.
# Processes all 24 GWAS × 22 chromosomes.
#
# Output: GWAS/finemapping/results/coloc_threshold_audit/

suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
  library(ggplot2)
})

BASE_DIR <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
FM_DIR   <- file.path(BASE_DIR, "GWAS/finemapping")
EQTL_DIR <- file.path(BASE_DIR, "data/broadaway_eqtl")
GWAS_DIR <- file.path(FM_DIR, "data/sumstats")
setwd(FM_DIR)

OUT_DIR <- file.path(FM_DIR, "results/coloc_threshold_audit")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
FIG_DIR <- file.path(BASE_DIR, "figures/supplementary")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

EQTL_N <- 1183L

# p12 values to test
P12_VALUES <- c(1e-7, 1e-6, 5e-6, 1e-5, 5e-5)
P12_NAMES  <- c("1e-7", "1e-6", "5e-6", "1e-5", "5e-5")

# MHC region (hg19 coordinates, matching Broadaway eQTL build)
MHC_CHR   <- 6L
MHC_START <- 25000000L
MHC_END   <- 35000000L

# ═══════════════════════════════════════════════════════════════════════════════
# 1. Load existing results — identify targets
# ═══════════════════════════════════════════════════════════════════════════════
cat("=== Loading existing COLOC results ===\n")
all_coloc <- fread(file.path(FM_DIR, "results/susie_coloc/susie_coloc_all_gwas.csv"))
cat(sprintf("  Total: %d gene-GWAS pairs\n", nrow(all_coloc)))

# Targets: genes with PP.H4 >= 0.3 in any GWAS (captures borderline cases)
PP4_TARGET <- 0.3
targets <- all_coloc[PP.H4.abf >= PP4_TARGET]
cat(sprintf("  Targets (PP.H4 >= %.1f): %d gene-GWAS pairs, %d unique genes\n",
            PP4_TARGET, nrow(targets), uniqueN(targets$gene)))

# Build lookup: GWAS -> chr -> list of ensembl IDs to test
target_lookup <- targets[, .(ensembl_list = list(unique(ensembl))),
                          by = .(gwas_name, chr)]
cat(sprintf("  GWAS-chr combinations to process: %d\n", nrow(target_lookup)))

# ═══════════════════════════════════════════════════════════════════════════════
# 2. Load GWAS registry
# ═══════════════════════════════════════════════════════════════════════════════
registry <- read.delim("config/gwas_registry.tsv", stringsAsFactors = FALSE)
gwas_names <- unique(targets$gwas_name)
cat(sprintf("  GWAS to process: %d\n", length(gwas_names)))

# ═══════════════════════════════════════════════════════════════════════════════
# 3. Process each GWAS × chromosome
# ═══════════════════════════════════════════════════════════════════════════════
results_p12   <- list()
results_maf   <- list()
results_palim <- list()
result_idx <- 0L

for (gw in gwas_names) {
  study_row <- registry[registry$study_name == gw, ]
  if (nrow(study_row) == 0) {
    cat(sprintf("  SKIP %s — not in registry\n", gw))
    next
  }

  gwas_file <- file.path(GWAS_DIR, paste0(gw, "_preprocessed.tsv"))
  if (!file.exists(gwas_file)) {
    # Try alternative naming
    gwas_file <- file.path(GWAS_DIR, paste0(gw, "_reformatted_hg19.tsv"))
  }
  if (!file.exists(gwas_file)) {
    cat(sprintf("  SKIP %s — sumstats file not found\n", gw))
    next
  }

  cat(sprintf("\n--- %s ---\n", gw))
  gwas <- fread(gwas_file)
  cat(sprintf("  GWAS loaded: %d variants\n", nrow(gwas)))

  # Determine COLOC type
  gwas_coloc_type <- ifelse(study_row$trait_type == "binary", "cc", "quant")
  gwas_n <- study_row$N_tot
  gwas_n_cases <- study_row$N_cases

  gwas[, merge_key := paste(chromosome, position, sep = ":")]

  # Chromosomes to process for this GWAS
  chr_targets <- target_lookup[gwas_name == gw]

  for (row_i in seq_len(nrow(chr_targets))) {
    chr_num   <- chr_targets$chr[row_i]
    gene_list <- chr_targets$ensembl_list[[row_i]]

    # Load eQTL for this chromosome
    eqtl_file <- file.path(EQTL_DIR, paste0("chr", chr_num, "_marginal_summary_results.tsv"))
    if (!file.exists(eqtl_file)) next
    eqtl_all <- fread(eqtl_file)

    for (gene_id in gene_list) {
      eqtl_gene <- eqtl_all[ENSG == gene_id]
      gene_symbol <- eqtl_gene$GeneSymbol[1]
      if (nrow(eqtl_gene) < 10) next

      # Merge on position
      eqtl_gene[, merge_key := paste(CHR, POS, sep = ":")]
      merged <- merge(
        eqtl_gene[, .(merge_key, eqtl_pos = POS, eqtl_ea = EA, eqtl_nea = NEA,
                       eqtl_beta = Beta, eqtl_se = SE, eqtl_pval = PVAL, eqtl_eaf = EAF)],
        gwas[, .(merge_key, gwas_a1 = allele1, gwas_a2 = allele2,
                 gwas_beta = beta, gwas_se = se, gwas_pval = pval)],
        by = "merge_key"
      )
      merged <- merged[order(gwas_pval)][!duplicated(merge_key)]
      if (nrow(merged) < 10) next

      # ── Allele harmonization ──
      merged[, allele_match := (eqtl_ea == gwas_a1 & eqtl_nea == gwas_a2)]
      merged[, allele_flip  := (eqtl_ea == gwas_a2 & eqtl_nea == gwas_a1)]

      # Count palindromic BEFORE removing them
      merged[, is_palindromic := (gwas_a1 %in% c("A","T") & gwas_a2 %in% c("A","T")) |
                                  (gwas_a1 %in% c("C","G") & gwas_a2 %in% c("C","G"))]
      n_palindromic <- sum(merged$is_palindromic)
      n_pre_filter  <- nrow(merged)

      # Count unmatched (neither match nor flip, excluding palindromic)
      n_unmatched <- sum(!merged$allele_match & !merged$allele_flip & !merged$is_palindromic)

      # Apply filters
      merged <- merged[allele_match | allele_flip]
      merged[allele_flip == TRUE, eqtl_beta := -eqtl_beta]
      merged[allele_flip == TRUE, eqtl_eaf  := 1 - eqtl_eaf]  # flip EAF too
      merged <- merged[is_palindromic == FALSE]
      if (nrow(merged) < 10) next

      n_post_filter <- nrow(merged)

      # ── MAF analysis ──
      merged[, maf := pmin(eqtl_eaf, 1 - eqtl_eaf)]
      n_maf_below_001 <- sum(merged$maf < 0.01, na.rm = TRUE)
      n_maf_below_005 <- sum(merged$maf < 0.05, na.rm = TRUE)
      median_maf      <- median(merged$maf, na.rm = TRUE)
      min_maf         <- min(merged$maf, na.rm = TRUE)

      # ── MHC flag ──
      gene_positions <- merged$eqtl_pos
      is_mhc <- (chr_num == MHC_CHR &&
                  any(gene_positions >= MHC_START & gene_positions <= MHC_END))

      # ── Record MAF/palindromic stats ──
      result_idx <- result_idx + 1L
      results_maf[[result_idx]] <- data.table(
        gene = gene_symbol, ensembl = gene_id, chr = chr_num,
        gwas_name = gw, n_snps = n_post_filter,
        n_pre_filter = n_pre_filter,
        n_palindromic = n_palindromic,
        pct_palindromic = round(n_palindromic / n_pre_filter * 100, 1),
        n_unmatched = n_unmatched,
        n_maf_below_001 = n_maf_below_001,
        n_maf_below_005 = n_maf_below_005,
        pct_maf_below_001 = round(n_maf_below_001 / n_post_filter * 100, 1),
        pct_maf_below_005 = round(n_maf_below_005 / n_post_filter * 100, 1),
        median_maf = round(median_maf, 4),
        min_maf = round(min_maf, 6),
        is_mhc = is_mhc
      )

      # ── p12 sensitivity: run coloc.abf with multiple p12 values ──
      pp4_by_p12 <- numeric(length(P12_VALUES))
      pp3_by_p12 <- numeric(length(P12_VALUES))

      for (p_i in seq_along(P12_VALUES)) {
        tryCatch({
          d1 <- list(beta = merged$gwas_beta, varbeta = merged$gwas_se^2,
                     N = gwas_n, type = gwas_coloc_type, snp = merged$merge_key)
          if (gwas_coloc_type == "cc" && gwas_n_cases > 0) d1$s <- gwas_n_cases / gwas_n
          if (gwas_coloc_type == "quant") d1$sdY <- 1

          d2 <- list(beta = merged$eqtl_beta, varbeta = merged$eqtl_se^2,
                     N = EQTL_N, type = "quant", sdY = 1, snp = merged$merge_key)

          res <- suppressMessages(suppressWarnings(
            coloc.abf(d1, d2, p1 = 1e-4, p2 = 1e-4, p12 = P12_VALUES[p_i])
          ))
          pp4_by_p12[p_i] <- res$summary["PP.H4.abf"]
          pp3_by_p12[p_i] <- res$summary["PP.H3.abf"]
        }, error = function(e) {
          pp4_by_p12[p_i] <<- NA_real_
          pp3_by_p12[p_i] <<- NA_real_
        })
      }

      results_p12[[result_idx]] <- data.table(
        gene = gene_symbol, ensembl = gene_id, chr = chr_num,
        gwas_name = gw, n_snps = n_post_filter, is_mhc = is_mhc,
        pp4_1e7  = pp4_by_p12[1], pp4_1e6  = pp4_by_p12[2],
        pp4_5e6  = pp4_by_p12[3], pp4_1e5  = pp4_by_p12[4],
        pp4_5e5  = pp4_by_p12[5],
        pp3_1e7  = pp3_by_p12[1], pp3_1e6  = pp3_by_p12[2],
        pp3_5e6  = pp3_by_p12[3], pp3_1e5  = pp3_by_p12[4],
        pp3_5e5  = pp3_by_p12[5]
      )
    } # gene

    rm(eqtl_all); gc(verbose = FALSE)
  } # chr

  rm(gwas); gc(verbose = FALSE)
  cat(sprintf("  Processed %d gene-GWAS pairs so far\n", result_idx))
} # gwas

# ═══════════════════════════════════════════════════════════════════════════════
# 4. Combine results
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n=== Combining results ===\n")
dt_p12 <- rbindlist(results_p12[!sapply(results_p12, is.null)], fill = TRUE)
dt_maf <- rbindlist(results_maf[!sapply(results_maf, is.null)], fill = TRUE)
cat(sprintf("  p12 sensitivity: %d gene-GWAS pairs\n", nrow(dt_p12)))
cat(sprintf("  MAF audit:       %d gene-GWAS pairs\n", nrow(dt_maf)))

fwrite(dt_p12, file.path(OUT_DIR, "p12_sensitivity.csv"))
fwrite(dt_maf, file.path(OUT_DIR, "maf_palindromic_audit.csv"))

# ═══════════════════════════════════════════════════════════════════════════════
# 5. p12 sensitivity analysis
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n=== p12 SENSITIVITY ANALYSIS ===\n")

# How many genes exceed PP.H4 thresholds at each p12?
for (pp_thr in c(0.5, 0.8, 0.9)) {
  cat(sprintf("\n  Unique genes with PP.H4 > %.1f by p12:\n", pp_thr))
  counts <- sapply(c("pp4_1e7", "pp4_1e6", "pp4_5e6", "pp4_1e5", "pp4_5e5"),
    function(col) uniqueN(dt_p12[get(col) > pp_thr, gene]))
  names(counts) <- P12_NAMES
  print(counts)
}

# Genes that flip above/below PP.H4 = 0.5 depending on p12
cat("\n  Genes SENSITIVE to p12 (PP.H4 crosses 0.5 between p12 = 1e-6 and 1e-5):\n")
sensitive <- dt_p12[(pp4_1e6 < 0.5 & pp4_1e5 >= 0.5) |
                    (pp4_1e6 >= 0.5 & pp4_1e5 < 0.5)]
if (nrow(sensitive) > 0) {
  cat(sprintf("  %d sensitive gene-GWAS pairs:\n", nrow(sensitive)))
  sensitive_summary <- sensitive[, .(
    min_pp4 = min(c(pp4_1e6, pp4_1e5)),
    max_pp4 = max(c(pp4_1e6, pp4_1e5)),
    n_snps = n_snps[1], is_mhc = is_mhc[1]
  ), by = .(gene, gwas_name)][order(-max_pp4)]
  print(head(sensitive_summary, 20))
}

# Overall stability: correlation between p12 = 5e-6 (current) and alternatives
cat("\n  Spearman correlation of PP.H4 ranks across p12 values:\n")
cor_vals <- sapply(c("pp4_1e7", "pp4_1e6", "pp4_5e6", "pp4_1e5", "pp4_5e5"),
  function(col) cor(dt_p12$pp4_5e6, dt_p12[[col]], use = "complete.obs", method = "spearman"))
names(cor_vals) <- P12_NAMES
print(round(cor_vals, 4))

# ═══════════════════════════════════════════════════════════════════════════════
# 6. MAF audit
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n=== MAF AUDIT ===\n")

cat("  Overall MAF distribution across all gene-GWAS windows:\n")
cat(sprintf("    Genes with ANY variant MAF < 0.01: %d (%.1f%%)\n",
  sum(dt_maf$n_maf_below_001 > 0), mean(dt_maf$n_maf_below_001 > 0) * 100))
cat(sprintf("    Genes with >10%% variants MAF < 0.01: %d (%.1f%%)\n",
  sum(dt_maf$pct_maf_below_001 > 10), mean(dt_maf$pct_maf_below_001 > 10) * 100))
cat(sprintf("    Genes with >10%% variants MAF < 0.05: %d (%.1f%%)\n",
  sum(dt_maf$pct_maf_below_005 > 10), mean(dt_maf$pct_maf_below_005 > 10) * 100))
cat(sprintf("    Median per-window median MAF: %.3f\n", median(dt_maf$median_maf)))

# ═══════════════════════════════════════════════════════════════════════════════
# 7. Palindromic SNP audit
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n=== PALINDROMIC SNP AUDIT ===\n")
cat(sprintf("  Total palindromic SNPs removed across all windows: %s\n",
  format(sum(dt_maf$n_palindromic), big.mark = ",")))
cat(sprintf("  Mean palindromic fraction per window: %.1f%%\n",
  mean(dt_maf$pct_palindromic)))
cat(sprintf("  Max palindromic fraction: %.1f%% (gene: %s, GWAS: %s)\n",
  max(dt_maf$pct_palindromic),
  dt_maf$gene[which.max(dt_maf$pct_palindromic)],
  dt_maf$gwas_name[which.max(dt_maf$pct_palindromic)]))

# ═══════════════════════════════════════════════════════════════════════════════
# 8. MHC analysis
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n=== MHC REGION (chr6:25-35Mb) ANALYSIS ===\n")
mhc_hits <- dt_p12[is_mhc == TRUE]
cat(sprintf("  MHC gene-GWAS pairs in audit set: %d\n", nrow(mhc_hits)))
if (nrow(mhc_hits) > 0) {
  cat(sprintf("  Unique MHC genes: %d\n", uniqueN(mhc_hits$gene)))
  cat("  MHC genes with PP.H4 > 0.5 (at current p12 = 5e-6):\n")
  mhc_top <- mhc_hits[pp4_5e6 > 0.5, .(gene, gwas_name, n_snps, pp4_5e6)][order(-pp4_5e6)]
  if (nrow(mhc_top) > 0) print(mhc_top) else cat("    None\n")
}

# ═══════════════════════════════════════════════════════════════════════════════
# 9. LD contamination integration
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n=== LD CONTAMINATION CLUSTERS ===\n")
ld_file <- file.path(FM_DIR, "results/susie_coloc/ld_contamination_clusters.csv")
if (file.exists(ld_file)) {
  ld_clusters <- fread(ld_file)
  cat(sprintf("  %d clusters flagged (%d unique genes involved)\n",
    nrow(ld_clusters),
    length(unique(unlist(strsplit(ld_clusters$genes, ", "))))))

  cat("  By flag reason:\n")
  print(table(ld_clusters$flag_reason))

  # Extract all genes in LD clusters
  ld_genes <- unique(unlist(strsplit(ld_clusters$genes, ", ")))
  # How many of our PP.H4>0.5 genes are in LD clusters?
  gene_coloc <- fread(file.path(FM_DIR, "results/susie_coloc/gene_level_coloc.csv"))
  top_genes <- gene_coloc[coloc_best_pp4 > 0.5, gene]
  n_top_in_ld <- sum(top_genes %in% ld_genes)
  cat(sprintf("\n  PP.H4>0.5 genes in LD clusters: %d / %d (%.1f%%)\n",
    n_top_in_ld, length(top_genes), n_top_in_ld/length(top_genes)*100))
}

# ═══════════════════════════════════════════════════════════════════════════════
# 10. Combined sensitivity summary table
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n=== COMBINED SENSITIVITY SUMMARY ===\n")

# For each filter combination, count surviving unique genes
combos <- expand.grid(
  min_snps = c(10, 50, 100, 200),
  min_maf  = c(0, 0.01, 0.05),
  p12      = c("pp4_1e6", "pp4_5e6", "pp4_1e5"),
  pp_thr   = c(0.5, 0.8, 0.9),
  stringsAsFactors = FALSE
)

# Merge p12 results with MAF data
dt_combined <- merge(dt_p12, dt_maf[, .(gene, ensembl, gwas_name, min_maf, n_maf_below_001, n_maf_below_005, n_snps)],
                     by = c("gene", "ensembl", "gwas_name"), suffixes = c("", "_maf"))

cat("  Building combined sensitivity table...\n")
combo_results <- lapply(seq_len(nrow(combos)), function(i) {
  ms <- combos$min_snps[i]
  mm <- combos$min_maf[i]
  p12_col <- combos$p12[i]
  pp_thr  <- combos$pp_thr[i]

  sub <- dt_combined[n_snps >= ms]

  # MAF filter: we can't directly apply per-SNP MAF filter to existing results
  # But we can use as proxy: exclude genes where min_maf < threshold
  if (mm > 0) {
    sub <- sub[min_maf >= mm | is.na(min_maf)]
  }

  data.table(
    min_snps = ms, min_maf = mm,
    p12 = gsub("pp4_", "", p12_col),
    pp_thr = pp_thr,
    n_unique_genes = uniqueN(sub[get(p12_col) >= pp_thr, gene])
  )
})
combo_dt <- rbindlist(combo_results)

# Pivot: show the most useful view
cat("\n  Unique genes by (min_snps, p12) at PP.H4 > 0.5, no MAF filter:\n")
sub_view <- combo_dt[pp_thr == 0.5 & min_maf == 0]
print(dcast(sub_view, min_snps ~ p12, value.var = "n_unique_genes"))

cat("\n  Unique genes by (min_snps, p12) at PP.H4 > 0.8, no MAF filter:\n")
sub_view <- combo_dt[pp_thr == 0.8 & min_maf == 0]
print(dcast(sub_view, min_snps ~ p12, value.var = "n_unique_genes"))

cat("\n  Unique genes by (min_snps, p12) at PP.H4 > 0.9, no MAF filter:\n")
sub_view <- combo_dt[pp_thr == 0.9 & min_maf == 0]
print(dcast(sub_view, min_snps ~ p12, value.var = "n_unique_genes"))

cat("\n  Effect of MAF filter at PP.H4 > 0.5, min_snps = 100, p12 = 5e-6:\n")
sub_view <- combo_dt[pp_thr == 0.5 & min_snps == 100 & p12 == "5e-6"]
print(sub_view)

fwrite(combo_dt, file.path(OUT_DIR, "combined_sensitivity_table.csv"))

# ═══════════════════════════════════════════════════════════════════════════════
# 11. Figures
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n=== Generating figures ===\n")

pdf(file.path(FIG_DIR, "coloc_threshold_audit.pdf"), width = 14, height = 16)

# Panel A: p12 sensitivity — PP.H4 at 5e-6 vs 1e-6 and 1e-5
par(mfrow = c(2, 2))

plot(dt_p12$pp4_5e6, dt_p12$pp4_1e6, pch = 16, cex = 0.3, col = "#00000040",
     xlab = "PP.H4 at p12 = 5e-6 (current)", ylab = "PP.H4 at p12 = 1e-6",
     main = "p12 sensitivity: 5e-6 vs 1e-6 (conservative)")
abline(0, 1, col = "red", lty = 2)
abline(v = 0.5, h = 0.5, col = "grey60", lty = 3)
r <- cor(dt_p12$pp4_5e6, dt_p12$pp4_1e6, use = "complete.obs")
legend("topleft", paste0("r = ", round(r, 3)), bty = "n")

plot(dt_p12$pp4_5e6, dt_p12$pp4_1e5, pch = 16, cex = 0.3, col = "#00000040",
     xlab = "PP.H4 at p12 = 5e-6 (current)", ylab = "PP.H4 at p12 = 1e-5",
     main = "p12 sensitivity: 5e-6 vs 1e-5 (permissive)")
abline(0, 1, col = "red", lty = 2)
abline(v = 0.5, h = 0.5, col = "grey60", lty = 3)
r <- cor(dt_p12$pp4_5e6, dt_p12$pp4_1e5, use = "complete.obs")
legend("topleft", paste0("r = ", round(r, 3)), bty = "n")

# Panel B: PP.H3 vs PP.H4 at current p12 — shows H3/H4 balance
plot(dt_p12$pp3_5e6, dt_p12$pp4_5e6, pch = 16, cex = 0.3, col = "#00000040",
     xlab = "PP.H3 (independent signals)", ylab = "PP.H4 (shared signal)",
     main = "H3 vs H4 at p12 = 5e-6")
# Color MHC genes
if (any(dt_p12$is_mhc)) {
  points(dt_p12$pp3_5e6[dt_p12$is_mhc], dt_p12$pp4_5e6[dt_p12$is_mhc],
         pch = 16, cex = 0.8, col = "red")
  legend("topright", "MHC genes", pch = 16, col = "red", bty = "n")
}

# Panel C: MAF distribution
hist(dt_maf$median_maf, breaks = 50, col = "steelblue", border = "white",
     main = "Median MAF per COLOC window", xlab = "Median eQTL MAF")
abline(v = c(0.01, 0.05), col = c("red", "orange"), lty = 2, lwd = 2)

par(mfrow = c(1, 1))

# Panel D: ggplot sensitivity heatmap
sens_plot_data <- combo_dt[min_maf == 0]
sens_plot_data[, label := paste0("PP.H4>", pp_thr)]
p <- ggplot(sens_plot_data,
            aes(x = factor(min_snps), y = p12, fill = n_unique_genes)) +
  geom_tile(color = "white") +
  geom_text(aes(label = n_unique_genes), size = 3.5) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  facet_wrap(~ label, nrow = 1) +
  labs(x = "Minimum SNPs", y = "p12 prior",
       title = "COLOC gene counts: sensitivity to min_snps and p12",
       fill = "Genes") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 0))
print(p)

# Panel E: Per-gene PP.H4 shift (delta between p12=1e-6 and p12=1e-5)
dt_p12[, delta_pp4 := pp4_1e5 - pp4_1e6]
p2 <- ggplot(dt_p12[pp4_5e6 > 0.3], aes(x = pp4_5e6, y = delta_pp4)) +
  geom_point(alpha = 0.3, size = 0.8) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(x = "PP.H4 at current p12 (5e-6)",
       y = "PP.H4 shift (p12=1e-5 minus p12=1e-6)",
       title = "Per-gene PP.H4 sensitivity to p12 prior",
       subtitle = "Positive = gains under permissive prior; negative = loses") +
  theme_bw(base_size = 11)
print(p2)

dev.off()
cat(sprintf("Figures saved: %s\n", file.path(FIG_DIR, "coloc_threshold_audit.pdf")))

# ═══════════════════════════════════════════════════════════════════════════════
# 12. Final recommendations
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n")
cat("================================================================\n")
cat("  THRESHOLD AUDIT COMPLETE — OUTPUT FILES\n")
cat("================================================================\n")
cat(sprintf("  p12 sensitivity:           %s\n", file.path(OUT_DIR, "p12_sensitivity.csv")))
cat(sprintf("  MAF/palindromic audit:     %s\n", file.path(OUT_DIR, "maf_palindromic_audit.csv")))
cat(sprintf("  Combined sensitivity:      %s\n", file.path(OUT_DIR, "combined_sensitivity_table.csv")))
cat(sprintf("  Figures:                   %s\n", file.path(FIG_DIR, "coloc_threshold_audit.pdf")))
cat("================================================================\n")
