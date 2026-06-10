#!/usr/bin/env Rscript
# 35t_susie_abf_comparison.R
# ---------------------------------------------------------------------------
# SuSiE-COLOC vs ABF-only COLOC Comparison
#
# Compares colocalization posterior probabilities from SuSiE fine-mapping
# (multi-causal-variant model) against standard ABF-only COLOC across 6 GWAS.
# Identifies multi-signal loci uniquely resolvable by SuSiE, quantifies
# concordance at the PP.H4 > 0.5 threshold, and reports SuSiE fallback rates.
#
# Inputs:
#   - SuSiE results: RNA-seq/results/causal_inference/susie_broadaway_{gwas}/susie_coloc_results.csv
#   - ABF results:   RNA-seq/results/causal_inference/{abf_dir}/coloc_results.csv
#
# Outputs (to RNA-seq/results/causal_inference/susie_comparison/):
#   - susie_abf_comparison.csv   — per-gene, per-GWAS matched PP.H4 values
#   - susie_abf_summary.csv      — per-GWAS concordance & summary statistics
#   - susie_multi_signal_genes.csv — genes with n_cs_pairs > 1 (SuSiE-unique)
#
# No figures — see figS_susie_comparison.R for visualization.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

# Source shared gene-symbol mapper
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

CAUSAL  <- file.path(BASE, "RNA-seq/results/causal_inference")
OUTDIR  <- file.path(CAUSAL, "susie_comparison")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

cat("=== Script 35t: SuSiE vs ABF COLOC Comparison ===\n")
cat("Start time:", format(Sys.time()), "\n\n")

# ==============================================================================
# 1. Define GWAS registry — maps GWAS name to SuSiE and ABF directories
# ==============================================================================
gwas_registry <- data.table(
  gwas_name  = c("UKBB_ALT",  "UKBB_AST",  "UKBB_GGT",
                 "FINNGEN_NAFLD", "FINNGEN_NASH", "FINNGEN_HCC"),
  susie_dir  = c("susie_broadaway_ukbb_alt",  "susie_broadaway_ukbb_ast",  "susie_broadaway_ukbb_ggt",
                 "susie_broadaway_finngen_nafld", "susie_broadaway_finngen_nash", "susie_broadaway_finngen_hcc"),
  abf_dir    = c("broadaway_ukbb",  "broadaway_ukbb_ast",  "broadaway_ukbb_ggt",
                 "finngen_nafld", "finngen_nash", "finngen_hcc")
)

cat("--- Step 1: Loading results for", nrow(gwas_registry), "GWAS ---\n")

# ==============================================================================
# 2. Load SuSiE and ABF results for each GWAS
# ==============================================================================
all_comparisons <- list()
all_summaries   <- list()
all_multisignal <- list()

for (i in seq_len(nrow(gwas_registry))) {
  gname     <- gwas_registry$gwas_name[i]
  susie_path <- file.path(CAUSAL, gwas_registry$susie_dir[i], "susie_coloc_results.csv")
  abf_path   <- file.path(CAUSAL, gwas_registry$abf_dir[i], "coloc_results.csv")

  cat("\n  [", gname, "]\n")

  # --- Check SuSiE file ---
  if (!file.exists(susie_path)) {
    cat("    SuSiE results not found:", susie_path, "\n")
    cat("    Skipping.\n")
    next
  }
  susie <- fread(susie_path)
  if (nrow(susie) == 0) {
    cat("    SuSiE results empty. Skipping.\n")
    next
  }

  # Validate expected columns
  required_susie <- c("gene", "ensembl", "PP.H4.abf", "PP.H4.susie", "n_cs_pairs", "method")
  missing_cols <- setdiff(required_susie, names(susie))
  if (length(missing_cols) > 0) {
    cat("    SuSiE file missing columns:", paste(missing_cols, collapse = ", "), "\n")
    cat("    Skipping.\n")
    next
  }

  cat("    SuSiE: ", nrow(susie), " genes (", sum(susie$method == "susie"), " SuSiE, ",
      sum(susie$method == "abf_fallback"), " ABF fallback)\n", sep = "")

  # --- Check ABF file ---
  if (!file.exists(abf_path)) {
    cat("    ABF results not found:", abf_path, "\n")
    cat("    Will compare internal ABF vs SuSiE only (no original ABF cross-check).\n")
    abf <- NULL
  } else {
    abf <- fread(abf_path)
    cat("    ABF:   ", nrow(abf), " genes\n", sep = "")
  }

  # ============================================================================
  # 2a. Per-gene comparison table
  # ============================================================================
  comp <- susie[, .(gene, ensembl,
                    PP.H4.abf.susie_script = PP.H4.abf,
                    PP.H4.susie            = PP.H4.susie,
                    n_cs_pairs             = n_cs_pairs,
                    method                 = method)]

  # n_snps_triple may not exist in all SuSiE outputs
  if ("n_snps_triple" %in% names(susie)) {
    comp[, n_snps_triple := susie$n_snps_triple]
  } else {
    comp[, n_snps_triple := NA_integer_]
  }

  # Merge with original ABF PP.H4 if available
  if (!is.null(abf) && nrow(abf) > 0) {
    abf_sub <- abf[, .(gene, PP.H4.abf.original = PP.H4)]
    comp <- merge(comp, abf_sub, by = "gene", all.x = TRUE)
  } else {
    comp[, PP.H4.abf.original := NA_real_]
  }

  comp[, gwas := gname]

  # Classification columns
  comp[, abf_sig   := PP.H4.abf.susie_script > 0.5]
  comp[, susie_sig := PP.H4.susie > 0.5]
  comp[, concordant := abf_sig == susie_sig]
  comp[, category := fcase(
    abf_sig & susie_sig,   "Both_sig",
    !abf_sig & !susie_sig, "Both_nonsig",
    abf_sig & !susie_sig,  "ABF_only",
    !abf_sig & susie_sig,  "SuSiE_only"
  )]

  all_comparisons[[gname]] <- comp

  # ============================================================================
  # 2b. Per-GWAS summary statistics
  # ============================================================================
  n_total       <- nrow(comp)
  n_susie_ok    <- sum(comp$method == "susie")
  n_abf_fallback <- sum(comp$method == "abf_fallback")
  fallback_rate <- n_abf_fallback / n_total

  # Spearman rho — only for genes where SuSiE succeeded
  susie_ok <- comp[method == "susie"]
  if (nrow(susie_ok) >= 3) {
    rho_test <- cor.test(susie_ok$PP.H4.abf.susie_script,
                         susie_ok$PP.H4.susie,
                         method = "spearman", exact = FALSE)
    spearman_rho <- rho_test$estimate
    spearman_p   <- rho_test$p.value
  } else {
    spearman_rho <- NA_real_
    spearman_p   <- NA_real_
  }

  # Concordance at PP.H4 > 0.5
  concordance_pct <- 100 * mean(comp$concordant, na.rm = TRUE)

  # Category counts
  n_both_sig    <- sum(comp$category == "Both_sig", na.rm = TRUE)
  n_both_nonsig <- sum(comp$category == "Both_nonsig", na.rm = TRUE)
  n_abf_only    <- sum(comp$category == "ABF_only", na.rm = TRUE)
  n_susie_only  <- sum(comp$category == "SuSiE_only", na.rm = TRUE)

  # Multi-signal loci
  n_multi_signal <- sum(comp$n_cs_pairs > 1, na.rm = TRUE)

  # ABF original vs SuSiE-script ABF correlation (sanity check)
  if (sum(!is.na(comp$PP.H4.abf.original)) >= 3) {
    abf_check <- cor.test(comp[!is.na(PP.H4.abf.original)]$PP.H4.abf.susie_script,
                          comp[!is.na(PP.H4.abf.original)]$PP.H4.abf.original,
                          method = "spearman", exact = FALSE)
    abf_abf_rho <- abf_check$estimate
  } else {
    abf_abf_rho <- NA_real_
  }

  summary_row <- data.table(
    gwas             = gname,
    n_genes          = n_total,
    n_susie_ok       = n_susie_ok,
    n_abf_fallback   = n_abf_fallback,
    fallback_rate    = round(fallback_rate, 3),
    spearman_rho     = round(spearman_rho, 4),
    spearman_p       = spearman_p,
    concordance_pct  = round(concordance_pct, 1),
    n_both_sig       = n_both_sig,
    n_both_nonsig    = n_both_nonsig,
    n_abf_only       = n_abf_only,
    n_susie_only     = n_susie_only,
    n_multi_signal   = n_multi_signal,
    abf_abf_rho      = round(abf_abf_rho, 4)
  )

  all_summaries[[gname]] <- summary_row

  cat("    Spearman rho (ABF vs SuSiE):   ", round(spearman_rho, 4), "\n")
  cat("    Concordance at PP.H4 > 0.5:    ", round(concordance_pct, 1), "%\n")
  cat("    Multi-signal loci (n_cs > 1):  ", n_multi_signal, "\n")
  cat("    Fallback rate:                 ", round(100 * fallback_rate, 1), "%\n")
  cat("    Categories — Both sig:", n_both_sig, "| Both nonsig:", n_both_nonsig,
      "| ABF-only:", n_abf_only, "| SuSiE-only:", n_susie_only, "\n")

  # ============================================================================
  # 2c. Multi-signal genes
  # ============================================================================
  multi <- comp[n_cs_pairs > 1]
  if (nrow(multi) > 0) {
    all_multisignal[[gname]] <- multi
  }
}

# ==============================================================================
# 3. Combine and save outputs
# ==============================================================================
cat("\n--- Step 3: Combining and saving outputs ---\n")

if (length(all_comparisons) == 0) {
  cat("WARNING: No SuSiE results found for any GWAS. No outputs written.\n")
  cat("Expected files at: ", CAUSAL, "/susie_broadaway_*/susie_coloc_results.csv\n")
  cat("End time:", format(Sys.time()), "\n")
  quit(status = 0)
}

# 3a. Per-gene comparison
comparison_dt <- rbindlist(all_comparisons, use.names = TRUE, fill = TRUE)
comparison_dt <- add_symbols(comparison_dt, gene_col = "ensembl")
comparison_dt <- comparison_dt[order(gwas, -PP.H4.susie)]
# Alias to the column name the figure (figS_susie_comparison.R) expects:
#   PP.H4.abf — the ABF posterior computed inside the SuSiE script on the same
#   merged SNP set (PP.H4.abf.susie_script), the apples-to-apples comparator
#   against PP.H4.susie. PP.H4.abf.susie_script / .original are retained above.
comparison_dt[, PP.H4.abf := PP.H4.abf.susie_script]
fwrite(comparison_dt, file.path(OUTDIR, "susie_abf_comparison.csv"))
cat("  Saved susie_abf_comparison.csv:", nrow(comparison_dt), "rows across",
    uniqueN(comparison_dt$gwas), "GWAS\n")

# 3b. Per-GWAS summary
summary_dt <- rbindlist(all_summaries, use.names = TRUE)
fwrite(summary_dt, file.path(OUTDIR, "susie_abf_summary.csv"))
cat("  Saved susie_abf_summary.csv:", nrow(summary_dt), "GWAS\n")

# Print summary table
cat("\n  Per-GWAS summary:\n")
cat("  ", paste(rep("-", 110), collapse = ""), "\n")
cat(sprintf("  %-16s %6s %6s %8s %8s %8s %6s %6s %6s %6s\n",
            "GWAS", "Genes", "SuSiE", "Fallback", "Rho", "Concord", "Both", "ABF+", "SuSiE+", "Multi"))
cat("  ", paste(rep("-", 110), collapse = ""), "\n")
for (j in seq_len(nrow(summary_dt))) {
  r <- summary_dt[j]
  cat(sprintf("  %-16s %6d %6d %7.1f%% %8.4f %7.1f%% %6d %6d %6d %6d\n",
              r$gwas, r$n_genes, r$n_susie_ok,
              100 * r$fallback_rate, r$spearman_rho,
              r$concordance_pct, r$n_both_sig,
              r$n_abf_only, r$n_susie_only, r$n_multi_signal))
}
cat("  ", paste(rep("-", 110), collapse = ""), "\n")

# 3c. Multi-signal genes
if (length(all_multisignal) > 0) {
  multi_dt <- rbindlist(all_multisignal, use.names = TRUE, fill = TRUE)
  multi_dt <- add_symbols(multi_dt, gene_col = "ensembl")
  multi_dt <- multi_dt[order(-n_cs_pairs, gwas)]
  fwrite(multi_dt, file.path(OUTDIR, "susie_multi_signal_genes.csv"))
  cat("  Saved susie_multi_signal_genes.csv:", nrow(multi_dt), "entries (",
      uniqueN(multi_dt$gene), "unique genes)\n")

  # Report top multi-signal genes
  cat("\n  Top multi-signal genes (n_cs_pairs > 1):\n")
  top_multi <- multi_dt[, .(max_cs = max(n_cs_pairs),
                            n_gwas = uniqueN(gwas),
                            gwas_list = paste(unique(gwas), collapse = ";")),
                        by = .(gene, symbol)]
  top_multi <- top_multi[order(-max_cs, -n_gwas)]
  for (k in seq_len(min(20, nrow(top_multi)))) {
    r <- top_multi[k]
    cat(sprintf("    %s (%s): max %d CS pairs across %d GWAS [%s]\n",
                r$symbol, r$gene, r$max_cs, r$n_gwas, r$gwas_list))
  }
} else {
  cat("  No multi-signal genes found across any GWAS.\n")
  # Write empty file for downstream consistency
  fwrite(data.table(gene = character(), ensembl = character(), gwas = character(),
                    n_cs_pairs = integer(), PP.H4.susie = numeric(), symbol = character()),
         file.path(OUTDIR, "susie_multi_signal_genes.csv"))
}

# ==============================================================================
# 4. Cross-GWAS aggregation
# ==============================================================================
cat("\n--- Step 4: Cross-GWAS aggregation ---\n")

# Genes that gain or lose significance with SuSiE
gained <- comparison_dt[category == "SuSiE_only", .(gene, symbol, gwas, PP.H4.abf.susie_script, PP.H4.susie, n_cs_pairs)]
lost   <- comparison_dt[category == "ABF_only",   .(gene, symbol, gwas, PP.H4.abf.susie_script, PP.H4.susie, n_cs_pairs)]

cat("  Genes gaining significance with SuSiE (SuSiE-only):", nrow(gained), "\n")
if (nrow(gained) > 0) {
  cat("    Top 10:\n")
  gained_top <- gained[order(-PP.H4.susie)][seq_len(min(10, nrow(gained)))]
  for (k in seq_len(nrow(gained_top))) {
    r <- gained_top[k]
    cat(sprintf("      %s | %s | ABF=%.3f → SuSiE=%.3f | CS=%d\n",
                r$symbol, r$gwas, r$PP.H4.abf.susie_script, r$PP.H4.susie,
                r$n_cs_pairs))
  }
}

cat("  Genes losing significance with SuSiE (ABF-only):", nrow(lost), "\n")
if (nrow(lost) > 0) {
  cat("    Top 10:\n")
  lost_top <- lost[order(-PP.H4.abf.susie_script)][seq_len(min(10, nrow(lost)))]
  for (k in seq_len(nrow(lost_top))) {
    r <- lost_top[k]
    cat(sprintf("      %s | %s | ABF=%.3f → SuSiE=%.3f | CS=%d\n",
                r$symbol, r$gwas, r$PP.H4.abf.susie_script, r$PP.H4.susie,
                r$n_cs_pairs))
  }
}

# ==============================================================================
# 5. Aggregate summary
# ==============================================================================
cat("\n--- Step 5: Aggregate summary ---\n")

agg <- summary_dt[, .(
  n_gwas           = .N,
  total_genes      = sum(n_genes),
  total_susie_ok   = sum(n_susie_ok),
  mean_fallback    = round(mean(fallback_rate), 3),
  mean_rho         = round(mean(spearman_rho, na.rm = TRUE), 4),
  mean_concordance = round(mean(concordance_pct, na.rm = TRUE), 1),
  total_both_sig   = sum(n_both_sig),
  total_abf_only   = sum(n_abf_only),
  total_susie_only = sum(n_susie_only),
  total_multi      = sum(n_multi_signal)
)]

cat("  GWAS analyzed:            ", agg$n_gwas, "\n")
cat("  Total gene-GWAS pairs:    ", agg$total_genes, "\n")
cat("  SuSiE success:            ", agg$total_susie_ok,
    " (", round(100 * agg$total_susie_ok / agg$total_genes, 1), "%)\n", sep = "")
cat("  Mean fallback rate:       ", round(100 * agg$mean_fallback, 1), "%\n")
cat("  Mean Spearman rho:        ", agg$mean_rho, "\n")
cat("  Mean concordance (>0.5):  ", agg$mean_concordance, "%\n")
cat("  Total colocalizing (both):", agg$total_both_sig, "\n")
cat("  ABF-only signal:          ", agg$total_abf_only, "\n")
cat("  SuSiE-only signal:        ", agg$total_susie_only, "\n")
cat("  Multi-signal loci:        ", agg$total_multi, "\n")

cat("\nEnd time:", format(Sys.time()), "\n")
cat("Outputs in:", OUTDIR, "\n")
cat("=== Done ===\n")
