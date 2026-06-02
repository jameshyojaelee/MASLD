#!/usr/bin/env Rscript
# 50e_eaf_concordance_check.R
# ---------------------------------------------------------------------------
# Cross-population EAF concordance diagnostic
# Compares European eQTL EAF (Broadaway) vs non-European GWAS EAF (Pan-UKBB)
# to quantify allele frequency divergence between populations.
#
# No COLOC is run — this is a pure merge + correlation diagnostic.
# ---------------------------------------------------------------------------

library(data.table)

cat("=== Script 50e: EAF Concordance Check ===\n")
cat("Start time:", format(Sys.time()), "\n\n")

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

EQTL_DIR  <- file.path(BASE, "data/broadaway_eqtl")
OUTDIR    <- file.path(BASE, "RNA-seq/results/causal_inference/eaf_concordance")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# GWAS files to test (AFR and CSA)
gwas_configs <- list(
  list(name = "PanUKBB_AFR_ALT", ancestry = "AFR",
       file = file.path(BASE, "GWAS/MR_Data/PanUKBB/PanUKBB_AFR_ALT_harmonised_hg38.tsv.gz")),
  list(name = "PanUKBB_CSA_ALT", ancestry = "CSA",
       file = file.path(BASE, "GWAS/MR_Data/PanUKBB/PanUKBB_CSA_ALT_harmonised_hg38.tsv.gz"))
)

# Also test EUR-EUR as positive control
eur_gwas <- file.path(BASE, "GWAS/MR_Data/GCST90019492_UKBB_ALT_harmonised.tsv.gz")
if (file.exists(eur_gwas)) {
  gwas_configs <- c(gwas_configs, list(
    list(name = "UKBB_EUR_ALT", ancestry = "EUR", file = eur_gwas)
  ))
}

# ==============================================================================
# Load all Broadaway eQTL data
# ==============================================================================
cat("--- Loading Broadaway eQTLs ---\n")
eqtl_files <- list.files(EQTL_DIR, pattern = "^chr[0-9]+_marginal_summary_results\\.tsv$",
                          full.names = TRUE)
cat("  Found", length(eqtl_files), "chromosome files\n")

eqtl_all <- rbindlist(lapply(eqtl_files, function(f) {
  dt <- fread(f, select = c("CHR", "POS", "EA", "NEA", "EAF", "GeneSymbol"))
  dt
}))
cat("  Total eQTL rows:", nrow(eqtl_all), "\n")

# Create merge key (chr:pos)
eqtl_all[, merge_key := paste0(CHR, ":", POS)]
# Deduplicate by merge_key (keep first occurrence)
eqtl_all <- eqtl_all[!duplicated(merge_key)]
cat("  Unique variants:", nrow(eqtl_all), "\n")

# ==============================================================================
# Compare EAF for each GWAS
# ==============================================================================
results <- list()

for (cfg in gwas_configs) {
  cat("\n--- Processing:", cfg$name, "(", cfg$ancestry, ") ---\n")

  if (!file.exists(cfg$file)) {
    cat("  WARNING: File not found:", cfg$file, "\n")
    next
  }

  gwas <- fread(cfg$file)
  cat("  Loaded GWAS:", nrow(gwas), "variants\n")

  # Standardize EAF column name
  if ("effect_allele_frequency" %in% names(gwas) && !"EAF" %in% names(gwas)) {
    setnames(gwas, "effect_allele_frequency", "EAF")
  }
  if (!"EAF" %in% names(gwas)) {
    cat("  WARNING: No EAF column in GWAS\n")
    next
  }

  # Create merge key
  if ("base_pair_location" %in% names(gwas)) {
    gwas[, merge_key := paste0(chromosome, ":", base_pair_location)]
  } else if ("pos_hg38" %in% names(gwas)) {
    gwas[, merge_key := paste0(chromosome, ":", pos_hg38)]
  } else {
    cat("  WARNING: Cannot identify position column\n")
    next
  }

  # Merge
  merged <- merge(eqtl_all, gwas, by = "merge_key", suffixes = c(".eqtl", ".gwas"))
  cat("  Merged variants:", nrow(merged), "\n")

  # Identify EAF columns
  eaf_eqtl <- if ("EAF.eqtl" %in% names(merged)) "EAF.eqtl" else "EAF"
  eaf_gwas <- if ("EAF.gwas" %in% names(merged)) "EAF.gwas" else NULL

  if (is.null(eaf_gwas)) {
    cat("  WARNING: EAF columns did not collide during merge\n")
    next
  }

  # Harmonize alleles before comparing EAF
  # eQTL has EA/NEA, GWAS has effect_allele/other_allele
  # When EA != effect_allele, EAF needs to be flipped (1 - EAF)
  ea_eqtl <- if ("EA" %in% names(merged)) "EA" else if ("EA.eqtl" %in% names(merged)) "EA.eqtl" else NULL
  ea_gwas <- if ("effect_allele" %in% names(merged)) "effect_allele" else if ("effect_allele.gwas" %in% names(merged)) "effect_allele.gwas" else NULL
  oa_gwas <- if ("other_allele" %in% names(merged)) "other_allele" else if ("other_allele.gwas" %in% names(merged)) "other_allele.gwas" else NULL

  if (!is.null(ea_eqtl) && !is.null(ea_gwas)) {
    # Classify allele matching
    merged[, allele_match := "none"]
    merged[toupper(get(ea_eqtl)) == toupper(get(ea_gwas)), allele_match := "direct"]
    if (!is.null(oa_gwas)) {
      merged[allele_match == "none" & toupper(get(ea_eqtl)) == toupper(get(oa_gwas)),
             allele_match := "flipped"]
    }
    # Harmonize GWAS EAF: flip when alleles are swapped
    merged[, eaf_gwas_harmonized := get(eaf_gwas)]
    merged[allele_match == "flipped", eaf_gwas_harmonized := 1 - get(eaf_gwas)]

    n_direct <- sum(merged$allele_match == "direct")
    n_flipped <- sum(merged$allele_match == "flipped")
    n_none <- sum(merged$allele_match == "none")
    cat("  Allele matching: direct =", n_direct, ", flipped =", n_flipped,
        ", unmatched =", n_none, "\n")
  } else {
    cat("  WARNING: Cannot find allele columns for harmonization\n")
    merged[, eaf_gwas_harmonized := get(eaf_gwas)]
    merged[, allele_match := "unknown"]
  }

  # Filter to non-missing EAF and matched alleles
  comp <- merged[allele_match %in% c("direct", "flipped") &
                   !is.na(get(eaf_eqtl)) & !is.na(eaf_gwas_harmonized) &
                   get(eaf_eqtl) > 0 & get(eaf_eqtl) < 1 &
                   eaf_gwas_harmonized > 0 & eaf_gwas_harmonized < 1]
  cat("  Variants with valid harmonized EAF:", nrow(comp), "\n")

  if (nrow(comp) < 50) {
    cat("  WARNING: Too few variants for correlation\n")
    next
  }

  eaf_e <- comp[[eaf_eqtl]]
  eaf_g <- comp[["eaf_gwas_harmonized"]]

  # Compute metrics
  pearson_r <- cor(eaf_e, eaf_g, method = "pearson")
  spearman_rho <- cor(eaf_e, eaf_g, method = "spearman")
  eaf_diff <- abs(eaf_e - eaf_g)
  n_large_diff <- sum(eaf_diff > 0.3)
  pct_large_diff <- 100 * n_large_diff / nrow(comp)
  median_diff <- median(eaf_diff)
  mean_diff <- mean(eaf_diff)

  cat("  Pearson r:", round(pearson_r, 4), "\n")
  cat("  Spearman rho:", round(spearman_rho, 4), "\n")
  cat("  Median |EAF diff|:", round(median_diff, 4), "\n")
  cat("  |EAF diff| > 0.3:", n_large_diff, "(", round(pct_large_diff, 1), "%)\n")

  results[[cfg$name]] <- data.table(
    source = cfg$name,
    ancestry = cfg$ancestry,
    n_variants = nrow(comp),
    pearson_r = round(pearson_r, 4),
    spearman_rho = round(spearman_rho, 4),
    median_eaf_diff = round(median_diff, 4),
    mean_eaf_diff = round(mean_diff, 4),
    n_large_diff_0.3 = n_large_diff,
    pct_large_diff_0.3 = round(pct_large_diff, 2)
  )

  # Plot EAF scatter
  pdf(file.path(OUTDIR, paste0("eaf_scatter_", cfg$name, ".pdf")), width = 6, height = 6)
  plot(eaf_e, eaf_g,
       pch = 16, cex = 0.3, col = rgb(0, 0, 0, 0.1),
       xlab = "EAF (European eQTL, Broadaway)",
       ylab = paste0("EAF (", cfg$ancestry, " GWAS)"),
       main = paste0("EAF Concordance: EUR eQTL vs ", cfg$ancestry, " GWAS\n",
                     "r = ", round(pearson_r, 3),
                     ", N = ", nrow(comp),
                     ", |diff|>0.3: ", n_large_diff))
  abline(0, 1, col = "red", lwd = 2)
  abline(0.3, 1, col = "grey", lty = 2)
  abline(-0.3, 1, col = "grey", lty = 2)
  dev.off()
  cat("  Saved: eaf_scatter_", cfg$name, ".pdf\n")
}

# ==============================================================================
# Save summary
# ==============================================================================
if (length(results) > 0) {
  summary_dt <- rbindlist(results)
  fwrite(summary_dt, file.path(OUTDIR, "eaf_concordance_summary.csv"))
  cat("\n--- Summary ---\n")
  print(summary_dt)
  cat("\nSaved: eaf_concordance_summary.csv\n")
} else {
  cat("\nWARNING: No results to save\n")
}

cat("\nFinished:", format(Sys.time()), "\n")
