#!/usr/bin/env Rscript
# 63b_otters_concordance.R
# ---------------------------------------------------------------------------
# Compare OTTERS-Broadaway TWAS vs S-PrediXcan-GTEx TWAS
#
# Input:
#   - RNA-seq/results/causal_inference/otters_broadaway/{gwas}/otters_twas_combined.csv
#   - RNA-seq/results/causal_inference/{gwas}/twas_spredixcan_liver.csv
#
# Output:
#   - RNA-seq/results/causal_inference/otters_broadaway/concordance/
#       - otters_vs_spredixcan_concordance.csv
#       - concordance_summary.txt
#       - scatter_z.pdf
#       - venn_significant.pdf
#
# Usage:
#   Rscript RNA-seq/63b_otters_concordance.R
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
setwd(BASE)

OUT_DIR <- file.path(BASE, "RNA-seq/results/causal_inference/otters_broadaway/concordance")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# GWAS to compare
GWAS_LIST <- c("ghodsian", "ukbb_alt", "ukbb_ast", "ukbb_ggt")

cat("=== OTTERS vs S-PrediXcan Concordance Analysis ===\n")
cat("Start:", format(Sys.time()), "\n\n")

# Load gene annotation for symbol mapping
anno <- fread("data/broadaway_eqtl/otters_format/gene_anno.txt")
gene_map <- setNames(anno$GeneName, anno$TargetID)

all_concordance <- list()

for (gwas_name in GWAS_LIST) {

  cat("--- GWAS:", gwas_name, "---\n")

  # Load S-PrediXcan results (GTEx)
  spredixcan_file <- file.path(BASE, "RNA-seq/results/causal_inference",
                                gwas_name, "twas_spredixcan_liver.csv")
  if (!file.exists(spredixcan_file)) {
    cat("  S-PrediXcan file not found, skipping\n")
    next
  }
  spx <- fread(spredixcan_file)
  setnames(spx, c("gene", "zscore", "pvalue", "fdr"),
           c("ensembl_id", "spx_z", "spx_pval", "spx_fdr"),
           skip_absent = TRUE)

  # Load OTTERS results
  otters_file <- file.path(BASE, "RNA-seq/results/causal_inference/otters_broadaway",
                            gwas_name, "otters_twas_combined.csv")
  if (!file.exists(otters_file)) {
    cat("  OTTERS file not found, skipping\n")
    next
  }
  ott <- fread(otters_file)

  # Standardize OTTERS column names
  # OTTERS testing.py outputs: CHROM GeneStart GeneEnd TargetID n_snps FUSION_Z FUSION_PVAL
  if ("TargetID" %in% names(ott) && !"ensembl_id" %in% names(ott)) {
    setnames(ott, "TargetID", "ensembl_id")
  }
  if ("FUSION_Z" %in% names(ott) && !"otters_z" %in% names(ott)) {
    setnames(ott, "FUSION_Z", "otters_z")
  }
  if ("FUSION_PVAL" %in% names(ott) && !"otters_pval" %in% names(ott)) {
    setnames(ott, "FUSION_PVAL", "otters_pval")
  }
  # Handle alternate column names from our Python combiner (63c)
  z_candidates <- c("Z", "best_z", "Zscore", "TWAS_Z")
  p_candidates <- c("P", "acat_pval", "best_pval", "pvalue")
  for (zc in z_candidates) {
    if (zc %in% names(ott) && !"otters_z" %in% names(ott)) {
      setnames(ott, zc, "otters_z"); break
    }
  }
  for (pc in p_candidates) {
    if (pc %in% names(ott) && !"otters_pval" %in% names(ott)) {
      setnames(ott, pc, "otters_pval"); break
    }
  }

  cat("  S-PrediXcan genes:", nrow(spx), "\n")
  cat("  OTTERS genes:", nrow(ott), "\n")

  # Merge on Ensembl ID
  if (!"ensembl_id" %in% names(ott)) {
    cat("  WARNING: No ensembl_id column in OTTERS output. Columns:", paste(names(ott), collapse=", "), "\n")
    next
  }

  # Deduplicate OTTERS (keep best p-value per gene if multiple models)
  if ("otters_pval" %in% names(ott)) {
    ott <- ott[order(otters_pval)][!duplicated(ensembl_id)]
  }

  merged <- merge(
    spx[, .(ensembl_id, spx_z, spx_pval, spx_fdr)],
    ott[, intersect(names(ott), c("ensembl_id", "otters_z", "otters_pval", "n_snps", "fdr")),
        with = FALSE],
    by = "ensembl_id"
  )

  if ("fdr" %in% names(merged) && !"otters_fdr" %in% names(merged)) {
    setnames(merged, "fdr", "otters_fdr")
  }

  cat("  Overlapping genes:", nrow(merged), "\n")

  if (nrow(merged) < 10) {
    cat("  Too few overlapping genes for concordance analysis\n")
    next
  }

  # Filter to valid z-scores
  valid <- merged[is.finite(spx_z) & is.finite(otters_z)]
  cat("  Valid z-scores:", nrow(valid), "\n")

  # Concordance metrics
  spearman_rho <- cor(valid$spx_z, valid$otters_z, method = "spearman")
  pearson_r    <- cor(valid$spx_z, valid$otters_z, method = "pearson")

  # Direction concordance
  dir_concord <- mean(sign(valid$spx_z) == sign(valid$otters_z), na.rm = TRUE)

  # Significance overlap
  spx_sig <- valid$spx_pval < 0.05
  ott_sig <- valid$otters_pval < 0.05
  both_sig <- spx_sig & ott_sig
  either_sig <- spx_sig | ott_sig
  jaccard <- sum(both_sig) / max(sum(either_sig), 1)

  # FDR overlap
  spx_fdr_sig <- if ("spx_fdr" %in% names(valid)) valid$spx_fdr < 0.05 else rep(FALSE, nrow(valid))
  ott_fdr_sig <- if ("otters_fdr" %in% names(valid)) valid$otters_fdr < 0.05 else rep(FALSE, nrow(valid))

  cat(sprintf("  Spearman rho: %.3f\n", spearman_rho))
  cat(sprintf("  Pearson r: %.3f\n", pearson_r))
  cat(sprintf("  Direction concordance: %.1f%%\n", dir_concord * 100))
  cat(sprintf("  S-PrediXcan p<0.05: %d | OTTERS p<0.05: %d | Both: %d | Jaccard: %.3f\n",
              sum(spx_sig), sum(ott_sig), sum(both_sig), jaccard))
  cat(sprintf("  S-PrediXcan FDR<0.05: %d | OTTERS FDR<0.05: %d\n",
              sum(spx_fdr_sig), sum(ott_fdr_sig)))

  # OTTERS-only and S-PrediXcan-only significant genes
  otters_only_genes <- nrow(ott) - nrow(merged)
  cat(sprintf("  OTTERS-only genes (not in S-PrediXcan): %d\n", otters_only_genes))

  # Save per-GWAS concordance
  valid$gwas <- gwas_name
  valid$gene_symbol <- gene_map[valid$ensembl_id]
  all_concordance[[gwas_name]] <- valid

  # Summary row
  summary_row <- data.table(
    gwas = gwas_name,
    n_spredixcan = nrow(spx),
    n_otters = nrow(ott),
    n_overlap = nrow(valid),
    n_otters_only = otters_only_genes,
    spearman_rho = spearman_rho,
    pearson_r = pearson_r,
    direction_concordance = dir_concord,
    spx_p05 = sum(spx_sig),
    otters_p05 = sum(ott_sig),
    both_p05 = sum(both_sig),
    jaccard_p05 = jaccard,
    spx_fdr05 = sum(spx_fdr_sig),
    otters_fdr05 = sum(ott_fdr_sig)
  )

  if (!exists("summary_dt")) {
    summary_dt <- summary_row
  } else {
    summary_dt <- rbind(summary_dt, summary_row)
  }

  # ---- Scatter plot ----
  p <- ggplot(valid, aes(x = spx_z, y = otters_z)) +
    geom_point(alpha = 0.3, size = 0.8, color = "grey40") +
    geom_point(data = valid[spx_sig | ott_sig],
               alpha = 0.6, size = 1.2, color = "steelblue") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.5) +
    labs(
      title = paste0("OTTERS-Broadaway vs S-PrediXcan-GTEx: ", toupper(gwas_name)),
      subtitle = sprintf("N=%d genes | rho=%.3f | direction=%.1f%%",
                          nrow(valid), spearman_rho, dir_concord * 100),
      x = "S-PrediXcan z-score (GTEx Liver, N=208)",
      y = "OTTERS z-score (Broadaway Liver, N=1,183)"
    ) +
    theme_minimal(base_size = 11) +
    coord_equal()

  ggsave(file.path(OUT_DIR, paste0("scatter_z_", gwas_name, ".pdf")),
         p, width = 6, height = 6)
  cat("  Saved scatter plot\n")
}

# ---- Save combined results ----
if (length(all_concordance) > 0) {
  combined <- rbindlist(all_concordance, fill = TRUE)
  fwrite(combined, file.path(OUT_DIR, "otters_vs_spredixcan_concordance.csv"))
  cat("\nSaved concordance data:", nrow(combined), "rows\n")
}

if (exists("summary_dt")) {
  fwrite(summary_dt, file.path(OUT_DIR, "concordance_summary.csv"))
  cat("\n=== Concordance Summary ===\n")
  print(summary_dt)

  # Write text summary
  sink(file.path(OUT_DIR, "concordance_summary.txt"))
  cat("OTTERS-Broadaway vs S-PrediXcan-GTEx Concordance\n")
  cat("================================================\n\n")
  for (i in seq_len(nrow(summary_dt))) {
    row <- summary_dt[i]
    cat(sprintf("GWAS: %s\n", row$gwas))
    cat(sprintf("  S-PrediXcan genes: %d (GTEx Liver, N=208)\n", row$n_spredixcan))
    cat(sprintf("  OTTERS genes: %d (Broadaway Liver, N=1,183)\n", row$n_otters))
    cat(sprintf("  Overlapping genes: %d\n", row$n_overlap))
    cat(sprintf("  OTTERS-only genes: %d\n", row$n_otters_only))
    cat(sprintf("  Spearman rho: %.3f\n", row$spearman_rho))
    cat(sprintf("  Direction concordance: %.1f%%\n", row$direction_concordance * 100))
    cat(sprintf("  S-PrediXcan p<0.05: %d | OTTERS p<0.05: %d | Both: %d\n",
                row$spx_p05, row$otters_p05, row$both_p05))
    cat(sprintf("  Jaccard (p<0.05): %.3f\n", row$jaccard_p05))
    cat(sprintf("  S-PrediXcan FDR<0.05: %d | OTTERS FDR<0.05: %d\n",
                row$spx_fdr05, row$otters_fdr05))
    cat("\n")
  }
  sink()
  cat("\nSaved summary to concordance_summary.txt\n")
}

cat("\n=== Done:", format(Sys.time()), "===\n")
