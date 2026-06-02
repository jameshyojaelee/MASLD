#!/usr/bin/env Rscript
# 36b_pseudobulk_bulk_concordance.R
# ---------------------------------------------------------------------------
# Baseline concordance: hepatocyte pseudobulk DE vs bulk dream mega-analysis
#
# Measures current pseudobulk-vs-bulk correlation BEFORE any pipeline changes,
# so subsequent improvements (adding GSE136103, re-running DE) can be quantified.
#
# Metrics:
#   - Spearman rho + Pearson r of LFC (all shared genes)
#   - Direction concordance (% same-sign LFC among co-significant genes)
#   - Jaccard index (both padj < 0.1)
#   - Conserved stratification (rho within CC genes)
#   - Hepatocyte-intrinsic enrichment (Fisher's test)
#
# Inputs:
#   - Analysis/SingleCell/results_gpu_v2/pseudobulk_de/Hepatocytes_de.csv
#   - RNA-seq/Human/.../results/integration/dream_results.csv
#   - RNA-seq/Human/.../results/gene_annotation/human_ensg_to_symbol.tsv
#   - RNA-seq/results/multi_evidence/multi_evidence_atlas.csv (Conserved)
#   - RNA-seq/results/causal_inference/deconv_attribution_scores.csv (optional)
#
# Outputs (RNA-seq/results/pseudobulk/):
#   - pseudobulk_vs_bulk_concordance.csv     (per-gene merged table)
#   - pseudobulk_bulk_summary.csv            (aggregate metrics)
#   - pseudobulk_bulk_concordance.pdf        (scatter plot)
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

cat("=== Script 36b: Pseudobulk-Bulk Concordance (Baseline) ===\n")
cat("Start time:", format(Sys.time()), "\n\n")

# ==============================================================================
# Configuration
# ==============================================================================
BASE_DIR <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

RESULTS_DIR <- file.path(BASE_DIR, "RNA-seq/results/pseudobulk")
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(RESULTS_DIR, "figures"), recursive = TRUE, showWarnings = FALSE)

PB_FILE <- file.path(BASE_DIR,
  "Analysis/SingleCell/results_gpu_v2/pseudobulk_de/Hepatocytes_de.csv")
DREAM_FILE <- file.path(BASE_DIR,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/dream_results.csv")
GENE_CACHE <- file.path(BASE_DIR,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/gene_annotation/human_ensg_to_symbol.tsv")
ATLAS_FILE <- file.path(BASE_DIR,
  "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
ATTRIB_FILE <- file.path(BASE_DIR,
  "RNA-seq/results/causal_inference/deconv_attribution_scores.csv")

# Thresholds
PADJ_THRESH <- 0.1

# ==============================================================================
# 1. Load data
# ==============================================================================
cat("--- Step 1: Loading data ---\n")

stopifnot(file.exists(PB_FILE))
stopifnot(file.exists(DREAM_FILE))
stopifnot(file.exists(GENE_CACHE))

pb <- fread(PB_FILE)
cat("  Pseudobulk hepatocyte DE:", nrow(pb), "genes\n")
cat("  Contrast:", unique(pb$contrast), "\n")

dream <- fread(DREAM_FILE)
cat("  Bulk dream results:", nrow(dream), "genes\n")

ann <- fread(GENE_CACHE)
cat("  Gene annotation:", nrow(ann), "entries\n")

# ==============================================================================
# 2. Gene ID harmonization
# ==============================================================================
cat("\n--- Step 2: Gene ID harmonization ---\n")

# Dream has versioned ENSG IDs (e.g., ENSG00000310526.1)
# Strip version to get base ENSG
dream[, ensembl_base := sub("\\.\\d+$", "", gene)]

# Pseudobulk has a mix of ENSG IDs (no version) and gene symbols
pb[, is_ensg := grepl("^ENSG\\d+$", gene)]
n_ensg <- sum(pb$is_ensg)
n_symbol <- sum(!pb$is_ensg)
cat("  Pseudobulk gene IDs: ", n_ensg, " ENSG + ", n_symbol, " symbols\n")

# Build annotation lookup: symbol -> ensembl_base (deduplicated, prefer protein_coding)
ann_dedup <- ann[symbol != "" & !is.na(symbol)]
ann_dedup[, priority := ifelse(gene_type == "protein_coding", 1L, 2L)]
setorder(ann_dedup, symbol, priority)
ann_dedup <- ann_dedup[!duplicated(symbol)]
symbol_to_ensg <- ann_dedup[, .(symbol, ensembl_base = gene_base)]

# Map pseudobulk to ensembl_base
# Route 1: ENSG IDs -> use directly
pb_ensg <- pb[is_ensg == TRUE, .(ensembl_base = gene, pb_logFC = logFC,
                                   pb_padj = padj, pb_tstat = t_stat)]

# Route 2: Symbols -> map via annotation
pb_sym <- merge(
  pb[is_ensg == FALSE, .(symbol = gene, pb_logFC = logFC, pb_padj = padj, pb_tstat = t_stat)],
  symbol_to_ensg,
  by = "symbol"
)[, .(ensembl_base, pb_logFC, pb_padj, pb_tstat)]

# Combine
pb_mapped <- rbind(pb_ensg, pb_sym)
cat("  Mapped pseudobulk genes:", nrow(pb_mapped), "(",
    nrow(pb_ensg), "ENSG +", nrow(pb_sym), "via symbol)\n")

# Deduplicate (keep first occurrence if any ENSG collision)
pb_mapped <- pb_mapped[!duplicated(ensembl_base)]

# ==============================================================================
# 3. Merge pseudobulk with dream
# ==============================================================================
cat("\n--- Step 3: Merging pseudobulk x dream ---\n")

concordance <- merge(
  pb_mapped,
  dream[, .(ensembl_base, bulk_logFC = logFC, bulk_padj = padj, bulk_tstat = t)],
  by = "ensembl_base"
)
cat("  Shared genes:", nrow(concordance), "\n")

if (nrow(concordance) == 0) {
  cat("  ERROR: No shared genes. Check ID mapping.\n")
  quit(save = "no", status = 1)
}

# Add gene symbols via annotation
concordance <- merge(concordance,
  ann_dedup[, .(ensembl_base = gene_base, symbol)],
  by = "ensembl_base", all.x = TRUE)

# ==============================================================================
# 4. Compute concordance metrics
# ==============================================================================
cat("\n--- Step 4: Concordance metrics ---\n")

# Full-genome LFC correlation
r_all <- cor(concordance$pb_logFC, concordance$bulk_logFC,
             use = "complete.obs", method = "pearson")
rho_all <- cor(concordance$pb_logFC, concordance$bulk_logFC,
               use = "complete.obs", method = "spearman")
cat("  All shared genes (n=", nrow(concordance), "):\n")
cat("    Pearson r:  ", round(r_all, 4), "\n")
cat("    Spearman rho:", round(rho_all, 4), "\n")

# Significance overlap
pb_sig <- concordance[pb_padj < PADJ_THRESH]
bulk_sig <- concordance[bulk_padj < PADJ_THRESH]
both_sig <- concordance[pb_padj < PADJ_THRESH & bulk_padj < PADJ_THRESH]

n_pb_sig <- nrow(pb_sig)
n_bulk_sig <- nrow(bulk_sig)
n_both_sig <- nrow(both_sig)

jaccard <- n_both_sig / (n_pb_sig + n_bulk_sig - n_both_sig)
cat("\n  Significance (padj < ", PADJ_THRESH, "):\n")
cat("    Pseudobulk DEGs:", n_pb_sig, "\n")
cat("    Bulk dream DEGs:", n_bulk_sig, "\n")
cat("    Overlap:", n_both_sig, "\n")
cat("    Jaccard:", round(jaccard, 4), "\n")

# Direction concordance among co-significant genes
if (n_both_sig > 0) {
  dir_conc <- mean(sign(both_sig$pb_logFC) == sign(both_sig$bulk_logFC))
  cat("    Direction concordance (co-sig):", round(dir_conc * 100, 1), "%\n")
} else {
  dir_conc <- NA_real_
  cat("    Direction concordance: N/A (no co-significant genes)\n")
}

# Correlation among significant genes only
if (n_both_sig >= 10) {
  r_sig <- cor(both_sig$pb_logFC, both_sig$bulk_logFC,
               use = "complete.obs", method = "pearson")
  rho_sig <- cor(both_sig$pb_logFC, both_sig$bulk_logFC,
                 use = "complete.obs", method = "spearman")
  cat("    LFC Pearson r (co-sig):", round(r_sig, 4), "\n")
  cat("    LFC Spearman rho (co-sig):", round(rho_sig, 4), "\n")
} else {
  r_sig <- rho_sig <- NA_real_
}

# ==============================================================================
# 5. Conserved stratification
# ==============================================================================
cat("\n--- Step 5: Conserved stratification ---\n")

rho_cc <- r_cc <- n_cc <- NA_real_
if (file.exists(ATLAS_FILE)) {
  atlas <- fread(ATLAS_FILE, select = c("ensembl_id", "human_symbol", "is_conserved"))
  cc_genes <- atlas[is_conserved == TRUE, ensembl_id]
  cat("  Conserved genes in atlas:", length(cc_genes), "\n")

  concordance[, is_cc := ensembl_base %in% cc_genes]
  cc_data <- concordance[is_cc == TRUE]
  n_cc <- nrow(cc_data)
  cat("  Conserved in concordance table:", n_cc, "\n")

  if (n_cc >= 10) {
    r_cc <- cor(cc_data$pb_logFC, cc_data$bulk_logFC,
                use = "complete.obs", method = "pearson")
    rho_cc <- cor(cc_data$pb_logFC, cc_data$bulk_logFC,
                  use = "complete.obs", method = "spearman")
    cat("    Pearson r (CC):", round(r_cc, 4), "\n")
    cat("    Spearman rho (CC):", round(rho_cc, 4), "\n")

    # Direction concordance among CC genes that are sig in both
    cc_both <- cc_data[pb_padj < PADJ_THRESH & bulk_padj < PADJ_THRESH]
    if (nrow(cc_both) > 0) {
      cc_dir <- mean(sign(cc_both$pb_logFC) == sign(cc_both$bulk_logFC))
      cat("    Direction concordance (CC co-sig, n=",
          nrow(cc_both), "):", round(cc_dir * 100, 1), "%\n")
    }
  }
} else {
  cat("  Atlas file not found, skipping CC stratification.\n")
  concordance[, is_cc := FALSE]
}

# ==============================================================================
# 6. Hepatocyte-intrinsic enrichment (optional)
# ==============================================================================
cat("\n--- Step 6: Hepatocyte-intrinsic enrichment ---\n")

fisher_or <- fisher_p <- NA_real_
if (file.exists(ATTRIB_FILE)) {
  attrib <- fread(ATTRIB_FILE)
  if ("attribution_class" %in% names(attrib)) {
    hep_intrinsic <- attrib[attribution_class == "Hepatocyte_intrinsic", gene]
    cat("  Hepatocyte-intrinsic genes:", length(hep_intrinsic), "\n")

    # Map to ensembl_base for matching
    hep_ensg <- hep_intrinsic[grepl("^ENSG", hep_intrinsic)]
    hep_sym_mapped <- symbol_to_ensg[symbol %in% hep_intrinsic, ensembl_base]
    hep_all_ensg <- unique(c(hep_ensg, hep_sym_mapped))

    pb_sig_genes <- concordance[pb_padj < PADJ_THRESH, ensembl_base]
    all_tested <- concordance$ensembl_base

    a <- sum(pb_sig_genes %in% hep_all_ensg)
    b <- length(pb_sig_genes) - a
    c <- sum(all_tested %in% hep_all_ensg) - a
    d <- length(all_tested) - a - b - c
    if (all(c(a, b, c, d) >= 0)) {
      ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2))
      fisher_or <- as.numeric(ft$estimate)
      fisher_p <- ft$p.value
      cat("  Pseudobulk DEGs in hepatocyte-intrinsic:", a, "/", length(pb_sig_genes), "\n")
      cat("  Fisher OR:", round(fisher_or, 2), ", p =", format(fisher_p, digits = 3), "\n")
    }
  }
} else {
  cat("  Attribution file not found, skipping.\n")
}

# ==============================================================================
# 7. Save per-gene concordance table
# ==============================================================================
cat("\n--- Step 7: Saving results ---\n")

# Add significance flags
concordance[, pb_sig := pb_padj < PADJ_THRESH]
concordance[, bulk_sig := bulk_padj < PADJ_THRESH]
concordance[, both_sig := pb_sig & bulk_sig]
concordance[, direction_concordant := sign(pb_logFC) == sign(bulk_logFC)]

fwrite(concordance, file.path(RESULTS_DIR, "pseudobulk_vs_bulk_concordance.csv"))
cat("  Saved: pseudobulk_vs_bulk_concordance.csv (", nrow(concordance), " genes)\n")

# Summary metrics table
summary_dt <- data.table(
  metric = c("n_shared_genes", "pearson_r_all", "spearman_rho_all",
             "n_pb_sig", "n_bulk_sig", "n_both_sig", "jaccard",
             "direction_concordance_cosig",
             "pearson_r_cosig", "spearman_rho_cosig",
             "n_conserved", "pearson_r_cc", "spearman_rho_cc",
             "fisher_or_hep_intrinsic", "fisher_p_hep_intrinsic",
             "padj_threshold", "timestamp"),
  value = c(nrow(concordance), r_all, rho_all,
            n_pb_sig, n_bulk_sig, n_both_sig, jaccard,
            dir_conc,
            r_sig, rho_sig,
            n_cc, r_cc, rho_cc,
            fisher_or, fisher_p,
            PADJ_THRESH, as.character(Sys.time()))
)
fwrite(summary_dt, file.path(RESULTS_DIR, "pseudobulk_bulk_summary.csv"))
cat("  Saved: pseudobulk_bulk_summary.csv\n")

# ==============================================================================
# 8. Scatter plot
# ==============================================================================
cat("\n--- Step 8: Scatter plot ---\n")

concordance[, sig_class := fifelse(
  pb_sig & bulk_sig, "Both significant",
  fifelse(pb_sig, "Pseudobulk only",
  fifelse(bulk_sig, "Bulk only", "Neither"))
)]

# Color palette
sig_colors <- c(
  "Both significant" = "#D7191C",
  "Pseudobulk only" = "#2C7BB6",
  "Bulk only" = "#FDAE61",
  "Neither" = "grey80"
)

p <- ggplot(concordance, aes(x = bulk_logFC, y = pb_logFC)) +
  geom_point(data = concordance[sig_class == "Neither"],
             aes(color = sig_class), alpha = 0.1, size = 0.3) +
  geom_point(data = concordance[sig_class == "Bulk only"],
             aes(color = sig_class), alpha = 0.3, size = 0.5) +
  geom_point(data = concordance[sig_class == "Pseudobulk only"],
             aes(color = sig_class), alpha = 0.4, size = 0.5) +
  geom_point(data = concordance[sig_class == "Both significant"],
             aes(color = sig_class), alpha = 0.6, size = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.3) +
  scale_color_manual(values = sig_colors) +
  labs(
    title = "Hepatocyte Pseudobulk vs Bulk Dream DE Concordance",
    subtitle = paste0(
      "n = ", format(nrow(concordance), big.mark = ","),
      " shared genes | r = ", round(r_all, 3),
      " | rho = ", round(rho_all, 3),
      ifelse(!is.na(dir_conc),
             paste0(" | dir.conc = ", round(dir_conc * 100, 1), "%"),
             "")
    ),
    x = "Bulk Dream log2FC (10-cohort mega-analysis)",
    y = "Hepatocyte Pseudobulk log2FC (scVI integration)",
    color = paste0("Significance (padj < ", PADJ_THRESH, ")")
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, color = "grey30"),
    legend.position = "bottom",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8)
  ) +
  coord_cartesian(xlim = quantile(concordance$bulk_logFC, c(0.001, 0.999)),
                  ylim = quantile(concordance$pb_logFC, c(0.001, 0.999)))

pdf(file.path(RESULTS_DIR, "figures/pseudobulk_bulk_concordance.pdf"),
    width = 7, height = 6.5)
print(p)
dev.off()
cat("  Saved: figures/pseudobulk_bulk_concordance.pdf\n")

# Conserved overlay plot (if available)
if (any(concordance$is_cc)) {
  p_cc <- ggplot(concordance, aes(x = bulk_logFC, y = pb_logFC)) +
    geom_point(alpha = 0.05, size = 0.3, color = "grey80") +
    geom_point(data = concordance[is_cc == TRUE],
               aes(color = "Conserved"), alpha = 0.6, size = 1.0) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
    scale_color_manual(values = c("Conserved" = "#D7191C")) +
    labs(
      title = "Conserved Genes: Pseudobulk vs Bulk",
      subtitle = paste0(
        "CC genes: n = ", sum(concordance$is_cc),
        " | r = ", round(r_cc, 3),
        " | rho = ", round(rho_cc, 3)
      ),
      x = "Bulk Dream log2FC",
      y = "Pseudobulk log2FC",
      color = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom") +
    coord_cartesian(xlim = quantile(concordance$bulk_logFC, c(0.001, 0.999)),
                    ylim = quantile(concordance$pb_logFC, c(0.001, 0.999)))

  pdf(file.path(RESULTS_DIR, "figures/pseudobulk_bulk_cc_overlay.pdf"),
      width = 7, height = 6.5)
  print(p_cc)
  dev.off()
  cat("  Saved: figures/pseudobulk_bulk_cc_overlay.pdf\n")
}

# ==============================================================================
# Summary
# ==============================================================================
cat("\n=== Script 36b: Concordance Summary ===\n")
cat("  Shared genes:     ", nrow(concordance), "\n")
cat("  Spearman rho:     ", round(rho_all, 4), "\n")
cat("  Pearson r:        ", round(r_all, 4), "\n")
cat("  Jaccard:          ", round(jaccard, 4), "\n")
cat("  Dir. concordance: ", ifelse(!is.na(dir_conc), paste0(round(dir_conc*100,1), "%"), "N/A"), "\n")
if (!is.na(rho_cc)) {
  cat("  CC Spearman rho:  ", round(rho_cc, 4), "\n")
}
cat("\n  SUCCESS CRITERION: rho > 0.3 =>",
    ifelse(!is.na(rho_all) && rho_all > 0.3, "PASSED", "BELOW THRESHOLD"), "\n")
cat("\nEnd time:", format(Sys.time()), "\n")
