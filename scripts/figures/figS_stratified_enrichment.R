#!/usr/bin/env Rscript
# figS_stratified_enrichment.R — Script 214
# Supplementary figure: GWAS-RNA stratified enrichment analysis
# 4 panels: (a) Stratified QQ, (b) Gene-set forest, (c) PIP distribution, (d) DEG enrichment
# Output: figures/supplementary/figS04_coloc/figS_stratified_enrichment.pdf (170mm x 170mm)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

outdir <- FIGS04_DIR
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(FIGS04_DIR, "panels"), showWarnings = FALSE, recursive = TRUE)

datadir <- file.path(BASE, "RNA-seq/results/gwas_rna_integration")

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
qq_data   <- fread(file.path(datadir, "stratified_qq_data.csv"))
gs_enrich <- fread(file.path(datadir, "geneset_enrichment_results.csv"))
deg_strat <- fread(file.path(datadir, "deg_stratum_coloc_summary.csv"))
# Use full credible sets (9.4K variants, full PIP range) instead of high-PIP only
pip_data  <- fread(file.path(BASE, "GWAS/finemapping/results/credible_sets.csv"),
                   select = c("variant_id", "recommended_pip"))

cat("Loaded:", nrow(qq_data), "QQ rows,", nrow(gs_enrich), "gene sets,",
    nrow(deg_strat), "DEG strata,", nrow(pip_data), "credible set variants\n")

# ---------------------------------------------------------------------------
# Panel A: Stratified QQ-like plot — cumulative COLOC PP.H4 by DEG stratum
# ---------------------------------------------------------------------------
# Strata colors
strata_colors <- c(
  "Bottom 50%"   = masld_colors$down,
  "Middle"       = "#7B1FA2",
  "Top 10% DEGs" = masld_colors$up
)

# Order strata factor
qq_data[, stratum := factor(stratum, levels = c("Top 10% DEGs", "Middle", "Bottom 50%"))]

pA <- ggplot(qq_data, aes(x = expected_rank, y = observed_coloc, color = stratum)) +
  geom_line(linewidth = 0.4, alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50",
              linewidth = 0.3) +
  scale_color_manual(values = strata_colors, name = "DEG stratum") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Expected rank (uniform quantile)",
       y = "COLOC PP.H4",
       title = "Stratified QQ: COLOC by DEG strength") +
  theme_masld(base_size = 7) +
  theme(legend.position = c(0.25, 0.85),
        legend.background = element_blank(),
        legend.key.size = unit(0.2, "cm"),
        legend.text = element_text(size = 5.5))

# ---------------------------------------------------------------------------
# Panel B: Gene-set enrichment forest plot (17 gene sets)
# ---------------------------------------------------------------------------
# Compute 95% CI for log(OR) then exponentiate
gs_enrich[, log_or := log(fisher_or)]
gs_enrich[, log_or_se := abs(log_or / qnorm(fisher_p / 2, lower.tail = FALSE))]
# Handle p=0 or very small p: cap SE
gs_enrich[is.infinite(log_or_se) | is.nan(log_or_se), log_or_se := 0.5]
gs_enrich[, or_lo := exp(log_or - 1.96 * log_or_se)]
gs_enrich[, or_hi := exp(log_or + 1.96 * log_or_se)]
gs_enrich[, sig := fdr_coloc < 0.05]

# Clean gene-set names for display
gs_enrich[, gene_set_label := gsub("_", " ", gene_set)]
gs_enrich[, gene_set_label := gsub("sc ", "sc: ", gene_set_label)]

# Order by fold enrichment
gs_enrich[, gene_set_label := factor(gene_set_label,
                                      levels = gene_set_label[order(fold_enrichment)])]

pB <- ggplot(gs_enrich, aes(x = fisher_or, y = gene_set_label)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray60", linewidth = 0.3) +
  geom_errorbarh(aes(xmin = or_lo, xmax = or_hi), height = 0.3,
                 color = "gray40", linewidth = 0.3) +
  geom_point(aes(fill = sig), shape = 21, size = 1.5, stroke = 0.3) +
  scale_fill_manual(values = c("TRUE" = masld_colors$up, "FALSE" = masld_colors$ns),
                    labels = c("TRUE" = "FDR < 0.05", "FALSE" = "NS"),
                    name = "Significance") +
  scale_x_log10() +
  labs(x = "Fisher OR (COLOC PP.H4 > 0.5)", y = NULL,
       title = "Gene-set enrichment for COLOC") +
  theme_masld(base_size = 7) +
  theme(legend.position = c(0.8, 0.15),
        legend.background = element_blank(),
        legend.key.size = unit(0.2, "cm"),
        axis.text.y = element_text(size = 5.5))

# ---------------------------------------------------------------------------
# Panel C: PIP distribution — deduplicated unique variants (max PIP across studies)
# ---------------------------------------------------------------------------
pip_valid <- pip_data[!is.na(recommended_pip)]
# Deduplicate: one row per variant, keeping max PIP across GWAS studies
pip_uniq <- pip_valid[, .(recommended_pip = max(recommended_pip),
                          n_studies = .N), by = variant_id]

cat("Panel C: ", nrow(pip_valid), " CS entries -> ", nrow(pip_uniq), " unique variants\n")

# Granular bins across full 0-1 range
pip_uniq[, pip_bin := cut(recommended_pip,
                          breaks = c(0, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.5, 0.7, 0.9, 0.95, 1.0),
                          include.lowest = TRUE)]

pip_bin_counts <- pip_uniq[!is.na(pip_bin), .N, by = pip_bin][order(pip_bin)]

# Color gradient: low PIP = grey, high PIP = magenta
pip_bin_counts[, midpoint := c(0.0005, 0.003, 0.0075, 0.03, 0.075, 0.15, 0.35, 0.6, 0.8, 0.925, 0.975)]

pC <- ggplot(pip_bin_counts, aes(x = pip_bin, y = N, fill = midpoint)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = format(N, big.mark = ",")),
            vjust = -0.3, size = 1.8, color = "black") +
  scale_fill_gradient(low = masld_colors$ns, high = masld_colors$up,
                      name = "PIP", guide = "none") +
  scale_y_log10(expand = expansion(mult = c(0, 0.15)),
                labels = scales::comma) +
  labs(x = "Recommended PIP (max across studies)", y = "Unique variants (log scale)",
       title = sprintf("Fine-mapped variant PIP (%s unique variants in 95%% CS)",
                       format(nrow(pip_uniq), big.mark = ","))) +
  theme_masld(base_size = 7) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5.5))

# ---------------------------------------------------------------------------
# Panel D: DEG enrichment — observed vs expected COLOC overlap
# ---------------------------------------------------------------------------
# Use the All_DEGs row from geneset enrichment
all_degs_row <- gs_enrich[gene_set == "All_DEGs"]

# Build bar data: observed fold enrichment vs 1.0 expected
bar_dt <- data.table(
  label = c("Expected (null)", "Observed DEGs"),
  fold  = c(1.0, all_degs_row$fold_enrichment),
  type  = c("Expected", "Observed")
)
bar_dt[, label := factor(label, levels = c("Expected (null)", "Observed DEGs"))]

pval_label <- sprintf("Wilcoxon P = %s\nFisher OR = %.2f",
                       format(all_degs_row$wilcox_p_coloc, digits = 3, scientific = TRUE),
                       all_degs_row$fisher_or)

pD <- ggplot(bar_dt, aes(x = label, y = fold, fill = type)) +
  geom_col(width = 0.5) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray60", linewidth = 0.3) +
  annotate("text", x = 1.5, y = max(bar_dt$fold) * 1.1,
           label = pval_label, size = 2, color = "gray20", hjust = 0.5) +
  scale_fill_manual(values = c("Expected" = masld_colors$ns,
                                "Observed" = masld_colors$up),
                    guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(x = NULL, y = "Fold enrichment (COLOC in DEGs vs non-DEGs)",
       title = "DEG enrichment for COLOC signals") +
  theme_masld(base_size = 7)

# ---------------------------------------------------------------------------
# Assemble
# ---------------------------------------------------------------------------
combined <- (pA | pB) / (pC | pD) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "plain"))

outfile <- file.path(outdir, "figS_stratified_enrichment.pdf")
save_fig(combined, outfile, width = 170 / 25.4, height = 170 / 25.4)

cat("Saved:", outfile, "\n")
cat("Done.\n")
