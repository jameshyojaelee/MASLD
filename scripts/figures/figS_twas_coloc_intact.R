#!/usr/bin/env Rscript
# figS_twas_coloc_intact.R — Supplementary: TWAS × COLOC concordance + INTACT validation
#
# 4 panels:
# (a) TWAS z-score vs COLOC PP.H4 scatter — shows which genes have both vs one signal
# (b) INTACT Bayesian filtering — raw TWAS ranks vs INTACT-filtered ranks
# (c) Cell-type proportion × COLOC gene expression (hepatocytes)
# (d) Cell-type proportion × COLOC gene expression (macrophages)
#
# Output: figures/supplementary/figS04_coloc/figS_twas_coloc_intact.pdf

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

outdir <- FIGS04_DIR
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(FIGS04_DIR, "panels"), showWarnings = FALSE, recursive = TRUE)

# ===========================================================================
# Load data
# ===========================================================================
cat("Loading data...\n")

# COLOC gene-level
coloc <- fread(file.path(BASE,
  "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv"))
coloc <- coloc[gene != "" & !is.na(gene)]

# TWAS
twas <- fread(file.path(BASE,
  "RNA-seq/results/causal_inference/twas_multi_gwas_combined.csv"))
twas_best <- twas[, .(twas_z = zscore[which.min(pvalue)],
                       twas_p = min(pvalue)),
                  by = .(gene = gene_name)]

# INTACT
intact <- fread(file.path(BASE,
  "RNA-seq/results/gwas_rna_integration/intact_scores.csv"))

# Proportion × COLOC
prop <- fread(file.path(BASE,
  "RNA-seq/results/gwas_rna_integration/proportion_coloc_correlation.csv"))

# Merge TWAS + COLOC
tc <- merge(
  coloc[, .(gene, coloc_pp4 = coloc_best_pp4, coloc_n = coloc_n_gwas_h4_05)],
  twas_best,
  by = "gene"
)
cat("  TWAS × COLOC:", nrow(tc), "genes\n")

# Add INTACT
tc <- merge(tc, intact[, .(gene, multi_intact_score, intact_score_ct)],
            by = "gene", all.x = TRUE)

# ===========================================================================
# Panel A: TWAS z-score vs COLOC PP.H4 scatter
# ===========================================================================
cat("Panel A: TWAS vs COLOC scatter...\n")

# Classify genes
tc[, category := fcase(
  coloc_pp4 > 0.9 & abs(twas_z) > 4, "Both COLOC + TWAS",
  coloc_pp4 > 0.9, "COLOC only",
  abs(twas_z) > 4 & twas_p < 0.05, "TWAS only",
  default = "Neither"
)]

cat("  Categories:\n")
print(tc[, .N, by = category][order(-N)])

# Labels for key genes
priority <- c("PNPLA3", "TM6SF2", "HSD17B13", "MBOAT7", "GCKR",
              "THRB", "NR1H4", "PPARA", "MARC1", "GPAM",
              "FASN", "SCD", "CYP7A1", "FABP1", "AKR1B10")
# Also label top both-signal genes
top_both <- tc[category == "Both COLOC + TWAS"][order(-coloc_pp4 * abs(twas_z))][1:10]$gene
label_genes <- unique(c(priority, top_both))
tc[, show_label := gene %in% label_genes & category != "Neither"]

# Spearman correlation
rho <- cor(tc$coloc_pp4, abs(tc$twas_z), method = "spearman", use = "complete.obs")

p_a <- ggplot(tc, aes(x = abs(twas_z), y = coloc_pp4)) +
  geom_point(data = tc[category == "Neither"],
             color = masld_colors$ns, size = 0.3, alpha = 0.2) +
  geom_point(data = tc[category == "TWAS only"],
             color = masld_colors$twas, size = 0.8, alpha = 0.7) +
  geom_point(data = tc[category == "COLOC only"],
             color = masld_colors$gwas, size = 0.8, alpha = 0.7) +
  geom_point(data = tc[category == "Both COLOC + TWAS"],
             color = masld_colors$up, size = 1.2, alpha = 0.9) +
  geom_hline(yintercept = 0.9, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_vline(xintercept = 4, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_text_repel(
    data = tc[show_label == TRUE],
    aes(label = gene),
    size = 1.8, fontface = "italic",
    max.overlaps = 20, segment.size = 0.15, seed = 42
  ) +
  annotate("text", x = max(abs(tc$twas_z), na.rm = TRUE) * 0.7,
           y = 0.95, size = 2, hjust = 0, color = "grey30",
           label = sprintf("rho = %.3f\nBoth: %d\nCOLOC only: %d\nTWAS only: %d",
                           rho,
                           sum(tc$category == "Both COLOC + TWAS"),
                           sum(tc$category == "COLOC only"),
                           sum(tc$category == "TWAS only"))) +
  scale_x_continuous(limits = c(0, NA)) +
  labs(x = "|Bulk TWAS z-score| (S-PrediXcan)",
       y = "ABF COLOC PP.H4 (best across 24 GWAS)",
       title = "Bulk TWAS vs COLOC concordance") +
  theme_masld(base_size = 7)

# ===========================================================================
# Panel B: INTACT Bayesian filtering effect
# ===========================================================================
cat("Panel B: INTACT filtering...\n")

# Compare raw |TWAS z| rank vs INTACT score rank
# Shows how INTACT re-ranks genes by incorporating COLOC
tc_intact <- tc[!is.na(multi_intact_score)]
tc_intact[, twas_rank := frank(-abs(twas_z))]
tc_intact[, intact_rank := frank(-multi_intact_score)]
tc_intact[, rank_change := twas_rank - intact_rank]  # positive = promoted by INTACT

# Classify: promoted (COLOC confirms TWAS), demoted (TWAS without COLOC = LD artifact)
tc_intact[, intact_effect := fcase(
  multi_intact_score > 0.5 & abs(twas_z) > 2, "Confirmed (TWAS + COLOC)",
  abs(twas_z) > 4 & multi_intact_score < 0.1, "Filtered (TWAS artifact)",
  multi_intact_score > 0.5 & abs(twas_z) < 2, "COLOC-driven",
  default = "Unchanged"
)]

cat("  INTACT effects:\n")
print(tc_intact[, .N, by = intact_effect][order(-N)])

# Labels for promoted/filtered genes
tc_intact[, show_label_b := gene %in% priority |
            (intact_effect == "Filtered (TWAS artifact)" & twas_rank <= 30) |
            (intact_effect == "Confirmed (TWAS + COLOC)" & intact_rank <= 20)]

p_b <- ggplot(tc_intact, aes(x = twas_rank, y = intact_rank)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.3) +
  geom_point(data = tc_intact[intact_effect == "Unchanged"],
             color = masld_colors$ns, size = 0.3, alpha = 0.2) +
  geom_point(data = tc_intact[intact_effect == "Confirmed (TWAS + COLOC)"],
             color = masld_colors$conserved, size = 1, alpha = 0.8) +
  geom_point(data = tc_intact[intact_effect == "Filtered (TWAS artifact)"],
             color = masld_colors$discordant, size = 1, alpha = 0.8) +
  geom_point(data = tc_intact[intact_effect == "COLOC-driven"],
             color = masld_colors$gwas, size = 0.8, alpha = 0.7) +
  geom_text_repel(
    data = tc_intact[show_label_b == TRUE],
    aes(label = gene),
    size = 1.6, fontface = "italic",
    max.overlaps = 20, segment.size = 0.15, seed = 42
  ) +
  annotate("text", x = max(tc_intact$twas_rank) * 0.6,
           y = max(tc_intact$intact_rank) * 0.05,
           size = 2, hjust = 0, color = "grey30",
           label = sprintf("Confirmed: %d\nFiltered: %d\nCOLOC-driven: %d",
                           sum(tc_intact$intact_effect == "Confirmed (TWAS + COLOC)"),
                           sum(tc_intact$intact_effect == "Filtered (TWAS artifact)"),
                           sum(tc_intact$intact_effect == "COLOC-driven"))) +
  labs(x = "scTWAS rank (by |z-score|)",
       y = "Multi-INTACT rank (bulk + scTWAS × COLOC)",
       title = "Multi-INTACT re-ranking: COLOC confirmation vs LD filtering") +
  theme_masld(base_size = 7)

# ===========================================================================
# Panel C: Proportion × COLOC — Hepatocytes
# ===========================================================================
cat("Panel C: Hepatocyte proportion correlation...\n")

hep_prop <- prop[cell_type == "Hepatocytes" & !is.na(coloc_pp4)]

if (nrow(hep_prop) > 0) {
  rho_hep <- cor(hep_prop$coloc_pp4, abs(hep_prop$rho), method = "spearman",
                 use = "complete.obs")

  p_c <- ggplot(hep_prop, aes(x = coloc_pp4, y = rho)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.3) +
    geom_point(aes(color = fdr < 0.05), size = 0.6, alpha = 0.5) +
    geom_smooth(method = "lm", color = masld_colors$deg, linewidth = 0.5,
                se = TRUE, alpha = 0.2) +
    scale_color_manual(values = c("TRUE" = masld_colors$up, "FALSE" = masld_colors$ns),
                       labels = c("TRUE" = "FDR < 0.05", "FALSE" = "NS"), name = NULL) +
    annotate("text", x = 0.8, y = max(hep_prop$rho, na.rm = TRUE) * 0.9,
             label = sprintf("n = %d genes\nrho = %.3f", nrow(hep_prop), rho_hep),
             size = 2, hjust = 0.5, color = "grey30") +
    labs(x = "COLOC PP.H4",
         y = "Spearman rho\n(expression ~ hepatocyte proportion)",
         title = "Hepatocyte proportion correlation") +
    theme_masld(base_size = 7) +
    theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
          legend.background = element_blank(), legend.key.size = unit(0.2, "cm"))
} else {
  p_c <- ggplot() + theme_void() + labs(title = "No hepatocyte data")
}

# ===========================================================================
# Panel D: Proportion × COLOC — Macrophages
# ===========================================================================
cat("Panel D: Macrophage proportion correlation...\n")

mac_prop <- prop[cell_type == "Macrophages" & !is.na(coloc_pp4)]

if (nrow(mac_prop) > 0) {
  rho_mac <- cor(mac_prop$coloc_pp4, abs(mac_prop$rho), method = "spearman",
                 use = "complete.obs")

  p_d <- ggplot(mac_prop, aes(x = coloc_pp4, y = rho)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.3) +
    geom_point(aes(color = fdr < 0.05), size = 0.6, alpha = 0.5) +
    geom_smooth(method = "lm", color = masld_colors$up, linewidth = 0.5,
                se = TRUE, alpha = 0.2) +
    scale_color_manual(values = c("TRUE" = masld_colors$up, "FALSE" = masld_colors$ns),
                       labels = c("TRUE" = "FDR < 0.05", "FALSE" = "NS"), name = NULL) +
    annotate("text", x = 0.8, y = max(mac_prop$rho, na.rm = TRUE) * 0.9,
             label = sprintf("n = %d genes\nrho = %.3f", nrow(mac_prop), rho_mac),
             size = 2, hjust = 0.5, color = "grey30") +
    labs(x = "COLOC PP.H4",
         y = "Spearman rho\n(expression ~ macrophage proportion)",
         title = "Macrophage proportion correlation") +
    theme_masld(base_size = 7) +
    theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
          legend.background = element_blank(), legend.key.size = unit(0.2, "cm"))
} else {
  p_d <- ggplot() + theme_void() + labs(title = "No macrophage data")
}

# ===========================================================================
# Combine and save
# ===========================================================================
cat("Assembling...\n")

p_combined <- (p_a | p_b) / (p_c | p_d) +
  plot_annotation(tag_levels = list(c("a", "b", "c", "d"))) &
  theme(plot.tag = element_text(size = 8, face = "bold"))

ggsave(file.path(outdir, "figS_twas_coloc_intact.pdf"),
       p_combined, width = 170, height = 140, units = "mm", device = cairo_pdf)
cat("Saved: figS_twas_coloc_intact.pdf\n")
cat("Done.\n")
