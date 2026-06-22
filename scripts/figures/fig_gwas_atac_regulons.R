#!/usr/bin/env Rscript
# fig_gwas_atac_regulons.R
# MASLD fine-mapped variants disrupting TF binding motifs at disease regulons.
#
# Panels:
#   A. Scatter  — all PIP>=0.2 variant-TF pairs: max_pip vs |alleleDiff|,
#                 colored by disease-regulon membership; top hits labeled
#                 (TF @ nearest gene)
#   B. Bar      — # of credible variants (PIP>=0.2) disrupting each
#                 disease-regulon TF, with per-TF variant lists (nearest gene)
#   C. Heatmap  — high-confidence (PIP>=0.8) variants × disrupted TFs;
#                 annotated with the variant's assigned / nearest gene

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(scales)
})

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUT_DIR <- FIGS05_DIR   # epigenomic + spatial supplementary
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
cat("Output directory:", OUT_DIR, "\n")

PANEL_THEME <- theme_masld(base_size = 8) +
  theme(
    plot.title       = element_text(size = 8.5, face = "bold", hjust = 0,
                                    margin = margin(b = 3)),
    plot.subtitle    = element_text(size = 6.8, colour = "#555555",
                                    margin = margin(b = 5)),
    axis.title       = element_text(size = 7.5),
    axis.text        = element_text(size = 7),
    legend.text      = element_text(size = 6.5),
    legend.title     = element_text(size = 7),
    legend.key.size  = unit(0.32, "cm"),
    panel.grid.minor = element_blank()
  )

REG_COL <- "#C2185B"   # disease-regulon TF hits
OTH_COL <- "#BDBDBD"   # non-regulon TFs
HI_COL  <- "#1565C0"   # max PIP highlight

# =============================================================================
# Load + annotate data
# =============================================================================
cat("\n-- Loading motif disruption + variant annotation --\n")

motif <- fread(file.path(BASE,
  "GWAS/finemapping/results/gwas_atac/motif_disruption_scores.csv"))
ms <- unique(motif[, .(SNP_id, tf_name, alleleDiff, max_pip,
                        motif_in_disease_regulon, activity_padj)])

va <- fread(file.path(BASE,
  "GWAS/finemapping/results/gwas_atac/gwas_atac_variant_annotation.csv"))
# one record per variant (prefer rows with nearest_gene annotated)
va[, has_gene := !is.na(nearest_gene) & nearest_gene != ""]
setorder(va, -has_gene, -max_pip)
va_uniq <- va[, .SD[1], by = variant_id][,
  .(variant_id, nearest_gene, linked_gene, coloc_best_pp4, coloc_best_gwas)]

merged <- merge(ms, va_uniq, by.x = "SNP_id", by.y = "variant_id", all.x = TRUE)

# keep credible (GWAS-supported) disruptors
credible <- merged[max_pip >= 0.2]
cat("  PIP >= 0.2 variant-TF pairs:", nrow(credible), "\n")
cat("  Unique variants :", length(unique(credible$SNP_id)), "\n")
cat("  Unique TFs      :", length(unique(credible$tf_name)), "\n")
cat("  Disease-regulon TF hits:",
    length(unique(credible[motif_in_disease_regulon == TRUE]$tf_name)),
    "\n")

# Helper: clean variant locus label "chrN:pos"
short_locus <- function(snp_id) gsub("^(\\d+):(\\d+):.+$",
                                      "chr\\1:\\2", snp_id)

credible[, locus       := short_locus(SNP_id)]
credible[, gene_label  := ifelse(!is.na(nearest_gene) & nearest_gene != "",
                                  nearest_gene, "intergenic")]
credible[, abs_diff    := abs(alleleDiff)]
credible[, regulon_lab := ifelse(motif_in_disease_regulon,
                                  "SCENIC+ MASLD disease-regulon TF",
                                  "Other TF")]
credible[, regulon_lab := factor(regulon_lab,
                                  levels = c("SCENIC+ MASLD disease-regulon TF",
                                             "Other TF"))]

# =============================================================================
# Panel A: Scatter — PIP vs motif disruption magnitude
# =============================================================================
cat("\n-- Panel A: PIP vs disruption scatter --\n")

# label the disease-regulon hits + top-PIP non-regulon hits
credible[, label := ""]
credible[motif_in_disease_regulon == TRUE,
         label := paste0(tf_name, " @ ", gene_label)]
# top 5 non-regulon hits (highest PIP * abs_diff)
credible[, score := max_pip * abs_diff]
top_novel_idx <- credible[motif_in_disease_regulon == FALSE
                          ][order(-score)]$SNP_id[1:6]
credible[motif_in_disease_regulon == FALSE & SNP_id %in% top_novel_idx
         & !duplicated(paste(SNP_id, tf_name)),
         label := paste0(tf_name, " @ ", gene_label)]
# dedupe labels (one label per disease-regulon point)
credible[duplicated(paste(label, regulon_lab)), label := ""]

pA <- ggplot(credible, aes(x = max_pip, y = abs_diff)) +
  # PIP threshold reference
  geom_vline(xintercept = 0.8, linetype = "dashed",
             linewidth = 0.35, colour = "#888888") +
  geom_vline(xintercept = 0.5, linetype = "dotted",
             linewidth = 0.3,  colour = "#AAAAAA") +
  # non-regulon (background)
  geom_point(data = credible[motif_in_disease_regulon == FALSE],
             aes(size = abs_diff),
             colour = OTH_COL, alpha = 0.55, shape = 16) +
  # disease-regulon (foreground)
  geom_point(data = credible[motif_in_disease_regulon == TRUE],
             aes(size = abs_diff),
             colour = REG_COL, alpha = 0.95, shape = 16) +
  geom_text_repel(aes(label = label),
                  colour = "black",
                  size = 2.3, segment.size = 0.25,
                  min.segment.length = 0.1,
                  box.padding = 0.3,
                  max.overlaps = Inf, seed = 42) +
  annotate("text", x = 0.8,  y = max(credible$abs_diff) * 1.05,
           label = "PIP 0.8", hjust = -0.1, size = 2.3,
           colour = "#666666") +
  annotate("text", x = 0.5,  y = max(credible$abs_diff) * 1.05,
           label = "PIP 0.5", hjust = -0.1, size = 2.3,
           colour = "#888888") +
  scale_size_continuous(range = c(1.2, 3.8), guide = "none") +
  scale_colour_manual(values = c(`SCENIC+ MASLD disease-regulon TF` = REG_COL,
                                  `Other TF`                        = "#555555"),
                      guide = "none") +
  scale_x_continuous(limits = c(0.18, 1.05),
                     breaks = c(0.2, 0.5, 0.8, 1.0)) +
  labs(x = "GWAS max PIP", y = "|motif alleleDiff|",
       title = "Credible MASLD variants: predicted TF-motif disruption (motifbreakR)",
       subtitle = sprintf(paste0(
         "%d variant-TF pairs (PIP >= 0.2); magenta = SCENIC+ MASLD ",
         "disease-regulon TFs (n = %d pairs, %d TFs)"),
         nrow(credible),
         sum(credible$motif_in_disease_regulon),
         length(unique(credible[motif_in_disease_regulon==TRUE]$tf_name)))) +
  PANEL_THEME

ggsave(file.path(OUT_DIR, "figS_gwas_atac_regulons_panelA_pip_vs_disruption.pdf"), pA,
       width = 7.0, height = 4.5, useDingbats = FALSE)
cat("  Saved panelA\n")

# =============================================================================
# Panel B: High-PIP variant × TF heatmap (PIP >= 0.8)
# =============================================================================
cat("\n-- Panel B: high-PIP variant × TF heatmap --\n")

hi <- credible[max_pip >= 0.8]
cat("  PIP >= 0.8 pairs (all TFs):", nrow(hi), "\n")

# Collapse: keep only TFs that are either (a) in disease regulon
# or (b) hit by >= 2 high-PIP variants. This drops the long tail of
# TFs disrupted by a single variant and makes the matrix readable.
tf_stats <- hi[, .(n_vars = uniqueN(SNP_id),
                   regulon = any(motif_in_disease_regulon)),
               by = tf_name]
keep_tfs <- tf_stats[n_vars >= 2 | regulon == TRUE]$tf_name
hi_f <- hi[tf_name %in% keep_tfs]
cat("  Collapsed to", length(keep_tfs),
    "TFs (disease-regulon OR hit >=2 variants):", nrow(hi_f), "pairs\n")

# distance to TSS context (for y-axis annotation)
hi_dist <- hi_f[, .(d = max(0, na.rm = TRUE))]   # placeholder; use gwas_atac_variant_annotation
d_stats <- va_uniq[variant_id %in% unique(hi_f$SNP_id)]
# join distance_to_tss from original va
va_src <- fread(file.path(BASE,
  "GWAS/finemapping/results/gwas_atac/gwas_atac_variant_annotation.csv"))
d_per_var <- unique(va_src[variant_id %in% unique(hi_f$SNP_id),
                            .(SNP_id = variant_id,
                              distance_to_tss = distance_to_tss)])
d_per_var <- d_per_var[, .(d_tss = min(distance_to_tss, na.rm = TRUE)),
                       by = SNP_id]
max_d_kb <- round(max(d_per_var$d_tss, na.rm = TRUE) / 1000, 1)
cat("  Max nearest-gene distance across these variants:",
    max_d_kb, "kb\n")

# variant labels, ordered by PIP
hi_f[, var_label := paste0(locus, "  (", gene_label, ")")]
var_order <- hi_f[, .(mp = max(max_pip)), by = var_label][order(-mp)]$var_label
hi_f[, var_label := factor(var_label, levels = rev(var_order))]

# TF order: disease-regulon first, then by # variants hit (desc)
tf_order_hi <- tf_stats[tf_name %in% keep_tfs
                        ][order(-regulon, -n_vars)]$tf_name
hi_f[, tf_name := factor(tf_name, levels = tf_order_hi)]

# Add missing (variant, TF) combos as NA so we see blank tiles where
# no motif disruption was detected — makes the density explicit.
grid_all <- CJ(var_label = levels(hi_f$var_label),
               tf_name   = levels(hi_f$tf_name))
hi_grid <- merge(grid_all,
                 hi_f[, .(var_label, tf_name, alleleDiff,
                          motif_in_disease_regulon)],
                 by = c("var_label","tf_name"), all.x = TRUE)
hi_grid[, var_label := factor(var_label, levels = levels(hi_f$var_label))]
hi_grid[, tf_name   := factor(tf_name,   levels = levels(hi_f$tf_name))]

regulon_tfs <- tf_stats[regulon == TRUE]$tf_name

pC <- ggplot(hi_grid, aes(x = tf_name, y = var_label)) +
  # full grid background
  geom_tile(fill = "#FAFAFA", colour = "white", linewidth = 0.6) +
  # disrupted tiles
  geom_tile(data = hi_grid[!is.na(alleleDiff)],
            aes(fill = alleleDiff),
            colour = "white", linewidth = 0.6) +
  scale_fill_gradient2(
    low = "#1565C0", mid = "#FFFFFF", high = "#C2185B",
    midpoint = 0, limits = c(-2.5, 2.5),
    breaks = c(-2, -1, 0, 1, 2), name = "alleleDiff\n(alt − ref)"
  ) +
  scale_x_discrete(position = "top") +
  labs(
    x = NULL,
    y = sprintf("Variant (nearest gene within %.0f kb)", ceiling(max_d_kb)),
    title = "High-confidence MASLD variants (PIP >= 0.8): predicted motif disruption (motifbreakR)",
    subtitle = sprintf(paste0(
      "%d variants × %d TFs (disease-regulon + TFs hit by >=2 variants); ",
      "bold = MASLD disease-regulon TF"),
      length(unique(hi_f$var_label)),
      length(unique(hi_f$tf_name)))
  ) +
  PANEL_THEME +
  theme(
    axis.text.x = element_text(angle = 55, hjust = 0, vjust = 0,
                                face = ifelse(levels(hi_f$tf_name)
                                              %in% regulon_tfs,
                                              "bold", "plain"),
                                colour = "black",
                                size = 6.8),
    axis.text.y = element_text(family = "mono", size = 6.5),
    panel.grid  = element_blank(),
    legend.position = "right"
  )

ggsave(file.path(OUT_DIR, "figS_gwas_atac_regulons_panelB_high_pip_heatmap.pdf"), pC,
       width = 6.5, height = 4.0, useDingbats = FALSE)
cat("  Saved panelB\n")

# =============================================================================
# Save the underlying summary table (TSV) for the supplementary methods
# =============================================================================
out_tsv <- file.path(OUT_DIR, "figS_gwas_atac_regulons_hits.tsv")
fwrite(credible[order(-max_pip),
                 .(SNP_id, locus, tf_name, nearest_gene, linked_gene,
                   max_pip, alleleDiff, abs_diff, motif_in_disease_regulon,
                   activity_padj, coloc_best_pp4, coloc_best_gwas)],
       out_tsv, sep = "\t")
cat("\nWrote supporting table:", out_tsv, "\n")

# =============================================================================
# Composite layout (A top, B middle, C bottom)
# =============================================================================
cat("\n-- Assembling composite --\n")

composite <- (pA / pC) +
  plot_layout(heights = c(1.2, 1)) +
  plot_annotation(
    title = "MASLD fine-mapped variants disrupt disease-regulon TF motifs",
    subtitle = paste0(
      "PIP-credible variants (>= 0.2) × hepatocyte scATAC peaks × motifbreakR · ",
      "annotated with GWAS fine-mapping max PIP and nearest gene"
    ),
    theme = theme(
      plot.title    = element_text(size = 10, face = "bold"),
      plot.subtitle = element_text(size = 7, colour = "#444444",
                                   margin = margin(b = 6))
    )
  )

ggsave(file.path(OUT_DIR, "figS_gwas_atac_regulons.pdf"), composite,
       width = 10, height = 8, useDingbats = FALSE)

# clean up deprecated panel files from prior runs
for (f in c("figS_gwas_atac_regulons_panelC_high_pip_heatmap.pdf",
            "figS_gwas_atac_regulons_panelB_tf_ranking.pdf")) {
  old <- file.path(OUT_DIR, f)
  if (file.exists(old)) file.remove(old)
}

cat("\nDone. Outputs in:", OUT_DIR, "\n")
cat("  figS_gwas_atac_regulons_panelA_pip_vs_disruption.pdf\n")
cat("  figS_gwas_atac_regulons_panelB_high_pip_heatmap.pdf\n")
cat("  figS_gwas_atac_regulons.pdf                     (composite)\n")
cat("  figS_gwas_atac_regulons_hits.tsv                (data table)\n")
