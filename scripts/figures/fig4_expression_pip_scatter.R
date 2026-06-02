#!/usr/bin/env Rscript
# fig4_expression_pip_scatter.R — Expression × COLOC/INTACT Quadrant Scatter
#
# Centerpiece figure proving evidence source orthogonality:
# x = dream logFC, y = COLOC PP.H4 (version A) or INTACT score (version B)
# Marginal density plots showing near-zero correlation (MI annotation)
# Quadrant counts: genetic+expression convergent vs single-source
#
# Output: figures/main/fig4_causal_architecture/panel_expression_pip_scatter.pdf
#
# SLURM: cpu partition, 4 CPUs, 16GB, 48h

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
dir.create(file.path(FIG3_DIR, "panels"), showWarnings = FALSE, recursive = TRUE)
outdir <- file.path(FIG3_DIR, "panels")

# ===========================================================================
# 1. Load data
# ===========================================================================
cat("Loading data...\n")

dream <- fread(file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/dream_results.csv"))

coloc <- fread(file.path(BASE,
  "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv"))
coloc <- coloc[gene != "" & !is.na(gene)]

# Map ENSEMBL → symbol
ensembl_map <- coloc[gene != "" & ensembl != "", .(ensembl, symbol = gene)]
ensembl_map[, ensembl_base := sub("\\.\\d+$", "", ensembl)]
dream[, ensembl_base := sub("\\.\\d+$", "", gene)]
dream <- merge(dream, ensembl_map[, .(ensembl_base, symbol)],
               by = "ensembl_base", all.x = TRUE)
dream[!is.na(symbol), gene_symbol := symbol]
dream[is.na(symbol), gene_symbol := gene]

# INTACT scores
intact_file <- file.path(BASE, "RNA-seq/results/gwas_rna_integration/intact_scores.csv")
intact <- if (file.exists(intact_file)) fread(intact_file) else data.table(gene = character())

# sc-eQTL cell-type assignments
sceqtl_file <- file.path(BASE,
  "RNA-seq/results/causal_inference/sceqtl/sceqtl_multi_celltype_summary.csv")
sceqtl <- if (file.exists(sceqtl_file)) fread(sceqtl_file) else data.table()

# ===========================================================================
# 2. Build plot data
# ===========================================================================
cat("Building scatter data...\n")

# Merge dream + COLOC
scatter <- merge(
  dream[, .(gene_symbol, logFC, padj, t)],
  coloc[, .(gene_symbol = gene, coloc_pp4 = coloc_best_pp4)],
  by = "gene_symbol"
)
cat("  Genes with both expression + COLOC:", nrow(scatter), "\n")

# Add INTACT
if (nrow(intact) > 0) {
  scatter <- merge(scatter,
    intact[, .(gene_symbol = gene, intact_bulk = multi_intact_score,
               intact_ct = multi_intact_score)],
    by = "gene_symbol", all.x = TRUE)
}

# Add cell-type assignment
if (nrow(sceqtl) > 0) {
  scatter <- merge(scatter,
    sceqtl[, .(gene_symbol = gene, best_cell_type)],
    by = "gene_symbol", all.x = TRUE)
}

# ===========================================================================
# 3. Compute MI (mutual information) — confirms orthogonality
# ===========================================================================
cat("Computing mutual information...\n")

# Discretize into bins for MI calculation
scatter[, logFC_bin := cut(logFC, breaks = 20)]
scatter[, coloc_bin := cut(coloc_pp4, breaks = 20)]

# MI = sum over bins of p(x,y) * log(p(x,y) / (p(x)*p(y)))
joint <- scatter[!is.na(logFC_bin) & !is.na(coloc_bin), .N,
                 by = .(logFC_bin, coloc_bin)]
joint[, p_xy := N / sum(N)]
marginal_x <- scatter[!is.na(logFC_bin), .N, by = logFC_bin]
marginal_x[, p_x := N / sum(N)]
marginal_y <- scatter[!is.na(coloc_bin), .N, by = coloc_bin]
marginal_y[, p_y := N / sum(N)]

joint <- merge(joint, marginal_x[, .(logFC_bin, p_x)], by = "logFC_bin")
joint <- merge(joint, marginal_y[, .(coloc_bin, p_y)], by = "coloc_bin")
joint[, mi_term := p_xy * log2(p_xy / (p_x * p_y))]
mi_value <- sum(joint$mi_term, na.rm = TRUE)

spearman_rho <- cor(scatter$logFC, scatter$coloc_pp4, method = "spearman",
                    use = "complete.obs")

cat(sprintf("  MI = %.4f bits, Spearman rho = %.4f\n", mi_value, spearman_rho))

# ===========================================================================
# 4. Define quadrants and labels
# ===========================================================================
# Thresholds
deg_thresh <- 0.5  # |logFC| > 0.5
coloc_thresh <- 0.9  # PP.H4 > 0.9 (stringent)

scatter[, quadrant := fcase(
  coloc_pp4 > coloc_thresh & logFC > deg_thresh, "Genetic + Up",
  coloc_pp4 > coloc_thresh & logFC < -deg_thresh, "Genetic + Down",
  coloc_pp4 > coloc_thresh & abs(logFC) <= deg_thresh, "Genetic only",
  coloc_pp4 <= coloc_thresh & padj < 0.1, "Expression only",
  default = "Neither"
)]

cat("Quadrant counts:\n")
print(scatter[, .N, by = quadrant][order(-N)])

# Label selection — keep sparse to avoid clutter at high PP.H4
# Only label genes that are in the convergent quadrants (Genetic + Up/Down)
priority <- c("PNPLA3", "TM6SF2", "HSD17B13", "MBOAT7", "GCKR",
              "THRB", "NR1H4", "MARC1")

# Top convergent: high PP.H4 AND strong |logFC| — only label the best
top_convergent <- scatter[coloc_pp4 > 0.9 & abs(logFC) > 0.5][
  order(-coloc_pp4 * abs(logFC))][1:8]$gene_symbol

label_genes <- unique(c(priority, top_convergent))
# Only label if gene is above the COLOC threshold or is a priority gene
scatter[, show_label := gene_symbol %in% label_genes &
          (coloc_pp4 > 0.5 | gene_symbol %in% priority)]

# ===========================================================================
# 5. Version A: Expression × COLOC PP.H4
# ===========================================================================
cat("Generating Version A (COLOC PP.H4)...\n")

# Color by quadrant — maximally distinct colors per category
scatter[, point_color := quadrant]
color_scale <- scale_color_manual(
  values = c("Genetic + Up" = "#D32F2F",       # strong red
             "Genetic + Down" = "#1565C0",      # deep blue
             "Genetic only" = "#FF9800",         # orange
             "Expression only" = "#4CAF50",      # green
             "Neither" = "#E0E0E0"),             # light grey
  name = "Category")

# Main scatter
p_main <- ggplot(scatter, aes(x = logFC, y = coloc_pp4)) +
  geom_point(aes(color = point_color),
             size = 0.5, alpha = 0.5) +
  color_scale +
  # Threshold lines
  geom_hline(yintercept = coloc_thresh, linetype = "dashed",
             color = "grey50", linewidth = 0.3) +
  geom_vline(xintercept = c(-deg_thresh, deg_thresh), linetype = "dashed",
             color = "grey50", linewidth = 0.3) +
  # Labels — sparse, well-spaced
  geom_text_repel(
    data = scatter[show_label == TRUE],
    aes(label = gene_symbol),
    size = 2, fontface = "italic",
    max.overlaps = 15, segment.size = 0.2,
    min.segment.length = 0.1, box.padding = 0.5,
    point.padding = 0.3, force = 3, force_pull = 0.5,
    seed = 42
  ) +
  labs(
    x = expression("Integrated log"[2]*"FC"),
    y = "COLOC PP.H4"
  ) +
  theme_masld(base_size = 7) +
  theme(legend.position = "right",
        legend.key.size = unit(0.25, "cm"))

# Marginal density: top (logFC distribution by COLOC status)
p_top <- ggplot(scatter, aes(x = logFC, fill = coloc_pp4 > coloc_thresh)) +
  geom_density(alpha = 0.6, linewidth = 0.3) +
  scale_fill_manual(values = c("TRUE" = masld_colors$up, "FALSE" = masld_colors$ns),
                    labels = c("TRUE" = paste0("PP.H4>", coloc_thresh), "FALSE" = "Other"),
                    name = NULL) +
  theme_masld(base_size = 7) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.line.x = element_blank(),
        legend.position = c(0.02, 0.9), legend.justification = c(0, 1),
        legend.key.size = unit(0.2, "cm"), legend.text = element_text(size = 5),
        legend.background = element_blank(),
        plot.margin = margin(2, 2, 0, 2))

# Marginal right: proportion of DEGs at each PP.H4 bin (not raw counts)
# Raw histogram fails because DEG counts are dwarfed by non-DEGs at low PP.H4
scatter[, is_deg := padj < 0.1 & abs(logFC) > 0.5]
scatter[, pp4_bin := cut(coloc_pp4, breaks = seq(0, 1, 0.05), include.lowest = TRUE)]
pp4_prop <- scatter[!is.na(pp4_bin), .(
  prop_deg = mean(is_deg, na.rm = TRUE),
  n_total = .N
), by = pp4_bin]
pp4_prop[, pp4_mid := seq(0.025, 0.975, 0.05)[as.integer(pp4_bin)]]

p_right <- ggplot(pp4_prop, aes(x = pp4_mid, y = prop_deg)) +
  geom_col(fill = "#4CAF50", width = 0.045, alpha = 0.8) +
  geom_hline(yintercept = mean(scatter$is_deg, na.rm = TRUE),
             linetype = "dashed", color = "grey50", linewidth = 0.3) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.05))) +
  coord_flip() +
  labs(y = "% DEGs") +
  theme_masld(base_size = 7) +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        axis.title.x = element_text(size = 5.5),
        plot.margin = margin(2, 2, 2, 0))

# Assemble with patchwork — wider marginals (2:5 ratio instead of 1:4)
p_a <- p_top + plot_spacer() + p_main + p_right +
  plot_layout(ncol = 2, nrow = 2,
              widths = c(5, 2), heights = c(2, 5))

# ===========================================================================
# 6. Left panel: SuSiE-COLOC PP.H4 × logFC (fine-mapping-based COLOC)
#    Right panel: INTACT score × logFC (Bayesian-filtered TWAS × COLOC)
# ===========================================================================
# Load SuSiE-COLOC for left panel (different from ABF COLOC used above)
atlas_file <- file.path(BASE, "RNA-seq/results/gwas_rna_integration/gwas_rna_scrna_integrated.csv")
if (file.exists(atlas_file)) {
  atlas_cols <- fread(atlas_file, select = c("gene", "coloc_susie_best_pp4"))
  scatter <- merge(scatter, atlas_cols[, .(gene_symbol = gene, susie_pp4 = coloc_susie_best_pp4)],
                   by = "gene_symbol", all.x = TRUE)
}

has_susie <- "susie_pp4" %in% names(scatter) && sum(!is.na(scatter$susie_pp4)) > 100
has_intact <- "intact_ct" %in% names(scatter) && sum(!is.na(scatter$intact_ct)) > 50

if (has_susie && has_intact) {
  cat("Generating combined: SuSiE-COLOC (left) + INTACT (right)...\n")

  # Left: SuSiE-COLOC PP.H4 × logFC
  scatter_susie <- scatter[!is.na(susie_pp4)]
  scatter_susie[, susie_cat := fcase(
    susie_pp4 > 0.9 & logFC > deg_thresh, "Genetic + Up",
    susie_pp4 > 0.9 & logFC < -deg_thresh, "Genetic + Down",
    susie_pp4 > 0.9, "Genetic only",
    padj < 0.1, "Expression only",
    default = "Neither"
  )]

  p_susie <- ggplot(scatter_susie, aes(x = logFC, y = susie_pp4, color = susie_cat)) +
    geom_point(size = 0.5, alpha = 0.5) +
    scale_color_manual(
      values = c("Genetic + Up" = "#D32F2F", "Genetic + Down" = "#1565C0",
                 "Genetic only" = "#FF9800", "Expression only" = "#4CAF50",
                 "Neither" = "#E0E0E0"),
      name = "Category") +
    geom_hline(yintercept = 0.9, linetype = "dashed", color = "grey50", linewidth = 0.3) +
    geom_vline(xintercept = c(-deg_thresh, deg_thresh), linetype = "dashed",
               color = "grey50", linewidth = 0.3) +
    geom_text_repel(
      data = scatter_susie[gene_symbol %in% label_genes & susie_pp4 > 0.5],
      aes(label = gene_symbol), color = "black",
      size = 2, fontface = "italic",
      max.overlaps = 12, segment.size = 0.2,
      box.padding = 0.5, force = 3, seed = 42
    ) +
    labs(x = expression("Integrated log"[2]*"FC"),
         y = "SuSiE-COLOC PP.H4 (fine-mapping)") +
    theme_masld(base_size = 7) +
    theme(legend.position = "bottom", legend.key.size = unit(0.25, "cm"),
          legend.text = element_text(size = 5))

  # Right: INTACT score × logFC
  scatter_intact <- scatter[!is.na(intact_ct)]
  scatter_intact[, intact_cat := fcase(
    intact_ct > 0.5 & logFC > deg_thresh, "Confirmed + Up",
    intact_ct > 0.5 & logFC < -deg_thresh, "Confirmed + Down",
    intact_ct > 0.5, "Confirmed only",
    padj < 0.1, "Expression only",
    default = "Neither"
  )]

  p_intact <- ggplot(scatter_intact, aes(x = logFC, y = intact_ct, color = intact_cat)) +
    geom_point(size = 0.5, alpha = 0.5) +
    scale_color_manual(
      values = c("Confirmed + Up" = "#D32F2F", "Confirmed + Down" = "#1565C0",
                 "Confirmed only" = "#FF9800", "Expression only" = "#4CAF50",
                 "Neither" = "#E0E0E0"),
      name = "Category") +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50", linewidth = 0.3) +
    geom_vline(xintercept = c(-deg_thresh, deg_thresh), linetype = "dashed",
               color = "grey50", linewidth = 0.3) +
    geom_text_repel(
      data = scatter_intact[gene_symbol %in% label_genes & intact_ct > 0.3],
      aes(label = gene_symbol), color = "black",
      size = 2, fontface = "italic",
      max.overlaps = 12, segment.size = 0.2,
      box.padding = 0.5, force = 3, seed = 42
    ) +
    labs(x = expression("Integrated log"[2]*"FC"),
         y = "Multi-INTACT score (bulk TWAS + scTWAS × COLOC)") +
    theme_masld(base_size = 7) +
    theme(legend.position = "bottom", legend.key.size = unit(0.25, "cm"),
          legend.text = element_text(size = 5))

  # Combined side-by-side
  p_combined <- p_susie + p_intact +
    plot_layout(ncol = 2, guides = "collect") &
    theme(legend.position = "bottom")

  ggsave(file.path(outdir, "panel_expression_pip_scatter.pdf"),
         p_combined, width = 180, height = 100, units = "mm", device = cairo_pdf)
  cat("  Saved SuSiE-COLOC + INTACT scatter\n")
} else {
  cat("  Insufficient data for combined scatter\n")
}

# Save Version A standalone
ggsave(file.path(outdir, "panel_expression_coloc_scatter.pdf"),
       p_a, width = 120, height = 110, units = "mm", device = cairo_pdf)

cat("Done.\n")
