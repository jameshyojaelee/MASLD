#!/usr/bin/env Rscript
# 214_spatial_figures.R — Spatial-stratified COLOC figures (Module C)
#
# Panel A: Periportal vs pericentral COLOC enrichment (forest plot with OR + 95% CI)
# Panel B: Spatial SVG × COLOC overlap (UpSet-style bar chart)
#
# Inputs:
#   - RNA-seq/results/stratified_causal/spatial_coloc_enrichment.csv  (from Script 213)
#   - RNA-seq/results/stratified_causal/zone_causal_scores.csv        (from Script 213)
#   - GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv
#   - Analysis/Spatial/results/svg/differential_svgs.csv
#
# Outputs:
#   - figures/supplementary/figS04_coloc/figS_spatial_enrichment.pdf
#   - figures/supplementary/figS04_coloc/figS_spatial_venn.pdf
#   - figures/supplementary/figS04_coloc/figS_spatial_combined.pdf
#
# SLURM: cpu partition, 4 CPUs, 32GB, 48h

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(grid)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

# Source publication theme
source(file.path(BASE, "scripts/figures/publication_theme.R"))

figdir <- file.path(BASE, "figures/supplementary/figS04_coloc")
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

# ===========================================================================
# 1. Load Script 213 outputs
# ===========================================================================
cat("Loading Script 213 outputs...\n")
enrichment <- fread(file.path(BASE, "RNA-seq/results/stratified_causal/spatial_coloc_enrichment.csv"))
zone_scores <- fread(file.path(BASE, "RNA-seq/results/stratified_causal/zone_causal_scores.csv"))

cat("  Enrichment tests:", nrow(enrichment), "\n")
cat("  Zone scores:", nrow(zone_scores), "\n")

# ===========================================================================
# 2. Panel A: Forest plot — zonation COLOC enrichment
# ===========================================================================
cat("\nGenerating Panel A: Forest plot...\n")

# Select the primary PP.H4 > 0.5 tests across all modules
forest_data <- enrichment[grepl("PP4>0.5", test_label) | module == "GWAS_ATAC_zonation"]

# Clean labels for display
forest_data[, label := test_label]
forest_data[, label := gsub("COLOC_PP4>0.5 . ", "", label)]
forest_data[, label := gsub("Hep_ATAC_peaks . ", "", label)]
forest_data[, label := gsub("_", " ", label)]

# Cap extreme OR and CI for plotting
forest_data[, or_plot := pmin(odds_ratio, 50)]
forest_data[, ci_upper_plot := pmin(ci_upper, 100)]
# Handle Inf/0 in CI
forest_data[is.infinite(ci_upper_plot) | is.na(ci_upper_plot), ci_upper_plot := 50]
forest_data[is.na(ci_lower), ci_lower := 0]

# Significance indicator
forest_data[, sig_label := fcase(
  padj < 0.001, "***",
  padj < 0.01,  "**",
  padj < 0.05,  "*",
  default = "ns"
)]

# Clean module labels for faceting
forest_data[, module_label := fcase(
  module == "COLOC_zonation", "COLOC x Zonation",
  module == "SVG_COLOC", "SVG x COLOC",
  module == "GWAS_ATAC_zonation", "GWAS-ATAC x Zonation"
)]

forest_data[, module_label := factor(module_label,
  levels = c("COLOC x Zonation", "SVG x COLOC", "GWAS-ATAC x Zonation"))]

pA <- ggplot(forest_data, aes(x = or_plot, y = reorder(label, or_plot))) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper_plot),
                 height = 0.2, colour = "grey40") +
  geom_point(aes(colour = padj < 0.05), size = 3) +
  geom_text(aes(label = paste0("n=", n_overlap, " ", sig_label)),
            hjust = -0.15, size = 2.8, colour = "grey30") +
  scale_colour_manual(values = c("TRUE" = masld_colors$up, "FALSE" = masld_colors$ns),
                      labels = c("TRUE" = "padj < 0.05", "FALSE" = "ns"),
                      name = "Significance") +
  scale_x_log10() +
  facet_grid(module_label ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Odds Ratio (log scale)", y = NULL,
       title = "Spatial and zonation enrichment of COLOC genes") +
  theme_masld() +
  theme(strip.text.y = element_text(angle = 0, hjust = 0, size = 8),
        legend.position = "bottom")

pdf(file.path(figdir, "figS_spatial_enrichment.pdf"), width = 7, height = 6)
print(pA)
dev.off()
cat("  Saved figS_spatial_enrichment.pdf\n")

# ===========================================================================
# 3. Panel B: SVG × COLOC overlap (UpSet-style bar chart)
# ===========================================================================
cat("\nGenerating Panel B: SVG x COLOC overlap...\n")

# Load raw data for set computation
coloc <- fread(file.path(BASE, "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv"))
coloc <- coloc[gene != "" & !is.na(gene)]
coloc_genes <- coloc[coloc_best_pp4 > 0.5, gene]

diff_svgs <- fread(file.path(BASE, "Analysis/Spatial/results/svg/differential_svgs.csv"))
if (names(diff_svgs)[1] == "V1") setnames(diff_svgs, "V1", "gene")
svg_genes <- diff_svgs[(svg_healthy == TRUE | svg_masld == TRUE), gene]

# DEGs from zonation table
deg_genes <- zone_scores[is_deg == TRUE, human_symbol]

# Compute 3-way set intersections
universe <- unique(c(coloc$gene, diff_svgs$gene, zone_scores$human_symbol))

set_deg <- intersect(deg_genes, universe)
set_coloc <- intersect(coloc_genes, universe)
set_svg <- intersect(svg_genes, universe)

# Build membership for each gene in universe
membership <- data.table(
  gene = universe,
  DEG = universe %in% set_deg,
  COLOC = universe %in% set_coloc,
  SVG = universe %in% set_svg
)

# Create intersection labels
membership[, combo := paste0(
  fifelse(DEG, "DEG", ""),
  fifelse(DEG & (COLOC | SVG), "+", ""),
  fifelse(COLOC, "COLOC", ""),
  fifelse((DEG | COLOC) & SVG, "+", ""),
  fifelse(SVG, "SVG", "")
)]

# Cleaner approach: enumerate all 7 non-empty intersections
membership[, group := fcase(
  DEG & COLOC & SVG,     "DEG + COLOC + SVG",
  DEG & COLOC & !SVG,    "DEG + COLOC",
  DEG & !COLOC & SVG,    "DEG + SVG",
  !DEG & COLOC & SVG,    "COLOC + SVG",
  DEG & !COLOC & !SVG,   "DEG only",
  !DEG & COLOC & !SVG,   "COLOC only",
  !DEG & !COLOC & SVG,   "SVG only",
  default = NA_character_
)]

# Count per group (exclude genes in no set)
group_counts <- membership[!is.na(group), .N, by = group]
setorder(group_counts, -N)

# Colour by number of sets in the intersection
group_counts[, n_sets := nchar(gsub("[^+]", "", group)) + 1L]
group_counts[grepl("only", group), n_sets := 1L]

set_fill <- c("1" = "#BDBDBD", "2" = "#42A5F5", "3" = masld_colors$up)

cat("  Set intersection sizes:\n")
print(group_counts)

pB <- ggplot(group_counts, aes(x = reorder(group, N), y = N, fill = factor(n_sets))) +
  geom_col(width = 0.7) +
  geom_text(aes(label = N), hjust = -0.1, size = 3) +
  coord_flip() +
  scale_fill_manual(values = set_fill, name = "Sets",
                    labels = c("1" = "Single", "2" = "Two-way", "3" = "Three-way")) +
  labs(x = NULL, y = "Number of genes",
       title = "DEG / COLOC / SVG overlap") +
  theme_masld() +
  theme(legend.position = "bottom")

pdf(file.path(figdir, "figS_spatial_venn.pdf"), width = 5.5, height = 5)
print(pB)
dev.off()
cat("  Saved figS_spatial_venn.pdf\n")

# ===========================================================================
# 4. Combined figure
# ===========================================================================
cat("\nGenerating combined figure...\n")

combined <- pA + pB +
  plot_layout(widths = c(3, 2)) +
  plot_annotation(tag_levels = "A",
                  theme = theme(plot.title = element_text(size = 14, face = "bold")))

pdf(file.path(figdir, "figS_spatial_combined.pdf"), width = 12, height = 6)
print(combined)
dev.off()
cat("  Saved figS_spatial_combined.pdf\n")

# ===========================================================================
# 5. Summary
# ===========================================================================
cat("\n=== Script 214 complete ===\n")
cat("Outputs:\n")
cat("  ", file.path(figdir, "figS_spatial_enrichment.pdf"), "\n")
cat("  ", file.path(figdir, "figS_spatial_venn.pdf"), "\n")
cat("  ", file.path(figdir, "figS_spatial_combined.pdf"), "\n")
