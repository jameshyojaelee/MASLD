#!/usr/bin/env Rscript
# fig4_gwas_volcano.R — GWAS-annotated DEG volcano plot
#
# Panel for Fig 4 (Causal Genetic Architecture):
# Dream DEGs as grey background, COLOC genes highlighted in color,
# shape = high-PIP fine-mapped variant nearby, size = evidence sources
#
# Output: figures/main/fig4_causal_architecture/panel_gwas_volcano.pdf
#
# SLURM: cpu partition, 4 CPUs, 16GB, 48h

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
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

# Dream DEGs (ENSEMBL IDs)
dream <- fread(file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/dream_results.csv"))

# COLOC gene-level
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

# High-PIP variants mapped to genes
hip <- fread(file.path(BASE, "GWAS/finemapping/results/high_pip_variants.csv"))

# Multi-evidence atlas for evidence count
atlas_file <- file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
atlas <- fread(atlas_file, select = c("human_symbol", "n_coloc_sources"))

# INTACT scores (if available)
intact_file <- file.path(BASE, "RNA-seq/results/gwas_rna_integration/intact_scores.csv")
intact <- if (file.exists(intact_file)) fread(intact_file) else data.table(gene = character())

# ===========================================================================
# 2. Build plot data
# ===========================================================================
cat("Building plot data...\n")

# Start with dream results
plot_dt <- dream[, .(
  gene_symbol,
  logFC,
  neglog10p = pmin(-log10(padj + 1e-300), 50)  # cap at 50
)]

# Add COLOC PP.H4
plot_dt <- merge(plot_dt,
  coloc[, .(gene_symbol = gene, coloc_pp4 = coloc_best_pp4,
            coloc_n_gwas = coloc_n_gwas_h4_05)],
  by = "gene_symbol", all.x = TRUE)

# Add evidence count
plot_dt <- merge(plot_dt,
  atlas[, .(gene_symbol = human_symbol, n_sources = n_coloc_sources)],
  by = "gene_symbol", all.x = TRUE)

# Add INTACT
if (nrow(intact) > 0) {
  plot_dt <- merge(plot_dt,
    intact[, .(gene_symbol = gene, multi_intact_score)],
    by = "gene_symbol", all.x = TRUE)
}

# Classify genes
plot_dt[, category := fcase(
  coloc_pp4 > 0.9 & logFC > 0, "COLOC_up",
  coloc_pp4 > 0.9 & logFC < 0, "COLOC_down",
  default = "Background"
)]

cat("  Categories:\n")
print(plot_dt[, .N, by = category][order(-N)])

# Clamp extreme values
plot_dt[logFC > 5, logFC := 5]
plot_dt[logFC < -5, logFC := -5]

# ===========================================================================
# 3. Label selection
# ===========================================================================
# Priority genes for labeling
priority_genes <- c(
  # Known MASLD causal
  "PNPLA3", "TM6SF2", "HSD17B13", "MBOAT7", "GCKR", "MARC1",
  # Drug targets
  "THRB", "NR1H4", "PPARA", "PPARG",
  # Top COLOC + DEG
  "CELSR2", "CDK18", "AKR1B10", "LPL", "TREM2",
  "FASN", "SCD", "CYP7A1", "FABP1"
)

# Top COLOC genes not in priority list
top_coloc <- plot_dt[category %in% c("COLOC_up", "COLOC_down")][
  order(-coloc_pp4)][1:30]$gene_symbol

label_genes <- unique(c(priority_genes, top_coloc))
plot_dt[, show_label := gene_symbol %in% label_genes & category != "Background"]

# Limit to 25 labels
if (sum(plot_dt$show_label) > 25) {
  labeled <- plot_dt[show_label == TRUE][order(-coloc_pp4)]
  plot_dt[, show_label := FALSE]
  plot_dt[gene_symbol %in% labeled$gene_symbol[1:25], show_label := TRUE]
}

# ===========================================================================
# 4. Plot
# ===========================================================================
cat("Generating volcano plot...\n")

# Order so COLOC genes are drawn on top
plot_dt[, plot_order := fifelse(category == "Background", 0, 1)]
setorder(plot_dt, plot_order)

# Size mapping: cap at 5 sources
plot_dt[, size_val := pmin(n_sources, 5)]
plot_dt[is.na(size_val), size_val := 0]

# Use mapped color aesthetic for proper legend
plot_dt[, category_label := fcase(
  category == "COLOC_up", "COLOC + Upregulated",
  category == "COLOC_down", "COLOC + Downregulated",
  default = "Not colocalized"
)]
plot_dt[, category_label := factor(category_label,
  levels = c("COLOC + Upregulated", "COLOC + Downregulated", "Not colocalized"))]

p <- ggplot(plot_dt, aes(x = logFC, y = neglog10p)) +
  # Background genes (not in legend — draw first)
  geom_point(data = plot_dt[category == "Background"],
             color = "#E0E0E0", size = 0.3, alpha = 0.3, show.legend = FALSE) +
  # COLOC genes with mapped color + legend
  geom_point(data = plot_dt[category != "Background"],
             aes(color = category_label, size = size_val), alpha = 0.85) +
  scale_color_manual(
    values = c("COLOC + Upregulated" = "#D32F2F",
               "COLOC + Downregulated" = "#1565C0",
               "Not colocalized" = "#E0E0E0"),
    name = "COLOC status",
    drop = TRUE
  ) +
  # Significance lines
  geom_hline(yintercept = -log10(0.1), linetype = "dashed",
             color = "grey50", linewidth = 0.3) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed",
             color = "grey50", linewidth = 0.3) +
  # Labels
  geom_text_repel(
    data = plot_dt[show_label == TRUE],
    aes(label = gene_symbol),
    color = "black",
    size = 2, fontface = "italic",
    max.overlaps = 30, segment.size = 0.2,
    min.segment.length = 0, box.padding = 0.3,
    seed = 42
  ) +
  scale_size_continuous(
    range = c(1, 3.5),
    breaks = c(1, 3, 5),
    name = "COLOC\nsources"
  ) +
  labs(
    x = expression("Integrated log"[2]*"FC (Disease vs Control)"),
    y = expression("-log"[10]*"(adjusted "*italic(P)*"-value)")
  ) +
  theme_masld(base_size = 7) +
  theme(legend.position = c(0.02, 0.98),
        legend.justification = c(0, 1),
        legend.background = element_rect(fill = "white", color = NA, linewidth = 0))

# Add annotation: COLOC gene counts
n_coloc_up <- sum(plot_dt$category == "COLOC_up")
n_coloc_dn <- sum(plot_dt$category == "COLOC_down")
n_coloc_sug <- sum(plot_dt$category == "COLOC_suggestive")
p <- p + annotate("text", x = max(plot_dt$logFC) * 0.7, y = max(plot_dt$neglog10p) * 0.95,
                  label = sprintf("COLOC PP.H4 > 0.9: %d\n(up: %d, down: %d)",
                                  n_coloc_up + n_coloc_dn, n_coloc_up, n_coloc_dn),
                  size = 2, hjust = 0, vjust = 1, color = "grey30")

# ===========================================================================
# 5. Save
# ===========================================================================
outfile <- file.path(outdir, "panel_gwas_volcano.pdf")
ggsave(outfile, p, width = 120, height = 85, units = "mm", device = cairo_pdf)
cat("Saved:", outfile, "\n")

cat("Done.\n")
