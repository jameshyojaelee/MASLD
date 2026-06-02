#!/usr/bin/env Rscript
# figS_gwas_celltype.R
# Figure for Analysis F2 — GWAS variant × cell-type architecture.
#
# Panels:
#   A — Risk gene × cell-type heat (variant density per cell type peaks)
#   B — Risk genes with highest credible-set variant count + max PIP
#   C — Cell-type risk gene load
#   D — Ligand/receptor role among risk genes
#
# Output: figures/supplementary/figS_celltype_biology/figS_F2_gwas_celltype.pdf

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})

source("scripts/figures/publication_theme.R")
source("scripts/figures/load_figure_data.R")

IN <- file.path(BASE, "RNA-seq/results/gwas_celltype")
roles <- fread(file.path(IN, "risk_gene_celltype_roles.csv"))
cov   <- fread(file.path(IN, "risk_gene_celltype_coverage.csv"))
ct_load <- fread(file.path(IN, "celltype_risk_gene_load.csv"))

# Panel A — variant heatmap
GWAS_ATAC <- file.path(BASE, "GWAS/finemapping/results/gwas_atac")
vann <- fread(file.path(GWAS_ATAC, "gwas_atac_variant_annotation.csv"))
vann[, atac_gene := fifelse(linked_gene != "" & !is.na(linked_gene),
                             linked_gene, nearest_gene)]
vann_risk <- vann[atac_gene %in% roles$atac_gene]
heat <- vann_risk[, .N, by = .(atac_gene, cell_type)]
heat[, atac_gene := factor(atac_gene, levels = sort(unique(atac_gene)))]

pA <- ggplot(heat, aes(cell_type, atac_gene, fill = N)) +
  geom_tile(color = "white", linewidth = 0.2) +
  geom_text(aes(label = N), size = 2.1) +
  scale_fill_gradient(low = "grey95", high = "#C0392B", name = "N variants\nin peak") +
  labs(x = "Cell type", y = NULL,
       title = "Credible-set variants per risk gene × cell type") +
  theme_masld() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 6.5),
        axis.text.y = element_text(size = 6.5))

# Panel B — n variants + PIP
topB <- cov[n_variants_in_peak >= 1][order(-max_pip, -n_variants_in_peak)][1:min(15, .N)]
topB[, atac_gene := factor(atac_gene, levels = rev(atac_gene))]
pB <- ggplot(topB, aes(n_variants_in_peak, atac_gene, fill = max_pip)) +
  geom_col(width = 0.7, color = "white") +
  geom_text(aes(label = sprintf("max_pip=%.2f", max_pip)),
            hjust = -0.1, size = 2.2) +
  scale_fill_gradient(low = "#FED976", high = "#C0392B", name = "max PIP") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.3))) +
  labs(x = "N credible-set variants in peak", y = NULL,
       title = "Risk genes ranked by variant density") +
  theme_masld() +
  theme(axis.text.y = element_text(size = 6.5))

# Panel C — cell-type risk gene load
ct_load[, cell_type := factor(cell_type, levels = rev(cell_type))]
pC <- ggplot(ct_load, aes(n_risk_genes, cell_type)) +
  geom_col(fill = "#4472C4", width = 0.65, color = "white") +
  geom_text(aes(label = n_risk_genes), hjust = -0.2, size = 2.4) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Distinct risk genes", y = NULL,
       title = "Cell-type risk gene load") +
  theme_masld() +
  theme(axis.text.y = element_text(size = 6.5))

# Panel D — ligand/receptor roles
role_dt <- data.table(
  role = c("Ligand","Receptor","Neither"),
  N = c(sum(roles$is_ligand, na.rm = TRUE),
        sum(roles$is_receptor, na.rm = TRUE),
        sum(!roles$is_ligand & !roles$is_receptor, na.rm = TRUE))
)
pD <- ggplot(role_dt, aes(role, N, fill = role)) +
  geom_col(width = 0.6, color = "white") +
  geom_text(aes(label = N), vjust = -0.3, size = 2.8) +
  scale_fill_brewer(palette = "Set1", guide = "none") +
  labs(x = NULL, y = "N risk genes",
       title = "Risk gene LIANA role (ligand/receptor/neither)") +
  theme_masld()

fig <- (pA | pB) / (pC | pD) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 8, face = "bold"))

out_path <- file.path(FIGS_CELLTYPE_DIR, "figS_F2_gwas_celltype.pdf")
ggsave(out_path, fig, width = 14, height = 11)
message("Saved: ", out_path)
