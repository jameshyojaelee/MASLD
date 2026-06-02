#!/usr/bin/env Rscript
# figS_mediation.R -- Supplementary: COLOC-DEG direction QC
#
# 2 panels (original panels a-b were MR-based; dropped 2026-04-22 when MR
# was ditched from the paper):
# (c) Proportion x COLOC: hepatocyte rho vs coloc_pp4
# (d) COLOC PP.H4 for up vs down DEGs
#
# Output: figures/supplementary/figS04_coloc/figS_mediation.pdf

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
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
# ---------------------------------------------------------------------------
# MR mediation results (mr_mediation_results.csv) retired alongside Script
# 203 on 2026-04-22 when MR was ditched from the paper. Panels (a) and (b)
# previously backed by that file have been removed.
# ===========================================================================
prop <- fread(file.path(BASE, "RNA-seq/results/gwas_rna_integration/proportion_coloc_correlation.csv"))
geneset <- fread(file.path(BASE, "RNA-seq/results/gwas_rna_integration/geneset_enrichment_results.csv"))

cat("Proportion correlations:", nrow(prop), "tests\n")

# ===========================================================================
# (c) Proportion x COLOC: hepatocyte correlation
# ===========================================================================
hep_prop <- prop[cell_type == "Hepatocytes"]
if (nrow(hep_prop) > 0) {
  p_c <- ggplot(hep_prop, aes(x = coloc_pp4, y = rho)) +
    geom_point(aes(color = fdr < 0.05), size = 0.8, alpha = 0.5) +
    geom_smooth(method = "lm", color = masld_colors$deg, linewidth = 0.5, se = FALSE) +
    scale_color_manual(values = c("TRUE" = masld_colors$up, "FALSE" = masld_colors$ns),
                       labels = c("TRUE" = "FDR < 0.05", "FALSE" = "ns"), name = NULL) +
    labs(x = "COLOC PP.H4", y = "Spearman rho\n(expression x hepatocyte proportion)",
         title = "Proportion-expression correlation") +
    theme_masld(base_size = 7)
} else {
  p_c <- ggplot() + theme_void() + labs(title = "No hepatocyte data")
}

# ===========================================================================
# (d) COLOC PP.H4 by DEG direction
# ===========================================================================
up_enrich <- geneset[gene_set == "DEGs_Up"]
dn_enrich <- geneset[gene_set == "DEGs_Down"]

dir_dt <- data.table(
  direction = c("Upregulated", "Downregulated"),
  fold_enrichment = c(up_enrich$fold_enrichment, dn_enrich$fold_enrichment),
  p_value = c(up_enrich$wilcox_p_coloc, dn_enrich$wilcox_p_coloc),
  n_genes = c(up_enrich$n_in_set, dn_enrich$n_in_set)
)

p_d <- ggplot(dir_dt, aes(x = direction, y = fold_enrichment, fill = direction)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_text(aes(label = sprintf("p=%.1e\nn=%s", p_value, format(n_genes, big.mark = ","))),
            vjust = -0.3, size = 2) +
  scale_fill_manual(values = c("Upregulated" = masld_colors$up,
                               "Downregulated" = masld_colors$down)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(x = NULL, y = "Fold enrichment\n(COLOC in DEGs vs background)",
       title = "GWAS enrichment by DEG direction") +
  theme_masld(base_size = 7) +
  theme(legend.position = "none")

# ===========================================================================
# Combine and save (single row of the two surviving COLOC panels)
# ===========================================================================
p_combined <- (p_c | p_d) +
  plot_annotation(tag_levels = list(c("c", "d")))

ggsave(file.path(outdir, "figS_mediation.pdf"),
       p_combined, width = 170, height = 70, units = "mm", device = cairo_pdf)
cat("Saved: figS_mediation.pdf\n")
cat("Done.\n")
