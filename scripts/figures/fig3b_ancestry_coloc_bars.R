#!/usr/bin/env Rscript
# Standalone renderer for Fig 3 Panel B (per-ancestry COLOC eGene counts).
# Mirrors the Panel 3b block in fig3_regulatory_architecture_v2.R; isolated
# so the panel can be regenerated quickly without rebuilding the full figure.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

PANEL_DIR <- file.path(FIG3_DIR, "panels")
dir.create(PANEL_DIR, showWarnings = FALSE, recursive = TRUE)

EUR_17 <- c(
  "2019_31311600_NAFLD_EUR", "2020_32298765_NAFLD_EUR",
  "2021_34128465_PDFF_EUR",  "2021_34841290_NAFLD_EUR",
  "2021_34957434_PDFF_EUR",  "2022_36402844_PDFF_EUR",
  "2023_36280732_NAFLD_deCode_EUR", "2023_36280732_NAFLD_Intermountain_EUR",
  "2023_36280732_NAFLD_UKBB_EUR", "FinnGen_HCC", "FinnGen_NAFLD", "FinnGen_NASH",
  "Ghouse_Cirrhosis", "Ghouse_HCC", "UKBB_ALT", "UKBB_AST", "UKBB_GGT"
)
BBJ_5 <- c("2020_32514122_Cirrhosis_EAS", "2020_32514122_HCC_EAS",
           "BBJ_ALT", "BBJ_AST", "BBJ_GGT")
AFR_3 <- c("PanUKBB_AFR_ALT", "PanUKBB_AFR_AST", "PanUKBB_AFR_GGT")
SAS_3 <- c("PanUKBB_CSA_ALT", "PanUKBB_CSA_AST", "PanUKBB_CSA_GGT")

ancestry_for_gwas <- function(g) {
  fcase(
    g %in% EUR_17, "EUR",
    g %in% BBJ_5,  "EAS",
    g %in% AFR_3,  "AFR",
    g %in% SAS_3,  "SAS",
    default       = NA_character_
  )
}

sc <- fread(file.path(BASE,
  "GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv"))
sc[, ancestry := ancestry_for_gwas(gwas_name)]
sc <- sc[!is.na(ancestry)]
sc[, pp4_best := pmax(PP.H4.susie, PP.H4.abf, na.rm = TRUE)]
sc[is.infinite(pp4_best), pp4_best := NA_real_]

per_gene_anc <- sc[!is.na(pp4_best),
                   .(max_pp4 = max(pp4_best, na.rm = TRUE)),
                   by = .(gene, ancestry)]

thresholds      <- c(0.5, 0.8, 0.9)
ancestry_levels <- c("EUR", "EAS", "AFR", "SAS")
ancestry_n_gwas <- c(EUR = 17L, EAS = 5L, AFR = 3L, SAS = 3L)

bar_dt <- rbindlist(lapply(thresholds, function(thr) {
  per_gene_anc[, .(n_genes = sum(max_pp4 > thr)),
               by = ancestry][, threshold := thr]
}))
bar_dt[, ancestry := factor(ancestry, levels = ancestry_levels)]
bar_dt[, threshold_label := factor(sprintf("PP4 > %.1f", threshold),
                                    levels = sprintf("PP4 > %.1f", thresholds))]
bar_dt[, n_gwas := ancestry_n_gwas[as.character(ancestry)]]
bar_dt[, ancestry_label := factor(
  sprintf("%s (n=%d GWAS)", ancestry, n_gwas),
  levels = sprintf("%s (n=%d GWAS)", ancestry_levels,
                   ancestry_n_gwas[ancestry_levels]))]

ancestry_colors_3b <- c(
  "EUR" = "#C9265E",
  "EAS" = "#F4511E",
  "AFR" = "#00695C",
  "SAS" = "#7B1FA2"
)
threshold_alphas <- c("PP4 > 0.5" = 0.45,
                       "PP4 > 0.8" = 0.75,
                       "PP4 > 0.9" = 1.00)

p3b <- ggplot(bar_dt, aes(x = ancestry_label, y = n_genes,
                           fill = ancestry, alpha = threshold_label)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75,
           color = "white", linewidth = 0.2) +
  geom_text(aes(label = n_genes),
            position = position_dodge(width = 0.8),
            vjust = -0.3, size = 1.9, color = "gray25") +
  scale_fill_manual(values = ancestry_colors_3b, guide = "none") +
  scale_alpha_manual(values = threshold_alphas, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL,
       y = "Colocalised eGenes (best PP.H4 across ancestry's GWAS)",
       title = "Per-ancestry COLOC support",
       subtitle = "Best SuSiE PP.H4 per gene; ABF fallback when SuSiE did not converge") +
  theme_masld() +
  theme(legend.position = "top",
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 6),
        plot.title = element_text(size = 8, face = "bold"),
        plot.subtitle = element_text(size = 6, color = "gray35"),
        axis.text.x = element_text(size = 6.5))

out_pdf <- file.path(PANEL_DIR, "fig3b.pdf")
out_csv <- file.path(FIG3_DIR, "fig3b_ancestry_coloc_counts.csv")
save_fig(p3b, out_pdf, width = fig_half_width, height = 3.0)
fwrite(bar_dt[, .(ancestry, n_gwas, threshold, n_genes)], out_csv)
cat("[fig3b] Wrote:", out_pdf, "\n")
cat("[fig3b] Wrote:", out_csv, "\n")
print(bar_dt[, .(ancestry, threshold, n_genes)])
