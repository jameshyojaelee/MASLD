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

# Ancestry is derived from the GWAS registry via gwas_ancestry() (load_figure_data.R,
# single source of truth) — NOT a hardcoded legacy list. The former inline
# ancestry_for_gwas() enumerated only the 23-GWAS legacy portfolio (14 EUR + 3 BBJ +
# 3 Pan-UKBB AFR + 3 Pan-UKBB CSA) and returned default=NA for everything else, which
# silently DROPPED every MVP stratum and had no AMR bin. The registry helper covers the
# full 50-GWAS MVP-expanded portfolio across all 5 ancestries (EUR/AFR/AMR/EAS/SAS).
sc <- fread(file.path(BASE,
  "GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv"))
sc[, ancestry := as.character(gwas_ancestry(gwas_name))]
sc[, pp4_best := pmax(PP.H4.susie, PP.H4.abf, na.rm = TRUE)]
sc[is.infinite(pp4_best), pp4_best := NA_real_]

per_gene_anc <- sc[!is.na(pp4_best),
                   .(max_pp4 = max(pp4_best, na.rm = TRUE)),
                   by = .(gene, ancestry)]

thresholds      <- c(0.5, 0.8, 0.9)
ancestry_levels <- GWAS_ANCESTRY_LEVELS   # EUR, AFR, AMR, EAS, SAS (canonical order)
# n GWAS per ancestry, computed from the portfolio actually present in the COLOC file
# (registry-driven; auto-tracks the 50-GWAS MVP-expanded portfolio incl. AMR, so the
# denominator labels can never drift from the data).
ancestry_n_gwas <- {
  m <- unique(sc[, .(gwas_name, ancestry)])[, .N, by = ancestry]
  setNames(m$N, m$ancestry)[ancestry_levels]
}

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

# Canonical 5-ancestry palette from load_figure_data.R (EUR/AFR/AMR/EAS/SAS; incl. AMR),
# shared with the other Fig 2 ancestry panels for consistency.
ancestry_colors_3b <- ANCESTRY_COLORS
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

out_pdf <- file.path(PANEL_DIR, "ancestry_coloc_counts.pdf")
out_csv <- file.path(FIG3_DIR, "ancestry_coloc_counts.csv")
save_fig(p3b, out_pdf, width = fig_half_width, height = 3.0)
fwrite(bar_dt[, .(ancestry, n_gwas, threshold, n_genes)], out_csv)
cat("[fig3b] Wrote:", out_pdf, "\n")
cat("[fig3b] Wrote:", out_csv, "\n")
print(bar_dt[, .(ancestry, threshold, n_genes)])
