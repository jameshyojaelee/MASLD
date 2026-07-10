#!/usr/bin/env Rscript
# figS_coloc_method_concordance.R  (2026-06-17)  — Fig 2 supplement (reviewer armor)
# Fine-mapping + colocalization conclusions are robust to method and LD-reference
# choice. Three facets:
#   A. SuSiE-coloc vs coloc.abf PP.H4 (why SuSiE-coloc is primary; methods cross
#      over: ABF nominates more at >0.5, SuSiE-coloc more confident at >0.9).
#   B. SuSiE vs CARMA fine-mapping PIP concordance (dual-method credible sets).
#   C. PP.H4 concordance across LD reference panels (1KG / TOP-LD / UKBB).
# All numbers computed from disk.
#
# Out: figures/main/fig2_genetics/panels/coloc_method_concordance.pdf (+ source CSV)
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
PANEL_DIR <- file.path(FIG3_DIR, "panels")

# ---- A. ABF vs SuSiE-coloc PP.H4 (gene level) ----
gl <- fread(file.path(BASE, "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv"))
gl[, abf := suppressWarnings(as.numeric(coloc_best_pp4))]
gl[, su  := suppressWarnings(as.numeric(coloc_best_susie_pp4))]
ab <- gl[!is.na(abf) & !is.na(su)]
rA <- cor(ab$abf, ab$su)
nA_susie <- sum(gl$su > 0.5, na.rm = TRUE); nA_abf <- sum(gl$abf > 0.5, na.rm = TRUE)
pA <- ggplot(ab, aes(abf, su)) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.3, color = "grey70") +
  geom_hline(yintercept = 0.5, linetype = "22", linewidth = 0.25, color = "grey60") +
  geom_vline(xintercept = 0.5, linetype = "22", linewidth = 0.25, color = "grey60") +
  rasterize_layer(geom_point(color = "#1565C0", size = 0.45, alpha = 0.35, shape = 16)) +
  annotate("text", x = 0.04, y = 0.97, hjust = 0, size = GEOM_TEXT_6PT, color = "grey20",
           label = sprintf("r = %.2f\nSuSiE>0.5: %d\nABF>0.5: %d", rA, nA_susie, nA_abf)) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  labs(x = "coloc.abf PP.H4", y = "SuSiE-coloc PP.H4") +
  theme_masld(base_size = 6)

# ---- B. SuSiE vs CARMA PIP ----
cf <- fread(file.path(BASE, "GWAS/finemapping/results/combined_finemapping.csv"),
            select = c("susie_pip", "carma_pip", "both_in_cs", "either_in_cs", "carma_outlier"))
both <- cf[!is.na(susie_pip) & !is.na(carma_pip)]
rB <- cor(both$susie_pip, both$carma_pip)
fracCS <- sum(cf$both_in_cs, na.rm = TRUE) / sum(cf$either_in_cs, na.rm = TRUE)
outl <- 100 * mean(cf$carma_outlier, na.rm = TRUE)
pB <- ggplot(both[susie_pip > 0.01 | carma_pip > 0.01], aes(carma_pip, susie_pip)) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.3, color = "grey70") +
  rasterize_layer(geom_point(color = "#00695C", size = 0.4, alpha = 0.25, shape = 16)) +
  annotate("text", x = 0.04, y = 0.97, hjust = 0, size = GEOM_TEXT_6PT, color = "grey20",
           label = sprintf("r = %.2f\nboth-in-CS: %.0f%%\nCARMA outlier: %.1f%%",
                           rB, 100 * fracCS, outl)) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  labs(x = "CARMA PIP", y = "SuSiE PIP") +
  theme_masld(base_size = 6)

# ---- C. LD-reference PP.H4 concordance ----
ld <- fread(file.path(BASE, "GWAS/finemapping/results/ld_panel_comparison/ld_panel_concordance.csv"))
rfun <- function(a, b) { ok <- !is.na(ld[[a]]) & !is.na(ld[[b]]); cor(ld[[a]][ok], ld[[b]][ok]) }
ldd <- data.table(
  pair = c("1KG vs\nTOP-LD", "1KG vs\nUKBB"),
  r = c(rfun("PP.H4.susie_1KG", "PP.H4.susie_TOP-LD"),
        rfun("PP.H4.susie_1KG", "PP.H4.susie_UKBB_v1")))
pC <- ggplot(ldd, aes(pair, r)) +
  geom_col(fill = "#9575CD", width = 0.6) +
  geom_text(aes(label = sprintf("%.2f", r)), vjust = -0.4, size = GEOM_TEXT_6PT, fontface = "plain") +
  scale_y_continuous(limits = c(0, 1.05), breaks = c(0, 0.5, 1), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "PP.H4 concordance (r)") +
  theme_masld(base_size = 6) +
  theme(axis.text.x = element_text(size = 6))

p <- (pA | pB | pC) + plot_layout(widths = c(1, 1, 0.8))
message("[caption] Colocalization is robust to method and LD-reference choice. A: coloc.abf vs SuSiE-coloc method. B: SuSiE vs CARMA fine-mapping. C: LD reference panel.")

save_fig(p, file.path(PANEL_DIR, "FigS2F_coloc_method_concordance.pdf"),
         width = fig_full_width * 0.66, height = 2.4)

fwrite(data.table(
  metric = c("ABF_vs_SuSiE_r", "SuSiE>0.5", "ABF>0.5", "SuSiExCARMA_r",
             "both_in_CS_frac", "CARMA_outlier_pct", "LD_1KG_TOPLD_r", "LD_1KG_UKBB_r"),
  value = round(c(rA, nA_susie, nA_abf, rB, fracCS, outl/100, ldd$r[1], ldd$r[2]), 4)),
  file.path(PANEL_DIR, "FigS2F_coloc_method_concordance_source.csv"))
cat(sprintf("[figS method] ABFvsSuSiE r=%.2f | SuSiExCARMA r=%.2f outl=%.1f%% | LD r=%.2f/%.2f\n",
            rA, rB, outl, ldd$r[1], ldd$r[2]))
