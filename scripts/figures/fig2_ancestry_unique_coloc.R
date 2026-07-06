#!/usr/bin/env Rscript
# fig2_ancestry_unique_coloc.R  (2026-06-25)  — Fig 2G
# Ancestry SPECIFICITY of colocalizing genes: how many MASLD effector genes are
# discovered ONLY in non-European ancestries (and would be missed by a EUR-only
# analysis). Purely GWAS x eQTL colocalization — no RNA-seq disease signal, so it
# is safe to show before the transcriptomic atlas is introduced in Fig 3.
#
# Two linked elements in one panel (the EUR-only count is ~75x the smallest non-EUR
# bin, so a single shared scale crushes the breakdown):
#   TOP  : a 100%-stacked proportion strip of all colocalizing genes
#          (EUR-only / shared / non-EUR-unique) — the whole-cohort context.
#   BELOW: the non-EUR-UNIQUE genes broken out by ancestry on their OWN scale
#          (EAS / AFR / SAS / AMR / multiple non-EUR), now fully readable.
#
# Colocalization = best PP.H4 > 0.5 (SuSiE with ABF fallback — SuSiE cannot converge
# on the sparse non-EUR PanUKBB/BBJ GWAS, so ABF is the only available coloc there).
# Reuses the ancestry map + pp4_best logic + palette from fig3b_ancestry_coloc_bars.R
# (which shows per-ancestry TOTALS — a different cut; this panel is the specificity bins).
#
# Out: figures/main/fig2_genetics/panels/Fig2G_ancestry_unique_coloc.pdf (+ source CSV)
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
PANEL_DIR <- file.path(FIG3_DIR, "panels")

# ---- ancestry map: registry-driven (gwas_ancestry() in load_figure_data.R) -----
# Retires the hardcoded name-sets, which had no AMR bin and routed every MVP stratum
# into EUR once the portfolio grew 23 -> 50 GWAS. NOTE (2026-07-05): this .R is a
# LEGACY DUPLICATE of fig2_ancestry_unique_coloc.py (both write the same
# Fig2G_ancestry_unique_coloc.pdf). The .py family is CANONICAL (it also produces the
# GWS/suggestive companions). Do NOT run this .R over the .py — the ancestry fix here
# only guarantees it can't emit a MISLABELED PDF if someone does.
sc <- fread(file.path(BASE, "GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv"))
sc[, ancestry := as.character(gwas_ancestry(gwas_name))]
sc <- sc[!is.na(ancestry)]
sc[, pp4_best := pmax(PP.H4.susie, PP.H4.abf, na.rm = TRUE)]   # best of SuSiE/ABF (user-chosen)
sc[is.infinite(pp4_best), pp4_best := NA_real_]
sig <- sc[!is.na(pp4_best) & pp4_best > 0.5]

# ---- per-gene ancestry set -> specificity bin ------------------------------
g <- sig[, .(has_eur = any(ancestry == "EUR"),
             n_nonEur = uniqueN(ancestry[ancestry != "EUR"]),
             ancs = paste(sort(unique(ancestry)), collapse = "+")), by = gene]
g[, bin := fcase(
  ancs == "EUR",            "EUR only",
  has_eur & n_nonEur >= 1,  "EUR + non-EUR (shared)",
  ancs == "EAS",            "EAS only",
  ancs == "AFR",            "AFR only",
  ancs == "SAS",            "SAS only",
  ancs == "AMR",            "AMR only",
  !has_eur & n_nonEur >= 2, "multiple non-EUR",
  default = NA_character_)]
cnt <- function(b) g[bin == b, .N]
n_eur <- cnt("EUR only"); n_shared <- cnt("EUR + non-EUR (shared)")
eas <- cnt("EAS only"); afr <- cnt("AFR only"); sas <- cnt("SAS only")
amr <- cnt("AMR only"); multi <- cnt("multiple non-EUR")
n_unique <- eas + afr + sas + amr + multi; n_total <- nrow(g)

# ---- TOP: 100%-stacked proportion strip (whole-cohort context) -------------
strip <- data.table(
  cat = factor(c("EUR-only","EUR + non-EUR (shared)","non-EUR-unique"),
               levels = c("EUR-only","EUR + non-EUR (shared)","non-EUR-unique")),
  n   = c(n_eur, n_shared, n_unique))
strip[, pct := 100 * n / sum(n)]
strip[, txt := fifelse(cat == "EUR-only", "grey15", "white")]
p_strip <- ggplot(strip, aes(x = n, y = "all", fill = cat)) +
  geom_col(position = "fill", width = 0.55, color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.0f%%", pct), color = cat),
            position = position_fill(vjust = 0.5), fontface = "bold", size = 2.7) +
  scale_fill_manual(values = c("EUR-only" = "#BDBDBD",
                               "EUR + non-EUR (shared)" = "#78909C",
                               "non-EUR-unique" = "#283593"), name = NULL) +
  scale_color_manual(values = c("EUR-only" = "grey15",
                                "EUR + non-EUR (shared)" = "white",
                                "non-EUR-unique" = "white"), guide = "none") +
  scale_x_continuous(expand = expansion(0)) +
  labs(x = NULL, y = NULL) +
  theme_masld() + theme_pub() +
  theme(legend.position = "top", legend.text = element_text(size = 7, color = "black"),
        legend.key.size = unit(0.30, "cm"),
        axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank())

# ---- BELOW: non-EUR-unique genes by ancestry, own scale --------------------
brk <- data.table(
  ancb = factor(c("EAS only","AFR only","SAS only","AMR only","multiple non-EUR"),
                levels = rev(c("EAS only","AFR only","SAS only","AMR only","multiple non-EUR"))),
  n    = c(eas, afr, sas, amr, multi))
# ancestry hexes mirror ANCESTRY_COLORS (load_figure_data.R): EAS red / AFR green /
# SAS purple / AMR orange; multi = neutral grey.
acol <- c("EAS only" = "#C44E52", "AFR only" = "#55A868", "SAS only" = "#8172B3",
          "AMR only" = "#DD8452", "multiple non-EUR" = "#455A64")
p_break <- ggplot(brk, aes(x = n, y = ancb, fill = ancb)) +
  geom_col(width = 0.66) +
  geom_text(aes(label = n), hjust = -0.3, size = 2.9, fontface = "bold", color = "black") +
  scale_fill_manual(values = acol, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.14))) +
  labs(x = sprintf("non-EUR-unique colocalizing genes (n=%d of %d)", n_unique, n_total), y = NULL) +
  theme_masld() + theme_pub() +
  theme(axis.text = element_text(color = "black"),
        axis.title.x = element_text(color = "black"))

p <- p_strip / p_break + plot_layout(heights = c(1, 2.4)) +
  plot_annotation(title = "Colocalization by ancestry specificity",
                  theme = theme(plot.title = element_text(size = 9, face = "bold")))

save_fig(p, file.path(PANEL_DIR, "Fig2G_ancestry_unique_coloc.pdf"),
         width = fig_col_width * 1.15, height = 2.9)

fwrite(g[, .N, by = bin][order(-N)], file.path(PANEL_DIR, "Fig2G_ancestry_unique_coloc_source.csv"))
cat(sprintf("[fig2G/.R legacy] wrote Fig2G_ancestry_unique_coloc.pdf | %d colocalizing genes; %d non-EUR-unique (EAS %d/AFR %d/SAS %d/AMR %d/multi %d)\n",
            n_total, n_unique, eas, afr, sas, amr, multi))
message(sprintf("CAPTION (Fig2G): Ancestry specificity of %d colocalizing MASLD effector genes (best PP.H4 > 0.5; SuSiE with ABF fallback; 50-GWAS portfolio incl. MVP). %d (%.0f%%) colocalize ONLY in non-European ancestries (EAS %d / AFR %d / SAS %d / AMR %d / multiple non-EUR %d) and would be missed by a European-only analysis; %d (%.0f%%) are shared with EUR and %d (%.0f%%) are EUR-only.",
                n_total, n_unique, 100*n_unique/n_total, eas, afr, sas, amr, multi,
                n_shared, 100*n_shared/n_total, n_eur, 100*n_eur/n_total))
