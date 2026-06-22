#!/usr/bin/env Rscript
# KEY MESSAGE (HONEST NEGATIVE — caQTL does NOT functionally validate the
# motif-disruption predictions): in Currin 2025 138-donor bulk-liver caQTL,
# motif-disruption direction agrees with caQTL effect direction in only
# 1,201 / 2,595 tested variants = 46.3% (i.e. 53.7% ANTI-concordant; binomial
# p = 1.6e-4 vs the 50% null — significantly BELOW chance, not above). The
# per-TF concordance distribution is centred exactly on the 50% null (median
# 50.0% across 328 TFs with n>=4). A few disease-regulon TFs sit above 50% on
# small n (HNF4A 83% n=12, THRB 80% n=10, RORA 67% n=9) but are balanced by an
# equal number below (ESR1 33%, FOXA1/FOXA2/ETS2 25%). Only 2 variants reach
# full 4-way concordance (motif+eQTL+DA+caQTL) and BOTH fall in ZNF701 (not a
# disease gene); 0 disease-regulon TF has a 4-way variant. This is why the
# variant->motif predictions are reported as PUTATIVE and the epigenetic panels
# were CUT from main Fig 4. (Prior version of this script claimed "80-83%
# concordant / 47 four-way" — both FALSE vs the data; corrected 2026-06-19.)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUT_PDF <- file.path(FIGS05_DIR, "figS05_scatac_caqtl_concordance.pdf")

GA_DIR <- file.path(BASE, "GWAS/finemapping/results/gwas_atac")
var_dt <- fread(file.path(GA_DIR, "currin_caqtl_concordance.csv"))
tf_dt  <- fread(file.path(GA_DIR, "currin_caqtl_per_tf.csv"))

# ---- Color palette (Liang style) -----------------------------------------
MAG  <- "#C9265E"
BLUE <- "#1565C0"
GRAY <- "#9E9E9E"
DARK <- "#37474F"

# ---- Panel a: variant-level motif x caQTL effect scatter -----------------
# Honest scatter: if motif disruption validated caQTL effects, points would
# cluster in the concordant (++ / --) quadrants. They DON'T — 46.3% concordant,
# i.e. an even (slightly anti-concordant) split. No cherry-picked TF labels.
v <- var_dt[currin_caqtl_tested == TRUE & !is.na(alleleDiff) & !is.na(currin_caqtl_beta)]
v[, concordant := concord_motif_caqtl == TRUE]
n_conc <- sum(v$concordant, na.rm = TRUE); n_tot <- nrow(v)
pct_conc <- 100 * n_conc / n_tot
binom_p  <- binom.test(n_conc, n_tot, 0.5)$p.value
v[, conc_lab := factor(fifelse(concordant, "Concordant", "Discordant"),
                       levels = c("Concordant", "Discordant"))]
setorder(v, -conc_lab)  # plot concordant on top, neither dominates

p_a <- ggplot(v, aes(x = alleleDiff, y = currin_caqtl_beta)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "gray70", linewidth = 0.25) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "gray70", linewidth = 0.25) +
  geom_point(aes(colour = conc_lab), shape = 16, size = 0.7, alpha = 0.55) +
  annotate("text", x = -max(abs(v$alleleDiff)) * 0.98,
           y = max(v$currin_caqtl_beta, na.rm = TRUE) * 0.96,
           hjust = 0, vjust = 1, size = PUB_GEOM_TEXT, colour = "black",
           label = sprintf("%.1f%% concordant\n(n = %s; p = %.1e vs 50%%)",
                           pct_conc, format(n_tot, big.mark = ","), binom_p)) +
  scale_colour_manual(values = c("Concordant" = BLUE, "Discordant" = GRAY), name = NULL) +
  scale_x_continuous(name = expression(bold("Motif disruption (alleleDiff)")),
                     limits = c(-1, 1) * max(abs(v$alleleDiff)),
                     expand = expansion(mult = 0.04)) +
  scale_y_continuous(name = expression(bold(paste("caQTL ", beta))),
                     expand = expansion(mult = 0.05)) +
  labs(tag = "a", colour = NULL) +
  theme_masld(base_size = 7) +
  theme_pub() +
  theme(panel.grid = element_blank(),
        legend.position = "bottom",
        legend.key.size = unit(0.25, "cm"),
        plot.tag = element_text(size = 11, face = "bold"),
        plot.tag.position = c(0.02, 0.97))

# ---- Panel b: per-TF concordance distribution (NO lollipop) ---------------
# Histogram of per-TF motif x caQTL concordance for all adequately-tested TFs.
# It is centred on the 50% null (median 50.0%) — no global enrichment. The few
# disease-regulon TFs of interest are marked as a rug to show they STRADDLE the
# null (some above, some below), not a coherent validation.
tfh <- tf_dt[n_caqtl_tested >= 4 & !is.na(pct_motif_caqtl_agree)]
med_tf <- median(tfh$pct_motif_caqtl_agree)
dz_rug <- tf_dt[any_disease_regulon == TRUE & n_caqtl_tested >= 5 &
                !is.na(pct_motif_caqtl_agree),
                .(tf_name, pct_motif_caqtl_agree, n_caqtl_tested)]

p_b <- ggplot(tfh, aes(x = pct_motif_caqtl_agree)) +
  geom_histogram(binwidth = 10, boundary = 0, fill = GRAY,
                 colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = 50, linetype = "dashed", colour = DARK, linewidth = 0.4) +
  geom_rug(data = dz_rug, aes(x = pct_motif_caqtl_agree),
           colour = MAG, linewidth = 0.5, length = unit(0.06, "npc"), inherit.aes = FALSE) +
  annotate("text", x = 50, y = Inf, vjust = 1.4, hjust = -0.06,
           size = PUB_GEOM_TEXT, colour = "black",
           label = sprintf("median %.0f%%\n(%d TFs)", med_tf, nrow(tfh))) +
  scale_x_continuous(name = expression(bold("Motif x caQTL concordance (%) per TF")),
                     breaks = c(0, 25, 50, 75, 100)) +
  scale_y_continuous(name = "TFs", expand = expansion(mult = c(0, 0.08))) +
  labs(tag = "b") +
  theme_masld(base_size = 7) +
  theme_pub() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        plot.tag = element_text(size = 11, face = "bold"),
        plot.tag.position = c(0.02, 0.97))

# ---- Compose (2 panels; the old one-bar "panel c" was a single ZNF701 datum,
#      folded into the caption per no-single-point-panel rule) ---------------
combo <- p_a + p_b + plot_layout(widths = c(2.6, 2.4))

ggsave(OUT_PDF, combo, width = 6.2, height = 2.9, device = cairo_pdf)
cat(sprintf("[saved] %s\n", OUT_PDF))

# ---- Caption / provenance to stdout (NOT on the panel) -------------------
n_4way <- sum(v$concord_4of4 == TRUE, na.rm = TRUE)
tf_4way <- tf_dt[n_4of4 >= 1]$tf_name
message(strrep("=", 78))
message("FIGURE LEGEND (paste into manuscript; stats live here, not on the panel)")
message(strrep("=", 78))
message(sprintf(
"Currin 2025 bulk-liver caQTL (138 donors) does NOT functionally validate the
motif-disruption predictions. (a) Motif-disruption direction (alleleDiff) vs
caQTL effect direction for %s tested fine-mapped variants: only %d are
direction-concordant (%.1f%%; binomial p = %.1e vs the 50%% null -- significantly
BELOW chance). (b) Per-TF concordance is centred exactly on the 50%% null
(median %.0f%% across %d TFs with n>=4); the disease-regulon TFs of interest
(magenta rug) STRADDLE it -- HNF4A 83%% (n=12), THRB 80%% (n=10), RORA 67%% (n=9)
above, balanced by ESR1 33%%, FOXA1/FOXA2/ETS2 25%% below.",
  format(n_tot, big.mark = ","), n_conc, pct_conc, binom_p, med_tf, nrow(tfh)))
message(sprintf(
"Only %d variant(s) reach full 4-way concordance (motif+eQTL+DA+caQTL), both in %s
-- NOT a disease gene; 0 disease-regulon TF has a 4-way variant. The variant->motif
disruption layer is therefore reported as PUTATIVE throughout.",
  n_4way, paste(tf_4way, collapse = ", ")))
message(strrep("=", 78))
cat(sprintf("  variant-level concordance: %d/%d = %.1f%% (p=%.2e); 4-way variants: %d (%s)\n",
            n_conc, n_tot, pct_conc, binom_p, n_4way, paste(tf_4way, collapse = ", ")))
