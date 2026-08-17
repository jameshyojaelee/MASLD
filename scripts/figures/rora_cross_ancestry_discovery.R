#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# Fig 4 RORA — the flagship cross-ancestry DISCOVERY. RORA is absent from the
# canonical NAFLD GWAS catalog (PNPLA3/TM6SF2/MBOAT7/GCKR/HSD17B13), yet its
# liver-cis-eQTL colocalizes with liver-enzyme GWAS at SuSiE PP.H4 > 0.99 in
# THREE independent ancestries (EUR/EAS/SAS), and the lead credible-set variants
# converge on the SAME ~5 kb region of chr15 — a shared causal signal, not three
# coincidences. Nominated de novo by the unsupervised score (no drug labels),
# then independently validated by entry into Phase 1 (TB-840 RORα agonist).
#
# DESIGN (2026-06-19): a LOCUS-CONCORDANCE diagram — 3 ancestry tracks on a shared
# chr15 genomic axis, each marking its lead credible-set variant (point shaded by
# fine-mapping PIP) with the COLOC PP.H4 annotated; the two lead-variant positions
# drawn as guide lines to show convergence. NO lollipop, NO saturated bars
# (see memory/feedback-no-lollipop); PIP is the varying, informative axis.
#
# HARD GATE (anti-FADS2): COLOC + top_snp + PIP read ONLY from the canonical
#   GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv (top_snp,
#   top_snp_PP, PP.H4.susie). bulk from canonical_deg_results.csv (C2). NEVER the
#   atlas *_coloc_pp4 convenience columns.
#
# Output: figures/main/fig5_molecular_context/fig4c_rora_cross_ancestry_discovery.pdf
# Env:    rnaseq
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

# ── RETIRED 2026-07-02 ────────────────────────────────────────────────────────
# This standalone panel was CUT from Fig 4 and its PDF deleted. RORA is still
# represented in Fig 4 by fig4c_atac_rora_motif (the regulatory-mechanism panel)
# and in the supplement by figS_rora_case_study. This script no longer emits.
# Body below retained for provenance only; remove the guard to regenerate.
message("[RETIRED 2026-07-02] fig4c_rora_cross_ancestry_discovery was cut from Fig 4; not regenerated.")
quit(save = "no", status = 0)

# ── Cross-ancestry COLOC + lead credible-set variant per ancestry ────────────
co <- fread(file.path(BASE, "GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv"))
rora <- co[gene == "RORA"]
anc_of <- function(g) fifelse(grepl("BBJ", g), "EAS",
                       fifelse(grepl("PanUKBB_CSA", g), "SAS",
                       fifelse(grepl("PanUKBB_AFR", g), "AFR", "EUR")))
rora[, ancestry := anc_of(gwas_name)]
rora[, trait := fifelse(grepl("GGT", gwas_name), "GGT",
                fifelse(grepl("ALT", gwas_name), "ALT",
                fifelse(grepl("AST", gwas_name), "AST", "other")))]
# GGT is the trait that colocalizes across all three ancestries — use it as the
# common cross-ancestry comparison; take the best SuSiE-grade row per ancestry.
ggt <- rora[trait == "GGT" & !is.na(PP.H4.susie) & PP.H4.susie > 0.9]
ggt[, pos := as.numeric(sub(".*:", "", top_snp))]
lead <- ggt[, .SD[which.max(PP.H4.susie)], by = ancestry][
  ancestry %in% c("EUR", "EAS", "SAS")]
lead[, anc_lab := factor(c(EUR = "European (UKBB)", EAS = "East Asian (BBJ)",
                           SAS = "South Asian (Pan-UKBB)")[ancestry],
       levels = c("South Asian (Pan-UKBB)", "East Asian (BBJ)", "European (UKBB)"))]
lead[, y := as.integer(anc_lab)]
lead[, pip := as.numeric(top_snp_PP)]

# ── bulk + mouse + drug footer values ────────────────────────────────────────
deg <- fread(file.path(INT_RESULTS, "canonical_deg_results.csv"))
bulk_lfc <- deg[symbol == "RORA", logFC][1]; bulk_padj <- deg[symbol == "RORA", padj][1]

# ── Plot: chr15 locus axis, 3 ancestry tracks, lead variant shaded by PIP ────
lead_pos <- sort(unique(lead$pos))                       # the (≤2) lead variants
xpad <- 4000
xr <- range(lead$pos) + c(-xpad, xpad * 3.2)   # extra right pad so right-anchored PP.H4 labels clear the EUR variant's PIP label

p <- ggplot(lead) +
  # converging lead-variant guide lines
  geom_vline(xintercept = lead_pos, linetype = "22", color = "grey75", linewidth = 0.3) +
  # per-ancestry locus track
  geom_segment(aes(x = xr[1], xend = xr[2], y = y, yend = y), color = "grey88", linewidth = 0.5) +
  # lead credible-set variant, shaded by PIP
  geom_point(aes(x = pos, y = y, fill = pip), shape = 21, size = 4.2, stroke = 0.3, color = "white") +
  geom_text(aes(x = pos, y = y + 0.26, label = sprintf("PIP %.2f", pip)),
            size = GEOM_TEXT_6PT, color = "black") +
  # COLOC PP.H4 annotated at the right (BLACK text — no colored fonts)
  geom_text(aes(x = xr[2], y = y, label = sprintf("PP.H4 %.3f", PP.H4.susie)),
            hjust = 1, vjust = -0.8, size = GEOM_TEXT_6PT, color = "black") +
  scale_fill_gradientn(colours = colorRampPalette(c("#F4A6C2", "#C9265E", "#7A1140"))(64),
                       limits = c(0.4, 1), name = "PIP",
                       guide = guide_colorbar(barheight = 2.2, barwidth = 0.35, ticks = FALSE)) +
  scale_y_continuous(breaks = lead$y, labels = lead$anc_lab,
                     limits = c(0.62, 3.55), expand = c(0, 0)) +
  scale_x_continuous(labels = function(z) sprintf("%.3f", z / 1e6),
                     breaks = lead_pos) +
  labs(x = "chr15 position (Mb)", y = NULL) +
  theme_masld_compact() +
  theme(axis.text.y = element_text(size = 6, color = "black"),
        legend.position = "right", legend.direction = "vertical",
        legend.key.height = unit(0.28, "cm"), legend.key.width = unit(0.22, "cm"),
        panel.grid = element_blank(), plot.margin = margin(3, 6, 3, 3))

out <- file.path(FIG4_DIR, "panels", "fig4c_rora_cross_ancestry_discovery.pdf")
save_fig(p, out, width = fig_half_width - 0.15, height = 1.7)
message("Saved: ", out)
print(lead[, .(ancestry, trait, top_snp, pos, pip = round(pip,3), PP.H4.susie = round(PP.H4.susie,3))])
message(sprintf("[rora] lead variants %s bp apart; bulk logFC %.2f padj %.1e", diff(range(lead$pos)), bulk_lfc, bulk_padj))
