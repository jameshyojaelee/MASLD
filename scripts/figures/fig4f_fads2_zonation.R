#!/usr/bin/env Rscript
# ============================================================================
# ⛔ SUPERSEDED / RETIRED 2026-06-18 — DO NOT USE.
# This script is WRONG: it claims FADS2 is East-Asian-specific by reading the
# atlas `coloc_susie_best_pp4` column, which is method-inconsistent/corrupted
# (see memory/reference-coloc-canonical-source-pitfall). The canonical
# susie_coloc_all_gwas.csv shows FADS2 is CROSS-ANCESTRY: EUR UKBB-GGT SuSiE
# 0.909 + EAS BBJ-ALT 0.948 (+ BBJ-AST 0.901); EAS BBJ-GGT = 0.003. The correct
# panel is `fads2_cross_ancestry.R` (-> _supp/fads2_cross_ancestry.pdf). Output
# below was repointed to a _superseded_ name so this can never overwrite the
# canonical panel. Kept only for provenance.
# ============================================================================
# KEY MESSAGE (RETIRED — INCORRECT): FADS2 colocalises with East-Asian (BBJ) liver-enzyme GWAS
# (SuSiE PP.H4 = 0.948, BBJ ALT) but is invisible to European, African, and
# South-Asian scans — cross-ancestry GWAS reveals disease biology a single-
# ancestry analysis misses.
#
# One compact panel: best SuSiE-COLOC PP.H4 per ancestry. Values are read from
# the atlas (so they trace to disk); the East-Asian bar is the only one to clear
# the 0.5 threshold. Supplementary cross-ancestry case for Figure 4.
#
# (A prior zonation-landscape scatter was dropped — it highlighted CYP3A4, not
#  FADS2, duplicating fig4f_cyp3a4_zonation.R. The cross-ancestry contrast is the
#  FADS2 story.)
#
# Output: figures/main/fig4_validation/fads2_cross_ancestry.pdf
# Env:    rnaseq

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

# ── Data ──────────────────────────────────────────────────────────────────────
atlas <- read.csv(file.path(BASE,
  "Analysis/Spatial/results/integration/multi_evidence_atlas_with_spatial.csv"),
  stringsAsFactors = FALSE)
f <- atlas %>% filter(human_symbol == "FADS2")
stopifnot(nrow(f) == 1)

# Best liver-enzyme PP.H4 per ancestry (max across ALT/AST/GGT/PDFF where present).
best <- function(...) {
  v <- suppressWarnings(as.numeric(c(...)))
  if (all(is.na(v))) return(NA_real_)
  max(v, na.rm = TRUE)
}
eur <- best(f$ukbb_alt_coloc_pp4, f$ast_coloc_pp4, f$ggt_coloc_pp4,
            f$pdff_coloc_pp4, f$coloc_best_pp4_polyfun)
eas <- best(f$coloc_susie_best_pp4)                       # 0.948 (BBJ ALT, SuSiE)
sas <- best(f$panukbb_csa_alt_coloc_pp4, f$panukbb_csa_ast_coloc_pp4,
            f$panukbb_csa_ggt_coloc_pp4)
afr <- best(f$panukbb_afr_alt_coloc_pp4, f$panukbb_afr_ast_coloc_pp4,
            f$panukbb_afr_ggt_coloc_pp4)

anc <- data.frame(
  ancestry = c("European\n(UKBB)", "East Asian\n(BBJ)",
               "South Asian\n(Pan-UKBB)", "African\n(Pan-UKBB)"),
  pp4      = c(eur, eas, sas, afr),
  is_eas   = c(FALSE, TRUE, FALSE, FALSE)
)
anc$ancestry <- factor(anc$ancestry, levels = anc$ancestry)
anc$fill     <- ifelse(anc$is_eas, masld_colors[["up"]], "#BDBDBD")

message(sprintf("[fads2] best PP.H4 — EUR %.3f | EAS %.3f | SAS %.3f | AFR %.3f",
                eur, eas, sas, afr))

# ── Panel: cross-ancestry COLOC ──────────────────────────────────────────────
p <- ggplot(anc, aes(x = ancestry, y = pp4, fill = fill)) +
  geom_col(width = 0.62) +
  geom_hline(yintercept = 0.5, linewidth = 0.3, linetype = "dashed", color = "gray60") +
  scale_fill_identity() +
  scale_y_continuous(limits = c(0, 1.02), breaks = c(0, 0.5, 1.0),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "SuSiE-COLOC PP.H4") +
  theme_masld() + theme_pub() +
  theme(axis.text.x = element_text(size = PUB_AXIS_TEXT))

out <- file.path(FIG4_DIR, "_supp", "_superseded_fig4f_fads2_zonation.pdf")
pdf(out, width = fig_half_width, height = fig_half_width * 0.72, useDingbats = FALSE)
print(p)
dev.off()
message("Saved: ", out)
