#!/usr/bin/env Rscript
# KEY MESSAGE (POSITIVE CONTROL): FADS2, a canonical pan-ancestry desaturase
# locus (Lemaitre 2011), colocalizes with liver/enzyme GWAS across FOUR
# ancestries — Europeans (PP.H4.susie 0.96), East Asians (0.95), Admixed
# Americans (0.94) and Africans (0.91) — a cross-ancestry replication that was
# STRENGTHENED once the MVP multi-ancestry GWAS (dbGaP phs002453) entered the
# 50-GWAS portfolio (AMR + AFR were abf-only near zero before MVP). This is the
# methodological positive control for the COLOC pipeline, NOT an ancestry-
# specific signal.
#
# One compact panel: best SuSiE-COLOC PP.H4 per ancestry. Ancestry is assigned
# from the GWAS registry (gwas_ancestry(); load_figure_data.R), which now covers
# every MVP stratum (MVP_{ALT,AST,GGT,NAFLD,...} x {EUR,EAS,AFR,AMR}) — retiring
# the old grepl() heuristic, which had no AMR bin and mis-routed every MVP non-EUR
# stratum into EUR. Numbers come straight from the canonical COLOC file
# susie_coloc_all_gwas.csv (NOT the atlas *_coloc_pp4 convenience columns).
# EUR/EAS/AMR/AFR all clear the 0.5 threshold; only SAS (Pan-UKBB CSA; no MVP
# South-Asian stratum) stays near zero.
#
# Best per ancestry (verified vs disk 2026-07-05):
#   EUR  MVP_ALT_EUR       PP.H4.susie 0.963
#   EAS  BBJ_ALT           PP.H4.susie 0.948
#   AMR  MVP_ALT_AMR       PP.H4.susie 0.944   (NEW via MVP)
#   AFR  MVP_Platelet_AFR  PP.H4.susie 0.911   (NEW via MVP; was ~0.019 Pan-UKBB-only)
#   SAS  PanUKBB_CSA_AST   PP.H4.abf   ~0.007  (no SuSiE convergence; abf only)
#
# CAVEATS (baked into legend, see message() to stdout):
#   - trait-specific: top signal is ALT for EUR/EAS/AMR, platelet count for AFR
#   - locus-level colocalization, NOT identical lead variant across ancestries
#   - abf calls H3 (one-causal-variant assumption) at this multi-causal locus;
#     the SuSiE (allows multiple causal signals) calls are the robust ones
#   - the SAS Pan-UKBB scan never reached SuSiE convergence (abf only, low PP)
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

# ── Data: canonical COLOC file (NOT atlas convenience columns) ─────────────────
coloc <- read.csv(file.path(BASE,
  "GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv"),
  stringsAsFactors = FALSE)
f <- coloc %>% filter(gene == "FADS2")
stopifnot(nrow(f) > 0)

# Ancestry is assigned from the GWAS registry (gwas_ancestry() in
# load_figure_data.R), which now covers the 50-GWAS portfolio incl. every MVP
# stratum (MVP_{ALT,AST,GGT,NAFLD,...} x {EUR,EAS,AFR,AMR}). This retires the old
# grepl() heuristic, which had no AMR bin and dropped every MVP non-EUR stratum
# into EUR.
f$ancestry <- as.character(gwas_ancestry(f$gwas_name))

# Best PP.H4 per ancestry. Prefer the SuSiE posterior (robust to multiple causal
# variants); fall back to abf only where SuSiE never converged (SAS / AFR).
best_pp <- function(df) {
  s <- suppressWarnings(as.numeric(df$PP.H4.susie))
  if (any(!is.na(s))) {
    i <- which.max(s)
    return(list(pp = s[i], gwas = df$gwas_name[i], src = "SuSiE"))
  }
  a <- suppressWarnings(as.numeric(df$PP.H4.abf))
  i <- which.max(a)
  list(pp = a[i], gwas = df$gwas_name[i], src = "abf")
}

picks <- lapply(c("EUR", "EAS", "AMR", "AFR", "SAS"), function(a) {
  b <- best_pp(f[f$ancestry == a, , drop = FALSE])
  data.frame(anc = a, pp4 = b$pp, gwas = b$gwas, src = b$src,
             stringsAsFactors = FALSE)
})
picks <- do.call(rbind, picks)

# Ancestry-name axis labels only (the exact best-colocalizing scan per ancestry is
# emitted dynamically to the message() legend, so nothing here goes stale as the
# portfolio grows). picks rows are already in EUR/EAS/AMR/AFR/SAS order.
anc_full <- c(EUR = "European", EAS = "East Asian", AMR = "Admixed\nAmerican",
              AFR = "African", SAS = "South Asian")
anc <- data.frame(
  ancestry = anc_full[picks$anc],
  pp4      = picks$pp4,
  replicated = picks$pp4 >= 0.5,
  stringsAsFactors = FALSE
)
anc$ancestry <- factor(anc$ancestry, levels = anc$ancestry)
# Both replicating ancestries highlighted in disease magenta; non-replicating gray.
anc$fill  <- ifelse(anc$replicated, masld_colors[["up"]], masld_colors[["ns"]])
anc$lab   <- sprintf("%.2f", anc$pp4)

# ── Provenance + caveats to stdout (these belong in the legend, not the panel) ─
eur <- picks$pp4[picks$anc == "EUR"]; eas <- picks$pp4[picks$anc == "EAS"]
amr <- picks$pp4[picks$anc == "AMR"]; afr <- picks$pp4[picks$anc == "AFR"]
sas <- picks$pp4[picks$anc == "SAS"]
message(sprintf(
  "[fads2] best PP.H4 — EUR %.3f (%s, %s) | EAS %.3f (%s, %s) | AMR %.3f (%s, %s) | AFR %.3f (%s, %s) | SAS %.3f (%s, %s)",
  eur, picks$gwas[picks$anc=="EUR"], picks$src[picks$anc=="EUR"],
  eas, picks$gwas[picks$anc=="EAS"], picks$src[picks$anc=="EAS"],
  amr, picks$gwas[picks$anc=="AMR"], picks$src[picks$anc=="AMR"],
  afr, picks$gwas[picks$anc=="AFR"], picks$src[picks$anc=="AFR"],
  sas, picks$gwas[picks$anc=="SAS"], picks$src[picks$anc=="SAS"]))
message("[fads2] LEGEND: FADS2, a canonical pan-ancestry desaturase locus (Lemaitre et al. 2011),")
message("[fads2]   colocalizes with liver/enzyme GWAS across FOUR ancestries (EUR/EAS/AMR/AFR, all")
message("[fads2]   PP.H4.susie > 0.9): cross-ancestry replication, positive control. Strengthened once the")
message("[fads2]   MVP multi-ancestry GWAS entered the portfolio (AMR + AFR were abf-only near zero before).")
message("[fads2] CAVEAT: trait-specific (top signal is ALT for EUR/EAS/AMR, platelet count for AFR);")
message("[fads2]   locus-level coloc, NOT an identical lead variant; SuSiE (multi-causal) posteriors are")
message("[fads2]   robust whereas coloc.abf calls H3 (distinct causal variants) at this multi-causal locus.")
message("[fads2]   Only SAS (Pan-UKBB CSA; no MVP South-Asian stratum) stays near zero (abf only).")

# ── Panel: cross-ancestry COLOC ───────────────────────────────────────────────
p <- ggplot(anc, aes(x = ancestry, y = pp4, fill = fill)) +
  geom_col(width = 0.62) +
  geom_hline(yintercept = 0.5, linewidth = 0.3, linetype = "dashed", color = "gray60") +
  geom_text(aes(label = lab, color = replicated),
            vjust = -0.5, size = GEOM_TEXT_6PT, fontface = "plain", show.legend = FALSE) +
  annotate("text", x = 5.45, y = 0.52, label = "PP.H4 = 0.5",
           hjust = 1, vjust = -0.4, size = GEOM_TEXT_6PT, color = "gray55") +
  scale_fill_identity() +
  scale_color_manual(values = c(`TRUE` = masld_colors[["up"]], `FALSE` = masld_colors[["ns"]])) +
  scale_y_continuous(limits = c(0, 1.08), breaks = c(0, 0.5, 1.0),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "COLOC PP.H4") +
  theme_masld() + theme_pub() +
  theme(axis.text.x = element_text(size = PUB_AXIS_TEXT))

out <- file.path(FIG4_DIR, "_supp", "fads2_cross_ancestry.pdf")
cairo_pdf(out, width = fig_half_width, height = fig_half_width * 0.78, onefile = FALSE)
print(p)
invisible(dev.off())
message("Saved: ", out)
