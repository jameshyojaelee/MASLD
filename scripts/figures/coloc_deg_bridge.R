# KEY MESSAGE: Most colocalized (genetically causal) genes carry NO disease transcriptional change, but a few drug-relevant genes sit at both extremes — THRB and RORA are both highly colocalized AND transcriptionally suppressed, while PPARG/PPARA are dysregulated yet have no GWAS-eQTL support. A Fig3->Fig2 bridge panel.
##############################################################################
# coloc_deg_bridge.R
# Fig S3 (coloc-DEG bridge; demoted from main 3E) — DEG-centric bridge between the RNA-seq chapter (Fig 3) and the
# genetics chapter (Fig 2). One point per gene: colocalization posterior
# (y = PP.H4) vs transcriptomic effect (x = bulk log2FC). Non-DEGs faded;
# DEGs emphasized; only the named drug-relevant genes are labeled.
#
# Output (PDF only): FIG2_DIR/panels/figs3_coloc_deg_bridge.pdf   (SUPP; was main 3E, demoted to supp
#   — the disconnect panel — gene-level LFC vs PP.H4
#   is the statistically honest disconnect; the 7-lineage cell-type scatter was
#   retired, its correlation being non-significant (p=0.35, hepatocyte-outlier-driven).
#   FIG2_DIR resolves to figures/main/fig3_RNAseq — known back-compat misnomer.)
# Sidecar: FIG2_DIR/panels/data/coloc_deg_bridge.csv
##############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(scales)
  library(ggrastr)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

PANEL_DIR <- file.path(FIG2_DIR, "panels")
DATA_DIR  <- file.path(PANEL_DIR, "data")
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Data ────────────────────────────────────────────────────────────────────
atlas <- fread(
  file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"),
  select = c("human_symbol", "bulk_logFC", "bulk_padj", "bulk_treat_fdr"))
setnames(atlas, "human_symbol", "gene")

# Colocalization posterior from the CANONICAL gene-level source (NOT the atlas
# coloc_susie_best_pp4 column, which the 2026-06-14 mega-review flagged as
# corrupted — a few PP.H4 values > 1 from ABF-fallback mis-sorting). This is the
# same source Fig 2 uses; SuSiE-preferred with ABF fallback, all values <= 1.
glc <- fread(file.path(BASE, "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv"),
             select = c("gene", "coloc_best_pp4", "coloc_best_susie_ancestry",
                        "coloc_susie_headline_cross_anc"))
glc <- glc[gene != "" & !is.na(gene)][
  , .(pp4            = max(coloc_best_pp4, na.rm = TRUE),
      # cross-ancestry provenance: headline PP.H4 driven by a NON-EUR GWAS. The eQTL
      # panel is EUR, so a non-EUR-GWAS headline holds a lower evidentiary bar. Named
      # anchor genes carrying this flag are marked "*" in the plot (caption below).
      cross_anc      = any(coloc_susie_headline_cross_anc %in% TRUE),
      susie_ancestry = coloc_best_susie_ancestry[1]), by = gene]
atlas <- merge(atlas, glc, by = "gene", all.x = TRUE)
atlas[is.na(pp4), pp4 := 0]
atlas[is.na(cross_anc), cross_anc := FALSE]

# Plottable universe: genes with a measured transcriptomic effect.
atlas <- atlas[!is.na(bulk_logFC)]

# Primary-DEG flag (Tier-1, TREAT canonical: treat_fdr<0.05 at lfc=0.25; the
# effect-size floor is folded into the test, so there is no separate |logFC| filter).
atlas[, is_deg := is_canonical_deg(atlas)]

# Named genes (manuscript). THRB/RORA = high PP.H4 + concordant suppression
# (THRB = resmetirom, the approved MASH target). FASN/SCD/FGF21 = leading MASLD
# drug targets (denifanstat / aramchol / efruxifermin) that are Tier-1 DEGs yet
# sit near PP.H4 = 0 (no common-variant colocalization support).
# (The earlier draft named GLP1R/PPARG/PPARA: GLP1R carries no bulk_logFC — not
#  expressed/tested in the hepatic bulk atlas; PPARG is not dysregulated
#  (logFC 0.02, padj 0.67); PPARA is sub-threshold. All replaced.)
# CYP2C19/GSTM1/IGFBP2 = strongly SUPPRESSED hepatic genes (drug metabolism /
#  glutathione detox / IGF axis) that, like the induced drug targets, are Tier-1
#  DEGs with no colocalization support — they anchor the down-regulated, low-PP.H4
#  quadrant so the "dysregulated yet not genetically causal" message is labeled on
#  both sides of zero, not only the up-regulated side.
# Labels restricted 2026-07-02 to genes discussed in the other main figures (Fig 2
# genetics / Fig 3 RNA-seq / Fig 4 validation), spanning the disconnect plane.
key_high      <- c("THRB", "RORA")          # colocalized + suppressed (Fig 2/3/4)
key_conv      <- c("HKDC1")                 # colocalized + induced, convergent hero (Fig 3/4)
key_zero_up   <- c("FASN", "SCD", "FGF21", "AKR1B10", "SERPINE1")  # strongly dysregulated, PP.H4~0 (Fig 3/4)
key_zero_down <- character(0)
key_zero      <- c(key_zero_up, key_zero_down)
key_causal_deg <- c("CYP3A4")   # colocalized + suppressed (Fig 4 zonation)
key_all       <- c(key_high, key_conv, key_causal_deg, key_zero)

atlas[, klass := fcase(
  gene %in% key_high,           "causal_supported",
  gene %in% key_conv,           "convergent",
  gene %in% key_causal_deg,     "causal_deg",
  gene %in% key_zero,           "drug_unsupported",
  is_deg == TRUE,               "deg",
  default                        = "nondeg"
)]

# Highlight class: ONLY the labeled anchor genes get colour (magenta = colocalized
# PP.H4>0.5 anchors; blue = dysregulated drug targets with no coloc). Everything
# else is uniform gray background.
atlas[, hl := fcase(
  klass %in% c("causal_supported", "convergent", "causal_deg"), "coloc",
  klass == "drug_unsupported", "drug",
  default = NA_character_)]
lab_dt <- atlas[gene %in% key_all]
# Mark anchors whose SuSiE COLOC headline is non-EUR (cross-ancestry) with a trailing
# "*" — the eQTL panel is EUR, so those hold a lower evidentiary bar.
lab_dt[, lab := fifelse(cross_anc == TRUE, paste0(gene, "*"), gene)]

# Coloc enrichment — computed LIVE from the plotted data (never hardcoded, so it
# tracks the canonical DEG definition). "Enrichment" = mean PP.H4 in Tier-1 DEGs
# relative to non-DEGs; significance = one-sided Wilcoxon (DEGs > non-DEGs). Under
# the TREAT canonical the effect is modest (~1.15x) yet highly significant, which
# is exactly the figure's message: the coloc and transcription axes are largely
# decoupled, with only a few genes (THRB/RORA/HKDC1) at both extremes.
.sup_map <- c("0"="⁰","1"="¹","2"="²","3"="³","4"="⁴",
              "5"="⁵","6"="⁶","7"="⁷","8"="⁸","9"="⁹","-"="⁻")
.fmt_p <- function(p) {
  if (!is.finite(p) || p <= 0) return("< 1×10⁻³⁰⁰")
  e <- floor(log10(p)); m <- round(p / 10^e)
  if (m >= 10) { m <- 1L; e <- e + 1L }
  esup <- paste(.sup_map[strsplit(as.character(e), "")[[1]]], collapse = "")
  sprintf("%d×10%s", m, esup)
}
enr_deg <- atlas[is_deg == TRUE,  pp4]
enr_non <- atlas[is_deg == FALSE, pp4]
enr_ratio <- mean(enr_deg) / mean(enr_non)
enr_w <- suppressWarnings(wilcox.test(enr_deg, enr_non, alternative = "greater"))
enr_lab <- sprintf("%.2f× coloc enrichment (Wilcoxon p = %s)",
                   enr_ratio, .fmt_p(enr_w$p.value))
cat(sprintf("[coloc] coloc enrichment: ratio=%.3fx  Wilcoxon(greater) p=%.3g\n",
            enr_ratio, enr_w$p.value))

# ── Colors (only the labeled anchor genes are coloured) ──────────────────────
COL_COLOC <- masld_colors$mash       # "#C9265E" magenta — colocalized anchors
COL_DRUG  <- "#1565C0"               # blue — dysregulated drug targets, no coloc
COL_BG    <- masld_colors$control    # "#9E9E9E" gray — uniform background

# ── Plot (stripped: gray cloud + coloured labeled anchors; uniform dot sizes; no
#    on-panel stats, no PP.H4 tick label) ──────────────────────────────────────
p <- ggplot() +
  geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.3, color = "gray55") +
  geom_vline(xintercept = 0, linewidth = 0.25, color = "gray80") +
  # all genes — ONE uniform gray background layer (rasterized)
  rasterize(
    geom_point(data = atlas, aes(x = bulk_logFC, y = pp4),
               color = COL_BG, size = 0.35, alpha = 0.28, shape = 16), dpi = 600) +
  # labeled anchor genes — uniform size, coloured by role, white halo
  geom_point(data = atlas[!is.na(hl)], aes(x = bulk_logFC, y = pp4, color = hl),
             size = 1.7, shape = 16) +
  geom_point(data = atlas[!is.na(hl)], aes(x = bulk_logFC, y = pp4),
             color = "white", size = 1.7, shape = 1, stroke = 0.3) +
  geom_text_repel(
    data = lab_dt, aes(x = bulk_logFC, y = pp4, label = lab),
    color = "black", size = 6 / ggplot2::.pt, fontface = "italic",
    box.padding = 0.45, point.padding = 0.25,
    segment.size = 0.2, segment.color = "gray60",
    min.segment.length = 0, max.overlaps = Inf,
    seed = 42, force = 7, bg.color = "white", bg.r = 0.12, show.legend = FALSE) +
  scale_color_manual(values = c(coloc = COL_COLOC, drug = COL_DRUG), guide = "none") +
  scale_x_continuous(expand = expansion(mult = 0.06)) +
  scale_y_continuous(limits = c(-0.01, 1.02), breaks = c(0, 0.5, 1.0)) +
  labs(x = expression(log[2]~"fold change (MASLD vs. control)"),
       y = "Colocalization posterior (PP.H4)") +
  theme_masld_compact()

# ── Save ────────────────────────────────────────────────────────────────────
# save_fig() writes PDF via cairo_pdf (vector, no dingbats) — PDF only.
out_pdf <- file.path(PANEL_DIR, "figs3_coloc_deg_bridge.pdf")
save_fig(p, out_pdf, width = fig_half_width, height = 3.0)
cat("[coloc] Saved:", out_pdf, "\n")

# Caption note for cross-ancestry provenance flagging.
.xanc <- lab_dt[cross_anc == TRUE, gene]
message("[coloc caption] genes marked * have a non-EUR (cross-ancestry) COLOC headline ",
        "(EUR eQTL — interpret with caution). Flagged anchors: ",
        if (length(.xanc)) paste(.xanc, collapse = ", ") else "none")

# ── Sidecar data CSV ─────────────────────────────────────────────────────────
csv_dt <- atlas[, .(gene,
                    bulk_logFC = round(bulk_logFC, 4),
                    bulk_padj  = signif(bulk_padj, 4),
                    PP_H4      = round(pp4, 4),
                    cross_anc,
                    susie_ancestry,
                    is_deg,
                    class      = klass)][order(-PP_H4)]
out_csv <- file.path(DATA_DIR, "coloc_deg_bridge.csv")
fwrite(csv_dt, out_csv)
cat("[coloc] Sidecar:", out_csv, "\n")

# ── Sanity check (printed to log) ────────────────────────────────────────────
cat("\n[coloc] Named-gene PP.H4 check:\n")
print(lab_dt[order(-pp4), .(gene, bulk_logFC = round(bulk_logFC, 3),
                            PP_H4 = round(pp4, 4), klass)])
cat(sprintf("[coloc] n total plotted = %d ; n DEGs = %d ; n PP.H4>0.5 = %d\n",
            nrow(atlas), sum(atlas$is_deg), sum(atlas$pp4 > 0.5)))
