# KEY MESSAGE: Most colocalized (genetically causal) genes carry NO disease transcriptional change, but a few drug-relevant genes sit at both extremes — THRB and RORA are both highly colocalized AND transcriptionally suppressed, while PPARG/PPARA are dysregulated yet have no GWAS-eQTL support. A Fig3->Fig2 bridge panel.
##############################################################################
# coloc_deg_bridge.R
# Fig 3E — DEG-centric bridge between the RNA-seq chapter (Fig 3) and the
# genetics chapter (Fig 2). One point per gene: colocalization posterior
# (y = PP.H4) vs transcriptomic effect (x = bulk log2FC). Non-DEGs faded;
# DEGs emphasized; only the named drug-relevant genes are labeled.
#
# Output (PDF only): FIG2_DIR/panels/figs3g_coloc_deg_bridge.pdf
#   (FIG2_DIR resolves to figures/main/fig3_RNAseq — known back-compat misnomer)
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
  select = c("human_symbol", "bulk_logFC", "bulk_padj",
             "coloc_susie_best_pp4", "coloc_abf_best_pp4"))
setnames(atlas, "human_symbol", "gene")

# Colocalization posterior: prefer SuSiE (matches manuscript THRB = 0.9999),
# fall back to ABF; untested genes -> 0.
atlas[, pp4 := fcoalesce(coloc_susie_best_pp4, coloc_abf_best_pp4)]
atlas[is.na(pp4), pp4 := 0]

# Plottable universe: genes with a measured transcriptomic effect.
atlas <- atlas[!is.na(bulk_logFC)]

# Primary-DEG flag (Tier-1 thresholds).
atlas[, is_deg := !is.na(bulk_padj) & bulk_padj < 0.05 & abs(bulk_logFC) > 0.5]

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
key_high      <- c("THRB", "RORA")          # high PP.H4 + transcriptionally SUPPRESSED
key_conv      <- c("HKDC1")                 # high PP.H4 + INDUCED — convergent anchor (hero of Fig 3I)
key_zero_up   <- c("FASN", "SCD", "FGF21")  # INDUCED drug targets, no genetic support
key_zero_down <- c("CYP2C19", "GSTM1", "IGFBP2")  # SUPPRESSED hepatic genes, no genetic support
key_zero      <- c(key_zero_up, key_zero_down)
# Additional genes that are BOTH highly colocalized (PP.H4 > 0.5) AND Tier-1
# DEGs — i.e. genetically causal AND transcriptionally dysregulated, beyond the
# named drug/convergent anchors. Chosen for biological recognizability and
# spatial spread (so labels don't pile up): SORT1/F2RL1/ANGPTL7 induced,
# CETP/CUX2 suppressed. Same magenta "causal" coloring as the anchors, smaller
# points so THRB/RORA/HKDC1 stay the dominant heroes.
key_causal_deg <- c("SORT1", "F2RL1", "ANGPTL7",   # induced (right side)
                    "CETP", "CUX2")                # suppressed (left side)
key_all       <- c(key_high, key_conv, key_causal_deg, key_zero)

atlas[, klass := fcase(
  gene %in% key_high,           "causal_supported",
  gene %in% key_conv,           "convergent",
  gene %in% key_causal_deg,     "causal_deg",
  gene %in% key_zero,           "drug_unsupported",
  is_deg == TRUE,               "deg",
  default                        = "nondeg"
)]

lab_dt <- atlas[gene %in% key_all]

# Manuscript-stated enrichment (text, not recomputed here): primary DEGs show a
# modest 1.24x enrichment for colocalization over non-DEGs (Wilcoxon p = 1.4e-38).

# ── Colors ──────────────────────────────────────────────────────────────────
COL_HIGH   <- masld_colors$mash      # "#C9265E" deep magenta — causal + suppressed
COL_ZERO   <- "#1565C0"              # deep blue — dysregulated, no genetic support
COL_DEG    <- "#C9265E"              # DEGs emphasized in magenta
COL_NONDEG <- masld_colors$control   # "#9E9E9E" neutral gray — faded background

# ── Plot ────────────────────────────────────────────────────────────────────
p <- ggplot() +
  # PP.H4 > 0.5 colocalization band
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.5, ymax = 1.04,
           fill = "#FBE9EF", alpha = 0.4) +
  geom_hline(yintercept = 0.5, linetype = "dashed",
             linewidth = 0.3, color = "gray55") +
  geom_vline(xintercept = 0, linewidth = 0.25, color = "gray70") +
  # Non-DEG cloud — faded, rasterized
  rasterize(
    geom_point(data = atlas[klass == "nondeg"],
               aes(x = bulk_logFC, y = pp4),
               color = COL_NONDEG, size = 0.3, alpha = 0.22, shape = 16),
    dpi = 600) +
  # DEGs — emphasized
  rasterize(
    geom_point(data = atlas[klass == "deg"],
               aes(x = bulk_logFC, y = pp4),
               color = COL_DEG, size = 0.65, alpha = 0.5, shape = 16),
    dpi = 600) +
  # Additional causal DEGs (high PP.H4 + Tier-1 DEG) — magenta, smaller than the
  # named heroes so THRB/RORA/HKDC1 still dominate; white halo to read above cloud
  geom_point(data = atlas[klass == "causal_deg"],
             aes(x = bulk_logFC, y = pp4),
             color = COL_HIGH, size = 1.9, shape = 16) +
  geom_point(data = atlas[klass == "causal_deg"],
             aes(x = bulk_logFC, y = pp4),
             color = "white", size = 1.9, shape = 1, stroke = 0.35) +
  # Named genes — solid fill + white halo so they read above the cloud
  geom_point(data = atlas[klass == "drug_unsupported"],
             aes(x = bulk_logFC, y = pp4),
             color = COL_ZERO, size = 2.6, shape = 16) +
  geom_point(data = atlas[klass == "drug_unsupported"],
             aes(x = bulk_logFC, y = pp4),
             color = "white", size = 2.6, shape = 1, stroke = 0.4) +
  geom_point(data = atlas[klass %in% c("causal_supported", "convergent")],
             aes(x = bulk_logFC, y = pp4),
             color = COL_HIGH, size = 3.0, shape = 16) +
  geom_point(data = atlas[klass %in% c("causal_supported", "convergent")],
             aes(x = bulk_logFC, y = pp4),
             color = "white", size = 3.0, shape = 1, stroke = 0.4) +
  geom_text_repel(
    data = lab_dt,
    aes(x = bulk_logFC, y = pp4, label = gene,
        color = klass),
    size = 2.55, fontface = "italic",
    box.padding = 0.7, point.padding = 0.4,
    segment.size = 0.25, segment.color = "gray55",
    min.segment.length = 0, max.overlaps = Inf,
    seed = 42, force = 9, force_pull = 0.4,
    max.iter = 100000, max.time = 3,
    bg.color = "white", bg.r = 0.12,
    show.legend = FALSE) +
  # PP.H4 = 0.5 reference label
  annotate("text", x = Inf, y = 0.515, label = "PP.H4 = 0.5",
           hjust = 1.05, size = 1.9, color = "gray50", fontface = "italic") +
  # Single terse enrichment stat, top-right (all other detail lives in caption).
  annotate("text", x = Inf, y = 1.045,
           label = "1.24× coloc enrichment (p = 1.4×10⁻³⁸)",
           hjust = 1.04, vjust = 1, size = 1.9, color = "gray45") +
  scale_color_manual(values = c(causal_supported = COL_HIGH,
                                convergent        = COL_HIGH,
                                causal_deg        = COL_HIGH,
                                drug_unsupported  = COL_ZERO),
                     guide = "none") +
  scale_x_continuous(expand = expansion(mult = 0.06)) +
  scale_y_continuous(limits = c(-0.01, 1.06),
                     breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
  labs(title = "Causal colocalization vs. disease transcription",
       x = expression(log[2]~"fold change (MASLD vs. control)"),
       y = "Colocalization posterior (PP.H4)") +
  theme_masld() +
  theme(plot.title = element_text(size = 7.5, face = "plain",
                                  margin = margin(b = 4)))

# ── Save ────────────────────────────────────────────────────────────────────
# save_fig() writes PDF via cairo_pdf (vector, no dingbats) — PDF only.
out_pdf <- file.path(PANEL_DIR, "figs3g_coloc_deg_bridge.pdf")
save_fig(p, out_pdf, width = fig_half_width + 0.7, height = 4.1)
cat("[fig3E] Saved:", out_pdf, "\n")

# ── Sidecar data CSV ─────────────────────────────────────────────────────────
csv_dt <- atlas[, .(gene,
                    bulk_logFC = round(bulk_logFC, 4),
                    bulk_padj  = signif(bulk_padj, 4),
                    PP_H4      = round(pp4, 4),
                    is_deg,
                    class      = klass)][order(-PP_H4)]
out_csv <- file.path(DATA_DIR, "coloc_deg_bridge.csv")
fwrite(csv_dt, out_csv)
cat("[fig3E] Sidecar:", out_csv, "\n")

# ── Sanity check (printed to log) ────────────────────────────────────────────
cat("\n[fig3E] Named-gene PP.H4 check:\n")
print(lab_dt[order(-pp4), .(gene, bulk_logFC = round(bulk_logFC, 3),
                            PP_H4 = round(pp4, 4), klass)])
cat(sprintf("[fig3E] n total plotted = %d ; n DEGs = %d ; n PP.H4>0.5 = %d\n",
            nrow(atlas), sum(atlas$is_deg), sum(atlas$pp4 > 0.5)))
