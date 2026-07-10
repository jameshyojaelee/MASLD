#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# Supplementary panels for the GLP-1RA section (mechanism-led; see Fig 5e for the
# main receptor-axis panel). Each panel is its own PDF (individual-panels rule).
#
#  S1 fibrogenic reversal : mouse-meta DISEASE logFC (x) vs semaglutide TREATMENT
#     logFC (y) for the 35-gene core fibrogenic program. Disease-up + treatment-down
#     = lower-right = reversal. Colour = reverses weight-independently in the
#     weight-stable CDA-HFD model (single model; per RevC scope).
#  S2 reversal breadth   : distribution of the continuous reversal score
#     (fraction_corrected) across disease-significant convergence Tier-1 targets —
#     most are corrected toward normal; a small residue moves with disease.
#
# NO lollipop; text BLACK; gene symbols italic; no titles (captions via message()).
# Output: figures/supplementary/figS_glp1ra/*.pdf   Env: rnaseq
# ─────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(ggrepel)
})
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
pdf_device <- if (capabilities("cairo")) cairo_pdf else grDevices::pdf
OUT <- file.path(FIG_SUPP, "figS_glp1ra"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
GD <- file.path(BASE, "RNA-seq/results/glp1ra")

## ── S1 fibrogenic reversal scatter (3-way category, all 35 genes incl. FN1) ───
# Reversal is defined relative to EACH gene's disease direction, so both disease-up
# genes (driven DOWN by drug, lower-right) and the single disease-down gene FN1
# (driven UP, upper-left) are shown. Colour = 3 honest classes, not the old binary
# that conflated "weight-dependent reverser" with "non-reverser".
fb <- fread(file.path(GD, "glp1ra_fibrogenic_reversal.csv"))          # all 35, no filter
# Colour = significantly reversed in the SAME pooled meta the axes show (axes+colour consistent).
# The CDA-HFD/non-obesity result (27/35) is stated in the caption/text, not encoded as colour.
fb[, revcat := fifelse(meta_reverses, "Reversed by semaglutide", "Not reversed")]
fb[, revcat := factor(revcat, levels = c("Reversed by semaglutide", "Not reversed"))]
p1 <- ggplot(fb, aes(disease_meta_lfc, treatment_meta_lfc)) +
  geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey70") +
  geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey70") +
  geom_point(aes(colour = revcat), size = 1.1) +
  geom_text_repel(aes(label = human_symbol), fontface = "italic", size = 1.6,
                  segment.size = 0.15, max.overlaps = 45, min.segment.length = 0) +
  scale_colour_manual(values = c("Reversed by semaglutide" = "#08519c",
                                 "Not reversed"            = "#9E9E9E"),
                      name = NULL) +
  labs(x = "Disease log2FC (MASLD vs control, mouse liver 5-model meta)",
       y = "Semaglutide log2FC (treated vs vehicle, meta)") +
  theme_masld() + theme(legend.position = "top")
ggsave(file.path(OUT, "figS_glp1ra_fibrogenic_reversal.pdf"), p1,
       width = 3.6, height = 3.5, device = pdf_device)

## ── S2 reversal breadth (continuous score) ───────────────────────────────────
part <- fread(file.path(GD, "glp1ra_target_reversal_partition.csv"))
br <- part[conv_tier1 == TRUE & mouse_disease_sig == TRUE & is.finite(fraction_corrected)]
med <- median(br$fraction_corrected, na.rm = TRUE)
p2 <- ggplot(br, aes(fraction_corrected)) +
  geom_histogram(binwidth = 0.1, fill = "#6baed6", colour = "white", linewidth = 0.1) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey40") +
  geom_vline(xintercept = med, linewidth = 0.3, linetype = "dashed", colour = "#08306b") +
  labs(x = "Fraction of disease dysregulation corrected by semaglutide",
       y = "Convergence Tier-1 targets") +
  theme_masld()
ggsave(file.path(OUT, "figS_glp1ra_reversal_breadth.pdf"), p2,
       width = 3.2, height = 2.2, device = pdf_device)

message("Supp panels written to ", OUT)
message("S1 CAPTION: Core fibrogenic program (35 genes incl. the disease-DOWN gene FN1). ",
        "Both axes = pooled 5-model mouse-liver meta; colour = significantly reversed in that ",
        "SAME meta (direction-aware: disease-up driven down = lower-right; FN1 driven up = upper-left). ",
        "32/35 reverse (blue); 3 not reversed (grey: MMP9/CCN2/MMP14). Of the 32, 27 also reverse ",
        "in the non-obesity CDA-HFD model (antifibrotic effect not restricted to obesity-driven disease). ",
        "Weight-independence of the HUMAN fibrosis benefit = trial mediation (25% weight-mediated; ",
        "Jara Nat Med 2025), not this figure.")
message("S2 CAPTION: Continuous reversal score across disease-significant convergence ",
        "Tier-1 targets; median fraction corrected = ", round(med, 2),
        "; most targets shift toward normal, a small residue (<0) moves with disease.")
