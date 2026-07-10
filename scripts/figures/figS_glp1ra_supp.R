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

## ── S1 fibrogenic reversal scatter ───────────────────────────────────────────
fb <- fread(file.path(GD, "glp1ra_fibrogenic_reversal.csv"))
fb <- fb[disease_meta_lfc > 0]                       # disease-up fibrogenic genes
fb[, wi := ifelse(cdahfd_reverses, "Weight-independent (CDA-HFD)", "Not / weight-dependent")]
p1 <- ggplot(fb, aes(disease_meta_lfc, treatment_meta_lfc)) +
  geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey60") +
  geom_point(aes(colour = wi), size = 1.1) +
  geom_text_repel(aes(label = human_symbol), fontface = "italic", size = 1.6,
                  segment.size = 0.15, max.overlaps = 40, min.segment.length = 0) +
  scale_colour_manual(values = c("Weight-independent (CDA-HFD)" = "#08519c",
                                 "Not / weight-dependent" = "#9E9E9E"), name = NULL) +
  labs(x = "MASLD disease log2FC (mouse liver, meta)",
       y = "Semaglutide treatment log2FC") +
  theme_masld() + theme(legend.position = "top")
ggsave(file.path(OUT, "figS_glp1ra_fibrogenic_reversal.pdf"), p1,
       width = 3.4, height = 3.2, device = pdf_device)

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
message("S1 CAPTION: Core fibrogenic program (35 genes) is disease-up and reversed by ",
        "semaglutide (32/35 at meta); 27/35 reverse in the weight-stable CDA-HFD model.")
message("S2 CAPTION: Continuous reversal score across disease-significant convergence ",
        "Tier-1 targets; median fraction corrected = ", round(med, 2),
        "; most targets shift toward normal, a small residue (<0) moves with disease.")
