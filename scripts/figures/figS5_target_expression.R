#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# Fig 5d COMPANION — hepatocyte expression on the calibration plane.
#
# Same axes as fig5c (x = genetic colocalization COLOC PP.H4; y = canonical raw
# logFC, disease vs control), but the point FILL encodes HEPATOCYTE EXPRESSION
# (single-cell mean CPM, log10 gradient) instead of drug-development stage.
# This answers "are the colocalizing / differentially-expressed targets actually
# expressed in hepatocytes?" — relevant because the in-vivo screen is hepatocyte-
# autonomous — directly on the evidence plane, rather than in a separate chart.
#
# Why a LOG colour scale: hepatocyte expression spans ~6 orders of magnitude
# across these genes (~0.002 CPM for LINC01561 to ~2,400 CPM for RORA), so a
# linear gradient would collapse almost everything; log10 spreads it cleanly.
# GLP1R / SLC5A2 (clinical MASLD targets, near-absent in hepatocytes) sit at the
# pale/low end; THRB / RORA / PNPLA3 / GPAM / MLIP at the dark/high end.
#
# House style: PDF only, all text BLACK, gene names italic, no title/subtitle/
# in-plot annotation sentences (descriptive line -> caption via message()).
# Env: rnaseq. Run from project root (lightweight render; login node OK).
#
# Output: figures/main/fig5_convergence/panels/figS5_target_expression.pdf
#         (reassigned from main Fig 5d to a supplementary panel, 2026-07-08;
#          kept in the same panels/ dir per the house layout)
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(viridisLite)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

# ── Data: the exact compact target set + values already plotted in fig5c ─────
#   coloc + logFC define the plane (identical to fig5c); hep_mean_cpm is the
#   hepatocyte single-cell mean CPM (our scRNA substrate); sig_class flags bulk-DEG
#   significance (kept for the audit CSV, not drawn — panel 5c carries it).
FIG5_DATA_DIR <- file.path(FIG5_DIR, "data")
dir.create(file.path(FIG5_DIR, "panels"), recursive = TRUE, showWarnings = FALSE)
dir.create(FIG5_DATA_DIR, recursive = TRUE, showWarnings = FALSE)
dt <- fread(file.path(FIG5_DATA_DIR, "drug_target_calibration.csv"),
            select = c("gene", "outcome", "coloc", "logFC", "treat_lfc", "treat_fdr",
                       "hep_mean_cpm", "hep_ratio", "sig_class", "hep_class",
                       "is_labeled", "is_plotted", "xj", "yj",
                       "label_nudge_x", "label_nudge_y"))
dt <- dt[is_plotted == TRUE]
dt[, plot_lfc := fifelse(is.na(logFC), 0, logFC)]   # GLP1R (no DEG) -> 0, as in 4h

FLOOR <- 0.01                                   # log-scale floor for genuine near-zeros
dt[, hep_plot := pmax(hep_mean_cpm, FLOOR)]

# ── Plot: COLOC × logFC plane, fill = log10 hepatocyte CPM ────────────────────
ymax <- max(abs(dt$plot_lfc), na.rm = TRUE)
ylim <- c(-ymax - 0.15, ymax + 0.15)
cpm_breaks <- c(0.01, 0.1, 1, 10, 100, 1000)

p <- ggplot(dt, aes(xj, yj)) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3) +
  geom_hline(yintercept = c(-0.25, 0.25), linetype = "22", color = "grey80", linewidth = 0.25) +
  geom_vline(xintercept = 0.5, linetype = "22", color = "grey75", linewidth = 0.3) +
  geom_point(aes(fill = hep_plot), shape = 21, color = "grey20",
             size = 2.2, stroke = 0.4) +
  geom_text_repel(data = dt[is_labeled == TRUE], aes(label = gene), color = "black",
                  fontface = "italic", size = GEOM_TEXT_6PT,
                  box.padding = 0.34, point.padding = 0.16, min.segment.length = 0,
                  force = 6, force_pull = 0.4,
                  segment.size = 0.18, segment.color = "grey65", seed = 1,
                  nudge_x = dt[is_labeled == TRUE]$label_nudge_x,
                  nudge_y = dt[is_labeled == TRUE]$label_nudge_y,
                  max.overlaps = Inf, max.time = 5, max.iter = 20000) +
  scale_fill_viridis_c(option = "mako", direction = -1, trans = "log10",
                       breaks = cpm_breaks,
                       labels = c("0", "0.1", "1", "10", "100", "1000"),
                       name = "Hepatocyte expression\n(scRNA mean CPM, log)",
                       guide = guide_colourbar(barwidth = 0.4, barheight = 3)) +
  scale_x_continuous(limits = c(-0.10, 1.20), breaks = c(0, 0.5, 1.0)) +
  scale_y_continuous(limits = ylim,
                     breaks = c(-1, -0.5, -0.25, 0, 0.25, 0.5, 1)) +
  labs(x = "coloc probabilities (PP.H4)",
       y = "Transcriptomic effect (logFC)") +
  theme_masld_compact() +
  theme(legend.position = "right",
        text = element_text(family = "Helvetica", size = 6, face = "plain"),
        axis.title = element_text(size = 6, face = "plain"),
        axis.text = element_text(size = 6, face = "plain"),
        legend.title = element_text(size = 6, face = "plain"),
        legend.text = element_text(size = 6, face = "plain"),
        plot.margin = margin(3, 4, 3, 3))

out <- file.path(FIG5_DIR, "panels", "figS5_target_expression.pdf")
save_fig(p, out, width = 3.4, height = 2.4)
message("Saved: ", out)
message("CAPTION: The fig5c calibration plane (x = COLOC PP.H4, y = canonical raw ",
        "logFC disease vs control) for the compact labelled target set, coloured by ",
        "HEPATOCYTE single-cell mean CPM (log10; human liver scRNA pseudobulk). ",
        "Expression spans ~6 orders of magnitude so the colour scale is log. ",
        "Near-zero hepatocyte expression floored at 0.01 CPM. Dashed y guides mark ",
        "the TREAT lfc=0.25 interval-null bound; significance is treat_fdr<0.05. ",
        "Drug-development stage and DEG significance are shown in panel 5c.")

# ── Audit CSV of plotted values ──────────────────────────────────────────────
fwrite(dt[order(-hep_mean_cpm),
          .(gene, outcome, coloc = round(coloc, 3), logFC = round(plot_lfc, 3),
            treat_lfc = round(treat_lfc, 3), treat_fdr = signif(treat_fdr, 3),
            hep_mean_cpm = round(hep_mean_cpm, 3), hep_ratio = round(hep_ratio, 3),
            hep_class, sig_class)],
       file.path(FIG5_DATA_DIR, "figS5_target_expression.csv"))

# ── Sanity print ─────────────────────────────────────────────────────────────
cat("==== fig5d expression-plane sanity (sorted by hepatocyte CPM) ====\n")
print(dt[order(-hep_mean_cpm),
         .(gene, coloc = round(coloc, 2), logFC = round(plot_lfc, 2),
           hep_cpm = round(hep_mean_cpm, 2), hep_class)])
cat(sprintf("\nPanel: %s\n", out))
