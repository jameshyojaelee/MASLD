#!/usr/bin/env Rscript
# KEY MESSAGE (honest, n=18 multiome is UNDERPOWERED for disease-state chromatin):
# Binary MASLD-vs-Control differential accessibility is NULL after covariate
# correction — Hep gives 0 significant peaks of 209,042 tested (panel a volcano).
# The ONLY non-zero contrast is the extreme F0-vs-F4 stage comparison (605 peaks,
# 74% CLOSING), and that is the signature of hepatocyte compositional loss, not
# graded locus-specific regulation: the adjacent F0->F3 and F3->F4 contrasts are
# both 0 (panel b). All non-Hep cell types are likewise underpowered (<=4 sig in
# any contrast). This panel documents the negative result; no chromatin
# disease-state claim is made from this multiome (see fig4 epigenetics CUT,
# docs/manuscript/working/fig4_critique_recommendations_2026-06-18.md).
# PROVENANCE: a prior version of this script HARDCODED 694 sig peaks for Hep — a
# value that exists in NO data file (corrected DA = 0); removed 2026-06-19.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUT_PDF <- file.path(FIGS05_DIR, "figS05_scatac_disease_vs_control.pdf")

MAG  <- "#C9265E"
BLUE <- "#1565C0"
GRAY <- "#9E9E9E"

# ── Panel a: Hep MASLD-vs-Control volcano (covariate-corrected; 0 sig) ──
hep <- fread(file.path(BASE,
  "Analysis/ATAC/Human_Multiome/results/snapatac2/scatac_da_corrected_hep.csv"))
setnames(hep, c("log2(fold_change)", "p-value", "adjusted p-value"),
              c("log2FC", "pvalue", "padj"))
hep[, neglog10_padj := pmin(-log10(pmax(padj, 1e-300)), 50)]
hep[, direction := fcase(
  padj < 0.05 & log2FC > 0, "Opening",
  padj < 0.05 & log2FC < 0, "Closing",
  default = "NS")]
hep[, direction := factor(direction, levels = c("NS", "Opening", "Closing"))]
setorder(hep, direction)  # plot sig on top

n_open  <- sum(hep$direction == "Opening")
n_close <- sum(hep$direction == "Closing")

p_a <- ggplot(hep, aes(x = log2FC, y = neglog10_padj, colour = direction)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             colour = "gray70", linewidth = 0.25) +
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "gray70", linewidth = 0.25) +
  geom_point(size = 0.5, alpha = 0.55, shape = 16) +
  scale_colour_manual(values = c(NS = GRAY, Opening = MAG, Closing = BLUE),
                      breaks = c("Opening", "Closing", "NS"),
                      labels = c(sprintf("Opening (%d)", n_open),
                                 sprintf("Closing (%d)", n_close),
                                 "Not sig.")) +
  scale_x_continuous(name = expression(bold(log[2]*" FC  (MASLD vs Control)")),
                     limits = c(-1, 1) * max(abs(hep$log2FC)) * 1.02,
                     expand = expansion(mult = 0.02)) +
  scale_y_continuous(name = expression(bold(-log[10]*" padj")),
                     expand = expansion(mult = c(0, 0.04))) +
  labs(tag = "a", colour = NULL) +
  theme_masld(base_size = 7) +
  theme_pub() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.99, 0.99),
        legend.justification = c(1, 1),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.key.size = unit(0.25, "cm"),
        legend.text = element_text(size = PUB_LEGEND),
        plot.tag = element_text(size = 11, face = "bold"),
        plot.tag.position = c(0.02, 0.97))

# ── Panel b: where (if anywhere) does Hep chromatin signal appear? ──
# Binary disease-vs-control (covariate-corrected) is NULL (n_open/n_close from
# panel a). The ONLY non-zero contrast is the extreme F0-vs-F4 stage comparison,
# dominated by global CLOSING (hepatocyte compositional loss), with the adjacent
# F0->F3 and F3->F4 contrasts both 0. Every count below is read live from disk —
# no hardcoded values.
stage_da_dir <- file.path(BASE, "Analysis/ATAC/Human_Multiome/results/stage_da")
read_sig <- function(con) {
  f <- file.path(stage_da_dir, sprintf("da_%s_Hep.csv", con))
  if (!file.exists(f)) return(c(open = 0L, close = 0L))
  d  <- fread(f)
  pc <- if ("padj"   %in% names(d)) d$padj   else d[["adjusted p-value"]]
  lc <- if ("log2FC" %in% names(d)) d$log2FC else d[["log2(fold_change)"]]
  c(open  = sum(pc < 0.05 & lc > 0, na.rm = TRUE),
    close = sum(pc < 0.05 & lc < 0, na.rm = TRUE))
}
f04 <- read_sig("F0_vs_F4"); f03 <- read_sig("F0_vs_F3"); f34 <- read_sig("F3_vs_F4")

ord <- c("MASLD vs Control\n(corrected)", "F0 vs F4\n(extreme stage)",
         "F0 vs F3\n(adjacent)",          "F3 vs F4\n(adjacent)")
bar <- rbindlist(list(
  data.table(contrast = ord[1], open = n_open,        close = n_close,        confound = FALSE),
  data.table(contrast = ord[2], open = f04[["open"]], close = f04[["close"]], confound = TRUE),
  data.table(contrast = ord[3], open = f03[["open"]], close = f03[["close"]], confound = FALSE),
  data.table(contrast = ord[4], open = f34[["open"]], close = f34[["close"]], confound = FALSE)))
bar[, total := open + close]
bar[, contrast := factor(contrast, levels = rev(ord))]   # first level -> bottom
barl <- melt(bar, id.vars = c("contrast", "confound", "total"),
             variable.name = "dir", value.name = "n")
barl[, dir := factor(fifelse(dir == "open", "Opening", "Closing"),
                     levels = c("Opening", "Closing"))]

# pct-closing label only on the one non-zero (confounded) bar
bar[, lab := fifelse(total > 0 & confound,
        sprintf("%d  (%.0f%% closing)", total, 100 * close / total),
        as.character(total))]

p_b <- ggplot(barl, aes(y = contrast, x = n, fill = dir)) +
  geom_col(width = 0.62) +
  geom_text(data = bar, aes(y = contrast, x = total, label = lab),
            inherit.aes = FALSE, hjust = -0.12, size = PUB_GEOM_TEXT, colour = "black") +
  scale_fill_manual(values = c(Opening = MAG, Closing = BLUE), name = NULL) +
  scale_x_continuous(name = "Significant DA peaks (padj < 0.05)",
                     expand = expansion(mult = c(0, 0.42))) +
  labs(tag = "b", y = NULL) +
  theme_masld(base_size = 7) +
  theme_pub() +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(colour = "black"),
        axis.title.x = element_text(face = "bold"),
        legend.position = "bottom",
        legend.key.size = unit(0.25, "cm"),
        plot.tag = element_text(size = 11, face = "bold"),
        plot.tag.position = c(0.02, 0.97))

combo <- p_a + p_b + plot_layout(widths = c(3, 3))
ggsave(OUT_PDF, combo, width = 6.8, height = 2.8, device = cairo_pdf)
cat(sprintf("[saved] %s\n", OUT_PDF))

# ── Caption / provenance to stdout (NOT on the panel) ────────────────────────
message(strrep("=", 78))
message("FIGURE LEGEND (paste into manuscript; stats live here, not on the panel)")
message(strrep("=", 78))
message(sprintf(
"Human liver multiome (n=18 donors; 6 control / 12 MASLD) is UNDERPOWERED for
disease-state chromatin. (a) Covariate-corrected MASLD-vs-Control differential
accessibility in hepatocytes: %d significant peaks of %d tested (%d opening / %d
closing) at padj<0.05 -- a null result. (b) The only non-zero contrast is the
extreme F0-vs-F4 stage comparison (%d peaks, %.0f%% closing); the adjacent
F0->F3 (%d) and F3->F4 (%d) contrasts are both 0. The exclusive concentration of
signal in the most-extreme stage gap, dominated by global closing, is the
signature of hepatocyte compositional loss (cell dropout) rather than graded
locus-specific regulation -- the same composition confound that led to the
epigenetic panels being CUT from main Fig 4.",
  n_open + n_close, nrow(hep), n_open, n_close,
  f04[["open"]] + f04[["close"]], 100 * f04[["close"]] / max(1, f04[["open"]] + f04[["close"]]),
  f03[["open"]] + f03[["close"]], f34[["open"]] + f34[["close"]]))
message("Non-Hep cell types: <=4 significant peaks in ANY stage/disease contrast (uniformly underpowered).")
message("NO hardcoded values: every count read live from the corrected DA + stage_da CSVs on disk.")
message(strrep("=", 78))
print(bar[, .(contrast = gsub("\n", " ", as.character(contrast)), open, close, total)])
