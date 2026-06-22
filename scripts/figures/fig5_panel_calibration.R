##############################################################################
# fig5_panel_calibration.R — Fig 5 PANEL 5d (NEW, 2026-06-10)
#
# Held-out calibration of the multi-evidence CONVERGENCE SCORE (Script 46d,
# the canonical evidence-weighted heuristic convergence ranking).
#
# What it shows: the convergence score ranks genes in three EXTERNAL held-out
# truth panels (Govaere consensus MASLD signature, NIDDK, Open Targets MASLD
# associations) better than chance, and improves over our OWN prior ranking
# method (Script 46b, retired). The held-out panels are TRUTH, not competitors.
#
# HONESTY GUARDRAIL (do NOT violate): 46b is OUR retired method; Govaere /
# NIDDK / Open Targets are held-out evaluation panels (TRUTH sets), NOT rival
# tools. Framing is strictly "convergence score predicts held-out panels and
# improves over our prior method." NEVER state or imply that the convergence
# score outperforms any of those external resources as a competing tool.
#
# Data: RNA-seq/results/multi_evidence/convergence_evidence_benchmark.csv
#   columns auroc_46d (convergence score) + auroc_46b (prior method) per panel.
#
# Output: $FIG5_DIR/panels/fig5d_calibration.pdf  (+ source CSV)
##############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUTDIR <- FIG5_DIR
PANDIR <- file.path(OUTDIR, "panels")
dir.create(PANDIR, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Load held-out benchmark
# ---------------------------------------------------------------------------
bench_f <- file.path(BASE, "RNA-seq/results/multi_evidence",
                     "convergence_evidence_benchmark.csv")
stopifnot(file.exists(bench_f))
bench <- fread(bench_f)

# Human-readable panel labels + ordering (best discrimination first).
# Govaere consensus signature is derived from the Govaere bulk cohort
# (GSE135251); display by accession. Keys (panel column values) unchanged.
panel_lab <- c(
  Govaere     = "GSE135251 consensus\nMASLD signature",
  NIDDK       = "NIDDK\nMASLD panel",
  OpenTargets = "Open Targets\nMASLD associations"
)
bench <- bench[panel %in% names(panel_lab)]
bench[, panel_label := panel_lab[panel]]

# Long format: one row per (panel, method).
long <- rbindlist(list(
  bench[, .(panel, panel_label, n_in_panel_atlas,
            method = "Convergence score", auroc = auroc_46d)],
  bench[, .(panel, panel_label, n_in_panel_atlas,
            method = "Prior method", auroc = auroc_46b)]
))

# Order panels by convergence-score AUROC (descending), methods score-first.
panel_order <- bench[order(-auroc_46d), panel_label]
long[, panel_label := factor(panel_label, levels = panel_order)]
long[, method := factor(method, levels = c("Convergence score", "Prior method"))]

# Per-panel label annotating held-out panel size (n positives in atlas).
panel_n <- unique(long[, .(panel_label, n_in_panel_atlas)])

cat("Held-out calibration (AUROC):\n")
print(bench[, .(panel, auroc_46d, auroc_46b)])

# ---------------------------------------------------------------------------
# Panel 5d — grouped bars: convergence score vs prior method per held-out panel
# Semantic colors: convergence score = disease magenta; prior method = gray.
# Chance line at 0.5.
# ---------------------------------------------------------------------------
method_cols <- c("Convergence score" = masld_colors$mash,
                 "Prior method"      = "#9E9E9E")  # control/neutral gray (invariant)

p5d <- ggplot(long, aes(x = panel_label, y = auroc,
                        fill = method, group = method)) +
  geom_hline(yintercept = 0.5, linetype = "dashed",
             linewidth = 0.3, color = "grey55") +
  geom_col(position = position_dodge(width = 0.72),
           width = 0.66, alpha = 0.95) +
  geom_text(aes(label = sprintf("%.2f", auroc)),
            position = position_dodge(width = 0.72),
            vjust = -0.4, size = 2.0, family = "Helvetica") +
  annotate("text", x = 3.42, y = 0.47, label = "chance (0.50)",
           hjust = 1, vjust = 1, size = 1.9, color = "grey45",
           fontface = "italic", family = "Helvetica") +
  scale_fill_manual(values = method_cols, name = NULL) +
  scale_y_continuous(limits = c(0, 1.0),
                     breaks = seq(0, 1, 0.25),
                     expand = expansion(mult = c(0, 0.08))) +
  # PI directive (2026-06-11): short single title line only. The framing
  # sentence ("convergence score predicts held-out panels and improves over
  # OUR prior method"; Govaere/NIDDK/Open Targets are held-out truth panels,
  # NOT competitors) lives in the figure caption, not on the panel.
  labs(x = NULL,
       y = "Held-out AUROC",
       title = "d",
       subtitle = "Held-out AUROC") +
  theme_masld(base_size = 7) +
  theme(plot.title = element_text(face = "bold", size = 10),
        plot.subtitle = element_text(size = 6.3, color = "grey30"),
        axis.text.x = element_text(size = 6.2, lineheight = 0.85),
        legend.position = c(0.99, 0.99),
        legend.justification = c(1, 1),
        legend.background = element_rect(fill = scales::alpha("white", 0.7),
                                         color = NA),
        legend.key.size = unit(0.28, "cm"),
        legend.text = element_text(size = 6))

ggsave(file.path(PANDIR, "fig5d_calibration.pdf"), p5d,
       width = 3.6, height = 3.3, device = cairo_pdf)
cat("Saved: panels/fig5d_calibration.pdf\n")

# Sidecar CSV for caption transparency.
out_csv <- bench[, .(panel, n_positives_in_atlas = n_in_panel_atlas,
                     auroc_convergence_score = auroc_46d,
                     auroc_prior_method = auroc_46b,
                     p_convergence_score = p_46d,
                     p_prior_method = p_46b)]
fwrite(out_csv, file.path(PANDIR, "fig5d_calibration_source.csv"))
cat("Saved: panels/fig5d_calibration_source.csv\n")

cat("\nDONE.\n")
