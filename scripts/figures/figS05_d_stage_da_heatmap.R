#!/usr/bin/env Rscript
# figS05_d_stage_da_heatmap.R
# COMPOSITION-CONFOUND DIAGNOSTIC (not a graded-regulation claim). Hepatocyte
# stage-stratified differential accessibility, opening (gaining) vs closing
# (losing) peaks per contrast. The signal is confined to the EXTREME F0<->F4
# endpoint gap (605 peaks, ~74% CLOSING) and the ordinal models that span it;
# the single-step ADJACENT contrasts (F0->F3, F3->F4, F2+F3->F4) are ~0. A
# signal that appears only across the widest stage gap and is dominated by
# global closing is the signature of hepatocyte compositional loss (cell
# dropout), NOT graded locus-specific regulation — the same confound that led
# to the epigenetic panels being CUT from main Fig 4. Counts read live from
# stage_da/*.csv (no hardcoding).

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

DA_DIR <- file.path(BASE, "Analysis/ATAC/Human_Multiome/results/stage_da")
PADJ   <- 0.05

cts   <- c("Hep")
# Contrasts ordered to expose the confound: single-step ADJACENT transitions
# first (these are ~0), then the EXTREME endpoint gap + ordinal models that span
# F0<->F4 (these carry the signal — driven by global closing = composition).
tests <- c("F0_vs_F3", "F3_vs_F4", "F2F3_vs_F4",
           "F0_vs_F2F3", "F0_vs_F4", "ordinal_F0_F3_F4", "stage_ordinal")
test_labels <- c("F0 vs F3",
                 "F3 vs F4",
                 "F2+F3 vs F4",
                 "F0 vs F2+F3",
                 "F0 vs F4",
                 "Ordinal F0/F3/F4",
                 "Ordinal (all stages)")

# ── load all DA results ────────────────────────────────────────────────────────
# IMPORTANT: always emit a row per contrast (0 when no sig peaks) so that the
# ADJACENT contrasts that are NULL still appear on the panel — that 0 IS the
# point (graded transitions show no signal). Never silently drop them.
rows <- rbindlist(lapply(seq_along(tests), function(ti) {
  rbindlist(lapply(cts, function(ct) {
    f <- file.path(DA_DIR, sprintf("da_%s_%s.csv", tests[ti], ct))
    d <- if (file.exists(f))
           tryCatch(fread(f, select = c("log2FC", "padj")), error = function(e) NULL)
         else NULL
    if (!is.null(d) && "padj" %in% names(d)) d <- d[!is.na(padj) & padj < PADJ] else d <- NULL
    data.table(cell_type = ct,
               contrast  = test_labels[ti],
               n_open  = if (!is.null(d)) sum(d$log2FC > 0) else 0L,
               n_close = if (!is.null(d)) sum(d$log2FC < 0) else 0L)
  }), fill = TRUE)
}), fill = TRUE)

cat(sprintf("[loaded] %d cell-type × contrast rows (incl. zero-signal contrasts)\n", nrow(rows)))
keep_cts <- cts   # Hep only; every contrast retained so adjacent-0 stays visible

# melt to long for bidirectional bars
bar <- rbindlist(list(
  rows[, .(cell_type, contrast, n =  n_open,  direction = "Opening")],
  rows[, .(cell_type, contrast, n = -n_close, direction = "Closing")]
))

# factor order: adjacent (left) -> extreme endpoint / ordinal (right)
bar[, contrast   := factor(contrast, levels = test_labels)]
bar[, cell_type  := factor(cell_type, levels = keep_cts)]
bar[, direction  := factor(direction, levels = c("Opening", "Closing"))]

# per-contrast total + %-closing (number, shown above each bar group)
tot <- rows[, .(total = n_open + n_close,
                pct_close = 100 * n_close / pmax(1, n_open + n_close)),
            by = .(cell_type, contrast)]
tot[, contrast := factor(contrast, levels = test_labels)]
tot[, lab := fifelse(total >= 20, sprintf("%d\n(%.0f%% closing)", total, pct_close),
                     as.character(total))]
# label sits just above the opening (positive) bar top for each contrast
tot[, ypos := rows$n_open[match(paste(cell_type, contrast),
              paste(rows$cell_type, rows$contrast))]]

# ── plot ──────────────────────────────────────────────────────────────────────
p <- ggplot(bar, aes(x = contrast, y = n, fill = direction)) +
  geom_col(width = 0.72) +
  geom_hline(yintercept = 0, linewidth = 0.35, color = "gray60") +
  geom_text(data = tot, aes(x = contrast, y = ypos, label = lab),
            inherit.aes = FALSE, vjust = -0.35, size = PUB_GEOM_TEXT,
            colour = "black", lineheight = 0.85) +
  facet_wrap(~ cell_type, nrow = 1, scales = "free_y") +
  scale_fill_manual(
    values = c(Opening = masld_colors$up,    # magenta — gaining accessibility
               Closing = masld_colors$down), # blue   — losing accessibility
    name = NULL
  ) +
  scale_y_continuous(labels = function(x) comma(abs(x)),
                     expand = expansion(mult = c(0.12, 0.28))) +
  labs(x = NULL, y = "DA peaks  (padj < 0.05)") +
  theme_masld(base_size = 10) +
  theme_pub() +
  theme(strip.text       = element_text(size = 9, face = "bold"),
        axis.text.x      = element_text(size = PUB_AXIS_TEXT, angle = 35, hjust = 1,
                                        colour = "black"),
        axis.title.y     = element_text(size = PUB_AXIS_TITLE, face = "bold"),
        legend.position  = "bottom",
        legend.key.size  = unit(0.3, "cm"),
        panel.grid       = element_blank())

out_pdf <- file.path(FIGS05_DIR, "figS05_d_stage_da_heatmap.pdf")
ggsave(out_pdf, p, width = 4.6, height = 3.0, device = cairo_pdf)
message("Wrote: ", out_pdf)

# ── Caption / provenance to stdout (NOT on the panel) ────────────────────────
adj  <- tot[contrast %in% c("F0 vs F3", "F3 vs F4", "F2+F3 vs F4")]
ext  <- tot[contrast == "F0 vs F4"]
message(strrep("=", 78))
message("FIGURE LEGEND (paste into manuscript; stats live here, not on the panel)")
message(strrep("=", 78))
message(sprintf(
"Hepatocyte stage-stratified differential accessibility is a COMPOSITION-CONFOUND
diagnostic, not evidence of graded locus-specific regulation. Significant DA
peaks (padj<0.05) concentrate in the extreme F0-vs-F4 endpoint gap (%d peaks,
%.0f%% CLOSING) while every single-step adjacent transition is ~0 (F0->F3, F3->F4,
F2+F3->F4: %s). A signal that emerges only across the widest stage gap and is
dominated by global closing is the expected signature of hepatocyte
compositional loss (cell dropout), not graded regulation. n=18 multiome.",
  ext$total[1], ext$pct_close[1],
  paste(sprintf("%s=%d", as.character(adj$contrast), adj$total), collapse = ", ")))
message(strrep("=", 78))

out_csv <- sub("\\.pdf$", ".csv", out_pdf)
fwrite(bar[, .(cell_type, contrast, direction, n = abs(n))], out_csv)
message("Wrote: ", out_csv)
