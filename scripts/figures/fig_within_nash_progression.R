##############################################################################
# Figure: Within-NASH Fibrosis Progression (C12) vs General Fibrosis (C3)
#
# Highlights the novel C12 contrast (early vs late NASH fibrosis) and how it
# compares to the general fibrosis contrast (C3), plus C11 NASH vs Control.
#
# Layout (3 panels side by side, ~12 x 4.5 inches):
#   (A) C12 vs C3 logFC scatter — shared, C12-only, C3-only DEGs
#   (B) C12-unique gene biology — top 20 C12-unique genes by |logFC|
#   (C) C11 volcano (NASH vs Control) — pure NASH initiation
#
# Output:
#   figures/supplementary/figS02_progression/fig_within_nash_progression.pdf
#   figures/supplementary/figS02_progression/panel_nash_scatter.pdf
#   figures/supplementary/figS02_progression/panel_nash_unique.pdf
#   figures/supplementary/figS02_progression/panel_nash_volcano.pdf
##############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(ggrepel)
  library(patchwork)
})

set.seed(42)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

source(file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/visualization/functions/theme_publication.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

PROG_DIR <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/progression")

OUT_DIR   <- FIGS02_DIR
PANEL_DIR <- file.path(OUT_DIR, "panels")
dir.create(PANEL_DIR, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
c12 <- fread(file.path(PROG_DIR, "c12_early_vs_late_nash_dream.csv"))
c3  <- fread(file.path(PROG_DIR, "c3_adv_vs_early_fib_dream.csv"))
c11 <- fread(file.path(PROG_DIR, "c11_nash_vs_ctrl_dream.csv"))

# Standardize column names (use symbol as primary label).
# C2 swap (2026-06-08): the c12/c3/c11 progression contrasts became limma-voom
# (gene,logFC,...,padj,...) and no longer carry a `symbol` column — derive it from
# the gene→symbol map so downstream symbol-based labelling/merging works.
c12 <- add_symbols(c12, "gene")
c3  <- add_symbols(c3,  "gene")
c11 <- add_symbols(c11, "gene")

# ---------------------------------------------------------------------------
# Panel A: C12 vs C3 logFC scatter
# ---------------------------------------------------------------------------
# Merge on gene (Ensembl ID)
merged <- merge(
  c12[, .(gene, logFC_c12 = logFC, padj_c12 = padj, symbol)],
  c3[,  .(gene, logFC_c3  = logFC, padj_c3  = padj)],
  by = "gene"
)

# Classify DEGs (padj < 0.1)
merged[, category := fcase(
  padj_c12 < 0.1 & padj_c3 < 0.1,  "Shared",
  padj_c12 < 0.1 & padj_c3 >= 0.1, "C12-only",
  padj_c12 >= 0.1 & padj_c3 < 0.1,  "C3-only",
  default = "NS"
)]

# Counts
n_shared   <- merged[category == "Shared",   .N]
n_c12_only <- merged[category == "C12-only", .N]
n_c3_only  <- merged[category == "C3-only",  .N]

# Correlation
rho <- cor(merged$logFC_c12, merged$logFC_c3, use = "complete.obs",
           method = "spearman")

# Key genes to label (shared markers)
shared_labels <- c("COL1A1", "ACTA2", "TIMP1", "PPARA", "NR1H4")

# Top 5 C12-unique genes by |logFC|
c12_unique <- merged[category == "C12-only"][order(-abs(logFC_c12))]
top5_c12 <- head(c12_unique$symbol, 5)
# Remove empty symbols
top5_c12 <- top5_c12[nzchar(top5_c12)]
if (length(top5_c12) < 5) {
  # Fill from remaining
  remaining <- c12_unique[nzchar(symbol) & !symbol %in% top5_c12]
  top5_c12 <- c(top5_c12, head(remaining$symbol,
                                5 - length(top5_c12)))
}

# Label column
merged[, label := fifelse(
  symbol %in% c(shared_labels, top5_c12), symbol, ""
)]

# Color palette
cat_colors <- c(
  "Shared"   = sanjana_colors[["Grey"]],
  "C12-only" = sanjana_colors[["Magenta"]],
  "C3-only"  = sanjana_colors[["Blue"]],
  "NS"       = "#E0E0E0"
)

# Annotation text
annot_text <- sprintf(
  "rho = %.3f\nShared: %s\nC12-only: %s\nC3-only: %s",
  rho,
  format(n_shared, big.mark = ","),
  format(n_c12_only, big.mark = ","),
  format(n_c3_only, big.mark = ",")
)

# Plot order: NS first, then shared, then C12-only and C3-only on top
merged[, plot_order := fcase(
  category == "NS",       1L,
  category == "Shared",   2L,
  category == "C3-only",  3L,
  category == "C12-only", 4L,
  default = 0L
)]
setorder(merged, plot_order)

p_a <- ggplot(merged, aes(x = logFC_c3, y = logFC_c12)) +
  geom_point(aes(color = category), size = 0.3, alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey40", linewidth = 0.3) +
  geom_text_repel(
    data = merged[nzchar(label)],
    aes(label = label, color = category),
    size = 2, fontface = "italic",
    max.overlaps = 20, segment.size = 0.2,
    min.segment.length = 0, seed = 42,
    box.padding = 0.3
  ) +
  scale_color_manual(values = cat_colors, name = "DEG status") +
  annotate("text", x = -Inf, y = Inf, label = annot_text,
           hjust = -0.05, vjust = 1.1, size = 2, color = "grey30") +
  labs(
    x = "logFC (C3: Advanced vs Early Fibrosis)",
    y = "logFC (C12: Within-NASH Fibrosis)",
    title = "A"
  ) +
  coord_cartesian(
    xlim = c(-1.5, 1.5),
    ylim = c(-1.5, 1.5)
  ) +
  theme_publication(base_size = 7) +
  theme(
    plot.title = element_text(size = 10, face = "plain"),
    legend.position = "bottom",
    legend.key.size = unit(0.5, "lines")
  ) +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1)))

# ---------------------------------------------------------------------------
# Panel B: C12-unique gene biology — top 20 by |logFC|
# ---------------------------------------------------------------------------
# Get all C12-unique DEGs
c12_unique_all <- merged[category == "C12-only" & nzchar(symbol)]
c12_unique_all[, abs_lfc := abs(logFC_c12)]
setorder(c12_unique_all, -abs_lfc)

top20 <- head(c12_unique_all, 20)

# Direction coloring
top20[, direction := fifelse(logFC_c12 > 0, "Up in late NASH", "Down in late NASH")]

# Factor for ordered lollipop
top20[, symbol := factor(symbol, levels = rev(top20$symbol))]

dir_colors <- c(
  "Up in late NASH"   = sanjana_colors[["Magenta"]],
  "Down in late NASH" = sanjana_colors[["Blue"]]
)

p_b <- ggplot(top20, aes(x = logFC_c12, y = symbol, color = direction)) +
  geom_segment(aes(x = 0, xend = logFC_c12, y = symbol, yend = symbol),
               linewidth = 0.4) +
  geom_point(size = 1.5) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey40") +
  scale_color_manual(values = dir_colors, name = "Direction") +
  labs(
    x = "logFC (C12: Within-NASH Fibrosis)",
    y = NULL,
    title = "B",
    subtitle = sprintf("Top 20 C12-unique DEGs (%s total)",
                        format(n_c12_only, big.mark = ","))
  ) +
  theme_publication(base_size = 7) +
  theme(
    plot.title = element_text(size = 10, face = "plain"),
    plot.subtitle = element_text(size = 6, face = "plain"),
    axis.text.y = element_text(face = "italic", size = 5.5),
    legend.position = "bottom",
    legend.key.size = unit(0.5, "lines"),
    panel.grid.major.y = element_blank()
  )

# ---------------------------------------------------------------------------
# Panel C: C11 volcano (NASH vs Control)
# ---------------------------------------------------------------------------
# Significance classification
c11[, sig_class := fcase(
  padj < 0.1 & logFC > 0,  "Up",
  padj < 0.1 & logFC < 0,  "Down",
  default = "NS"
)]

n_up   <- c11[sig_class == "Up",   .N]
n_down <- c11[sig_class == "Down", .N]

# Key NASH markers to highlight
nash_markers <- c("CYP2E1", "PPARA", "THRB", "COL1A1", "IL1B")
c11[, label := fifelse(symbol %in% nash_markers, symbol, "")]

# Clamp -log10 padj for display
c11[, neg_log10_padj := pmin(-log10(padj), 50)]

# Colors
vol_colors <- c(
  "Up"   = sanjana_colors[["Magenta"]],
  "Down" = sanjana_colors[["Blue"]],
  "NS"   = sanjana_colors[["Grey"]]
)

# Order for plotting
c11[, plot_order := fcase(
  sig_class == "NS",   1L,
  sig_class == "Up",   2L,
  sig_class == "Down", 2L,
  default = 0L
)]
setorder(c11, plot_order)

p_c <- ggplot(c11, aes(x = logFC, y = neg_log10_padj)) +
  geom_point(aes(color = sig_class), size = 0.3, alpha = 0.5) +
  geom_hline(yintercept = -log10(0.1), linetype = "dashed",
             color = "grey50", linewidth = 0.3) +
  geom_text_repel(
    data = c11[nzchar(label)],
    aes(label = label, color = sig_class),
    size = 2, fontface = "italic",
    max.overlaps = 20, segment.size = 0.2,
    min.segment.length = 0, seed = 42,
    box.padding = 0.4
  ) +
  scale_color_manual(values = vol_colors, name = "Significance") +
  annotate("text",
    x = max(c11$logFC, na.rm = TRUE) * 0.6,
    y = max(c11$neg_log10_padj, na.rm = TRUE) * 0.95,
    label = sprintf("Up: %s\nDown: %s",
                    format(n_up, big.mark = ","),
                    format(n_down, big.mark = ",")),
    hjust = 0, size = 2, color = "grey30"
  ) +
  labs(
    x = "logFC",
    y = expression(-log[10](p[adj])),
    title = "C",
    subtitle = "NASH vs Control (C11)"
  ) +
  theme_publication(base_size = 7) +
  theme(
    plot.title = element_text(size = 10, face = "plain"),
    plot.subtitle = element_text(size = 6, face = "plain"),
    legend.position = "bottom",
    legend.key.size = unit(0.5, "lines")
  ) +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1)))

# ---------------------------------------------------------------------------
# Assemble and save
# ---------------------------------------------------------------------------
combined <- p_a + p_b + p_c +
  plot_layout(ncol = 3, widths = c(1, 0.85, 1))

# Composite figure
save_pdf(combined,
         file.path(OUT_DIR, "fig_within_nash_progression.pdf"),
         width = 12, height = 4.5)

# Individual panels → panels/ subdirectory
save_pdf(p_a, file.path(PANEL_DIR, "panel_nash_scatter.pdf"),
         width = 4.5, height = 4.5)
save_pdf(p_b, file.path(PANEL_DIR, "panel_nash_unique.pdf"),
         width = 3.5, height = 4.5)
save_pdf(p_c, file.path(PANEL_DIR, "panel_nash_volcano.pdf"),
         width = 4.5, height = 4.5)

cat("Done. Output:\n")
cat("  ", file.path(OUT_DIR, "fig_within_nash_progression.pdf"), "\n")
cat("  ", file.path(PANEL_DIR, "panel_nash_scatter.pdf"), "\n")
cat("  ", file.path(PANEL_DIR, "panel_nash_unique.pdf"), "\n")
cat("  ", file.path(PANEL_DIR, "panel_nash_volcano.pdf"), "\n")

# ---------------------------------------------------------------------------
# Summary statistics
# ---------------------------------------------------------------------------
cat("\n--- Summary ---\n")
cat(sprintf("C12 total genes tested: %s\n", format(nrow(c12), big.mark = ",")))
cat(sprintf("C3 total genes tested: %s\n", format(nrow(c3), big.mark = ",")))
cat(sprintf("C11 total genes tested: %s\n", format(nrow(c11), big.mark = ",")))
cat(sprintf("Merged genes (C12 x C3): %s\n", format(nrow(merged), big.mark = ",")))
cat(sprintf("Shared DEGs: %s\n", format(n_shared, big.mark = ",")))
cat(sprintf("C12-only DEGs: %s\n", format(n_c12_only, big.mark = ",")))
cat(sprintf("C3-only DEGs: %s\n", format(n_c3_only, big.mark = ",")))
cat(sprintf("Spearman rho(C12, C3): %.3f\n", rho))
cat(sprintf("C11 Up: %s, Down: %s\n",
            format(n_up, big.mark = ","),
            format(n_down, big.mark = ",")))
