#!/usr/bin/env Rscript
# KEY MESSAGE: Binary MASLD-vs-Control DA recapitulates the F-stage-stratified
# signal — Hep gives 694 covariate-corrected sig peaks; other CTs underpowered.

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

# ── Panel a: Hep MASLD-vs-Control volcano (covariate-corrected, 694 sig) ──
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

# ── Panel b: Per-CT sig peaks comparison (F0 vs F4  vs  MASLD vs Control) ──
ct_short <- c("Hep", "Mac", "Fib", "Endo", "Chol")
ct_full  <- c("Hepatocytes", "Macrophages", "Fibroblasts",
              "Endothelial", "Cholangiocytes")

stage_da_dir <- file.path(BASE, "Analysis/ATAC/Human_Multiome/results/stage_da")
dc <- fread(file.path(BASE,
  "Analysis/ATAC/Human_Multiome/results/snapatac2/scatac_da_results.csv"))
setnames(dc, c("log2(fold_change)", "p-value", "adjusted p-value", "cell_type"),
              c("log2FC", "pvalue", "padj", "cell_type"))
dc_lookup <- c(Hepatocyte = "Hep", Macrophage = "Mac",
               Fibroblast = "Fib", Endothelial = "Endo",
               Cholangiocyte = "Chol")
dc[, ct := dc_lookup[cell_type]]

counts <- rbindlist(lapply(seq_along(ct_short), function(i) {
  ct_i <- ct_short[i]
  # F0 vs F4
  f <- file.path(stage_da_dir, sprintf("da_F0_vs_F4_%s.csv", ct_i))
  f04 <- if (file.exists(f)) sum(fread(f)$padj < 0.05, na.rm = TRUE) else 0L
  # MASLD vs Control (uncorrected, all CTs) — Hep replaced with corrected below
  dvc <- sum(dc[ct == ct_i]$padj < 0.05, na.rm = TRUE)
  data.table(cell_type = ct_full[i], short = ct_i,
             `F0 vs F4` = f04,
             `MASLD vs Control` = dvc)
}))
# Replace Hep MASLD vs Control with covariate-corrected count (694)
counts[short == "Hep", `MASLD vs Control` := 694L]

long <- melt(counts, id.vars = c("cell_type", "short"),
             variable.name = "contrast", value.name = "n")
long[, cell_type := factor(cell_type, levels = rev(ct_full))]
long[, contrast := factor(contrast, levels = c("F0 vs F4", "MASLD vs Control"))]

p_b <- ggplot(long, aes(y = cell_type, x = n, fill = contrast)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) +
  scale_fill_manual(values = c("F0 vs F4" = GRAY,
                               "MASLD vs Control" = MAG),
                    name = NULL) +
  scale_x_continuous(name = "Significant DA peaks",
                     expand = expansion(mult = c(0, 0.05))) +
  labs(tag = "b", y = NULL) +
  theme_masld(base_size = 7) +
  theme_pub() +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(face = "bold", colour = "black"),
        axis.title.x = element_text(face = "bold"),
        legend.position = "bottom",
        legend.key.size = unit(0.25, "cm"),
        plot.tag = element_text(size = 11, face = "bold"),
        plot.tag.position = c(0.02, 0.97))

combo <- p_a + p_b + plot_layout(widths = c(3, 2.5))
ggsave(OUT_PDF, combo, width = 6.4, height = 2.8, device = cairo_pdf)
cat(sprintf("[saved] %s\n", OUT_PDF))
cat(sprintf("  Hep MASLD vs Control (corrected): %d opening / %d closing\n",
            n_open, n_close))
print(counts)
