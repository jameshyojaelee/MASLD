#!/usr/bin/env Rscript
# figS_coloc_window_sensitivity.R
# Supplementary Figure: COLOC window size sensitivity
#
# 4 panels:
#   (a) Gene count bar chart: PP.H4>0.5/0.8/0.9 at 250kb, 500kb, 1Mb
#   (b) Retention rate: % of 1Mb PP.H4>0.5 genes surviving at narrower windows
#   (c) Scatter: PP.H4 at 500kb vs 1Mb (all gene-GWAS pairs with PP.H4>0.3)
#   (d) Scatter: PP.H4 at 250kb vs 1Mb
#
# Input:
#   GWAS/finemapping/results/coloc_threshold_audit/window_sensitivity.csv
# Output:
#   figures/supplementary/figS04_coloc/figS_coloc_window_sensitivity.pdf

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
  library(scales)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

IN_FILE  <- file.path(BASE, "GWAS/finemapping/results/coloc_threshold_audit/window_sensitivity.csv")
OUT_FILE <- file.path(FIGS04_DIR, "figS_coloc_window_sensitivity.pdf")
dir.create(FIGS04_DIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(IN_FILE)) stop("Run 06b_coloc_window_sensitivity.R first: ", IN_FILE)

dat <- fread(IN_FILE)
cat("Loaded:", nrow(dat), "rows\n")

# Window labels (ordered narrow -> wide)
win_labels <- c("250" = "±250 kb", "500" = "±500 kb", "1000" = "±1,000 kb (current)")
win_order  <- c("±250 kb", "±500 kb", "±1,000 kb (current)")
dat[, window_label := factor(win_labels[as.character(half_window_kb)], levels = win_order)]

# Gene-level max PP.H4 per window
gene_level <- dat[!is.na(PP.H4.abf),
                  .(max_pp4 = max(PP.H4.abf, na.rm = TRUE)),
                  by = .(gene, window_label)]

# ---------------------------------------------------------------------------
# Panel A: Gene count bar chart
# ---------------------------------------------------------------------------
counts <- rbindlist(lapply(c(0.5, 0.8, 0.9), function(thresh) {
  gene_level[, .(
    n_genes = sum(max_pp4 > thresh, na.rm = TRUE),
    threshold = paste0("PP.H4 > ", thresh)
  ), by = window_label]
}))
counts[, threshold := factor(threshold,
  levels = c("PP.H4 > 0.9", "PP.H4 > 0.8", "PP.H4 > 0.5"))]

pA <- ggplot(counts, aes(x = window_label, y = n_genes, fill = threshold)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_text(aes(label = n_genes),
            position = position_dodge(width = 0.75),
            vjust = -0.3, size = GEOM_TEXT_6PT, color = "black") +
  scale_fill_manual(
    values = c("PP.H4 > 0.5" = "#BBDEFB",
               "PP.H4 > 0.8" = "#1565C0",
               "PP.H4 > 0.9" = "#0D47A1"),
    name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "COLOC window (half-width)", y = "Unique genes") +
  theme_masld(base_size = 6) +
  theme(legend.position = "bottom",
        legend.key.size = unit(0.25, "cm"),
        axis.text.x = element_text(size = 6))

# ---------------------------------------------------------------------------
# Panel B: Retention of original 1Mb PP.H4>0.5 hits
# ---------------------------------------------------------------------------
orig_genes <- unique(dat[pp4_1mb > 0.5, gene])
n_orig <- length(orig_genes)

retention <- dat[gene %in% orig_genes & !is.na(PP.H4.abf),
                 .(max_pp4 = max(PP.H4.abf, na.rm = TRUE)),
                 by = .(gene, window_label)]
ret_summary <- retention[, .(
  retained_05 = sum(max_pp4 > 0.5) / n_orig * 100,
  retained_08 = sum(max_pp4 > 0.8) / n_orig * 100,
  retained_09 = sum(max_pp4 > 0.9) / n_orig * 100
), by = window_label]

ret_long <- melt(ret_summary, id.vars = "window_label",
                 variable.name = "threshold", value.name = "pct_retained")
ret_long[, threshold := factor(
  fcase(threshold == "retained_05", "PP.H4 > 0.5",
        threshold == "retained_08", "PP.H4 > 0.8",
        threshold == "retained_09", "PP.H4 > 0.9"),
  levels = c("PP.H4 > 0.5", "PP.H4 > 0.8", "PP.H4 > 0.9"))]

pB <- ggplot(ret_long, aes(x = window_label, y = pct_retained,
                            color = threshold, group = threshold)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2) +
  geom_text(aes(label = sprintf("%.0f%%", pct_retained)),
            vjust = -0.6, size = GEOM_TEXT_6PT, show.legend = FALSE) +
  scale_color_manual(
    values = c("PP.H4 > 0.5" = "#90CAF9",
               "PP.H4 > 0.8" = "#1565C0",
               "PP.H4 > 0.9" = "#0D47A1"),
    name = NULL) +
  scale_y_continuous(limits = c(0, 105), labels = function(x) paste0(x, "%")) +
  labs(x = "COLOC window (half-width)",
       y = paste0("% of 1Mb hits retained (n=", n_orig, ")")) +
  theme_masld(base_size = 6) +
  theme(legend.position = "bottom",
        legend.key.size = unit(0.25, "cm"),
        axis.text.x = element_text(size = 6))

# ---------------------------------------------------------------------------
# Panels C & D: Scatter plots (pair-level PP.H4)
# ---------------------------------------------------------------------------
make_scatter <- function(half_window_kb_val, window_label_str) {
  # Get max PP.H4 per gene per window
  w_dat <- dat[half_window_kb == half_window_kb_val & !is.na(PP.H4.abf),
               .(pp4_narrow = max(PP.H4.abf, na.rm = TRUE),
                 pp4_1mb    = max(pp4_1mb, na.rm = TRUE)),
               by = gene]

  # Pearson and Spearman correlation
  r_p <- round(cor(w_dat$pp4_narrow, w_dat$pp4_1mb, use = "complete.obs"), 3)
  r_s <- round(cor(w_dat$pp4_narrow, w_dat$pp4_1mb,
                   use = "complete.obs", method = "spearman"), 3)
  n   <- nrow(w_dat)

  # Colour by whether gene crosses PP.H4>0.5 threshold
  w_dat[, status := fcase(
    pp4_1mb > 0.5 & pp4_narrow > 0.5, "Both >0.5",
    pp4_1mb > 0.5 & pp4_narrow <= 0.5, "Lost at narrow",
    pp4_1mb <= 0.5 & pp4_narrow > 0.5, "Gained at narrow",
    default = "Both <0.5"
  )]
  w_dat[, status := factor(status,
    levels = c("Both >0.5", "Lost at narrow", "Gained at narrow", "Both <0.5"))]

  status_cols <- c(
    "Both >0.5"       = "#1565C0",
    "Lost at narrow"  = "#C62828",
    "Gained at narrow"= "#2E7D32",
    "Both <0.5"       = "#BDBDBD"
  )

  ggplot(w_dat, aes(x = pp4_1mb, y = pp4_narrow, color = status)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                color = "gray60", linewidth = 0.4) +
    geom_hline(yintercept = 0.5, linetype = "dotted", color = "gray50",
               linewidth = 0.3) +
    geom_vline(xintercept = 0.5, linetype = "dotted", color = "gray50",
               linewidth = 0.3) +
    geom_point(size = 0.8, alpha = 0.7) +
    scale_color_manual(values = status_cols, name = NULL,
                       guide = guide_legend(override.aes = list(size = 2))) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    annotate("text", x = 0.02, y = 0.97,
             label = sprintf("r = %.3f\nρ = %.3f\nn = %d", r_p, r_s, n),
             hjust = 0, vjust = 1, size = GEOM_TEXT_6PT, color = "black") +
    labs(x = "PP.H4 (±1,000 kb, current)",
         y = sprintf("PP.H4 (%s)", window_label_str)) +
    theme_masld(base_size = 6) +
    theme(legend.position = "bottom",
          legend.key.size = unit(0.2, "cm"),
          legend.text = element_text(size = 6))
}

pC <- make_scatter(500,  "±500 kb")
pD <- make_scatter(250,  "±250 kb")

# ---------------------------------------------------------------------------
# Assemble 2×2 figure
# ---------------------------------------------------------------------------
combined <- plot_grid(
  pA, pB, pC, pD,
  labels     = c("a", "b", "c", "d"),
  label_size = 9, label_fontface = "plain",
  ncol = 2,
  rel_heights = c(1, 1.1)
)

save_fig_tall(combined, OUT_FILE, width = fig_full_width, height = 7)
cat("Saved:", OUT_FILE, "\n")
message("[caption] a: gene counts by COLOC window size. b: retention of PP.H4>0.5 hits at narrower windows. c-d: PP.H4 at 500kb/250kb vs 1Mb (window: 500kb vs 1Mb / 250kb vs 1Mb).")

# ---------------------------------------------------------------------------
# Print summary to stdout (also appears in methods)
# ---------------------------------------------------------------------------
cat("\n=== Summary for methods ===\n")
for (w in c(250, 500, 1000)) {
  gl <- gene_level[window_label == win_labels[as.character(w)]]
  cat(sprintf("±%dkb: PP.H4>0.5: %d | >0.8: %d | >0.9: %d\n",
              w, gl[max_pp4>0.5,.N], gl[max_pp4>0.8,.N], gl[max_pp4>0.9,.N]))
}
cat(sprintf("\nRetention of %d original 1Mb PP.H4>0.5 genes:\n", n_orig))
print(ret_summary)
