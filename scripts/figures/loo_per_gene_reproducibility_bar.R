#!/usr/bin/env Rscript
# ============================================================================
# loo_per_gene_reproducibility_bar.R
# SUPP (robustness) — per-gene leave-one-cohort-out reproducibility
#
# Tier-1 integrated DEGs binned by the number of 5-fold leave-one-cohort-out
# refits in which each gene stays significant (padj<0.1). Converts the
# aggregate 81.6% mean recovery into a per-gene story: most DEGs survive ALL
# five LOO refits (~77%) and the vast majority survive >=4 (~94%).
#
# Output: <FIGS_ROBUST_DIR>/loo_per_gene_reproducibility_bar.pdf  (70 x 72 mm)
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

dir.create(FIGS_ROBUST_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_PDF <- file.path(FIGS_ROBUST_DIR, "loo_per_gene_reproducibility_bar.pdf")
OUT_CSV <- file.path(FIGS_ROBUST_DIR, "loo_per_gene_reproducibility_bar.csv")

INT <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")

# ---------------------------------------------------------------------------
# Per-gene LOO stability + Tier-1 DEG membership
# ---------------------------------------------------------------------------
stab <- fread(file.path(INT, "loo_cv_C2/loo_cv_C2_per_gene_stability.csv"))
deg  <- fread(file.path(INT, "canonical_deg_results.csv"))

deg[, is_deg := padj < 0.05 & abs(logFC) > 0.5]
deg_genes <- deg[is_deg == TRUE, gene]

x <- stab[gene %in% deg_genes]
n_deg <- nrow(x)
cat(sprintf("[info] Tier-1 DEGs with LOO stability records: %d (of %d Tier-1)\n",
            n_deg, length(deg_genes)))

# bin by folds-significant (clamp to 0..5; restrict to 5-fold-tested genes)
x[, nsig := pmin(n_folds_sig_01, 5L)]
bins <- x[, .(n = .N), by = nsig][order(-nsig)]
bins <- merge(data.table(nsig = 5:0), bins, by = "nsig", all.x = TRUE)
bins[is.na(n), n := 0]
bins[, frac := n / n_deg]

# hero cumulative fractions
all5  <- bins[nsig == 5, frac]
ge4   <- bins[nsig >= 4, sum(frac)]
ge3   <- bins[nsig >= 3, sum(frac)]
cat(sprintf("[hero] DEGs significant in ALL 5 LOO folds: %.1f%% (%d/%d)\n",
            100 * all5, bins[nsig == 5, n], n_deg))
cat(sprintf("[hero] DEGs significant in >=4 LOO folds: %.1f%%\n", 100 * ge4))
cat(sprintf("[hero] DEGs significant in >=3 LOO folds: %.1f%%\n", 100 * ge3))

fwrite(bins[, .(n_folds_sig = nsig, n_genes = n, fraction = frac)], OUT_CSV)

# ---------------------------------------------------------------------------
# Stacked single-bar: green->magenta gradient by reproducibility tier
# ---------------------------------------------------------------------------
plot_dt <- bins[nsig >= 1]   # drop 0-fold (genes lost in every refit ~ none)
plot_dt[, lab := factor(paste0(nsig, "/5"),
                        levels = paste0(5:1, "/5"))]
plot_dt[, bar := "Tier-1 DEGs"]

# darkest = most reproducible (all 5)
grad <- c("5/5" = "#1565C0", "4/5" = "#5B8FCB", "3/5" = "#9E9E9E",
          "2/5" = "#E0809F", "1/5" = "#C9265E")

plot_dt[, pct_lab := ifelse(frac >= 0.04, sprintf("%.1f%%", 100 * frac), "")]

p <- ggplot(plot_dt, aes(x = bar, y = frac, fill = lab)) +
  geom_col(width = 0.55, color = "white", linewidth = 0.4) +
  geom_text(aes(label = pct_lab), position = position_stack(vjust = 0.5),
            size = 2.2, color = "white", fontface = "bold") +
  # cumulative annotations
  annotate("segment", x = 1.32, xend = 1.32, y = 1 - all5, yend = 1,
           color = "#1565C0", linewidth = 0.5) +
  annotate("text", x = 1.40, y = 1 - all5 / 2,
           label = sprintf("all 5 folds\n%.1f%%", 100 * all5),
           hjust = 0, size = 2.1, color = "#1565C0", lineheight = 0.9) +
  annotate("segment", x = 1.62, xend = 1.62, y = 1 - ge4, yend = 1,
           color = "#5B8FCB", linewidth = 0.5) +
  annotate("text", x = 1.70, y = 1 - ge4 / 2,
           label = sprintf(">=4 folds\n%.1f%%", 100 * ge4),
           hjust = 0, size = 2.1, color = "#3A6FA8", lineheight = 0.9) +
  scale_fill_manual(values = grad, name = "LOO folds\nsignificant",
                    breaks = paste0(5:1, "/5")) +
  scale_y_continuous(name = "Fraction of Tier-1 DEGs",
                     labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(xlim = c(0.7, 2.2)) +
  labs(title = "Most integrated DEGs are reproducible per gene",
       subtitle = sprintf("%s Tier-1 DEGs, 5-fold leave-one-cohort-out refit (padj<0.1)",
                          format(n_deg, big.mark = ",")),
       x = NULL) +
  theme_masld(base_size = 7) +
  theme(
    plot.title      = element_text(size = 7.3, face = "bold", margin = margin(b = 2)),
    plot.subtitle   = element_text(size = 5.3, color = "#555555", margin = margin(b = 5)),
    axis.text.x     = element_blank(),
    axis.ticks.x    = element_blank(),
    legend.position = "left",
    legend.key.size = unit(0.2, "cm"),
    legend.text     = element_text(size = 5.5),
    legend.title    = element_text(size = 5.8)
  )

ggsave(OUT_PDF, p, width = 70 / 25.4, height = 72 / 25.4,
       units = "in", device = cairo_pdf)
cat(sprintf("[saved] %s\n", OUT_PDF))
