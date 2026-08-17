#!/usr/bin/env Rscript
# fig1g_deg_landscape_cv.R
# Panel A: DEG count heatmap (padj on y-axis, |LFC| on x-axis)
# Panel B: CV of N DEGs across padj thresholds at each LFC cutoff
# The canonical LFC cutoff is the first value where CV drops below 10%.

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})

source("scripts/figures/publication_theme.R")
source("scripts/figures/load_figure_data.R")  # provides FIG2_DIR (= fig3_RNAseq)

PROJ <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
DREAM <- file.path(PROJ, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/canonical_deg_results.csv")
OUTDIR <- file.path(PROJ, "figures/main/fig1_atlas_overview/panels")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
# fig1g_figs3c_deg_landscape_cv.pdf + fig1g_deg_landscape_data.csv relocated to the
# Figure-3 RNA-seq dir (FIG2_DIR = figures/main/fig3_RNAseq, back-compat constant name).
FIG2_PANEL_DIR <- file.path(FIG2_DIR, "panels")
dir.create(FIG2_PANEL_DIR, recursive = TRUE, showWarnings = FALSE)

d <- fread(DREAM)
cat("Loaded:", nrow(d), "genes\n")

lfcs <- seq(0.05, 1.0, by = 0.05)
padjs <- c(0.001, 0.005, 0.01, 0.025, 0.05, 0.1)

grid <- rbindlist(lapply(lfcs, function(lfc) {
  rbindlist(lapply(padjs, function(pa) {
    n <- sum(d$padj < pa & abs(d$logFC) > lfc, na.rm = TRUE)
    data.table(lfc = lfc, padj = pa, n_DEG = n)
  }))
}))

cv_dt <- grid[, .(mean_n = mean(n_DEG), sd_n = sd(n_DEG),
                   cv = sd(n_DEG) / mean(n_DEG) * 100), by = lfc]

# Find first LFC where CV drops below 10%
first_below_10 <- cv_dt[cv < 10, lfc][1]
cv_at_threshold <- cv_dt[lfc == first_below_10, cv]
n_at_threshold <- grid[lfc == first_below_10 & padj == 0.05, n_DEG]

cat("First LFC with CV < 10%:", first_below_10, "(CV =", round(cv_at_threshold, 1), "%)\n")
cat("DEGs at padj<0.05, |LFC|>", first_below_10, ":", n_at_threshold, "\n")

# Panel A: Heatmap
pA <- ggplot(grid, aes(x = factor(lfc), y = factor(padj), fill = n_DEG)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = formatC(n_DEG, format = "d", big.mark = ",")),
            size = GEOM_TEXT_6PT, color = "black") +
  scale_fill_gradient2(low = "white", mid = "#4393C3", high = "#2166AC",
                       midpoint = max(grid$n_DEG) / 2, name = "DEGs") +
  geom_vline(xintercept = which(levels(factor(lfcs)) == as.character(first_below_10)),
             linetype = "dashed", color = "red", linewidth = 0.6) +
  labs(x = "|log2FC| cutoff", y = "padj cutoff") +
  theme_minimal(base_size = 6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())

# Panel B: CV line with threshold annotation
pB <- ggplot(cv_dt, aes(x = lfc, y = cv)) +
  geom_line(linewidth = 0.8, color = "#2166AC") +
  geom_point(size = 1.5, color = "#2166AC") +
  geom_hline(yintercept = 10, linetype = "dashed", color = "#999999") +
  geom_vline(xintercept = first_below_10, linetype = "dashed", color = "red") +
  geom_point(data = cv_dt[lfc == first_below_10], size = 3.5, color = "red", shape = 1, stroke = 1.2) +
  annotate("text", x = first_below_10 + 0.05, y = cv_at_threshold + 3,
           label = sprintf("|LFC| = %.2f\nCV = %.1f%%\n%s DEGs (padj<0.05)",
                           first_below_10, cv_at_threshold,
                           formatC(n_at_threshold, big.mark = ",")),
           size = GEOM_TEXT_6PT, hjust = 0, color = "red") +
  annotate("text", x = 0.85, y = 11.5, label = "CV = 10%", size = GEOM_TEXT_6PT,
           color = "#999999") +
  scale_x_continuous(breaks = seq(0.1, 1.0, by = 0.1)) +
  labs(x = "|log2FC| cutoff",
       y = "CV of DEG count\nacross padj thresholds (%)") +
  theme_minimal(base_size = 6)

message("[caption] Panel A: DEG count by threshold. Panel B: threshold stability (CV of DEG count across padj thresholds).")
combined <- pA / pB + plot_layout(heights = c(2, 1))

out_pdf <- file.path(FIG2_PANEL_DIR, "figs3c_deg_landscape_cv.pdf")
ggsave(out_pdf, combined, width = 7, height = 6)
cat("Saved:", out_pdf, "\n")

out_csv <- file.path(FIG2_PANEL_DIR, "data", "deg_landscape_cv_data.csv")
fwrite(grid, out_csv)

out_cv <- file.path(OUTDIR, "fig1g_cv_data.csv")  # not relocated; stays in fig1
fwrite(cv_dt, out_cv)
cat("Data:", out_csv, out_cv, "\n")
