##############################################################################
# Supplementary Figure: DEG Threshold Sensitivity Landscape
# 3 panels: (a) DEG count heatmap (padj x |log2FC|),
#           (c) Marginal trade-off at padj = 0.05  [DEG count vs control recovery],
#           (d) Standard padj+logFC vs ashr framework concordance.
# Sensitivity sweep across the (padj, |log2FC|) grid. The canonical operating
# point (padj < 0.05, |LFC| > 0.5; LOO-CV stability + ~log2(1.5) convention)
# is highlighted for reference -- this panel does NOT claim that cell is
# "best" by any objective criterion. The cutoff is justified independently
# by the LOO-CV stability analysis
# (figures/main/fig1_atlas_overview/panels/fig1g.pdf).
##############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(patchwork)
  library(scales)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUT_DIR <- FIGS_SENS_DIR

# ═══════════════════════════════════════════════════════════════════════════
# Load data
# ═══════════════════════════════════════════════════════════════════════════
message("Loading dream results (raw padj + logFC)...")
dream <- fread(file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration",
  "dream_results_ashr.csv"),
  select = c("gene", "logFC", "padj", "symbol", "shrunk_logFC", "lfsr"))
message("  ", nrow(dream), " genes loaded")

# Positive controls
posctrl <- fread(file.path(BASE, "results/library/positive_control.csv"))
ctrl_symbols <- unique(posctrl$`Gene symbol`)
n_ctrl_total <- length(ctrl_symbols)
message("  ", n_ctrl_total, " positive control genes")

# How many controls are in dream at all?
n_ctrl_in_dream <- sum(ctrl_symbols %in% dream$symbol)
message("  ", n_ctrl_in_dream, " controls found in dream results")

# ═══════════════════════════════════════════════════════════════════════════
# Threshold sweep — padj x |logFC| (RAW, not shrunk)
# ═══════════════════════════════════════════════════════════════════════════
padj_vals <- c(1e-3, 0.005, 0.01, 0.025, 0.05, 0.1)
lfc_vals  <- c(0, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0, 1.5)

message("Running padj x |logFC| sweep (", length(padj_vals), " x ", length(lfc_vals), ")...")

sweep <- rbindlist(lapply(padj_vals, function(pa) {
  rbindlist(lapply(lfc_vals, function(lf) {
    hits <- dream[padj < pa & abs(logFC) >= lf]
    n_degs <- nrow(hits)
    n_ctrl <- sum(hits$symbol %in% ctrl_symbols)
    pct_ctrl <- 100 * n_ctrl / n_ctrl_in_dream  # % of recoverable controls
    data.table(padj_threshold = pa, lfc_threshold = lf,
               n_degs = n_degs, n_ctrl = n_ctrl, pct_ctrl = pct_ctrl)
  }))
}))

# Also do ashr sweep for concordance comparison
ashr_padj_vals <- c(0.01, 0.025, 0.05, 0.1)
ashr_lfc_vals  <- c(0, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0)

sweep_ashr <- rbindlist(lapply(ashr_padj_vals, function(pa) {
  rbindlist(lapply(ashr_lfc_vals, function(lf) {
    hits <- dream[lfsr < pa & abs(shrunk_logFC) >= lf]
    data.table(padj_threshold = pa, lfc_threshold = lf,
               n_degs = nrow(hits), framework = "ashr (lfsr + shrunk LFC)")
  }))
}))

# padj equivalent for comparison
sweep_padj_compare <- rbindlist(lapply(ashr_padj_vals, function(pa) {
  rbindlist(lapply(ashr_lfc_vals, function(lf) {
    hits <- dream[padj < pa & abs(logFC) >= lf]
    data.table(padj_threshold = pa, lfc_threshold = lf,
               n_degs = nrow(hits), framework = "Standard (padj + logFC)")
  }))
}))

# Format axis labels
sweep[, padj_label := factor(
  ifelse(padj_threshold < 0.001, format(padj_threshold, scientific = TRUE),
         as.character(padj_threshold)),
  levels = ifelse(padj_vals < 0.001, format(padj_vals, scientific = TRUE),
                  as.character(padj_vals))
)]
sweep[, lfc_label := factor(
  sprintf("%.1f", lfc_threshold),
  levels = sprintf("%.1f", lfc_vals)
)]

# DEG count labels (compact: 1.4k, 7.6k, etc.)
sweep[, deg_label := ifelse(n_degs >= 1000,
                            sprintf("%.1fk", n_degs / 1000),
                            as.character(n_degs))]

# ═══════════════════════════════════════════════════════════════════════════
# Panel (a): DEG Count Heatmap (|log2FC| on x, padj on y) + CV curve below
# ═══════════════════════════════════════════════════════════════════════════
message("Panel (a): DEG count heatmap + CV across padj...")

tier1_x_a <- which(levels(sweep$lfc_label) == sprintf("%.1f", 0.5))
n_padj_a  <- length(levels(sweep$padj_label))

p_a_hm <- ggplot(sweep, aes(x = lfc_label, y = padj_label, fill = n_degs)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = deg_label),
            size = 2, color = "gray15", fontface = "bold") +
  # Highlight canonical operating point (padj<0.05, |LFC|>0.5).
  # Note: this is a reference marker, not an "optimum" -- the canonical
  # cutoff is justified by LOO-CV stability (Fig 1g), not by this panel.
  annotate("rect",
           xmin = tier1_x_a - 0.5, xmax = tier1_x_a + 0.5,
           ymin = 0.5, ymax = n_padj_a + 0.5,
           color = "#FFB300", fill = NA, linewidth = 0.9) +
  scale_fill_gradient(
    low = "#F5F9FB", high = "#9CC2D6",
    trans = "log10",
    labels = comma,
    name = "DEGs"
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = "padj threshold",
       title = "DEG count across (padj, |log2FC|) sweep") +
  theme_masld() +
  theme(
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    legend.key.height = unit(0.8, "cm"),
    legend.key.width = unit(0.25, "cm")
  )

# CV of N DEGs across padj thresholds, at each |log2FC| cutoff
cv_lfc <- sweep[, .(cv_n_degs = sd(n_degs) / mean(n_degs)), by = .(lfc_threshold, lfc_label)][order(lfc_threshold)]
cv_min_a     <- min(cv_lfc$cv_n_degs, na.rm = TRUE)
cv_min_lfc_a <- cv_lfc[which.min(cv_n_degs), lfc_threshold]

p_a_curve <- ggplot(cv_lfc, aes(x = lfc_label, y = 100 * cv_n_degs, group = 1)) +
  geom_vline(xintercept = tier1_x_a,
             linetype = "dashed", color = "#FFB300", linewidth = 0.6) +
  geom_hline(yintercept = 100 * cv_min_a,
             linetype = "dotted", color = "gray50", linewidth = 0.4) +
  geom_line(color = "#37474F", linewidth = 0.7) +
  geom_point(color = "#37474F", size = 1.8) +
  scale_y_continuous(labels = function(v) paste0(v, "%")) +
  scale_x_discrete(expand = c(0, 0.5)) +
  labs(x = expression("|log"[2]*"FC| cutoff"),
       y = "CV of N DEGs across padj") +
  theme_masld() +
  theme(panel.grid.minor = element_blank())

p_a <- p_a_hm / p_a_curve + plot_layout(heights = c(2.4, 1))

# ═══════════════════════════════════════════════════════════════════════════
# Panel (c): Marginal Trade-off at padj = 0.05
# ═══════════════════════════════════════════════════════════════════════════
message("Panel (c): Marginal trade-off at padj = 0.05...")

marginal <- sweep[padj_threshold == 0.05]

# Scale factor for dual y-axis
max_degs <- max(marginal$n_degs)
scale_factor <- max_degs / 100

p_c <- ggplot(marginal, aes(x = lfc_threshold)) +
  # DEG count (bars)
  geom_col(aes(y = n_degs), fill = masld_colors$down, alpha = 0.7, width = 0.08) +
  # Control recovery (line + points)
  geom_line(aes(y = pct_ctrl * scale_factor), color = masld_colors$up, linewidth = 0.7) +
  geom_point(aes(y = pct_ctrl * scale_factor), color = masld_colors$up, size = 1.5) +
  # Canonical cutoff reference line + label
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "#FFB300", linewidth = 0.4) +
  annotate("text", x = 0.52, y = max_degs * 0.95,
           label = sweep[padj_threshold == 0.05 & lfc_threshold == 0.5,
                         sprintf("canonical |LFC|=0.5\n%s DEGs · %d%% ctrl",
                                 format(n_degs, big.mark = ","),
                                 round(pct_ctrl))],
           size = 2, color = "#B26A00", hjust = 0, vjust = 1, lineheight = 0.8) +
  scale_y_continuous(
    name = "Number of DEGs",
    labels = comma,
    sec.axis = sec_axis(~ . / scale_factor,
                        name = "% positive controls recovered",
                        labels = function(x) paste0(x, "%"))
  ) +
  scale_x_continuous(breaks = lfc_vals) +
  labs(x = "|log2FC| threshold (at padj < 0.05)",
       title = "DEG count vs control-recovery trade-off") +
  theme_masld() +
  theme(
    axis.title.y.right = element_text(color = masld_colors$up, angle = 90),
    axis.text.y.right = element_text(color = masld_colors$up)
  )

# ═══════════════════════════════════════════════════════════════════════════
# Panel (d): padj vs ashr framework concordance
# ═══════════════════════════════════════════════════════════════════════════
message("Panel (d): Standard vs ashr framework concordance...")

compare <- rbindlist(list(sweep_ashr, sweep_padj_compare), use.names = TRUE)
compare[, lfc_label := factor(sprintf("|LFC| >= %.1f", lfc_threshold))]

p_d <- ggplot(compare[padj_threshold == 0.05],
              aes(x = lfc_threshold, y = n_degs / 1000,
                  color = framework, shape = framework)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1.5) +
  scale_color_manual(values = c(
    "Standard (padj + logFC)" = masld_colors$down,
    "ashr (lfsr + shrunk LFC)" = masld_colors$up
  )) +
  labs(x = "|log2FC| or |shrunk LFC| threshold",
       y = "DEGs (thousands)",
       color = NULL, shape = NULL,
       title = "Standard vs ashr framework (sig < 0.05)") +
  theme_masld() +
  theme(legend.position = c(0.7, 0.8),
        legend.background = element_rect(fill = "white", color = NA),
        legend.text = element_text(size = 5))

# ═══════════════════════════════════════════════════════════════════════════
# Panel (e): Jaccard similarity of DEG sets at |LFC| >= 0.5 across padj
# ═══════════════════════════════════════════════════════════════════════════
message("Panel (e): Jaccard stability matrix at |LFC| >= 0.5...")

# Build gene sets at each padj threshold with fixed |LFC| >= 0.5
deg_sets_lfc05 <- lapply(padj_vals, function(pa) {
  dream[padj < pa & abs(logFC) >= 0.5, symbol]
})
names(deg_sets_lfc05) <- ifelse(
  padj_vals < 0.001,
  format(padj_vals, scientific = TRUE),
  as.character(padj_vals)
)

jaccard <- function(a, b) length(intersect(a, b)) / length(union(a, b))

# Full Jaccard matrix
jmat <- outer(seq_along(deg_sets_lfc05), seq_along(deg_sets_lfc05),
              Vectorize(function(i, j) jaccard(deg_sets_lfc05[[i]], deg_sets_lfc05[[j]])))
rownames(jmat) <- colnames(jmat) <- names(deg_sets_lfc05)

set_sizes <- vapply(deg_sets_lfc05, length, integer(1))

# Melt — lower triangle + diagonal only (upper triangle set to NA / blank)
jdt <- as.data.table(as.table(jmat), keep.rownames = FALSE)
setnames(jdt, c("padj_x", "padj_y", "jaccard"))
lvls <- names(deg_sets_lfc05)
jdt[, padj_x := factor(padj_x, levels = lvls)]
jdt[, padj_y := factor(padj_y, levels = rev(lvls))]

# row index > col index = lower triangle (after accounting for reversed y)
jdt[, xi := as.integer(padj_x)]
jdt[, yi := length(lvls) + 1L - as.integer(padj_y)]  # un-reverse for index compare
jdt[xi < yi, jaccard := NA_real_]                     # blank upper triangle

jdt[, is_diag := xi == yi]
jdt[, label := fcase(
  is_diag,  format(set_sizes[as.character(padj_x)], big.mark = ","),
  !is_diag & !is.na(jaccard), sprintf("%.3f", jaccard),
  default = ""
)]

p_e <- ggplot(jdt, aes(x = padj_x, y = padj_y, fill = jaccard)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = label, fontface = ifelse(is_diag, "bold", "plain")),
            size = 2.2, color = "gray10") +
  scale_fill_gradient(low = "#F5F9FB", high = "#4A90A4",
                      limits = c(0, 1), na.value = "white",
                      name = "Jaccard") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = "padj threshold", y = "padj threshold",
       title = "DEG-set stability at |log₂FC| ≥ 0.5",
       subtitle = "Jaccard similarity; diagonal = N DEGs") +
  theme_masld() +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        legend.key.height = unit(0.8, "cm"),
        legend.key.width  = unit(0.25, "cm"))

# ═══════════════════════════════════════════════════════════════════════════
# Save individual panels
# ═══════════════════════════════════════════════════════════════════════════
PANEL_DIR <- file.path(OUT_DIR, "panels")
dir.create(PANEL_DIR, showWarnings = FALSE, recursive = TRUE)

save_fig(p_a, file.path(PANEL_DIR, "deg_landscape_a_count.pdf"),
         width = fig_half_width, height = 4.2)
save_fig(p_c, file.path(PANEL_DIR, "deg_landscape_c_marginal.pdf"),
         width = fig_half_width, height = 2.5)
save_fig(p_d, file.path(PANEL_DIR, "deg_landscape_d_concordance.pdf"),
         width = fig_half_width, height = 2.5)
save_fig(p_e, file.path(PANEL_DIR, "deg_landscape_e_stability.pdf"),
         width = fig_half_width, height = fig_half_width)

# ═══════════════════════════════════════════════════════════════════════════
# Assemble composite figure
# ═══════════════════════════════════════════════════════════════════════════
message("Assembling composite figure...")
fig <- patchwork::wrap_elements(p_a) / (p_c | p_d | p_e)
fig <- auto_tag(fig)

OUT_PDF <- file.path(OUT_DIR, "figS_deg_threshold_landscape.pdf")
save_fig(fig, OUT_PDF, width = fig_full_width, height = 7.5)
message("Saved: ", OUT_PDF)

# ═══════════════════════════════════════════════════════════════════════════
# Save companion data
# ═══════════════════════════════════════════════════════════════════════════
CSV_OUT <- file.path(OUT_DIR, "figS_deg_threshold_landscape_data.csv")
fwrite(sweep[, .(padj_threshold, lfc_threshold, n_degs, n_ctrl, pct_ctrl)], CSV_OUT)
message("Saved companion CSV: ", CSV_OUT, " (", nrow(sweep), " rows)")

# Print summary table for quick reference
message("\n=== Key Threshold Combinations (padj < 0.05) ===")
marginal_summary <- sweep[padj_threshold == 0.05,
                          .(lfc_threshold, n_degs, n_ctrl, pct_ctrl = round(pct_ctrl, 1))]
print(marginal_summary)

message("\n=== figS_deg_threshold_landscape.R complete ===")
