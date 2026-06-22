#!/usr/bin/env Rscript
# ============================================================================
# fib_stage_vs_ctrl_degs.R
# Fig 2 panel — Fibrosis stage-specific up/down DEG counts vs STRICTLY
# HEALTHY controls (condition=="Control"; n=64), not vs all F0 samples.
# Mirrors fib_stage_degs.R but uses fibrosis_stage_vs_ctrl_dream.csv.
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

PANEL_DIR <- file.path(FIG2_DIR, "panels")
DATA_DIR  <- file.path(PANEL_DIR, "data")
dir.create(PANEL_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(DATA_DIR,  showWarnings = FALSE, recursive = TRUE)
OUT_PDF <- file.path(PANEL_DIR, "figs3d_fib_stage_vs_healthy_degs.pdf")

FIB_VS_CTRL <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/disease_signatures/fibrosis_stage_vs_ctrl_dream.csv")
SIZES_FILE  <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/disease_signatures/fibrosis_stage_vs_ctrl_sample_sizes.csv")

if (!file.exists(FIB_VS_CTRL)) {
  stop("fibrosis_stage_vs_ctrl_dream.csv not found. Run 14b_stage_vs_healthy_dream.R first.")
}

BASE_SIZE <- 7
LBL_SIZE  <- 7 / ggplot2::.pt
col_up    <- masld_colors$up
col_down  <- masld_colors$down

count_ud <- function(dt, padj_col, lfc_col = "logFC",
                     padj_thr = 0.05, lfc_thr = 0.5) {
  up   <- sum(dt[[padj_col]] < padj_thr & dt[[lfc_col]] >  lfc_thr, na.rm = TRUE)
  down <- sum(dt[[padj_col]] < padj_thr & dt[[lfc_col]] < -lfc_thr, na.rm = TRUE)
  data.table(up = up, down = down, total = up + down)
}

fib_de <- fread(FIB_VS_CTRL)
# contrast column has values like "F1_vs_Ctrl", "F2_vs_Ctrl", etc.
padj_col <- if ("padj" %in% names(fib_de)) "padj" else "adj.P.Val"
fib_counts <- fib_de[, count_ud(.SD, padj_col), by = contrast]

# Load sample sizes for subtitle annotation
sizes <- if (file.exists(SIZES_FILE)) fread(SIZES_FILE) else NULL
n_ctrl <- if (!is.null(sizes)) sizes[grp == "Ctrl", N] else "?"
n_label <- sprintf("Reference: %s healthy controls (condition='Control')", n_ctrl)

# Order stages F1-F4
stage_order <- paste0("F", 1:4, "_vs_Ctrl")
stage_labels <- paste0("F", 1:4)
fib_counts <- fib_counts[contrast %in% stage_order]
fib_counts[, contrast := factor(contrast, levels = stage_order, labels = stage_labels)]
setorder(fib_counts, contrast)

long <- melt(fib_counts, id.vars = "contrast",
             measure.vars = c("up", "down"),
             variable.name = "direction", value.name = "n")
long[, direction := factor(direction, levels = c("up", "down"),
                            labels = c("Up", "Down"))]
long[, signed_n  := ifelse(direction == "Up", n, -n)]
long[, lbl_vjust := ifelse(direction == "Up", -0.35, 1.15)]

panel <- ggplot(long, aes(x = contrast, y = signed_n, fill = direction)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey25") +
  geom_text(aes(label = comma(n), vjust = lbl_vjust),
            size = LBL_SIZE) +
  scale_fill_manual(values = c(Up = col_up, Down = col_down), name = NULL) +
  scale_y_continuous(labels = function(x) comma(abs(x)),
                     expand = expansion(mult = c(0.18, 0.18))) +
  labs(title = "Fibrosis stage DEGs vs healthy controls",
       subtitle = n_label,
       x = "Fibrosis stage (vs healthy controls)", y = "DEGs") +
  theme_masld(base_size = BASE_SIZE) +
  theme(plot.title    = element_text(size = BASE_SIZE),
        plot.subtitle = element_text(size = BASE_SIZE - 1, colour = "grey40"),
        axis.title    = element_text(size = BASE_SIZE),
        axis.text     = element_text(size = BASE_SIZE),
        legend.title  = element_text(size = BASE_SIZE),
        legend.text   = element_text(size = BASE_SIZE),
        plot.margin   = margin(3, 3, 3, 3),
        legend.position = "top",
        legend.justification = "left",
        axis.text.x = element_text(angle = 0, hjust = 0.5))

save_fig(panel, OUT_PDF,
         width = fig_half_width, height = 3.0)
fwrite(fib_counts, file.path(DATA_DIR, "fib_stage_vs_ctrl_degs_data.csv"))
cat(sprintf("[saved] %s\n", OUT_PDF))
