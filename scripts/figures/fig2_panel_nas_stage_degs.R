#!/usr/bin/env Rscript
# ============================================================================
# fig2_panel_nas_stage_degs.R
# Fig 2 panel — NAS stage-specific up/down DEG counts (vs NAS0).
# Restored from archive/legacy_pre_redesign_2026-05-17/fig2b.pdf (2026-05-19).
# Panel suffix dropped per convention; Illustrator handles labeling.
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
OUT_PDF <- file.path(PANEL_DIR, "fig2_panel_nas_stage_degs.pdf")

NAS_STAGE <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/disease_signatures/nas_score_dream.csv")

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

nas_de <- fread(NAS_STAGE)
nas_padj_col <- if ("padj" %in% names(nas_de)) "padj" else "adj.P.Val"
nas_counts <- nas_de[, count_ud(.SD, nas_padj_col), by = contrast]
nas_counts[, contrast := factor(contrast,
    levels = paste0("NAS", 1:7, "_vs_NAS0"),
    labels = paste0("NAS", 1:7))]
setorder(nas_counts, contrast)

long <- melt(nas_counts, id.vars = "contrast",
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
  labs(title = "NAS stage DEGs",
       x = "NAS stage (vs NAS0)", y = "DEGs") +
  theme_masld(base_size = BASE_SIZE) +
  theme(plot.title    = element_text(size = BASE_SIZE),
        axis.title    = element_text(size = BASE_SIZE),
        axis.text     = element_text(size = BASE_SIZE),
        legend.title  = element_text(size = BASE_SIZE),
        legend.text   = element_text(size = BASE_SIZE),
        plot.subtitle = element_blank(),
        plot.margin   = margin(3, 3, 3, 3),
        legend.position = "top",
        legend.justification = "left",
        axis.text.x = element_text(angle = 30, hjust = 1))

save_fig(panel, OUT_PDF,
         width = fig_half_width, height = 3.0)
fwrite(nas_counts, file.path(DATA_DIR, "fig2_panel_nas_stage_degs_data.csv"))
cat(sprintf("[saved] %s\n", OUT_PDF))
