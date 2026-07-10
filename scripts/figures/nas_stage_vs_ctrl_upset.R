#!/usr/bin/env Rscript
# ============================================================================
# nas_stage_vs_ctrl_upset.R
# Fig 2 panel — UpSet of NAS stage-specific DEG sets (NAS1-7 vs strictly
# healthy controls). Mirrors nas_stage_upset.R.
# DEG definition: padj<0.05, |logFC|>0.5. Top 11 intersections shown.
# ============================================================================

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

PANEL_DIR <- file.path(FIG2_DIR, "panels")
DATA_DIR  <- file.path(PANEL_DIR, "data")
dir.create(PANEL_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(DATA_DIR,  showWarnings = FALSE, recursive = TRUE)
OUT_PDF <- file.path(PANEL_DIR, "nas_stage_vs_healthy_upset.pdf")
OUT_CSV <- file.path(DATA_DIR,  "nas_stage_vs_healthy_upset_intersections.csv")

NAS_VS_CTRL <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/disease_signatures/nas_stage_vs_ctrl_dream.csv")

if (!file.exists(NAS_VS_CTRL)) {
  stop("nas_stage_vs_ctrl_dream.csv not found. Run 14b_stage_vs_healthy_dream.R first.")
}

PADJ_THR <- 0.05
LFC_THR  <- 0.5
TOP_N    <- 15L

BAR_FILL <- "#546E7A"

# ---------------------------------------------------------------------------
# Build per-stage DEG sets
# ---------------------------------------------------------------------------
de <- fread(NAS_VS_CTRL)
padj_col <- if ("padj" %in% names(de)) "padj" else "adj.P.Val"
de <- de[!is.na(get(padj_col)) & !is.na(logFC)]
de_sig <- de[get(padj_col) < PADJ_THR & abs(logFC) > LFC_THR, .(gene, contrast)]

stage_labels <- paste0("NAS", 1:7)
contrast_map <- setNames(stage_labels, paste0("NAS", 1:7, "_vs_Ctrl"))
de_sig[, stage := contrast_map[contrast]]
de_sig <- de_sig[!is.na(stage)]
# Drop any stages with 0 DEGs (they won't form meaningful intersections)
present_stages <- de_sig[, unique(stage)]
stage_labels   <- stage_labels[stage_labels %in% present_stages]

if (length(stage_labels) == 0) {
  stop("No significant DEGs found for any NAS stage vs control. Check analysis results.")
}

gene_sets <- dcast(unique(de_sig[, .(gene, stage)]),
                   gene ~ stage, fun.aggregate = length, fill = 0L)
for (s in stage_labels) {
  if (!s %in% names(gene_sets)) gene_sets[, (s) := 0L]
  gene_sets[, (paste0("in_", s)) := get(s) > 0L]
}
membership_cols <- paste0("in_", stage_labels)

gene_sets[, intersection_id := apply(.SD, 1,
  function(r) paste0(as.integer(r), collapse = "")),
  .SDcols = membership_cols]

intersections <- gene_sets[, .(n_genes = .N), by = intersection_id]
intersections[, c(membership_cols) := lapply(seq_along(membership_cols), function(i) {
  substr(intersection_id, i, i) == "1"
})]
intersections[, n_sets := rowSums(.SD), .SDcols = membership_cols]
intersections <- intersections[n_sets > 0]
setorder(intersections, -n_genes)
intersections <- head(intersections, TOP_N)
intersections[, ix_rank  := .I]
intersections[, ix_label := factor(ix_rank, levels = ix_rank)]

set_sizes <- data.table(
  stage   = stage_labels,
  n_genes = sapply(stage_labels, function(s) sum(gene_sets[[paste0("in_", s)]]))
)
set_sizes[, stage := factor(stage, levels = rev(stage_labels))]

# NAS color ramp: light orange → orange → purple (matches nas_stage_upset.R)
n_nas  <- length(stage_labels)
nas_palette <- colorRampPalette(c("#FFE0B2", "#F57F17", "#7B1FA2"))(n_nas)
nas_colors  <- setNames(nas_palette, stage_labels)

# ---------------------------------------------------------------------------
# Top: intersection-size bar
# ---------------------------------------------------------------------------
p_top <- ggplot(intersections, aes(x = ix_label, y = n_genes)) +
  geom_col(width = 0.75, fill = BAR_FILL, color = NA) +
  geom_text(aes(label = comma(n_genes)), vjust = -0.3, size = 1.9, color = "gray20") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18)), labels = comma) +
  labs(x = NULL, y = "DEGs (intersection size)") +
  theme_masld() +
  theme(axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin  = margin(2, 4, 0, 4))

# ---------------------------------------------------------------------------
# Bottom right: dot matrix
# ---------------------------------------------------------------------------
dot_dt <- melt(intersections[, c("ix_label", membership_cols), with = FALSE],
               id.vars = "ix_label",
               variable.name = "stage", value.name = "in_set")
dot_dt[, stage := factor(sub("^in_", "", stage), levels = rev(stage_labels))]

segs <- dot_dt[in_set == TRUE,
               .(ymin = min(as.integer(stage)),
                 ymax = max(as.integer(stage))),
               by = ix_label]
segs <- segs[ymin != ymax]

p_dots <- ggplot(dot_dt, aes(x = ix_label, y = stage)) +
  geom_point(aes(color = in_set), size = 2.0) +
  { if (nrow(segs)) geom_segment(data = segs,
      aes(x = ix_label, xend = ix_label, y = ymin, yend = ymax),
      inherit.aes = FALSE, color = "#212121", linewidth = 0.5) } +
  scale_color_manual(values = c("TRUE" = "#212121", "FALSE" = "#E0E0E0"),
                     guide = "none") +
  labs(x = NULL, y = NULL) +
  theme_masld() +
  theme(panel.grid    = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y  = element_text(face = "plain", size = 6.5),
        plot.margin  = margin(0, 4, 4, 4))

# ---------------------------------------------------------------------------
# Bottom left: set-size bar
# ---------------------------------------------------------------------------
p_left <- ggplot(set_sizes,
                 aes(y = stage, x = n_genes, fill = as.character(stage))) +
  geom_col(width = 0.7, color = NA) +
  geom_text(aes(label = comma(n_genes)), hjust = 1.1, size = 1.9, color = "white") +
  scale_fill_manual(values = nas_colors, guide = "none") +
  scale_x_reverse(expand = expansion(mult = c(0.05, 0)), labels = comma,
                  breaks = scales::breaks_pretty(n = 3)) +
  labs(x = "Set size", y = NULL) +
  theme_masld() +
  theme(panel.grid    = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x  = element_text(size = 5.5, angle = 35, hjust = 1),
        plot.margin  = margin(0, 1, 4, 4))

# ---------------------------------------------------------------------------
# Assemble
# ---------------------------------------------------------------------------
aligned <- align_plots(p_top, p_dots, align = "v", axis = "lr")

top_row    <- plot_grid(NULL, aligned[[1]], rel_widths = c(1.3, 5), nrow = 1)
bottom_row <- plot_grid(p_left, aligned[[2]], rel_widths = c(1.3, 5),
                        nrow = 1, align = "h", axis = "tb")
upset      <- plot_grid(top_row, bottom_row, ncol = 1, rel_heights = c(1.6, 1.0))

ggsave(OUT_PDF, upset,
       width = fig_full_width * 0.75, height = 3.4,
       device = cairo_pdf)

# Save intersection table with dynamic column names
out_cols <- c("intersection_id", membership_cols, "n_sets", "n_genes")
fwrite(intersections[, out_cols, with = FALSE], OUT_CSV)
cat(sprintf("[saved] %s\n", OUT_PDF))
cat(sprintf("[saved] %s\n", OUT_CSV))
