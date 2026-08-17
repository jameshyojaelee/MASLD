#!/usr/bin/env Rscript
# ============================================================================
# fib_stage_upset.R
# Fig 2 panel — UpSet of stage-specific DEG sets across F1/F2/F3/F4 (vs F0).
# DEG definition: analytical TREAT, fdr_treat<0.05 at lfc=0.25 (canonical 2026-06-29). Top 11 intersections shown.
# Alignment via cowplot::align_plots; single neutral bar color.
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
OUT_PDF <- file.path(PANEL_DIR, "fib_stage_upset.pdf")
OUT_CSV <- file.path(DATA_DIR,  "fib_stage_upset_intersections.csv")

FIB_STAGE <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/disease_signatures/fibrosis_stage_dream.csv")

# CANONICAL 2026-08-12: padj < 0.05 AND |log2FC| > 0.5, BH within each stage contrast.
# Supersedes the per-contrast analytical TREAT at lfc=0.25 (canonical 2026-06-29..2026-08-12).
CANON_LFC <- 0.50
FDR_THR   <- 0.05
TREAT_LFC <- 0.25  # retained: sensitivity arm only, via add_treat_fdr() below
TOP_N     <- 11L

BAR_FILL <- "#546E7A"   # single neutral color for intersection bars

# Analytical TREAT reconstruction (mirrors rebuild_cas13_library.R::add_treat_fdr,
# proven identical to limma::treat). The per-stage source CSV carries logFC + SE + t +
# P.Value but no treat_fdr, so we reconstruct the moderated-t TREAT statistic here:
# H0 is |true logFC| <= TREAT_LFC, the effect-size floor folded INTO the test (no
# separate |logFC| filter). df.total inferred per contrast from the (t, P.Value) pair.
infer_df_total <- function(dt) {
  pr <- dt[is.finite(t) & is.finite(P.Value) & P.Value > 0 & P.Value < 1 & abs(t) > 1e-6]
  idx <- unique(round(seq(1, nrow(pr), length.out = min(nrow(pr), 12))))
  median(vapply(idx, function(i)
    uniroot(function(df) 2 * pt(-abs(pr$t[i]), df = df) - pr$P.Value[i], c(0.1, 1e6))$root,
    numeric(1)))
}
add_treat_fdr <- function(dt, lfc, se_col = "SE") {
  out <- copy(dt)
  se <- if (se_col %in% names(out)) out[[se_col]] else abs(out$logFC / out$t)
  se[!is.finite(se) | se <= 0] <- NA_real_
  df_use <- infer_df_total(out)
  out[, p_treat := pt((abs(logFC) - lfc) / se, df = df_use, lower.tail = FALSE) +
                   pt((abs(logFC) + lfc) / se, df = df_use, lower.tail = FALSE)]
  out[, fdr_treat := p.adjust(p_treat, method = "BH")]
  out
}

# ---------------------------------------------------------------------------
# Build per-stage DEG sets
# ---------------------------------------------------------------------------
de <- fread(FIB_STAGE)
stopifnot(all(c("logFC", "SE", "t", "P.Value", "contrast") %in% names(de)))
de <- de[!is.na(logFC) & !is.na(SE)]
# Per-contrast analytical TREAT (each stage-vs-F0 contrast is its own test family)
# CANONICAL 2026-08-12: padj<0.05 & |log2FC|>0.5, per contrast family.
de[, padj := p.adjust(P.Value, method = "BH"), by = contrast]
de_sig <- de[padj < FDR_THR & abs(logFC) > CANON_LFC, .(gene, contrast)]

stage_labels <- paste0("F", 1:4)
contrast_to_stage <- setNames(stage_labels, paste0("F", 1:4, "_vs_F0"))
de_sig[, stage := contrast_to_stage[contrast]]
de_sig <- de_sig[!is.na(stage)]

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
intersections <- intersections[intersection_id != "0101"] # Remove F2-F4 only overlap (24 DEGs)
setorder(intersections, -n_genes)
intersections <- head(intersections, TOP_N)
intersections[, ix_rank := .I]
intersections[, ix_label := factor(ix_rank, levels = ix_rank)]

set_sizes <- data.table(
  stage   = stage_labels,
  n_genes = sapply(stage_labels, function(s) sum(gene_sets[[paste0("in_", s)]]))
)
set_sizes[, stage := factor(stage, levels = rev(stage_labels))]

stage_colors <- c(
  "F1" = unname(fibrosis_stage_colors["F1"]),
  "F2" = unname(fibrosis_stage_colors["F2"]),
  "F3" = unname(fibrosis_stage_colors["F3"]),
  "F4" = unname(fibrosis_stage_colors["F4"])
)

# ---------------------------------------------------------------------------
# Top: intersection-size bar (single neutral color)
# ---------------------------------------------------------------------------
p_top <- ggplot(intersections, aes(x = ix_label, y = n_genes)) +
  geom_col(width = 0.75, fill = BAR_FILL, color = NA) +
  geom_text(aes(label = comma(n_genes)), vjust = -0.3, size = GEOM_TEXT_6PT, color = "gray20") +
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
        axis.text.y  = element_text(face = "plain", size = 6),
        plot.margin  = margin(0, 4, 4, 4))

# ---------------------------------------------------------------------------
# Bottom left: set-size bar
# ---------------------------------------------------------------------------
p_left <- ggplot(set_sizes,
                 aes(y = stage, x = n_genes, fill = as.character(stage))) +
  geom_col(width = 0.7, color = NA) +
  geom_text(aes(label = comma(n_genes)), hjust = 1.1, size = GEOM_TEXT_6PT, color = "white") +
  scale_fill_manual(values = stage_colors, guide = "none") +
  scale_x_reverse(expand = expansion(mult = c(0.05, 0)), labels = comma,
                  breaks = scales::breaks_pretty(n = 3)) +
  labs(x = "Set size", y = NULL) +
  theme_masld() +
  theme(panel.grid    = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x  = element_text(size = 6, angle = 35, hjust = 1),
        plot.margin  = margin(0, 1, 4, 4))

# ---------------------------------------------------------------------------
# Align top bar and dot matrix vertically, then assemble with cowplot
# ---------------------------------------------------------------------------
aligned <- align_plots(p_top, p_dots, align = "v", axis = "lr")

top_row    <- plot_grid(NULL, aligned[[1]], rel_widths = c(1.3, 5), nrow = 1)
bottom_row <- plot_grid(p_left, aligned[[2]], rel_widths = c(1.3, 5), nrow = 1, align = "h", axis = "tb")
upset      <- plot_grid(top_row, bottom_row, ncol = 1, rel_heights = c(1.6, 1.0))


ggsave(OUT_PDF, upset,
       width = fig_full_width * 0.55, height = 3.2,
       device = cairo_pdf)
fwrite(intersections[, .(intersection_id, in_F1, in_F2, in_F3, in_F4,
                          n_sets, n_genes)], OUT_CSV)
cat(sprintf("[saved] %s\n", OUT_PDF))
cat(sprintf("[saved] %s\n", OUT_CSV))
