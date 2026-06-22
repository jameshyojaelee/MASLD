#!/usr/bin/env Rscript
# ============================================================================
# zonation_directional_polarity_xspecies.R
# Fig 3 (RNA-seq) panel — cross-species zonation directional polarity
#
# Shared-x two-track composite. For the SAME zoned DEGs in the SAME gene order:
#   TRACK 1 (top)    = human bulk disease logFC  (diverging lollipop)
#   TRACK 2 (bottom) = MOUSE ortholog meta logFC for those same genes
# Vertical alignment exposes per-gene human<->mouse zonation conservation
# (same sign) vs divergence (opposite sign). The lipogenic pericentral exemplars
# FASN/ACLY go human-UP but mouse-DOWN (species divergence); the periportal
# urea-cycle genes ASS1/ARG1/CPS1 go DOWN in both (conserved).
#
# Companion to (does NOT replace) zonation_directional_polarity.R, which is the
# single-species human-only panel.
#
# Output: figures/main/fig3_RNAseq/panels/zonation_directional_polarity_xspecies.pdf
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

PANEL_DIR <- file.path(FIG4_DIR, "_supp")          # moved to main Fig 4 (spatial/zonation home) 2026-06-18
DATA_DIR  <- file.path(PANEL_DIR, "data")
OUT_PDF   <- file.path(PANEL_DIR, "zonation_directional_polarity_xspecies.pdf")
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)

UP_COL   <- masld_colors$mash   # #C9265E Liang magenta (up)
DOWN_COL <- masld_colors$down   # #1565C0 Liang blue   (down)

# ---------------------------------------------------------------------------
# Load + filter to zoned DEGs that carry a mouse ortholog meta logFC
# ---------------------------------------------------------------------------
z <- fread(file.path(BASE, "RNA-seq/results/zonation/deg_zonation_classification.csv"))
z <- z[zonation_class %in% c("Pericentral", "Periportal") &
         (is_deg_up | is_deg_down) & !is.na(mouse_meta_logFC)]
z[, zone := factor(zonation_class, levels = c("Pericentral", "Periportal"))]

# Per-gene cross-species concordance (same sign of human vs mouse logFC)
z[, xspecies := ifelse(sign(bulk_logFC) == sign(mouse_meta_logFC),
                       "Conserved", "Divergent")]

# Shared gene order: within each zone, order by human bulk_logFC (descending).
# A single global factor used as the discrete x on BOTH tracks guarantees that
# column i is the SAME gene in the human and mouse tracks.
setorder(z, zone, -bulk_logFC)
z[, gene := factor(human_symbol, levels = human_symbol)]

# Long form: one row per (gene, species) for shared-x plotting downstream
hum <- z[, .(gene, zone, human_symbol, xspecies, is_conserved,
             value = bulk_logFC,
             direction = ifelse(bulk_logFC > 0, "Up", "Down"))]
mou <- z[, .(gene, zone, human_symbol, xspecies, is_conserved,
             value = mouse_meta_logFC,
             direction = ifelse(mouse_meta_logFC > 0, "Up", "Down"))]

# ---------------------------------------------------------------------------
# Genes to label (divergence + conserved exemplars)
# ---------------------------------------------------------------------------
lab_genes <- c("FASN", "ACLY",            # pericentral: human-up / mouse-down (divergent)
               "SCD", "ACACA",            # other lipogenic divergence if present
               "ASS1", "ARG1", "CPS1")    # periportal urea cycle: down in both (conserved)
hum[, lab := ifelse(human_symbol %in% lab_genes, human_symbol, NA_character_)]

# ---------------------------------------------------------------------------
# Hero numbers
# ---------------------------------------------------------------------------
n_total <- nrow(z)
n_cons  <- z[xspecies == "Conserved", .N]
n_div   <- z[xspecies == "Divergent", .N]
cat(sprintf("[hero] %d zoned DEGs with mouse ortholog logFC: %d conserved / %d divergent (by sign)\n",
            n_total, n_cons, n_div))
cat(sprintf("[hero] Pericentral n=%d  Periportal n=%d\n",
            z[zone == "Pericentral", .N], z[zone == "Periportal", .N]))

# Data sidecar
fwrite(z[, .(human_symbol, zonation_class, bulk_logFC, mouse_meta_logFC,
             xspecies, is_conserved,
             human_dir = ifelse(bulk_logFC > 0, "Up", "Down"),
             mouse_dir = ifelse(mouse_meta_logFC > 0, "Up", "Down"))],
       file.path(DATA_DIR, "zonation_directional_polarity_xspecies.csv"))

# ---------------------------------------------------------------------------
# Shared discrete-x geometry so the two tracks align column-for-column.
# Identical x expansion on both panels keeps gene positions equal.
# ---------------------------------------------------------------------------
X_EXPAND <- ggplot2::expansion(add = 0.6)

# Mark cross-species-divergent genes on the x axis (bold colour) so the eye can
# anchor the vertical sign-flips.
div_genes <- z[xspecies == "Divergent", as.character(gene)]
axis_cols <- ifelse(levels(z$gene) %in% div_genes, UP_COL, "gray30")
axis_face <- ifelse(levels(z$gene) %in% div_genes, "bold.italic", "italic")

# Common y range (symmetric) so up/down magnitudes read on the same visual scale
ymax <- max(abs(c(hum$value, mou$value))) * 1.10

base_track <- function(df, ylab, show_x) {
  p <- ggplot(df, aes(x = gene, y = value, color = direction)) +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "#9E9E9E") +
    geom_segment(aes(xend = gene, y = 0, yend = value), linewidth = 0.45) +
    geom_point(size = 1.4) +
    facet_grid(. ~ zone, scales = "free_x", space = "free_x") +
    scale_color_manual(values = c(Up = UP_COL, Down = DOWN_COL), name = NULL) +
    scale_x_discrete(expand = X_EXPAND) +
    scale_y_continuous(limits = c(-ymax, ymax),
                       expand = expansion(mult = c(0.02, 0.02))) +
    labs(x = NULL, y = ylab) +
    theme_masld(base_size = 7)
  if (show_x) {
    p <- p + theme(
      axis.text.x  = element_text(size = 5.2, angle = 90, hjust = 1, vjust = 0.5,
                                  face = axis_face, color = axis_cols),
      strip.text.x = element_blank(),          # zone strip only on top track
      panel.grid.major.x = element_blank(),
      legend.position = "none")
  } else {
    p <- p + theme(
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x  = element_blank(),
      strip.text.x = element_text(size = 6.8, face = "bold"),
      panel.grid.major.x = element_blank(),
      legend.position = "none")
  }
  p
}

# TRACK 1 — human (top), blanks the gene axis, carries the zone strip + labels
p_top <- base_track(hum, "Human\ndisease logFC", show_x = FALSE) +
  geom_text_repel(aes(label = lab), size = 2.0, fontface = "italic",
                  segment.size = 0.2, min.segment.length = 0,
                  max.overlaps = Inf, box.padding = 0.3,
                  na.rm = TRUE, show.legend = FALSE) +
  labs(title = "Human↔mouse zonated-DEG directional polarity") +
  theme(plot.title = element_text(size = 7.6, face = "bold", margin = margin(b = 4)))

# TRACK 2 — mouse (bottom), carries the gene axis
p_bottom <- base_track(mou, "Mouse ortholog\nmeta logFC", show_x = TRUE)

# ---------------------------------------------------------------------------
# Assemble: vertical stack, shared discrete x; only bottom track shows genes.
# ---------------------------------------------------------------------------
panel <- p_top / p_bottom +
  plot_layout(heights = c(1, 1.45))   # extra room for the rotated gene labels

ggsave(OUT_PDF, panel,
       width  = 120 / 25.4,
       height = 110 / 25.4,
       units  = "in",
       device = cairo_pdf)
cat(sprintf("[saved] %s\n", OUT_PDF))
