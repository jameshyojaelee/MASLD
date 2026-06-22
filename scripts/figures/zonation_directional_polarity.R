#!/usr/bin/env Rscript
# ============================================================================
# zonation_directional_polarity.R
# Fig 3 (RNA-seq) panel — directional polarity of zonated DEGs in MASLD
#
# Diverging lollipop, two grouped blocks (Pericentral vs Periportal). MASLD
# selectively REWIRES the pericentral zone (lipogenic/Wnt induced FASN/SCD/
# ACLY/AKR1B10/AXIN2 up; detox CYP1A2/2E1/3A4 down — bidirectional) while
# MONOLITHICALLY shutting down periportal identity (26/27 zoned periportal
# DEGs down).
#
# Output: figures/main/fig3_RNAseq/panels/zonation_directional_polarity.pdf
#   sized 90 x 120 mm
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

PANEL_DIR <- file.path(FIG4_DIR, "_supp")  # demoted to Fig 4 supplement (_supp) 2026-06-18 (panel reorg)
DATA_DIR  <- file.path(PANEL_DIR, "data")
OUT_PDF   <- file.path(PANEL_DIR, "zonation_directional_polarity.pdf")
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)

UP_COL   <- "#C9265E"  # Liang magenta
DOWN_COL <- "#1565C0"  # Liang blue

# ---------------------------------------------------------------------------
# Load + filter to zonated DEGs
# ---------------------------------------------------------------------------
z <- fread(file.path(BASE, "RNA-seq/results/zonation/deg_zonation_classification.csv"))
z <- z[zonation_class %in% c("Pericentral", "Periportal") & (is_deg_up | is_deg_down)]
z[, direction := ifelse(bulk_logFC > 0, "Up", "Down")]
z[, zone := factor(zonation_class, levels = c("Pericentral", "Periportal"))]

# Order within each block by logFC (descending) for a clean diverging fan
setorder(z, zone, -bulk_logFC)
z[, y := .I]                                  # global row index, top = highest pericentral up
z[, y := rev(y)]                              # flip so highest logFC sits at top

# ---------------------------------------------------------------------------
# Hero numbers
# ---------------------------------------------------------------------------
pc_up <- z[zone == "Pericentral" & direction == "Up", .N]
pc_dn <- z[zone == "Pericentral" & direction == "Down", .N]
pp_up <- z[zone == "Periportal"  & direction == "Up", .N]
pp_dn <- z[zone == "Periportal"  & direction == "Down", .N]
cat(sprintf("[hero] Pericentral: %d up / %d down (bidirectional)\n", pc_up, pc_dn))
cat(sprintf("[hero] Periportal:  %d up / %d down (near-uniform negative)\n", pp_up, pp_dn))

fwrite(z[, .(human_symbol, zonation_class, bulk_logFC, bulk_padj, direction, y)],
       file.path(DATA_DIR, "zonation_directional_polarity.csv"))

# ---------------------------------------------------------------------------
# Genes to label
# ---------------------------------------------------------------------------
lab_genes <- c("FASN", "SCD", "ACLY", "AKR1B10", "AXIN2",      # pericentral up
               "CYP2E1", "CYP1A2", "CYP3A4",                   # pericentral down
               "ASS1", "ARG1", "CPS1")                         # periportal down
z[, lab := ifelse(human_symbol %in% lab_genes, human_symbol, NA_character_)]

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
p <- ggplot(z, aes(x = bulk_logFC, y = y, color = direction)) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "#9E9E9E") +
  geom_segment(aes(x = 0, xend = bulk_logFC, y = y, yend = y), linewidth = 0.45) +
  geom_point(size = 1.4) +
  geom_text_repel(aes(label = lab), size = 2.0, fontface = "italic",
                  segment.size = 0.2, min.segment.length = 0,
                  max.overlaps = Inf, box.padding = 0.25,
                  na.rm = TRUE, show.legend = FALSE) +
  facet_grid(zone ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_color_manual(values = c(Up = UP_COL, Down = DOWN_COL), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  labs(x = "Bulk disease logFC", y = NULL,
       title = "Zonated DEG directional polarity") +
  theme_masld(base_size = 7) +
  theme(
    plot.title       = element_text(size = 7.3, face = "bold", margin = margin(b = 6)),
    axis.text.y      = element_blank(),
    axis.ticks.y     = element_blank(),
    panel.grid.major.y = element_blank(),
    strip.placement  = "outside",
    strip.text.y.left = element_text(size = 6.5, face = "bold", angle = 90),
    legend.position  = "top"
  )

ggsave(OUT_PDF, p,
       width  = 90 / 25.4,
       height = 120 / 25.4,
       units  = "in",
       device = cairo_pdf)
cat(sprintf("[saved] %s\n", OUT_PDF))
