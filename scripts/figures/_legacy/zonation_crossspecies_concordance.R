#!/usr/bin/env Rscript
# ============================================================================
# zonation_crossspecies_concordance.R
# Fig 3 SUPP (figS06_cross_species) — human vs mouse logFC of zonated DEGs
#
# Periportal nitrogen-handling/complement shutdown is CONSERVED human->mouse,
# but the induced pericentral lipogenic program (FASN/ACLY/ACACA) is
# paradoxically SUPPRESSED in mouse — a species-specific divergence of de
# novo lipogenesis zonation.
#
# Output: figures/supplementary/figS06_cross_species/
#         zonation_crossspecies_concordance.pdf  (90 x 90 mm)
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

OUT_DIR  <- file.path(FIG4_DIR, "_supp")  # demoted to Fig 4 supplement (_supp) 2026-06-18 (panel reorg)
DATA_DIR <- file.path(OUT_DIR, "data")
OUT_PDF  <- file.path(OUT_DIR, "zonation_crossspecies_concordance.pdf")
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)

PC_COL <- "#C9265E"  # Pericentral
PP_COL <- "#1565C0"  # Periportal

# ---------------------------------------------------------------------------
# Load + filter: zonated DEGs with a non-NA mouse meta logFC
# ---------------------------------------------------------------------------
z <- fread(file.path(BASE, "RNA-seq/results/zonation/deg_zonation_classification.csv"))
z <- z[zonation_class %in% c("Pericentral", "Periportal") &
         is_deg == TRUE & !is.na(mouse_meta_logFC)]
z[, zone := factor(zonation_class, levels = c("Pericentral", "Periportal"))]
z[, concordant := sign(bulk_logFC) == sign(mouse_meta_logFC)]

# ---------------------------------------------------------------------------
# Hero numbers
# ---------------------------------------------------------------------------
n_tot  <- nrow(z)
n_cons <- z[is_conserved == TRUE, .N]
n_disc <- z[is_conserved == FALSE, .N]
cat(sprintf("[hero] %d zonated DEGs plotted; %d conserved / %d discordant\n",
            n_tot, n_cons, n_disc))

fwrite(z[, .(human_symbol, zonation_class, bulk_logFC, mouse_meta_logFC,
             is_conserved, concordant)],
       file.path(DATA_DIR, "zonation_crossspecies_concordance.csv"))

# ---------------------------------------------------------------------------
# Genes to label: conserved periportal-down + discordant pericentral-lipogenic
# ---------------------------------------------------------------------------
lab_genes <- c("ASS1", "ARG1", "CPS1",           # conserved periportal nitrogen
               "FASN", "ACLY", "ACACA", "SCD")   # discordant pericentral DNL
z[, lab := ifelse(human_symbol %in% lab_genes, human_symbol, NA_character_)]

lims <- range(c(z$bulk_logFC, z$mouse_meta_logFC), na.rm = TRUE)

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
p <- ggplot(z, aes(x = bulk_logFC, y = mouse_meta_logFC, color = zone)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "#9E9E9E") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "#9E9E9E") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              linewidth = 0.3, color = "#BDBDBD") +
  geom_point(size = 1.6, alpha = 0.9) +
  geom_text_repel(aes(label = lab), size = 2.0, fontface = "italic",
                  segment.size = 0.2, min.segment.length = 0,
                  max.overlaps = Inf, box.padding = 0.3,
                  na.rm = TRUE, show.legend = FALSE) +
  scale_color_manual(values = c(Pericentral = PC_COL, Periportal = PP_COL),
                     name = NULL) +
  coord_equal(xlim = lims, ylim = lims) +
  labs(x = "Human bulk logFC", y = "Mouse meta logFC",
       title = "Cross-species zonation concordance") +
  theme_masld(base_size = 7) +
  theme(
    plot.title      = element_text(size = 7.3, face = "bold", margin = margin(b = 6)),
    legend.position = "top"
  )

ggsave(OUT_PDF, p,
       width  = 90 / 25.4,
       height = 90 / 25.4,
       units  = "in",
       device = cairo_pdf)
cat(sprintf("[saved] %s\n", OUT_PDF))
