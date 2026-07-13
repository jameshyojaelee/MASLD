#!/usr/bin/env Rscript
# ============================================================================
# cross_species_pathway_translatability_matrix.R
# Fig 3 SUPP (figS06_cross_species) — signed-NES pathway concordance heatmap
#
# Mouse diet models recapitulate human steatohepatitis/fibrosis pathway
# programs (signed NES rho 0.24-0.69) but INVERT the pure-steatosis program
# (nafl_specific rho negative in every diet).
#
# NOTE: overlaps the dedicated cross-species figure; placed here as a fig3
# SUPP panel.
#
# Output: figures/supplementary/figS06_cross_species/
#         cross_species_pathway_translatability_matrix.pdf  (90 x 70 mm)
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUT_DIR  <- FIGS06_DIR
DATA_DIR <- file.path(OUT_DIR, "data")
OUT_PDF  <- file.path(OUT_DIR, "cross_species_pathway_translatability_matrix.pdf")
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------
d <- fread(file.path(BASE,
  "Analysis/Cross_Species_Concordance/results/fgsea_pathway_concordance.csv"))

# Human signature ordering (readable: disease -> NASH -> fibrosis -> NAFL-specific)
sig_levels  <- c("disease_vs_ctrl", "nafl_vs_nash", "fibrosis", "nafl_specific")
sig_labels  <- c("Disease vs ctrl", "NAFL vs NASH", "Fibrosis", "NAFL-specific")
sig_levels  <- sig_levels[sig_levels %in% d$human_signature]
diet_levels <- c("HFD", "CDAHFD", "MCD", "FPC")
diet_levels <- diet_levels[diet_levels %in% d$diet]

d[, human_signature := factor(human_signature, levels = rev(sig_levels),
                              labels = rev(sig_labels[match(sig_levels, c("disease_vs_ctrl","nafl_vs_nash","fibrosis","nafl_specific"))]))]
d[, diet := factor(diet, levels = diet_levels)]

# ---------------------------------------------------------------------------
# Hero numbers
# ---------------------------------------------------------------------------
rng <- range(d$rho_NES, na.rm = TRUE)
neg_sig <- d[, .(all_neg = all(rho_NES < 0)), by = human_signature][all_neg == TRUE]
cat(sprintf("[hero] rho_NES range: %.3f to %.3f\n", rng[1], rng[2]))
cat(sprintf("[hero] negative across all diets: %s\n",
            paste(as.character(neg_sig$human_signature), collapse = ", ")))

fwrite(d, file.path(DATA_DIR, "cross_species_pathway_translatability_matrix.csv"))

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
lim <- max(abs(d$rho_NES), na.rm = TRUE)
p <- ggplot(d, aes(x = diet, y = human_signature, fill = rho_NES)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.2f", rho_NES)), size = 2.2) +
  scale_fill_gradient2(low = "#1565C0", mid = "#F5F5F5", high = "#C9265E",
                       midpoint = 0, limits = c(-lim, lim),
                       name = "NES rho") +
  labs(x = "Mouse diet model", y = NULL,
       title = "Cross-species pathway translatability") +
  coord_equal() +
  theme_masld(base_size = 7) +
  theme(
    plot.title      = element_text(size = 7.3, face = "bold", margin = margin(b = 6)),
    axis.text.x     = element_text(size = 6.5, face = "bold"),
    axis.text.y     = element_text(size = 6.5, face = "bold"),
    panel.grid      = element_blank(),
    legend.position = "right",
    legend.key.size = unit(0.3, "cm")
  )

ggsave(OUT_PDF, p,
       width  = 90 / 25.4,
       height = 70 / 25.4,
       units  = "in",
       device = cairo_pdf)
cat(sprintf("[saved] %s\n", OUT_PDF))
