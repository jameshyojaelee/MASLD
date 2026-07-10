##############################################################################
# Figure 4 panels (f)-(h): Spatial transcriptomics validation
#
# Narrative: discover -> localize -> integrate
#
# Two independent Visium datasets:
#   Guilliams et al. (Cell 2022) — fresh-frozen, 6,546 spots
#   Vu et al. (JHEP Reports 2025) — FFPE CytAssist, 17,512 spots
#
# Panels saved (g)-(h):
#   (g) Gene expression dynamics along pseudotime
#   (h) Spatial enrichment summary — dot plot
#
# Panel (f) Disease-emergent SVGs (Moran's I scatter) and the (f)-(h) composite
# were retired 2026-07-07 (superseded by the Fig4 spatial candidate gallery;
# see memory project-fig4-panelA-validation-wheel-2026-07-01) — p_f is still
# computed above but no longer saved; do not re-add its save_fig() call.
#
# Panels moved to supplementary (figS_spatial_validation.R):
#   Cross-dataset zonation concordance scatter
#   Periportal excess replicated — grouped bar
#   Fibrosis LR pairs spatial co-localization
##############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(scales)
  library(ggrepel)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
dir.create(file.path(FIG4_DIR, "panels"), showWarnings = FALSE, recursive = TRUE)

SPATIAL  <- file.path(BASE, "Analysis/Spatial/results")
VAL_DIR  <- file.path(SPATIAL, "validation_bulk")

# Dataset colors and labels
vu_col <- "#C9265E"
gu_col <- "#1565C0"
VU <- "Vu et al."
GU <- "Guilliams et al."
ds_pal <- setNames(c(gu_col, vu_col), c(GU, VU))

# Star helper
star <- function(p) {
  ifelse(p < 0.001, "***",
  ifelse(p < 0.01,  "**",
  ifelse(p < 0.05,  "*", "n.s.")))
}

# ==========================================================================
# (f) Disease-emergent SVGs: Moran's I scatter
# ==========================================================================
cat("Panel f: Disease-emergent SVGs...\n")

dsvg <- fread(file.path(SPATIAL, "svg/differential_svgs.csv"))

# Fix first column if unnamed
if (names(dsvg)[1] %in% c("V1", "")) setnames(dsvg, 1, "gene")
if (!("gene" %in% names(dsvg))) setnames(dsvg, 1, "gene")

# Coerce booleans
for (col in c("svg_healthy", "svg_masld")) {
  if (is.character(dsvg[[col]])) dsvg[, (col) := get(col) == "True"]
}

# Classify
dsvg[, display_cat := fcase(
  grepl("emergent", category, ignore.case = TRUE), "Disease-emergent",
  grepl("lost|resolved", category, ignore.case = TRUE), "Disease-resolved",
  default = "Stable SVG"
)]

# Only SVGs in at least one condition
dsvg_plot <- dsvg[svg_healthy == TRUE | svg_masld == TRUE]

cat_colors_f <- c(
  `Disease-emergent` = spatial_colors[["disease_emergent"]],
  `Disease-resolved` = spatial_colors[["disease_lost"]],
  `Stable SVG`       = spatial_colors[["stable_svg"]]
)

# Labels: top emergent + top resolved
top_em <- dsvg_plot[display_cat == "Disease-emergent"][order(-delta_I)][1:min(8, .N)]
top_rs <- dsvg_plot[display_cat == "Disease-resolved"][order(delta_I)][1:min(4, .N)]
top_lab_f <- rbind(top_em, top_rs)

# Count for legend
cat_n <- dsvg_plot[, .N, by = display_cat]

p_f <- ggplot(dsvg_plot[display_cat == "Stable SVG"],
              aes(x = morans_I_healthy, y = morans_I_masld)) +
  geom_point(color = spatial_colors[["stable_svg"]], size = 0.3, alpha = 0.3, shape = 16) +
  geom_point(data = dsvg_plot[display_cat != "Stable SVG"],
             aes(color = display_cat), size = 0.8, alpha = 0.7, shape = 16) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              linewidth = 0.3, color = "gray50") +
  geom_label_repel(data = top_lab_f, aes(label = gene, color = display_cat),
                   size = GEOM_TEXT_6PT, max.overlaps = 20,
                   label.padding = 0.07, box.padding = 0.25,
                   segment.size = 0.1, fill = alpha("white", 0.85),
                   show.legend = FALSE) +
  scale_color_manual(
    values = cat_colors_f, name = NULL,
    labels = setNames(
      paste0(cat_n$display_cat, " (n=", cat_n$N, ")"),
      cat_n$display_cat
    )
  ) +
  labs(x = expression("Moran's " * italic(I) * " (Healthy)"),
       y = expression("Moran's " * italic(I) * " (Steatotic)")) +
  theme_masld() +
  theme(legend.position = c(0.25, 0.88),
        legend.background = element_rect(fill = alpha("white", 0.9), color = NA),
        legend.key.size = unit(0.2, "cm")) +
  guides(color = guide_legend(override.aes = list(size = 1.5, alpha = 1)))

message("[caption] Disease creates spatially organized gene programs")
cat("  ", nrow(dsvg_plot), "SVGs plotted\n")

# ==========================================================================
# (g) Gene expression dynamics along disease pseudotime
# ==========================================================================
cat("Panel g: Pseudotime dynamics...\n")

traj <- fread(file.path(SPATIAL, "trajectory/gene_dynamics_along_trajectory.csv"))

# Fix unnamed index column
if (names(traj)[1] %in% c("V1", "")) setnames(traj, 1, "idx")

# Select genes for display: lipogenesis up, zonation marker, fibrosis, periportal down
show_genes <- c("FASN", "SCD", "CYP2E1", "COL1A1", "HAL")
traj_sub <- traj[gene %in% show_genes]
traj_sub[, gene := factor(gene, levels = show_genes)]

# Ribbon = spot-level dispersion (+/- SD), NOT a standard error / confidence
# band (F192). The trajectory CSV is pre-aggregated per gene x bin over Visium
# SPOTS (~134 spots/bin) of only ~5 GSE192741 donors, so spots are
# pseudoreplicates: dividing std_expr by sqrt(n_spots) would inflate the
# effective n ~25-130x and give an artificially tight CI. We plot the spot SD
# as a descriptive dispersion band instead. (A donor-level SE would require the
# per-donor trajectory, not present in this aggregated file.)
traj_sub[, sd_band := std_expr]

# Gene-specific colors
gene_pal <- c(
  FASN   = "#C2185B",
  SCD    = "#E91E63",
  CYP2E1 = "#7B1FA2",
  COL1A1 = "#880E4F",
  HAL    = "#1565C0"
)

p_g <- ggplot(traj_sub, aes(x = bin, y = mean_expr, color = gene, fill = gene)) +
  geom_ribbon(aes(ymin = mean_expr - sd_band, ymax = mean_expr + sd_band),
              alpha = 0.15, linewidth = 0) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 0.8, shape = 16) +
  scale_color_manual(values = gene_pal, name = NULL) +
  scale_fill_manual(values = gene_pal, name = NULL) +
  # Trajectory CSV has 20 pseudotime bins (0-19); the disease endpoint is bin 19,
  # not bin 9 (F191). Span the full data range so "Steatotic" sits at the end.
  scale_x_continuous(breaks = c(0, 9, 19),
                     labels = c("Healthy-like", "Transition", "Steatotic")) +
  labs(x = "Spatial pseudotime", y = "Mean expression") +
  theme_masld() +
  theme(legend.position = c(0.15, 0.85),
        legend.background = element_rect(fill = alpha("white", 0.9), color = NA),
        legend.key.size = unit(0.25, "cm"),
        axis.text.x = element_text(size = 6))

message("[caption] Gene dysregulation follows a spatial gradient")

# ==========================================================================
# (h) Spatial enrichment summary — dot plot
# ==========================================================================
cat("Panel h: Spatial enrichment...\n")

enrich   <- fread(file.path(VAL_DIR, "deg_spatial_enrichment.csv"))
# The on-disk file is conserved_core_spatial.csv with comparison=="Conserved_Core"
# (F190/F245); we relabel to "Conserved" in the select below so the line-183
# filter still matches.
cc_test  <- fread(file.path(VAL_DIR, "conserved_core_spatial.csv"))

# Combine (relabel Conserved_Core -> Conserved for the panel)
cc_fmt <- cc_test[, .(dataset, comparison = "Conserved",
                       fold_enrichment, mannwhitney_pval)]
enrich_slim <- enrich[, .(dataset, comparison, fold_enrichment, mannwhitney_pval)]
eall <- rbind(enrich_slim, cc_fmt, fill = TRUE)
# Match the "Strong DEGs" row by prefix so the panel renders regardless of the
# |LFC| threshold string carried in deg_spatial_enrichment.csv. NOTE: the on-disk
# enrichment is currently computed on the RAW |LFC|>0.5 DEG set (upstream spatial
# generator, not regenerated against the relaxed |shrunk_logFC|>0.3 Tier-1 set);
# the displayed label is the threshold-agnostic "Strong DEGs".
eall <- eall[grepl("^Strong DEGs", comparison) |
             comparison %in% c("Up-regulated DEGs", "Conserved")]

# Clean labels
eall[, label := fcase(
  grepl("^Strong DEGs", comparison), "Strong DEGs",
  comparison == "Up-regulated DEGs", "Up-regulated DEGs",
  comparison == "Conserved", "Conserved"
)]
eall[, label := factor(label,
                       levels = rev(c("Strong DEGs", "Up-regulated DEGs", "Conserved")))]
eall[, sig_star := star(mannwhitney_pval)]
eall[, nlp := -log10(pmax(mannwhitney_pval, 1e-30))]
eall[, sig := mannwhitney_pval < 0.05]

p_h <- ggplot(eall, aes(x = fold_enrichment, y = label)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray70", linewidth = 0.3) +
  geom_point(aes(color = dataset, size = nlp, shape = sig)) +
  geom_text(aes(label = sig_star, color = dataset),
            hjust = -0.5, vjust = 0.3, size = GEOM_TEXT_6PT, show.legend = FALSE) +
  # Display Guilliams Visium dataset by its GEO accession (data keys unchanged
  # so the color join against `dataset` still matches; Vu has no mapped accession).
  scale_color_manual(values = ds_pal, name = NULL,
                     labels = c(`Guilliams et al.` = "GSE192741",
                                `Vu et al.` = "Vu et al.")) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
  scale_size_continuous(range = c(1.5, 4),
                        name = expression(-log[10] ~ italic(p)),
                        breaks = c(5, 15)) +
  scale_x_continuous(breaks = c(1, 1.5, 2, 2.5, 3)) +
  labs(x = "Spatial enrichment (fold over non-DEGs)", y = NULL) +
  theme_masld() +
  theme(legend.position = "right", legend.key.size = unit(0.25, "cm"))

message("[caption] Cross-species conserved genes are most spatially structured")

# ==========================================================================
# Individual panels for Illustrator
# ==========================================================================
# NOTE: the composite fig4_spatial_panels.pdf and panels (f) f_emergent_svgs.pdf,
# (g) g_pseudotime.pdf, (h) h_enrichment.pdf were retired 2026-07-07 (superseded
# by the Fig4 spatial candidate gallery; see memory
# project-fig4-panelA-validation-wheel-2026-07-01) — do not re-add their
# save_fig() calls. p_f/p_g/p_h are still computed above and sourced in-memory
# by fig4_compact.R for the fig4_compact.pdf composite spatial row.

cat("Done (no standalone panel PDFs — p_f/p_g/p_h consumed by fig4_compact.R).\n")
