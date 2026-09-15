#!/usr/bin/env Rscript

# KEY MESSAGE: Cell-resolved MASLD modules have assay-specific regulatory,
# protein, and spatial outcomes; absence of support is not conflated with lack
# of assay coverage and no universal convergence score is constructed.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

BASE <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

FS <- 6
INK <- "#222222"
NEG <- "#2166AC"
POS <- "#B2182B"
ZERO <- "#F4F4F4"
UNTESTED <- "#BDBDBD"
LINEAGE <- c(
  hepatocytes = "#4C78A8",
  fibroblasts = "#E45756",
  macrophages = "#7A5195",
  cholangiocytes = "#59A14F"
)

theme_fig4 <- theme_minimal(base_size = FS, base_family = "Helvetica") +
  theme(
    text = element_text(size = FS, face = "plain", colour = INK),
    panel.grid = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_text(size = FS, face = "plain", colour = INK),
    axis.title = element_text(size = FS, face = "plain", colour = INK),
    legend.text = element_text(size = FS, face = "plain"),
    legend.title = element_text(size = FS, face = "plain"),
    strip.text = element_text(size = FS, face = "plain"),
    plot.title = element_text(size = FS, face = "plain", hjust = 0),
    plot.subtitle = element_text(size = FS, face = "plain", hjust = 0),
    legend.key = element_blank(),
    legend.background = element_blank(),
    plot.margin = margin(3, 3, 3, 3)
  )

wide_file <- file.path(PROGRAM_CONTEXT_DIR, "program_context_wide.tsv")
long_file <- file.path(PROGRAM_CONTEXT_DIR, "program_context_long.tsv")
if (!file.exists(wide_file) || !file.exists(long_file)) {
  stop("Run 05_build_program_context.R before drawing Figure 4 panels")
}
wide <- fread(wide_file)
long <- fread(long_file)
stopifnot(nrow(wide) == 22L, !anyDuplicated(wide$program_id))

abbr <- c(
  hepatocytes = "Hep",
  fibroblasts = "Fib",
  macrophages = "Mac",
  cholangiocytes = "Chol"
)
wide[, row_label := paste0(
  abbr[cell_type], " | ", program_name,
  fifelse(stability_fail, " [S]", "")
)]
row_levels <- rev(wide$row_label)
label_map <- setNames(wide$row_label, wide$program_id)

panel_filter <- tolower(trimws(Sys.getenv("FIG4_PANEL", "all")))
if (!panel_filter %in% c("all", "a", "b", "d", "e")) {
  stop("FIG4_PANEL must be one of: all, a, b, d, e")
}

save_pdf <- function(plot, filename, width, height) {
  panel_code <- sub("^fig4_([a-z]).*$", "\\1", filename)
  if (panel_filter != "all" && panel_code != panel_filter) {
    return(invisible(NULL))
  }
  out <- file.path(FIG4_DIR, filename)
  ggsave(
    out,
    plot,
    width = width,
    height = height,
    units = "in",
    device = "pdf",
    useDingbats = FALSE
  )
  message("[fig4] saved: ", out)
}

# -------------------------------------------------------------------------
# 4A: input separation. Genetics and the 22 Hotspot programs are inherited from
# Figures 2/3; protein process sets remain an assay-native parallel branch.
# -------------------------------------------------------------------------
boxes <- data.table(
  xmin = c(0.05, 0.05, 0.05, 2.05, 2.05, 2.05, 4.05),
  xmax = c(1.45, 1.45, 1.45, 3.45, 3.45, 3.45, 5.35),
  ymin = c(2.15, 1.15, 0.15, 2.15, 1.15, 0.15, 0.65),
  ymax = c(2.85, 1.85, 0.85, 2.85, 1.85, 0.85, 2.35),
  label = c(
    "Fig. 2\n35-study genetics",
    "Fig. 3\n117 -> 22 donor modules",
    "PXD051911-selected\n25-protein process rows",
    "Fine-mapped variants\nin open chromatin",
    "Frozen-module ATAC,\nDIA-MS, and spatial",
    "Selected-process liver\nDIA-MS display (4C)",
    "22-row module matrix\n(no modality sum)"
  ),
  fill = c("#E8EEF7", "#EEF5EA", "#FAEFE4", "#E8EEF7", "#EEF5EA", "#FAEFE4", "#F2F2F2")
)
arrows <- data.table(
  x = c(1.45, 1.45, 1.45, 3.45, 3.45),
  xend = c(2.05, 2.05, 2.05, 4.05, 4.05),
  y = c(2.50, 1.50, 0.50, 2.50, 1.50),
  yend = c(2.50, 1.50, 0.50, 1.95, 1.35)
)
p4a <- ggplot() +
  geom_rect(
    data = boxes,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = boxes$fill,
    colour = "#777777",
    linewidth = 0.25
  ) +
  geom_segment(
    data = arrows,
    aes(x = x, xend = xend, y = y, yend = yend),
    colour = "#666666",
    linewidth = 0.3,
    arrow = arrow(length = grid::unit(0.05, "in"), type = "closed")
  ) +
  geom_text(
    data = boxes,
    aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = label),
    size = FS / .pt,
    lineheight = 0.9,
    family = "Helvetica",
    fontface = "plain"
  ) +
  annotate(
    "text", x = 0.05, y = 3.05,
    label = "Fixed analysis inputs", hjust = 0,
    size = FS / .pt, family = "Helvetica"
  ) +
  annotate(
    "text", x = 2.05, y = 3.05,
    label = "Assay-native analyses", hjust = 0,
    size = FS / .pt, family = "Helvetica"
  ) +
  coord_cartesian(xlim = c(0, 5.4), ylim = c(0, 3.15), expand = FALSE, clip = "off") +
  theme_void(base_size = FS, base_family = "Helvetica") +
  theme(text = element_text(size = FS, face = "plain", colour = INK))
save_pdf(p4a, "fig4_a_input_firewall.pdf", 5.3, 1.55)

# -------------------------------------------------------------------------
# 4B: two-layer ATAC. Counts are descriptive promoter/open-peak loci; dynamic
# effects use donor pseudobulk and exact label permutations.
# -------------------------------------------------------------------------
static <- melt(
  wide,
  id.vars = c("program_id", "row_label", "cell_type", "display_order"),
  measure.vars = c("gwas_atac_direct_genes", "gwas_atac_enzyme_genes"),
  variable.name = "layer",
  value.name = "count"
)
static[, layer := factor(
  layer,
  levels = c("gwas_atac_direct_genes", "gwas_atac_enzyme_genes"),
  labels = c("Direct disease/PDFF", "Liver enzyme")
)]
static[, row_label := factor(row_label, levels = row_levels)]
p_static <- ggplot(static, aes(x = layer, y = row_label)) +
  geom_tile(aes(fill = sqrt(count)), colour = "white", linewidth = 0.25) +
  geom_text(aes(label = count), size = FS / .pt, family = "Helvetica") +
  scale_fill_gradient(low = "white", high = "#4C78A8", name = "sqrt(genes)") +
  labs(x = NULL, y = NULL) +
  theme_fig4 +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

dynamic <- rbindlist(list(
  wide[, .(
    program_id, row_label, display_order,
    cohort = "External GSE281367",
    effect = atac_external_effect,
    testable = atac_external_testable,
    robust = atac_external_robust
  )],
  wide[, .(
    program_id, row_label, display_order,
    cohort = "Internal GSE244832",
    effect = atac_internal_effect,
    testable = atac_internal_testable,
    robust = atac_internal_robust
  )]
), fill = TRUE)
dynamic[, status := fifelse(!testable | is.na(testable), "Untested", fifelse(robust, "Robust", "Indeterminate"))]
dynamic[, row_label := factor(row_label, levels = row_levels)]
dynamic[, cohort := factor(cohort, levels = c("External GSE281367", "Internal GSE244832"))]
atac_r <- cor(wide$atac_external_effect, wide$atac_internal_effect, use = "complete.obs")
atac_sign_agree <- sum(
  sign(wide$atac_external_effect) == sign(wide$atac_internal_effect),
  na.rm = TRUE
)
p_dynamic <- ggplot(dynamic, aes(x = cohort, y = row_label)) +
  geom_point(
    data = dynamic[status == "Indeterminate"],
    aes(size = abs(effect), fill = effect), shape = 21,
    colour = "#777777", stroke = 0.25, alpha = 0.65
  ) +
  geom_point(
    data = dynamic[status == "Robust"],
    aes(size = abs(effect), fill = effect), shape = 21, colour = INK, stroke = 0.3
  ) +
  geom_point(
    data = dynamic[status == "Untested"],
    shape = 4, colour = UNTESTED, stroke = 0.35, size = 1.4
  ) +
  scale_fill_gradient2(low = NEG, mid = ZERO, high = POS, midpoint = 0, name = "score difference") +
  scale_size_continuous(range = c(0.7, 2.4), name = "abs(effect)") +
  labs(x = NULL, y = NULL) +
  theme_fig4 +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )
p4b <- (p_static + p_dynamic + plot_layout(widths = c(1.15, 1.0))) +
  plot_annotation(
    title = "Two complementary ATAC tests",
    subtitle = paste0(
      "Left, static overlap: descriptive promoter counts (credible sets for 32/35 studies)\n",
      "Right, signed donor-level effects; no individual module passes FDR (cohort r=",
      sprintf("%.2f", atac_r), "; sign agreement ", atac_sign_agree, "/22)"
    ),
    theme = theme(
      plot.title = element_text(size = FS, face = "plain", family = "Helvetica", colour = INK),
      plot.subtitle = element_text(size = FS, face = "plain", family = "Helvetica", colour = INK)
    )
  )
save_pdf(p4b, "fig4_b_two_layer_atac.pdf", 5.4, 5.25)

# -------------------------------------------------------------------------
# 4D: raw versus composition-adjusted spatial organization, plus the Vu array
# series and the descriptive GSE disease delta.
# -------------------------------------------------------------------------
spatial_long <- rbindlist(list(
  wide[, .(program_id, row_label, display_order, metric = "GSE raw Moran Z", value = spatial_gse_raw_moran_z, testable = spatial_gse_testable, robust = FALSE)],
  wide[, .(program_id, row_label, display_order, metric = "GSE residual Moran Z", value = spatial_gse_residual_moran_z, testable = spatial_gse_testable, robust = spatial_gse_robust)],
  wide[, .(program_id, row_label, display_order, metric = "GSE + zonation Moran Z", value = spatial_gse_zonation_moran_z, testable = spatial_gse_testable, robust = FALSE)],
  wide[, .(program_id, row_label, display_order, metric = "Vu residual Moran Z", value = spatial_vu_residual_moran_z, testable = spatial_vu_testable, robust = spatial_vu_robust)],
  wide[, .(program_id, row_label, display_order, metric = "GSE disease delta", value = spatial_gse_disease_delta_descriptive, testable = spatial_gse_testable, robust = FALSE)]
), fill = TRUE)
spatial_long[, metric := factor(
  metric,
  levels = c("GSE raw Moran Z", "GSE residual Moran Z", "GSE + zonation Moran Z", "Vu residual Moran Z", "GSE disease delta")
)]
spatial_long[, row_label := factor(row_label, levels = row_levels)]
spatial_long[, status := fifelse(!testable | is.na(testable), "Untested", fifelse(robust, "Robust", "Tested/descriptive"))]
spatial_z <- spatial_long[metric != "GSE disease delta"]
spatial_delta <- spatial_long[metric == "GSE disease delta"]
lim_z <- quantile(abs(spatial_z$value), 0.98, na.rm = TRUE)
lim_delta <- max(abs(spatial_delta$value), na.rm = TRUE)
if (!is.finite(lim_z) || lim_z == 0) lim_z <- 1
if (!is.finite(lim_delta) || lim_delta == 0) lim_delta <- 1
p_spatial_z <- ggplot(spatial_z, aes(x = metric, y = row_label)) +
  geom_tile(aes(fill = pmax(pmin(value, lim_z), -lim_z)), colour = "white", linewidth = 0.2) +
  geom_point(data = spatial_z[status == "Robust"], shape = 21, fill = NA, colour = INK, size = 1.35, stroke = 0.4) +
  geom_point(data = spatial_z[status == "Untested"], shape = 4, colour = UNTESTED, size = 1.25, stroke = 0.35) +
  scale_fill_gradient2(low = NEG, mid = ZERO, high = POS, midpoint = 0, limits = c(-lim_z, lim_z), oob = squish, name = "matched-null Z") +
  labs(
    x = NULL, y = NULL
  ) +
  theme_fig4 +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")
p_spatial_delta <- ggplot(spatial_delta, aes(x = metric, y = row_label)) +
  geom_tile(aes(fill = pmax(pmin(value, lim_delta), -lim_delta)), colour = "white", linewidth = 0.2) +
  geom_point(data = spatial_delta[status == "Untested"], shape = 4, colour = UNTESTED, size = 1.25, stroke = 0.35) +
  scale_fill_gradient2(low = NEG, mid = ZERO, high = POS, midpoint = 0, limits = c(-lim_delta, lim_delta), oob = squish, name = "score difference") +
  labs(
    x = NULL, y = NULL
  ) +
  theme_fig4 +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )
p4d <- (p_spatial_z + p_spatial_delta + plot_layout(widths = c(2.8, 1.0))) +
  plot_annotation(
    title = "Spatial organization beyond matched lineage abundance",
    subtitle = paste0(
      "Lineage-specific matched-gene null; graphs split by connected tissue island\n",
      "GSE: 4 donors after H35 collapse; Vu: 10 arrays (biopsy map unavailable)\n",
      "Right: descriptive GSE delta"
    ),
    theme = theme(
      plot.title = element_text(size = FS, face = "plain", family = "Helvetica", colour = INK),
      plot.subtitle = element_text(size = FS, face = "plain", family = "Helvetica", colour = INK)
    )
  )
save_pdf(p4d, "fig4_d_spatial_organization.pdf", 5.4, 5.25)

# -------------------------------------------------------------------------
# 4E: compact all-program status matrix. Assay-native magnitudes remain in the
# source table and in the detailed panels; this synthesis shows only testability
# and prespecified robustness so incomparable assay units are not conflated.
# -------------------------------------------------------------------------
wide[, program_label := paste0(
  fifelse(fig3_direction == "up", "+ ", "- "),
  program_name,
  fifelse(!is.na(stability_fail) & stability_fail, " *", "")
)]
wide[, lineage_label := factor(
  cell_type,
  levels = c("hepatocytes", "fibroblasts", "macrophages", "cholangiocytes"),
  labels = c("Hepatocyte", "Fibroblast", "Macrophage", "Cholangiocyte")
)]

# A valid test whose robust rule does not pass is INDETERMINATE, never
# tested_negative: that state requires a prespecified adequate-negative rule
# that passed, and no such rule is defined for these assays.
status_from <- function(testable, robust) {
  fifelse(
    is.na(testable) | !testable,
    "untested",
    fifelse(!is.na(robust) & robust, "robust", "indeterminate")
  )
}

base_cols <- c(
  "program_id", "display_order", "lineage_label", "program_label"
)
open_cells <- wide[, c(
  mget(base_cols),
  list(
    column = "open_chromatin",
    count = fcoalesce(as.integer(gwas_atac_direct_genes), 0L) +
      fcoalesce(as.integer(gwas_atac_enzyme_genes), 0L)
  )
)]
open_cells[, status := fifelse(count > 0L, "overlap", "no_overlap")]

assay_cells <- function(column_name, testable_col, robust_col) {
  wide[, c(
    mget(base_cols),
    list(
      column = column_name,
      count = NA_integer_,
      status = status_from(get(testable_col), get(robust_col))
    )
  )]
}

matrix_e <- rbindlist(list(
  open_cells,
  assay_cells("atac_external", "atac_external_testable", "atac_external_robust"),
  assay_cells("atac_internal", "atac_internal_testable", "atac_internal_robust"),
  assay_cells("protein", "protein_testable", "protein_robust"),
  assay_cells("spatial_gse", "spatial_gse_testable", "spatial_gse_robust"),
  assay_cells("spatial_vu", "spatial_vu_testable", "spatial_vu_robust")
), use.names = TRUE)

column_levels <- c(
  "open_chromatin", "atac_external", "atac_internal",
  "protein", "spatial_gse", "spatial_vu"
)
matrix_e[, column := factor(column, levels = column_levels)]
matrix_e[, program_id := factor(program_id, levels = rev(wide$program_id))]

robust_fraction <- function(column_name) {
  d <- matrix_e[as.character(column) == column_name]
  sprintf(
    "%d/%d robust",
    sum(d$status == "robust"),
    sum(d$status != "untested")
  )
}
x_labels <- c(
  open_chromatin = sprintf(
    "Chromatin\nopen peaks\n%d/%d overlap",
    sum(open_cells$status == "overlap"), nrow(open_cells)
  ),
  atac_external = paste("ATAC", "GSE281367", robust_fraction("atac_external"), sep = "\n"),
  atac_internal = paste("ATAC", "GSE244832", robust_fraction("atac_internal"), sep = "\n"),
  protein = paste("Protein", "PXD051911", robust_fraction("protein"), sep = "\n"),
  spatial_gse = paste("Spatial", "GSE192741", robust_fraction("spatial_gse"), sep = "\n"),
  spatial_vu = paste("Spatial", "Vu et al.", robust_fraction("spatial_vu"), sep = "\n")
)
y_labels <- setNames(wide$program_label, wide$program_id)

p4e <- ggplot(matrix_e, aes(x = column, y = program_id)) +
  geom_vline(
    xintercept = c(1.5, 3.5, 4.5),
    colour = "#BDBDBD", linewidth = 0.2
  ) +
  geom_tile(
    width = 0.72, height = 0.72,
    fill = "white", colour = "#E1E1E1", linewidth = 0.25
  ) +
  geom_tile(
    data = matrix_e[status == "indeterminate"],
    width = 0.72, height = 0.72,
    fill = "#E2E2E2", colour = "#C8C8C8", linewidth = 0.25
  ) +
  geom_tile(
    data = matrix_e[status == "robust"],
    width = 0.72, height = 0.72,
    fill = "#2F6F73", colour = INK, linewidth = 0.3
  ) +
  geom_tile(
    data = matrix_e[status == "overlap"],
    width = 0.72, height = 0.72,
    fill = "#F1E2C2", colour = "#8C6D31", linewidth = 0.3
  ) +
  geom_text(
    data = matrix_e[status == "overlap"],
    aes(label = count),
    size = FS / .pt, family = "Helvetica", colour = INK
  ) +
  geom_text(
    data = matrix_e[status == "untested"],
    label = "x", size = FS / .pt, family = "Helvetica", colour = UNTESTED
  ) +
  facet_grid(
    rows = vars(lineage_label),
    scales = "free_y", space = "free_y", switch = "y"
  ) +
  scale_x_discrete(
    limits = column_levels, labels = x_labels,
    position = "top", drop = FALSE, expand = expansion(add = 0.25)
  ) +
  scale_y_discrete(labels = y_labels, drop = TRUE, expand = expansion(add = 0.18)) +
  labs(
    x = NULL, y = NULL,
    title = "Assay status for frozen Figure 3 programs",
    caption = paste0(
      "+/- = Figure 3 trajectory direction; * = stability flag. ",
      "Teal = robust; gray = indeterminate (valid test, robust rule not met); ",
      "number = open-peak genes; ",
      "x = not testable.\n",
      "Open-peak overlaps are confined to liver-enzyme loci; direct disease/PDFF overlap = 0/22."
    )
  ) +
  theme_fig4 +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, lineheight = 0.92),
    axis.text.y = element_text(lineheight = 0.92),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0, hjust = 1, vjust = 0.5),
    panel.spacing.y = grid::unit(0.06, "in"),
    legend.position = "none",
    plot.caption = element_text(size = FS, face = "plain", hjust = 0, lineheight = 0.9),
    plot.caption.position = "plot",
    plot.margin = margin(3, 4, 3, 3)
  )
save_pdf(p4e, "fig4_e_program_context_matrix.pdf", 5.8, 5.15)
