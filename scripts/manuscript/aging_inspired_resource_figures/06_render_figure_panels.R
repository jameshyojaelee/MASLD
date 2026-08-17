#!/usr/bin/env Rscript
# Render individual editable-vector panels from the validated all-117 contracts.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(gridExtra)
})

ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
source(file.path(ROOT, "scripts/figures/publication_theme.R"))
OUT <- Sys.getenv("FIG_CAND_ROOT")
if (!nzchar(OUT)) stop("FIG_CAND_ROOT is required")
SOURCE <- file.path(OUT, "source_tables")
PANELS <- file.path(OUT, "panels")
PROOFS <- file.path(OUT, "proofs")
dir.create(PANELS, recursive = TRUE, showWarnings = FALSE)
dir.create(PROOFS, recursive = TRUE, showWarnings = FALSE)
set.seed(20260815)
grDevices::pdf.options(useDingbats = FALSE)

save_vector <- function(plot, path, width, height) {
  if (file.exists(path)) stop("Refusing to overwrite: ", path)
  ggsave(path, plot, width = width, height = height,
         device = grDevices::cairo_pdf, limitsize = FALSE)
}

stack_with_title <- function(plot, title, caption = NULL) {
  pieces <- list(grid::textGrob(
    title, x = 0, hjust = 0,
    gp = grid::gpar(fontfamily = "Helvetica", fontsize = 6)
  ), plot)
  heights <- c(0.07, 0.86)
  if (!is.null(caption)) {
    pieces <- c(pieces, list(grid::textGrob(
      caption, x = 0, hjust = 0,
      gp = grid::gpar(fontfamily = "Helvetica", fontsize = 6)
    )))
    heights <- c(heights, 0.07)
  }
  arrangeGrob(grobs = pieces, ncol = 1, heights = heights)
}

lineage_colors <- c(
  cholangiocytes = "#518DC9",
  fibroblasts = "#E69F00",
  hepatocytes = "#0072B2",
  macrophages = "#C9265E",
  tcells = "#7B1FA2",
  "not BH q<0.05" = "#9E9E9E"
)
state_colors <- c(
  supported = "#00796B",
  descriptive = "#4C78A8",
  sensitive = "#E69F00",
  indeterminate = "#9E9E9E"
)
plain_theme <- theme_masld() + theme_pub() +
  theme(plot.title.position = "plot", plot.caption.position = "plot")

# 4A: derive, interpret, stress-test, transport.
# KEY MESSAGE: Fixed biological objects are interpreted, stress-tested, and transported.
n4a <- fread(file.path(SOURCE, "fig4a_analysis_logic_nodes.tsv"))
e4a <- fread(file.path(SOURCE, "fig4a_analysis_logic_edges.tsv"))
e4a <- merge(e4a, n4a[, .(from = node_id, x, y)], by = "from")
setnames(e4a, c("x", "y"), c("x_from", "y_from"))
e4a <- merge(e4a, n4a[, .(to = node_id, x, y)], by = "to")
setnames(e4a, c("x", "y"), c("x_to", "y_to"))
p4a <- ggplot() +
  geom_segment(
    data = e4a,
    aes(x = x_from, y = y_from, xend = x_to, yend = y_to),
    linewidth = 0.35, colour = "#5A5A5A",
    arrow = grid::arrow(length = grid::unit(1.4, "mm"), type = "closed")
  ) +
  geom_label(
    data = n4a,
    aes(x = x, y = y, label = paste0(title, "\n", detail)),
    family = "Helvetica", size = GEOM_TEXT_6PT, label.size = 0.25,
    label.padding = grid::unit(1.2, "mm"), fill = "white", colour = "black",
    lineheight = 0.95
  ) +
  annotate("text", x = 4.5, y = 1.42,
           label = "Identity stress test", family = "Helvetica",
           size = GEOM_TEXT_6PT, colour = "black") +
  coord_cartesian(xlim = c(-0.65, 6.65), ylim = c(-0.42, 1.55), clip = "off") +
  labs(
    title = "Fixed biological objects are interpreted, stress-tested, and transported",
    caption = "GSE189600 remains visible as uncorrected passthrough in the full-universe sensitivity."
  ) +
  theme_void(base_family = "Helvetica", base_size = 6) +
  theme(
    text = element_text(family = "Helvetica", size = 6, face = "plain", colour = "black"),
    plot.title = element_text(size = 6, face = "plain"),
    plot.caption = element_text(size = 6, face = "plain", hjust = 0),
    plot.margin = margin(5, 5, 5, 5)
  )
save_vector(p4a, file.path(PANELS, "fig4a_analysis_logic.pdf"), 7.1, 2.15)

# 4C: all 117 programs; only the prespecified examples and controls are labelled.
# KEY MESSAGE: All 117 frozen programs occupy a continuous disease-effect landscape.
d4c <- fread(file.path(SOURCE, "fig4c_all117_disease_skyline.tsv"))
d4c[, evidence_plot := pmin(family_evidence, 12)]
d4c[, point_group := fifelse(disease_qvalue < 0.05, cell_type, "not BH q<0.05")]
d4c[, cell_type := factor(cell_type,
                          levels = c("hepatocytes", "cholangiocytes", "fibroblasts", "macrophages", "tcells"),
                          labels = c("Hepatocytes", "Cholangiocytes", "Fibroblasts", "Macrophages", "T cells"))]
lab4c <- d4c[label_requested == TRUE & nzchar(direct_label)]
p4c <- ggplot(d4c, aes(disease_beta, evidence_plot)) +
  geom_hline(yintercept = -log10(0.05), linewidth = 0.25, linetype = 2, colour = "#9E9E9E") +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "#BDBDBD") +
  geom_point(aes(colour = point_group), size = 1.45, alpha = 0.9) +
  geom_text(
    data = lab4c,
    aes(label = direct_label, hjust = ifelse(disease_beta >= 0, 1.08, -0.08)),
    family = "Helvetica", size = GEOM_TEXT_6PT, colour = "black", check_overlap = TRUE
  ) +
  facet_wrap(~cell_type, nrow = 1, scales = "free_x") +
  scale_colour_manual(values = lineage_colors, drop = FALSE) +
  scale_y_continuous(breaks = c(0, 2, 4, 8, 12), labels = c("0", "2", "4", "8", "≥12")) +
  coord_cartesian(clip = "off") +
  labs(
    title = "All 117 frozen programs occupy a continuous disease-effect landscape",
    x = "Disease association β (Hotspot score per stage ordinal)",
    y = expression(-log[10]~"BH q (capped)"), colour = NULL
  ) +
  plain_theme +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.margin = margin(3, 10, 3, 3))
save_vector(p4c, file.path(PANELS, "fig4c_all117_disease_skyline.pdf"), 7.1, 2.55)

# 4D: correction transport plus dataset status in the same panel contract.
# KEY MESSAGE: Ambient correction changes some disease effects more than others.
d4d <- fread(file.path(SOURCE, "fig4d_ambient_effect_transport.tsv"))
d4d[, state_class := fifelse(grepl("retained|reversal_supported", evidence_state), "supported",
                      fifelse(grepl("ambient_sensitive", evidence_state), "sensitive", "indeterminate"))]
lab4d <- d4d[label_requested == TRUE & nzchar(direct_label)]
lims4d <- range(c(d4d$raw_beta, d4d$corrected_beta), finite = TRUE)
pad4d <- diff(lims4d) * 0.08
lims4d <- lims4d + c(-pad4d, pad4d)
p4d_main <- ggplot(d4d, aes(raw_beta, corrected_beta)) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.3, colour = "#9E9E9E") +
  geom_segment(
    data = lab4d,
    aes(x = raw_beta, y = raw_beta, xend = raw_beta, yend = corrected_beta),
    inherit.aes = FALSE, colour = "black", linewidth = 0.35,
    arrow = grid::arrow(length = grid::unit(1.1, "mm"), type = "closed")
  ) +
  geom_point(data = d4d[is.finite(corrected_beta)],
             aes(colour = state_class, shape = state_class), size = 1.55) +
  geom_point(data = d4d[!is.finite(corrected_beta)],
             aes(x = raw_beta, y = lims4d[1]), inherit.aes = FALSE,
             shape = 4, size = 1.7, stroke = 0.45, colour = "#7B1FA2") +
  annotate("text", x = lims4d[1], y = lims4d[1], hjust = 0, vjust = -0.8,
           label = paste0(sum(!is.finite(d4d$corrected_beta)),
                          " T-cell programs untestable"),
           family = "Helvetica", size = GEOM_TEXT_6PT, colour = "black") +
  geom_text(
    data = lab4d, aes(label = direct_label), family = "Helvetica",
    size = GEOM_TEXT_6PT, hjust = -0.08, vjust = 0.2, colour = "black",
    check_overlap = TRUE
  ) +
  scale_colour_manual(values = state_colors[c("supported", "sensitive", "indeterminate")]) +
  scale_shape_manual(values = c(supported = 16, sensitive = 17, indeterminate = 1)) +
  coord_equal(xlim = lims4d, ylim = lims4d, clip = "off") +
  labs(
    title = "Ambient correction changes some disease effects more than others",
    x = "Raw disease β", y = "Ambient-corrected disease β",
    colour = "Program state", shape = "Program state"
  ) +
  plain_theme +
  theme(legend.position = "right", plot.margin = margin(3, 14, 2, 3))

q4d <- fread(file.path(SOURCE, "fig4d_dataset_status_strip.tsv"))
q4d[, status_label := fifelse(correction_status == "corrected", "corrected", "uncorrected passthrough")]
q4d[, dataset := factor(dataset, levels = dataset)]
p4d_strip <- ggplot(q4d, aes(dataset, 1, fill = status_label)) +
  geom_tile(colour = "white", linewidth = 0.4, height = 0.65) +
  geom_text(aes(label = ifelse(dataset == "GSE189600", "GSE189600\npassthrough", as.character(dataset))),
            family = "Helvetica", size = GEOM_TEXT_6PT, colour = "black", lineheight = 0.9) +
  scale_fill_manual(values = c(corrected = "#8FC9BD", `uncorrected passthrough` = "#9E9E9E")) +
  labs(x = NULL, y = NULL, fill = "Dataset status") +
  theme_void(base_family = "Helvetica", base_size = 6) +
  theme(text = element_text(size = 6, family = "Helvetica", face = "plain", colour = "black"),
        legend.position = "bottom", legend.text = element_text(size = 6),
        legend.title = element_text(size = 6), plot.margin = margin(0, 3, 3, 3))
p4d <- arrangeGrob(p4d_main, p4d_strip, ncol = 1, heights = c(4.2, 1.0))
save_vector(p4d, file.path(PANELS, "fig4d_ambient_effect_transport.pdf"), 7.1, 4.2)

# 4E or S4E: the entire family follows the predeclared two-program promotion gate.
# KEY MESSAGE: Same-atlas contrasts localize IGFBP7 and BICC1 expression descriptively.
d4e <- fread(file.path(SOURCE, "fig4e_same_atlas_lineage_specificity.tsv"))
d4e[, comparison_lineage := factor(
  comparison_lineage,
  levels = rev(c("Hepatocytes", "Fibroblasts", "Cholangiocytes", "Endothelial cells", "Macrophages", "T cells"))
)]
d4e[, evidence := fifelse(comparison_lineage == "Hepatocytes", "reference",
                   fifelse(qvalue < 0.05, "BH q<0.05", "not BH q<0.05"))]
p4e <- ggplot(d4e, aes(beta, comparison_lineage)) +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "#9E9E9E") +
  geom_segment(aes(x = ci_low, xend = ci_high, yend = comparison_lineage), linewidth = 0.4) +
  geom_point(aes(fill = evidence), shape = 21, size = 2.1, stroke = 0.35) +
  facet_wrap(~program_name, nrow = 1, scales = "free_x") +
  scale_fill_manual(values = c(reference = "#9E9E9E", `BH q<0.05` = "#C9265E", `not BH q<0.05` = "white")) +
  labs(
    title = "Same-atlas descriptive localization",
    subtitle = "Dataset-adjusted donor-paired lineage-minus-hepatocyte contrasts; HC3 is a displayed robustness check",
    x = "Standardized fixed-weight program-score difference", y = NULL, fill = NULL,
    caption = paste0("Family destination: ", unique(d4e$destination), ". This analysis does not establish transcript origin or cell autonomy.")
  ) +
  plain_theme + theme(legend.position = "bottom")
name4e <- if (unique(d4e$destination) == "Figure_4E") {
  "fig4e_same_atlas_lineage_specificity.pdf"
} else {
  "figS4e_same_atlas_lineage_specificity.pdf"
}
save_vector(p4e, file.path(PANELS, name4e), 7.1, 2.65)

# 4F: separate axes prevent any visual conversion between single-cell and bulk units.
# KEY MESSAGE: Disease association and tissue-state transport align without combining units.
sc4f <- fread(file.path(SOURCE, "fig4f_singlecell_estimates.tsv"))
sc4f[, estimate_type := factor(estimate_type, levels = c("corrected", "raw"), labels = c("Ambient-corrected", "Raw"))]
p4f_sc <- ggplot(sc4f, aes(beta, estimate_type, colour = estimate_type)) +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "#9E9E9E") +
  geom_segment(aes(x = ci_low, xend = ci_high, yend = estimate_type), linewidth = 0.45) +
  geom_point(size = 2) +
  facet_wrap(~program_name, ncol = 1, scales = "free_x") +
  scale_colour_manual(values = c(`Ambient-corrected` = "#C9265E", Raw = "#4C78A8")) +
  labs(title = "Single-cell association", x = "Hotspot-score β per stage ordinal", y = NULL, colour = NULL) +
  plain_theme + theme(legend.position = "bottom")
bulk4f <- fread(file.path(SOURCE, "fig4f_bulk_stage_transport.tsv"))
bulk4f[, stage := factor(stage, levels = c("F1", "F2", "F3", "F4"))]
p4f_bulk <- ggplot(bulk4f, aes(stage, effect, group = program_uid)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "#9E9E9E") +
  geom_line(linewidth = 0.5, colour = "#4C78A8") +
  geom_point(size = 1.8, colour = "#4C78A8") +
  facet_wrap(~program_name, ncol = 1, scales = "free_y") +
  labs(title = "Bulk tissue-state transport", x = "Fibrosis stage versus F0",
       y = "L1-weighted member-gene log2FC") +
  plain_theme
p4f_body <- arrangeGrob(p4f_sc, p4f_bulk, nrow = 1, widths = c(1, 1.25))
p4f <- stack_with_title(
  p4f_body,
  "Disease association and tissue-state transport are aligned without combining their units",
  "Bulk points are fixed member-gene projections, not program-level significance tests."
)
save_vector(p4f, file.path(PANELS, "fig4f_tissue_state_transport.pdf"), 7.1, 4.0)

# 5F: representative maps paired with observed organization and its matched-gene
# null, without new statistics.
# KEY MESSAGE: Disease association and matched-null tissue organization are different claims.
d5f <- fread(file.path(SOURCE, "fig5f_spatial_map_and_matched_null.tsv"))
d5f[, dataset_label := factor(dataset_label, levels = rev(unique(dataset_label)))]
d5f[, state_class := fifelse(within_source_call == "supported", "supported", "indeterminate")]
p5f_null <- ggplot(d5f, aes(residual_moran_i, dataset_label)) +
  geom_segment(aes(x = q050, xend = q950, yend = dataset_label), colour = "#9E9E9E", linewidth = 2.0) +
  geom_point(aes(x = q500), shape = 23, size = 1.8, fill = "white", colour = "black") +
  geom_segment(aes(x = q500, xend = residual_moran_i, yend = dataset_label), linewidth = 0.35, linetype = 2) +
  geom_point(aes(colour = state_class), size = 2.2) +
  geom_text(aes(label = paste0("q=", formatC(primary_qvalue_spatial, format = "g", digits = 2))),
            family = "Helvetica", size = GEOM_TEXT_6PT, hjust = -0.1, colour = "black") +
  facet_wrap(~program_short, nrow = 1, scales = "free_x") +
  scale_colour_manual(values = state_colors[c("supported", "indeterminate")]) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Matched-null calibration",
    x = "Residual Moran I", y = NULL, colour = "Evidence state",
    caption = "Gray interval: matched-gene null 5th–95th percentiles; diamond: null median; colored point: observed."
  ) +
  plain_theme + theme(legend.position = "bottom", plot.margin = margin(3, 12, 3, 3))
map_paths <- unique(d5f$map_panel_path)
if (length(map_paths) != 1L || !file.exists(map_paths)) {
  stop("One existing representative-map panel is required")
}
map_png <- tempfile("fig5f_spatial_program_maps_", fileext = ".png")
pdftoppm <- Sys.which("pdftoppm")
if (!nzchar(pdftoppm)) stop("pdftoppm is required to embed the retained map panel")
map_stem <- sub("\\.png$", "", map_png)
status <- system2(
  pdftoppm,
  c("-f", "1", "-singlefile", "-png", "-r", "450", map_paths, map_stem)
)
if (status != 0L || !file.exists(map_png)) stop("Representative-map rasterization failed")
map_grob <- grid::rasterGrob(png::readPNG(map_png), interpolate = TRUE)
unlink(map_png)
p5f_map <- ggplot() +
  annotation_custom(map_grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  theme_void()
p5f_body <- arrangeGrob(p5f_map, p5f_null, nrow = 1,
                        widths = c(0.72, 1.28))
p5f <- stack_with_title(
  p5f_body,
  "Disease association and tissue organization are different claims",
  "Representative maps show abundance-adjusted scores; matched-null statistics quantify program-specific organization, not lineage origin or a spatial disease contrast."
)
save_vector(p5f, file.path(PANELS, "fig5f_spatial_maps_and_matched_null.pdf"), 7.1, 3.0)

# 6C: three assay-native evidence ribbons. Point color is state, never a score.
# KEY MESSAGE: Recurring examples retain native evidence states and route alternatives to experiments.
d6 <- fread(file.path(SOURCE, "fig6c_evidence_nodes.tsv"))
d6[, state_class := fifelse(grepl("supported|retained", state), "supported",
                     fifelse(grepl("ambient_sensitive|attenuation", state), "sensitive",
                     fifelse(grepl("descriptive|directional|mixed", state), "descriptive", "indeterminate")))]
unit_map <- c(
  "Inherited shared signal" = "PP.H4 (pending)",
  "Five-cohort RNA remodeling" = "cohort log2FC",
  "Cross-sectional histologic stage" = "adjacent-stage log2FC",
  "Liver protein decrease" = "protein log2FC",
  "Adjusted protein covariation" = "partial Spearman rho",
  "Donor-level disease association" = "Hotspot score / stage",
  "Ambient recalibration" = "Hotspot score / stage",
  "Same-atlas localization" = "paired standardized score",
  "Bulk tissue-state transport" = "member-gene log2FC",
  "Matched-null tissue organization" = "residual Moran I"
)
d6[, unit_short := unname(unit_map[node])]
d6[, label := vapply(
  paste(node, display, paste0("[", unit_short, "]"), sep = "\n"),
  function(x) paste(strwrap(x, width = 24), collapse = "\n"), character(1)
)]
experiments <- unique(d6[, .(example, next_experiment)])[, .(
  n_experiments = .N,
  experiment_label = if (.N == 1L) next_experiment[1] else
    paste0(.N, " unresolved alternatives route to ", .N,
           " discriminating experiments in panel 6D")
), by = example]
experiments[, `:=`(
  node_order = 6,
  node = "Discriminating experiment",
  state_class = "indeterminate",
  label = vapply(paste0("Next: ", experiment_label),
                 function(x) paste(strwrap(x, width = 27), collapse = "\n"), character(1))
)]
examples <- c("GNMT–MAT1A–CYP2C19", "ECM/IGFBP7", "Ductular-injury/BICC1")
d6[, example := factor(example, levels = rev(examples))]
experiments[, example := factor(example, levels = rev(examples))]
p6c <- ggplot(d6, aes(node_order, example)) +
  geom_segment(data = unique(d6[, .(example)]),
               aes(x = 1, xend = 6, y = example, yend = example),
               inherit.aes = FALSE, linewidth = 0.45, colour = "#BDBDBD") +
  geom_point(aes(colour = state_class, shape = state_class), size = 2.5) +
  geom_point(data = experiments, aes(node_order, example), inherit.aes = FALSE,
             shape = 21, size = 2.7, stroke = 0.5, fill = "white", colour = "black") +
  geom_text(aes(label = label), family = "Helvetica", size = GEOM_TEXT_6PT,
            vjust = -0.35, lineheight = 0.9, colour = "black") +
  geom_text(data = experiments, aes(node_order, example, label = label), inherit.aes = FALSE,
            family = "Helvetica", size = GEOM_TEXT_6PT, vjust = 1.3,
            lineheight = 0.9, colour = "black") +
  scale_colour_manual(values = state_colors) +
  scale_shape_manual(values = c(supported = 16, sensitive = 17,
                                descriptive = 15, indeterminate = 1)) +
  scale_x_continuous(breaks = 1:6, labels = NULL, limits = c(0.7, 6.3), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0.55, 3.5), clip = "off") +
  labs(
    title = "Recurring examples retain assay-native evidence states and route uncertainty to an experiment",
    x = NULL, y = NULL, colour = "Assay-native state", shape = "Assay-native state",
    caption = "No combined score, vote count, or rank is calculated. Full effect units, provenance, state reasons, and alternatives are in the source fingerprint table."
  ) +
  theme_classic(base_family = "Helvetica", base_size = 6) +
  theme(
    text = element_text(size = 6, family = "Helvetica", face = "plain", colour = "black"),
    axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6, colour = "black"),
    axis.title = element_blank(), legend.position = "bottom",
    legend.text = element_text(size = 6), legend.title = element_text(size = 6),
    plot.title = element_text(size = 6, face = "plain"),
    plot.caption = element_text(size = 6, face = "plain", hjust = 0),
    plot.margin = margin(14, 8, 10, 8)
  )
save_vector(p6c, file.path(PANELS, "fig6c_evidence_ribbons.pdf"), 12.0, 5.1)

# Composite visual proofs are review aids, not assembly artifacts.
proof4_body <- arrangeGrob(
  p4a, p4c, arrangeGrob(p4d_main, p4e, nrow = 1), p4f,
  ncol = 1, heights = c(0.75, 1, 1.25, 1.35)
)
proof4 <- stack_with_title(
  proof4_body, "Figure 4 aging-inspired rebuild: composite visual proof"
)
save_vector(proof4, file.path(PROOFS, "figure4_composite_proof.pdf"), 12.0, 13.0)
proof56_body <- arrangeGrob(p5f, p6c, ncol = 1, heights = c(1, 1.8))
proof56 <- stack_with_title(
  proof56_body, "Figures 5–6 aging-inspired rebuild: composite visual proof"
)
save_vector(proof56, file.path(PROOFS, "figure5_6_composite_proof.pdf"), 12.0, 8.5)

writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
cat("Rendered 7 individual panels and 2 composite proofs in", OUT, "\n")
