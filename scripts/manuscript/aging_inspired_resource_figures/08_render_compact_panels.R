#!/usr/bin/env Rscript
# KEY MESSAGE: Fixed donor-level programs can be stress-tested for ambient
# sensitivity, localized descriptively, and transported to tissue without
# collapsing assay-native effects into one score.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
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

compact_theme <- theme_masld_compact() +
  theme(
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    plot.caption = element_blank(),
    plot.margin = margin(2, 2, 2, 2),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0)
  )

program_colors <- c(
  "ECM/IGFBP7" = masld_colors$up,
  "Ductular-injury/BICC1" = "#D18B00"
)
state_colors <- c(
  supported = "#00796B",
  descriptive = "#4C78A8",
  sensitive = "#D18B00",
  indeterminate = "#9E9E9E"
)
lineage_colors <- c(
  hepatocytes = ct_palette[["Hepatocytes"]],
  fibroblasts = ct_palette[["Fibroblasts"]],
  macrophages = ct_palette[["Macrophages"]],
  cholangiocytes = ct_palette[["Cholangiocytes"]],
  tcells = ct_palette[["T cells"]]
)

pdf_raster_grob <- function(path, dpi = 450, crop = TRUE) {
  if (!file.exists(path)) stop("Missing PDF: ", path)
  pdftoppm <- Sys.which("pdftoppm")
  if (!nzchar(pdftoppm)) stop("pdftoppm is required")
  png_path <- tempfile("pdf_panel_", fileext = ".png")
  stem <- sub("\\.png$", "", png_path)
  status <- system2(pdftoppm, c("-f", "1", "-singlefile", "-png", "-r",
                                as.character(dpi), path, stem))
  if (status != 0L || !file.exists(png_path)) stop("PDF rasterization failed")
  img <- png::readPNG(png_path)
  unlink(png_path)
  if (crop) {
    rgb <- img[, , seq_len(min(3L, dim(img)[3L])), drop = FALSE]
    nonwhite <- apply(rgb, c(1, 2), min) < 0.985
    rows <- which(rowSums(nonwhite) > 0)
    cols <- which(colSums(nonwhite) > 0)
    if (length(rows) && length(cols)) {
      pad <- max(2L, round(min(dim(img)[1:2]) * 0.006))
      r1 <- max(1L, min(rows) - pad)
      r2 <- min(dim(img)[1], max(rows) + pad)
      c1 <- max(1L, min(cols) - pad)
      c2 <- min(dim(img)[2], max(cols) + pad)
      img <- img[r1:r2, c1:c2, , drop = FALSE]
    }
  }
  grid::rasterGrob(img, interpolate = TRUE)
}

grob_panel <- function(grob) {
  ggplot() +
    annotation_custom(grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    theme_void()
}

# 4A. Five compact visual steps; the caption carries the model details.
n4a <- data.table(
  x = 1:5,
  glyph = c("117", "β", "Δ", "cell", "tissue"),
  label = c("fixed programs", "donor effects", "ambient shift",
            "paired lineage", "transport"),
  colour = c("#4C78A8", masld_colors$up, "#D18B00", "#00796B", "#7B1FA2")
)
e4a <- data.table(x = 1:4, xend = 2:5)
p4a <- ggplot() +
  geom_segment(
    data = e4a, aes(x = x + 0.22, xend = xend - 0.22, y = 0.62, yend = 0.62),
    linewidth = 0.35, colour = "#9E9E9E",
    arrow = grid::arrow(length = grid::unit(1.2, "mm"), type = "closed")
  ) +
  geom_point(data = n4a, aes(x, 0.62, fill = colour), shape = 21, size = 10,
             stroke = 0.5, colour = "black") +
  geom_text(data = n4a, aes(x, 0.62, label = glyph), size = GEOM_TEXT_6PT,
            family = "Helvetica", colour = "black") +
  geom_text(data = n4a, aes(x, 0.22, label = label), size = GEOM_TEXT_6PT,
            family = "Helvetica", colour = "black") +
  scale_fill_identity() +
  coord_cartesian(xlim = c(0.65, 5.35), ylim = c(0.02, 0.95), clip = "off") +
  theme_void(base_family = "Helvetica", base_size = 6) +
  theme(text = element_text(size = 6, family = "Helvetica", face = "plain"),
        plot.margin = margin(2, 2, 2, 2))
save_vector(p4a, file.path(PANELS, "fig4a_analysis_logic.pdf"), 3.0, 1.25)

# 4C. Complete disease landscape as a compact donor-effect forest.
d4c <- fread(file.path(SOURCE, "fig4c_all117_disease_skyline.tsv"))
ct_order <- c("hepatocytes", "fibroblasts", "macrophages", "cholangiocytes", "tcells")
ct_labels <- c(hepatocytes = "Hepatocyte", fibroblasts = "Fibroblast",
               macrophages = "Macrophage", cholangiocytes = "Cholangiocyte",
               tcells = "T cell")
d4c[, cell_type := factor(cell_type, levels = ct_order)]
d4c[, y := as.numeric(cell_type) + seq(-0.28, 0.28, length.out = .N), by = cell_type]
d4c[, role := fifelse(grepl("^hero_", display_role), "hero",
               fifelse(!is.na(display_role) & nzchar(display_role), "control", "background"))]
hero4c <- d4c[role == "hero"]
control4c <- d4c[role == "control"]
labels4c <- rbindlist(list(hero4c, control4c), use.names = TRUE, fill = TRUE)
labels4c[, short_label := fifelse(grepl("IGFBP7", direct_label), "IGFBP7",
                           fifelse(grepl("BICC1", direct_label), "BICC1",
                                   sub(" control$", "", direct_label)))]
p4c <- ggplot(d4c, aes(disease_beta, y)) +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "#BDBDBD") +
  geom_segment(aes(x = raw_ci_low, xend = raw_ci_high, yend = y, colour = cell_type),
               linewidth = 0.25, alpha = 0.48) +
  geom_point(aes(fill = cell_type), shape = 21, size = 1.15,
             stroke = 0.15, colour = "white") +
  geom_point(data = control4c, shape = 1, size = 2.15, stroke = 0.45,
             colour = "black") +
  geom_point(data = hero4c, shape = 21, size = 2.4, stroke = 0.55,
             fill = "white", colour = "black") +
  geom_text_repel(
    data = labels4c, aes(label = short_label),
    size = GEOM_TEXT_6PT, family = "Helvetica", colour = "black",
    seed = 20260815, direction = "both", force = 1.5,
    box.padding = 0.15, point.padding = 0.15, min.segment.length = 0,
    segment.size = 0.2, segment.color = "#6F6F6F", max.overlaps = Inf,
    ylim = c(0.64, 1.38)
  ) +
  scale_colour_manual(values = lineage_colors, guide = "none") +
  scale_fill_manual(values = lineage_colors, guide = "none") +
  scale_y_continuous(breaks = seq_along(ct_order), labels = ct_labels[ct_order]) +
  labs(x = "Stage effect (Hotspot score / ordinal stage)", y = NULL) +
  coord_cartesian(clip = "off") + compact_theme +
  theme(axis.text.y = element_text(size = 6), plot.margin = margin(2, 8, 2, 2))
save_vector(p4c, file.path(PANELS, "fig4c_all117_disease_skyline.pdf"), 7.1, 2.45)

# 4D. Raw-to-corrected transport, with only the two examples and four controls
# carrying arrows or labels.
d4d <- fread(file.path(SOURCE, "fig4d_ambient_effect_transport.tsv"))
d4d[, role := fifelse(grepl("^hero_", display_role), "hero",
               fifelse(!is.na(display_role) & nzchar(display_role), "control", "background"))]
finite4d <- d4d[is.finite(corrected_beta)]
nt4d <- d4d[!is.finite(corrected_beta)]
hero4d <- finite4d[role == "hero"]
control4d <- finite4d[role == "control"]
lims4d <- range(c(finite4d$raw_beta, finite4d$corrected_beta), finite = TRUE)
lims4d <- lims4d + c(-0.05, 0.05) * diff(lims4d)
p4d <- ggplot(finite4d, aes(raw_beta, corrected_beta)) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.3, colour = "#9E9E9E") +
  geom_point(shape = 21, size = 1.2, stroke = 0.2, fill = "white",
             colour = "#9E9E9E") +
  geom_segment(
    data = control4d,
    aes(x = raw_beta, y = raw_beta, xend = raw_beta, yend = corrected_beta),
    inherit.aes = FALSE, linewidth = 0.35, colour = "#595959",
    arrow = grid::arrow(length = grid::unit(0.9, "mm"), type = "closed")
  ) +
  geom_segment(
    data = hero4d,
    aes(x = raw_beta, y = raw_beta, xend = raw_beta, yend = corrected_beta,
        colour = direct_label),
    inherit.aes = FALSE, linewidth = 0.7,
    arrow = grid::arrow(length = grid::unit(1.1, "mm"), type = "closed")
  ) +
  geom_point(data = hero4d, aes(colour = direct_label), size = 2.1) +
  geom_text(data = hero4d, aes(label = ifelse(grepl("IGFBP7", direct_label), "IGFBP7", "BICC1")),
            size = GEOM_TEXT_6PT, family = "Helvetica", hjust = -0.12,
            vjust = 0.25, colour = "black") +
  geom_text(data = control4d,
            aes(label = sub(" control$", "", direct_label)),
            size = GEOM_TEXT_6PT, family = "Helvetica", hjust = -0.12,
            vjust = 0.25, colour = "black", check_overlap = TRUE) +
  annotate("text", x = lims4d[1], y = lims4d[2], hjust = 0, vjust = 1,
           label = "5 corrected\n2 passthrough*\n*GSE189600, Liver Atlas",
           size = GEOM_TEXT_6PT,
           family = "Helvetica", colour = "black", lineheight = 0.9) +
  annotate("text", x = lims4d[1], y = lims4d[1], hjust = 0, vjust = -0.25,
           label = "11 T-cell programs NT", size = GEOM_TEXT_6PT,
           family = "Helvetica", colour = "black") +
  scale_colour_manual(values = program_colors, guide = "none") +
  coord_equal(xlim = lims4d, ylim = lims4d, clip = "off") +
  labs(x = "Raw stage effect", y = "Corrected stage effect") + compact_theme +
  theme(plot.margin = margin(3, 9, 3, 3))
save_vector(p4d, file.path(PANELS, "fig4d_ambient_effect_transport.pdf"), 2.25, 2.45)

# 4E. Two stacked mini-forests; target lineage is the only saturated mark.
d4e <- fread(file.path(SOURCE, "fig4e_same_atlas_lineage_specificity.tsv"))
lineage_order <- rev(c("Hepatocytes", "Fibroblasts", "Cholangiocytes",
                       "Endothelial cells", "Macrophages", "T cells"))
d4e[, comparison_lineage := factor(comparison_lineage, levels = lineage_order)]
d4e[, target := comparison_lineage == target_lineage]
d4e[, reference := comparison_lineage == "Hepatocytes"]
p4e <- ggplot(d4e, aes(beta, comparison_lineage)) +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "#BDBDBD") +
  geom_segment(aes(x = ci_low, xend = ci_high, yend = comparison_lineage),
               linewidth = 0.4, colour = "#6F6F6F") +
  geom_point(data = d4e[target == FALSE & reference == FALSE], shape = 21,
             size = 1.65, stroke = 0.3, fill = "white", colour = "#6F6F6F") +
  geom_point(data = d4e[reference == TRUE], shape = 21, size = 1.65,
             stroke = 0.3, fill = "#9E9E9E", colour = "black") +
  geom_point(data = d4e[target == TRUE], aes(fill = program_name), shape = 21,
             size = 2.25, stroke = 0.4, colour = "black") +
  facet_wrap(~program_name, ncol = 1, scales = "free_x") +
  scale_fill_manual(values = program_colors, guide = "none") +
  labs(x = "Lineage − hepatocyte score", y = NULL) + compact_theme +
  theme(strip.text = element_text(size = 6), axis.text.y = element_text(size = 6),
        panel.spacing.y = grid::unit(1.5, "mm"))
save_vector(p4e, file.path(PANELS, "fig4e_same_atlas_lineage_specificity.pdf"), 2.75, 2.45)

# 4F. Single-cell correction and bulk transport remain in separate aligned axes.
sc4f <- fread(file.path(SOURCE, "fig4f_singlecell_estimates.tsv"))
sc4f[, estimate_type := factor(estimate_type, levels = c("raw", "corrected"))]
sc4f[, program_short := fifelse(grepl("IGFBP7", program_name), "IGFBP7", "BICC1")]
wide4f <- dcast(sc4f, program_name ~ estimate_type, value.var = "beta")
wide4f[, program_short := fifelse(grepl("IGFBP7", program_name), "IGFBP7", "BICC1")]
p4f_sc <- ggplot(sc4f, aes(beta, program_short, colour = program_name)) +
  geom_segment(
    data = wide4f,
    aes(x = raw, xend = corrected, y = program_short, yend = program_short,
        colour = program_name),
    inherit.aes = FALSE, linewidth = 0.65,
    arrow = grid::arrow(length = grid::unit(1.1, "mm"), type = "closed")
  ) +
  geom_point(data = sc4f[estimate_type == "raw"], size = 2.0) +
  geom_point(data = sc4f[estimate_type == "corrected"], shape = 21, size = 2.1,
             stroke = 0.45, fill = "white") +
  scale_colour_manual(values = program_colors, guide = "none") +
  labs(x = "Single-cell effect", y = NULL) + compact_theme +
  theme(axis.text.y = element_text(size = 6))
bulk4f <- fread(file.path(SOURCE, "fig4f_bulk_stage_transport.tsv"))
bulk4f[, stage_num := match(stage, c("F1", "F2", "F3", "F4"))]
lab4f <- bulk4f[stage == "F4"]
p4f_bulk <- ggplot(bulk4f, aes(stage_num, effect, colour = program_name,
                               group = program_name)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "#BDBDBD") +
  geom_line(linewidth = 0.65) + geom_point(size = 1.55) +
  geom_text(data = lab4f,
            aes(label = ifelse(grepl("IGFBP7", program_name), "IGFBP7", "BICC1")),
            size = GEOM_TEXT_6PT, family = "Helvetica", hjust = 1.12,
            colour = "black") +
  scale_colour_manual(values = program_colors, guide = "none") +
  scale_x_continuous(breaks = 1:4, labels = c("F1", "F2", "F3", "F4"),
                     limits = c(0.9, 4.7)) +
  labs(x = "Fibrosis stage vs F0", y = "Bulk member-gene log2FC") + compact_theme +
  theme(plot.margin = margin(2, 8, 2, 2))
p4f <- arrangeGrob(p4f_sc, p4f_bulk, nrow = 1, widths = c(1.05, 1.25))
save_vector(p4f, file.path(PANELS, "fig4f_tissue_state_transport.pdf"), 3.2, 2.45)

# 5F. Retained maps are cropped to their artwork and paired directly with the
# matched-null interval, without a second title, subtitle, caption, or legend.
d5f <- fread(file.path(SOURCE, "fig5f_spatial_map_and_matched_null.tsv"))
map_paths <- unique(d5f$map_panel_path)
if (length(map_paths) != 1L) stop("Expected one retained map panel")
p5f_map <- grob_panel(pdf_raster_grob(map_paths, dpi = 500, crop = TRUE))
d5f[, program_short := factor(program_short, levels = c("BICC1", "IGFBP7"))]
d5f[, dataset_label := factor(dataset_label, levels = rev(c("GSE", "Vu*")))]
d5f[, state := fifelse(within_source_call == "supported", "supported", "indeterminate")]
p5f_null <- ggplot(d5f, aes(residual_moran_i, dataset_label)) +
  geom_segment(aes(x = q050, xend = q950, yend = dataset_label),
               linewidth = 2.2, colour = "#BDBDBD") +
  geom_point(aes(x = q500), shape = 23, size = 1.65, fill = "white",
             colour = "black") +
  geom_segment(aes(x = q500, xend = residual_moran_i, yend = dataset_label),
               linewidth = 0.3, linetype = 2, colour = "#6F6F6F") +
  geom_point(aes(fill = state), shape = 21, size = 2.1, stroke = 0.35,
             colour = "black") +
  geom_text(aes(label = paste0("q=", formatC(primary_qvalue_spatial,
                                              format = "g", digits = 2))),
            size = GEOM_TEXT_6PT, family = "Helvetica", hjust = 0.5,
            vjust = -0.85,
            colour = "black") +
  facet_grid(. ~ program_short, scales = "free_x") +
  scale_fill_manual(values = state_colors[c("supported", "indeterminate")],
                    guide = "none") +
  labs(x = "Residual Moran I  (gray: matched-gene 5–95%)", y = NULL) +
  coord_cartesian(clip = "off") + compact_theme +
  theme(plot.margin = margin(2, 10, 2, 2), strip.text = element_text(size = 6))
p5f <- arrangeGrob(p5f_map, p5f_null, nrow = 1, widths = c(0.9, 1.25))
save_vector(p5f, file.path(PANELS, "fig5f_spatial_maps_and_matched_null.pdf"), 4.6, 2.55)

# 6C. Two compact five-node ribbons. Text is limited to assay and one native
# value; qualifications and full units remain in the exact source table.
d6 <- fread(file.path(SOURCE, "fig6c_evidence_nodes.tsv"))
d6[, block := fifelse(grepl("^GNMT", example), "genetic_state", "program")]
d6[, x := node_order]
d6[, y := fifelse(grepl("^GNMT", example), 3,
           fifelse(example == "ECM/IGFBP7", 2, 1))]
d6[, state_class := fifelse(
  state %in% c("supported_multicohort_remodeling", "supported_stage_association",
               "supported_selection_conditioned"), "supported",
  fifelse(grepl("attenuation", state), "sensitive",
  fifelse(grepl("descriptive|directional|mixed", state), "descriptive",
          "indeterminate")))]
d6[state == "descriptive_localization_supported", state_class := "descriptive"]
d6[example == "ECM/IGFBP7" & node_order == 5, state_class := "supported"]
d6[, assay_short := fifelse(block == "genetic_state",
  c("Genetics", "RNA", "Stage", "Protein", "Covariation")[node_order],
  c("Disease", "Ambient", "Lineage", "Bulk", "Spatial")[node_order])]
d6[, value_short := ""]
d6[example == "GNMT–MAT1A–CYP2C19", value_short := c(
  "pending", "5/5 · 4/5 · 5/5 ↓", "4/4 ↓", "3 proteins ↓", "ρ < 0"
)]
d6[example == "ECM/IGFBP7", value_short := c(
  "β 0.40", "β 0.13", "fib β 2.27", "F4 0.76", "supported"
)]
d6[example == "Ductular-injury/BICC1", value_short := c(
  "β 0.40", "β 0.39", "chol β 1.57", "F4 1.10", "indeterminate"
)]
d6[, label_y := y + fifelse(y == 2, 0.24, 0.22)]
line6 <- unique(d6[, .(y, block)])
p6c <- ggplot(d6, aes(x, y)) +
  geom_segment(data = line6, aes(x = 1, xend = 5, y = y, yend = y),
               inherit.aes = FALSE, linewidth = 0.45, colour = "#BDBDBD") +
  geom_point(data = d6[state_class == "indeterminate"], shape = 21, size = 2.5,
             stroke = 0.45, fill = "white", colour = "#9E9E9E") +
  geom_point(data = d6[state_class != "indeterminate"],
             aes(fill = state_class, shape = state_class), size = 2.5,
             stroke = 0.35, colour = "black") +
  geom_text(aes(y = label_y, label = assay_short), size = GEOM_TEXT_6PT,
            family = "Helvetica", colour = "black") +
  geom_text(aes(y = y - 0.22, label = value_short), size = GEOM_TEXT_6PT,
            family = "Helvetica", colour = "black") +
  annotate("text", x = 0.68, y = 3, label = "GNMT state",
           hjust = 1, size = GEOM_TEXT_6PT, family = "Helvetica") +
  annotate("text", x = 0.68, y = 2, label = "ECM/IGFBP7",
           hjust = 1, size = GEOM_TEXT_6PT, family = "Helvetica") +
  annotate("text", x = 0.68, y = 1, label = "BICC1",
           hjust = 1, size = GEOM_TEXT_6PT, family = "Helvetica") +
  scale_fill_manual(values = state_colors, guide = "none") +
  scale_shape_manual(values = c(supported = 21, descriptive = 22, sensitive = 24),
                     guide = "none") +
  coord_cartesian(xlim = c(0.65, 5.25), ylim = c(0.55, 3.45), clip = "off") +
  theme_void(base_family = "Helvetica", base_size = 6) +
  theme(text = element_text(size = 6, family = "Helvetica", face = "plain"),
        plot.margin = margin(4, 5, 4, 36))
save_vector(p6c, file.path(PANELS, "fig6c_evidence_ribbons.pdf"), 5.0, 2.55)

# 6D. Every unresolved alternative has exactly one dot in an experiment-class
# matrix; the long-form wording stays in the source table.
d6d <- fread(file.path(SOURCE, "fig6d_experiment_router.tsv"))
d6d[, route := fifelse(grepl("Allele-aware", next_experiment), "Allele edit",
                fifelse(grepl("Cell-type", next_experiment), "Cell perturb",
                fifelse(grepl("Prospective", next_experiment), "Longitudinal",
                fifelse(grepl("Independent liver proteomics", next_experiment), "Proteomics",
                fifelse(grepl("methionine-cycle", next_experiment), "Flux phenotype",
                fifelse(grepl("Fibroblast-restricted", next_experiment), "Fibroblast tissue",
                        "Cholangiocyte tissue"))))))]
d6d[, node_short := fifelse(node == "Inherited shared signal", "Signal",
                     fifelse(node == "Five-cohort RNA remodeling", "RNA",
                     fifelse(node == "Cross-sectional histologic stage", "Stage",
                     fifelse(node == "Liver protein decrease", "Protein",
                     fifelse(node == "Adjusted protein covariation", "Covariation",
                     fifelse(node == "Donor-level disease association", "Disease",
                     fifelse(node == "Ambient recalibration", "Ambient",
                     fifelse(node == "Same-atlas localization", "Lineage",
                     fifelse(node == "Bulk tissue-state transport", "Bulk", "Spatial")))))))))]
example_short <- c(
  "GNMT–MAT1A–CYP2C19" = "GNMT state",
  "ECM/IGFBP7" = "IGFBP7",
  "Ductular-injury/BICC1" = "BICC1"
)
d6d[, row_label := paste0(example_short[example], " · ", node_short)]
example_order <- c("GNMT–MAT1A–CYP2C19", "ECM/IGFBP7", "Ductular-injury/BICC1")
node_order6d <- c("Signal", "RNA", "Stage", "Protein", "Covariation",
                  "Disease", "Ambient", "Lineage", "Bulk", "Spatial")
d6d[, sort_group := match(example, example_order)]
d6d[, sort_node := match(node_short, node_order6d)]
row_levels <- d6d[order(sort_group, sort_node), unique(row_label)]
d6d[, row_label := factor(row_label, levels = rev(row_levels))]
route_levels <- c("Allele edit", "Cell perturb", "Longitudinal", "Proteomics",
                  "Flux phenotype", "Fibroblast tissue", "Cholangiocyte tissue")
d6d[, route := factor(route, levels = route_levels)]
route_colors <- c("GNMT–MAT1A–CYP2C19" = "#00796B",
                  "ECM/IGFBP7" = masld_colors$up,
                  "Ductular-injury/BICC1" = "#D18B00")
p6d <- ggplot(d6d, aes(route, row_label)) +
  geom_vline(xintercept = seq_along(route_levels), linewidth = 0.2,
             colour = "#E5E5E5") +
  geom_hline(yintercept = c(5.5, 10.5), linewidth = 0.3,
             colour = "#BDBDBD") +
  geom_point(aes(fill = example), shape = 21, size = 2.35,
             stroke = 0.35, colour = "black") +
  scale_fill_manual(values = route_colors, guide = "none") +
  labs(x = NULL, y = NULL) + compact_theme +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6), axis.ticks = element_blank(),
        axis.line = element_blank(), plot.margin = margin(2, 2, 2, 2))
save_vector(p6d, file.path(PANELS, "fig6d_experiment_router.pdf"), 5.0, 3.05)

# Compact review proofs use the intended manuscript-scale panel proportions.
umap_path <- file.path(
  ROOT, "figures/candidates/pi-figure-redesign-2026-08-13-v8/figure4/panels",
  "fig4b_scrna_umap_embeddable.pdf"
)
p4b <- grob_panel(pdf_raster_grob(umap_path, dpi = 450, crop = TRUE))
proof4_top <- arrangeGrob(p4a, p4b, nrow = 1, widths = c(3.0, 3.1))
proof4_bottom <- arrangeGrob(p4d, p4e, p4f, nrow = 1,
                            widths = c(2.25, 2.75, 3.2))
proof4 <- arrangeGrob(proof4_top, p4c, proof4_bottom, ncol = 1,
                      heights = c(1.35, 2.45, 2.45))
save_vector(proof4, file.path(PROOFS, "figure4_composite_proof.pdf"), 7.1, 6.4)
proof56 <- arrangeGrob(p5f, p6c, p6d, ncol = 1, heights = c(2.55, 2.55, 3.05))
save_vector(proof56, file.path(PROOFS, "figure5_6_composite_proof.pdf"), 5.2, 8.25)

writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
cat("Rendered 8 compact panels and 2 manuscript-scale proofs in", OUT, "\n")
