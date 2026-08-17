#!/usr/bin/env Rscript

# Figure 5E multimodal program summary.
#
# Selection is data-driven: retain a frozen Figure 3 program when it has at
# least one credible-set variant linked to a lineage-matched open promoter, at
# least one robust MASLD-only protein-histology association, or robust spatial
# organization in both cohorts. Each column is summarized at the program level:
# distinct open-promoter program genes at fine-mapped liver-enzyme loci, a
# bootstrap interval for protein association with an outcome-only histologic-
# burden axis, and a nested shared-core/range spatial bar. All tests and BH
# adjustments remain those from the full frozen-program universe; display
# filtering never triggers re-testing. No cross-modality score or ranking is
# constructed.
# KEY MESSAGE: Frozen programs occupy distinct regulatory, histologic-protein, and replicated spatial contexts without a cross-modality score.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})
grDevices::pdf.options(useDingbats = FALSE)

BASE <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

FS <- 6
FAM <- "Helvetica"
INK <- "#222222"
NEG <- "#2166AC"
POS <- "#B2182B"
MID <- "#FFFFFF"
GOLD <- "#C79A3E"
SPATIAL_SHARED <- "#3D86A8"
SPATIAL_RANGE <- "#C8E0EB"
CT_COLS <- c(
  hepatocytes = unname(ct_palette["Hepatocytes"]),
  fibroblasts = unname(ct_palette["Fibroblasts"]),
  macrophages = unname(ct_palette["Macrophages"]),
  cholangiocytes = unname(ct_palette["Cholangiocytes"])
)
CT_ABBR <- c(
  hepatocytes = "Hep", fibroblasts = "Fib",
  macrophages = "Mac", cholangiocytes = "Chol"
)

burden_file <- file.path(
  PROGRAM_CONTEXT_DIR, "proteomics", "module_protein_histology_burden_masld.tsv"
)
wide_file <- file.path(PROGRAM_CONTEXT_DIR, "program_context_wide.tsv")
spatial_file <- file.path(
  PROGRAM_CONTEXT_DIR, "spatial", "spatial_program_results.tsv"
)
inputs <- c(burden_file, wide_file, spatial_file)
if (any(!file.exists(inputs))) {
  stop("Missing candidate input(s): ", paste(inputs[!file.exists(inputs)], collapse = ", "))
}

burden <- fread(burden_file)
wide <- fread(wide_file)
spatial <- fread(spatial_file)
stopifnot(
  nrow(wide) == 22L,
  !anyDuplicated(wide$program_id),
  nrow(burden) == 22L,
  !anyDuplicated(burden$program_id),
  nrow(spatial) == 22L * 2L
)
stopifnot(
  all(grepl("rank-residual partial Spearman", burden$inference_method, fixed = TRUE)),
  all(burden[testable == TRUE, n_permutations] >= 999L)
)

# Selection requires a concrete, assay-native result rather than a hand-picked
# program list. A spatial program must be robust in both cohorts to qualify.
burden_any <- burden[, .(
  histology_burden_robust = robust %in% TRUE
), by = program_id]
spatial_any <- spatial[, .(
  spatial_shared = all(testable == TRUE & robust == TRUE)
), by = program_id]
selected <- merge(wide, burden_any, by = "program_id", all.x = TRUE, sort = FALSE)
selected <- merge(selected, spatial_any, by = "program_id", all.x = TRUE, sort = FALSE)
stopifnot(all(fcoalesce(as.integer(selected$gwas_atac_direct_genes), 0L) == 0L))
selected[, open_promoter_genes :=
  fcoalesce(as.integer(gwas_atac_enzyme_genes), 0L)]
selected[is.na(histology_burden_robust), histology_burden_robust := FALSE]
selected[is.na(spatial_shared), spatial_shared := FALSE]
selected <- selected[
  open_promoter_genes > 0L | histology_burden_robust == TRUE |
    spatial_shared == TRUE
]
setorder(selected, display_order)
selected[, y := .N:1L]
stopifnot(nrow(selected) > 0L, uniqueN(selected$program_id) == nrow(selected))

selection_sidecar <- selected[, .(
  program_id, display_order, cell_type, module, program_name, stability_fail,
  open_promoter_genes, histology_burden_robust, spatial_shared,
  displayed_n = .N,
  frozen_universe_n = nrow(wide),
  selection_rule = paste0(
    "open_promoter_gene OR robust_histology_burden OR ",
    "robust_spatial_in_both"
  )
)]
fwrite(
  selection_sidecar,
  file.path(FIG5_CONTEXT_DIR, "panels", "data", "fig5e_multimodal_program_summary_selection.tsv"),
  sep = "\t",
  quote = FALSE
)

selected[, short_name := sub(" \\([^()]+\\)$", "", program_name)]
selected[short_name == "Fatty-acid / peroxisomal metab",
         short_name := "FA / peroxisomal metabolism"]
selected[short_name == "Glucocorticoid resp",
         short_name := "Glucocorticoid response"]
selected[short_name == "Amino-acid metabolism",
         short_name := "AA metabolism"]
selected[short_name == "Hepatocyte stress",
         short_name := "Stress"]
selected[short_name == "NRF2 antioxidant",
         short_name := "NRF2"]
selected[short_name == "Insulin response",
         short_name := "Insulin"]
selected[short_name == "Glucocorticoid response",
         short_name := "Glucocorticoid"]
selected[short_name == "Glucocorticoid / stress",
         short_name := "Gluc/stress"]
selected[short_name == "Endocytosis / TLR2",
         short_name := "Endocytosis"]
selected[, row_label := paste0(module, "  ", short_name)]
selected[, label_col := INK]

n <- nrow(selected)
YMIN <- -1.45
YMAX <- n + 2.05
groups <- selected[, .(
  ymin = min(y) - 0.46,
  ymax = max(y) + 0.46,
  cy = mean(y),
  colour = CT_COLS[.BY[[1]]],
  label = CT_ABBR[.BY[[1]]]
), by = cell_type]
breaks_y <- selected[, which(cell_type[-.N] != cell_type[-1L])]
breaks_y <- if (length(breaks_y)) selected$y[breaks_y] - 0.5 else numeric()

base_void <- theme_void(base_family = FAM, base_size = FS) +
  theme(
    text = element_text(family = FAM, size = FS, face = "plain", colour = INK),
    plot.margin = margin(1, 1, 1, 1),
    legend.position = "none"
  )

separator_layers <- function(xmin, xmax) {
  if (!length(breaks_y)) return(list())
  list(geom_segment(
    data = data.table(y = breaks_y),
    aes(x = xmin, xend = xmax, y = y, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.22,
    linetype = "dashed",
    colour = "#9A9A9A"
  ))
}

# ── Shared program labels + cell-type brackets ──────────────────────────────
label_brackets <- rbind(
  groups[, .(x = 0.48, xend = 0.48, y = ymin + 0.12, yend = ymax - 0.12, colour)],
  groups[, .(x = 0.48, xend = 0.57, y = ymin + 0.12, yend = ymin + 0.12, colour)],
  groups[, .(x = 0.48, xend = 0.57, y = ymax - 0.12, yend = ymax - 0.12, colour)]
)

p_labels <- ggplot() +
  geom_segment(
    data = label_brackets,
    aes(x = x, xend = xend, y = y, yend = yend),
    colour = label_brackets$colour,
    linewidth = 0.35 / .pt,
    lineend = "round"
  ) +
  geom_text(
    data = groups,
    aes(x = 0.42, y = cy, label = label),
    hjust = 1,
    size = FS / .pt,
    family = FAM,
    colour = INK
  ) +
  geom_text(
    data = selected,
    aes(x = 0.64, y = y, label = row_label),
    hjust = 0,
    size = FS / .pt,
    family = FAM,
    colour = selected$label_col
  ) +
  separator_layers(0.03, 2.62) +
  coord_cartesian(xlim = c(0.02, 2.62), ylim = c(YMIN, YMAX), expand = FALSE) +
  base_void

# ── Fine-mapped liver-enzyme loci: linked program-gene count ────────────────
# Variant identities and traits remain in Panel 4D. Here the unit matches the
# row: the number of distinct program genes whose promoter contains at least
# one fine-mapped liver-enzyme variant in a lineage-matched open ATAC peak.
genetic <- selected[, .(
  program_id, y, n_linked_genes = open_promoter_genes
)]
GEN_MAX <- max(3L, max(genetic$n_linked_genes, na.rm = TRUE))
x_ticks_gen <- data.table(x = 0:GEN_MAX, label = as.character(0:GEN_MAX))

p_genetic <- ggplot() +
  geom_rect(
    data = genetic,
    aes(xmin = 0, xmax = GEN_MAX, ymin = y - 0.16, ymax = y + 0.16),
    fill = "#F1F1F1",
    colour = NA
  ) +
  geom_rect(
    data = genetic[n_linked_genes > 0L],
    aes(xmin = 0, xmax = n_linked_genes, ymin = y - 0.16, ymax = y + 0.16),
    fill = GOLD,
    colour = NA
  ) +
  annotate(
    "text", x = 0, y = n + 1.63,
    label = "Open-promoter n",
    hjust = 0,
    size = FS / .pt,
    family = FAM,
    colour = INK
  ) +
  annotate("segment", x = 0, xend = GEN_MAX, y = -0.18, yend = -0.18,
           linewidth = 0.25, colour = INK) +
  geom_segment(
    data = x_ticks_gen,
    aes(x = x, xend = x, y = -0.18, yend = -0.32),
    linewidth = 0.25,
    colour = INK
  ) +
  geom_text(
    data = x_ticks_gen,
    aes(x = x, y = -0.58, label = label),
    size = FS / .pt,
    family = FAM,
    colour = INK
  ) +
  separator_layers(-0.08, GEN_MAX + 0.08) +
  coord_cartesian(
    xlim = c(-0.08, GEN_MAX + 0.08), ylim = c(YMIN, YMAX), expand = FALSE
  ) +
  base_void

# ── Protein score ↔ histologic-burden forest intervals ─────────────────────
# A single program-level estimate replaces five feature wedges. The burden
# axis is PC1 of steatosis, ballooning, inflammation, and fibrosis; NAS is
# excluded because it is a composite of the first three. Intervals are patient-
# bootstrap 95% CIs with resampling stratified by acquisition batch.
hist_plot <- merge(
  burden,
  selected[, .(program_id, y)],
  by = "program_id",
  all.y = TRUE,
  sort = FALSE
)
hist_plot[, estimate_col := fifelse(rho >= 0, POS, NEG)]
HIST_LIM <- max(
  0.6,
  ceiling(max(abs(c(hist_plot$ci_low, hist_plot$ci_high)), na.rm = TRUE) * 10) / 10
)
HIST_LIM <- min(1, HIST_LIM)
hist_ticks <- data.table(x = pretty(c(-HIST_LIM, HIST_LIM), n = 4L))
hist_ticks <- hist_ticks[x >= -HIST_LIM & x <= HIST_LIM]

p_hist <- ggplot() +
  annotate(
    "segment", x = 0, xend = 0, y = 0.48, yend = n + 0.48,
    linewidth = 0.2,
    colour = "#B8B8B8"
  ) +
  geom_segment(
    data = hist_plot[testable == TRUE],
    aes(x = ci_low, xend = ci_high, y = y, yend = y),
    colour = alpha(hist_plot[testable == TRUE]$estimate_col, 0.72),
    linewidth = 0.48,
    lineend = "round"
  ) +
  geom_segment(
    data = hist_plot[testable == TRUE],
    aes(x = rho, xend = rho, y = y - 0.17, yend = y + 0.17),
    colour = hist_plot[testable == TRUE]$estimate_col,
    linewidth = 0.65,
    lineend = "round"
  ) +
  annotate(
    "text", x = 0, y = n + 1.63,
    label = "Histology partial ρ",
    hjust = 0.5,
    size = FS / .pt,
    family = FAM,
    colour = INK
  ) +
  annotate("segment", x = -HIST_LIM, xend = HIST_LIM,
           y = -0.18, yend = -0.18, linewidth = 0.25, colour = INK) +
  geom_segment(
    data = hist_ticks,
    aes(x = x, xend = x, y = -0.18, yend = -0.32),
    linewidth = 0.25,
    colour = INK
  ) +
  geom_text(
    data = hist_ticks,
    aes(x = x, y = -0.58, label = sprintf("%g", x)),
    size = FS / .pt,
    family = FAM,
    colour = INK
  ) +
  separator_layers(-HIST_LIM, HIST_LIM) +
  coord_cartesian(
    xlim = c(-HIST_LIM, HIST_LIM), ylim = c(YMIN, YMAX), expand = FALSE
  ) +
  base_void

# ── Spatial excess over matched genes: shared core + cohort range ───────────
# The pale segment is the interval between the two cohort estimates. A dark
# segment from zero to the smaller common effect is drawn only when both cohort
# tests pass their native robustness gate. Thus no average or maximum can hide
# discordance, and the two studies occupy one aligned program-level track.
spatial_long <- copy(spatial[program_id %in% selected$program_id])
spatial_long[, excess_moran := residual_moran_i - residual_null_mean]
spatial_gse <- spatial_long[dataset == "GSE192741", .(
  program_id, gse_effect = excess_moran,
  gse_testable = testable, gse_robust = robust
)]
spatial_vu <- spatial_long[dataset == "Vu_et_al_2025", .(
  program_id, vu_effect = excess_moran,
  vu_testable = testable, vu_robust = robust
)]
spatial_plot <- merge(spatial_gse, spatial_vu, by = "program_id", all = TRUE)
spatial_plot <- merge(
  selected[, .(program_id, y)], spatial_plot,
  by = "program_id", all.x = TRUE, sort = FALSE
)
spatial_plot[, `:=`(
  range_low = pmin(gse_effect, vu_effect),
  range_high = pmax(gse_effect, vu_effect),
  shared_robust = gse_testable == TRUE & vu_testable == TRUE &
    gse_robust == TRUE & vu_robust == TRUE &
    sign(gse_effect) == sign(vu_effect),
  shared_effect = fifelse(
    gse_effect > 0 & vu_effect > 0,
    pmin(gse_effect, vu_effect),
    fifelse(gse_effect < 0 & vu_effect < 0, pmax(gse_effect, vu_effect), 0)
  )
)]
spatial_plot[, `:=`(
  shared_low = pmin(0, shared_effect),
  shared_high = pmax(0, shared_effect)
)]

SP_MIN <- min(
  0,
  floor(min(spatial_plot$range_low, na.rm = TRUE) * 100) / 100
)
SP_MAX <- max(
  0.01,
  ceiling(max(spatial_plot$range_high, na.rm = TRUE) * 100) / 100
)
spatial_ticks <- data.table(x = pretty(c(SP_MIN, SP_MAX), n = 4L))
spatial_ticks <- spatial_ticks[x >= SP_MIN & x <= SP_MAX]
key_span <- SP_MAX - SP_MIN

p_spatial <- ggplot() +
  annotate(
    "segment", x = 0, xend = 0, y = 0.48, yend = n + 0.48,
    linewidth = 0.2,
    colour = "#B8B8B8"
  ) +
  geom_rect(
    data = spatial_plot[gse_testable == TRUE & vu_testable == TRUE],
    aes(xmin = range_low, xmax = range_high,
        ymin = y - 0.17, ymax = y + 0.17),
    fill = SPATIAL_RANGE,
    colour = NA
  ) +
  geom_rect(
    data = spatial_plot[shared_robust == TRUE],
    aes(xmin = shared_low, xmax = shared_high,
        ymin = y - 0.17, ymax = y + 0.17),
    fill = SPATIAL_SHARED,
    colour = NA
  ) +
  annotate("text", x = SP_MIN, y = n + 1.63,
           label = "Spatial Moran excess",
           hjust = 0, size = FS / .pt, family = FAM, colour = INK) +
  annotate(
    "rect",
    xmin = c(SP_MIN, SP_MIN + 0.48 * key_span),
    xmax = c(SP_MIN + 0.09 * key_span, SP_MIN + 0.57 * key_span),
    ymin = n + 0.64, ymax = n + 0.86,
    fill = c(SPATIAL_SHARED, SPATIAL_RANGE), colour = NA
  ) +
  annotate(
    "text",
    x = c(SP_MIN + 0.12 * key_span, SP_MIN + 0.60 * key_span),
    y = n + 0.75,
    label = c("shared", "range"),
    hjust = 0,
    size = FS / .pt,
    family = FAM,
    colour = INK
  ) +
  annotate("segment", x = SP_MIN, xend = SP_MAX,
           y = -0.18, yend = -0.18, linewidth = 0.25, colour = INK) +
  geom_segment(
    data = spatial_ticks,
    aes(x = x, xend = x, y = -0.18, yend = -0.32),
    linewidth = 0.25,
    colour = INK
  ) +
  geom_text(
    data = spatial_ticks,
    aes(x = x, y = -0.58, label = sprintf("%g", x)),
    size = FS / .pt,
    family = FAM,
    colour = INK
  ) +
  separator_layers(SP_MIN, SP_MAX) +
  coord_cartesian(xlim = c(SP_MIN, SP_MAX), ylim = c(YMIN, YMAX), expand = FALSE) +
  base_void

candidate <- p_labels + p_genetic + p_hist + p_spatial +
  plot_layout(widths = c(1.52, 0.94, 1.28, 1.30))

out <- file.path(FIG5_CONTEXT_DIR, "panels", "fig5e_multimodal_program_summary.pdf")
ggsave(
  out,
  candidate,
  # Compact full-panel width: redundant cell-type prefixes are carried by the
  # brackets, while shorter track labels preserve the fixed 6 pt typography.
  width = 3.63,
  height = 1.88,
  units = "in",
  device = grDevices::cairo_pdf
)
message(
  "[fig5e summary] saved: ", out,
  "; selected ", n,
  "/22 programs (open promoter OR robust histology OR shared spatial)"
)
message(
  "CAPTION (Fig. 5E): Assay-native context for the displayed subset of the 22 frozen Figure 3 programs. ",
  "Rows are included when they contain at least one lineage-matched open-promoter gene at a fine-mapped ",
  "liver-enzyme locus, pass the batch-stratified Freedman-Lane histology-burden test, or show robust spatial ",
  "organization in both cohorts. The direct MASLD/PDFF ",
  "credible-set subset contributes zero open-promoter hits. Histology estimates are rank-residual partial ",
  "Spearman correlations with patient-bootstrap intervals. Spatial bars show the cohort range and a dark ",
  "shared core only when both native matched-null tests pass. Display filtering never triggers re-testing, ",
  "and no cross-modal score or rank is calculated."
)
