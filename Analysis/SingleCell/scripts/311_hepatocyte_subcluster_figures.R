#!/usr/bin/env Rscript
# 311: Hepatocyte Subclustering — Supplementary Figure Panel
#
# Groups 43 Leiden subtypes into ~5 interpretable meta-subtypes based on
# disease enrichment class + fibrosis stage correlation, then generates
# individual panel PDFs for flexible arrangement.
#
# Meta-subtype grouping:
#   Disease-Progressor: MASLD-enriched, stage rho > 0.25 (sig)
#   Disease-Associated:   MASLD-enriched, 0.1 < rho <= 0.25 (sig)
#   Disease-Neutral:  MASLD-enriched, rho <= 0.1 or n.s.
#   Healthy:          Healthy_enriched
#   Neutral:          neutral enrichment class
#
# Input: hepatocyte_subtypes/ CSVs from 309-310
# Output:
#   panels/panelA_umap_meta.pdf, panelB_umap_condition.pdf, ...
#   figS_hepatocyte_subtypes.pdf (combined)
#   meta_subtype_mapping.csv
#
# Environment: rnaseq (CPU)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

SUB_DIR  <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes")
OUT_DIR  <- FIGS_HEPSUB_DIR
PANEL_DIR <- file.path(OUT_DIR, "panels")
dir.create(PANEL_DIR, recursive = TRUE, showWarnings = FALSE)

# ═══════════════════════════════════════════════════════════════════════════
# Load data
# ═══════════════════════════════════════════════════════════════════════════
message("Loading data...")
meta       <- fread(file.path(SUB_DIR, "hepatocyte_subtype_metadata.csv"))
axis_scores <- fread(file.path(SUB_DIR, "subtype_axis_scores.csv"))
enrichment <- fread(file.path(SUB_DIR, "subtype_disease_enrichment.csv"))
markers    <- fread(file.path(SUB_DIR, "subtype_markers.csv"))

# ═══════════════════════════════════════════════════════════════════════════
# Meta-subtype grouping
# ═══════════════════════════════════════════════════════════════════════════
message("Computing meta-subtypes...")

enr <- enrichment[, .(subtype, enrichment_class, stage_spearman_rho, stage_spearman_p)]

enr[, meta_subtype := fcase(
  enrichment_class == "Healthy_enriched", "Healthy",
  enrichment_class == "neutral",          "Neutral",
  # MASLD-enriched: split by stage correlation
  enrichment_class == "MASLD_enriched" &
    stage_spearman_rho > 0.25 & stage_spearman_p < 0.05, "Disease-Progressor",
  enrichment_class == "MASLD_enriched" &
    stage_spearman_rho > 0.10 & stage_spearman_p < 0.05, "Disease-Associated",
  default = "Disease-Neutral"
)]

# Save the mapping
mapping <- enr[, .(subtype, enrichment_class, stage_spearman_rho, stage_spearman_p, meta_subtype)]
fwrite(mapping, file.path(SUB_DIR, "meta_subtype_mapping.csv"))
message("Meta-subtype mapping:")
print(mapping[, .N, by = meta_subtype][order(-N)])

# Merge onto cell metadata
meta[, subtype_chr := as.character(hepatocyte_subtype)]
mapping[, subtype_chr := as.character(subtype)]
meta <- merge(meta, mapping[, .(subtype_chr, meta_subtype)],
              by = "subtype_chr", all.x = TRUE)

# Ordered factor for consistent display
meta_levels <- c("Disease-Progressor", "Disease-Associated", "Disease-Neutral", "Neutral", "Healthy")
meta[, meta_subtype := factor(meta_subtype, levels = meta_levels)]

# ═══════════════════════════════════════════════════════════════════════════
# Color palette — 5 meta-subtypes
# ═══════════════════════════════════════════════════════════════════════════
meta_colors <- c(
  "Disease-Progressor" = "#880E4F",  # Dark magenta (most disease-associated)
  "Disease-Associated"   = "#E91E63",  # Bright magenta
  "Disease-Neutral"     = "#F48FB1",
  "Neutral"          = "#BDBDBD",  # Soft pink
  "Healthy"          = "#1565C0"   # Deep blue
)

disease_stage_colors <- c(
  Healthy = "#E3F2FD", Steatosis = "#90CAF9",
  Steatohepatitis = "#AB47BC", Cirrhosis = "#6A1B9A"
)

# Subsample for UMAP plotting
set.seed(42)
n_plot <- min(100000, nrow(meta))
plot_meta <- meta[sample(.N, n_plot)]

# ═══════════════════════════════════════════════════════════════════════════
# Panel A: UMAP by meta-subtype
# ═══════════════════════════════════════════════════════════════════════════
message("Panel A: UMAP by meta-subtype")
pA <- ggplot(plot_meta, aes(x = UMAP_1, y = UMAP_2, color = meta_subtype)) +
  geom_point(size = 0.05, alpha = 0.3) +
  scale_color_manual(values = meta_colors) +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  labs(title = "Hepatocyte subtypes", color = "Meta-subtype") +
  theme_masld() +
  theme(legend.position = "right",
        legend.key.size = unit(0.4, "cm"))

save_fig(pA, file.path(PANEL_DIR, "panelA_umap_meta.pdf"),
         width = fig_half_width + 1, height = 3.5)

# ═══════════════════════════════════════════════════════════════════════════
# Panel B: UMAP split by condition
# ═══════════════════════════════════════════════════════════════════════════
message("Panel B: UMAP split by condition")
plot_meta[, condition_plot := fifelse(
  condition %in% c("MASLD", "NAFLD", "NASH", "Cirrhotic"), "MASLD", condition)]
plot_meta_cond <- plot_meta[condition_plot %in% c("Healthy", "MASLD")]

pB <- ggplot(plot_meta_cond, aes(x = UMAP_1, y = UMAP_2, color = meta_subtype)) +
  geom_point(size = 0.05, alpha = 0.3) +
  scale_color_manual(values = meta_colors) +
  facet_wrap(~condition_plot) +
  labs(title = "By condition") +
  theme_masld() +
  theme(legend.position = "none")

save_fig(pB, file.path(PANEL_DIR, "panelB_umap_condition.pdf"),
         width = fig_full_width, height = 3.5)

# ═══════════════════════════════════════════════════════════════════════════
# Panel C: Stacked barplot — proportions per condition
# ═══════════════════════════════════════════════════════════════════════════
message("Panel C: Proportion barplot")
prop_dt <- meta[condition %in% c("Healthy", "MASLD", "NAFLD", "NASH", "Cirrhotic")]
prop_dt[, condition_plot := fifelse(
  condition %in% c("MASLD", "NAFLD", "NASH", "Cirrhotic"), "MASLD", condition)]
prop_dt <- prop_dt[condition_plot %in% c("Healthy", "MASLD")]

sample_props <- prop_dt[, .N, by = .(sample, condition_plot, meta_subtype)]
sample_props[, prop := N / sum(N), by = sample]
mean_props <- sample_props[, .(mean_prop = mean(prop)), by = .(condition_plot, meta_subtype)]

pC <- ggplot(mean_props, aes(x = condition_plot, y = mean_prop, fill = meta_subtype)) +
  geom_bar(stat = "identity", position = "stack", width = 0.75) +
  scale_fill_manual(values = meta_colors) +
  labs(title = "Subtype proportions", x = "", y = "Mean proportion", fill = "Meta-subtype") +
  theme_masld()

save_fig(pC, file.path(PANEL_DIR, "panelC_proportions.pdf"),
         width = 2.2, height = 3.5)

# ═══════════════════════════════════════════════════════════════════════════
# Panel D: Boxplot — proportion across fibrosis stages
# ═══════════════════════════════════════════════════════════════════════════
message("Panel D: Disease stage proportion")
meta_stage <- meta[!is.na(disease_stage_coarse) & disease_stage_coarse != ""]
stage_props <- meta_stage[, .N, by = .(sample, disease_stage_coarse, meta_subtype)]
stage_props[, prop := N / sum(N), by = sample]

stage_props[, disease_stage_coarse := factor(
  disease_stage_coarse,
  levels = c("Healthy", "Steatosis", "Steatohepatitis", "Cirrhosis")
)]

pD <- ggplot(stage_props, aes(x = disease_stage_coarse, y = prop, fill = disease_stage_coarse)) +
  geom_boxplot(outlier.size = 0.3) +
  scale_fill_manual(values = disease_stage_colors, guide = "none") +
  facet_wrap(~meta_subtype, nrow = 1, scales = "free_y") +
  labs(title = "Meta-subtype proportion by disease stage",
       x = "Disease stage", y = "Proportion") +
  theme_masld() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

save_fig(pD, file.path(PANEL_DIR, "panelD_stage_proportions.pdf"),
         width = fig_full_width, height = 3)

# ═══════════════════════════════════════════════════════════════════════════
# Panel E: Dot plot — top markers per meta-subtype
# ═══════════════════════════════════════════════════════════════════════════
message("Panel E: Marker dot plot")
markers[, subtype_chr := as.character(subtype)]
markers <- merge(markers, mapping[, .(subtype_chr, meta_subtype)],
                 by = "subtype_chr", all.x = TRUE)

# Top 3 unique genes per meta-subtype by highest logfoldchange
top_meta <- markers[, {
  # Sort by logfoldchanges descending, pick top unique genes
  .SD[order(-logfoldchanges)][!duplicated(names)][1:min(5, .N)]
}, by = meta_subtype]

# Remove NA rows
top_meta <- top_meta[!is.na(names)]

pE <- ggplot(top_meta, aes(x = meta_subtype, y = reorder(names, logfoldchanges),
                           size = pct_nz_group, color = logfoldchanges)) +
  geom_point() +
  scale_color_gradient2(low = masld_colors$down, mid = "white", high = masld_colors$up,
                        midpoint = 0) +
  scale_size_continuous(range = c(0.5, 4)) +
  labs(title = "Top marker genes per meta-subtype", x = "", y = "",
       color = "logFC", size = "% expressed") +
  theme_masld() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

save_fig(pE, file.path(PANEL_DIR, "panelE_markers.pdf"),
         width = fig_half_width + 1, height = 4)

# ═══════════════════════════════════════════════════════════════════════════
# Panel F: Heatmap — axis scores (meta-subtype-level means)
# ═══════════════════════════════════════════════════════════════════════════
message("Panel F: Axis score heatmap")

# Compute mean axis scores per meta-subtype
axis_scores[, subtype_chr := as.character(subtype)]
axis_scores <- merge(axis_scores, mapping[, .(subtype_chr, meta_subtype)],
                     by = "subtype_chr", all.x = TRUE)

axis_cols <- grep("^axis_", names(axis_scores), value = TRUE)
axis_meta <- axis_scores[, lapply(.SD, mean), by = meta_subtype, .SDcols = axis_cols]

axis_mat <- as.matrix(axis_meta[, ..axis_cols])
rownames(axis_mat) <- axis_meta$meta_subtype
colnames(axis_mat) <- gsub("axis_", "", colnames(axis_mat))

# Scale columns
axis_scaled <- scale(axis_mat)

# Row order
row_order <- match(meta_levels, rownames(axis_scaled))
row_order <- row_order[!is.na(row_order)]
axis_scaled <- axis_scaled[row_order, , drop = FALSE]

pF_hm <- Heatmap(
  axis_scaled,
  name = "Z-score",
  col = colorRamp2(c(-2, 0, 2), c(masld_colors$down, "white", masld_colors$up)),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_title = "Meta-subtype",
  column_title = "Biological axis",
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  column_names_rot = 45,
  heatmap_width = unit(8, "cm"),
  heatmap_height = unit(4, "cm")
)

pdf(file.path(PANEL_DIR, "panelF_axis_heatmap.pdf"), width = 5, height = 3.5)
draw(pF_hm)
dev.off()

# ═══════════════════════════════════════════════════════════════════════════
# Panel G: Spatial ssGSEA overlay — best meta-subtype
# ═══════════════════════════════════════════════════════════════════════════
message("Panel G: Spatial overlay")
spatial_scores_file <- file.path(SUB_DIR, "crossmodal/spatial/subtype_spatial_scores.csv")
pG <- NULL
if (file.exists(spatial_scores_file)) {
  sp_scores <- fread(spatial_scores_file)

  # Read spatial coordinates
  library(rhdf5)
  sp_h5 <- file.path(BASE, "Analysis/Spatial/results/preprocessed/merged_spatial.h5ad")
  if (file.exists(sp_h5)) {
    sp_coords <- h5read(sp_h5, "obsm/spatial")
    if (nrow(sp_coords) == 2) sp_coords <- t(sp_coords)
    sp_scores[, x := sp_coords[1:nrow(sp_scores), 1]]
    sp_scores[, y := sp_coords[1:nrow(sp_scores), 2]]

    # Average scores per meta-subtype: map Hep_N columns to meta-subtypes
    score_cols <- grep("^Hep_", names(sp_scores), value = TRUE)

    if (length(score_cols) > 0) {
      # Map each Hep_N column to its meta-subtype
      col_to_meta <- mapping[, .(col = paste0("Hep_", subtype), meta_subtype)]

      for (ms in meta_levels) {
        ms_cols <- col_to_meta[meta_subtype == ms, col]
        ms_cols <- ms_cols[ms_cols %in% score_cols]
        if (length(ms_cols) > 0) {
          sp_scores[, (ms) := rowMeans(.SD, na.rm = TRUE), .SDcols = ms_cols]
        }
      }

      # Show Disease-Progressor (most interesting)
      best_col <- "Disease-Progressor"
      if (best_col %in% names(sp_scores)) {
        pG <- ggplot(sp_scores, aes(x = x, y = -y, color = get(best_col))) +
          geom_point(size = 0.3) +
          scale_color_gradient2(low = masld_colors$down, mid = "white", high = masld_colors$up,
                                midpoint = median(sp_scores[[best_col]], na.rm = TRUE)) +
          labs(title = "Spatial: Disease-Progressor signature", color = "ssGSEA") +
          coord_fixed() +
          theme_masld() +
          theme(axis.text = element_blank(), axis.ticks = element_blank(),
                axis.title = element_blank())

        save_fig(pG, file.path(PANEL_DIR, "panelG_spatial.pdf"),
                 width = fig_half_width, height = 3.5)
      }
    }
  }
}

# ═══════════════════════════════════════════════════════════════════════════
# Panel H: Plasma proteomics — fibrosis association per meta-subtype
# ═══════════════════════════════════════════════════════════════════════════
message("Panel H: Plasma fibrosis association")
plasma_file <- file.path(SUB_DIR, "crossmodal/plasma/subtype_plasma_fibrosis_association.csv")
pH <- NULL
if (file.exists(plasma_file)) {
  plasma_assoc <- fread(plasma_file)
  # Map Hep_N to meta-subtype
  plasma_assoc[, subtype_num := as.integer(gsub("Hep_", "", subtype))]
  plasma_assoc <- merge(plasma_assoc,
                        mapping[, .(subtype, meta_subtype)],
                        by.x = "subtype_num", by.y = "subtype", all.x = TRUE)

  # Aggregate to meta-subtype level: weighted mean rho
  plasma_meta <- plasma_assoc[, .(
    mean_rho = mean(spearman_rho, na.rm = TRUE),
    min_pval = min(spearman_pval, na.rm = TRUE),
    n_subtypes = .N
  ), by = meta_subtype]
  plasma_meta[, meta_subtype := factor(meta_subtype, levels = meta_levels)]
  plasma_meta[, sig := fifelse(min_pval < 0.05, "*", "")]

  pH <- ggplot(plasma_meta[!is.na(meta_subtype)],
               aes(x = meta_subtype, y = mean_rho, fill = meta_subtype)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = sig), vjust = -0.5, size = 5) +
    scale_fill_manual(values = meta_colors, guide = "none") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    labs(title = "Plasma signature vs fibrosis stage",
         x = "", y = "Mean Spearman rho") +
    theme_masld() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  save_fig(pH, file.path(PANEL_DIR, "panelH_plasma_fibrosis.pdf"),
           width = 2.6, height = 3)
}

# ═══════════════════════════════════════════════════════════════════════════
# Panel I: Enrichment summary — OR + rho per fine-grained subtype
# ═══════════════════════════════════════════════════════════════════════════
message("Panel I: Enrichment summary")
enr_plot <- copy(enrichment)
enr_plot <- merge(enr_plot, mapping[, .(subtype, meta_subtype)],
                  by = "subtype", all.x = TRUE)
enr_plot[, meta_subtype := factor(meta_subtype, levels = meta_levels)]
enr_plot[, log10_or := log10(pmax(odds_ratio, 1e-3))]  # Floor for Healthy (OR~0)

pI <- ggplot(enr_plot, aes(x = log10_or, y = stage_spearman_rho,
                           color = meta_subtype, size = n_cells / 1000)) +
  geom_point(alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  scale_color_manual(values = meta_colors) +
  scale_size_continuous(range = c(1, 6), name = "Cells (K)") +
  labs(title = "Subtype disease enrichment",
       x = "log10(OR) MASLD vs Healthy",
       y = "Fibrosis stage rho",
       color = "Meta-subtype") +
  theme_masld()

save_fig(pI, file.path(PANEL_DIR, "panelI_enrichment_scatter.pdf"),
         width = fig_half_width + 1, height = 3.5)

# ═══════════════════════════════════════════════════════════════════════════
# Combined figure (2 pages)
# ═══════════════════════════════════════════════════════════════════════════
message("Assembling combined figure...")
out_path <- file.path(OUT_DIR, "figS_hepatocyte_subtypes.pdf")

pdf(out_path, width = 14, height = 16)

# Page 1: ggplot panels
p_page1 <- (pA | pB) /
  (pC | pI) /
  pD /
  pE +
  plot_annotation(tag_levels = list(c("A", "B", "C", "I", "D", "E"))) +
  plot_layout(heights = c(1, 1, 0.8, 1.2))
print(p_page1)

# Page 2: Heatmap + spatial + plasma
grid::grid.newpage()
grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 2)))

grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(pF_hm, newpage = FALSE)
grid::upViewport()

if (!is.null(pG)) {
  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  print(pG, newpage = FALSE)
  grid::upViewport()
}

if (!is.null(pH)) {
  grid::pushViewport(grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
  print(pH, newpage = FALSE)
  grid::upViewport()
}

grid::upViewport()

dev.off()
message("Saved combined figure: ", out_path)
message("Saved individual panels: ", PANEL_DIR)
message("\nMeta-subtype summary:")
print(meta[, .N, by = meta_subtype][order(-N)])
message("\n311 complete.")
