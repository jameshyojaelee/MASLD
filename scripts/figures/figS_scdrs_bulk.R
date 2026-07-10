#!/usr/bin/env Rscript
# ============================================================================
# figS_scdrs_bulk.R  (v2, 2026-05-12 overhaul)
#
# Bulk-DEG-anchored scDRS — question-driven supplementary figure.
# Each panel asks a biological question in its title and states the answer
# (computed from the data) in its subtitle.
#
#   A — Q1. Which cell types harbor F-transition genetic risk?
#       Ranked CT forest plot (12 CT × 4 F-transitions). Z-score per signed
#       group-analysis; FDR<0.05 dots filled.
#
#   B — Q2. Do disease-state cells activate F-transition signatures more
#       than control cells of the same type?
#       Per-cell scDRS norm_score by disease_stage_coarse × cell_type, for
#       the top-4 CTs. Violin + median crossbar. Wilcoxon p (Healthy vs
#       Steatohepatitis) annotated in each facet.
#
#   C — Q3. Does activation climb monotonically with disease stage?
#       Median norm_score (± IQR ribbon) per disease stage × CT × transition.
#       Line plot per CT, faceted by transition.
#
#   D — Visual. Where does the signal land inside the macrophage manifold?
#       Macrophage UMAP × per-cell scDRS norm_score per F-transition.
#
# Inputs:
#   Analysis/SingleCell/results_gpu_v2/disease_signatures/
#     scdrs_celltype_enrichment.csv          (Panel A — bulk-DEG CT enrichment)
#     scdrs_cell_scores.parquet              (Panels B/C/D — per-cell scores)
#     scdrs_obs_metadata.csv                 (Panels B/C — disease_stage_coarse)
#   Analysis/SingleCell/results_gpu_v2/pseudotime/Macrophages_metadata.csv
#                                            (Panel D — UMAP coords)
#
# Outputs:
#   figures/supplementary/figS_scdrs/figS_scdrs_bulk.pdf
#   figures/supplementary/figS_scdrs/panels/figS_scdrs_bulk_{A..D}_*.pdf
#   figures/supplementary/figS_scdrs/panel_data/figS_scdrs_bulk_panel{A..D}_data.csv
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUT_DIR    <- FIGS_SCDRS_DIR
DATA_DIR   <- FIGS_SCDRS_DATA_DIR
PANELS_DIR <- file.path(FIGS_SCDRS_DIR, "panels")
dir.create(PANELS_DIR, recursive = TRUE, showWarnings = FALSE)

SIG_DIR    <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/disease_signatures")
PREP_DIR   <- file.path(SIG_DIR, "figure_prep")  # produced by figS_scdrs_prep.py
CT_ENRICH  <- file.path(SIG_DIR, "scdrs_celltype_enrichment.csv")
BULK_CS    <- file.path(PREP_DIR, "bulk_cell_scores_top4ct_signed.csv.gz")
BULK_MED   <- file.path(PREP_DIR, "bulk_cell_scores_median_by_stage.csv")
MAC_UMAP   <- file.path(PREP_DIR, "macrophage_umap_signed_scores.csv.gz")

# Sanity check: complain loudly if prep CSVs are absent
for (f in c(BULK_CS, BULK_MED, MAC_UMAP)) {
  if (!file.exists(f))
    stop(sprintf("Missing %s — run scripts/figures/figS_scdrs_prep.py first.", f))
}

log_msg <- function(...) cat("[figS_scdrs_bulk] ", ..., "\n", sep = "")
save_panel <- function(p, slug, width, height) {
  out <- file.path(PANELS_DIR, paste0("figS_scdrs_bulk_", slug, ".pdf"))
  save_fig(p, out, width = width, height = height)
  log_msg("  panel saved: ", out)
}

BASE_SIZE <- 7
LBL_PT    <- 7
LBL_SIZE  <- LBL_PT / ggplot2::.pt

theme_scdrs <- function() {
  theme_masld(base_size = BASE_SIZE) +
    theme(
      plot.title    = element_text(size = BASE_SIZE, face = "plain"),
      plot.subtitle = element_text(size = BASE_SIZE - 1, colour = "grey25",
                                   lineheight = 1.05),
      axis.text.x   = element_text(angle = 45, hjust = 1),
      legend.key.size = unit(0.3, "cm"),
      panel.grid    = element_blank(),
      strip.background = element_rect(fill = "grey96", colour = NA),
      strip.text       = element_text(size = BASE_SIZE - 1, face = "plain")
    )
}

# ---- Cell-type aggregation (fine -> umbrella) ----
umbrella_map <- c(
  "Hepatocytes"             = "Hepatocytes",
  "Cholangiocytes"          = "Cholangiocytes",
  "Endothelial cells"       = "Endothelial cells",
  "Fibroblasts"             = "Fibroblasts",
  "Macrophages"             = "Macrophages",
  "Mono+mono derived cells" = "Monocytes",
  "Mono+mono_derived"       = "Monocytes",
  "cDC1s"                   = "Dendritic cells",
  "cDC2s"                   = "Dendritic cells",
  "pDCs"                    = "Dendritic cells",
  "Mig.cDCs"                = "Dendritic cells",
  "Neutrophils"             = "Granulocytes",
  "Basophils"               = "Granulocytes",
  "T cells"                 = "T cells",
  "B cells"                 = "B cells",
  "Plasma cells"            = "Plasma cells",
  "Resident NK"             = "NK cells",
  "Circulating NK/NKT"      = "NK cells"
)
umbrella_palette <- c(
  "Hepatocytes"       = "#1B5E20",
  "Cholangiocytes"    = "#00BCD4",
  "Endothelial cells" = "#7CB342",
  "Fibroblasts"       = "#FF6F00",
  "Macrophages"       = "#212121",
  "Monocytes"         = "#EC407A",
  "Dendritic cells"   = "#6A1B9A",
  "Granulocytes"      = "#F9A825",
  "T cells"           = "#1A237E",
  "B cells"           = "#5D4037",
  "Plasma cells"      = "#8D6E63",
  "NK cells"          = "#00897B"
)
umbrella_order <- c(
  "Hepatocytes", "Cholangiocytes", "Endothelial cells", "Fibroblasts",
  "Macrophages", "Monocytes", "Dendritic cells", "Granulocytes",
  "T cells", "B cells", "Plasma cells", "NK cells"
)
TOP_CT <- c("Hepatocytes", "Macrophages", "Fibroblasts", "Cholangiocytes")

trans_levels <- c("F1_vs_F0", "F2_vs_F1", "F3_vs_F2", "F4_vs_F3")
trans_labels <- c(
  "F1_vs_F0" = "F0→F1",
  "F2_vs_F1" = "F1→F2",
  "F3_vs_F2" = "F2→F3",
  "F4_vs_F3" = "F3→F4"
)

# Disease-stage axis (drop "Unknown"; order = ordinal severity)
STAGE_ORDER  <- c("Healthy", "Steatosis", "Steatohepatitis")
STAGE_LABELS <- c("Healthy" = "Healthy", "Steatosis" = "Steatosis",
                  "Steatohepatitis" = "Steatohepatitis")
STAGE_PALETTE <- c(
  "Healthy"          = "#9E9E9E",   # neutral gray per project convention
  "Steatosis"        = "#F48FB1",
  "Steatohepatitis"  = "#C2185B"
)

# ----------------------------------------------------------------------------
# Panel A — Q1. Which cell types harbor F-transition genetic risk?
# ----------------------------------------------------------------------------
log_msg("Panel A: cell-type forest (Q1)")
ct <- fread(CT_ENRICH)[mode == "signed"]
ct[, umbrella := umbrella_map[as.character(cell_type)]]
ct <- ct[!is.na(umbrella)]

# Aggregate fine -> umbrella (weighted by n_cell). Combine assoc_mcz into
# a single representative z by Stouffer-style weighting on signed z.
A_data <- ct[, .(
  z      = weighted.mean(assoc_mcz, w = n_cell, na.rm = TRUE),
  mcp    = min(assoc_mcp, na.rm = TRUE),
  n_sig  = sum(`n_fdr_0.05`, na.rm = TRUE),
  n_cell = sum(n_cell, na.rm = TRUE)
), by = .(umbrella, transition)]
A_data[, transition := factor(transition, levels = trans_levels,
                              labels = trans_labels[trans_levels])]
A_data[, sig := mcp < 0.05]
# Rank CTs by max |z| across transitions for the forest order
ct_rank <- A_data[, .(rank_z = max(abs(z), na.rm = TRUE)), by = umbrella][
  order(-rank_z)]
A_data[, umbrella := factor(umbrella, levels = rev(ct_rank$umbrella))]
fwrite(A_data, file.path(DATA_DIR, "figS_scdrs_bulk_panelA_data.csv"))

# Findings sentence (Stouffer-weighted top CT for each transition)
top_per_trans <- A_data[sig == TRUE,
  .SD[which.max(z)], by = transition]
finding_A <- if (nrow(top_per_trans) > 0) {
  paste0("Top FDR<0.05 cell type per transition: ",
         paste(sprintf("%s = %s (z=%.1f)",
                       top_per_trans$transition,
                       top_per_trans$umbrella,
                       top_per_trans$z),
               collapse = "; "))
} else "No cell type reached FDR<0.05 in any transition."

panel_A <- ggplot(A_data, aes(x = z, y = umbrella, colour = umbrella)) +
  geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_segment(aes(x = 0, xend = z, yend = umbrella),
               linewidth = 0.4, alpha = 0.6) +
  geom_point(aes(size = -log10(pmax(mcp, 1e-4)),
                 shape = sig, fill = umbrella),
             stroke = 0.4) +
  facet_wrap(~ transition, nrow = 1) +
  scale_colour_manual(values = umbrella_palette, guide = "none") +
  scale_fill_manual(values = umbrella_palette, guide = "none") +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 19),
                     labels = c("FDR≥0.05", "FDR<0.05"),
                     name = NULL) +
  scale_size_continuous(range = c(0.8, 3.5),
                        name = "-log10\n(group mcp)",
                        breaks = c(1, 2, 3)) +
  labs(x = "scDRS group z-score (Stouffer-weighted across fine CTs)", y = NULL,
       title = "A — Q1. Which cell types harbor F-transition genetic risk?",
       subtitle = finding_A) +
  theme_scdrs() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "right",
        legend.box = "vertical")

save_panel(panel_A, "A_Q1_celltype_forest", width = 7.5, height = 3.6)

# ----------------------------------------------------------------------------
# Load prepped per-cell scores (subsampled, joined with disease stage)
# ----------------------------------------------------------------------------
log_msg("Loading prepped per-cell tables (figS_scdrs_prep.py outputs)")
cs <- fread(BULK_CS)
cs[, cell_type   := factor(cell_type, levels = TOP_CT)]
cs[, disease_stage := factor(disease_stage, levels = STAGE_ORDER)]
cs[, transition  := factor(transition, levels = trans_levels,
                           labels = trans_labels[trans_levels])]
log_msg("  per-cell rows: ", nrow(cs))

# ----------------------------------------------------------------------------
# Panel B — Q2. Do disease-state cells activate signatures more than control?
# ----------------------------------------------------------------------------
log_msg("Panel B: disease-state violins (Q2)")

# Two-sided Wilcoxon (direction-agnostic); delta = sign of biology.
# Positive delta = scDRS goes UP from healthy to steatohepatitis (activation).
# Negative delta = scDRS goes DOWN (depletion / identity loss).
pvals_B <- cs[, {
  h <- norm_score[disease_stage == "Healthy"]
  sh <- norm_score[disease_stage == "Steatohepatitis"]
  if (length(h) >= 5 & length(sh) >= 5) {
    w <- suppressWarnings(wilcox.test(sh, h))
    .(p = w$p.value,
      delta = median(sh) - median(h),
      n_h  = length(h),
      n_sh = length(sh))
  } else {
    .(p = NA_real_, delta = NA_real_, n_h = length(h), n_sh = length(sh))
  }
}, by = .(cell_type, transition)]
pvals_B[, p_label := ifelse(is.na(p), "n.s.",
                            ifelse(p < 1e-300, "p<1e-300",
                                   ifelse(p < 0.001, sprintf("p=%.1g", p),
                                          sprintf("p=%.3f", p))))]
pvals_B[, dir_label := ifelse(is.na(delta), "·",
                              ifelse(delta > 0, "↑disease",
                                     ifelse(delta < 0, "↓disease", "·")))]
pvals_B[, sig := !is.na(p) & p < 0.05]

# cs is already subsampled by the prep step (cap MAX_PER_BIN per bin)
cs_B <- cs
log_msg("  per-bin cells for violin/strip: ", nrow(cs_B))

fwrite(pvals_B, file.path(DATA_DIR, "figS_scdrs_bulk_panelB_data.csv"))

# Direction-aware finding (compact summary; per-facet Δ shown in plot).
sig_pv <- pvals_B[!is.na(p) & p < 0.05]
up_fac   <- sig_pv[delta > 0]
down_fac <- sig_pv[delta < 0]
standout <- if (nrow(up_fac) > 0) {
  u <- up_fac[which.max(delta)]
  sprintf("%s @ %s is the sole ↑disease facet (Δ=%+.2f)",
          u$cell_type, u$transition, u$delta)
} else "no facets show ↑disease"
finding_B <- sprintf(
  "%d/%d facets p<0.05 — %d ↑disease, %d ↓disease. %s.",
  nrow(sig_pv), nrow(pvals_B), nrow(up_fac), nrow(down_fac), standout
)

# Color the facet border by direction (red = ↑disease, blue = ↓disease).
panel_B <- ggplot(cs_B,
                  aes(x = disease_stage, y = norm_score, fill = disease_stage)) +
  geom_violin(scale = "width", width = 0.85, linewidth = 0.2,
              draw_quantiles = c(0.5), trim = TRUE, alpha = 0.85) +
  geom_text(data = pvals_B,
            aes(x = 2, y = Inf,
                label = paste0(dir_label, "\n", p_label),
                colour = ifelse(!sig, "grey60",
                                ifelse(delta > 0, "#B40426", "#3B4CC0"))),
            inherit.aes = FALSE,
            size = LBL_SIZE * 0.85, vjust = 1.15, lineheight = 0.9) +
  scale_colour_identity() +
  facet_grid(cell_type ~ transition, switch = "y") +
  scale_fill_manual(values = STAGE_PALETTE, name = "Disease stage") +
  scale_y_continuous(limits = c(-3, 4.5),
                     oob = scales::squish, expand = c(0, 0.05)) +
  labs(x = NULL, y = "per-cell scDRS norm_score",
       title = "B — Q2. How does per-cell scDRS shift from healthy to disease-state cells of the same type?",
       subtitle = finding_B) +
  theme_scdrs() +
  theme(legend.position = "bottom",
        strip.placement = "outside",
        axis.text.x  = element_text(angle = 0, hjust = 0.5, size = BASE_SIZE - 2))

save_panel(panel_B, "B_Q2_disease_state_violins", width = 7.5, height = 5.0)

# ----------------------------------------------------------------------------
# Panel C — Q3. Does activation climb monotonically with disease stage?
# ----------------------------------------------------------------------------
log_msg("Panel C: stage trajectory (Q3)")

# Use the full-population median table (no subsampling) for the trajectory.
C_data <- fread(BULK_MED)
setnames(C_data, c("median", "n"), c("median_z", "n_cells"))
C_data[, cell_type    := factor(cell_type, levels = TOP_CT)]
C_data[, disease_stage := factor(disease_stage, levels = STAGE_ORDER)]
C_data[, transition   := factor(transition, levels = trans_levels,
                                labels = trans_labels[trans_levels])]
fwrite(C_data, file.path(DATA_DIR, "figS_scdrs_bulk_panelC_data.csv"))

# Slope (Healthy median -> Steatohepatitis median) per CT × transition.
# Sort by |slope| and label as activating (↑) vs depleting (↓).
slope_C <- C_data[, {
  h  <- median_z[disease_stage == "Healthy"]
  sh <- median_z[disease_stage == "Steatohepatitis"]
  .(slope = if (length(h) & length(sh)) sh - h else NA_real_)
}, by = .(cell_type, transition)]
top_up   <- head(slope_C[!is.na(slope) & slope > 0][order(-slope)], 1)
top_down <- head(slope_C[!is.na(slope) & slope < 0][order(slope)], 1)
fmt_one <- function(dt) {
  if (nrow(dt) == 0) return("(none)")
  sprintf("%s @ %s (Δ=%+.2f)", dt$cell_type, dt$transition, dt$slope)
}
finding_C <- paste0(
  "Top ↑activator: ", fmt_one(top_up),
  "  ·  Top ↓depleter: ", fmt_one(top_down)
)

panel_C <- ggplot(C_data,
                  aes(x = disease_stage, y = median_z,
                      group = cell_type, colour = cell_type)) +
  geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_ribbon(aes(ymin = q25, ymax = q75, fill = cell_type),
              alpha = 0.18, colour = NA) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.8) +
  facet_wrap(~ transition, nrow = 1) +
  scale_colour_manual(values = umbrella_palette[TOP_CT], name = "cell type") +
  scale_fill_manual(values   = umbrella_palette[TOP_CT], guide = "none") +
  labs(x = "Disease stage", y = "median per-cell scDRS (IQR ribbon)",
       title = "C — Q3. How does per-cell scDRS shift across disease stages (Healthy → Steatosis → Steatohepatitis)?",
       subtitle = finding_C) +
  theme_scdrs() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "right")

save_panel(panel_C, "C_Q3_stage_trajectory", width = 7.5, height = 3.4)

# ----------------------------------------------------------------------------
# Panel D — Visual. Macrophage UMAP × per-cell scDRS norm_score
# ----------------------------------------------------------------------------
log_msg("Panel D: macrophage UMAP per F-transition")

D_join <- fread(MAC_UMAP)
D_join[, transition_lab := factor(transition, levels = trans_levels,
                                  labels = trans_labels[trans_levels])]
log_msg("  macrophage UMAP rows: ", nrow(D_join))

fwrite(D_join[, .(cell_id, UMAP_1, UMAP_2,
                  transition = transition_lab, norm_score)],
       file.path(DATA_DIR, "figS_scdrs_bulk_panelD_data.csv"))

zmax_D <- quantile(abs(D_join$norm_score), 0.99, na.rm = TRUE)
D_join[, norm_clip := pmax(pmin(norm_score, zmax_D), -zmax_D)]
frac_pos <- D_join[, .(frac_pos = mean(norm_score > 0, na.rm = TRUE)),
                   by = transition_lab]
finding_D <- paste0("Fraction of macrophages with norm_score > 0 per transition: ",
                    paste(sprintf("%s = %.0f%%",
                                  frac_pos$transition_lab,
                                  100 * frac_pos$frac_pos),
                          collapse = "; "))

panel_D <- ggplot(D_join, aes(x = UMAP_1, y = UMAP_2, colour = norm_clip)) +
  rasterize_layer(geom_point(size = 0.05, stroke = 0, alpha = 0.55), dpi = 300) +
  facet_wrap(~ transition_lab, ncol = 4) +
  scale_colour_gradient2(low = masld_colors$down, mid = "grey88",
                         high = masld_colors$up, midpoint = 0,
                         limits = c(-zmax_D, zmax_D),
                         name = "scDRS\nnorm_score",
                         guide = guide_colorbar(barwidth = unit(2.4, "cm"),
                                                barheight = unit(0.25, "cm"))) +
  labs(x = "UMAP 1", y = "UMAP 2",
       title = "D — Where does the F-transition signal land inside macrophages?",
       subtitle = finding_D) +
  theme_scdrs() +
  theme(legend.position = "bottom",
        axis.text       = element_blank(),
        axis.ticks      = element_blank())

save_panel(panel_D, "D_macrophage_umap_per_transition", width = 7.5, height = 3.4)

# ----------------------------------------------------------------------------
# Compose
# ----------------------------------------------------------------------------
log_msg("Composing figS_scdrs_bulk.pdf")
full <- panel_A / panel_B / panel_C / panel_D +
  plot_layout(heights = c(0.85, 1.30, 0.85, 0.95))

out_pdf <- file.path(OUT_DIR, "figS_scdrs_bulk.pdf")
save_fig(full, out_pdf, width = 7.5, height = 14)
log_msg("Wrote ", out_pdf)
log_msg("Done.")
