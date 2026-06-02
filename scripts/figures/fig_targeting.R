#!/usr/bin/env Rscript
# fig_targeting.R — 4-panel targeting figure (TAS heatmap, heterogeneity,
#                   therapeutic modality, personalized targets)
# Output: figures/supplementary/figS10_prediction/fig_targeting.pdf

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(data.table)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

PROG  <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/progression")
STAGE <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/staging_classifier")
OUT   <- FIGS10_DIR

# ── Transition color palette (blue → red gradient) ──────────────────────────
transition_cols <- c(
  "F0_to_F1" = "#2196F3",  # blue

"F1_to_F2" = "#7B1FA2",  # violet
  "F2_to_F3" = "#E91E63",  # magenta
  "F3_to_F4" = "#B71C1C",  # dark red
  "quiescent" = "#BDBDBD",
  "unassigned" = "#E0E0E0"
)

# Modality palette
modality_cols <- c(
  "anti_fibrotic"      = "#C2185B",
  "anti_inflammatory"  = "#1565C0",
  "metabolic"          = "#F57F17",
  "vascular"           = "#00695C",
  "other"              = "#BDBDBD"
)

modality_labels <- c(
  "anti_fibrotic"      = "Anti-fibrotic",
  "anti_inflammatory"  = "Anti-inflammatory",
  "metabolic"          = "Metabolic",
  "vascular"           = "Vascular",
  "other"              = "Other"
)

# ── Load data ────────────────────────────────────────────────────────────────
tas  <- fread(file.path(PROG, "transition_activity_scores.csv"))
pt   <- fread(file.path(PROG, "consensus_pseudotime.csv"))
subs <- fread(file.path(PROG, "transition_subtype_assignments.csv"))
ppc  <- fread(file.path(PROG, "patient_pathway_class.csv"))
top10 <- fread(file.path(PROG, "patient_top10_targets.csv"))
meta <- fread(file.path(STAGE, "modeling_metadata.csv"))

# ══════════════════════════════════════════════════════════════════════════════
# (a) TAS heatmap: z-scored TAS × samples ordered by pseudotime
# ══════════════════════════════════════════════════════════════════════════════
pa <- tryCatch({
  # Merge TAS with pseudotime and fibrosis stage
  tas_cols <- grep("^TAS_F", names(tas), value = TRUE)  # fibrosis transitions only
  hm <- merge(tas[, c("sample_id", tas_cols), with = FALSE],
              pt[, .(sample_id, pseudotime_consensus)], by = "sample_id")
  hm <- merge(hm, meta[, .(sample_id, fibrosis_stage)], by = "sample_id",
              all.x = TRUE)

  # Z-score each TAS column
  for (col in tas_cols) {
    mu <- mean(hm[[col]], na.rm = TRUE)
    sd_val <- sd(hm[[col]], na.rm = TRUE)
    if (sd_val > 0) set(hm, j = col, value = (hm[[col]] - mu) / sd_val)
  }

  # Order by pseudotime
  hm <- hm[order(pseudotime_consensus)]
  hm[, row_idx := .I]

  # Map fibrosis stage
  hm[, fib_label := fifelse(is.na(fibrosis_stage) | fibrosis_stage < 0,
                             "NA",
                             paste0("F", fibrosis_stage))]

  # Pivot long for geom_tile
  hm_long <- melt(hm, id.vars = c("sample_id", "row_idx", "fib_label",
                                    "pseudotime_consensus"),
                  measure.vars = tas_cols, variable.name = "transition",
                  value.name = "z_score")

  # Clean transition labels
  hm_long[, transition := gsub("TAS_", "", transition)]
  hm_long[, transition := gsub("_to_", " -> ", transition)]

  # Clamp extreme z-scores for visual clarity
  hm_long[, z_clamped := pmin(pmax(z_score, -3), 3)]

  # Transition ordering
  trans_order <- c("F0 -> F1", "F1 -> F2", "F2 -> F3", "F3 -> F4")
  hm_long[, transition := factor(transition, levels = trans_order)]

  # Fibrosis stage annotation (build bands)
  fib_bands <- hm[, .(ymin = min(row_idx), ymax = max(row_idx)),
                   by = fib_label]
  fib_bands <- fib_bands[fib_label != "NA"]
  fib_bands[, ymid := (ymin + ymax) / 2]
  fib_bands[, fib_label := factor(fib_label,
    levels = c("F0", "F1", "F2", "F3", "F4"))]

  # Build heatmap with ggrastr if available
  tile_layer <- geom_tile(aes(x = transition, y = row_idx, fill = z_clamped))
  if (requireNamespace("ggrastr", quietly = TRUE)) {
    tile_layer <- ggrastr::rasterise(tile_layer, dpi = 300)
  }

  p <- ggplot(hm_long) +
    tile_layer +
    scale_fill_gradient2(low = "#1565C0", mid = "white", high = "#C2185B",
                         midpoint = 0, limits = c(-3, 3),
                         name = "TAS (z)") +
    # Fibrosis stage bands on left margin
    geom_rect(data = fib_bands,
              aes(xmin = 0.3, xmax = 0.5, ymin = ymin, ymax = ymax),
              fill = fibrosis_stage_colors[as.character(fib_bands$fib_label)],
              inherit.aes = FALSE) +
    geom_text(data = fib_bands,
              aes(x = 0.4, y = ymid, label = fib_label),
              size = 1.8, inherit.aes = FALSE) +
    scale_x_discrete(expand = expansion(add = c(0.5, 0.2))) +
    labs(x = NULL, y = "Patients (ordered by pseudotime)") +
    theme_masld(base_size = 7) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
          legend.key.width = unit(0.25, "cm"),
          legend.key.height = unit(0.5, "cm"))
  p
}, error = function(e) {
  message("Panel (a) error: ", conditionMessage(e))
  placeholder("Panel a: TAS heatmap\n(error)")
})

# ══════════════════════════════════════════════════════════════════════════════
# (b) Within-stage heterogeneity: stacked bar of dominant transition per stage
# ══════════════════════════════════════════════════════════════════════════════
pb <- tryCatch({
  # Merge subtype assignments with fibrosis stage from metadata
  het <- merge(subs, meta[, .(sample_id, fibrosis_stage)], by = "sample_id",
               all.x = TRUE)

  # Use fibrosis_stage from metadata (more complete); fallback to subs column
  het[, fib := fifelse(!is.na(fibrosis_stage.y) & fibrosis_stage.y >= 0,
                       fibrosis_stage.y,
                       as.numeric(fibrosis_stage.x))]
  het <- het[!is.na(fib) & fib >= 0]
  het[, fib_label := paste0("F", fib)]

  # Classify dominant transition
  het[, dom := fib_dominant]
  het[dom == "" | is.na(dom), dom := "unassigned"]

  # Compute proportions
  het_prop <- het[, .N, by = .(fib_label, dom)]
  het_prop[, total := sum(N), by = fib_label]
  het_prop[, prop := N / total]

  # Order
  het_prop[, fib_label := factor(fib_label,
    levels = c("F0", "F1", "F2", "F3", "F4"))]
  het_prop[, dom := factor(dom,
    levels = c("F0_to_F1", "F1_to_F2", "F2_to_F3", "F3_to_F4",
               "quiescent", "unassigned"))]

  p <- ggplot(het_prop, aes(x = fib_label, y = prop, fill = dom)) +
    geom_col(width = 0.7, color = "white", linewidth = 0.2) +
    scale_fill_manual(values = transition_cols, name = "Dominant\ntransition",
                      drop = FALSE) +
    scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.02))) +
    labs(x = "Fibrosis stage", y = "Proportion of patients") +
    theme_masld(base_size = 7) +
    theme(legend.key.size = unit(0.25, "cm"))
  p
}, error = function(e) {
  message("Panel (b) error: ", conditionMessage(e))
  placeholder("Panel b: Heterogeneity\n(error)")
})

# ══════════════════════════════════════════════════════════════════════════════
# (c) Therapeutic modality distribution by fibrosis stage
# ══════════════════════════════════════════════════════════════════════════════
pc <- tryCatch({
  mod <- merge(ppc, meta[, .(sample_id, fibrosis_stage)], by = "sample_id",
               all.x = TRUE)
  mod <- mod[!is.na(fibrosis_stage) & fibrosis_stage >= 0]
  mod[, fib_label := paste0("F", fibrosis_stage)]

  # Compute proportions per stage
  mod_prop <- mod[, .N, by = .(fib_label, dominant_modality)]
  mod_prop[, total := sum(N), by = fib_label]
  mod_prop[, prop := N / total]

  # Relabel
  mod_prop[, modality := modality_labels[dominant_modality]]
  mod_prop[is.na(modality), modality := dominant_modality]

  mod_prop[, fib_label := factor(fib_label,
    levels = c("F0", "F1", "F2", "F3", "F4"))]

  # Ordered factor for modalities
  mod_order <- c("Anti-fibrotic", "Anti-inflammatory", "Metabolic",
                 "Vascular", "Other")
  mod_prop[, modality := factor(modality, levels = mod_order)]

  p <- ggplot(mod_prop, aes(x = fib_label, y = prop, fill = modality)) +
    geom_col(width = 0.7, color = "white", linewidth = 0.2) +
    scale_fill_manual(
      values = setNames(modality_cols, modality_labels),
      name = "Therapeutic\nmodality", drop = FALSE) +
    scale_y_continuous(labels = scales::percent,
                       expand = expansion(mult = c(0, 0.02))) +
    labs(x = "Fibrosis stage", y = "Proportion of patients") +
    theme_masld(base_size = 7) +
    theme(legend.key.size = unit(0.25, "cm"))
  p
}, error = function(e) {
  message("Panel (c) error: ", conditionMessage(e))
  placeholder("Panel c: Modality\n(error)")
})

# ══════════════════════════════════════════════════════════════════════════════
# (d) Top personalized targets dot plot (frequency × mean score per transition)
# ══════════════════════════════════════════════════════════════════════════════
pd <- tryCatch({
  # Keep only fibrosis transitions
  tg <- top10[grepl("^F[0-3]_to_F[1-4]$", assigned_transition)]

  # Top 5 most frequent genes overall
  gene_freq <- tg[, .N, by = gene_symbol][order(-N)]
  top5 <- head(gene_freq$gene_symbol, 5)
  tg5 <- tg[gene_symbol %in% top5]

  # Per-gene per-transition: frequency in top-10 lists and mean score
  dot_dat <- tg5[, .(freq = .N, mean_score = mean(personalized_score, na.rm = TRUE)),
                 by = .(gene_symbol, assigned_transition)]

  # Clean labels
  dot_dat[, transition := gsub("_to_", " -> ", assigned_transition)]
  trans_order <- c("F0 -> F1", "F1 -> F2", "F2 -> F3", "F3 -> F4")
  dot_dat[, transition := factor(transition, levels = trans_order)]

  # Gene order by total frequency
  gene_order <- tg5[, .N, by = gene_symbol][order(N)]$gene_symbol
  dot_dat[, gene_symbol := factor(gene_symbol, levels = gene_order)]

  p <- ggplot(dot_dat, aes(x = transition, y = gene_symbol)) +
    geom_point(aes(size = freq, color = mean_score)) +
    scale_size_continuous(range = c(1, 5), name = "Patients\n(top-10)") +
    scale_color_gradient(low = "#42A5F5", high = "#C2185B",
                         name = "Mean\nscore") +
    labs(x = "Fibrosis transition", y = NULL) +
    theme_masld(base_size = 7) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
          panel.grid.major = element_line(color = "grey92", linewidth = 0.2))
  p
}, error = function(e) {
  message("Panel (d) error: ", conditionMessage(e))
  placeholder("Panel d: Top targets\n(error)")
})

# ══════════════════════════════════════════════════════════════════════════════
# Assemble and save
# ══════════════════════════════════════════════════════════════════════════════
combined <- (pa + pb) / (pc + pd) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 9, face = "bold"))

width_in  <- 180 / 25.4
height_in <- 120 / 25.4

pdf_dev <- if (capabilities("cairo")) cairo_pdf else grDevices::pdf
ggsave(file.path(OUT, "fig_targeting.pdf"), combined,
       width = width_in, height = height_in, device = pdf_dev)
message("Saved: ", file.path(OUT, "fig_targeting.pdf"))

# Individual panels
ggsave(file.path(OUT, "panel_targeting_a.pdf"), pa,
       width = width_in / 2, height = height_in / 2, device = pdf_dev)
ggsave(file.path(OUT, "panel_targeting_b.pdf"), pb,
       width = width_in / 2, height = height_in / 2, device = pdf_dev)
ggsave(file.path(OUT, "panel_targeting_c.pdf"), pc,
       width = width_in / 2, height = height_in / 2, device = pdf_dev)
ggsave(file.path(OUT, "panel_targeting_d.pdf"), pd,
       width = width_in / 2, height = height_in / 2, device = pdf_dev)
message("Individual panels saved to: ", OUT)

message("Done.")
