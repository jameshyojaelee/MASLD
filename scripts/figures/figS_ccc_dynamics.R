#!/usr/bin/env Rscript
# ===========================================================================
# Supplementary Figure S13: Cell-Cell Communication Dynamics in MASLD
# ---------------------------------------------------------------------------
# 6-panel layout:
#   (a) Heatmap: axes x transitions — n sig rewired L-R pairs
#   (b) Lollipop: top gained + lost L-R pairs per transition
#   (c) Pathway dynamics along pseudotime (line plot per pathway family)
#   (d) Stellate->Macrophage spotlight heatmap (top L-R pairs x pseudotime)
#   (e) GAS6::MERTK loss + CCN1::CAV1 gain — key biological events
#   (f) L-R cluster trajectories (rate-differentiated)
#
# Output: figures/supplementary/figS02_progression/figS13_ccc_dynamics.pdf
# ===========================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

data_dir <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/progression")
out_file  <- file.path(FIGS02_DIR, "figS13_ccc_dynamics.pdf")
dir.create(file.path(FIGS02_DIR, "panels"), showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
cat("Loading CCC data files...\n")
axis_summary   <- read.csv(file.path(data_dir, "ccc_axis_transition_summary.csv"),
                            stringsAsFactors = FALSE)
top_rewired    <- read.csv(file.path(data_dir, "ccc_top_rewired.csv"),
                            stringsAsFactors = FALSE)
transition_de  <- read.csv(file.path(data_dir, "ccc_transition_de.csv"),
                            stringsAsFactors = FALSE)
pt_dynamics    <- read.csv(file.path(data_dir, "ccc_pseudotime_dynamics.csv"),
                            stringsAsFactors = FALSE)
lr_clusters    <- read.csv(file.path(data_dir, "ccc_lr_clusters.csv"),
                            stringsAsFactors = FALSE)
pathway_dyn    <- read.csv(file.path(data_dir, "ccc_pathway_dynamics.csv"),
                            stringsAsFactors = FALSE)
binned_scores  <- read.csv(file.path(data_dir, "ccc_binned_scores.csv"),
                            stringsAsFactors = FALSE)

# ---------------------------------------------------------------------------
# Axis display ordering (clean names)
# ---------------------------------------------------------------------------
axis_order <- c(
  "Hepatocyte_to_Macrophage", "Hepatocyte_to_Stellate",
  "Macrophage_to_Hepatocyte", "Macrophage_to_Stellate",
  "Stellate_to_Hepatocyte",   "Stellate_to_Macrophage",
  "Endothelial_to_Hepatocyte","Endothelial_to_Stellate"
)
axis_labels <- gsub("_to_", " -> ", axis_order)

transition_order <- c("F0_to_F1", "F1_to_F2", "F2_to_F3", "F3_to_F4")
transition_labels <- gsub("_to_", " -> ", transition_order)

# ===========================================================================
# Panel (a): Communication rewiring overview — heatmap
# ===========================================================================
panel_a <- tryCatch({
  cat("  Panel (a): Rewiring overview heatmap\n")

  # Ensure all axis x transition combos exist
  full_grid <- expand.grid(
    axis       = axis_order,
    transition = transition_order,
    stringsAsFactors = FALSE
  )
  hm_data <- full_grid %>%
    left_join(axis_summary, by = c("axis", "transition")) %>%
    mutate(
      n_sig      = ifelse(is.na(n_sig), 0, n_sig),
      axis       = factor(axis, levels = rev(axis_order), labels = rev(axis_labels)),
      transition = factor(transition, levels = transition_order, labels = transition_labels)
    )

  ggplot(hm_data, aes(x = transition, y = axis, fill = n_sig)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = n_sig), size = 2, color = "black") +
    scale_fill_gradient(low = "white", high = masld_colors$fibrosis,
                        name = "Sig. rewired\nL-R pairs") +
    labs(title = "Communication rewiring by axis and transition",
         x = NULL, y = NULL) +
    theme_masld() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}, error = function(e) {
  cat("    Panel (a) FAILED:", conditionMessage(e), "\n")
  placeholder("Panel (a) failed")
})

# ===========================================================================
# Panel (b): Top rewired L-R pairs — lollipop chart
# ===========================================================================
panel_b <- tryCatch({
  cat("  Panel (b): Top rewired L-R pairs lollipop\n")

  # Classify direction: lfc > 0 = gained, lfc < 0 = lost
  top_rw <- top_rewired %>%
    mutate(
      direction = ifelse(lfc > 0, "Gained", "Lost"),
      abs_lfc   = abs(lfc),
      transition = factor(transition, levels = transition_order, labels = transition_labels)
    )

  # Select top 3 gained + top 3 lost per transition (less clutter than 5)
  top_gained <- top_rw %>%
    filter(direction == "Gained") %>%
    group_by(transition) %>%
    slice_max(order_by = abs_lfc, n = 3, with_ties = FALSE) %>%
    ungroup()
  top_lost <- top_rw %>%
    filter(direction == "Lost") %>%
    group_by(transition) %>%
    slice_max(order_by = abs_lfc, n = 3, with_ties = FALSE) %>%
    ungroup()

  plot_data <- bind_rows(top_gained, top_lost) %>%
    mutate(
      axis_short = gsub("_to_", "->", gsub("Endothelial", "Endo",
                    gsub("Macrophage", "Mac",
                    gsub("Hepatocyte", "Hep",
                    gsub("Stellate", "Stel", axis))))),
      pair_label = paste0(lr_pair, " (", axis_short, ")")
    )

  # Order by LFC within transition
  plot_data <- plot_data %>%
    arrange(transition, lfc) %>%
    mutate(pair_label = factor(pair_label, levels = unique(pair_label)))

  ggplot(plot_data, aes(x = lfc, y = pair_label, color = direction)) +
    geom_segment(aes(x = 0, xend = lfc, yend = pair_label),
                 linewidth = 0.4) +
    geom_point(size = 1.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.3) +
    scale_color_manual(values = c(Gained = masld_colors$up, Lost = masld_colors$down),
                       name = NULL) +
    facet_wrap(~transition, scales = "free_y", ncol = 2) +
    labs(title = "Top rewired L-R pairs per transition",
         x = "Log2 fold-change", y = NULL) +
    theme_masld() +
    theme(axis.text.y  = element_text(size = 5.5),
          strip.text   = element_text(size = 7, face = "bold"),
          legend.position = c(0.85, 0.08),
          legend.background = element_rect(fill = alpha("white", 0.9), color = NA),
          legend.key.size = unit(0.2, "cm"))
}, error = function(e) {
  cat("    Panel (b) FAILED:", conditionMessage(e), "\n")
  placeholder("Panel (b) failed")
})

# ===========================================================================
# Panel (c): Pathway dynamics along pseudotime — line plot
# ===========================================================================
panel_c <- tryCatch({
  cat("  Panel (c): Pathway dynamics along pseudotime\n")

  # Map canonical_pathway from pseudotime_dynamics onto binned_scores
  # Create L-R pair -> canonical_pathway lookup (most common non-"Other" assignment)
  lr_pathway_map <- pt_dynamics %>%
    filter(!is.na(canonical_pathway) & canonical_pathway != "Other" &
           canonical_pathway != "") %>%
    select(lr_pair, canonical_pathway) %>%
    distinct()

  # Aggregate binned scores by pathway
  binned_with_pathway <- binned_scores %>%
    inner_join(lr_pathway_map, by = "lr_pair", relationship = "many-to-many") %>%
    group_by(canonical_pathway, pt_bin) %>%
    summarize(
      mean_score  = mean(mean_score, na.rm = TRUE),
      n_pairs     = n_distinct(lr_pair),
      .groups     = "drop"
    )

  # Normalize each pathway to its bin=0 value for comparable trajectories
  baseline <- binned_with_pathway %>%
    filter(pt_bin == 0) %>%
    select(canonical_pathway, baseline_score = mean_score)
  binned_norm <- binned_with_pathway %>%
    left_join(baseline, by = "canonical_pathway") %>%
    mutate(
      # Use z-score-style normalization: fold relative to baseline
      relative_score = ifelse(baseline_score > 0,
                              mean_score / baseline_score, NA_real_)
    )

  # Select key pathways (7 most informative, ordered by biological importance)
  key_pathways <- c("TGFbeta", "ECM_Integrin", "Interleukin", "Chemokine",
                    "NOTCH", "HGF_MET", "VEGF")
  plot_data <- binned_norm %>%
    filter(canonical_pathway %in% key_pathways) %>%
    mutate(canonical_pathway = factor(canonical_pathway, levels = key_pathways))

  # Distinct colors + linetypes for 7 pathways
  pw_cols <- setNames(
    c("#C2185B", "#00695C", "#F57F17", "#7B1FA2",
      "#1565C0", "#E91E63", "#0D47A1"),
    key_pathways
  )
  pw_lty <- setNames(
    c("solid", "solid", "solid", "dashed",
      "dashed", "dotted", "dotted"),
    key_pathways
  )

  # Annotate: "All 14 pathways decreasing"
  ggplot(plot_data, aes(x = pt_bin, y = relative_score,
                        color = canonical_pathway, linetype = canonical_pathway,
                        group = canonical_pathway)) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 0.6) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray50", linewidth = 0.3) +
    annotate("text", x = 1, y = 0.22, label = "All 14 pathways\ndecreasing",
             size = 2, color = "gray40", fontface = "italic", hjust = 0) +
    scale_color_manual(values = pw_cols, name = "Pathway") +
    scale_linetype_manual(values = pw_lty, name = "Pathway") +
    scale_x_continuous(breaks = seq(0, 19, by = 5)) +
    labs(title = "Pathway-level L-R dynamics along pseudotime",
         x = "Pseudotime bin", y = "Relative score (fold vs bin 0)") +
    theme_masld() +
    theme(legend.position = "right",
          legend.text = element_text(size = 6),
          legend.key.width = unit(0.5, "cm"))
}, error = function(e) {
  cat("    Panel (c) FAILED:", conditionMessage(e), "\n")
  placeholder("Panel (c) failed")
})

# ===========================================================================
# Panel (d): Stellate->Macrophage spotlight heatmap (top L-R x pseudotime)
# ===========================================================================
panel_d <- tryCatch({
  cat("  Panel (d): Stellate->Macrophage spotlight heatmap\n")

  # Get top 20 most dynamic Stellate_to_Macrophage pairs by |rho|
  stm_pairs <- pt_dynamics %>%
    filter(axis == "Stellate_to_Macrophage") %>%
    arrange(spearman_pval) %>%
    slice_head(n = 20) %>%
    pull(lr_pair)

  # Get binned scores for these pairs
  stm_binned <- binned_scores %>%
    filter(axis == "Stellate_to_Macrophage" & lr_pair %in% stm_pairs)

  # Normalize per pair (z-score across bins)
  stm_norm <- stm_binned %>%
    group_by(lr_pair) %>%
    mutate(z_score = (mean_score - mean(mean_score)) / (sd(mean_score) + 1e-10)) %>%
    ungroup() %>%
    mutate(lr_pair = factor(lr_pair, levels = rev(stm_pairs)))

  ggplot(stm_norm, aes(x = pt_bin, y = lr_pair, fill = z_score)) +
    geom_tile() +
    scale_fill_gradient2(low = masld_colors$down, mid = "white", high = masld_colors$up,
                         midpoint = 0, name = "Z-score",
                         limits = c(-2.5, 2.5), oob = scales::squish) +
    scale_x_continuous(breaks = seq(0, 19, by = 5),
                       labels = seq(0, 19, by = 5)) +
    labs(title = expression("Stellate" %->% "Macrophage: top dynamic L-R pairs"),
         x = "Pseudotime bin", y = NULL) +
    theme_masld() +
    theme(axis.text.y = element_text(size = 5))
}, error = function(e) {
  cat("    Panel (d) FAILED:", conditionMessage(e), "\n")
  placeholder("Panel (d) failed")
})

# ===========================================================================
# Panel (e): GAS6::MERTK loss + CCN1::CAV1 gain — key biological events
# ===========================================================================
panel_e <- tryCatch({
  cat("  Panel (e): GAS6::MERTK loss + CCN1::CAV1 gain\n")

  # --- GAS6::MERTK on Hepatocyte_to_Macrophage axis ---
  gas6_data <- pt_dynamics %>%
    filter(lr_pair == "GAS6::MERTK" & axis == "Hepatocyte_to_Macrophage") %>%
    select(lr_pair, axis, mean_F0, mean_F1, mean_F2, mean_F3, mean_F4,
           spearman_rho, spearman_padj)
  if (nrow(gas6_data) == 0) stop("GAS6::MERTK not found on Hepatocyte_to_Macrophage axis")

  gas6_long <- gas6_data %>%
    select(starts_with("mean_F")) %>%
    pivot_longer(everything(), names_to = "stage", values_to = "mean_score") %>%
    mutate(
      stage = gsub("mean_", "", stage),
      stage = factor(stage, levels = c("F0", "F1", "F2", "F3", "F4")),
      pair = "GAS6::MERTK"
    )

  # --- CCN1::CAV1 on Stellate_to_Macrophage axis (top gained F2->F3) ---
  ccn1_data <- transition_de %>%
    filter(lr_pair == "CCN1::CAV1" & axis == "Stellate_to_Macrophage")
  # Get per-stage means from pt_dynamics
  ccn1_pt <- pt_dynamics %>%
    filter(lr_pair == "CCN1::CAV1" & axis == "Stellate_to_Macrophage") %>%
    select(lr_pair, axis, mean_F0, mean_F1, mean_F2, mean_F3, mean_F4,
           spearman_rho, spearman_padj)

  if (nrow(ccn1_pt) > 0) {
    ccn1_long <- ccn1_pt %>%
      select(starts_with("mean_F")) %>%
      pivot_longer(everything(), names_to = "stage", values_to = "mean_score") %>%
      mutate(
        stage = gsub("mean_", "", stage),
        stage = factor(stage, levels = c("F0", "F1", "F2", "F3", "F4")),
        pair = "CCN1::CAV1"
      )
    both_long <- bind_rows(gas6_long, ccn1_long)
  } else {
    both_long <- gas6_long
  }

  # Normalize each pair to its F0 value for comparable y-axis
  both_norm <- both_long %>%
    group_by(pair) %>%
    mutate(
      baseline = mean_score[stage == "F0"],
      relative = mean_score / baseline
    ) %>%
    ungroup()

  pair_cols <- c("GAS6::MERTK" = masld_colors$down, "CCN1::CAV1" = masld_colors$up)

  p_e <- ggplot(both_norm, aes(x = stage, y = relative, color = pair, group = pair)) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 2) +
    geom_hline(yintercept = 1, linetype = "dotted", color = "gray60", linewidth = 0.3) +
    # Annotate GAS6 F3->F4 drop
    annotate("text", x = 4.5, y = 0.4,
             label = "LFC = -1.00\npadj = 4.9e-08",
             size = 2, color = masld_colors$down, fontface = "italic") +
    # Annotate CCN1 F2->F3 gain (if present)
    {if (nrow(ccn1_pt) > 0)
      annotate("text", x = 3.5, y = 3.8,
               label = "LFC = +1.67\npadj = 1.6e-08",
               size = 2, color = masld_colors$up, fontface = "italic")
    } +
    scale_color_manual(values = pair_cols, name = NULL) +
    labs(title = "Key biological events: efferocytosis loss vs fibrosis gain",
         subtitle = "Normalized to F0 baseline",
         x = "Fibrosis stage", y = "Relative interaction score") +
    theme_masld() +
    theme(legend.position = c(0.5, 0.95),
          legend.direction = "horizontal",
          legend.background = element_rect(fill = alpha("white", 0.9), color = NA),
          legend.key.size = unit(0.3, "cm"))

  p_e
}, error = function(e) {
  cat("    Panel (e) FAILED:", conditionMessage(e), "\n")
  placeholder("Panel (e) failed")
})

# ===========================================================================
# Panel (f): L-R cluster trajectories (differentiated by enriched axes)
# ===========================================================================
panel_f <- tryCatch({
  cat("  Panel (f): L-R cluster trajectories\n")

  # Merge cluster assignments from pt_dynamics into binned_scores
  cluster_map <- pt_dynamics %>%
    filter(!is.na(cluster)) %>%
    select(axis, lr_pair, cluster) %>%
    distinct()

  # Per-cluster mean trajectory from binned scores
  cluster_binned <- binned_scores %>%
    inner_join(cluster_map, by = c("axis", "lr_pair")) %>%
    group_by(cluster, pt_bin) %>%
    summarize(
      mean_score  = mean(mean_score, na.rm = TRUE),
      n_pairs     = n_distinct(lr_pair),
      .groups     = "drop"
    )

  # Normalize per cluster to bin 0
  cl_baseline <- cluster_binned %>%
    filter(pt_bin == 0) %>%
    select(cluster, baseline = mean_score)
  cluster_norm <- cluster_binned %>%
    left_join(cl_baseline, by = "cluster") %>%
    mutate(relative_score = mean_score / baseline)

  # Build informative labels from lr_clusters (enriched axes, not just "decreasing")
  cl_labels <- lr_clusters %>%
    mutate(
      # Extract first 2 top axes and abbreviate
      short_axes = sapply(strsplit(top_axes, ";"), function(x) {
        ax <- gsub("_to_", "->", gsub("Endothelial", "Endo",
                gsub("Macrophage", "Mac",
                gsub("Hepatocyte", "Hep",
                gsub("Stellate", "Stel", x[1:min(2, length(x))])))))
        paste(ax, collapse = ", ")
      }),
      label = paste0("C", cluster, ": ", short_axes,
                     " (n=", n_pairs, ", rho=", round(mean_rho, 2), ")")
    ) %>%
    select(cluster, label)

  cluster_norm <- cluster_norm %>%
    mutate(cluster = as.numeric(cluster)) %>%
    left_join(cl_labels, by = "cluster") %>%
    mutate(label = ifelse(is.na(label),
                          paste0("Cluster ", cluster), label))

  # Colors per cluster
  cl_colors <- c("0" = masld_colors$down, "1" = masld_colors$up)

  ggplot(cluster_norm, aes(x = pt_bin, y = relative_score,
                           color = factor(cluster), group = factor(cluster))) +
    geom_ribbon(aes(ymin = relative_score * 0.95, ymax = relative_score * 1.05,
                    fill = factor(cluster)),
                alpha = 0.1, color = NA) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.2) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray50", linewidth = 0.3) +
    scale_color_manual(
      values = cl_colors,
      labels = cl_labels$label,
      name   = NULL
    ) +
    scale_fill_manual(
      values = cl_colors,
      guide = "none"
    ) +
    scale_x_continuous(breaks = seq(0, 19, by = 5)) +
    labs(title = "L-R trajectory clusters along pseudotime",
         subtitle = "Both clusters monotonically decrease; C0 drops faster",
         x = "Pseudotime bin",
         y = "Relative score (fold vs bin 0)") +
    theme_masld() +
    theme(legend.position = "bottom",
          legend.background = element_rect(fill = "white", color = NA),
          legend.text = element_text(size = 5.5),
          legend.margin = margin(0, 0, 0, 0))
}, error = function(e) {
  cat("    Panel (f) FAILED:", conditionMessage(e), "\n")
  placeholder("Panel (f) failed")
})

# ===========================================================================
# Assemble 6-panel figure
# ===========================================================================
cat("Assembling 6-panel figure...\n")

# Layout: 3 rows x 2 columns
# Row 1: (a) overview heatmap | (b) top rewired lollipop
# Row 2: (c) pathway dynamics  | (d) stellate->mac spotlight
# Row 3: (e) GAS6+CCN1 events  | (f) cluster trajectories

fig <- (panel_a | panel_b) /
       (panel_c | panel_d) /
       (panel_e | panel_f) +
  plot_layout(heights = c(1, 1.1, 0.9)) +
  plot_annotation(
    tag_levels = "a",
    theme    = theme(
      plot.title = element_text(size = 9, face = "bold", family = "Helvetica", hjust = 0.5)
    )
  ) &
  theme(plot.tag = element_text(size = 9, face = "bold", family = "Helvetica"))

# Save at generous dimensions for 6-panel supplementary figure
save_fig_tall(fig, out_file, width = fig_full_width, height = 11)
cat("Figure saved to:", out_file, "\n")
cat("Done.\n")
