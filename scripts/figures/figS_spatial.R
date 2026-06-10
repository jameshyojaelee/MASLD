#!/usr/bin/env Rscript
# ==========================================================================
# Supplementary Figure: Spatial transcriptomics validation of MASLD programs
# 6-panel layout (3 rows x 2 columns):
#   a: Hepatocyte-intrinsic spatial validation scatter (bulk vs spatial)
#   b: Cell type proportions per sample (cell2location)
#   c: Spatial enrichment forest plot (gene sets among SVGs)
#   d: Differential SVGs — Moran's I healthy vs steatotic
#   e: DEG zonation classification (bar chart)
#   f: Disease-emergent L-R communication (dot plot)
# ==========================================================================

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

OUT <- file.path(FIGS05_DIR, "figS_spatial.pdf")
dir.create(file.path(FIGS05_DIR, "panels"), showWarnings = FALSE, recursive = TRUE)

# Spatial results directory
SPATIAL <- file.path(BASE, "Analysis/Spatial/results")

# cell2location column prefix to strip
C2L_PREFIX <- "q05cell_abundance_w_sf_means_per_cluster_mu_fg_"

# Condition colors
condition_colors <- c(Healthy = masld_colors$control, Steatotic = masld_colors$up)

# Cell type palette (ordered by expected abundance)
celltype_colors <- c(
  Hepatocytes           = "#0D47A1",
  Fibroblasts           = "#C2185B",
  `Endothelial cells`   = "#7B1FA2",
  Cholangiocytes        = "#00695C",
  `T cells`             = "#42A5F5",
  Macrophages           = "#E91E63",
  Neutrophils           = "#F48FB1",
  `B cells`             = "#AD1457",
  `Mono+mono derived cells` = "#880E4F",
  `Resident NK`         = "#1A237E",
  `Circulating NK/NKT`  = "#F57F17",
  `Plasma cells`        = "#64B5F6",
  cDC1s                 = "#90CAF9",
  cDC2s                 = "#BDBDBD",
  pDCs                  = "#757575",
  Basophils             = "#E0E0E0",
  Other                 = "#9E9E9E"
)

# Pre-initialize placeholders
p_a <- placeholder("Panel a: hep_intrinsic_validation.csv not found")
p_b <- placeholder("Panel b: spatial_cell_type_proportions.csv not found")
p_c <- placeholder("Panel c: spatial_enrichment_tests.csv not found")
p_d <- placeholder("Panel d: differential_svgs.csv not found")
p_e <- placeholder("Panel e: deg_zonation_classification.csv not found")
p_f <- placeholder("Panel f: differential_lr_pairs.csv not found")

# ==========================================================================
# Panel a: Hepatocyte-intrinsic spatial validation scatter
# X = bulk attribution score, Y = spatial fold change, color = validated
# ==========================================================================
hep_val_path <- file.path(SPATIAL, "cell2location/hep_intrinsic_validation.csv")
if (file.exists(hep_val_path)) {
  hep_val <- fread(hep_val_path)

  n_total <- nrow(hep_val)
  n_valid <- sum(hep_val$validated == TRUE | hep_val$validated == "True")
  pct_valid <- round(n_valid / n_total * 100, 1)

  # Compute Spearman correlation for annotation
  cor_test <- cor.test(hep_val$attribution_raw, hep_val$spatial_fc,
                       method = "spearman", exact = FALSE)

  # Convert validated to factor for color mapping
  hep_val[, validated := factor(
    fifelse(validated == TRUE | validated == "True", "Validated", "Not validated"),
    levels = c("Validated", "Not validated")
  )]

  # Label top validated genes (highest spatial_fc)
  top_genes <- hep_val[validated == "Validated"][order(-spatial_fc)][1:min(10, sum(hep_val$validated == "Validated"))]

  val_colors <- c(Validated = masld_colors$hep_intrinsic, `Not validated` = masld_colors$ns)

  p_a <- ggplot(hep_val, aes(x = attribution_raw, y = spatial_fc, color = validated)) +
    geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.3, color = "gray50") +
    rasterize_layer(
      geom_point(size = 0.8, alpha = 0.6, shape = 16)
    ) +
    geom_smooth(data = hep_val, aes(x = attribution_raw, y = spatial_fc),
                method = "lm", se = TRUE, linewidth = 0.4,
                color = "gray30", fill = "gray80", alpha = 0.3,
                inherit.aes = FALSE) +
    geom_label_repel(data = top_genes, aes(label = symbol),
                     size = 1.6, max.overlaps = 20,
                     label.padding = 0.12, box.padding = 0.4,
                     segment.size = 0.15, show.legend = FALSE) +
    scale_color_manual(values = val_colors,
                       name = NULL,
                       labels = c(
                         Validated = paste0("Validated (n=", n_valid, ")"),
                         `Not validated` = paste0("Not validated (n=", n_total - n_valid, ")")
                       )) +
    annotate("text", x = Inf, y = Inf,
             label = paste0("rho = ", sprintf("%.3f", cor_test$estimate),
                            ", P = ", sprintf("%.1e", cor_test$p.value)),
             hjust = 1.1, vjust = 1.5, size = 2.2) +
    annotate("text", x = Inf, y = Inf,
             label = paste0(pct_valid, "% validated"),
             hjust = 1.1, vjust = 3.2, size = 2.2, fontface = "bold") +
    labs(x = "Bulk attribution score",
         y = "Spatial FC (hep / non-hep spots)",
         title = "Hepatocyte-intrinsic spatial validation") +
    theme_masld() +
    guides(color = guide_legend(override.aes = list(size = 1.5, alpha = 1)))
}

# ==========================================================================
# Panel b: Cell type proportions per sample (stacked bar)
# ==========================================================================
prop_path <- file.path(SPATIAL, "cell2location/spatial_cell_type_proportions.csv")
if (file.exists(prop_path)) {
  props <- fread(prop_path)

  # Strip c2l prefix from column names
  old_names <- names(props)
  new_names <- gsub(C2L_PREFIX, "", old_names)
  setnames(props, old_names, new_names)

  # Assign condition from the canonical GSE192741 registry (F228/F247). The
  # previous map was INVERTED and named non-existent samples: canonical labels
  # (config/spatial_datasets.yaml + zonation_scores.csv) are
  #   Healthy  = {JBO018, JBO022}
  #   Steatotic= {JBO014, JBO015, JBO019}
  # (JBO016/JBO020 do not exist in this dataset). Define condition off the
  # Healthy set positively so unseen IDs do not silently default to Steatotic.
  props[, condition := fifelse(sample_id %in% c("JBO018", "JBO022"),
                               "Healthy", "Steatotic")]

  # Melt to long format
  ct_cols <- setdiff(names(props), c("sample_id", "condition"))
  props_long <- melt(props, id.vars = c("sample_id", "condition"),
                     measure.vars = ct_cols,
                     variable.name = "cell_type", value.name = "proportion")

  # Collapse rare cell types (<2% mean proportion) into "Other"
  ct_means <- props_long[, .(mean_prop = mean(proportion)), by = cell_type]
  top_cts <- ct_means[mean_prop >= 0.02]$cell_type
  props_long[!cell_type %in% top_cts, cell_type := "Other"]
  props_long <- props_long[, .(proportion = sum(proportion)), by = .(sample_id, condition, cell_type)]

  # Order cell types by mean proportion (largest first)
  ct_order <- props_long[, .(mean_prop = mean(proportion)), by = cell_type][order(-mean_prop)]$cell_type
  # Ensure "Other" is at bottom of stack
  ct_order <- c(setdiff(ct_order, "Other"), "Other")
  props_long[, cell_type := factor(cell_type, levels = rev(ct_order))]

  # Order samples: Healthy first, then Steatotic
  sample_order <- props[order(condition, sample_id)]$sample_id
  props_long[, sample_id := factor(sample_id, levels = sample_order)]

  p_b <- ggplot(props_long, aes(x = sample_id, y = proportion, fill = cell_type)) +
    geom_col(width = 0.75, color = "white", linewidth = 0.15) +
    scale_fill_manual(values = celltype_colors, name = "Cell type") +
    scale_y_continuous(labels = percent_format(), expand = expansion(mult = c(0, 0.02))) +
    labs(x = NULL, y = "Proportion",
         title = "Cell type composition (cell2location)") +
    theme_masld() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 6),
          legend.key.size = unit(0.25, "cm"),
          legend.text = element_text(size = 5)) +
    guides(fill = guide_legend(ncol = 1, reverse = TRUE))

  # Add condition bracket annotation (below x-axis). Samples are ordered
  # Healthy-first (n=2: JBO018/JBO022) then Steatotic (n=3: JBO014/JBO015/JBO019)
  # by the order(condition, sample_id) sort above, so the bracket split is at 2.5
  # (was 3.5 under the inverted 3-Healthy/2-Steatotic assumption; F228/F247).
  p_b <- p_b +
    annotate("segment", x = 0.5, xend = 2.5, y = -0.04, yend = -0.04,
             color = masld_colors$control, linewidth = 0.6) +
    annotate("text", x = 1.5, y = -0.065, label = "Healthy",
             color = masld_colors$control, size = 2, fontface = "bold") +
    annotate("segment", x = 2.5, xend = 5.5, y = -0.04, yend = -0.04,
             color = masld_colors$up, linewidth = 0.6) +
    annotate("text", x = 4, y = -0.065, label = "Steatotic",
             color = masld_colors$up, size = 2, fontface = "bold") +
    coord_cartesian(clip = "off") +
    theme(plot.margin = margin(5, 5, 20, 5))
}

# ==========================================================================
# Panel c: SVG enrichment forest plot (odds ratios + 95% CI)
# ==========================================================================
enrich_path <- file.path(SPATIAL, "integration/spatial_enrichment_tests.csv")
if (file.exists(enrich_path)) {
  enrich <- fread(enrich_path)

  # Filter to SVG-level tests (not zonation subtests)
  # Keep rows where gene_set is one of the main sets
  main_sets <- c("Dream_DEGs", "Conserved", "Hepatocyte_intrinsic",
                 "Sex_dimorphic_DEGs", "Druggable_genes", "Top5pct_multi_evidence")
  enrich_svg <- enrich[gene_set %in% main_sets]

  if (nrow(enrich_svg) > 0) {
    # 95% CI for log(OR) via Woolf's method on the REAL 2x2 cell counts (F202/F256).
    # The previous code back-solved SE from the BH-adjusted p-value
    # (z = qnorm(1 - padj/2); SE = |log OR|/z), which is circular: that CI always
    # brackets OR=1 exactly at the significance boundary by construction and is
    # neither a Wald nor a Fisher interval. The producer (06_integration.py:381-386)
    # builds the 2x2 from set-membership; we reconstruct the three cells the CSV
    # carries: a = n_overlap (svg & set); b = n_svg - n_overlap (svg not in set);
    # c = n_set - n_overlap (set not svg). The 4th cell d (universe - svg - set +
    # overlap) is not persisted; d ~ thousands so 1/d is negligible and dropping it
    # only widens the interval (conservative). This is a genuine Woolf SE from the
    # cell counts, not a value reconstructed from the p-value.
    enrich_svg[, log_or := log(odds_ratio)]
    enrich_svg[, a_cell := n_overlap]
    enrich_svg[, b_cell := pmax(n_svg - n_overlap, 0)]
    enrich_svg[, c_cell := pmax(n_set  - n_overlap, 0)]
    # Woolf SE = sqrt(1/a + 1/b + 1/c [+ 1/d]); omit the negligible 1/d term.
    enrich_svg[, se_woolf := sqrt(1 / pmax(a_cell, 0.5) +
                                  1 / pmax(b_cell, 0.5) +
                                  1 / pmax(c_cell, 0.5))]
    enrich_svg[, ci_lo := exp(log_or - 1.96 * se_woolf)]
    enrich_svg[, ci_hi := exp(log_or + 1.96 * se_woolf)]

    # Significance label
    enrich_svg[, sig := fifelse(fisher_padj_bh < 0.05, "FDR < 0.05", "n.s.")]

    # Clean display names
    enrich_svg[, display := gsub("_", " ", gene_set)]
    enrich_svg[, display := factor(display, levels = rev(display))]

    sig_colors <- c(`FDR < 0.05` = masld_colors$up, n.s. = masld_colors$ns)

    p_c <- ggplot(enrich_svg, aes(x = odds_ratio, y = display, color = sig)) +
      geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.3, color = "gray50") +
      geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi), width = 0.2, linewidth = 0.4,
                    orientation = "y") +
      geom_point(size = 2, shape = 16) +
      geom_text(aes(label = paste0("n=", n_overlap)),
                hjust = -0.3, vjust = -0.6, size = 1.8, show.legend = FALSE) +
      scale_color_manual(values = sig_colors, name = NULL) +
      scale_x_log10() +
      labs(x = "Odds ratio (log scale)",
           y = NULL,
           title = "Gene set enrichment among SVGs") +
      theme_masld() +
      theme(legend.position = c(0.8, 0.15),
            legend.background = element_rect(fill = "white", color = NA))
  }
}

# ==========================================================================
# Panel d: Differential SVGs — Moran's I scatter
# ==========================================================================
dsvg_path <- file.path(SPATIAL, "svg/differential_svgs.csv")
if (file.exists(dsvg_path)) {
  dsvg <- fread(dsvg_path)

  # Fix column: first column is gene name (unnamed row index in CSV)
  if (names(dsvg)[1] == "V1" || names(dsvg)[1] == "") {
    setnames(dsvg, 1, "gene")
  } else if (!("gene" %in% names(dsvg))) {
    setnames(dsvg, 1, "gene")
  }

  # Ensure boolean columns are logical
  for (col in c("svg_healthy", "svg_masld")) {
    if (is.character(dsvg[[col]])) dsvg[, (col) := get(col) == "True"]
  }

  # Categorize (actual CSV values: disease_emergent_SVG, disease_lost_SVG, stable)
  dsvg[, display_cat := fcase(
    grepl("emergent", category, ignore.case = TRUE), "Disease-emergent",
    grepl("lost|resolved", category, ignore.case = TRUE), "Disease-resolved",
    category == "stable",  "Stable SVG",
    default = "Not SVG"
  )]

  # Keep only genes that are SVG in at least one condition
  dsvg_plot <- dsvg[svg_healthy == TRUE | svg_masld == TRUE]

  # Order: Not SVG / Stable at bottom, emergent/resolved on top
  dsvg_plot[, plot_order := fifelse(display_cat %in% c("Disease-emergent", "Disease-resolved"), 1L, 0L)]
  setorder(dsvg_plot, plot_order)

  cat_colors <- c(
    `Disease-emergent` = masld_colors$up,
    `Disease-resolved` = masld_colors$down,
    `Stable SVG`       = masld_colors$ns
  )

  # Label top disease-emergent genes (largest delta_I)
  top_emergent <- dsvg_plot[display_cat == "Disease-emergent"][order(-delta_I)][1:min(8, .N)]
  top_resolved <- dsvg_plot[display_cat == "Disease-resolved"][order(delta_I)][1:min(5, .N)]
  top_labels_d <- rbind(top_emergent, top_resolved)

  # Count categories for legend
  cat_n <- dsvg_plot[, .N, by = display_cat]

  p_d <- ggplot(dsvg_plot, aes(x = morans_I_healthy, y = morans_I_masld, color = display_cat)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                linewidth = 0.3, color = "gray50") +
    rasterize_layer(
      geom_point(size = 0.5, alpha = 0.5, shape = 16)
    ) +
    geom_label_repel(data = top_labels_d, aes(label = gene),
                     size = 1.5, max.overlaps = 30,
                     label.padding = 0.1, box.padding = 0.3,
                     segment.size = 0.12, show.legend = FALSE) +
    scale_color_manual(
      values = cat_colors, name = NULL,
      labels = setNames(
        paste0(cat_n$display_cat, " (n=", cat_n$N, ")"),
        cat_n$display_cat
      )
    ) +
    labs(x = expression("Moran's " * italic(I) * " (Healthy)"),
         y = expression("Moran's " * italic(I) * " (Steatotic)"),
         title = "Spatially variable gene dynamics") +
    theme_masld() +
    guides(color = guide_legend(override.aes = list(size = 1.5, alpha = 1)))
}

# ==========================================================================
# Panel e: DEG zonation classification (horizontal bar + Conserved)
# ==========================================================================
zon_path <- file.path(SPATIAL, "zonation/deg_zonation_classification.csv")
if (file.exists(zon_path)) {
  zon <- fread(zon_path)
  n_zon <- nrow(zon)

  # Load Conserved gene list for cross-reference
  conc <- load_concordance_atlas()
  if (!is.null(conc)) {
    core_genes <- conc[primary_category == "Conserved", human_symbol]
    zon[, is_core := gene %in% core_genes]
  } else {
    zon[, is_core := FALSE]
  }

  # Summary counts
  zon_summary <- zon[, .(All = .N, `Conserved` = sum(is_core)), by = zonation_class]
  zon_melt <- melt(zon_summary, id.vars = "zonation_class",
                   variable.name = "subset", value.name = "count")

  # Use clean display names (no newlines — let ggplot2 wrap)
  zon_melt[, zonation_class := factor(zonation_class,
           levels = rev(c("Pan-lobular", "Periportal-enriched", "Pericentral-enriched")))]

  zon_class_colors <- c(
    `Pan-lobular`          = masld_colors$ns,
    `Periportal-enriched`  = masld_colors$down,
    `Pericentral-enriched` = masld_colors$up
  )

  subset_alpha <- c(All = 1, `Conserved` = 0.5)

  p_e <- ggplot(zon_melt[subset == "All"],
                aes(x = count, y = zonation_class, fill = zonation_class)) +
    geom_col(width = 0.6) +
    geom_col(data = zon_melt[subset == "Conserved"],
             aes(x = count, y = zonation_class),
             width = 0.6, fill = masld_colors$conserved, alpha = 0.8) +
    geom_text(aes(label = paste0(comma(count), " (",
                                  sprintf("%.1f", count / n_zon * 100), "%)")),
              hjust = -0.05, size = 2, color = "gray20") +
    geom_text(data = zon_melt[subset == "Conserved" & count > 0],
              aes(x = count, label = paste0(count, " Core")),
              hjust = -0.05, size = 1.8, color = masld_colors$conserved,
              fontface = "bold") +
    scale_fill_manual(values = zon_class_colors, guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.35))) +
    labs(x = paste0("DEGs (N=", comma(n_zon), " classified)"),
         y = NULL,
         title = "Zonation classification of MASLD DEGs") +
    theme_masld()
}

# ==========================================================================
# Panel f: Disease-emergent L-R communication (dot plot)
# ==========================================================================
dlr_path <- file.path(SPATIAL, "communication/differential_lr_pairs.csv")
if (file.exists(dlr_path)) {
  dlr <- fread(dlr_path)

  # Focus on disease-emergent pairs with highest effect sizes
  emergent <- dlr[category == "disease_emergent"]

  if (nrow(emergent) > 0) {
    # Clean lr_pair (remove tuple formatting)
    emergent[, lr_clean := gsub("^\\('|'\\)$", "", lr_pair)]
    emergent[, lr_clean := gsub("', '", " \u2192 ", lr_clean)]

    # Create axis label: source -> target: lr_pair
    emergent[, display := paste0(source, " \u2192 ", target, ": ", lr_clean)]

    # Top 20 by delta_mean_expr
    top_lr <- emergent[order(-delta_mean_expr)][1:min(20, .N)]
    top_lr[, display := factor(display, levels = rev(display))]

    pair_colors <- c(autocrine = "#7B1FA2", paracrine = masld_colors$up)

    p_f <- ggplot(top_lr, aes(x = delta_mean_expr, y = display,
                               color = pair_type, size = mean_expr_disease)) +
      geom_segment(aes(x = 0, xend = delta_mean_expr, y = display, yend = display),
                   linewidth = 0.3, color = "gray70") +
      geom_point(shape = 16) +
      scale_color_manual(values = pair_colors, name = "Interaction") +
      scale_size_continuous(range = c(1, 3.5), name = "Mean expr\n(disease)") +
      labs(x = expression(Delta * " mean expression"),
           y = NULL,
           title = "Top disease-emergent L-R pairs") +
      theme_masld() +
      theme(axis.text.y = element_text(size = 4.5))
  }
}

# ==========================================================================
# Assemble: 3 rows x 2 columns
# ==========================================================================
fig_spatial <- (p_a | p_b) / (p_c | p_d) / (p_e | p_f) +
  plot_annotation(
    title = "Spatial transcriptomics validation (Guilliams et al.)",
    tag_levels = "a"
  ) &
  theme(plot.tag = element_text(size = 8, face = "bold"))

save_fig_tall(fig_spatial, OUT, height = 10)
message("Spatial supplementary figure saved to ", OUT)
