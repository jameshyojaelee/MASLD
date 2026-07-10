#!/usr/bin/env Rscript
# ==========================================================================
# Supplementary Figure 7: Extended Cell-Type Composition & Attribution
# 5-panel layout (3 rows):
#   Row 1: a (Hep+Mac fractions) | b (attribution examples)
#   Row 2: c (all non-parenchymal cell types, Control vs Disease)
#   Row 3: d (stacked bar composition) | e (per-dataset hepatocyte fractions)
# ==========================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(scales)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUT <- file.path(FIGS03_DIR, "supp_fig7_extended_celltype.pdf")

BP_DIR <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration",
                    "results/deconvolution/bayesprism")

# Pre-initialize placeholders
p_a <- placeholder("Panel a: cell-type fractions not found")
p_b <- placeholder("Panel b: data unavailable")
p_c <- placeholder("Panel c: full cell-type data not found")
p_d <- placeholder("Panel d: composition data not found")
p_e <- placeholder("Panel e: per-dataset data not found")

# ---------------------------------------------------------------------------
# Load full 16-cell-type BayesPrism proportions from per-cohort files
# ---------------------------------------------------------------------------
DECONV_DIR <- file.path(BASE, "Analysis/Deconvolution/results")
# 9 cohorts in the cohort presentation; PRJNA512027 (Gerhard 2018) excluded
# for L0/S0 library-prep batch confound with diagnosis.
datasets_10 <- c("GSE126848", "GSE130970", "GSE135251", "GSE162694", "GSE167523",
                  "GSE174478", "GSE193066", "GSE213621", "GSE240729")

full_bp <- rbindlist(lapply(datasets_10, function(ds) {
  f <- file.path(DECONV_DIR, ds, paste0(ds, "_bayesprism_proportions.tsv"))
  if (!file.exists(f)) return(NULL)
  # R write.table with row.names=TRUE: header has N cols, data has N+1 (rowname + values)
  raw <- read.table(f, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
  mat <- as.data.table(raw)
  mat[, sample_id := rownames(raw)]
  mat[, dataset := ds]
  mat
}), fill = TRUE)

# Load deconvolution data
deconv <- load_deconv_attribution()

# ==========================================================================
# Panel a: Cell-type fraction distributions across disease stages
# ==========================================================================
prop_path <- file.path(BP_DIR, "unified_bayesprism_proportions.csv")
meta_path <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration",
                       "metadata/unified_metadata.csv")

if (file.exists(prop_path) && file.exists(meta_path)) {
  props <- fread(prop_path)
  meta  <- fread(meta_path, select = c("sample_id", "group_binary", "condition"))

  props <- merge(props, meta, by = "sample_id", all.x = TRUE)
  props <- props[!is.na(group_binary)]

  hep_col <- intersect(c("music_Hepatocytes", "bp_Hepatocytes"), names(props))[1]
  mac_col <- intersect(c("music_Macrophages", "bp_Macrophages"), names(props))[1]

  if (!is.na(hep_col) && !is.na(mac_col)) {
    plot_props <- melt(props,
                       id.vars = c("sample_id", "dataset", "group_binary"),
                       measure.vars = c(hep_col, mac_col),
                       variable.name = "cell_type",
                       value.name = "proportion")
    plot_props[, cell_type := gsub("^music_|^bp_", "", cell_type)]
    plot_props[, group_binary := factor(group_binary, levels = c("Control", "Disease"))]

    disease_fill <- c(Control = masld_colors$control, Disease = masld_colors$masld)

    p_a <- ggplot(plot_props, aes(x = group_binary, y = proportion, fill = group_binary)) +
      geom_boxplot(outlier.size = 0.3, linewidth = 0.3, width = 0.6, alpha = 0.8) +
      geom_jitter(width = 0.15, size = 0.15, alpha = 0.2, color = "gray30") +
      facet_wrap(~ cell_type, scales = "free_y", nrow = 1) +
      scale_fill_manual(values = disease_fill, guide = "none") +
      scale_y_continuous(labels = percent_format()) +
      labs(x = NULL,
           y = "Estimated cell-type proportion",
           title = "MuSiC cell-type fractions") +
      theme_masld() +
      theme(
        strip.text  = element_text(size = 7, face = "plain"),
        plot.title  = element_text(size = 7, face = "plain"),
        axis.text.x = element_text(size = 6)
      )
  }
}

# ==========================================================================
# Panel b: Example genes — Intrinsic vs Composition-driven vs Unmasked
# ==========================================================================
if (!is.null(deconv)) {
  # Normalize columns
  if (!all(c("logFC_unadj", "logFC_adj", "category") %in% names(deconv))) {
    unadj_col <- intersect(c("logFC_unadj", "logFC_unadjusted", "logFC"), names(deconv))[1]
    adj_col   <- intersect(c("logFC_adj", "logFC_adjusted", "logFC_deconv"), names(deconv))[1]
    cat_col   <- intersect(c("category", "attribution_class", "class"), names(deconv))[1]
    if (!is.na(unadj_col)) setnames(deconv, unadj_col, "logFC_unadj")
    if (!is.na(adj_col))   setnames(deconv, adj_col,   "logFC_adj")
    if (!is.na(cat_col))   setnames(deconv, cat_col,   "category")
  }

  has_cols <- all(c("logFC_unadj", "logFC_adj", "category") %in% names(deconv))

  if (has_cols) {
    sig_deconv_f <- deconv[category != "Not_significant" & !grepl("^ENS(G|MUSG)", symbol)]

    if (nrow(sig_deconv_f) > 0) {
      pick_top <- function(cat) {
        sub <- sig_deconv_f[category == cat]
        if (nrow(sub) == 0) return(NULL)
        sub[which.max(abs(logFC_unadj - logFC_adj))]
      }

      examples <- rbindlist(lapply(
        c("Hepatocyte_intrinsic", "Composition_driven", "Unmasked"),
        pick_top
      ))

      if (nrow(examples) > 0) {
        bar_dt <- melt(examples[, .(symbol, category, logFC_unadj, logFC_adj)],
                       id.vars = c("symbol", "category"),
                       variable.name = "adjustment",
                       value.name = "logFC")
        bar_dt[, adjustment := fifelse(adjustment == "logFC_unadj",
                                       "Unadjusted", "Adjusted")]
        bar_dt[, adjustment := factor(adjustment, levels = c("Unadjusted", "Adjusted"))]

        bar_dt[, category_label := gsub("_", " ", category)]
        bar_dt[, facet_label := paste0(symbol, " (", category_label, ")")]
        facet_order <- bar_dt[!duplicated(facet_label)]$facet_label
        bar_dt[, facet_label := factor(facet_label, levels = facet_order)]

        adj_fill <- c(Unadjusted = "gray60", Adjusted = masld_colors$hep_intrinsic)

        p_b <- ggplot(bar_dt, aes(x = adjustment, y = logFC, fill = adjustment)) +
          geom_col(width = 0.6) +
          geom_hline(yintercept = 0, linewidth = 0.3, color = "gray30") +
          facet_wrap(~ facet_label, nrow = 1, scales = "free_y") +
          scale_fill_manual(values = adj_fill, name = NULL) +
          labs(x = NULL,
               y = expression("log"[2]*"FC"),
               title = "Effect of composition adjustment") +
          theme_masld() +
          theme(
            strip.text     = element_text(size = 6, face = "bold.italic"),
            plot.title     = element_text(size = 7, face = "plain"),
            axis.text.x    = element_text(angle = 30, hjust = 1, size = 6),
            legend.position = "none"
          )
      }
    }
  }
}

# ==========================================================================
# Panels c-e: Full BayesPrism cell-type deconvolution (16 types, 9 cohorts)
# ==========================================================================
meta_path2 <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration",
                        "metadata/unified_metadata.csv")

if (nrow(full_bp) > 0 && file.exists(meta_path2)) {
  meta2 <- fread(meta_path2, select = c("sample_id", "group_binary", "dataset"))
  full_bp <- merge(full_bp, meta2, by = c("sample_id", "dataset"), all.x = TRUE)
  full_bp <- full_bp[!is.na(group_binary)]

  # Identify cell-type columns (everything except metadata)
  meta_cols <- c("sample_id", "dataset", "group_binary")
  ct_cols <- setdiff(names(full_bp), meta_cols)

  cat("Full BayesPrism:", nrow(full_bp), "samples,", length(ct_cols), "cell types\n")

  # ========================================================================
  # Panel c: Non-parenchymal cell types (Control vs Disease boxplots)
  # Show the 6 most abundant non-hepatocyte types
  # ========================================================================
  ct_means <- sapply(ct_cols, function(x) mean(full_bp[[x]], na.rm = TRUE))
  # Exclude Hepatocytes (already in panel a), take top 6 by mean proportion
  ct_means_npc <- ct_means[!grepl("Hepatocyte", names(ct_means))]
  top_npc <- names(sort(ct_means_npc, decreasing = TRUE))[1:min(6, length(ct_means_npc))]

  npc_long <- melt(full_bp, id.vars = meta_cols, measure.vars = top_npc,
                   variable.name = "cell_type", value.name = "proportion")
  npc_long[, group_binary := factor(group_binary, levels = c("Control", "Disease"))]
  # Clean cell type names for display
  npc_long[, cell_type := factor(cell_type, levels = top_npc)]

  disease_fill <- c(Control = masld_colors$control, Disease = masld_colors$masld)

  p_c <- ggplot(npc_long, aes(x = group_binary, y = proportion, fill = group_binary)) +
    geom_boxplot(outlier.size = 0.3, linewidth = 0.3, width = 0.6, alpha = 0.8) +
    geom_jitter(width = 0.15, size = 0.1, alpha = 0.15, color = "gray30") +
    facet_wrap(~ cell_type, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = disease_fill, guide = "none") +
    scale_y_continuous(labels = percent_format()) +
    labs(x = NULL,
         y = "Estimated cell-type proportion",
         title = "Non-parenchymal cell-type fractions (BayesPrism)") +
    theme_masld() +
    theme(
      strip.text  = element_text(size = 6, face = "plain"),
      plot.title  = element_text(size = 7, face = "plain"),
      axis.text.x = element_text(size = 6)
    )

  # ========================================================================
  # Panel d: Stacked bar — mean composition by condition
  # ========================================================================
  comp_long <- melt(full_bp, id.vars = meta_cols, measure.vars = ct_cols,
                    variable.name = "cell_type", value.name = "proportion")
  comp_mean <- comp_long[, .(mean_prop = mean(proportion, na.rm = TRUE)),
                         by = .(group_binary, cell_type)]
  # Order cell types by overall abundance (Hepatocytes on bottom)
  ct_order <- names(sort(ct_means, decreasing = FALSE))
  comp_mean[, cell_type := factor(cell_type, levels = ct_order)]
  comp_mean[, group_binary := factor(group_binary, levels = c("Control", "Disease"))]

  # Use ct_palette from publication_theme.R, fill missing with grays
  fill_vals <- ct_palette[ct_cols]
  missing <- is.na(fill_vals)
  if (any(missing)) {
    gray_seq <- colorRampPalette(c("#78909C", "#CFD8DC"))(sum(missing))
    fill_vals[missing] <- gray_seq
  }
  names(fill_vals) <- ct_cols

  p_d <- ggplot(comp_mean, aes(x = group_binary, y = mean_prop, fill = cell_type)) +
    geom_col(width = 0.7, color = "white", linewidth = 0.2) +
    scale_fill_manual(values = fill_vals, name = "Cell type") +
    scale_y_continuous(labels = percent_format(), expand = expansion(mult = c(0, 0.02))) +
    labs(x = NULL, y = "Mean proportion",
         title = "Mean cell-type composition") +
    theme_masld() +
    theme(
      plot.title     = element_text(size = 7, face = "plain"),
      legend.text    = element_text(size = 5),
      legend.title   = element_text(size = 6),
      legend.key.size = unit(0.25, "cm"),
      axis.text.x    = element_text(size = 6)
    ) +
    guides(fill = guide_legend(ncol = 2, reverse = TRUE))

  # ========================================================================
  # Panel e: Per-dataset hepatocyte fraction (Control vs Disease)
  # ========================================================================
  hep_col_full <- intersect(c("Hepatocytes"), ct_cols)
  if (length(hep_col_full) > 0) {
    hep_ds <- full_bp[, .(sample_id, dataset, group_binary, proportion = get(hep_col_full))]
    hep_ds[, group_binary := factor(group_binary, levels = c("Control", "Disease"))]
    # Order datasets by median hepatocyte fraction
    ds_order <- hep_ds[, .(med = median(proportion, na.rm = TRUE)), by = dataset][order(-med)]$dataset
    hep_ds[, dataset := factor(dataset, levels = ds_order)]

    p_e <- ggplot(hep_ds, aes(x = dataset, y = proportion, fill = group_binary)) +
      geom_boxplot(outlier.size = 0.3, linewidth = 0.3, width = 0.7, alpha = 0.8,
                   position = position_dodge(width = 0.8)) +
      scale_fill_manual(values = disease_fill, name = "Group") +
      scale_y_continuous(labels = percent_format()) +
      labs(x = NULL, y = "Hepatocyte proportion",
           title = "Per-dataset hepatocyte fraction") +
      theme_masld() +
      theme(
        plot.title     = element_text(size = 7, face = "plain"),
        axis.text.x    = element_text(size = 5, angle = 45, hjust = 1),
        legend.position = "bottom",
        legend.text     = element_text(size = 6),
        legend.title    = element_text(size = 6)
      )
  }
}

# ==========================================================================
# Assemble: 3-row layout
#   Row 1: a (Hep+Mac) | b (attribution examples)
#   Row 2: c (non-parenchymal cell types, full width)
#   Row 3: d (stacked bars) | e (per-dataset hepatocyte)
# ==========================================================================
row1 <- p_a + p_b + plot_layout(widths = c(1, 1.5))
row2 <- p_c
row3 <- p_d + p_e + plot_layout(widths = c(1, 2))

figS7 <- row1 / row2 / row3 +
  plot_layout(heights = c(1, 1, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "plain"))

save_fig_tall(figS7, OUT, width = fig_full_width, height = 10)
message("Supp Fig 7 (Extended Cell-Type) saved to ", OUT)
