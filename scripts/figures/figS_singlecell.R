suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(scales)
  library(dplyr)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUT <- file.path(FIGS03_DIR, "figS_singlecell.pdf")

scvi_dir <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2")
h5ad_file <- file.path(scvi_dir, "integrated_atlas.h5ad")
scvi_ready <- file.exists(h5ad_file)

pb_de_dir <- file.path(scvi_dir, "pseudobulk_de")

# ============================================================================
# Cell-type color palette (17 types, visually distinct, colorblind-friendly)
# ============================================================================
ct_palette <- c(
  "Hepatocytes"          = "#0D47A1",
  "Cholangiocytes"       = "#1565C0",
  "Endothelial cells"    = "#2E7D32",
  "Fibroblasts"          = "#F57F17",
  "Macrophages"          = "#C2185B",
  "Mono+mono derived cells" = "#E91E63",
  "T cells"              = "#7B1FA2",
  "B cells"              = "#9C27B0",
  "Resident NK"          = "#00695C",
  "Circulating NK/NKT"   = "#00897B",
  "Plasma cells"         = "#5D4037",
  "Neutrophils"          = "#FF6F00",
  "cDC1s"                = "#AD1457",
  "cDC2s"                = "#D81B60",
  "pDCs"                 = "#6A1B9A",
  "Mig.cDCs"             = "#AB47BC",
  "Basophils"            = "#78909C"
)

if (scvi_ready) {
  message("scVI results found: ", h5ad_file)

  # ==========================================================================
  # Load h5ad metadata via rhdf5
  # ==========================================================================
  library(rhdf5)

  umap_raw <- h5read(h5ad_file, "obsm/X_umap")
  umap_coords <- data.frame(UMAP_1 = umap_raw[1, ], UMAP_2 = umap_raw[2, ])

  items <- h5ls(h5ad_file)

  read_obs_col <- function(key) {
    group_path <- paste0("/obs/", key)
    obs_groups <- subset(items, group == "/obs" & otype == "H5I_GROUP")$name
    if (key %in% obs_groups) {
      codes <- h5read(h5ad_file, paste0(group_path, "/codes"))
      cats  <- h5read(h5ad_file, paste0(group_path, "/categories"))
      return(cats[codes + 1L])
    }
    vals <- h5read(h5ad_file, group_path)
    if (is.list(vals) && "codes" %in% names(vals))
      return(vals$categories[vals$codes + 1L])
    as.vector(vals)
  }

  obs_meta <- data.frame(
    cell_type = tryCatch(read_obs_col("cell_type"), error = function(e) NA),
    dataset   = tryCatch(read_obs_col("dataset"),   error = function(e) NA),
    sample    = tryCatch(read_obs_col("sample"),     error = function(e) NA),
    condition = tryCatch(read_obs_col("condition"),  error = function(e) NA),
    preparation_method = tryCatch(read_obs_col("preparation_method"), error = function(e) NA),
    condition_harmonized = tryCatch(read_obs_col("condition_harmonized"), error = function(e) NA)
  )

  # Collapse condition into simplified disease status
  # If condition_harmonized is available (from add_sample_metadata.py), use it directly
  if (!all(is.na(obs_meta$condition_harmonized))) {
    obs_meta$disease_status <- dplyr::case_when(
      obs_meta$condition_harmonized == "Healthy" ~ "Control",
      obs_meta$condition_harmonized == "MASLD"   ~ "MASLD",
      TRUE ~ "Unknown"
    )
    message("  Using condition_harmonized from h5ad metadata")
  } else {
    # Fallback: manual condition mapping (pre-add_sample_metadata.py)
    obs_meta$disease_status <- dplyr::case_when(
      obs_meta$condition == "Healthy"  ~ "Control",
      obs_meta$condition == "Mixed"    ~ "Control",  # Liver_Atlas healthy donors
      obs_meta$condition %in% c("MASLD", "NAFLD", "NASH") ~ "MASLD",
      TRUE ~ "Unknown"
    )
    message("  WARNING: condition_harmonized not found, using fallback mapping")
  }

  plot_df <- cbind(umap_coords, obs_meta)
  n_total_cells <- nrow(plot_df)
  message("  Loaded ", n_total_cells, " human cells")

  # Subsample for UMAP rendering
  if (nrow(plot_df) > 100000) {
    set.seed(42)
    plot_df_sub <- plot_df[sample(nrow(plot_df), 100000), ]
  } else {
    plot_df_sub <- plot_df
  }

  # ========================================================================
  # Panel A: UMAP by cell type
  # ========================================================================
  message("Panel A: UMAP...")

  # Shuffle to avoid overplotting bias
  plot_df_sub <- plot_df_sub[sample(nrow(plot_df_sub)), ]

  p_a <- ggplot(plot_df_sub, aes(x = UMAP_1, y = UMAP_2, color = cell_type)) +
    rasterize_layer(geom_point(size = 0.02, alpha = 0.4, stroke = 0, shape = 16)) +
    scale_color_manual(values = ct_palette, na.value = "grey80") +
    theme_masld() +
    theme(legend.position = "right",
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.key.size = unit(0.2, "cm"),
          legend.spacing.y = unit(0.05, "cm")) +
    guides(color = guide_legend(override.aes = list(size = 1.5, alpha = 1),
                                ncol = 1)) +
    labs(x = "UMAP 1", y = "UMAP 2", color = NULL,
         title = paste0("Cell types (", format(n_total_cells, big.mark = ","), " cells)"))

  # ========================================================================
  # Panel C: Cell-type composition shifts (Grouped Boxplots: MASLD vs Control)
  # ========================================================================
  message("Panel C: Proportions...")
  props_file <- file.path(scvi_dir, "cell_type_proportions.csv")

  if (file.exists(props_file)) {
    props_raw <- fread(props_file)

    # Check if updated proportions file has preparation_method column
    has_prep <- "preparation_method" %in% names(props_raw)
    has_cond_harm <- "condition_harmonized" %in% names(props_raw)

    # Identify metadata vs cell-type proportion columns
    meta_cols <- intersect(names(props_raw),
                           c("sample", "dataset", "condition", "preparation_method",
                             "condition_harmonized"))
    ct_cols <- setdiff(names(props_raw), meta_cols)

    if (has_prep) {
      # Filter to unsorted/nuclei/CD45- samples only (exclude sorted for proportions)
      de_preps <- c("nuclei", "cd45_negative", "unsorted")
      props_raw <- props_raw[preparation_method %in% de_preps]
      message("  Filtered to ", nrow(props_raw), " unsorted/nuclei/CD45- samples for proportions")
    }

    if (has_cond_harm) {
      # Use condition_harmonized for disease status
      props_raw[, disease_status := fifelse(condition_harmonized == "Healthy", "Control",
                                   fifelse(condition_harmonized == "MASLD", "MASLD", "Unknown"))]
    } else {
      # Fallback: merge from obs_meta
      sample_status <- unique(obs_meta[, c("sample", "disease_status")])
      props_raw <- merge(props_raw, sample_status, by = "sample", all.x = TRUE)
    }

    # Filter to only Control and MASLD samples
    props_filt <- props_raw[disease_status %in% c("Control", "MASLD")]

    # Calc row sums for remaining valid cell types
    props_filt[, row_sum := rowSums(.SD, na.rm=TRUE), .SDcols = ct_cols]
    props <- props_filt[row_sum > 0]
    props[, row_sum := NULL]

    # Melt to long format
    id_vars <- intersect(names(props), c("sample", "dataset", "disease_status"))
    props_long <- melt(props, id.vars = id_vars,
                       measure.vars = ct_cols,
                       variable.name = "cell_type", value.name = "proportion")
    
    # Reorder cell types by descending median proportion (across all to get a stable order)
    ct_medians <- props_long[, .(med_prop = median(proportion)), by = cell_type][order(-med_prop)]
    
    # Filter to focus on the top 6 most abundant cell types for clarity
    top_6_cts <- ct_medians$cell_type[1:6]
    props_long <- props_long[cell_type %in% top_6_cts]
    props_long[, cell_type := factor(cell_type, levels = top_6_cts)]

    # Explicitly set MASLD/Control palette and order
    props_long[, disease_status := factor(disease_status, levels = c("Control", "MASLD"))]
    ds_palette <- c("Control" = "#2E7D32", "MASLD" = "#C62828")

    p_c <- ggplot(props_long, aes(x = disease_status, y = proportion, fill = disease_status)) +
      geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.5, alpha = 0.8, linewidth = 0.3) +
      facet_wrap(~ cell_type, scales = "free_y", ncol = 3) +
      scale_fill_manual(values = ds_palette) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
      theme_masld() +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
            axis.text.y = element_text(size = 6),
            strip.text = element_text(size = 7, face = "plain", margin = margin(b=4, t=4))) +
      labs(x = NULL, y = "Proportion of cells", title = "Cell-type abundance")
      
  } else {
    p_c <- placeholder("Cell proportions\n(file not found)")
  }

  # ========================================================================
  # Panel B: UMAP by disease status
  # ========================================================================
  message("Panel B: UMAP by disease status...")

  disease_palette <- c(
    "Control"       = "#2E7D32",
    "MASLD"         = "#C62828"
  )

  # Filter out Mixed/Unknown for cleaner visualization
  plot_df_disease <- plot_df_sub[plot_df_sub$disease_status %in% c("Control", "MASLD"), ]
  plot_df_disease$disease_status <- factor(
    plot_df_disease$disease_status,
    levels = c("Control", "MASLD")
  )
  plot_df_disease <- plot_df_disease[order(plot_df_disease$disease_status), ]

  n_ctrl <- sum(obs_meta$disease_status == "Control")
  n_masld <- sum(obs_meta$disease_status == "MASLD")

  p_b <- ggplot(plot_df_disease, aes(x = UMAP_1, y = UMAP_2, color = disease_status)) +
    rasterize_layer(geom_point(size = 0.02, alpha = 0.4, stroke = 0, shape = 16)) +
    scale_color_manual(values = disease_palette, na.value = "grey90") +
    theme_masld() +
    theme(legend.position = "right",
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.key.size = unit(0.3, "cm")) +
    guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
    labs(x = "UMAP 1", y = "UMAP 2", color = NULL,
         title = paste0("Disease status (",
                        format(n_ctrl, big.mark = ","), " ctrl, ",
                        format(n_masld, big.mark = ","), " MASLD)"))

  # ========================================================================
  # Panel D: DEG counts per cell type (horizontal bar)
  # ========================================================================
  message("Panel D: DEG summary...")
  de_files <- list.files(pb_de_dir, pattern = "_de\\.csv$", full.names = TRUE)

  if (length(de_files) > 0) {
    deg_summary <- rbindlist(lapply(de_files, function(f) {
      ct <- gsub("_de\\.csv$", "", basename(f))
      dt <- fread(f)
      n_up <- sum(dt$padj < 0.05 & dt$logFC > 0, na.rm = TRUE)
      n_down <- sum(dt$padj < 0.05 & dt$logFC < 0, na.rm = TRUE)
      data.table(cell_type = ct, direction = c("Up", "Down"), n = c(n_up, n_down))
    }))

    # Clean up duplicate cell types (e.g., "T cells" vs "T_cells")
    deg_summary[, ct_clean := gsub("_", " ", gsub("\\+", "+", cell_type))]
    # Deduplicate: keep the version with more DEGs
    deg_summary[, total := sum(n), by = .(ct_clean, direction)]
    deg_summary <- deg_summary[, .SD[which.max(n)], by = .(ct_clean, direction)]

    # Order by total DEGs
    ct_totals <- deg_summary[, .(total = sum(n)), by = ct_clean][order(total)]
    deg_summary[, ct_clean := factor(ct_clean, levels = ct_totals$ct_clean)]
    deg_summary[direction == "Down", n := -n]

    p_d <- ggplot(deg_summary, aes(x = n, y = ct_clean, fill = direction)) +
      geom_bar(stat = "identity", width = 0.7) +
      scale_fill_manual(values = c(Up = masld_colors$up, Down = masld_colors$down)) +
      geom_vline(xintercept = 0, linewidth = 0.3) +
      scale_x_continuous(labels = function(x) format(abs(x), big.mark = ",")) +
      theme_masld() +
      theme(legend.position = "bottom",
            legend.key.size = unit(0.2, "cm")) +
      labs(x = "Number of DEGs (padj < 0.05)", y = NULL, fill = NULL,
           title = "Pseudobulk DE per cell type")
  } else {
    p_d <- placeholder("DEG counts per cell type\n(run pseudobulk_de.R first)")
  }

} else {
  message("scVI results not yet available at: ", h5ad_file)
  p_a <- placeholder("Panel a: UMAP by cell type\n(awaiting integration)")
  p_b <- placeholder("Panel b: UMAP by disease status\n(awaiting integration)")
  p_c <- placeholder("Panel c: Cell-type proportions\n(awaiting integration)")
  p_d <- placeholder("Panel d: Pseudobulk DE per cell type\n(awaiting integration)")
}

# ============================================================================
# Assemble: 4-panel layout
# Top row: A and B get equal width
# Bottom row: C gets 2x the width of D, and the row gets 2x the height
# ============================================================================
message("Assembling Figure...")

design <- "
AAABBB
CCCCDD
CCCCDD
"

fig_out <- p_a + p_b + p_c + p_d +
  plot_layout(design = design) +
  plot_annotation(
    tag_levels = "a"
  ) &
  theme(plot.tag = element_text(size = 8, face = "plain"))

save_fig_tall(fig_out, OUT, height = 8)
message("Figure saved to ", OUT)
