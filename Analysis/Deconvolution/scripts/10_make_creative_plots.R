args <- commandArgs(trailingOnly = TRUE)
results_dir <- if (length(args) >= 1) args[[1]] else "RNA-seq/deconvolution/results"
bulk_dir <- if (length(args) >= 2) args[[2]] else "RNA-seq/deconvolution/bulk"
method <- if (length(args) >= 3) args[[3]] else "music" # music or bayesprism

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

# Configure PDF device
pdf.options(useDingbats = FALSE)

palette_base <- c(
  "#D81B60", # magenta
  "#EC407A", # pink
  "#FF5DA2", # hot pink
  "#FB8C00", # orange
  "#FFA726", # light orange
  "#FF7043", # coral
  "#FDD835", # gold
  "#AB47BC", # purple-pink
  "#8E24AA", # deep magenta
  "#F06292"  # soft pink
)

make_palette <- function(n) {
  if (n <= length(palette_base)) {
    return(palette_base[seq_len(n)])
  }
  return(colorRampPalette(palette_base)(n))
}

theme_creative <- function() {
  theme_minimal(base_family = "Helvetica", base_size = 7) +
    theme(
      plot.title = element_text(face = "bold", size = 8, color = "black"),
      plot.subtitle = element_text(size = 7, color = "black"),
      axis.title = element_text(size = 8, color = "black"),
      axis.text = element_text(size = 6, color = "black"),
      panel.grid.major = element_line(color = "#F3E5F5"),
      panel.grid.minor = element_blank(),
      legend.title = element_text(size = 7, color = "black"),
      legend.text = element_text(size = 6, color = "black")
    )
}

read_props <- function(prop_path, meta_path) {
  props <- read.delim(prop_path, check.names = FALSE, row.names = 1)
  props <- as.matrix(props)
  storage.mode(props) <- "numeric"

  meta <- read.delim(meta_path, check.names = FALSE, row.names = 1)
  meta$sample_id <- rownames(meta)

  row_hit <- intersect(rownames(props), meta$sample_id)
  col_hit <- intersect(colnames(props), meta$sample_id)

  if (length(row_hit) > 0) {
    props_samples <- props[row_hit, , drop = FALSE]
  } else if (length(col_hit) > 0) {
    props_samples <- t(props)
    props_samples <- props_samples[intersect(rownames(props_samples), meta$sample_id), , drop = FALSE]
  } else {
    stop("No overlapping samples between props and metadata")
  }

  meta <- meta[match(rownames(props_samples), meta$sample_id), , drop = FALSE]
  list(props = props_samples, meta = meta)
}

pick_group <- function(meta) {
  candidates <- c(
    "diet",
    "disease",
    "group_in_paper",
    "Stage",
    "stage",
    "fibrosis_stage",
    "sex",
    "week"
  )
  candidates <- candidates[candidates %in% colnames(meta)]
  if (length(candidates) == 0) return(NULL)
  for (col in candidates) {
    n_unique <- length(unique(meta[[col]]))
    if (n_unique >= 2 && n_unique <= 8) return(col)
  }
  return(candidates[[1]])
}

save_plot <- function(plot, path_base, width = 8, height = 5) {
  ggsave(paste0(path_base, ".pdf"), plot = plot, width = width, height = height)
}

# Dataset Labels
dataset_labels <- c(
  "GSE130970" = "Human (Hoang)",
  "GSE135251" = "Human (Govaere)",
  "inhouse_MCD" = "Cas13 mouse MCD",
  "GSE156918" = "Mouse MCD (Paquette)",
  "GSE205974" = "Mouse MCD (Yue)"
)

# Helper to get display name
get_display_name <- function(dataset) {
  if (dataset %in% names(dataset_labels)) {
    return(dataset_labels[[dataset]])
  }
  return(dataset)
}

make_dataset_plots <- function(dataset) {
  display_name <- get_display_name(dataset)
  if (method == "music") {
    prop_path <- file.path(results_dir, dataset, paste0(dataset, "_music_prop_weighted.tsv"))
  } else {
    prop_path <- file.path(results_dir, dataset, paste0(dataset, "_bayesprism_proportions.tsv"))
  }
  
  meta_path <- file.path(bulk_dir, dataset, paste0(dataset, "_metadata.tsv"))
  out_dir <- file.path(results_dir, dataset, "qc", paste0("creative_", method))

  if (!file.exists(prop_path) || !file.exists(meta_path)) {
    message("Skipping ", dataset, ": missing inputs")
    return(NULL)
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  data <- read_props(prop_path, meta_path)
  props <- data$props
  meta <- data$meta

  celltypes <- colnames(props)
  pal <- make_palette(length(celltypes))
  names(pal) <- celltypes

  # Long format
  df_long <- as.data.frame(props)
  df_long$sample_id <- rownames(df_long)
  df_long <- pivot_longer(df_long, cols = -sample_id, names_to = "celltype", values_to = "fraction")

  # Stacked bar sorted by dominant celltype
  dominant <- df_long %>%
    group_by(sample_id) %>%
    summarize(top_fraction = max(fraction, na.rm = TRUE), .groups = "drop")
  order_samples <- dominant %>% arrange(desc(top_fraction)) %>% pull(sample_id)
  df_long$sample_id <- factor(df_long$sample_id, levels = order_samples)

  stacked <- ggplot(df_long, aes(x = sample_id, y = fraction, fill = celltype)) +
    geom_col(width = 0.9) +
    coord_cartesian(ylim = c(0, 1)) +
    scale_fill_manual(values = pal) +
    labs(title = paste0(display_name, ": Cell type composition"), x = "Sample", y = "Fraction") +
    theme_creative() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6))
  save_plot(stacked, file.path(out_dir, paste0(dataset, "_", method, "_creative_stacked_sorted")), width = max(6, length(order_samples) * 0.15), height = 4)

  # Heatmap with clustering
  sample_dist <- dist(props)
  sample_order <- hclust(sample_dist)$order
  cell_dist <- dist(t(props))
  cell_order <- hclust(cell_dist)$order

  props_ord <- props[sample_order, cell_order, drop = FALSE]
  df_heat <- as.data.frame(props_ord)
  df_heat$sample_id <- rownames(df_heat)
  df_heat <- pivot_longer(df_heat, cols = -sample_id, names_to = "celltype", values_to = "fraction")
  df_heat$celltype <- factor(df_heat$celltype, levels = colnames(props_ord))
  df_heat$sample_id <- factor(df_heat$sample_id, levels = rownames(props_ord))

  heatmap <- ggplot(df_heat, aes(x = sample_id, y = celltype, fill = fraction)) +
    geom_tile() +
    scale_fill_gradient(low = "#FCE4EC", high = "#FF6F00") +
    labs(title = paste0(display_name, ": Clustered heatmap"), x = "Sample", y = "Cell type") +
    theme_creative() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  save_plot(heatmap, file.path(out_dir, paste0(dataset, "_", method, "_creative_heatmap_clustered")), width = max(6, nrow(props_ord) * 0.1), height = max(3, ncol(props_ord) * 0.25))

  # PCA of proportions (drop zero-variance columns)
  var_cols <- apply(props, 2, var, na.rm = TRUE)
  props_pca <- props[, var_cols > 0, drop = FALSE]
  pca_df <- NULL
  if (ncol(props_pca) >= 2) {
    pca <- prcomp(props_pca, scale. = TRUE)
    pca_df <- as.data.frame(pca$x[, 1:2, drop = FALSE])
    pca_df$sample_id <- rownames(props)
  }
  group_col <- pick_group(meta)
  if (!is.null(pca_df)) {
    if (!is.null(group_col)) {
      pca_df$group <- as.factor(meta[[group_col]])
      pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
        geom_point(size = 2.2, alpha = 0.85) +
        scale_color_manual(values = make_palette(length(levels(pca_df$group)))) +
        labs(title = paste0(display_name, ": PCA of proportions"), subtitle = paste0("Grouped by ", group_col)) +
        theme_creative()
    } else {
      pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
        geom_point(size = 2.2, color = "#D81B60", alpha = 0.85) +
        labs(title = paste0(display_name, ": PCA of proportions")) +
        theme_creative()
    }
    save_plot(pca_plot, file.path(out_dir, paste0(dataset, "_", method, "_creative_pca")), width = 5, height = 4)
  } else {
    message("Skipping PCA for ", dataset, ": not enough variable cell types")
  }

  # Boxplot by group (if available)
  if (!is.null(group_col)) {
    df_long <- df_long %>% left_join(meta[, c("sample_id", group_col)], by = "sample_id")
    df_long$group <- as.factor(df_long[[group_col]])

    boxplot <- ggplot(df_long, aes(x = group, y = fraction, fill = celltype)) +
      geom_boxplot(outlier.size = 0.6, alpha = 0.85) +
      scale_fill_manual(values = pal) +
      facet_wrap(~celltype, scales = "free_y") +
      labs(title = paste0(display_name, ": Cell type distribution by ", group_col), x = group_col, y = "Fraction") +
      theme_creative() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
    save_plot(boxplot, file.path(out_dir, paste0(dataset, "_", method, "_creative_boxplot_by_", group_col)), width = 8, height = 5)
  }

  # Radial mean composition
  mean_df <- df_long %>%
    group_by(celltype) %>%
    summarize(mean_fraction = mean(fraction, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(mean_fraction))
  mean_df$celltype <- factor(mean_df$celltype, levels = mean_df$celltype)

  radial <- ggplot(mean_df, aes(x = celltype, y = mean_fraction, fill = celltype)) +
    geom_col(width = 1) +
    coord_polar() +
    scale_fill_manual(values = pal) +
    labs(title = paste0(display_name, ": Mean composition (radial)"), x = NULL, y = NULL) +
    theme_creative() +
    theme(axis.text.x = element_text(size = 8), axis.text.y = element_blank(), axis.ticks = element_blank(), legend.position = "none")
  save_plot(radial, file.path(out_dir, paste0(dataset, "_", method, "_creative_radial")), width = 5, height = 5)

  # Celltype correlation heatmap
  cor_mat <- cor(props, use = "pairwise.complete.obs")
  cor_df <- as.data.frame(cor_mat)
  cor_df$celltype1 <- rownames(cor_df)
  cor_df <- pivot_longer(cor_df, cols = -celltype1, names_to = "celltype2", values_to = "correlation")

  corr <- ggplot(cor_df, aes(x = celltype1, y = celltype2, fill = correlation)) +
    geom_tile() +
    scale_fill_gradient2(low = "#F06292", mid = "white", high = "#FF8F00", midpoint = 0) +
    labs(title = paste0(display_name, ": Cell type correlation"), x = NULL, y = NULL) +
    theme_creative() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot(corr, file.path(out_dir, paste0(dataset, "_", method, "_creative_correlation")), width = 6, height = 5)

  invisible(TRUE)
}

# Dataset list from results directory
all_datasets <- list.dirs(results_dir, full.names = FALSE, recursive = FALSE)
all_datasets <- all_datasets[all_datasets != "summary_report.md"]

# Filter for known datasets
all_datasets <- all_datasets[all_datasets %in% names(dataset_labels)]

for (ds in all_datasets) {
  message("Plotting ", ds)
  make_dataset_plots(ds)
}

# Cross-dataset mean heatmap
means <- list()
for (ds in all_datasets) {
  if (method == "music") {
    prop_path <- file.path(results_dir, ds, paste0(ds, "_music_prop_weighted.tsv"))
  } else {
    prop_path <- file.path(results_dir, ds, paste0(ds, "_bayesprism_proportions.tsv"))
  }
  meta_path <- file.path(bulk_dir, ds, paste0(ds, "_metadata.tsv"))
  if (!file.exists(prop_path) || !file.exists(meta_path)) next
  data <- read_props(prop_path, meta_path)
  props <- data$props
  mean_vec <- colMeans(props, na.rm = TRUE)
  means[[ds]] <- mean_vec
}

if (length(means) > 0) {
  all_celltypes <- unique(unlist(lapply(means, names)))
  mat <- do.call(rbind, lapply(means, function(v) {
    v[all_celltypes][is.na(v[all_celltypes])] <- 0
    v[all_celltypes]
  }))
  rownames(mat) <- names(means)
  mat[is.na(mat)] <- 0

  df <- as.data.frame(mat)
  df$dataset_id <- rownames(mat)
  # Map to display names
  df$dataset <- sapply(df$dataset_id, get_display_name)
  df <- pivot_longer(df, cols = -c(dataset, dataset_id), names_to = "celltype", values_to = "mean_fraction")

  out_dir <- file.path(results_dir, "summary_plots")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  cross <- ggplot(df, aes(x = dataset, y = celltype, fill = mean_fraction)) +
    geom_tile() +
    scale_fill_gradient(low = "#FCE4EC", high = "#C2185B") +
    labs(title = "Mean cell type fractions across datasets", x = NULL, y = NULL) +
    theme_creative() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot(cross, file.path(out_dir, paste0(method, "_creative_cross_dataset_heatmap")), width = 6, height = 4)

  bubble <- ggplot(df, aes(x = dataset, y = celltype, size = mean_fraction, color = mean_fraction)) +
    geom_point(alpha = 0.85) +
    scale_color_gradient(low = "#FCE4EC", high = "#C2185B") +
    scale_size(range = c(1, 8)) +
    labs(title = "Cross-dataset mean fractions (bubble)", x = NULL, y = NULL) +
    theme_creative() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")
  save_plot(bubble, file.path(out_dir, paste0(method, "_creative_cross_dataset_bubble")), width = 6.5, height = 4.5)

  top_celltypes <- df %>%
    group_by(celltype) %>%
    summarize(global_mean = mean(mean_fraction, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(global_mean)) %>%
    slice_head(n = 6) %>%
    pull(celltype)

  df_top <- df %>% filter(celltype %in% top_celltypes)
  df_other <- df %>%
    group_by(dataset) %>%
    summarize(other = max(0, 1 - sum(mean_fraction[celltype %in% top_celltypes])), .groups = "drop") %>%
    mutate(celltype = "Other")
  df_stacked <- bind_rows(
    df_top,
    df_other %>% rename(mean_fraction = other)
  )

  stacked_pal <- make_palette(length(top_celltypes))
  names(stacked_pal) <- top_celltypes
  stacked_pal <- c(stacked_pal, Other = "#F8BBD0")

  stacked <- ggplot(df_stacked, aes(x = dataset, y = mean_fraction, fill = celltype)) +
    geom_col(position = "fill", width = 0.8) +
    scale_fill_manual(values = stacked_pal) +
    labs(title = "Top cell types across datasets (stacked)", x = NULL, y = "Fraction") +
    theme_creative() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot(stacked, file.path(out_dir, paste0(method, "_creative_cross_dataset_stacked_top")), width = 6.5, height = 4)

  lollipop_df <- df %>%
    group_by(dataset) %>%
    arrange(desc(mean_fraction)) %>%
    slice_head(n = 5) %>%
    ungroup()

  lollipop <- ggplot(lollipop_df, aes(x = mean_fraction, y = reorder(celltype, mean_fraction))) +
    geom_segment(aes(x = 0, xend = mean_fraction, yend = celltype), color = "#F48FB1", linewidth = 1) +
    geom_point(color = "#C2185B", size = 2.8) +
    facet_wrap(~dataset, scales = "free_y") +
    labs(title = "Top 5 cell types per dataset", x = "Mean fraction", y = NULL) +
    theme_creative()
  save_plot(lollipop, file.path(out_dir, paste0(method, "_creative_cross_dataset_lollipop")), width = 7.5, height = 5)
}
