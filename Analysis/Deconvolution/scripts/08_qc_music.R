args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript 08_qc_music.R <music_prop_tsv> <metadata_tsv> <output_dir> <dataset_name> [min_mean]",
       call. = FALSE)
}

music_prop_tsv <- args[[1]]
metadata_tsv <- args[[2]]
output_dir <- args[[3]]
dataset_name <- args[[4]]
min_mean <- if (length(args) >= 5) as.numeric(args[[5]]) else 0.01

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

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

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

props <- read.delim(music_prop_tsv, check.names = FALSE, row.names = 1)
props <- as.matrix(props)
if (nrow(props) == 0 || ncol(props) == 0) {
  stop("MuSiC proportions table is empty")
}

storage.mode(props) <- "numeric"

cell_means <- rowMeans(props, na.rm = TRUE)
keep <- cell_means >= min_mean
if (any(!keep)) {
  other <- colSums(props[!keep, , drop = FALSE], na.rm = TRUE)
  props <- rbind(props[keep, , drop = FALSE], Other = other)
  cell_means <- rowMeans(props, na.rm = TRUE)
}

props <- props[order(cell_means, decreasing = TRUE), , drop = FALSE]

meta <- read.delim(metadata_tsv, check.names = FALSE, row.names = 1)
meta$sample_id <- rownames(meta)

samples <- intersect(colnames(props), meta$sample_id)
if (length(samples) == 0) {
  row_samples <- intersect(rownames(props), meta$sample_id)
  if (length(row_samples) > 0) {
    props <- t(props)
    samples <- intersect(colnames(props), meta$sample_id)
  }
}
if (length(samples) == 0) {
  stop("No overlapping samples between MuSiC proportions and metadata")
}
props <- props[, samples, drop = FALSE]
meta <- meta[match(samples, meta$sample_id), , drop = FALSE]

props_df <- as.data.frame(t(props))
props_df$sample_id <- rownames(props_df)

long <- tidyr::pivot_longer(
  props_df,
  cols = -sample_id,
  names_to = "celltype",
  values_to = "fraction"
)

cell_summary <- long %>%
  group_by(celltype) %>%
  summarize(
    mean_fraction = mean(fraction, na.rm = TRUE),
    median_fraction = median(fraction, na.rm = TRUE),
    sd_fraction = sd(fraction, na.rm = TRUE),
    pct_nonzero = mean(fraction > 0, na.rm = TRUE) * 100,
    pct_gt_01 = mean(fraction >= 0.01, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  arrange(desc(mean_fraction))

sample_summary <- long %>%
  group_by(sample_id) %>%
  summarize(
    total_fraction = sum(fraction, na.rm = TRUE),
    top_celltype = celltype[which.max(fraction)],
    top_fraction = max(fraction, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(meta, by = "sample_id")

plot_width <- min(20, max(5, length(samples) * 0.15))
heatmap_height <- min(20, max(3, nrow(props) * 0.25))

write.table(
  cell_summary,
  file.path(output_dir, paste0(dataset_name, "_celltype_summary.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  sample_summary,
  file.path(output_dir, paste0(dataset_name, "_sample_summary.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

stacked_plot <- ggplot(long, aes(x = sample_id, y = fraction, fill = celltype)) +
  geom_col(width = 0.9) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = paste0(get_display_name(dataset_name), " MuSiC cell type proportions"),
    x = "Sample",
    y = "Estimated fraction"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
    panel.grid.major.x = element_blank()
  )

ggsave(
  filename = file.path(output_dir, paste0(dataset_name, "_stacked_bar.png")),
  plot = stacked_plot,
  width = plot_width,
  height = 5,
  dpi = 300
)

ggsave(
  filename = file.path(output_dir, paste0(dataset_name, "_stacked_bar.pdf")),
  plot = stacked_plot,
  width = plot_width,
  height = 5
)

heatmap_plot <- ggplot(long, aes(x = sample_id, y = celltype, fill = fraction)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(
    title = paste0(get_display_name(dataset_name), " MuSiC proportion heatmap"),
    x = "Sample",
    y = "Cell type"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

ggsave(
  filename = file.path(output_dir, paste0(dataset_name, "_heatmap.png")),
  plot = heatmap_plot,
  width = plot_width,
  height = heatmap_height,
  dpi = 300
)

ggsave(
  filename = file.path(output_dir, paste0(dataset_name, "_heatmap.pdf")),
  plot = heatmap_plot,
  width = plot_width,
  height = heatmap_height
)

preferred_groups <- c(
  "diet",
  "disease",
  "group_in_paper",
  "Stage",
  "stage",
  "fibrosis_stage",
  "sex",
  "week"
)

available_groups <- preferred_groups[preferred_groups %in% colnames(meta)]
if (length(available_groups) > 0) {
  group_col <- available_groups[[1]]
  grouped <- long %>%
    left_join(meta[, c("sample_id", group_col), drop = FALSE], by = "sample_id") %>%
    group_by(.data[[group_col]], celltype) %>%
    summarize(mean_fraction = mean(fraction, na.rm = TRUE), .groups = "drop")

  grouped[[group_col]] <- as.factor(grouped[[group_col]])

  group_plot <- ggplot(
    grouped,
    aes(x = .data[[group_col]], y = mean_fraction, fill = celltype)
  ) +
    geom_col(width = 0.7) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
      title = paste0(get_display_name(dataset_name), " mean proportions by ", group_col),
      x = group_col,
      y = "Mean fraction"
    ) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(
    filename = file.path(output_dir, paste0(dataset_name, "_group_means_by_", group_col, ".png")),
    plot = group_plot,
    width = 7,
    height = 4,
    dpi = 300
  )
  ggsave(
    filename = file.path(output_dir, paste0(dataset_name, "_group_means_by_", group_col, ".pdf")),
    plot = group_plot,
    width = 7,
    height = 4
  )
}

cat("QC outputs written to", output_dir, "\n")
