#!/usr/bin/env Rscript
#' Publication-quality TPM Distribution Plots (Unfiltered Version)
#'
#' Computes TPM directly from RAW count matrices (unfiltered) to ensure
#' fair comparison across all datasets. Uses cat_palette from publication_color_themes.R.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ggridges)
  library(patchwork)
  library(scales)
  library(data.table)
})

# =============================================================================
# Configuration
# =============================================================================

ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
source(file.path(ROOT, "publication_color_themes.R"))

OUTPUT_DIR <- file.path(ROOT, "RNA-seq", "TPM_analysis", "unfiltered")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Raw Data Paths (Unfiltered Count Matrices)
# =============================================================================

RAW_DATA_CONFIG <- list(
  # Mouse MCD datasets (featureCounts format)
  list(
    label = "Cas13 mouse MCD",
    species = "Mouse",
    type = "featurecounts",
    counts_path = file.path(ROOT, "RNA-seq/in-house_MCD_RNAseq/counts/featurecounts/gene_counts.txt"),
    metadata_path = file.path(ROOT, "RNA-seq/in-house_MCD_RNAseq/metadata/samples.tsv"),
    filter_col = "diet",
    filter_pattern = "mcd"
  ),
  list(
    label = "Mouse MCD (Paquette)",
    species = "Mouse",
    type = "featurecounts",
    counts_path = file.path(ROOT, "RNA-seq/other_MCD_RNAseq/GSE156918/counts/featurecounts/gene_counts.txt"),
    metadata_path = file.path(ROOT, "RNA-seq/other_MCD_RNAseq/GSE156918/metadata/samples.tsv"),
    filter_col = "diet",
    filter_pattern = "mcd"
  ),
  list(
    label = "Mouse MCD (Yue)",
    species = "Mouse",
    type = "featurecounts",
    counts_path = file.path(ROOT, "RNA-seq/other_MCD_RNAseq/GSE205974/counts/featurecounts/gene_counts.txt"),
    metadata_path = file.path(ROOT, "RNA-seq/other_MCD_RNAseq/GSE205974/metadata/samples.tsv"),
    filter_col = "diet",
    filter_pattern = "mcd"
  ),
  # Human MASLD datasets (gene counts matrix + individual counts for lengths)
  list(
    label = "Human (Hoang)",
    species = "Human",
    type = "gene_matrix",
    counts_path = file.path(ROOT, "RNA-seq/patient_RNAseq/results/GSE130970/counts/gene_counts_matrix.txt"),
    length_dir = file.path(ROOT, "RNA-seq/patient_RNAseq/results/GSE130970/counts/individual"),
    samplesheet = file.path(ROOT, "RNA-seq/patient_RNAseq/data/samplesheets/GSE130970_samplesheet.csv"),
    filter_col = "nafld_activity_score",
    filter_min = 4  # NAS >= 4 for "high"
  ),
  list(
    label = "Human (Govaere)",
    species = "Human",
    type = "gene_matrix",
    counts_path = file.path(ROOT, "RNA-seq/patient_RNAseq/results/GSE135251/counts/final_count_matrix_all_216_samples.txt"),
    # Use GSE130970's individual counts for lengths (same human gene annotation)
    length_dir = file.path(ROOT, "RNA-seq/patient_RNAseq/results/GSE130970/counts/individual"),
    samplesheet = file.path(ROOT, "RNA-seq/patient_RNAseq/data/samplesheets/GSE135251_samplesheet.csv"),
    filter_col = "disease",
    filter_exclude = "control"  # Exclude controls
  )
)

# =============================================================================
# TPM Computation Functions
# =============================================================================

read_featurecounts <- function(path) {
  df <- fread(path, sep = "\t", skip = "Geneid", header = TRUE)
  gene_col <- names(df)[1]
  length_col <- "Length"
  
  # Get metadata columns to exclude
  meta_cols <- c(gene_col, "Chr", "Start", "End", "Strand", "Length")
  sample_cols <- setdiff(names(df), meta_cols)
  
  gene_ids <- df[[gene_col]]
  lengths <- as.numeric(df[[length_col]])
  counts <- as.matrix(df[, ..sample_cols])
  rownames(counts) <- gene_ids
  
  list(gene_ids = gene_ids, lengths = lengths, counts = counts)
}

read_gene_counts_matrix <- function(path) {
  df <- fread(path, sep = "\t", header = TRUE)
  gene_col <- names(df)[1]
  gene_ids <- df[[gene_col]]
  
  # Extract count columns and convert to numeric matrix
  count_cols <- names(df)[-1]
  counts <- as.matrix(df[, ..count_cols])
  
  # Ensure numeric mode (data.table can produce character matrices)
  if (!is.numeric(counts)) {
    counts <- apply(counts, 2, as.numeric)
  }
  
  rownames(counts) <- gene_ids
  list(gene_ids = gene_ids, counts = counts)
}

read_lengths_from_individual <- function(dir_path) {
  files <- list.files(dir_path, pattern = "_counts\\.txt$", full.names = TRUE)
  if (length(files) == 0) stop("No individual count files found in ", dir_path)
  
  df <- fread(files[1], sep = "\t", skip = "Geneid", header = TRUE)
  gene_ids <- df[[1]]
  lengths <- as.numeric(df$Length)
  names(lengths) <- gene_ids
  lengths
}

extract_sample_id_from_path <- function(col) {
  # Extract sample ID from featureCounts path column
  basename <- basename(col)
  if (grepl("\\.Aligned", basename)) {
    return(sub("\\.Aligned.*", "", basename))
  }
  sub("\\..*", "", basename)
}

get_sample_filter <- function(config) {
  if (config$type == "featurecounts") {
    if (!file.exists(config$metadata_path)) return(NULL)
    meta <- fread(config$metadata_path, sep = "\t")
    if (!"sample_id" %in% names(meta)) return(NULL)
    if (!is.null(config$filter_col) && config$filter_col %in% names(meta)) {
      mask <- grepl(config$filter_pattern, meta[[config$filter_col]], ignore.case = TRUE)
      return(meta$sample_id[mask])
    }
    return(meta$sample_id)
  } else if (config$type == "gene_matrix") {
    if (!file.exists(config$samplesheet)) return(NULL)
    meta <- fread(config$samplesheet)
    sample_col <- if ("sample" %in% names(meta)) "sample" else if ("Run" %in% names(meta)) "Run" else names(meta)[1]
    
    if (!is.null(config$filter_min) && config$filter_col %in% names(meta)) {
      vals <- as.numeric(meta[[config$filter_col]])
      mask <- !is.na(vals) & vals >= config$filter_min
      return(meta[[sample_col]][mask])
    } else if (!is.null(config$filter_exclude) && config$filter_col %in% names(meta)) {
      mask <- !grepl(config$filter_exclude, meta[[config$filter_col]], ignore.case = TRUE)
      return(meta[[sample_col]][mask])
    }
    return(meta[[sample_col]])
  }
  NULL
}

subset_counts_by_samples <- function(counts, sample_ids, type) {
  if (is.null(sample_ids) || length(sample_ids) == 0) return(counts)
  
  if (type == "featurecounts") {
    col_ids <- sapply(colnames(counts), extract_sample_id_from_path)
    keep <- col_ids %in% sample_ids
  } else {
    keep <- colnames(counts) %in% sample_ids
  }
  
  if (sum(keep) == 0) return(counts)
  counts[, keep, drop = FALSE]
}

compute_tpm <- function(counts, lengths) {
  # Filter to shared genes with valid lengths
  valid_lengths <- !is.na(lengths) & lengths > 0
  
  if (is.null(names(lengths))) {
    # Unnamed lengths - match by position
    n_shared <- min(nrow(counts), length(lengths))
    valid_idx <- which(valid_lengths[1:n_shared])
    counts <- counts[valid_idx, , drop = FALSE]
    lengths_use <- lengths[valid_idx]
  } else {
    # Named lengths - match by gene name
    shared_genes <- intersect(rownames(counts), names(lengths)[valid_lengths])
    if (length(shared_genes) == 0) {
      stop("No shared genes between counts and lengths")
    }
    counts <- counts[shared_genes, , drop = FALSE]
    lengths_use <- lengths[shared_genes]
  }
  
  # Compute TPM
  length_kb <- lengths_use / 1000
  rpk <- counts / length_kb
  scale_factors <- colSums(rpk, na.rm = TRUE) / 1e6
  tpm <- t(t(rpk) / scale_factors)
  
  # Return mean TPM per gene
  rowMeans(tpm, na.rm = TRUE)
}

# =============================================================================
# Main Data Loading
# =============================================================================

load_all_tpm_data <- function() {
  all_data <- list()
  
  for (config in RAW_DATA_CONFIG) {
    message(sprintf("Processing: %s", config$label))
    
    if (!file.exists(config$counts_path)) {
      message(sprintf("  Skipping: counts file not found"))
      next
    }
    
    tryCatch({
      if (config$type == "featurecounts") {
        data <- read_featurecounts(config$counts_path)
        sample_ids <- get_sample_filter(config)
        counts <- subset_counts_by_samples(data$counts, sample_ids, config$type)
        tpm_mean <- compute_tpm(counts, data$lengths)
      } else if (config$type == "gene_matrix") {
        data <- read_gene_counts_matrix(config$counts_path)
        lengths <- read_lengths_from_individual(config$length_dir)
        sample_ids <- get_sample_filter(config)
        counts <- subset_counts_by_samples(data$counts, sample_ids, config$type)
        tpm_mean <- compute_tpm(counts, lengths)
      }
      
      tpm_vals <- as.numeric(tpm_mean)
      tpm_vals <- tpm_vals[!is.na(tpm_vals)]
      
      all_data[[config$label]] <- data.frame(
        dataset = config$label,
        species = config$species,
        tpm = tpm_vals,
        log2_tpm = log2(tpm_vals + 1),
        stringsAsFactors = FALSE
      )
      
      message(sprintf("  Loaded %d genes", length(tpm_vals)))
      
    }, error = function(e) {
      message(sprintf("  Error: %s", e$message))
    })
  }
  
  bind_rows(all_data)
}

# Compute sequencing depth statistics from raw count matrices
load_depth_stats <- function() {
  depth_data <- list()
  
  for (config in RAW_DATA_CONFIG) {
    message(sprintf("Computing depth stats: %s", config$label))
    
    if (!file.exists(config$counts_path)) {
      next
    }
    
    tryCatch({
      if (config$type == "featurecounts") {
        data <- read_featurecounts(config$counts_path)
        counts <- data$counts
      } else if (config$type == "gene_matrix") {
        data <- read_gene_counts_matrix(config$counts_path)
        counts <- data$counts
        # Ensure numeric for gene_matrix type (may be character from data.table)
        gene_names <- rownames(counts)
        counts <- apply(counts, 2, as.numeric)
        rownames(counts) <- gene_names
      }
      
      # Filter to relevant samples if specified
      sample_ids <- get_sample_filter(config)
      counts <- subset_counts_by_samples(counts, sample_ids, config$type)
      
      # Ensure counts is numeric matrix for calculations
      if (!is.numeric(counts)) {
        counts <- apply(counts, 2, as.numeric)
      }
      
      # Per-sample metrics
      library_sizes <- colSums(counts, na.rm = TRUE)              # Total counts per sample
      genes_detected <- colSums(counts > 0, na.rm = TRUE)         # Genes with > 0 counts per sample
      
      # Per-gene metrics (across all samples)
      mean_counts_per_gene <- rowMeans(counts)
      median_counts_per_gene <- median(mean_counts_per_gene, na.rm = TRUE)
      
      # Complexity: genes detected per million reads
      complexity_per_sample <- genes_detected / (library_sizes / 1e6)
      
      depth_data[[config$label]] <- data.frame(
        dataset = config$label,
        species = config$species,
        n_samples = ncol(counts),
        n_genes = nrow(counts),
        mean_library_size = mean(library_sizes),
        sd_library_size = sd(library_sizes),
        mean_genes_detected = mean(genes_detected),
        sd_genes_detected = sd(genes_detected),
        complexity_ratio = mean(complexity_per_sample),
        median_counts_per_gene = median_counts_per_gene,
        stringsAsFactors = FALSE
      )
      
      message(sprintf("  %d samples, %.1fM mean library size", 
                      ncol(counts), mean(library_sizes)/1e6))
      
    }, error = function(e) {
      message(sprintf("  Error: %s", e$message))
    })
  }
  
  bind_rows(depth_data)
}

# =============================================================================
# Theme and Plot Functions
# =============================================================================

theme_publication <- function(base_size = 14) {
  theme_minimal(base_size = base_size) +
    theme(
      text = element_text(family = "sans"),
      plot.title = element_text(size = base_size + 4, face = "bold", hjust = 0.5, margin = margin(b = 15)),
      plot.subtitle = element_text(size = base_size, hjust = 0.5, color = "grey40", margin = margin(b = 10)),
      axis.title = element_text(size = base_size, face = "bold"),
      axis.text = element_text(size = base_size - 2),
      legend.title = element_text(size = base_size, face = "bold"),
      legend.text = element_text(size = base_size - 2),
      legend.position = "bottom",
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      strip.text = element_text(size = base_size, face = "bold"),
      plot.margin = margin(20, 20, 20, 20),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}

compute_summary_stats <- function(data) {
  data %>%
    group_by(dataset, species) %>%
    summarize(
      n_genes = n(),
      median_tpm = median(tpm, na.rm = TRUE),
      mean_tpm = mean(tpm, na.rm = TRUE),
      q25_tpm = quantile(tpm, 0.25, na.rm = TRUE),
      q75_tpm = quantile(tpm, 0.75, na.rm = TRUE),
      pct_zero = mean(tpm == 0, na.rm = TRUE) * 100,
      pct_low = mean(tpm > 0 & tpm < 1, na.rm = TRUE) * 100,
      pct_med = mean(tpm >= 1 & tpm < 10, na.rm = TRUE) * 100,
      pct_high = mean(tpm >= 10, na.rm = TRUE) * 100,
      .groups = "drop"
    )
}

# Explicit dataset ordering: Human first, then Mouse
DATASET_ORDER <- c(
  "Human (Govaere)",
  "Human (Hoang)",
  "Cas13 mouse MCD",
  "Mouse MCD (Paquette)",
  "Mouse MCD (Yue)"
)

get_dataset_colors <- function(data) {
  # Manual mapping to match filtered/publication plots
  # Cas13 -> Yellow (#FFCD69) [1]
  # Paquette -> Teal (#2A8C7D) [2]
  # Yue -> Red (#ED5565) [3]
  # Hoang -> Green (#8BC163) [4]
  # Govaere -> Pink (#ED87BD) [5]
  
  manual_colors <- c(
    "Cas13 mouse MCD" = cat_palette[1],
    "Mouse MCD (Paquette)" = cat_palette[2],
    "Mouse MCD (Yue)" = cat_palette[3],
    "Human (Hoang)" = cat_palette[4],
    "Human (Govaere)" = cat_palette[5]
  )
  
  unique_datasets <- unique(data$dataset)
  colors <- manual_colors[unique_datasets]
  
  # Fallback for any missing keys
  if (any(is.na(colors))) {
    missing <- unique_datasets[is.na(colors)]
    warning(paste("Missing colors for:", paste(missing, collapse = ", ")))
    # Fill NAs with grey
    colors[is.na(colors)] <- "#888888"
  }
  
  colors
}

get_species_colors <- function() {
  c("Mouse" = cat_palette[2], "Human" = cat_palette[1])
}

# Ridgeline Plot
plot_ridgeline <- function(data, colors) {
  # Order datasets: Human first, then Mouse (reversed for ridgeline Y axis)
  ordered_levels <- intersect(DATASET_ORDER, unique(data$dataset))
  data$dataset <- factor(data$dataset, levels = rev(ordered_levels))
  
  # Cap values at 5 for visualization
  data$log2_tpm <- pmin(data$log2_tpm, 5)
  
  ggplot(data, aes(x = log2_tpm, y = dataset, fill = dataset)) +
    geom_density_ridges(
      alpha = 0.85,
      scale = 2.5,
      rel_min_height = 0.01,
      quantile_lines = TRUE,
      quantiles = 2,
      color = "white",
      linewidth = 0.8
    ) +
    geom_vline(xintercept = 5, linetype = "dashed", color = "grey50") +
    annotate("text", x = 4.8, y = length(unique(data$dataset)) + 0.5, 
             label = "Capped at 5", hjust = 1, vjust = 0, size = 4, fontface = "italic", color = "grey40") +
    scale_fill_manual(values = colors) +
    scale_x_continuous(
      name = expression(bold(log[2](TPM + 1))),
      limits = c(0, 5.2),
      breaks = seq(0, 5, 1),
      expand = c(0, 0)
    ) +
    labs(
      title = "TPM Distribution Across Datasets (Unfiltered)",
      subtitle = "Computed from raw count matrices - all genes included (Capped at 5)",
      y = NULL
    ) +
    theme_publication() +
    theme(
      legend.position = "none",
      axis.text.y = element_text(size = 12, face = "bold"),
      panel.grid.major.y = element_blank()
    )
}

# Improved Violin: Split zero vs non-zero (better for skewed data)
plot_violin_box <- function(data, colors) {
  # Order datasets: Human first, then Mouse
  ordered_levels <- intersect(DATASET_ORDER, unique(data$dataset))
  data$dataset <- factor(data$dataset, levels = ordered_levels)
  
  # For non-zero genes only
  nonzero_data <- data %>% filter(tpm > 0)
  nonzero_data$dataset <- factor(nonzero_data$dataset, levels = ordered_levels)
  
  # Create violin for non-zero genes with outlier points visible
  p <- ggplot(nonzero_data, aes(x = dataset, y = log2_tpm, fill = dataset)) +
    geom_violin(alpha = 0.7, color = "white", linewidth = 0.5, trim = FALSE, scale = "width") +
    geom_boxplot(width = 0.15, fill = "white", alpha = 0.9, color = "grey30", 
                 outlier.shape = 16, outlier.size = 0.8, outlier.alpha = 0.3) +
    stat_summary(fun = median, geom = "point", shape = 18, size = 3, color = "grey20") +
    scale_fill_manual(values = colors) +
    scale_y_continuous(
      name = expression(bold(log[2](TPM + 1))),
      limits = c(0, NA), breaks = seq(0, 16, 2)
    ) +
    labs(
      title = "Expression Distribution (Expressed Genes Only)",
      x = NULL
    ) +
    theme_publication() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 25, hjust = 1, size = 11, face = "bold"))
  
  p
}

# Additional plot: Zero vs Expressed bar comparison
plot_zero_comparison <- function(data, colors) {
  # Calculate zero and expressed counts
  summary_df <- data %>%
    group_by(dataset, species) %>%
    summarize(
      n_zero = sum(tpm == 0),
      n_expressed = sum(tpm > 0),
      pct_zero = mean(tpm == 0) * 100,
      .groups = "drop"
    ) %>%
    pivot_longer(
      cols = c(n_zero, n_expressed),
      names_to = "category",
      values_to = "count"
    ) %>%
    mutate(category = factor(category, levels = c("n_zero", "n_expressed"),
                             labels = c("Zero Expression", "Expressed")))
  
  ggplot(summary_df, aes(x = dataset, y = count / 1000, fill = category)) +
    geom_col(position = "dodge", width = 0.7, color = "white", linewidth = 0.3) +
    scale_fill_manual(
      values = c("Zero Expression" = cat_palette[18], "Expressed" = cat_palette[2]),
      name = "Expression Status"
    ) +
    scale_y_continuous(name = "Number of Genes (thousands)", expand = c(0, 0)) +
    labs(
      title = "Zero vs Expressed Gene Counts",
      subtitle = "Absolute gene counts per dataset",
      x = NULL
    ) +
    theme_publication() +
    theme(
      axis.text.x = element_text(angle = 25, hjust = 1, size = 11, face = "bold"),
      legend.position = "right"
    )
}

# Species Density
plot_species_density <- function(data) {
  species_colors <- get_species_colors()
  
  ggplot(data, aes(x = log2_tpm, fill = species, color = species)) +
    geom_density(alpha = 0.5, linewidth = 1.2) +
    scale_fill_manual(values = species_colors, name = "Species") +
    scale_color_manual(values = species_colors, name = "Species") +
    scale_x_continuous(name = expression(bold(log[2](TPM + 1))), limits = c(0, 12), breaks = seq(0, 12, 2)) +
    labs(title = "Expression Distribution by Species (Unfiltered)", subtitle = "Aggregated density: Mouse MCD vs Human MASLD", y = "Density") +
    theme_publication() +
    theme(legend.position = c(0.85, 0.85), legend.background = element_rect(fill = "white", color = "grey80", linewidth = 0.3))
}

# ECDF
plot_ecdf <- function(data, colors) {
  ggplot(data, aes(x = log2_tpm, color = dataset)) +
    stat_ecdf(linewidth = 1.2, alpha = 0.9) +
    scale_color_manual(values = colors, name = "Dataset") +
    scale_x_continuous(name = expression(bold(log[2](TPM + 1))), limits = c(0, 12), breaks = seq(0, 12, 2)) +
    scale_y_continuous(name = "Cumulative Proportion", labels = percent_format(), breaks = seq(0, 1, 0.2)) +
    labs(title = "Cumulative Expression Distribution (Unfiltered)", subtitle = "ECDF from raw count matrices") +
    theme_publication() +
    theme(legend.position = c(0.85, 0.35), legend.background = element_rect(fill = "white", color = "grey80", linewidth = 0.3))
}

# Expression Bins
plot_expression_bins <- function(summary_stats) {
  # Order datasets: Human first, then Mouse
  ordered_levels <- intersect(DATASET_ORDER, summary_stats$dataset)
  summary_stats$dataset <- factor(summary_stats$dataset, levels = ordered_levels)
  
  bins_df <- summary_stats %>%
    select(dataset, pct_zero, pct_low, pct_med, pct_high) %>%
    pivot_longer(cols = starts_with("pct_"), names_to = "bin", values_to = "percentage") %>%
    mutate(bin = factor(bin, levels = c("pct_zero", "pct_low", "pct_med", "pct_high"),
                        labels = c("Zero (TPM = 0)", "Low (0 < TPM < 1)", "Medium (1 ≤ TPM < 10)", "High (TPM ≥ 10)")))
  
  # Grey for Zero, then orange → magenta gradient for expressed genes
  bin_colors <- c(
    "Zero (TPM = 0)" = "#A0A0A0",        # Grey (no expression)
    "Low (0 < TPM < 1)" = "#F2A45E",     # Orange
    "Medium (1 ≤ TPM < 10)" = "#E07BB6", # Light pink/magenta
    "High (TPM ≥ 10)" = "#C23B75"        # Magenta (highest expression)
  )
  
  ggplot(bins_df, aes(x = dataset, y = percentage, fill = bin)) +
    geom_col(position = "stack", width = 0.7, color = "white", linewidth = 0.3) +
    scale_fill_manual(values = bin_colors, name = "Expression Level") +
    scale_y_continuous(name = "Percentage of Genes", labels = percent_format(scale = 1), expand = c(0, 0)) +
    labs(title = "Gene Expression Level Distribution (Unfiltered)", subtitle = "Proportion of genes in each expression category", x = NULL) +
    theme_publication() +
    theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 11, face = "bold"), legend.position = "right")
}

# Summary Heatmap
plot_summary_heatmap <- function(summary_stats) {
  metrics_df <- summary_stats %>%
    select(dataset, median_tpm, mean_tpm, pct_zero, pct_high) %>%
    pivot_longer(cols = -dataset, names_to = "metric", values_to = "value") %>%
    mutate(metric = factor(metric, levels = c("median_tpm", "mean_tpm", "pct_zero", "pct_high"),
                           labels = c("Median TPM", "Mean TPM", "% Zero Expression", "% High Expression"))) %>%
    group_by(metric) %>%
    mutate(scaled = (value - min(value)) / (max(value) - min(value) + 1e-9)) %>%
    ungroup()
  
  ggplot(metrics_df, aes(x = metric, y = dataset, fill = scaled)) +
    geom_tile(color = "white", linewidth = 1.5) +
    geom_text(aes(label = ifelse(grepl("%", metric), sprintf("%.1f%%", value), sprintf("%.2f", value))),
              size = 4, fontface = "bold", color = ifelse(metrics_df$scaled > 0.5, "white", "grey20")) +
    # Yellow → Orange → Magenta gradient for intuitive low-to-high scaling
    scale_fill_gradientn(colors = c("#FFDD78", "#F2A45E", "#E07BB6", "#C23B75"), guide = "none") +
    labs(title = "Expression Summary Metrics (Unfiltered)", subtitle = "Key statistics from raw count matrices", x = NULL, y = NULL) +
    theme_publication() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 11, face = "bold"),
          axis.text.y = element_text(size = 12, face = "bold"), panel.grid = element_blank())
}

# Sequencing Depth Comparison Plot
# This explains TPM differences by showing underlying library metrics
plot_sequencing_depth <- function(depth_stats, colors) {
  # Order datasets: Human first, then Mouse
  ordered_levels <- intersect(DATASET_ORDER, depth_stats$dataset)
  depth_stats$dataset <- factor(depth_stats$dataset, levels = ordered_levels)
  
  # Panel 1: Library Size (total counts per sample)
  p1 <- ggplot(depth_stats, aes(x = dataset, y = mean_library_size / 1e6, fill = dataset)) +
    geom_col(width = 0.7, color = "white", linewidth = 0.3) +
    geom_errorbar(aes(ymin = (mean_library_size - sd_library_size) / 1e6,
                      ymax = (mean_library_size + sd_library_size) / 1e6),
                  width = 0.2, color = "grey30") +
    scale_fill_manual(values = colors) +
    scale_y_continuous(name = "Mean Library Size\n(millions of reads)", expand = c(0, 0)) +
    labs(title = "A. Sequencing Depth", x = NULL) +
    theme_publication() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 35, hjust = 1, size = 9))
  
  # Panel 2: Genes Detected
  p2 <- ggplot(depth_stats, aes(x = dataset, y = mean_genes_detected / 1000, fill = dataset)) +
    geom_col(width = 0.7, color = "white", linewidth = 0.3) +
    geom_errorbar(aes(ymin = (mean_genes_detected - sd_genes_detected) / 1000,
                      ymax = (mean_genes_detected + sd_genes_detected) / 1000),
                  width = 0.2, color = "grey30") +
    scale_fill_manual(values = colors) +
    scale_y_continuous(name = "Mean Genes Detected\n(thousands)", expand = c(0, 0)) +
    labs(title = "B. Gene Detection", x = NULL) +
    theme_publication() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 35, hjust = 1, size = 9))
  
  # Panel 3: Library Complexity (genes per million reads)
  p3 <- ggplot(depth_stats, aes(x = dataset, y = complexity_ratio, fill = dataset)) +
    geom_col(width = 0.7, color = "white", linewidth = 0.3) +
    scale_fill_manual(values = colors) +
    scale_y_continuous(name = "Library Complexity\n(genes per M reads)", expand = c(0, 0)) +
    labs(title = "C. Library Complexity", x = NULL) +
    theme_publication() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 35, hjust = 1, size = 9))
  
  # Panel 4: Median counts per gene
  p4 <- ggplot(depth_stats, aes(x = dataset, y = median_counts_per_gene, fill = dataset)) +
    geom_col(width = 0.7, color = "white", linewidth = 0.3) +
    scale_fill_manual(values = colors) +
    scale_y_continuous(name = "Median Counts\nper Gene", expand = c(0, 0)) +
    labs(title = "D. Count Distribution", x = NULL) +
    theme_publication() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 35, hjust = 1, size = 9))
  
  (p1 | p2) / (p3 | p4) +
    plot_annotation(
      title = "Sequencing Depth and Library Metrics by Dataset",
      subtitle = "Explaining TPM distribution differences through underlying count statistics",
      theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
                    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "grey40"))
    )
}

# Multi-panel Summary
create_summary_figure <- function(data, summary_stats, colors) {
  p1 <- plot_ridgeline(data, colors)
  p2 <- plot_violin_box(data, colors)
  p3 <- plot_species_density(data)
  p4 <- plot_expression_bins(summary_stats)
  
  (p1 | p2) / (p3 | p4) +
    plot_annotation(
      title = "TPM Expression Distribution Analysis (Unfiltered)",
      subtitle = "Computed from raw count matrices - all genes included for fair comparison",
      theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
                    plot.subtitle = element_text(size = 14, hjust = 0.5, color = "grey40"))
    ) +
    plot_layout(heights = c(1, 1))
}

# =============================================================================
# Main Execution
# =============================================================================

main <- function() {
  message("=== TPM Distribution Plots (Unfiltered Raw Counts) ===\n")
  message("Loading and computing TPM from raw count matrices...")
  
  data <- load_all_tpm_data()
  
  if (nrow(data) == 0) stop("No data loaded. Check file paths.")
  
  message(sprintf("\nTotal: %d observations from %d datasets", nrow(data), length(unique(data$dataset))))
  
  summary_stats <- compute_summary_stats(data)
  colors <- get_dataset_colors(data)
  
  # Save stats
  stats_path <- file.path(OUTPUT_DIR, "tpm_distribution_unfiltered_summary_stats.tsv")
  write.table(summary_stats, stats_path, sep = "\t", row.names = FALSE, quote = FALSE)
  message(sprintf("\nSaved summary statistics to %s", stats_path))
  
  # Print summary
  message("\nDataset Summary:")
  print(as.data.frame(summary_stats[, c("dataset", "species", "n_genes", "pct_zero")]))
  
  message("\nGenerating plots...")
  
  pdf(file.path(OUTPUT_DIR, "tpm_distribution_unfiltered_ridgeline.pdf"), width = 10, height = 7)
  print(plot_ridgeline(data, colors)); dev.off()
  message("  - Ridgeline saved")
  
  pdf(file.path(OUTPUT_DIR, "tpm_distribution_unfiltered_violin.pdf"), width = 10, height = 7)
  print(plot_violin_box(data, colors)); dev.off()
  message("  - Violin saved")
  
  pdf(file.path(OUTPUT_DIR, "tpm_distribution_unfiltered_species.pdf"), width = 10, height = 6)
  print(plot_species_density(data)); dev.off()
  message("  - Species comparison saved")
  

  pdf(file.path(OUTPUT_DIR, "tpm_distribution_unfiltered_bins.pdf"), width = 10, height = 7)
  print(plot_expression_bins(summary_stats)); dev.off()
  message("  - Expression bins saved")
  
  pdf(file.path(OUTPUT_DIR, "tpm_distribution_unfiltered_zero_comparison.pdf"), width = 10, height = 7)
  print(plot_zero_comparison(data, colors)); dev.off()
  message("  - Zero comparison saved")

  
  pdf(file.path(OUTPUT_DIR, "tpm_distribution_unfiltered_heatmap.pdf"), width = 10, height = 6)
  print(plot_summary_heatmap(summary_stats)); dev.off()
  message("  - Summary heatmap saved")
  
  # Sequencing depth comparison plot
  message("\nComputing sequencing depth statistics...")
  depth_stats <- load_depth_stats()
  
  pdf(file.path(OUTPUT_DIR, "tpm_distribution_unfiltered_depth_metrics.pdf"), width = 12, height = 10)
  print(plot_sequencing_depth(depth_stats, colors)); dev.off()
  message("  - Sequencing depth metrics saved")
  
  # Save depth stats
  depth_path <- file.path(OUTPUT_DIR, "sequencing_depth_stats.tsv")
  write.table(depth_stats, depth_path, sep = "\t", row.names = FALSE, quote = FALSE)
  message(sprintf("  - Depth stats table saved to %s", depth_path))
  
  pdf(file.path(OUTPUT_DIR, "tpm_distribution_unfiltered_summary.pdf"), width = 16, height = 12)
  print(create_summary_figure(data, summary_stats, colors)); dev.off()
  message("  - Comprehensive summary saved")
  
  message("\n=== All plots saved to: ", OUTPUT_DIR, " ===")
  message("Done!")
}

if (sys.nframe() == 0) {
  main()
}
