#!/usr/bin/env Rscript

# Generic DESeq2 driver that supports multiple study-specific groupings.
# Usage:
#   Rscript deseq2_generic_analysis.R \
#       <count_matrix.tsv> \
#       <metadata.csv> \
#       <output_dir> \
#       <group_specs> \
#       [sample_column] \
#       [min_group_size]
# Example group specs:
#   "fibrosis_stage=fibrosis,steatosis_grade=steatosis,nas_score=nas"

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(matrixStats)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
})

has_plotly <- requireNamespace("plotly", quietly = TRUE)
has_htmlwidgets <- requireNamespace("htmlwidgets", quietly = TRUE)
enable_3d_pca <- FALSE

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript deseq2_generic_analysis.R <count_matrix> <metadata> <output_dir> <group_specs> [sample_column] [min_group_size] [covariates]")
}

count_matrix_file <- args[1]
metadata_file <- args[2]
output_dir <- args[3]
group_specs <- args[4]
sample_column <- if (length(args) >= 5) tolower(args[5]) else "run"
min_group_size <- if (length(args) >= 6) as.integer(args[6]) else 3
covariates <- if (length(args) >= 7) args[7] else ""
covariate_list <- covariates |> strsplit(",") |> unlist() |> trimws()
covariate_list <- covariate_list[covariate_list != ""]

if (!file.exists(count_matrix_file)) {
  stop("Count matrix not found: ", count_matrix_file)
}
if (!file.exists(metadata_file)) {
  stop("Metadata file not found: ", metadata_file)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading count matrix: ", count_matrix_file)
count_df <- read.delim(count_matrix_file, header = TRUE, check.names = FALSE, comment.char = "")

if (!"Geneid" %in% colnames(count_df) && !"GeneID" %in% colnames(count_df)) {
  stop("Count matrix must have a Geneid/GeneID column")
}

gene_col <- if ("Geneid" %in% colnames(count_df)) "Geneid" else "GeneID"
rownames(count_df) <- count_df[[gene_col]]
count_df[[gene_col]] <- NULL

# Drop helper columns produced by featureCounts if present
annotation_cols <- c("Chr", "Start", "End", "Strand", "Length")
annotation_cols <- annotation_cols[annotation_cols %in% colnames(count_df)]
if (length(annotation_cols) > 0) {
  count_df <- count_df[, setdiff(colnames(count_df), annotation_cols), drop = FALSE]
}

# Convert to numeric matrix and drop rows that are entirely NA (e.g. BAM path rows)
count_df[] <- lapply(count_df, function(col) {
  if (is.character(col)) {
    suppressWarnings(as.numeric(col))
  } else {
    col
  }
})
count_matrix <- as.matrix(count_df)
count_matrix <- round(count_matrix)
mode(count_matrix) <- "integer"
all_na <- apply(count_matrix, 1, function(x) all(is.na(x)))
if (any(all_na)) {
  count_matrix <- count_matrix[!all_na, , drop = FALSE]
}

sample_ids <- colnames(count_matrix)

message("Loaded count matrix with ", nrow(count_matrix), " genes and ", length(sample_ids), " samples")

message("Reading metadata: ", metadata_file)
metadata <- as.data.frame(readr::read_csv(metadata_file, show_col_types = FALSE))
colnames(metadata) <- make.names(tolower(colnames(metadata)))

if (!sample_column %in% colnames(metadata)) {
  stop("Sample column '", sample_column, "' not found in metadata. Available columns: ", paste(colnames(metadata), collapse = ", "))
}

metadata$sample_id <- as.character(metadata[[sample_column]])

common_samples <- intersect(sample_ids, metadata$sample_id)
if (length(common_samples) < 2) {
  stop("Fewer than 2 overlapping samples between counts and metadata")
}

# Reorder counts and metadata to the same sample order
count_matrix <- count_matrix[, common_samples, drop = FALSE]
metadata <- metadata[match(common_samples, metadata$sample_id), , drop = FALSE]
rownames(metadata) <- metadata$sample_id

message("Samples in common: ", length(common_samples))

# --- Custom Derived Columns for Strict Comparisons ---
if ("fibrosis_stage" %in% colnames(metadata)) {
  metadata$fibrosis_strict <- NA_character_
  
  # Logic for GSE130970 (using NAS)
  if ("nafld_activity_score" %in% colnames(metadata)) {
    nas <- suppressWarnings(as.numeric(metadata$nafld_activity_score))
    fib <- suppressWarnings(as.numeric(metadata$fibrosis_stage))
    
    # Healthy = F0 AND NAS=0
    is_healthy <- !is.na(nas) & nas == 0 & !is.na(fib) & fib == 0
    # Fibrosis = Stage >= 1
    is_fibrosis <- !is.na(fib) & fib >= 1
    
    metadata$fibrosis_strict[is_healthy] <- "Healthy"
    metadata$fibrosis_strict[is_fibrosis] <- "Fibrosis"
    
    # Explicitly calculate F0_Disease for verifying exclusion
    is_f0_disease <- !is.na(fib) & fib == 0 & !is.na(nas) & nas > 0
    if (any(is_f0_disease)) {
      message("Deriving fibrosis_strict: Excluding ", sum(is_f0_disease), " F0 samples with NAS > 0 (Early Disease) from Controls.")
    }
  } 
  # Logic for GSE135251 (using Stage)
  else if ("stage" %in% colnames(metadata)) {
    stg <- tolower(trimws(as.character(metadata$stage)))
    fib <- suppressWarnings(as.numeric(metadata$fibrosis_stage))
    
    # Healthy = Stage "control"
    is_healthy <- !is.na(stg) & stg == "control"
    # Fibrosis = Stage >= 1
    is_fibrosis <- !is.na(fib) & fib >= 1
    
    metadata$fibrosis_strict[is_healthy] <- "Healthy"
    metadata$fibrosis_strict[is_fibrosis] <- "Fibrosis"
    
    # Check exclusions
    is_f0_excluded <- !is.na(fib) & fib == 0 & !is_healthy
    if (any(is_f0_excluded)) {
      message("Deriving fibrosis_strict: Excluding ", sum(is_f0_excluded), " F0 samples not labeled 'control' from Controls.")
    }
  }
}
# -----------------------------------------------------

# No raw-count prefilter
count_matrix_filtered <- count_matrix
message("No raw-count prefilter applied; genes retained: ", nrow(count_matrix_filtered))

# Prepare base dataset for normalized counts
base_coldata <- DataFrame(row.names = colnames(count_matrix_filtered))
dds_base <- DESeqDataSetFromMatrix(
  countData = count_matrix_filtered,
  colData = base_coldata,
  design = ~1
)
dds_base <- estimateSizeFactors(dds_base)
normalized_counts <- counts(dds_base, normalized = TRUE)
write.csv(normalized_counts, file.path(output_dir, "deseq2_normalized_counts.csv"), quote = FALSE)

# Helper functions -----------------------------------------------------------

safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  ifelse(nchar(x) == 0, "level", x)
}

create_group_factor <- function(values, alias) {
  if (length(values) == 0 || all(is.na(values))) {
    return(NULL)
  }

  values_vec <- as.vector(values)
  names(values_vec) <- names(values)
  n <- length(values_vec)

  make_factor <- function(labels, levels) {
    factor(labels, levels = levels)
  }

  result_labels <- rep(NA_character_, n)
  level_order <- character()

  if (alias %in% c("nas_high")) {
    numeric_vals <- suppressWarnings(as.integer(values_vec))
    control <- numeric_vals == 0
    high <- numeric_vals >= 4
    result_labels[control] <- "NAS_0"
    result_labels[high] <- "NAS_4plus"
    level_order <- c("NAS_0", "NAS_4plus")
  } else if (alias %in% c("nas_low")) {
    numeric_vals <- suppressWarnings(as.integer(values_vec))
    control <- numeric_vals == 0
    low <- numeric_vals >= 1 & numeric_vals <= 3
    result_labels[control] <- "NAS_0"
    result_labels[low] <- "NAS_1to3"
    level_order <- c("NAS_0", "NAS_1to3")
  } else if (alias %in% c("ballooning", "ballooning_grade", "cytological_ballooning_grade")) {
    numeric_vals <- suppressWarnings(as.integer(values_vec))
    valid <- !is.na(numeric_vals)
    result_labels[valid] <- paste0("B", numeric_vals[valid])
    level_order <- paste0("B", sort(unique(numeric_vals[valid])))
  } else if (alias %in% c("masld", "nas_any")) {
    numeric_vals <- suppressWarnings(as.integer(values_vec))
    control <- numeric_vals == 0
    masld <- numeric_vals >= 1
    result_labels[control] <- "Control"
    result_labels[masld] <- "MASLD"
    level_order <- c("Control", "MASLD")
  } else if (alias %in% c("fibrosis", "fibrosis_stage", "fibrosis_any", "fibrosis_binary")) {
    numeric_vals <- suppressWarnings(as.integer(values_vec))
    result_labels[numeric_vals == 0] <- "F0"
    result_labels[numeric_vals >= 1] <- "F1to4"
    level_order <- c("F0", "F1to4")
  } else if (alias %in% c("ballooning_any", "ballooning")) {
    numeric_vals <- suppressWarnings(as.integer(values_vec))
    result_labels[numeric_vals == 0] <- "B0"
    result_labels[numeric_vals >= 1] <- "B1to2"
    level_order <- c("B0", "B1to2")
  } else if (alias %in% c("fibrosis_trim")) {
    numeric_vals <- suppressWarnings(as.integer(values_vec))
    result_labels[numeric_vals == 0] <- "F0"
    result_labels[numeric_vals == 2] <- "F2"
    result_labels[numeric_vals == 3] <- "F3"
    result_labels[numeric_vals == 4] <- "F4"
    level_order <- c("F0", "F2", "F3", "F4")
  } else if (alias %in% c("steatosis_any", "steatosis_binary")) {
    numeric_vals <- suppressWarnings(as.integer(values_vec))
    result_labels[numeric_vals == 0] <- "S0"
    result_labels[numeric_vals >= 1] <- "S1plus"
    level_order <- c("S0", "S1plus")
  } else if (alias %in% c("steatosis", "steatosis_grade")) {
    numeric_vals <- suppressWarnings(as.integer(values_vec))
    valid <- !is.na(numeric_vals)
    result_labels[valid] <- paste0("S", numeric_vals[valid])
    level_order <- paste0("S", sort(unique(numeric_vals[valid])))
  } else if (alias %in% c("steatosis_any")) {
    numeric_vals <- suppressWarnings(as.integer(values_vec))
    result_labels[numeric_vals == 0] <- "S0"
    result_labels[numeric_vals >= 1] <- "S1plus"
    level_order <- c("S0", "S1plus")
  } else if (alias %in% c("inflammation", "lobular_inflammation")) {
    numeric_vals <- suppressWarnings(as.integer(values_vec))
    valid <- !is.na(numeric_vals)
    result_labels[valid] <- paste0("L", numeric_vals[valid])
    level_order <- paste0("L", sort(unique(numeric_vals[valid])))
  } else if (alias %in% c("inflammation_any")) {
    numeric_vals <- suppressWarnings(as.integer(values_vec))
    result_labels[numeric_vals == 0] <- "L0"
    result_labels[numeric_vals >= 1] <- "L1to2"
    level_order <- c("L0", "L1to2")
  } else if (alias %in% c("nas", "nas_score", "nafld_activity_score")) {
    numeric_vals <- suppressWarnings(as.integer(values_vec))
    valid <- !is.na(numeric_vals)
    result_labels[valid] <- as.character(numeric_vals[valid])
    level_order <- as.character(sort(unique(numeric_vals[valid])))
  } else if (alias %in% c("disease")) {
    values_chr <- trimws(as.character(values_vec))
    values_chr[values_chr == ""] <- NA_character_
    result_labels <- values_chr
    level_order <- sort(unique(values_chr[!is.na(values_chr)]))
  } else if (alias %in% c("stage")) {
    values_chr <- tolower(trimws(as.character(values_vec)))
    values_chr[values_chr == ""] <- NA_character_
    result_labels <- values_chr
    desired_order <- c("control", "early", "moderate", "advanced", "cirrhosis", "severe")
    present <- desired_order[desired_order %in% result_labels]
    if (length(present) == 0) {
      level_order <- sort(unique(result_labels[!is.na(result_labels)]))
    } else {
      level_order <- present
    }
  } else if (is.numeric(values_vec)) {
    numeric_vals <- suppressWarnings(as.numeric(values_vec))
    valid <- !is.na(numeric_vals)
    result_labels[valid] <- as.character(numeric_vals[valid])
    level_order <- as.character(sort(unique(numeric_vals[valid])))
  } else if (alias %in% c("fibrosis_strict")) {
    values_chr <- as.character(values_vec)
    result_labels <- values_chr
    level_order <- c("Healthy", "Fibrosis")
  } else {
    values_chr <- as.character(values_vec)
    result_labels <- values_chr
    level_order <- sort(unique(values_chr[!is.na(values_chr)]))
  }

  if (length(level_order) > 0) {
    level_order <- level_order[level_order %in% result_labels]
  }
  out <- make_factor(result_labels, levels = level_order)
  if (!is.null(names(values))) {
    names(out) <- names(values)
  }
  out
}

group_color_map <- function(levels_vec) {
  control_levels <- c("control", "nas_0", "f0", "s0", "control_patients", "control patients", "healthy")
  control_levels <- unique(c(control_levels, tolower(control_levels)))
  lvl_lower <- tolower(levels_vec)
  colors <- rep(NA_character_, length(levels_vec))
  control_idx <- lvl_lower %in% control_levels
  if (any(control_idx)) {
    colors[control_idx] <- "#9c9c9c" # gray
    colors[!control_idx] <- "#ED87BD" # magenta/pink
  } else if (length(levels_vec) == 2) {
    colors[1] <- "#9c9c9c"
    colors[2] <- "#ED87BD"
  } else {
    colors <- RColorBrewer::brewer.pal(min(8, max(3, length(levels_vec))), "Set2")
  }
  names(colors) <- levels_vec
  colors
}

compute_pca_scores <- function(vsd, ntop = 500) {
  mat <- assay(vsd)
  if (nrow(mat) == 0 || ncol(mat) < 2) {
    return(NULL)
  }

  mat <- mat - rowMeans(mat)
  rv <- matrixStats::rowVars(mat)
  ntop <- min(ntop, length(rv))
  select <- order(rv, decreasing = TRUE)[seq_len(ntop)]
  pca <- prcomp(t(mat[select, , drop = FALSE]))

  percent_var <- (pca$sdev^2) / sum(pca$sdev^2) * 100
  pcs_to_keep <- min(3, ncol(pca$x))
  scores <- as.data.frame(pca$x[, seq_len(pcs_to_keep), drop = FALSE])
  colnames(scores) <- paste0("PC", seq_len(pcs_to_keep))
  scores$sample_id <- rownames(pca$x)

  list(scores = scores, percent_var = percent_var)
}

save_plotly_pdf_if_available <- function(plt, pdf_file, width = 900, height = 700) {
  if (is.null(pdf_file) || !has_plotly) {
    return(invisible(NULL))
  }
  if (!requireNamespace("plotly", quietly = TRUE)) {
    message("plotly not available; skipping PDF export for ", pdf_file)
    return(invisible(NULL))
  }

  saved <- FALSE

  if ("save_image" %in% getNamespaceExports("plotly")) {
    tryCatch({
      plotly::save_image(plt, pdf_file, format = "pdf", width = width, height = height)
      saved <- TRUE
    }, error = function(e) {
      message("Failed to export PDF via save_image for ", pdf_file, ": ", conditionMessage(e))
    })
  }

  if (!saved && "orca" %in% getNamespaceExports("plotly")) {
    tryCatch({
      plotly::orca(plt, pdf_file, width = width, height = height)
      saved <- TRUE
    }, error = function(e) {
      message("Failed to export PDF via orca for ", pdf_file, ": ", conditionMessage(e))
    })
  }

  if (!saved) {
    message("No plotly PDF export backend available; skipping ", pdf_file)
  }
}

save_3d_pca_plot <- function(pca_scores, percent_var, group_factor, title, out_file, pdf_file = NULL) {
  if (!has_plotly || !has_htmlwidgets) {
    message("plotly/htmlwidgets not available; skipping 3D PCA for ", title)
    return(invisible(NULL))
  }

  required_cols <- c("PC1", "PC2", "PC3", "sample_id")
  if (!all(required_cols %in% colnames(pca_scores))) {
    message("Not enough principal components for 3D PCA: ", title)
    return(invisible(NULL))
  }

  plot_df <- pca_scores
  plot_df$group <- group_factor[match(plot_df$sample_id, names(group_factor))]
  plot_df <- plot_df[!is.na(plot_df$group), , drop = FALSE]
  plot_df$group <- droplevels(plot_df$group)
  if (nrow(plot_df) < 3 || nlevels(plot_df$group) < 2) {
    message("Insufficient samples for 3D PCA: ", title)
    return(invisible(NULL))
  }

  # Guard against missing percent_var entries
  pv <- percent_var
  if (length(pv) < 3) {
    pv <- c(pv, rep(NA_real_, 3 - length(pv)))
  }

  plt <- plotly::plot_ly(
    plot_df,
    x = ~PC1,
    y = ~PC2,
    z = ~PC3,
    color = ~group,
    colors = RColorBrewer::brewer.pal(min(8, max(3, length(unique(plot_df$group)))), "Set1")
  ) %>%
    plotly::add_markers(size = 6, opacity = 0.85) %>%
    plotly::layout(
      title = title,
      scene = list(
        xaxis = list(title = paste0("PC1 (", round(pv[1], 1), "%)")),
        yaxis = list(title = paste0("PC2 (", round(pv[2], 1), "%)")),
        zaxis = list(title = paste0("PC3 (", round(pv[3], 1), "%)"))
      ),
      legend = list(itemsizing = "constant")
    )

  libdir <- file.path(dirname(out_file), paste0(tools::file_path_sans_ext(basename(out_file)), "_libs"))
  htmlwidgets::saveWidget(plt, out_file, selfcontained = FALSE, libdir = libdir)
  save_plotly_pdf_if_available(plt, pdf_file)
}

compute_separation_metrics <- function(pca_scores, group_factor, alias) {
  if (is.null(pca_scores) || !"sample_id" %in% colnames(pca_scores)) {
    return(NULL)
  }

  pc_cols <- intersect(colnames(pca_scores), c("PC1", "PC2", "PC3"))
  if (length(pc_cols) < 2) {
    return(NULL)
  }

  df <- pca_scores[, c(pc_cols, "sample_id"), drop = FALSE]
  df$group <- group_factor[match(df$sample_id, names(group_factor))]
  df <- df[!is.na(df$group), , drop = FALSE]
  df$group <- droplevels(df$group)

  if (nrow(df) < 3 || nlevels(df$group) < 2) {
    return(NULL)
  }

  centroids <- aggregate(df[, pc_cols, drop = FALSE], by = list(df$group), FUN = mean)
  between <- if (nrow(centroids) > 1) mean(dist(as.matrix(centroids[, pc_cols, drop = FALSE]))) else NA_real_

  within_vals <- sapply(levels(df$group), function(g) {
    pts <- as.matrix(df[df$group == g, pc_cols, drop = FALSE])
    if (nrow(pts) <= 1) return(NA_real_)
    mean(dist(pts))
  })
  within <- mean(within_vals, na.rm = TRUE)

  sep_ratio <- if (is.na(between) || is.na(within) || within == 0) NA_real_ else between / within

  data.frame(
    alias = alias,
    groups = paste(levels(df$group), collapse = "|"),
    n_samples = nrow(df),
    mean_between_distance = between,
    mean_within_distance = within,
    separation_ratio = sep_ratio,
    stringsAsFactors = FALSE
  )
}

plot_heatmap <- function(vsd, alias, group_factor, out_file) {
  mat <- assay(vsd)
  if (nrow(mat) == 0) {
    return(invisible(NULL))
  }
  top <- head(order(matrixStats::rowVars(mat), decreasing = TRUE), n = min(50, nrow(mat)))
  if (length(top) == 0) {
    return(invisible(NULL))
  }
  heatmap_mat <- mat[top, , drop = FALSE]
  heatmap_mat <- heatmap_mat - rowMeans(heatmap_mat)
  annot <- data.frame(Group = group_factor)
  rownames(annot) <- names(group_factor)
  group_cols <- group_color_map(levels(group_factor))
  pdf(out_file, width = 10, height = 8)
  pheatmap(
    heatmap_mat,
    annotation_col = annot,
    show_rownames = nrow(heatmap_mat) <= 40,
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(255),
    annotation_colors = list(Group = group_cols),
    border_color = NA,
    main = paste("Sample Heatmap -", alias)
  )
  dev.off()
}

generate_volcano <- function(res_df, comparison_name, out_file) {
  res_plot <- res_df %>%
    filter(!is.na(padj)) %>%
    mutate(
      neg_log10_padj = -log10(padj),
      significance = case_when(
        padj <= 0.05 & log2FoldChange >= 1 ~ "Up",
        padj <= 0.05 & log2FoldChange <= -1 ~ "Down",
        TRUE ~ "Not Significant"
      )
    )
  p <- ggplot(res_plot, aes(x = log2FoldChange, y = neg_log10_padj, color = significance)) +
    geom_point(alpha = 0.6, size = 1.2) +
    scale_color_manual(values = c("Up" = "#d73027", "Down" = "#4575b4", "Not Significant" = "grey70")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
    theme_bw() +
    labs(
      title = paste("Volcano Plot -", comparison_name),
      x = "log2 Fold Change",
      y = "-log10 adjusted p-value"
    )
  pdf(out_file, width = 8, height = 6)
  print(p)
  dev.off()
}

generate_comparison_outputs <- function(dds, vsd, group_factor, level_test, level_ref, base_dir) {
  comparison_name <- paste0(level_test, "_vs_", level_ref)
  safe_comp <- safe_name(comparison_name)
  comparison_dir <- file.path(base_dir, safe_comp)
  dir.create(comparison_dir, recursive = TRUE, showWarnings = FALSE)

  res <- results(dds, contrast = c("group", level_test, level_ref))
  
  # Apply LFC shrinkage for more conservative fold change estimates
  shrink_type <- if (requireNamespace("ashr", quietly = TRUE)) "ashr" else "normal"
  message(sprintf("    Applying LFC shrinkage (type: %s)...", shrink_type))
  res <- lfcShrink(dds, contrast = c("group", level_test, level_ref), 
                   res = res, type = shrink_type)
  
  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  res_df$ensembl_id <- sub("\\..*$", "", res_df$gene_id)
  res_df$gene_symbol <- NA_character_

  if (requireNamespace("AnnotationDbi", quietly = TRUE) &&
      requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    orgdb <- getExportedValue("org.Hs.eg.db", "org.Hs.eg.db")
    keys <- unique(res_df$ensembl_id)
    key_type <- if (all(grepl("^ENST", keys))) "ENSEMBLTRANS" else "ENSEMBL"
    columns <- if (key_type == "ENSEMBLTRANS") c("SYMBOL", "ENSEMBL") else c("SYMBOL")
    map_tbl <- suppressMessages(AnnotationDbi::select(
      orgdb,
      keys = keys,
      keytype = key_type,
      columns = columns
    ))
    if (!is.null(map_tbl) && nrow(map_tbl) > 0) {
      map_tbl <- map_tbl[!duplicated(map_tbl[[key_type]]), ]
      symbol_map <- setNames(map_tbl$SYMBOL, map_tbl[[key_type]])
      res_df$gene_symbol <- unname(symbol_map[res_df$ensembl_id])
      if ("ENSEMBL" %in% colnames(map_tbl) && !"parent_gene_id" %in% colnames(res_df)) {
        gene_map <- setNames(map_tbl$ENSEMBL, map_tbl[[key_type]])
        res_df$parent_gene_id <- unname(gene_map[res_df$ensembl_id])
      }
    }
  } else {
    warning("org.Hs.eg.db not available; gene_symbol column will be NA")
  }

  res_df <- dplyr::select(
    res_df,
    gene_id,
    ensembl_id,
    dplyr::any_of("parent_gene_id"),
    gene_symbol,
    baseMean,
    log2FoldChange,
    lfcSE,
    dplyr::any_of("stat"),
    pvalue,
    padj
  )
  res_df <- res_df[order(res_df$padj, na.last = NA), ]
  write.csv(res_df, file.path(comparison_dir, "differential_expression.csv"), row.names = FALSE, quote = FALSE)

  sig_df <- res_df %>%
    filter(!is.na(padj), padj <= 0.05, !is.na(log2FoldChange), abs(log2FoldChange) >= 1)

  upregulated_df <- sig_df %>%
    filter(log2FoldChange >= 1)
  write.csv(upregulated_df, file.path(comparison_dir, "significant_upregulated_genes.csv"), row.names = FALSE, quote = FALSE)

  pdf(file.path(comparison_dir, "ma_plot.pdf"), width = 8, height = 6)
  plotMA(res, main = paste("MA Plot:", comparison_name), ylim = c(-5, 5))
  dev.off()

  generate_volcano(res_df, comparison_name, file.path(comparison_dir, "volcano_plot.pdf"))

  top_genes <- head(res_df$gene_id, n = 50)
  top_genes <- top_genes[!is.na(top_genes)]
  top_genes <- intersect(top_genes, rownames(assay(vsd)))
  if (length(top_genes) > 2) {
    heatmap_mat <- assay(vsd)[top_genes, , drop = FALSE]
    heatmap_mat <- heatmap_mat - rowMeans(heatmap_mat)
    annot <- data.frame(Group = group_factor)
    rownames(annot) <- names(group_factor)
    annot <- annot[colnames(heatmap_mat), , drop = FALSE]
    group_cols <- group_color_map(levels(group_factor))
    gene_symbol_map <- setNames(res_df$gene_symbol, res_df$gene_id)
    rownames(heatmap_mat) <- ifelse(!is.na(gene_symbol_map[rownames(heatmap_mat)]),
                                    gene_symbol_map[rownames(heatmap_mat)],
                                    rownames(heatmap_mat))
    pdf(file.path(comparison_dir, "heatmap_top_de_genes.pdf"), width = 10, height = 8)
    pheatmap(
      heatmap_mat,
      annotation_col = annot,
      show_rownames = length(top_genes) <= 50,
      clustering_distance_rows = "correlation",
      clustering_distance_cols = "correlation",
      color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(255),
      annotation_colors = list(Group = group_cols),
      border_color = NA,
      main = paste("Top DE Genes:", comparison_name)
    )
    dev.off()
  }

  list(
    summary = data.frame(
      comparison = comparison_name,
      total_DE_genes = nrow(sig_df),
      upregulated = sum(sig_df$log2FoldChange > 0),
      downregulated = sum(sig_df$log2FoldChange < 0),
      significant_upregulated = nrow(upregulated_df),
      stringsAsFactors = FALSE
    )
  )
}

baseline_level_for_alias <- function(alias) {
  if (alias %in% c("fibrosis", "fibrosis_stage", "fibrosis_trim", "fibrosis_any", "fibrosis_binary")) {
    return("F0")
  }
  if (alias %in% c("steatosis", "steatosis_grade", "steatosis_any", "steatosis_binary")) {
    return("S0")
  }
  if (alias %in% c("stage")) {
    return("control")
  }
  if (alias %in% c("nas_high", "nas_low")) {
    return("NAS_0")
  }
  if (alias %in% c("masld", "nas_any")) {
    return("Control")
  }
  if (alias %in% c("ballooning", "ballooning_grade", "cytological_ballooning_grade", "ballooning_any")) {
    return("B0")
  }
  if (alias %in% c("inflammation_any")) {
    return("L0")
  }
  if (alias %in% c("fibrosis_strict")) {
    return("Healthy")
  }
  NULL
}

# Split group specifications -------------------------------------------------
group_entries <- strsplit(group_specs, ",")[[1]]
group_entries <- trimws(group_entries)
group_entries <- group_entries[group_entries != ""]

if (length(group_entries) == 0) {
  stop("No group specifications provided")
}

analysis_summaries <- list()

for (entry in group_entries) {
  parts <- strsplit(entry, "=")[[1]]
  column_name <- trimws(parts[1])
  alias <- if (length(parts) >= 2) trimws(parts[2]) else column_name
  column_name <- tolower(column_name)
  alias <- tolower(alias)

  if (!column_name %in% colnames(metadata)) {
    warning("Skipping grouping '", entry, "' because column '", column_name, "' is not in metadata")
    next
  }

  message("\n=== Processing grouping: ", column_name, " (alias: ", alias, ") ===")
  group_factor_all <- create_group_factor(metadata[[column_name]], alias)
  if (is.null(group_factor_all)) {
    warning("Grouping '", column_name, "' has only missing values; skipping")
    next
  }
  names(group_factor_all) <- rownames(metadata)

  valid_levels <- table(group_factor_all)
  eligible_levels <- names(valid_levels[valid_levels >= min_group_size])
  if (length(eligible_levels) < 2) {
    warning("Grouping '", column_name, "' has fewer than two levels with at least ", min_group_size, " samples; skipping")
    next
  }

  keep_samples <- names(group_factor_all) %in% rownames(metadata) &
    group_factor_all %in% eligible_levels
  selected_samples <- names(group_factor_all)[keep_samples]

  counts_sub <- count_matrix_filtered[, selected_samples, drop = FALSE]
  group_factor <- droplevels(group_factor_all[selected_samples])
  coldata <- DataFrame(group = group_factor, row.names = selected_samples)
  covariate_terms <- character()

  if (length(covariate_list) > 0) {
    for (covar in covariate_list) {
      covar_col <- tolower(make.names(covar))
      if (!covar_col %in% colnames(metadata)) {
        warning("Covariate '", covar, "' not found in metadata; skipping")
        next
      }
      vals <- metadata[selected_samples, covar_col, drop = TRUE]
      # Skip covariates with no variation
      if (all(is.na(vals)) || length(unique(na.omit(vals))) < 2) {
        warning("Covariate '", covar, "' has no variation; skipping")
        next
      }
      # Treat character/factor as factor; numeric stays numeric
      if (is.numeric(vals)) {
        coldata[[covar_col]] <- as.numeric(vals)
      } else {
        coldata[[covar_col]] <- factor(vals)
      }
      covariate_terms <- c(covariate_terms, covar_col)
    }
  }

  if (length(covariate_terms) > 0) {
    complete_idx <- complete.cases(as.data.frame(coldata[, covariate_terms, drop = FALSE]))
    if (!all(complete_idx)) {
      coldata <- coldata[complete_idx, , drop = FALSE]
      counts_sub <- counts_sub[, complete_idx, drop = FALSE]
      group_factor <- droplevels(group_factor[complete_idx])
    }
  }

  if (ncol(counts_sub) < 2 || nlevels(group_factor) < 2) {
    warning("Grouping '", alias, "' has fewer than two samples/levels after covariate filtering; skipping")
    next
  }

  design_terms <- c(covariate_terms, "group")
  design_formula <- as.formula(paste("~", paste(design_terms, collapse = " + ")))

  if (length(covariate_terms) > 0) {
    mm <- model.matrix(design_formula, data = coldata)
    mm_rank <- qr(mm)$rank
    while (mm_rank < ncol(mm) && length(covariate_terms) > 0) {
      warning("Design not full rank; dropping covariate ", tail(covariate_terms, 1))
      covariate_terms <- head(covariate_terms, -1)
      design_terms <- c(covariate_terms, "group")
      design_formula <- as.formula(paste("~", paste(design_terms, collapse = " + ")))
      mm <- model.matrix(design_formula, data = coldata)
      mm_rank <- qr(mm)$rank
    }
    if (mm_rank < ncol(mm)) {
      warning("Design still not full rank; proceeding with group only")
      design_terms <- "group"
      design_formula <- as.formula("~ group")
    }
  }

  dds <- DESeqDataSetFromMatrix(
    countData = counts_sub,
    colData = coldata,
    design = design_formula
  )
  dds <- DESeq(dds)
  vsd <- vst(dds, blind = FALSE)

  alias_dir <- file.path(output_dir, alias)
  dir.create(alias_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(alias_dir, "deseq2_results"), recursive = TRUE, showWarnings = FALSE)

  # Save normalized counts for this subset
  normalized_counts_subset <- counts(dds, normalized = TRUE)
  write.csv(normalized_counts_subset, file.path(alias_dir, "deseq2_normalized_counts.csv"), quote = FALSE)

  # PCA plots (2D and optional 3D)
  pca_res <- compute_pca_scores(vsd)
  if (!is.null(pca_res)) {
    pca_df <- pca_res$scores
    pca_df$group <- group_factor[match(pca_df$sample_id, names(group_factor))]
    pca_df <- pca_df[!is.na(pca_df$group), , drop = FALSE]
    percent_var <- round(pca_res$percent_var, 1)

    if (all(c("PC1", "PC2") %in% colnames(pca_df))) {
      pv1 <- if (length(percent_var) >= 1) percent_var[1] else NA_real_
      pv2 <- if (length(percent_var) >= 2) percent_var[2] else NA_real_
      pdf(file.path(alias_dir, "deseq2_pca_plot.pdf"), width = 8, height = 6)
      print(
        ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
          geom_point(size = 3, alpha = 0.8) +
          theme_bw(base_size = 14) +
          theme(
            plot.title = element_text(size = 16, face = "bold"),
            axis.title = element_text(size = 14),
            axis.text = element_text(size = 12),
            legend.title = element_text(size = 12),
            legend.text = element_text(size = 11)
          ) +
          labs(
            title = paste("PCA -", alias),
            x = paste0("PC1: ", pv1, "%"),
            y = paste0("PC2: ", pv2, "%")
          )
      )
      dev.off()
    }

    if (enable_3d_pca && alias %in% c("nas_high", "nas_low", "fibrosis", "fibrosis_any")) {
      save_3d_pca_plot(
        pca_scores = pca_df,
        percent_var = percent_var,
        group_factor = group_factor,
        title = paste("3D PCA -", alias),
        out_file = file.path(alias_dir, "deseq2_pca_plot_3d.html"),
        pdf_file = file.path(alias_dir, "deseq2_pca_plot_3d.pdf")
      )
    }
  }

  # Heatmap of top variable genes
  plot_heatmap(
    vsd = vsd,
    alias = alias,
    group_factor = group_factor,
    out_file = file.path(alias_dir, "deseq2_sample_heatmap.pdf")
  )

  # Pairwise comparisons
  group_levels <- levels(group_factor)
  comparison_stats <- list()
  if (length(group_levels) < 2) {
    warning("Grouping '", alias, "' has fewer than two levels after filtering; skipping comparisons")
    next
  }

  for (i in 1:(length(group_levels) - 1)) {
    for (j in (i + 1):length(group_levels)) {
      level1 <- group_levels[i]
      level2 <- group_levels[j]
      result <- generate_comparison_outputs(
        dds = dds,
        vsd = vsd,
        group_factor = group_factor,
        level_test = level2,
        level_ref = level1,
        base_dir = file.path(alias_dir, "deseq2_results")
      )
      if (!is.null(result$summary)) {
        comparison_stats[[length(comparison_stats) + 1]] <- result$summary
      }
    }
  }

  baseline_level <- baseline_level_for_alias(alias)
  if (!is.null(baseline_level) && baseline_level %in% group_levels) {
    baseline_dir <- file.path(alias_dir, "deseq2_results_control_vs")
    baseline_stats <- list()
    for (target_level in setdiff(group_levels, baseline_level)) {
      baseline_result <- generate_comparison_outputs(
        dds = dds,
        vsd = vsd,
        group_factor = group_factor,
        level_test = target_level,
        level_ref = baseline_level,
        base_dir = baseline_dir
      )
      if (!is.null(baseline_result$summary)) {
        baseline_stats[[length(baseline_stats) + 1]] <- baseline_result$summary
      }
    }
    if (length(baseline_stats) > 0) {
      baseline_summary <- bind_rows(baseline_stats)
      write.csv(baseline_summary, file.path(alias_dir, "deseq2_control_vs_summary.csv"), row.names = FALSE, quote = FALSE)
      analysis_summaries[[paste0(alias, "_control_vs")]] <- baseline_summary
    }
  }

  if (length(comparison_stats) > 0) {
    comparison_summary <- bind_rows(comparison_stats)
    write.csv(comparison_summary, file.path(alias_dir, "deseq2_comparison_summary.csv"), row.names = FALSE, quote = FALSE)
    analysis_summaries[[alias]] <- comparison_summary
  } else {
    warning("No comparisons were generated for grouping '", alias, "'")
  }
}

# Cross-analysis PCA comparisons on a shared PCA space -----------------------
comparison_dir <- file.path(output_dir, "comparison_plots", "pca")
dir.create(comparison_dir, recursive = TRUE, showWarnings = FALSE)

global_vsd <- vst(dds_base, blind = TRUE)
global_pca_res <- compute_pca_scores(global_vsd)

comparison_aliases <- c("nas_high", "nas_low", "masld", "fibrosis")
comparison_columns <- list(
  nas_high = c("nas_score", "nafld_activity_score"),
  nas_low = c("nas_score", "nafld_activity_score"),
  masld = c("nas_score", "nafld_activity_score"),
  fibrosis = c("fibrosis_stage")
)

if (!is.null(global_pca_res)) {
  global_pca_df <- global_pca_res$scores
  write.csv(global_pca_df, file.path(comparison_dir, "global_pca_scores.csv"), row.names = FALSE, quote = FALSE)
  percent_var_global <- round(global_pca_res$percent_var, 1)

  comparison_plot_data <- list()
  separation_results <- list()

  for (alias in comparison_aliases) {
    candidates <- comparison_columns[[alias]]
    column_name <- NULL
    for (cand in candidates) {
      if (cand %in% colnames(metadata)) {
        column_name <- cand
        break
      }
    }

    if (is.null(column_name)) {
      next
    }

    group_factor_all <- create_group_factor(metadata[[column_name]], alias)
    if (is.null(group_factor_all)) {
      next
    }
    names(group_factor_all) <- rownames(metadata)

    valid_levels <- table(group_factor_all)
    eligible_levels <- names(valid_levels[valid_levels >= min_group_size])
    if (length(eligible_levels) < 2) {
      next
    }

    keep_samples <- names(group_factor_all) %in% global_pca_df$sample_id & group_factor_all %in% eligible_levels
    sample_ids_alias <- names(group_factor_all)[keep_samples]
    if (length(sample_ids_alias) < 3) {
      next
    }

    plot_df <- global_pca_df[match(sample_ids_alias, global_pca_df$sample_id), ]
    plot_df$group <- droplevels(group_factor_all[sample_ids_alias])
    plot_df$alias <- alias
    comparison_plot_data[[length(comparison_plot_data) + 1]] <- plot_df

    pv1 <- if (length(percent_var_global) >= 1) percent_var_global[1] else NA_real_
    pv2 <- if (length(percent_var_global) >= 2) percent_var_global[2] else NA_real_

    pdf(file.path(comparison_dir, paste0("global_pca_pc1_pc2_", alias, ".pdf")), width = 8, height = 6)
    print(
      ggplot(plot_df, aes(x = PC1, y = PC2, color = group)) +
        geom_point(size = 2.5, alpha = 0.85) +
        theme_bw(base_size = 14) +
        theme(
          plot.title = element_text(size = 16, face = "bold"),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 11)
        ) +
        labs(
          title = paste("Global PCA (shared axes) -", alias),
          x = paste0("PC1: ", pv1, "%"),
          y = paste0("PC2: ", pv2, "%")
        )
    )
    dev.off()

    if (enable_3d_pca && alias %in% c("nas_high", "nas_low", "fibrosis")) {
      save_3d_pca_plot(
        pca_scores = plot_df,
        percent_var = percent_var_global,
        group_factor = group_factor_all,
        title = paste("Global 3D PCA -", alias),
        out_file = file.path(comparison_dir, paste0("global_pca_", alias, "_3d.html")),
        pdf_file = file.path(comparison_dir, paste0("global_pca_", alias, "_3d.pdf"))
      )
    }

    sep <- compute_separation_metrics(
      pca_scores = plot_df,
      group_factor = group_factor_all,
      alias = alias
    )
    if (!is.null(sep)) {
      separation_results[[length(separation_results) + 1]] <- sep
    }
  }

  if (length(comparison_plot_data) > 0) {
    combined_df <- bind_rows(comparison_plot_data)
    pdf(file.path(comparison_dir, "global_pca_pc1_pc2_facets.pdf"), width = 10, height = 8)
    print(
      ggplot(combined_df, aes(x = PC1, y = PC2, color = group)) +
        geom_point(size = 2, alpha = 0.8) +
        theme_bw(base_size = 14) +
        theme(
          strip.text = element_text(size = 14),
          plot.title = element_text(size = 16, face = "bold"),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 11)
        ) +
        facet_wrap(~alias, scales = "free") +
        labs(
          title = "Global PCA (shared axes) by analysis label",
          x = paste0("PC1: ", percent_var_global[1], "%"),
          y = paste0("PC2: ", percent_var_global[2], "%")
        )
    )
    dev.off()
  }

  if (length(separation_results) > 0) {
    separation_df <- bind_rows(separation_results)
    write.csv(separation_df, file.path(comparison_dir, "pca_separation_metrics.csv"), row.names = FALSE, quote = FALSE)
  }
}

# Save session info
writeLines(capture.output(sessionInfo()), file.path(output_dir, "deseq2_session_info.txt"))

message("\nDESeq2 analysis completed. Results written to: ", output_dir)
