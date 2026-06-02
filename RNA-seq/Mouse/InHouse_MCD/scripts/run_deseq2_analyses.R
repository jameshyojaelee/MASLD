#!/usr/bin/env Rscript
# Run multiple DESeq2 analyses for mouse liver diet study.

suppressPackageStartupMessages({
  library(DESeq2)
  library(tximport)
  library(readr)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(pheatmap)
  library(ggrepel)
})

# Unified color scheme (peach & magenta, matching QC palette)
diet_colors <- c("Control" = "#F2A45E", "MCD" = "#C23B75")
sex_shapes  <- c("Female" = 16, "Male" = 17)
padj_cutoff <- 0.1
overwrite_existing <- tolower(Sys.getenv("DESEQ2_OVERWRITE", "false")) %in% c("1", "true", "yes")
run_salmon_analyses <- tolower(Sys.getenv("RUN_SALMON_ANALYSES", "false")) %in% c("1", "true", "yes")
# Optional: specify coefficient name for LRT directionality (e.g., "week_num.dietMCD").
lrt_direction_coef <- Sys.getenv("LRT_DIRECTION_COEF", "")

shrink_type <- if (requireNamespace("ashr", quietly = TRUE)) "ashr" else "normal"
if (shrink_type == "normal") {
  message("Package 'ashr' not found; falling back to normal shrinkage.")
}

root <- normalizePath("./")
metadata_path <- file.path(root, "metadata", "samples.tsv")
counts_path <- file.path(root, "counts", "featurecounts", "gene_counts.txt")
salmon_dir <- file.path(root, "quant", "salmon")
report_dir <- file.path(root, "reports")
if (!dir.exists(report_dir)) dir.create(report_dir, recursive = TRUE)

samples <- read_tsv(metadata_path, col_types = cols()) %>%
  mutate(
    diet = factor(diet, levels = c("Control", "MCD")),
    sex = factor(sex, levels = c("Female", "Male")),
    week = factor(week),
    week_num = as.numeric(week)
  )
message("Loaded metadata for ", nrow(samples), " samples")
samples_df <- as.data.frame(samples)
rownames(samples_df) <- samples_df$sample_id

# Parse featureCounts output (skip comment lines starting with #).
raw_counts <- read_tsv(counts_path, comment = "#", col_types = cols())
message("Loaded featureCounts table with ", nrow(raw_counts), " genes")

# Harmonise featureCounts column names with sample identifiers returned by STAR.
bam_suffix <- ".Aligned.sortedByCoord.out.bam"
bam_cols <- colnames(raw_counts)[grepl(paste0(bam_suffix, "$"), colnames(raw_counts))]
bam_sample_ids <- str_replace(basename(bam_cols), paste0("\\", bam_suffix, "$"), "")
valid_bam <- bam_sample_ids %in% samples$sample_id
rename_map <- setNames(bam_sample_ids[valid_bam], bam_cols[valid_bam])
if (length(rename_map) > 0) {
  colnames(raw_counts)[match(names(rename_map), colnames(raw_counts))] <- rename_map
}
message("Renamed ", length(rename_map), " featureCounts columns to sample IDs")
count_matrix <- raw_counts %>%
  select(Geneid, any_of(samples$sample_id)) %>%
  column_to_rownames(var = "Geneid") %>%
  as.matrix()
message("Constructed count matrix with dimensions ", paste(dim(count_matrix), collapse = " x "))

stopifnot(all(samples$sample_id %in% colnames(count_matrix)))
count_matrix <- count_matrix[, samples$sample_id]

# Load gene annotations from GTF
gtf_path <- file.path(root, "reference", "raw", "gencode.vM33.annotation.gtf")
tx_annot <- NULL
if (file.exists(gtf_path)) {
  message("Loading gene annotations from GTF...")
  gtf_data <- read_tsv(gtf_path, comment = "#", col_names = FALSE, col_types = cols())
  gene_annot <- gtf_data %>%
    filter(X3 == "gene") %>%
    mutate(
      gene_id = str_match(X9, 'gene_id "([^"]+)"')[, 2],
      gene_name = str_match(X9, 'gene_name "([^"]+)"')[, 2],
      gene_type = str_match(X9, 'gene_type "([^"]+)"')[, 2]
    ) %>%
    select(gene = gene_id, gene_name, gene_type) %>%
    distinct()
  message("Loaded annotations for ", nrow(gene_annot), " genes")
  tx2gene <- gtf_data %>%
    filter(X3 == "transcript") %>%
    mutate(
      gene_id = str_match(X9, 'gene_id "([^"]+)"')[, 2],
      transcript_id = str_match(X9, 'transcript_id "([^"]+)"')[, 2]
    ) %>%
    filter(!is.na(gene_id), !is.na(transcript_id)) %>%
    select(transcript_id, gene_id) %>%
    distinct()
  message("Constructed tx2gene table with ", nrow(tx2gene), " transcript-gene pairs")
  tx_annot <- gtf_data %>%
    filter(X3 == "transcript") %>%
    mutate(
      gene_id = str_match(X9, 'gene_id "([^"]+)"')[, 2],
      transcript_id = str_match(X9, 'transcript_id "([^"]+)"')[, 2],
      gene_name = str_match(X9, 'gene_name "([^"]+)"')[, 2],
      gene_type = str_match(X9, 'gene_type "([^"]+)"')[, 2],
      transcript_type = str_match(X9, 'transcript_type "([^"]+)"')[, 2]
    ) %>%
    filter(!is.na(transcript_id)) %>%
    transmute(
      gene = transcript_id,
      transcript_id = transcript_id,
      gene_id = gene_id,
      gene_name = gene_name,
      gene_type = gene_type,
      transcript_type = transcript_type
    ) %>%
    distinct()
  message("Prepared transcript annotation table with ", nrow(tx_annot), " entries")
} else {
  message("GTF file not found, gene names will not be added to volcano plots")
  gene_annot <- NULL
  tx2gene <- NULL
  tx_annot <- NULL
}

txi_salmon_gene <- NULL
txi_salmon_transcripts <- NULL
if (run_salmon_analyses) {
  salmon_quant_paths <- file.path(salmon_dir, samples$sample_id, "quant.sf")
  names(salmon_quant_paths) <- samples$sample_id
  missing_salmon <- salmon_quant_paths[!file.exists(salmon_quant_paths)]
  if (length(missing_salmon) > 0) {
    message(
      "Skipping tximport: missing quant.sf for samples ",
      paste(names(missing_salmon), collapse = ", ")
    )
  } else if (is.null(tx2gene) || nrow(tx2gene) == 0) {
    message("Skipping tximport: tx2gene mapping unavailable.")
  } else {
    message("Running tximport to summarise Salmon transcript quantifications...")
    txi_salmon_gene <- tximport(
      salmon_quant_paths,
      type = "salmon",
      tx2gene = tx2gene,
      countsFromAbundance = "lengthScaledTPM"
    )
    message("tximport completed; gene-level counts available for ",
            nrow(txi_salmon_gene$counts), " genes")

    message("Running tximport in transcript-output mode for transcript-level analyses...")
    txi_salmon_transcripts <- tximport(
      salmon_quant_paths,
      type = "salmon",
      txOut = TRUE,
      countsFromAbundance = "lengthScaledTPM"
    )
    message("Transcript-level tximport completed; ",
            nrow(txi_salmon_transcripts$counts), " transcripts quantified")
  }
} else {
  message("Skipping Salmon tximport/analyses (RUN_SALMON_ANALYSES=false).")
}

# Helper to save normalized counts.
write_normalized_counts <- function(dds, outfile, gene_annotations = NULL) {
  norm_counts <- counts(dds, normalized = TRUE) %>%
    as.data.frame() %>%
    rownames_to_column("gene")

  if (!is.null(gene_annotations)) {
    norm_counts <- norm_counts %>%
      left_join(gene_annotations, by = "gene") %>%
      select(gene, gene_name, gene_type, everything())
  }

  write_csv(norm_counts, outfile)
  message("Saved normalized counts to ", outfile)
  invisible(norm_counts)
}

# Helper to save DESeq2 results table.
write_results <- function(res,
                          outfile,
                          gene_annotations = NULL,
                          upregulated_csv = NULL,
                          padj_threshold = padj_cutoff,
                          log2fc_threshold = 1,
                          lfc_column = "log2FoldChange") {
  res_tbl <- as.data.frame(res) %>%
    rownames_to_column("gene")

  if (!is.null(gene_annotations)) {
    res_tbl <- res_tbl %>%
      left_join(gene_annotations, by = "gene")
  }

  res_tbl <- res_tbl %>%
    arrange(padj)

  write_tsv(res_tbl, outfile)

  if (!is.null(upregulated_csv)) {
    if (!lfc_column %in% colnames(res_tbl)) {
      warning("Requested LFC column '", lfc_column, "' not found in results; skipping upregulated CSV.")
      return(invisible(res_tbl))
    }
    up_tbl <- res_tbl %>%
      filter(!is.na(padj), padj < padj_threshold,
             !is.na(.data[[lfc_column]]), .data[[lfc_column]] >= log2fc_threshold)
    write_csv(up_tbl, upregulated_csv)
  }

  invisible(res_tbl)
}

# Choose a coefficient to infer LRT directionality.
choose_lrt_coef <- function(dds, preferred = lrt_direction_coef) {
  available <- resultsNames(dds)
  if (!is.null(preferred) && nzchar(preferred) && preferred %in% available) {
    return(preferred)
  }
  # Prefer diet:week interaction terms when available.
  interaction_hits <- available[grepl("diet.*week|week.*diet", available)]
  if (length(interaction_hits) > 0) {
    return(interaction_hits[1])
  }
  # Fall back to main diet effect if present.
  diet_hits <- available[grepl("^diet|diet", available)]
  if (length(diet_hits) > 0) {
    return(diet_hits[1])
  }
  NA_character_
}

# Helper to subset tximport object to selected samples
subset_tximport <- function(txi_obj, sample_ids) {
  if (is.null(txi_obj)) {
    stop("txi_obj is NULL; cannot subset.")
  }
  missing_ids <- setdiff(sample_ids, colnames(txi_obj$counts))
  if (length(missing_ids) > 0) {
    stop("tximport object missing samples: ", paste(missing_ids, collapse = ", "))
  }
  txi_subset <- txi_obj
  matrix_slots <- intersect(
    names(txi_subset),
    c("counts", "abundance", "length", "variance", "varianceMat")
  )
  for (slot in matrix_slots) {
    if (!is.null(txi_subset[[slot]])) {
      txi_subset[[slot]] <- txi_subset[[slot]][, sample_ids, drop = FALSE]
    }
  }
  if (!is.null(txi_subset$infRep)) {
    txi_subset$infRep <- txi_subset$infRep[sample_ids]
  }
  txi_subset
}

# Helper to save MA plot and volcano plot with larger fonts and gene labels.
plot_contrast <- function(res, prefix, gene_annotations = NULL) {
  res_df <- as.data.frame(res) %>%
    rownames_to_column("gene")
  res_df$padj[is.na(res_df$padj)] <- 1
  res_df$log2FoldChange[is.na(res_df$log2FoldChange)] <- 0

  # Add gene names if annotations provided
  if (!is.null(gene_annotations)) {
    res_df <- res_df %>%
      left_join(gene_annotations %>% select(gene, gene_name), by = "gene")
  } else {
    res_df$gene_name <- res_df$gene
  }

  res_df <- res_df %>%
    mutate(
      signif = padj < 0.1,
      neg_log10_padj = -log10(padj),
      significance = case_when(
        padj < 0.1 & log2FoldChange > 1 ~ "Up (padj<0.1, log2FC>1)",
        padj < 0.1 & log2FoldChange < -1 ~ "Down (padj<0.1, log2FC<-1)",
        padj < 0.1 ~ "Significant (padj<0.1)",
        TRUE ~ "Not significant"
      )
    )

  # MA plot with larger fonts
  ma_plot <- ggplot(res_df, aes(baseMean, log2FoldChange, colour = signif)) +
    geom_point(alpha = 0.5, size = 2.5) +
    scale_x_log10() +
    scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "#C23B75")) +
    theme_bw(base_size = 20) +
    labs(title = "MA Plot", colour = "padj < 0.1",
         x = "Mean Expression", y = "log2 Fold Change") +
    theme(
      plot.title = element_text(size = 26, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 24, face = "bold"),
      axis.text = element_text(size = 18),
      legend.title = element_text(size = 22, face = "bold"),
      legend.text = element_text(size = 18),
      legend.position = "right",
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank()
    )
  ggsave(filename = paste0(prefix, "_MA.pdf"), plot = ma_plot, width = 10, height = 8, dpi = 300)
  ggsave(filename = paste0(prefix, "_MA.png"), plot = ma_plot, width = 10, height = 8, dpi = 300)

  # Volcano plot with larger fonts and top gene labels
  top_genes <- res_df %>%
    filter(!is.na(padj), abs(log2FoldChange) > 2) %>%
    arrange(padj) %>%
    head(15)

  volcano <- ggplot(res_df, aes(log2FoldChange, neg_log10_padj)) +
    geom_point(aes(color = significance), alpha = 0.6, size = 2.5) +
    geom_point(data = top_genes, aes(x = log2FoldChange, y = neg_log10_padj),
               color = "black", size = 4, shape = 21, fill = "yellow", stroke = 1.5) +
    geom_text_repel(data = top_genes,
                    aes(label = gene_name),
                    size = 6.5,
                    fontface = "bold",
                    box.padding = 1.0,
                    point.padding = 0.7,
                    segment.size = 0.6,
                    segment.color = "gray30",
                    max.overlaps = 20,
                    min.segment.length = 0,
                    force = 2) +
    geom_hline(yintercept = -log10(0.1), linetype = "dashed", color = "gray30", linewidth = 1.0) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray30", linewidth = 1.0) +
    scale_color_manual(
      values = c(
        "Up (padj<0.1, log2FC>1)" = "#C23B75",   # magenta
        "Down (padj<0.1, log2FC<-1)" = "#F2A45E", # peach
        "Significant (padj<0.1)" = "#E07BB6",    # light magenta
        "Not significant" = "gray70"
      ),
      breaks = c("Up (padj<0.1, log2FC>1)", "Down (padj<0.1, log2FC<-1)",
                 "Significant (padj<0.1)", "Not significant")
    ) +
    theme_bw(base_size = 20) +
    labs(title = "Volcano Plot",
         x = "log2 Fold Change",
         y = "-log10(adjusted p-value)",
         color = "Significance") +
    theme(
      plot.title = element_text(size = 26, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 24, face = "bold"),
      axis.text = element_text(size = 18),
      legend.title = element_text(size = 22, face = "bold"),
      legend.text = element_text(size = 18),
      legend.position = "bottom",
      legend.key.size = unit(1.5, "lines"),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank()
    ) +
    guides(color = guide_legend(override.aes = list(size = 5, alpha = 1)))
  ggsave(filename = paste0(prefix, "_volcano.pdf"), plot = volcano, width = 14, height = 10, dpi = 300)
  ggsave(filename = paste0(prefix, "_volcano.png"), plot = volcano, width = 14, height = 10, dpi = 300)
}

# Variance stabilized PCA plot per dataset with MUCH LARGER fonts.
plot_pca <- function(dds, outfile, sample_info) {
  vsd <- vst(dds, blind = TRUE)
  pca <- plotPCA(vsd, intgroup = c("diet", "sex", "week"), returnData = TRUE)
  percentVar <- round(100 * attr(pca, "percentVar"), 1)

  # Add sample names for labeling
  pca$sample_name <- sample_info$sample_id[match(rownames(pca), sample_info$sample_id)]

  p <- ggplot(pca, aes(PC1, PC2, colour = diet, shape = sex)) +
    geom_point(size = 10, alpha = 0.8, stroke = 2.5) +
    geom_text_repel(aes(label = sample_name),
                    size = 8,
                    box.padding = 1.0,
                    point.padding = 0.7,
                    segment.size = 0.5,
                    max.overlaps = 20,
                    show.legend = FALSE) +
    scale_color_manual(values = diet_colors) +
    scale_shape_manual(values = sex_shapes) +
    theme_bw(base_size = 22) +
    labs(
      title = "PCA Plot",
      x = paste0("PC1 (", percentVar[1], "% variance)"),
      y = paste0("PC2 (", percentVar[2], "% variance)"),
      color = "Diet",
      shape = "Sex"
    ) +
    theme(
      plot.title = element_text(size = 30, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 26, face = "bold"),
      axis.text = element_text(size = 20),
      legend.title = element_text(size = 24, face = "bold"),
      legend.text = element_text(size = 20),
      legend.position = "right",
      legend.key.size = unit(2.5, "lines"),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank(),
      plot.margin = margin(20, 20, 20, 20)
    )

  ggsave(outfile, plot = p, width = 14, height = 10, dpi = 300)
  ggsave(sub("\\.pdf$", ".png", outfile), plot = p, width = 14, height = 10, dpi = 300)

  # Also save PCA coordinates
  pca_coords <- pca %>%
    mutate(
      group = paste(diet, sex, week, sep = ":"),
      name = sample_name
    ) %>%
    select(PC1, PC2, group, diet, sex, week, name)

  outfile_tsv <- sub("\\.pdf$", "_coordinates.tsv", outfile)
  write_tsv(pca_coords, outfile_tsv)
}

run_analysis <- function(subset_samples, design_formula, output_dir, contrast = NULL,
                         test = c("Wald", "LRT"), reduced = NULL, name = NULL,
                         gene_annotations = NULL, count_matrix_source = count_matrix,
                         txi = NULL, overwrite = overwrite_existing) {
  test <- match.arg(test)
  if (dir.exists(output_dir) && !overwrite) {
    existing_files <- list.files(output_dir, include.dirs = FALSE)
    if (length(existing_files) > 0) {
      message("Skipping analysis in ", output_dir,
              " (existing outputs detected; set overwrite = TRUE or DESEQ2_OVERWRITE=1 to regenerate).")
      return(invisible(NULL))
    }
  }
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  sample_info <- samples %>% filter(sample_id %in% subset_samples)
  sample_info <- sample_info[match(subset_samples, sample_info$sample_id), ]
  sample_info <- as.data.frame(sample_info)
  rownames(sample_info) <- sample_info$sample_id
  message("Running analysis in ", output_dir, " | design: ", deparse(design_formula),
          " | samples: ", length(subset_samples))
  if (!is.null(txi)) {
    subset_txi <- subset_tximport(txi, subset_samples)
    dds <- DESeqDataSetFromTximport(
      txi = subset_txi,
      colData = sample_info,
      design = design_formula
    )
  } else {
    if (is.null(count_matrix_source)) {
      stop("count_matrix_source is NULL while txi is not provided.")
    }
    missing_samples <- setdiff(subset_samples, colnames(count_matrix_source))
    if (length(missing_samples) > 0) {
      stop("Count matrix missing samples: ", paste(missing_samples, collapse = ", "))
    }
    dds <- DESeqDataSetFromMatrix(
      countData = count_matrix_source[, subset_samples, drop = FALSE],
      colData = sample_info,
      design = design_formula
    )
  }
  message("No raw-count prefilter applied for DESeq2")
  if (test == "LRT" && is.null(reduced)) {
    stop("Reduced formula is required for LRT")
  }
  if (test == "LRT") {
    dds <- DESeq(dds, test = test, reduced = reduced)
  } else {
    dds <- DESeq(dds, test = test)
  }
  message("DESeq() completed")
  plot_pca(dds, file.path(output_dir, "PCA.pdf"), sample_info)
  if (test == "Wald") {
    if (is.null(contrast) && is.null(name)) {
      stop("Wald test requires contrast or name")
    }
    if (!is.null(contrast)) {
      res <- results(dds, contrast = contrast)
      res <- lfcShrink(dds, contrast = contrast, res = res, type = shrink_type)
    } else {
      res <- results(dds, name = name)
      res <- lfcShrink(dds, coef = name, type = shrink_type)
    }
    write_results(
      res,
      file.path(output_dir, "deseq2_results.tsv"),
      gene_annotations = gene_annotations,
      upregulated_csv = file.path(output_dir, "deseq2_upregulated.csv")
    )
    plot_contrast(res, file.path(output_dir, "deseq2"), gene_annotations)
    message("Saved Wald results (including upregulated CSV) to ", output_dir)
  } else {
    res <- results(dds)
    direction_coef <- choose_lrt_coef(dds)
    if (!is.na(direction_coef) && nzchar(direction_coef)) {
      direction_res <- results(dds, name = direction_coef)
      res$direction_log2FoldChange <- direction_res$log2FoldChange
      res$direction_coef <- direction_coef
      message("LRT directionality inferred from coefficient: ", direction_coef)
    } else {
      message("LRT directionality coefficient not found; direction columns will be NA.")
    }
    lfc_column <- if ("direction_log2FoldChange" %in% colnames(res)) {
      "direction_log2FoldChange"
    } else {
      "log2FoldChange"
    }
    write_results(
      res,
      file.path(output_dir, "lrt_results.tsv"),
      gene_annotations = gene_annotations,
      upregulated_csv = file.path(output_dir, "lrt_upregulated.csv"),
      lfc_column = lfc_column
    )
    message("Saved LRT results (including upregulated CSV) to ", output_dir)
  }
  saveRDS(dds, file.path(output_dir, "dds.rds"))
  message("Saved dds object for ", output_dir)
}

# Version 1: Diet effect per sex and combined
run_analysis(
  subset_samples = samples$sample_id[samples$sex == "Female"],
  design_formula = ~ week + diet,
  output_dir = file.path(root, "analysis_v1_weekPooled_female"),
  contrast = c("diet", "MCD", "Control"),
  test = "Wald",
  gene_annotations = gene_annot,
  overwrite = overwrite_existing
)

run_analysis(
  subset_samples = samples$sample_id[samples$sex == "Male"],
  design_formula = ~ week + diet,
  output_dir = file.path(root, "analysis_v1_weekPooled_male"),
  contrast = c("diet", "MCD", "Control"),
  test = "Wald",
  gene_annotations = gene_annot,
  overwrite = overwrite_existing
)

run_analysis(
  subset_samples = samples$sample_id,
  design_formula = ~ sex + week + diet,
  output_dir = file.path(root, "analysis_v1_weekPooled_combined"),
  contrast = c("diet", "MCD", "Control"),
  test = "Wald",
  gene_annotations = gene_annot,
  overwrite = overwrite_existing
)

# Version 2: Time-course per sex (diet x week interaction)
run_analysis(
  subset_samples = samples$sample_id[samples$sex == "Female"],
  design_formula = ~ week_num + diet + week_num:diet,
  output_dir = file.path(root, "analysis_v2_timeCourse_female"),
  test = "LRT",
  reduced = ~ week_num + diet,
  overwrite = overwrite_existing
)

run_analysis(
  subset_samples = samples$sample_id[samples$sex == "Male"],
  design_formula = ~ week_num + diet + week_num:diet,
  output_dir = file.path(root, "analysis_v2_timeCourse_male"),
  test = "LRT",
  reduced = ~ week_num + diet,
  overwrite = overwrite_existing
)

# Export normalized counts for all samples
message("Exporting DESeq2 normalized counts for all samples...")
dds_all <- DESeqDataSetFromMatrix(
  countData = count_matrix[, samples$sample_id],
  colData = samples_df[samples$sample_id, , drop = FALSE],
  design = ~ sex + week + diet
)
message("No raw-count prefilter applied for normalized counts export")
dds_all <- DESeq(dds_all)
norm_counts_feature <- file.path(root, "normalized_counts_all_samples.csv")
if (!file.exists(norm_counts_feature)) {
  write_normalized_counts(
    dds_all,
    norm_counts_feature,
    gene_annotations = gene_annot
  )
  message("Normalized counts exported successfully")
} else {
  message("Skipping normalized counts export (file already exists): ", norm_counts_feature)
}

# Targeted comparisons: Week-specific MCD cohorts vs all controls
control_samples <- samples %>%
  filter(diet == "Control") %>%
  pull(sample_id)

for (wk in c("1", "2", "3")) {
  week_mcd_samples <- samples %>%
    filter(diet == "MCD", week == wk) %>%
    pull(sample_id)
  target_subset <- samples %>%
    filter(sample_id %in% c(week_mcd_samples, control_samples)) %>%
    pull(sample_id)
  output_dir <- file.path(root, sprintf("analysis_v5_week%sMCD_vs_controls", wk))

  run_analysis(
    subset_samples = target_subset,
    design_formula = ~ sex + diet,
    output_dir = output_dir,
    contrast = c("diet", "MCD", "Control"),
    test = "Wald",
    gene_annotations = gene_annot,
    overwrite = overwrite_existing
  )
}

if (!is.null(txi_salmon_gene)) {
  message("Running Salmon-based week-specific MCD vs Control analyses...")
  for (wk in c("1", "2", "3")) {
    week_mcd_samples <- samples %>%
      filter(diet == "MCD", week == wk) %>%
      pull(sample_id)
    target_subset <- samples %>%
      filter(sample_id %in% c(week_mcd_samples, control_samples)) %>%
      pull(sample_id)
    output_dir <- file.path(root, sprintf("analysis_v5_salmon_week%sMCD_vs_controls", wk))

    run_analysis(
      subset_samples = target_subset,
      design_formula = ~ sex + diet,
      output_dir = output_dir,
      contrast = c("diet", "MCD", "Control"),
      test = "Wald",
      gene_annotations = gene_annot,
      txi = txi_salmon_gene,
      overwrite = overwrite_existing
    )
  }

  salmon_norm_out <- file.path(root, "normalized_counts_all_samples_salmon.csv")
  if (!file.exists(salmon_norm_out)) {
    txi_all_samples <- subset_tximport(txi_salmon_gene, samples$sample_id)
    dds_all_salmon <- DESeqDataSetFromTximport(
      txi = txi_all_samples,
      colData = samples_df[samples$sample_id, , drop = FALSE],
      design = ~ sex + week + diet
    )
    message("No raw-count prefilter applied for salmon normalized counts export")
    dds_all_salmon <- DESeq(dds_all_salmon)
    write_normalized_counts(
      dds_all_salmon,
      salmon_norm_out,
      gene_annotations = gene_annot
    )
    message("Salmon-based normalized counts exported successfully")
  } else {
    message("Skipping Salmon-based normalized counts export (file already exists): ",
            salmon_norm_out)
  }
} else if (run_salmon_analyses) {
  message("Salmon-based analyses were skipped because tximport results were unavailable.")
} else {
  message("Salmon-based analyses disabled (RUN_SALMON_ANALYSES=false).")
}

if (!is.null(txi_salmon_transcripts)) {
  message("Running Salmon transcript-level week-specific MCD vs Control analyses...")
  for (wk in c("1", "2", "3")) {
    week_mcd_samples <- samples %>%
      filter(diet == "MCD", week == wk) %>%
      pull(sample_id)
    target_subset <- samples %>%
      filter(sample_id %in% c(week_mcd_samples, control_samples)) %>%
      pull(sample_id)
    output_dir <- file.path(root, sprintf("analysis_v5_salmon_transcript_week%sMCD_vs_controls", wk))

    run_analysis(
      subset_samples = target_subset,
      design_formula = ~ sex + diet,
      output_dir = output_dir,
      contrast = c("diet", "MCD", "Control"),
      test = "Wald",
      gene_annotations = tx_annot,
      txi = txi_salmon_transcripts,
      overwrite = overwrite_existing
    )
  }

  salmon_tx_norm_out <- file.path(root, "normalized_counts_all_samples_salmon_transcripts.csv")
  if (!file.exists(salmon_tx_norm_out)) {
    txi_all_transcripts <- subset_tximport(txi_salmon_transcripts, samples$sample_id)
    dds_all_transcripts <- DESeqDataSetFromTximport(
      txi = txi_all_transcripts,
      colData = samples_df[samples$sample_id, , drop = FALSE],
      design = ~ sex + week + diet
    )
    message("No raw-count prefilter applied for transcript normalized counts export")
    dds_all_transcripts <- DESeq(dds_all_transcripts)
    write_normalized_counts(
      dds_all_transcripts,
      salmon_tx_norm_out,
      gene_annotations = tx_annot
    )
    message("Salmon transcript-level normalized counts exported successfully")
  } else {
    message("Skipping Salmon transcript-level normalized counts export (file already exists): ",
            salmon_tx_norm_out)
  }
} else if (run_salmon_analyses) {
  message("Salmon transcript-level analyses were skipped because tximport results were unavailable.")
} else {
  message("Salmon transcript-level analyses disabled (RUN_SALMON_ANALYSES=false).")
}

# Version 3: Time-course combined sexes with interaction term
run_analysis(
  subset_samples = samples$sample_id,
  design_formula = ~ sex + week_num + diet + week_num:diet,
  output_dir = file.path(root, "analysis_v3_timeCourse_combined"),
  test = "LRT",
  reduced = ~ sex + week_num + diet,
  overwrite = overwrite_existing
)

# Version 4: Sex differences within each diet (Wald)
run_analysis(
  subset_samples = samples$sample_id[samples$diet == "Control"],
  design_formula = ~ week + sex,
  output_dir = file.path(root, "analysis_v4_sexDifferences_control"),
  contrast = c("sex", "Male", "Female"),
  test = "Wald",
  gene_annotations = gene_annot,
  overwrite = overwrite_existing
)

run_analysis(
  subset_samples = samples$sample_id[samples$diet == "MCD"],
  design_formula = ~ week + sex,
  output_dir = file.path(root, "analysis_v4_sexDifferences_MCD"),
  contrast = c("sex", "Male", "Female"),
  test = "Wald",
  gene_annotations = gene_annot,
  overwrite = overwrite_existing
)

message("All DESeq2 analyses completed.")
