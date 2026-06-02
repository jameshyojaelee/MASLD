#!/usr/bin/env Rscript

# Meta-analysis across DESeq2 result sets (v1–v5) with extensive visualisations.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(scales)
  library(forcats)
})

root <- normalizePath("./")
output_dir <- file.path(root, "meta_analysis")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Harmonised metadata for each DESeq2 result set that we want to compare.
analysis_info <- tribble(
  ~analysis_id,        ~version, ~category,           ~analysis_label,                     ~path,                                                                                      ~test, ~description,
  "v1_female",         "v1",     "Diet effect",       "Female MCD vs Control",             file.path(root, "analysis_v1_weekPooled_female",     "deseq2_results.tsv"),                               "Wald", "Week-pooled female cohort; diet effect (MCD vs Control).",
  "v1_male",           "v1",     "Diet effect",       "Male MCD vs Control",               file.path(root, "analysis_v1_weekPooled_male",       "deseq2_results.tsv"),                               "Wald", "Week-pooled male cohort; diet effect (MCD vs Control).",
  "v1_combined",       "v1",     "Diet effect",       "Combined-sex MCD vs Control",       file.path(root, "analysis_v1_weekPooled_combined",   "deseq2_results.tsv"),                               "Wald", "Week-pooled sexes combined; diet effect (MCD vs Control).",
  "v2_female",         "v2",     "Diet x week LRT",   "Female diet + week LRT",            file.path(root, "analysis_v2_timeCourse_female",     "lrt_results_annotated.tsv"),                       "LRT",  "Female-only time-course; LRT for diet + week effects.",
  "v2_male",           "v2",     "Diet x week LRT",   "Male diet + week LRT",              file.path(root, "analysis_v2_timeCourse_male",       "lrt_results_annotated.tsv"),                       "LRT",  "Male-only time-course; LRT for diet + week effects.",
  "v3_combined",       "v3",     "Diet x week LRT",   "Combined diet-week interaction LRT",file.path(root, "analysis_v3_timeCourse_combined",   "lrt_results_annotated.tsv"),                       "LRT",  "Combined sexes; LRT for diet x week interaction.",
  "v4_control",        "v4",     "Sex difference",    "Control diet Male vs Female",       file.path(root, "analysis_v4_sexDifferences_control","deseq2_results.tsv"),                               "Wald", "Control-fed animals; male vs female contrast.",
  "v4_MCD",            "v4",     "Sex difference",    "MCD diet Male vs Female",           file.path(root, "analysis_v4_sexDifferences_MCD",    "deseq2_results.tsv"),                               "Wald", "MCD-fed animals; male vs female contrast.",
  "v5_week3",          "v5",     "Focused contrast",  "Week 3 MCD vs pooled Controls",     file.path(root, "analysis_v5_week3MCD_vs_controls",  "deseq2_results.tsv"),                               "Wald", "Week 3 MCD animals vs pooled controls (sex adjusted)."
)

missing_files <- analysis_info %>%
  filter(!file.exists(path))

if (nrow(missing_files) > 0) {
  stop("Missing expected result files:\n", paste(missing_files$path, collapse = "\n"))
}

analysis_order <- analysis_info$analysis_id
label_order <- analysis_info$analysis_label

# Helper to load each result table and tag with metadata.
read_analysis_table <- function(path, analysis_id, analysis_label, version, category, test, description) {
  df <- read_tsv(path, col_types = cols())

  required_cols <- c("gene", "log2FoldChange", "padj")
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop("Result file ", path, " is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  if (!"gene_name" %in% colnames(df)) {
    df$gene_name <- NA_character_
  }
  if (!"gene_type" %in% colnames(df)) {
    df$gene_type <- NA_character_
  }

  df %>%
    mutate(
      gene = as.character(gene),
      gene_name = as.character(gene_name),
      gene_type = as.character(gene_type),
      gene_label = if_else(!is.na(gene_name) & gene_name != "", gene_name, gene),
      gene_display = if_else(!is.na(gene_name) & gene_name != "", paste0(gene_name, " (", gene, ")"), gene),
      analysis_id = analysis_id,
      analysis_label = analysis_label,
      version = version,
      category = category,
      test = test,
      description = description
    ) %>%
    relocate(analysis_id, analysis_label, version, category, test, gene, gene_label, gene_display, gene_name, gene_type)
}

all_results <- pmap_dfr(
  analysis_info,
  ~ read_analysis_table(..5, ..1, ..4, ..2, ..3, ..6, ..7)
) %>%
  mutate(
    padj = as.numeric(padj),
    log2FoldChange = as.numeric(log2FoldChange),
    baseMean = if ("baseMean" %in% names(.)) as.numeric(baseMean) else NA_real_,
    stat = if ("stat" %in% names(.)) as.numeric(stat) else NA_real_,
    pvalue = if ("pvalue" %in% names(.)) as.numeric(pvalue) else NA_real_
  )

all_results <- all_results %>%
  mutate(
    significant = !is.na(padj) & padj <= 0.05,
    direction = case_when(
      significant & !is.na(log2FoldChange) & log2FoldChange > 0 ~ "Up",
      significant & !is.na(log2FoldChange) & log2FoldChange < 0 ~ "Down",
      TRUE ~ "NS"
    ),
    abs_log2fc = abs(log2FoldChange),
    neg_log10_p = if_else(!is.na(padj) & padj > 0, -log10(padj), NA_real_),
    analysis_label = factor(analysis_label, levels = label_order),
    analysis_id = factor(analysis_id, levels = analysis_order)
  )

# Utility to save plots as png + pdf siblings.
save_plot_dual <- function(plot, filename_base, width = 10, height = 6) {
  ggsave(paste0(filename_base, ".png"), plot = plot, width = width, height = height, dpi = 300)
  ggsave(paste0(filename_base, ".pdf"), plot = plot, width = width, height = height)
}

# Helper to make filesystem-friendly identifiers.
slugify <- function(x) {
  x %>%
    stringr::str_to_lower() %>%
    stringr::str_replace_all("[^a-z0-9]+", "_") %>%
    stringr::str_replace_all("_+", "_") %>%
    stringr::str_replace_all("^_|_$", "")
}

safe_max <- function(x) {
  if (length(x) == 0 || all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
}

safe_min <- function(x) {
  if (length(x) == 0 || all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
}

first_non_missing <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) NA_character_ else x[1]
}

# Summary counts per analysis.
counts_summary <- all_results %>%
  group_by(analysis_id, analysis_label, version, category, test, description) %>%
  summarise(
    n_up = sum(direction == "Up", na.rm = TRUE),
    n_down = sum(direction == "Down", na.rm = TRUE),
    n_sig = n_up + n_down,
    median_log2fc_up = if_else(n_up > 0, median(log2FoldChange[direction == "Up"], na.rm = TRUE), NA_real_),
    median_log2fc_down = if_else(n_down > 0, median(log2FoldChange[direction == "Down"], na.rm = TRUE), NA_real_),
    .groups = "drop"
  ) %>%
  arrange(match(analysis_id, analysis_order))

write_csv(counts_summary, file.path(output_dir, "analysis_significant_counts.csv"))

counts_long <- counts_summary %>%
  select(analysis_label, n_up, n_down) %>%
  pivot_longer(cols = c(n_up, n_down), names_to = "direction", values_to = "count") %>%
  mutate(
    direction = factor(direction, levels = c("n_up", "n_down"), labels = c("Up", "Down")),
    analysis_label = factor(analysis_label, levels = label_order)
  )

palette_direction <- c("Up" = "#C23B75", "Down" = "#3075B1")

p_counts <- ggplot(counts_long, aes(x = analysis_label, y = count, fill = direction)) +
  geom_col(position = "stack", width = 0.75) +
  coord_flip() +
  scale_fill_manual(values = palette_direction) +
  theme_minimal(base_size = 14) +
  labs(
    x = NULL,
    y = "Number of genes (padj <= 0.05)",
    fill = "Direction",
    title = "Significant genes per analysis",
    subtitle = "Padj <= 0.05; direction determined by shrunken log2 fold-change sign"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.position = "top",
    plot.margin = margin(18, 26, 18, 18)
  )

save_plot_dual(p_counts, file.path(output_dir, "significant_counts_bar"), width = 10, height = 6)

# Build helper for pairwise set overlaps (Jaccard indices, intersections, etc).
extract_sets <- function(direction_value, analysis_subset) {
  analysis_subset <- as.character(analysis_subset)
  base <- tibble(analysis_id = analysis_subset)
  observed <- all_results %>%
    filter(direction == direction_value, analysis_id %in% analysis_subset) %>%
    group_by(analysis_id) %>%
    summarise(
      genes = list(unique(gene)),
      count = n_distinct(gene),
      .groups = "drop"
    )

  base %>%
    left_join(observed, by = "analysis_id") %>%
    mutate(
      genes = map(genes, ~ if (is.null(.x)) character(0) else .x),
      count = replace_na(count, 0)
    )
}

compute_pairwise_overlap <- function(set_tbl, label_levels) {
  ids <- as.character(set_tbl$analysis_id)
  gene_map <- setNames(set_tbl$genes, ids)
  expand_grid(analysis_a = ids, analysis_b = ids) %>%
    mutate(
      genes_a = map(analysis_a, ~ gene_map[[as.character(.x)]]),
      genes_b = map(analysis_b, ~ gene_map[[as.character(.x)]]),
      intersection = map2_int(genes_a, genes_b, ~ length(intersect(.x, .y))),
      union = map2_int(genes_a, genes_b, ~ length(union(.x, .y))),
      jaccard = if_else(union == 0, 0, intersection / union)
    ) %>%
    select(-genes_a, -genes_b) %>%
    left_join(analysis_info %>% select(analysis_a = analysis_id, analysis_a_label = analysis_label),
              by = "analysis_a") %>%
    left_join(analysis_info %>% select(analysis_b = analysis_id, analysis_b_label = analysis_label),
              by = "analysis_b") %>%
    mutate(
      analysis_a_label = factor(analysis_a_label, levels = label_levels),
      analysis_b_label = factor(analysis_b_label, levels = label_levels)
    )
}

plot_overlap_heatmap <- function(overlap_tbl, title_text, filename_base) {
  heat <- ggplot(overlap_tbl, aes(x = analysis_b_label, y = analysis_a_label, fill = jaccard)) +
    geom_tile(color = "white", size = 0.3) +
    geom_text(aes(label = scales::percent(jaccard, accuracy = 1)), size = 3, color = "black") +
    scale_fill_gradient2(low = "#F2F0F7", mid = "#9E9AC8", high = "#54278F", midpoint = 0.35, limits = c(0, 1),
                         labels = percent_format(accuracy = 1), name = "Jaccard") +
    theme_minimal(base_size = 13) +
    labs(
      x = NULL,
      y = NULL,
      title = title_text,
      subtitle = "Jaccard index of gene sets (intersection / union)"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 11),
      legend.position = "top",
      panel.grid = element_blank(),
      plot.margin = margin(20, 24, 20, 20)
    )
  save_plot_dual(heat, file.path(output_dir, filename_base), width = 9, height = 7)
}

# Diet-focused subset (MCD vs Control contrasts only).
diet_ids <- analysis_info %>%
  filter(category %in% c("Diet effect", "Focused contrast")) %>%
  pull(analysis_id)
diet_ids <- analysis_order[analysis_order %in% diet_ids]

if (length(diet_ids) >= 2) {
  diet_labels <- analysis_info %>%
    filter(analysis_id %in% diet_ids) %>%
    arrange(match(analysis_id, diet_ids)) %>%
    pull(analysis_label)

  up_overlap_diet <- extract_sets("Up", diet_ids) %>%
    compute_pairwise_overlap(label_levels = diet_labels)
  down_overlap_diet <- extract_sets("Down", diet_ids) %>%
    compute_pairwise_overlap(label_levels = diet_labels)

  write_csv(up_overlap_diet, file.path(output_dir, "pairwise_jaccard_up_diet.csv"))
  write_csv(down_overlap_diet, file.path(output_dir, "pairwise_jaccard_down_diet.csv"))

  plot_overlap_heatmap(up_overlap_diet, "Upregulated overlap (diet contrasts only)", "upregulated_jaccard_heatmap_diet_only")
  plot_overlap_heatmap(down_overlap_diet, "Downregulated overlap (diet contrasts only)", "downregulated_jaccard_heatmap_diet_only")
}

# Sex-specific subset (Male vs Female within diets).
sex_ids <- c("v4_control", "v4_MCD")
sex_ids <- sex_ids[sex_ids %in% analysis_order]

if (length(sex_ids) >= 2) {
  sex_labels <- analysis_info %>%
    filter(analysis_id %in% sex_ids) %>%
    arrange(match(analysis_id, sex_ids)) %>%
    pull(analysis_label)

  up_overlap_sex <- extract_sets("Up", sex_ids) %>%
    compute_pairwise_overlap(label_levels = sex_labels)
  down_overlap_sex <- extract_sets("Down", sex_ids) %>%
    compute_pairwise_overlap(label_levels = sex_labels)

  write_csv(up_overlap_sex, file.path(output_dir, "pairwise_jaccard_up_sex.csv"))
  write_csv(down_overlap_sex, file.path(output_dir, "pairwise_jaccard_down_sex.csv"))

  plot_overlap_heatmap(up_overlap_sex, "Upregulated overlap (sex contrasts only)", "upregulated_jaccard_heatmap_sex_only")
  plot_overlap_heatmap(down_overlap_sex, "Downregulated overlap (sex contrasts only)", "downregulated_jaccard_heatmap_sex_only")
}

compute_logfc_correlation <- function(analysis_subset, label_levels, filename_base, title_text, csv_name) {
  sig_gene_ids <- all_results %>%
    filter(direction != "NS", analysis_id %in% analysis_subset) %>%
    pull(gene) %>%
    unique()

  logfc_matrix <- all_results %>%
    filter(gene %in% sig_gene_ids, analysis_id %in% analysis_subset) %>%
    select(gene, analysis_id, log2FoldChange) %>%
    mutate(analysis_id = as.character(analysis_id)) %>%
    pivot_wider(names_from = analysis_id, values_from = log2FoldChange) %>%
    select(all_of(as.character(analysis_subset)))

  if (nrow(logfc_matrix) > 1) {
    cor_mat <- cor(logfc_matrix, use = "pairwise.complete.obs", method = "spearman")
    cor_df <- as.data.frame(as.table(cor_mat), stringsAsFactors = FALSE) %>%
      as_tibble() %>%
      rename(analysis_a = Var1, analysis_b = Var2, correlation = Freq) %>%
      mutate(
        analysis_a = as.character(analysis_a),
        analysis_b = as.character(analysis_b)
      ) %>%
      left_join(analysis_info %>% select(analysis_a = analysis_id, analysis_a_label = analysis_label),
                by = "analysis_a") %>%
      left_join(analysis_info %>% select(analysis_b = analysis_id, analysis_b_label = analysis_label),
                by = "analysis_b") %>%
      mutate(
        analysis_a_label = factor(analysis_a_label, levels = label_levels),
        analysis_b_label = factor(analysis_b_label, levels = label_levels)
      )

    write_csv(cor_df, file.path(output_dir, csv_name))

    cor_plot <- ggplot(cor_df, aes(x = analysis_b_label, y = analysis_a_label, fill = correlation)) +
      geom_tile(color = "white", size = 0.3) +
      geom_text(aes(label = scales::number(correlation, accuracy = 0.01)), size = 3) +
      scale_fill_gradient2(low = "#0571B0", mid = "white", high = "#CA0020", midpoint = 0,
                           limits = c(-1, 1), oob = scales::squish, name = "Spearman") +
      theme_minimal(base_size = 13) +
      labs(
        x = NULL,
        y = NULL,
        title = title_text,
        subtitle = "Spearman correlation computed over genes significant in >=1 contrast within group"
      ) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(size = 11),
        legend.position = "top",
        panel.grid = element_blank()
      )

    save_plot_dual(cor_plot, file.path(output_dir, filename_base), width = 9, height = 7)
  }
}

if (length(diet_ids) >= 2) {
  compute_logfc_correlation(
    analysis_subset = diet_ids,
    label_levels = diet_labels,
    filename_base = "log2fc_correlation_heatmap_diet",
    title_text = "Correlation of diet-driven log2 fold-changes",
    csv_name = "pairwise_log2fc_correlations_diet.csv"
  )
}

if (length(sex_ids) >= 2) {
  compute_logfc_correlation(
    analysis_subset = sex_ids,
    label_levels = sex_labels,
    filename_base = "log2fc_correlation_heatmap_sex",
    title_text = "Correlation of sex bias log2 fold-changes",
    csv_name = "pairwise_log2fc_correlations_sex.csv"
  )
}

# Gene-level summary (how often a gene is significant and direction consistency).
gene_summary <- all_results %>%
  group_by(gene, gene_label, gene_display, gene_type) %>%
  summarise(
    analyses_up = sum(direction == "Up"),
    analyses_down = sum(direction == "Down"),
    analyses_total = analyses_up + analyses_down,
    prop_up = if_else(analyses_total > 0, analyses_up / analyses_total, NA_real_),
    max_log2fc = safe_max(log2FoldChange[direction != "NS"]),
    min_log2fc = safe_min(log2FoldChange[direction != "NS"]),
    strongest_padj = safe_min(padj[direction != "NS"]),
    .groups = "drop"
  ) %>%
  mutate(
    gene_display = if_else(is.na(gene_display) | gene_display == "", gene, gene_display),
    strongest_padj = if_else(is.infinite(strongest_padj), NA_real_, strongest_padj)
  ) %>%
  arrange(desc(analyses_total), desc(analyses_up), gene_label)

write_csv(gene_summary, file.path(output_dir, "gene_significance_counts.csv"))

# Tile plot showing direction of change for genes repeatedly significant.
top_gene_count <- 40
top_genes <- gene_summary %>%
  filter(analyses_total > 0) %>%
  slice_head(n = top_gene_count) %>%
  pull(gene)

direction_matrix <- all_results %>%
  filter(gene %in% top_genes) %>%
  group_by(gene, gene_display, analysis_label) %>%
  summarise(
    value = case_when(
      any(direction == "Up") ~ 1,
      any(direction == "Down") ~ -1,
      TRUE ~ 0
    ),
    .groups = "drop"
  ) %>%
  mutate(
    analysis_label = factor(analysis_label, levels = label_order)
  )

gene_order <- gene_summary %>%
  filter(gene %in% top_genes) %>%
  mutate(gene_display = if_else(is.na(gene_display) | gene_display == "", gene, gene_display)) %>%
  arrange(desc(analyses_total), desc(analyses_up)) %>%
  pull(gene_display)

direction_matrix <- direction_matrix %>%
  mutate(
    gene_display = factor(gene_summary$gene_display[match(gene, gene_summary$gene)], levels = rev(unique(gene_order)))
  )

palette_direction_tile <- c("-1" = "#2166AC", "0" = "#F5F5F5", "1" = "#B2182B")

direction_plot <- ggplot(direction_matrix, aes(x = analysis_label, y = gene_display, fill = factor(value))) +
  geom_tile(color = "white", size = 0.2) +
  scale_fill_manual(
    values = palette_direction_tile,
    breaks = c("-1", "0", "1"),
    labels = c("Down", "Not significant", "Up"),
    name = "Direction"
  ) +
  theme_minimal(base_size = 11) +
  labs(
    x = NULL,
    y = NULL,
    title = paste0("Direction of change for top ", top_gene_count, " recurrent genes"),
    subtitle = "Genes ordered by number of analyses with padj <= 0.05"
  ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 8),
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 11),
      legend.position = "top",
      panel.grid = element_blank(),
      plot.margin = margin(20, 26, 20, 20)
    )

save_plot_dual(direction_plot, file.path(output_dir, "top_recurrent_genes_heatmap"), width = 10, height = 12)

# Diet-focused recurrent gene heatmap (>=100 genes if available).
if (length(diet_ids) >= 2) {
  diet_results <- all_results %>%
    filter(analysis_id %in% diet_ids)

  gene_summary_diet <- diet_results %>%
    group_by(gene, gene_label, gene_display, gene_type) %>%
    summarise(
      analyses_up = sum(direction == "Up"),
      analyses_down = sum(direction == "Down"),
      analyses_total = analyses_up + analyses_down,
      prop_up = if_else(analyses_total > 0, analyses_up / analyses_total, NA_real_),
      max_log2fc = safe_max(log2FoldChange[direction != "NS"]),
      min_log2fc = safe_min(log2FoldChange[direction != "NS"]),
      strongest_padj = safe_min(padj[direction != "NS"]),
      .groups = "drop"
    ) %>%
    mutate(
      gene_display = if_else(is.na(gene_display) | gene_display == "", gene, gene_display),
      strongest_padj = if_else(is.infinite(strongest_padj), NA_real_, strongest_padj)
    ) %>%
    arrange(desc(analyses_total), desc(analyses_up), gene_label)

  write_csv(gene_summary_diet, file.path(output_dir, "gene_significance_counts_diet.csv"))

  top_gene_count_diet <- 100
  top_genes_diet <- gene_summary_diet %>%
    filter(analyses_total > 0) %>%
    slice_head(n = top_gene_count_diet) %>%
    pull(gene)

  if (length(top_genes_diet) > 0) {
    direction_matrix_diet <- diet_results %>%
      filter(gene %in% top_genes_diet) %>%
      group_by(gene, gene_display, analysis_label) %>%
      summarise(
        value = case_when(
          any(direction == "Up") ~ 1,
          any(direction == "Down") ~ -1,
          TRUE ~ 0
        ),
        .groups = "drop"
      ) %>%
      mutate(
        analysis_label = factor(analysis_label, levels = diet_labels)
      )

    gene_order_diet <- gene_summary_diet %>%
      filter(gene %in% top_genes_diet) %>%
      mutate(gene_display = if_else(is.na(gene_display) | gene_display == "", gene, gene_display)) %>%
      arrange(desc(analyses_total), desc(analyses_up)) %>%
      pull(gene_display)

    direction_matrix_diet <- direction_matrix_diet %>%
      mutate(
        gene_display = factor(gene_summary_diet$gene_display[match(gene, gene_summary_diet$gene)], levels = rev(unique(gene_order_diet)))
      )

    direction_plot_diet <- ggplot(direction_matrix_diet, aes(x = analysis_label, y = gene_display, fill = factor(value))) +
      geom_tile(color = "white", size = 0.2) +
      scale_fill_manual(
        values = palette_direction_tile,
        breaks = c("-1", "0", "1"),
        labels = c("Down", "Not significant", "Up"),
        name = "Direction"
      ) +
      theme_minimal(base_size = 11) +
      labs(
        x = NULL,
        y = NULL,
        title = "Top 100 recurrent genes (diet contrasts)",
        subtitle = "Genes ordered by number of diet analyses with padj <= 0.05"
      ) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 6.5),
        plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(size = 11),
        legend.position = "top",
        panel.grid = element_blank(),
        plot.margin = margin(20, 28, 20, 20)
      )

    save_plot_dual(direction_plot_diet, file.path(output_dir, "top_recurrent_genes_heatmap_diet_only"), width = 10, height = 14)
  }
}

# Diet-specific uniqueness profiles.
if (length(diet_ids) >= 2) {
  diet_presence <- diet_results %>%
    filter(direction != "NS") %>%
    select(gene, gene_label, gene_type, analysis_id, analysis_label, direction, log2FoldChange, padj) %>%
    distinct(gene, analysis_label, .keep_all = TRUE)

  diet_gene_participation <- diet_presence %>%
    group_by(gene) %>%
    summarise(
      n_analyses = n(),
      gene_label = gene_label[!is.na(gene_label) & gene_label != "" ][1],
      .groups = "drop"
    )

  unique_genes <- diet_gene_participation %>%
    filter(n_analyses == 1) %>%
    pull(gene)

  diet_specificity_counts <- diet_presence %>%
    mutate(category = if_else(gene %in% unique_genes, "Unique to this contrast", "Shared across contrasts")) %>%
    group_by(analysis_label, category) %>%
    summarise(gene_count = n_distinct(gene), .groups = "drop") %>%
    mutate(
      analysis_label = factor(analysis_label, levels = diet_labels),
      category = factor(category, levels = c("Shared across contrasts", "Unique to this contrast"))
    )

  diet_specificity_plot <- ggplot(diet_specificity_counts, aes(x = analysis_label, y = gene_count, fill = category)) +
    geom_col(width = 0.7) +
    coord_flip() +
    scale_fill_manual(values = c("Shared across contrasts" = "#B3CDE3", "Unique to this contrast" = "#FBB4AE")) +
    theme_minimal(base_size = 13) +
    labs(
      x = NULL,
      y = "Number of significant genes",
      fill = NULL,
      title = "Diet-specific gene signatures",
      subtitle = "Counts of genes significant uniquely in each diet contrast versus shared across diet contrasts"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 11),
      axis.text.y = element_text(size = 11),
      legend.position = "top",
      plot.margin = margin(22, 32, 22, 22)
    )

  save_plot_dual(diet_specificity_plot, file.path(output_dir, "diet_specificity_bar"), width = 10, height = 6.5)

  diet_unique_details <- diet_presence %>%
    filter(gene %in% unique_genes) %>%
    arrange(analysis_label, padj) %>%
    mutate(
      gene_label = if_else(!is.na(gene_label) & gene_label != "", gene_label, gene),
      effect_direction = if_else(direction == "Up", "Higher in MCD", "Lower in MCD")
    ) %>%
    select(
      analysis_label,
      gene,
      gene_label,
      gene_type,
      log2FoldChange,
      padj,
      effect_direction
    )

  write_csv(diet_unique_details, file.path(output_dir, "diet_unique_gene_signatures.csv"))

  diet_consistency <- diet_presence %>%
    group_by(gene) %>%
    summarise(
      n_analyses = n(),
      mean_log2fc = mean(log2FoldChange, na.rm = TRUE),
      sd_log2fc = sd(log2FoldChange, na.rm = TRUE),
      min_padj = min(padj, na.rm = TRUE),
      gene_label = first_non_missing(gene_label),
      gene_type = first_non_missing(gene_type),
      .groups = "drop"
    ) %>%
    mutate(
      sd_log2fc = replace_na(sd_log2fc, 0),
      neg_log10_p = if_else(is.finite(min_padj) & min_padj > 0, -log10(min_padj), NA_real_),
      neg_log10_p = replace_na(neg_log10_p, 0),
      gene_label = if_else(is.na(gene_label) | gene_label == "", gene, gene_label)
    )

  write_csv(diet_consistency, file.path(output_dir, "diet_gene_consistency_summary.csv"))

  top_consistent <- diet_consistency %>%
    arrange(desc(n_analyses), desc(neg_log10_p)) %>%
    slice_head(n = 15)

  diet_consistency_plot <- ggplot(diet_consistency, aes(x = n_analyses, y = mean_log2fc)) +
    geom_point(aes(size = neg_log10_p, colour = sd_log2fc), alpha = 0.7) +
    geom_text_repel(
      data = top_consistent,
      aes(label = gene_label),
      size = 3,
      box.padding = 0.35,
      point.padding = 0.2,
      show.legend = FALSE
    ) +
    scale_colour_gradient(low = "#7fc97f", high = "#d73027", name = "SD log2FC") +
    scale_size_continuous(name = "-log10(min padj)", range = c(2, 9)) +
    theme_minimal(base_size = 13) +
    labs(
      x = "Number of diet contrasts with padj <= 0.05",
      y = "Mean log2 fold-change (MCD vs Control)",
      title = "Diet contrast consistency landscape",
      subtitle = "Bubble size scales with statistical strength; colour reflects variability across contrasts"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 11),
      legend.position = "right",
      plot.margin = margin(24, 32, 24, 24)
    )

  save_plot_dual(diet_consistency_plot, file.path(output_dir, "diet_consistency_scatter"), width = 11, height = 7)

  shared_genes <- diet_consistency %>%
    filter(n_analyses == length(diet_ids)) %>%
    arrange(desc(abs(mean_log2fc))) %>%
    mutate(gene_label = if_else(is.na(gene_label) | gene_label == "", gene, gene_label)) %>%
    slice_head(n = min(40, nrow(.)))

  if (nrow(shared_genes) > 0) {
    shared_gene_levels <- rev(shared_genes$gene_label)

    shared_matrix <- diet_results %>%
      filter(gene %in% shared_genes$gene) %>%
      mutate(gene_label = if_else(!is.na(gene_label) & gene_label != "", gene_label, gene)) %>%
      group_by(gene_label, analysis_label) %>%
      summarise(
        mean_log2fc = mean(log2FoldChange, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        analysis_label = factor(analysis_label, levels = diet_labels),
        gene_label = factor(gene_label, levels = shared_gene_levels)
      )

    shared_heatmap <- ggplot(shared_matrix, aes(x = analysis_label, y = gene_label, fill = mean_log2fc)) +
      geom_tile(color = "white", size = 0.25) +
      scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, name = "Mean log2FC") +
      theme_minimal(base_size = 12) +
      labs(
        x = NULL,
        y = NULL,
        title = "Genes consistently altered across all diet contrasts",
        subtitle = "Top 40 by absolute mean log2 fold-change (MCD vs Control)"
      ) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(size = 11),
        legend.position = "top",
        plot.margin = margin(24, 30, 24, 24)
      )

    save_plot_dual(shared_heatmap, file.path(output_dir, "diet_shared_genes_heatmap"), width = 10, height = 9)

    write_csv(shared_matrix, file.path(output_dir, "diet_shared_genes_log2fc.csv"))
  }
}

# Sex-difference transition analysis (Control vs MCD).
sex_ids_transition <- sex_ids
sex_results <- all_results %>%
  filter(analysis_id %in% sex_ids_transition)

if (length(sex_ids_transition) >= 2 && nrow(sex_results) > 0) {
  sex_labels <- analysis_info %>%
    filter(analysis_id %in% sex_ids_transition) %>%
    arrange(match(analysis_id, sex_ids_transition)) %>%
    pull(analysis_label)

  sex_direction_map <- sex_results %>%
    mutate(
      direction_simple = case_when(
        direction == "Up" ~ "Male > Female",
        direction == "Down" ~ "Female > Male",
        TRUE ~ "Not significant"
      ),
      analysis_label = factor(analysis_label, levels = sex_labels)
    )

  sex_matrix <- sex_direction_map %>%
    select(gene, gene_label, analysis_label, direction_simple, log2FoldChange, padj) %>%
    distinct(gene, analysis_label, .keep_all = TRUE) %>%
    mutate(gene_label = if_else(!is.na(gene_label) & gene_label != "", gene_label, gene)) %>%
    pivot_wider(
      names_from = analysis_label,
      values_from = direction_simple,
      values_fill = "Not significant"
    )

  transition_counts <- sex_matrix %>%
    mutate(
      control_state = .data[[sex_labels[1]]],
      mcd_state = .data[[sex_labels[2]]]
    ) %>%
    count(control_state, mcd_state, name = "gene_count")

  transition_plot <- ggplot(transition_counts, aes(x = mcd_state, y = control_state, fill = gene_count)) +
    geom_tile(color = "white", size = 0.4) +
    geom_text(aes(label = gene_count), fontface = "bold", size = 4) +
    scale_fill_gradient(low = "#deebf7", high = "#3182bd", name = "Genes") +
    theme_minimal(base_size = 13) +
    labs(
      x = paste0(sex_labels[2], " classification"),
      y = paste0(sex_labels[1], " classification"),
      title = "Sex bias transitions between diets",
      subtitle = "Counts of genes by male/female bias in control vs MCD comparisons"
    ) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1),
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 11),
      legend.position = "top",
      plot.margin = margin(24, 28, 24, 24)
    )

  save_plot_dual(transition_plot, file.path(output_dir, "sex_bias_transition_heatmap"), width = 8.5, height = 7)

  sex_categories <- sex_matrix %>%
    mutate(
      category = case_when(
        .data[[sex_labels[1]]] == "Male > Female" & .data[[sex_labels[2]]] == "Male > Female" ~ "Male-biased in both diets",
        .data[[sex_labels[1]]] == "Female > Male" & .data[[sex_labels[2]]] == "Female > Male" ~ "Female-biased in both diets",
        .data[[sex_labels[1]]] == "Male > Female" & .data[[sex_labels[2]]] == "Female > Male" ~ paste0("Switches to female bias in MCD"),
        .data[[sex_labels[1]]] == "Female > Male" & .data[[sex_labels[2]]] == "Male > Female" ~ paste0("Switches to male bias in MCD"),
        .data[[sex_labels[1]]] != "Not significant" & .data[[sex_labels[2]]] == "Not significant" ~ paste0("Only significant in ", sex_labels[1]),
        .data[[sex_labels[1]]] == "Not significant" & .data[[sex_labels[2]]] != "Not significant" ~ paste0("Only significant in ", sex_labels[2]),
        TRUE ~ "Not significant in either"
      )
    ) %>%
    select(gene, gene_label, category)

  write_csv(sex_categories, file.path(output_dir, "sex_bias_transition_categories.csv"))

  sex_category_counts <- sex_categories %>%
    count(category, name = "gene_count") %>%
    mutate(category = forcats::fct_reorder(category, gene_count))

  sex_category_plot <- ggplot(sex_category_counts, aes(x = category, y = gene_count)) +
    geom_col(fill = "#9ecae1", width = 0.7) +
    coord_flip() +
    theme_minimal(base_size = 13) +
    labs(
      x = NULL,
      y = "Number of genes",
      title = "How sex bias classifications change between diets",
      subtitle = "Categories derived from control vs MCD male/female contrasts"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 11),
      axis.text.y = element_text(size = 10),
      plot.margin = margin(24, 28, 24, 24)
    )

  save_plot_dual(sex_category_plot, file.path(output_dir, "sex_bias_category_counts"), width = 9, height = 6)

  sex_density_data <- sex_direction_map %>%
    mutate(
      bias = case_when(
        direction == "Up" ~ "Male > Female",
        direction == "Down" ~ "Female > Male",
        TRUE ~ "Not significant"
      ),
      analysis_label = str_wrap(as.character(analysis_label), width = 40)
    )

  sex_density_plot <- ggplot(sex_density_data, aes(x = log2FoldChange, colour = analysis_label, fill = analysis_label)) +
    geom_density(alpha = 0.18, linewidth = 1.1) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey30") +
    theme_minimal(base_size = 13) +
    labs(
      x = "log2 fold-change (Male vs Female)",
      y = "Density",
      colour = NULL,
      fill = NULL,
      title = "Distribution of sex bias across diets",
      subtitle = "Positive values indicate male-biased expression; negative values female-biased"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 11),
      legend.position = "top",
      plot.margin = margin(24, 30, 24, 24)
    )

  save_plot_dual(sex_density_plot, file.path(output_dir, "sex_bias_density"), width = 9.5, height = 6.5)
}

# Focused pairwise scatter plots for key comparisons.
pairwise_comparisons <- tribble(
  ~analysis_x,     ~analysis_y,
  "v1_female",     "v1_male",
  "v1_combined",   "v5_week3",
  "v1_female",     "v5_week3",
  "v1_male",       "v5_week3",
  "v4_control",    "v4_MCD"
)

create_scatter <- function(analysis_x, analysis_y) {
  info_x <- analysis_info %>% filter(analysis_id == analysis_x)
  info_y <- analysis_info %>% filter(analysis_id == analysis_y)
  label_x <- info_x$analysis_label[[1]]
  label_y <- info_y$analysis_label[[1]]
  sig_x_label <- str_wrap(paste("Significant in", label_x), width = 40)
  sig_y_label <- str_wrap(paste("Significant in", label_y), width = 40)
  title_text <- str_wrap(paste("Comparison of log2FC:", label_x, "vs", label_y), width = 70)
  x_axis_label <- str_wrap(paste0("log2 fold-change (", label_x, ")"), width = 45)
  y_axis_label <- str_wrap(paste0("log2 fold-change (", label_y, ")"), width = 45)

  df_x <- all_results %>%
    filter(analysis_id == analysis_x) %>%
    transmute(
      gene,
      gene_label_x = gene_label,
      log2fc_x = log2FoldChange,
      padj_x = padj,
      direction_x = direction,
      gene_display_x = gene_display
    )

  df_y <- all_results %>%
    filter(analysis_id == analysis_y) %>%
    transmute(
      gene,
      gene_label_y = gene_label,
      log2fc_y = log2FoldChange,
      padj_y = padj,
      direction_y = direction,
      gene_display_y = gene_display
    )

  combined <- full_join(df_x, df_y, by = "gene") %>%
    mutate(
      label_text = coalesce(gene_label_x, gene_label_y, gene),
      category = case_when(
        direction_x %in% c("Up", "Down") & direction_y %in% c("Up", "Down") & direction_x == direction_y ~ paste(direction_x, "in both"),
        direction_x %in% c("Up", "Down") & direction_y %in% c("Up", "Down") & direction_x != direction_y ~ "Opposite direction",
        direction_x %in% c("Up", "Down") & (is.na(direction_y) | direction_y == "NS") ~ sig_x_label,
        direction_y %in% c("Up", "Down") & (is.na(direction_x) | direction_x == "NS") ~ sig_y_label,
      TRUE ~ "Not significant"
    ),
    category = factor(
      category,
      levels = c(
        "Up in both",
        "Down in both",
        "Opposite direction",
        sig_x_label,
        sig_y_label,
        "Not significant"
      )
    )
  )

  x_limits <- range(combined$log2fc_x, na.rm = TRUE)
  y_limits <- range(combined$log2fc_y, na.rm = TRUE)
  if (any(is.finite(x_limits))) {
    x_limits <- grDevices::extendrange(x_limits, f = 0.15)
  } else {
    x_limits <- c(-1, 1)
  }
  if (any(is.finite(y_limits))) {
    y_limits <- grDevices::extendrange(y_limits, f = 0.15)
  } else {
    y_limits <- c(-1, 1)
  }

  highlight <- combined %>%
    filter(category %in% c("Up in both", "Down in both", "Opposite direction")) %>%
    arrange(desc(abs(coalesce(log2fc_x, 0)) + abs(coalesce(log2fc_y, 0)))) %>%
    slice_head(n = 15)

  palette_scatter <- c("#C23B75", "#3075B1", "#F2A45E", "#7B3294", "#A6CEE3", "grey80")
  names(palette_scatter) <- c(
    "Up in both",
    "Down in both",
    "Opposite direction",
    sig_x_label,
    sig_y_label,
    "Not significant"
  )

  scatter <- ggplot(combined, aes(x = log2fc_x, y = log2fc_y, color = category)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey80") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey80") +
    geom_abline(slope = 1, intercept = 0, color = "grey60", linetype = "dotted") +
    geom_point(alpha = 0.6, size = 2.2) +
    geom_text_repel(
      data = highlight,
      aes(label = label_text),
      size = 3,
      max.overlaps = Inf,
      box.padding = 0.4,
      point.padding = 0.2,
      show.legend = FALSE
    ) +
    scale_color_manual(values = palette_scatter, guide = guide_legend(override.aes = list(size = 4, alpha = 1))) +
    theme_minimal(base_size = 13) +
    labs(
      x = x_axis_label,
      y = y_axis_label,
      color = NULL,
      title = title_text,
      subtitle = "Colored by significance classification (padj <= 0.05)"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 11),
      legend.position = "top",
      plot.margin = margin(24, 34, 24, 24)
    ) +
    coord_cartesian(xlim = x_limits, ylim = y_limits, expand = TRUE, clip = "off")

  filename_base <- file.path(
    output_dir,
    paste0("scatter_", slugify(info_x$analysis_id), "_vs_", slugify(info_y$analysis_id))
  )
  save_plot_dual(scatter, filename_base, width = 11, height = 8.5)
}

walk2(pairwise_comparisons$analysis_x, pairwise_comparisons$analysis_y, create_scatter)

# Export master table of per-analysis statistics (wide format for log2FC if needed downstream).
wide_logfc <- all_results %>%
  select(gene, gene_label, gene_display, analysis_id, log2FoldChange, padj, direction) %>%
  mutate(analysis_id = as.character(analysis_id)) %>%
  pivot_wider(
    names_from = analysis_id,
    values_from = c(log2FoldChange, padj, direction),
    names_glue = "{analysis_id}_{.value}"
  )

write_tsv(wide_logfc, file.path(output_dir, "log2fc_padj_matrix.tsv"))

cat("Meta-analysis completed. Outputs saved to:", output_dir, "\n")
