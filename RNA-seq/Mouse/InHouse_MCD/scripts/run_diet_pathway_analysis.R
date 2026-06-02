#!/usr/bin/env Rscript

# Pathway enrichment for diet-consistent responders.

suppressPackageStartupMessages({
  if (!requireNamespace("gprofiler2", quietly = TRUE)) {
    install.packages("gprofiler2", repos = "https://cloud.r-project.org")
  }
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(forcats)
  library(gprofiler2)
})

root <- normalizePath("./")
input_path <- file.path(root, "meta_analysis", "diet_consistency_top_hits.csv")
output_dir <- file.path(root, "meta_analysis")

if (!file.exists(input_path)) {
  stop("Input file not found: ", input_path, ". Run the meta-analysis pipeline first.")
}

genes_tbl <- read_csv(input_path, show_col_types = FALSE) %>%
  mutate(gene_label = if_else(is.na(gene_label) | gene_label == "", gene, gene_label))

gene_list <- genes_tbl %>%
  pull(gene_label) %>%
  unique() %>%
  na.omit()

if (length(gene_list) == 0) {
  stop("No genes available for enrichment after filtering.")
}

message("Running g:Profiler enrichment on ", length(gene_list), " genes...")
gost_res <- tryCatch(
  gost(
    query = gene_list,
    organism = "mmusculus",
    ordered_query = TRUE,
    multi_query = FALSE,
    correction_method = "fdr",
    sources = c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC", "WP", "CORUM")
  ),
  error = function(e) {
    stop("g:Profiler request failed: ", e$message)
  }
)

if (is.null(gost_res) || nrow(gost_res$result) == 0) {
  warning("No enrichment terms returned.")
  quit(status = 0)
}

enrichment_tbl <- gost_res$result %>%
  transmute(
    source,
    term_id,
    term_name,
    adjusted_p_value = p_value,
    term_size,
    intersection_size,
    query_size,
    precision = intersection_size / term_size,
    recall = intersection_size / query_size
  ) %>%
  arrange(adjusted_p_value, desc(intersection_size))

output_tsv <- file.path(output_dir, "diet_consistency_top_hits_enrichment.tsv")
write_tsv(enrichment_tbl, output_tsv)
message("Saved enrichment table to ", output_tsv)

top_terms <- enrichment_tbl %>%
  filter(adjusted_p_value <= 0.1) %>%
  slice_head(n = 15)

if (nrow(top_terms) > 0) {
  plot_data <- top_terms %>%
    mutate(
      term_display = str_wrap(paste0(source, ": ", term_name), width = 60),
      term_display = fct_reorder(term_display, adjusted_p_value, .desc = TRUE)
    )

  plot <- ggplot(plot_data, aes(x = term_display, y = -log10(adjusted_p_value), fill = source)) +
    geom_col(width = 0.7) +
    coord_flip() +
    theme_minimal(base_size = 13) +
    labs(
      x = NULL,
      y = "-log10(FDR)",
      title = "Enrichment of stable MCD responders",
      subtitle = "Top pathways among genes significant in all diet contrasts",
      fill = "Source"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 11),
      plot.margin = margin(24, 24, 24, 24)
    )

  ggsave(
    filename = file.path(output_dir, "diet_consistency_top_hits_enrichment.png"),
    plot = plot, width = 10, height = 7, dpi = 300
  )
  ggsave(
    filename = file.path(output_dir, "diet_consistency_top_hits_enrichment.pdf"),
    plot = plot, width = 10, height = 7
  )
  message("Saved enrichment plot to meta_analysis/diet_consistency_top_hits_enrichment.[png|pdf]")
} else {
  message("No terms passed FDR ≤ 0.1; enrichment plot skipped.")
}

message("Pathway enrichment completed.")
