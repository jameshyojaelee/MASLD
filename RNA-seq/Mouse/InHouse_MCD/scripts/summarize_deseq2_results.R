#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(purrr)
  library(glue)
  library(DESeq2)
  library(data.table)
})

root <- normalizePath("./")
metadata_path <- file.path(root, "metadata", "samples.tsv")
counts_path <- file.path(root, "counts", "featurecounts", "gene_counts.txt")
gtf_path <- file.path(root, "reference", "raw", "gencode.vM33.annotation.gtf")
analysis_dirs <- tribble(
  ~analysis_id, ~directory, ~result_file, ~test_type, ~description,
  "v1_female", "analysis_v1_weekPooled_female", "deseq2_results.tsv", "wald", "Diet effect (MCD vs Control) – females",
  "v1_male", "analysis_v1_weekPooled_male", "deseq2_results.tsv", "wald", "Diet effect (MCD vs Control) – males",
  "v1_combined", "analysis_v1_weekPooled_combined", "deseq2_results.tsv", "wald", "Diet effect (MCD vs Control) – pooled sexes",
  "v2_female", "analysis_v2_timeCourse_female", "lrt_results.tsv", "lrt", "Diet × week (time-course) – females",
  "v2_male", "analysis_v2_timeCourse_male", "lrt_results.tsv", "lrt", "Diet × week (time-course) – males",
  "v3_combined", "analysis_v3_timeCourse_combined", "lrt_results.tsv", "lrt", "Diet × week – pooled sexes",
  "v4_control", "analysis_v4_sexDifferences_control", "deseq2_results.tsv", "wald", "Sex effect (Male vs Female) – Control diet",
  "v4_mcd", "analysis_v4_sexDifferences_MCD", "deseq2_results.tsv", "wald", "Sex effect (Male vs Female) – MCD diet",
  "v4_interaction", "analysis_v4_sexInteraction", "deseq2_results.tsv", "wald", "Diet:Sex interaction (Male MCD vs Female Control baseline)",
  "v5_week1", "analysis_v5_week1MCD_vs_controls", "deseq2_results.tsv", "wald", "Diet effect (MCD vs Control) – week 1 MCD cohort",
  "v5_week2", "analysis_v5_week2MCD_vs_controls", "deseq2_results.tsv", "wald", "Diet effect (MCD vs Control) – week 2 MCD cohort",
  "v5_week3", "analysis_v5_week3MCD_vs_controls", "deseq2_results.tsv", "wald", "Diet effect (MCD vs Control) – week 3 MCD cohort"
)

if (!file.exists(gtf_path)) {
  stop("GTF file not found at ", gtf_path)
}

gtf_cols <- c("seqname","source","feature","start","end","score","strand","frame","attribute")
message("Parsing GTF annotations…")
gtf <- fread(gtf_path, sep = "\t", header = FALSE, col.names = gtf_cols, quote = "",
             showProgress = FALSE, data.table = FALSE)
gene_rows <- gtf %>% filter(feature == "gene")
extract_attr <- function(attr, key) {
  m <- str_match(attr, glue('{key} "([^"]+)"'))
  m[, 2]
}
gene_map <- gene_rows %>%
  transmute(
    gene_id = extract_attr(attribute, "gene_id"),
    gene_name = extract_attr(attribute, "gene_name"),
    gene_type = extract_attr(attribute, "gene_type")
  ) %>%
  distinct()
rm(gtf, gene_rows)
gc()
message("Loaded annotations for ", nrow(gene_map), " genes")

summary_rows <- list()
summary_tables <- list()

for (i in seq_len(nrow(analysis_dirs))) {
  entry <- analysis_dirs[i, ]
  dir_path <- file.path(root, entry$directory)
  result_path <- file.path(dir_path, entry$result_file)
  dds_path <- file.path(dir_path, "dds.rds")
  if (!file.exists(result_path)) {
    warning("Missing result file for ", entry$analysis_id, " at ", result_path)
    next
  }
  message("Processing ", entry$analysis_id, " (", entry$description, ")")
  res <- read_tsv(result_path, show_col_types = FALSE)
  lfc_col <- if ("direction_log2FoldChange" %in% colnames(res)) {
    "direction_log2FoldChange"
  } else {
    "log2FoldChange"
  }
  res <- res %>%
    select(-any_of(c("gene_name", "gene_type"))) %>%
    left_join(gene_map, by = c("gene" = "gene_id")) %>%
    relocate(gene_name, gene_type, .after = gene)
  annotated_path <- file.path(dir_path, str_replace(entry$result_file, ".tsv$", "_annotated.tsv"))
  write_tsv(res, annotated_path)

  is_signif <- res$padj < 0.1 & !is.na(res$padj)
  effect_filter <- rep(FALSE, nrow(res))
  if (entry$test_type == "wald") {
    effect_filter <- is_signif & abs(res[[lfc_col]]) >= 1 & res$baseMean >= 10
  } else {
    effect_filter <- is_signif & abs(res[[lfc_col]]) >= 1
  }
  res_filtered <- res[effect_filter, ]
  filtered_path <- file.path(dir_path, str_replace(entry$result_file, ".tsv$", "_filtered.tsv"))
  write_tsv(res_filtered, filtered_path)

  sig_count <- sum(is_signif, na.rm = TRUE)
  filtered_count <- nrow(res_filtered)
  up_count <- if (entry$test_type == "wald") sum(res_filtered$log2FoldChange > 0, na.rm = TRUE) else NA_integer_
  down_count <- if (entry$test_type == "wald") sum(res_filtered$log2FoldChange < 0, na.rm = TRUE) else NA_integer_

  top_up <- res %>%
    filter(is_signif) %>%
    arrange(desc(.data[[lfc_col]])) %>%
    slice_head(n = 5) %>%
    mutate(direction = "up")
  top_down <- res %>%
    filter(is_signif) %>%
    arrange(.data[[lfc_col]]) %>%
    slice_head(n = 5) %>%
    mutate(direction = "down")
  top_genes <- bind_rows(top_up, top_down) %>%
    mutate(analysis_id = entry$analysis_id)
  summary_tables[[entry$analysis_id]] <- top_genes

  if (file.exists(dds_path)) {
    dds <- readRDS(dds_path)
    vsd <- vst(dds, blind = TRUE)
    pca_data <- plotPCA(vsd, intgroup = c("diet", "sex", "week"), returnData = TRUE)
    percentVar <- round(100 * attr(pca_data, "percentVar"), 1)
    pca_out <- pca_data %>% as_tibble() %>%
      mutate(analysis_id = entry$analysis_id)
    write_tsv(pca_out, file.path(dir_path, "pca_coordinates.tsv"))
  } else {
    percentVar <- c(NA_real_, NA_real_)
  }

summary_rows[[entry$analysis_id]] <- tibble(
    analysis_id = entry$analysis_id,
    description = entry$description,
    test_type = entry$test_type,
    directory = entry$directory,
    total_genes = nrow(res),
    significant_padj_lt_0_1 = sig_count,
    effect_size_filtered = filtered_count,
    filtered_up = up_count,
    filtered_down = down_count,
    pc1_percent = percentVar[1],
    pc2_percent = percentVar[2],
    result_file = entry$result_file,
    annotated_results = basename(annotated_path),
    filtered_results = basename(filtered_path)
  )
}

summary_table <- bind_rows(summary_rows)
top_gene_table <- bind_rows(summary_tables) %>%
  select(analysis_id, direction, gene, gene_name, gene_type, log2FoldChange, padj, baseMean)

reports_dir <- file.path(root, "reports")
if (!dir.exists(reports_dir)) dir.create(reports_dir, recursive = TRUE)
summary_path <- file.path(reports_dir, "deseq2_summary.tsv")
top_gene_path <- file.path(reports_dir, "deseq2_top_genes.tsv")
write_tsv(summary_table, summary_path)
write_tsv(top_gene_table, top_gene_path)
message("Wrote summary table to ", summary_path)
message("Wrote top gene table to ", top_gene_path)
