#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(DESeq2)
  library(readr)
  library(dplyr)
  library(tibble)
  library(stringr)
})

root <- normalizePath(".")
metadata_path <- file.path(root, "metadata", "samples.tsv")
counts_path <- file.path(root, "counts", "featurecounts", "gene_counts.txt")
out_dir <- file.path(root, "analysis_mcd_vs_control")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

samples <- read_tsv(metadata_path, col_types = cols()) %>%
  mutate(
    diet = factor(diet, levels = c("Control", "MCD"))
  )
samples_df <- as.data.frame(samples)
rownames(samples_df) <- samples_df$sample_id

raw_counts <- read_tsv(counts_path, comment = "#", col_types = cols())

bam_paths <- file.path(
  root,
  "alignments",
  "star",
  samples$sample_id,
  paste0(samples$sample_id, ".Aligned.sortedByCoord.out.bam")
)
path_to_sample <- setNames(samples$sample_id, bam_paths)
sample_cols <- intersect(colnames(raw_counts), names(path_to_sample))
colnames(raw_counts)[match(sample_cols, colnames(raw_counts))] <- path_to_sample[sample_cols]

count_matrix <- raw_counts %>%
  select(Geneid, any_of(samples$sample_id)) %>%
  column_to_rownames(var = "Geneid") %>%
  as.matrix()

stopifnot(all(samples$sample_id %in% colnames(count_matrix)))
count_matrix <- count_matrix[, samples$sample_id]

design_formula <- ~ diet
if ("genotype" %in% colnames(samples_df)) {
  geno_levels <- samples_df$genotype[!is.na(samples_df$genotype) & samples_df$genotype != ""]
  if (length(unique(geno_levels)) > 1) {
    samples_df$genotype <- factor(samples_df$genotype)
    design_formula <- ~ genotype + diet
  }
}

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = samples_df,
  design = design_formula
)
dds <- DESeq(dds)
res <- results(dds, contrast = c("diet", "MCD", "Control"))

res_tbl <- as.data.frame(res) %>%
  rownames_to_column("gene")

gtf_path <- file.path(root, "reference", "raw", "gencode.vM33.annotation.gtf")
if (file.exists(gtf_path)) {
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
  res_tbl <- res_tbl %>%
    left_join(gene_annot, by = "gene") %>%
    select(gene, gene_name, gene_type, everything())
}

write_tsv(res_tbl, file.path(out_dir, "deseq2_mcd_vs_control.tsv"))

norm_counts <- counts(dds, normalized = TRUE) %>%
  as.data.frame() %>%
  rownames_to_column("gene")
write_csv(norm_counts, file.path(out_dir, "normalized_counts.csv"))

writeLines(
  paste0("design: ", as.character(design_formula)),
  con = file.path(out_dir, "design.txt")
)
