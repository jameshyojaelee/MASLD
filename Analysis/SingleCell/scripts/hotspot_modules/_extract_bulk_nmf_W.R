#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(NMF); library(readr); library(tibble); library(dplyr) })
ROOT <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
cache <- readRDS(file.path(ROOT, "RNA-seq/results/subtypes/nmf_results_cache_clean.rds"))
# cache is a list: $nmf_results (keyed by k), $metrics, $mat, $mat_nn
nm <- cache$nmf_results[["6"]]           # canonical k=6 NMFfitX1 object
W <- NMF::basis(nm)                      # gene x k matrix (5000 x 6)
colnames(W) <- paste0("P", seq_len(ncol(W)))
W_df <- as_tibble(W, rownames = "ensembl_versioned")

# Map Ensembl versioned IDs to gene symbols so the namespace matches
# hotspot module_genes.tsv (which uses symbols).
biotypes <- read_tsv(file.path(ROOT, "data/gencode_v49_gene_metadata.tsv.gz"),
                     show_col_types = FALSE) |>
  select(ensembl_versioned = gene_id, gene = gene_name) |>
  distinct(ensembl_versioned, .keep_all = TRUE)

W_sym <- W_df |>
  left_join(biotypes, by = "ensembl_versioned") |>
  filter(!is.na(gene), nchar(gene) > 0) |>
  group_by(gene) |>
  summarise(across(starts_with("P"), max), .groups = "drop") |>
  relocate(gene)

write_tsv(W_sym, file.path(ROOT, "data/reference_signatures/bulk_nmf_k6_W.tsv"))
message(sprintf("Wrote %dx%d W matrix (gene-symbol indexed)", nrow(W_sym), ncol(W_sym) - 1))
