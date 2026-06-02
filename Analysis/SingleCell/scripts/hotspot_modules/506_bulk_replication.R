#!/usr/bin/env Rscript
# Score Hotspot modules in bulk dream meta-analysis.
# bulk_score[s, m] = sum_g (logFC_meta[g] * hotspot_weight[g, m]) over top-200 genes per module
suppressPackageStartupMessages({ library(dplyr); library(readr); library(tidyr) })

ROOT <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
RES <- file.path(ROOT, "Analysis/SingleCell/results_gpu_v2/hotspot_modules")
DREAM <- file.path(ROOT,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/dream_results_ashr.csv")

dream <- read_csv(DREAM, show_col_types = FALSE) |>
  select(gene = symbol, dream_logFC = logFC, dream_padj = padj) |>
  filter(!is.na(gene), nchar(gene) > 0)

cell_types <- c("global","hepatocytes","macrophages","fibroblasts",
                "endothelial_cells","cholangiocytes","tcells")

agg_rows <- list()
for (ct in cell_types) {
  mg_path <- file.path(RES, ct, "module_genes.tsv")
  if (!file.exists(mg_path)) { message(sprintf("[WARN] skipping %s (missing %s)", ct, mg_path)); next }
  mg <- read_tsv(mg_path, show_col_types = FALSE)
  per_mod <- mg |>
    group_by(module) |>
    slice_max(order_by = weight, n = 200, with_ties = FALSE) |>
    ungroup() |>
    inner_join(dream, by = "gene")
  bulk_score <- per_mod |>
    group_by(module) |>
    summarise(
      bulk_module_score = sum(dream_logFC * weight),
      n_genes_overlapping = n(),
      mean_dream_logFC = mean(dream_logFC, na.rm = TRUE),
      pct_concordant_up = mean(dream_logFC > 0, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(cell_type = ct, .before = 1)
  agg_rows[[ct]] <- bulk_score
}
out <- bind_rows(agg_rows)
write_tsv(out, file.path(RES, "bulk_replication.tsv"))
cat(sprintf("Wrote %d module bulk-replication rows\n", nrow(out)))
