#!/usr/bin/env Rscript
# Score Hotspot modules in the canonical (C2 limma-voom-qw) bulk DEG analysis.
# bulk_score[s, m] = sum_g (logFC_canon[g] * hotspot_weight[g, m]) over top-200 genes per module
# C2 swap 2026-06-08: repointed dream_results_ashr.csv -> canonical_deg_results.csv
# (limma-voom quality-weighted C2), matching 506b_bulk_replication_nocirr.R.
suppressPackageStartupMessages({ library(dplyr); library(readr); library(tidyr) })

ROOT <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
RES <- file.path(ROOT, "Analysis/SingleCell/results_gpu_v2/hotspot_modules")
CANON <- file.path(ROOT,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/canonical_deg_results.csv")

dream <- read_csv(CANON, show_col_types = FALSE) |>
  select(gene = symbol, bulk_logFC = logFC, bulk_padj = padj) |>
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
      bulk_module_score = sum(bulk_logFC * weight),
      n_genes_overlapping = n(),
      mean_bulk_logFC = mean(bulk_logFC, na.rm = TRUE),
      pct_concordant_up = mean(bulk_logFC > 0, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(cell_type = ct, .before = 1)
  agg_rows[[ct]] <- bulk_score
}
out <- bind_rows(agg_rows)
write_tsv(out, file.path(RES, "bulk_replication.tsv"))
cat(sprintf("Wrote %d module bulk-replication rows\n", nrow(out)))
