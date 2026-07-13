#!/usr/bin/env Rscript
# CIRRHOSIS-EXCLUDED hepatocyte bulk replication (remediation R1).
#
# Mirrors 506_bulk_replication.R, but:
#   (1) reads ONLY the cirrhosis-excluded hepatocyte modules
#       (hepatocytes_nocirrhosis/module_genes.tsv), and
#   (2) replicates against the C2 CANONICAL bulk logFC
#       (canonical_deg_results.csv; limma-voom quality-weighted C2),
#       NOT dream_results_ashr.csv.
#
# bulk_score[m] = sum_g (logFC_canon[g] * hotspot_weight[g, m]) over top-200
#                 genes per module (by within-module local-corr weight).
#
# Output: hepatocytes_nocirrhosis/bulk_replication_nocirrhosis.tsv
# Never touches the canonical bulk_replication.tsv.
suppressPackageStartupMessages({ library(dplyr); library(readr); library(tidyr) })

ROOT <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
RES  <- file.path(ROOT, "Analysis/SingleCell/results_gpu_v2/hotspot_modules")
NOCIRR_DIR <- file.path(RES, "hepatocytes_nocirrhosis")
CELL_TYPE  <- "hepatocytes"

# C2 canonical bulk DEG (limma-voom QW, ~ dataset + inferred_sex + group_binary).
CANON <- file.path(ROOT,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/canonical_deg_results.csv")

bulk <- read_csv(CANON, show_col_types = FALSE) |>
  select(gene = symbol, bulk_logFC = logFC, bulk_padj = padj) |>
  filter(!is.na(gene), nchar(gene) > 0)

mg_path <- file.path(NOCIRR_DIR, "module_genes.tsv")
if (!file.exists(mg_path)) {
  stop(sprintf("Missing %s — run 501 with --exclude-stage Cirrhosis first.", mg_path))
}
mg <- read_tsv(mg_path, show_col_types = FALSE)

per_mod <- mg |>
  group_by(module) |>
  slice_max(order_by = weight, n = 200, with_ties = FALSE) |>
  ungroup() |>
  inner_join(bulk, by = "gene")

out <- per_mod |>
  group_by(module) |>
  summarise(
    bulk_module_score   = sum(bulk_logFC * weight),
    n_genes_overlapping = n(),
    mean_bulk_logFC     = mean(bulk_logFC, na.rm = TRUE),
    pct_concordant_up   = mean(bulk_logFC > 0, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(cell_type = CELL_TYPE, .before = 1)

dir.create(NOCIRR_DIR, recursive = TRUE, showWarnings = FALSE)
out_path <- file.path(NOCIRR_DIR, "bulk_replication_nocirrhosis.tsv")
write_tsv(out, out_path)
cat(sprintf("Wrote %d module bulk-replication rows to %s (canonical C2 logFC)\n",
            nrow(out), out_path))
