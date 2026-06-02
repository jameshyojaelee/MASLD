#!/usr/bin/env Rscript
# Build all_modules.tsv (canonical module catalog) and append 6 columns to
# multi_evidence_atlas.csv. Reads ALL outputs of stages 504/505/506/507/loo.
suppressPackageStartupMessages({ library(dplyr); library(readr); library(tidyr); library(purrr) })

ROOT <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
RES <- file.path(ROOT, "Analysis/SingleCell/results_gpu_v2/hotspot_modules")
ATLAS <- file.path(ROOT, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")

# ---- Build canonical all_modules.tsv ----
cell_types <- c("global","hepatocytes","macrophages","fibroblasts",
                "endothelial_cells","cholangiocytes","tcells")

module_genes_long <- map_dfr(cell_types, function(ct) {
  p <- file.path(RES, ct, "module_genes.tsv")
  if (!file.exists(p)) { message(sprintf("[WARN] skipping %s (missing %s)", ct, p)); return(tibble()) }
  read_tsv(p, show_col_types = FALSE) |>
    mutate(cell_type = ct, .before = 1)
})

novelty <- read_tsv(file.path(RES, "novelty_matches.tsv"), show_col_types = FALSE)
loo <- if (file.exists(file.path(RES, "loo_stability.tsv")))
  read_tsv(file.path(RES, "loo_stability.tsv"), show_col_types = FALSE) else
  tibble(cell_type = character(), module = integer(), stability_score = numeric(), stability_fail = logical())
pheno <- read_tsv(file.path(RES, "phenotype_correlations.tsv"), show_col_types = FALSE)
bulk <- read_tsv(file.path(RES, "bulk_replication.tsv"), show_col_types = FALSE)
val  <- if (file.exists(file.path(RES, "olink_visium_validation.tsv")))
  read_tsv(file.path(RES, "olink_visium_validation.tsv"), show_col_types = FALSE) else
  tibble(cell_type = character(), module = integer())

# Normalize the module column to integer across all inputs.
# 505 (which reads donor_scores_all.tsv) inherits 502's namespaced
# "<celltype>__<int>" form; 504, 506 use per-CT module_genes.tsv where
# module is plain int. Strip "<celltype>__" prefix if present.
strip_module_ns <- function(df) {
  if (!"module" %in% names(df)) return(df)
  if (is.character(df$module)) {
    # Remove everything up to and including the LAST "__"; handles
    # cell-type names that themselves contain underscores
    # (e.g. "endothelial_cells__1" -> "1").
    df$module <- as.integer(sub("^.*__", "", df$module))
  } else {
    df$module <- as.integer(df$module)
  }
  df
}
novelty <- strip_module_ns(novelty)
loo <- strip_module_ns(loo)
pheno <- strip_module_ns(pheno)
bulk <- strip_module_ns(bulk)
val <- strip_module_ns(val)

# Reshape pheno: one row per (ct, module) with disease_stage and F_stage hits
pheno_wide <- pheno |>
  filter(axis %in% c("disease_stage", "F_stage", "F_stage_documented_only", "NAS")) |>
  pivot_wider(id_cols = c(cell_type, module),
              names_from = axis, values_from = c(beta, q),
              names_glue = "{axis}_{.value}")

# Tag progression_module: significant in disease_stage AND F_stage with concordant direction
modules_summary <- pheno_wide |>
  mutate(progression_module = (
    !is.na(disease_stage_q) & disease_stage_q < 0.05 &
    !is.na(F_stage_q) & F_stage_q < 0.05 &
    sign(disease_stage_beta) == sign(F_stage_beta)
  )) |>
  mutate(activity_module = (
    !is.na(disease_stage_q) & disease_stage_q < 0.05 &
    (!is.na(F_stage_q) & F_stage_q >= 0.10)
  ))

all_modules <- novelty |>
  left_join(loo, by = c("cell_type", "module")) |>
  left_join(modules_summary, by = c("cell_type", "module")) |>
  left_join(bulk, by = c("cell_type", "module")) |>
  left_join(val, by = c("cell_type", "module")) |>
  mutate(bulk_replicated = !is.na(bulk_module_score) &
                            sign(bulk_module_score) == sign(disease_stage_beta))

write_tsv(all_modules, file.path(RES, "all_modules.tsv"))
cat(sprintf("Wrote all_modules.tsv with %d rows\n", nrow(all_modules)))

# ---- Build per-gene summary then merge into atlas ----
gene_summary <- module_genes_long |>
  group_by(gene) |>
  summarise(
    hotspot_n_modules = n(),
    .groups = "drop"
  )

# Top module per gene (highest weight across all CT runs)
top_module <- module_genes_long |>
  group_by(gene) |>
  slice_max(weight, n = 1, with_ties = FALSE) |>
  transmute(gene,
            hotspot_top_module = sprintf("%s__%d", cell_type, module),
            hotspot_top_module_weight = weight,
            cell_type, module)

# Annotate top_module with progression/novelty/bulk_replicated flags
top_module <- top_module |>
  left_join(all_modules |>
              transmute(cell_type, module,
                        is_novel, progression_module, bulk_replicated),
            by = c("cell_type", "module")) |>
  transmute(gene, hotspot_top_module, hotspot_top_module_weight,
            hotspot_novel_member = is_novel,
            hotspot_progression_member = progression_module,
            hotspot_bulk_replicated = bulk_replicated)

per_gene <- gene_summary |> left_join(top_module, by = "gene")

# Append to atlas
atlas <- read_csv(ATLAS, show_col_types = FALSE)
new_cols <- c("hotspot_n_modules", "hotspot_top_module", "hotspot_top_module_weight",
              "hotspot_novel_member", "hotspot_progression_member", "hotspot_bulk_replicated")
existing <- intersect(new_cols, names(atlas))
if (length(existing) > 0) {
  message("Removing existing hotspot columns before re-add: ", paste(existing, collapse = ", "))
  atlas <- atlas |> select(-all_of(existing))
}

# Auto-detect the gene join column (multi_evidence_atlas.csv uses human_symbol)
gene_col <- intersect(c("human_symbol", "gene_symbol", "gene_name", "gene"), names(atlas))[1]
if (is.na(gene_col)) stop("Could not find gene column in multi_evidence_atlas.csv")
atlas_out <- atlas |> left_join(per_gene, by = setNames("gene", gene_col))

write_csv(atlas_out, ATLAS)
cat(sprintf("Atlas %d genes x %d cols (was %d); +%d hotspot columns\n",
            nrow(atlas_out), ncol(atlas_out), ncol(atlas), length(new_cols)))
