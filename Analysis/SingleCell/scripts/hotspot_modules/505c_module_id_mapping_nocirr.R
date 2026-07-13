#!/usr/bin/env Rscript
# POST-HOC module-ID mapping for the cirrhosis-excluded hepatocyte re-run (R1).
#
# Hotspot module integer ids are arbitrary and NOT stable across runs, so the
# nocirrhosis modules do not carry the same numbers as the canonical
# hep-19/20/24/26/27 named in the manuscript. This maps each nocirrhosis module
# to its best-matching canonical hepatocyte module by gene-membership Jaccard
# (reciprocal-best-match flagged), then attaches the refreshed phenotype
# beta/q (505b) and bulk-replication score (506b) so the manuscript's NAMED
# modules can be reported with their cirrhosis-excluded statistics.
#
# RUN ONLY AFTER the chain (17133007 -> 17133008) completes, i.e. once
# hepatocytes_nocirrhosis/{module_genes.tsv,phenotype_correlations_nocirrhosis.tsv,
# bulk_replication_nocirrhosis.tsv} all exist.
#
# Output:
#   hepatocytes_nocirrhosis/module_id_map_to_canonical.tsv
#       nocirr_module, best_canonical_module, jaccard, n_shared, n_nocirr, n_canon,
#       reciprocal_best, + refreshed disease_stage/F_stage beta & q + bulk score
suppressPackageStartupMessages({ library(dplyr); library(readr); library(tidyr); library(purrr) })

ROOT <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
RES  <- file.path(ROOT, "Analysis/SingleCell/results_gpu_v2/hotspot_modules")
CANON_DIR  <- file.path(RES, "hepatocytes")
NOCIRR_DIR <- file.path(RES, "hepatocytes_nocirrhosis")

# Manuscript-named canonical hepatocyte modules (for convenience flagging).
NAMED_CANONICAL <- c(19, 20, 24, 26, 27)

stopifnot(file.exists(file.path(CANON_DIR,  "module_genes.tsv")))
stopifnot(file.exists(file.path(NOCIRR_DIR, "module_genes.tsv")))

canon  <- read_tsv(file.path(CANON_DIR,  "module_genes.tsv"), show_col_types = FALSE)
nocirr <- read_tsv(file.path(NOCIRR_DIR, "module_genes.tsv"), show_col_types = FALSE)

canon_sets  <- split(canon$gene,  as.integer(canon$module))
nocirr_sets <- split(nocirr$gene, as.integer(nocirr$module))

jaccard <- function(a, b) {
  i <- length(intersect(a, b)); u <- length(union(a, b))
  if (u == 0) return(0); i / u
}

# Full Jaccard matrix (nocirr x canonical)
grid <- expand_grid(nocirr_module = as.integer(names(nocirr_sets)),
                    canonical_module = as.integer(names(canon_sets))) |>
  mutate(
    n_shared = map2_int(nocirr_module, canonical_module,
                        ~ length(intersect(nocirr_sets[[as.character(.x)]],
                                           canon_sets[[as.character(.y)]]))),
    jaccard  = map2_dbl(nocirr_module, canonical_module,
                        ~ jaccard(nocirr_sets[[as.character(.x)]],
                                  canon_sets[[as.character(.y)]])),
    n_nocirr = map_int(nocirr_module,   ~ length(nocirr_sets[[as.character(.x)]])),
    n_canon  = map_int(canonical_module, ~ length(canon_sets[[as.character(.x)]]))
  )

# Best canonical match per nocirr module + reciprocal-best check
best_for_nocirr <- grid |> group_by(nocirr_module) |>
  slice_max(jaccard, n = 1, with_ties = FALSE) |> ungroup()
best_for_canon  <- grid |> group_by(canonical_module) |>
  slice_max(jaccard, n = 1, with_ties = FALSE) |> ungroup() |>
  transmute(canonical_module, canon_best_nocirr = nocirr_module)

mapping <- best_for_nocirr |>
  rename(best_canonical_module = canonical_module) |>
  left_join(best_for_canon, by = c("best_canonical_module" = "canonical_module")) |>
  mutate(reciprocal_best = nocirr_module == canon_best_nocirr,
         canonical_is_named = best_canonical_module %in% NAMED_CANONICAL) |>
  select(nocirr_module, best_canonical_module, jaccard, n_shared,
         n_nocirr, n_canon, reciprocal_best, canonical_is_named)

# Attach refreshed phenotype beta/q (disease_stage + F_stage) and bulk score.
pheno_path <- file.path(NOCIRR_DIR, "phenotype_correlations_nocirrhosis.tsv")
if (file.exists(pheno_path)) {
  pheno <- read_tsv(pheno_path, show_col_types = FALSE) |>
    filter(axis %in% c("disease_stage", "F_stage")) |>
    mutate(module = as.integer(module)) |>
    pivot_wider(id_cols = module, names_from = axis,
                values_from = c(beta, q, n_donors),
                names_glue = "{axis}_{.value}")
  mapping <- mapping |> left_join(pheno, by = c("nocirr_module" = "module"))
}
bulk_path <- file.path(NOCIRR_DIR, "bulk_replication_nocirrhosis.tsv")
if (file.exists(bulk_path)) {
  bulk <- read_tsv(bulk_path, show_col_types = FALSE) |>
    mutate(module = as.integer(module)) |>
    select(module, bulk_module_score, n_genes_overlapping, pct_concordant_up)
  mapping <- mapping |> left_join(bulk, by = c("nocirr_module" = "module"))
}

out_path <- file.path(NOCIRR_DIR, "module_id_map_to_canonical.tsv")
write_tsv(mapping, out_path)
cat(sprintf("Wrote %d nocirr->canonical module mappings to %s\n", nrow(mapping), out_path))
cat("\nMatches to manuscript-named canonical modules (19/20/24/26/27):\n")
print(mapping |> filter(canonical_is_named) |>
      select(nocirr_module, best_canonical_module, jaccard, n_shared,
             reciprocal_best, disease_stage_beta, disease_stage_q,
             F_stage_beta, F_stage_q))
