#!/usr/bin/env Rscript
# CIRRHOSIS-EXCLUDED hepatocyte phenotype correlation (remediation R1).
#
# Mirrors 505_phenotype_correlation.R, but reads the SEPARATE
# hepatocytes_nocirrhosis/donor_scores.tsv (produced by 501 with
# --exclude-stage Cirrhosis), runs the SAME mixed-effects regression
#   lmer(score ~ axis + (1|dataset))
# for the disease_stage and F_stage axes (plus F_stage_documented_only, NAS
# for parity), BH-adjusts within axis, and writes
#   hepatocytes_nocirrhosis/phenotype_correlations_nocirrhosis.tsv
#
# This script NEVER touches the canonical hepatocytes/ outputs nor the
# pooled phenotype_correlations.tsv.
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(lme4); library(lmerTest)
})

ROOT <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
RES  <- file.path(ROOT, "Analysis/SingleCell/results_gpu_v2/hotspot_modules")
NOCIRR_DIR <- file.path(RES, "hepatocytes_nocirrhosis")
CELL_TYPE  <- "hepatocytes"

donor_meta <- read_tsv(
  file.path(ROOT, "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/donor_metadata_extended.tsv"),
  show_col_types = FALSE) |>
  mutate(disease_stage_ordinal = case_when(
    disease_stage_coarse == "Healthy"         ~ 0,
    disease_stage_coarse == "Steatosis"       ~ 1,
    disease_stage_coarse == "Steatohepatitis" ~ 2,
    disease_stage_coarse == "Cirrhosis"       ~ 3,
    TRUE ~ NA_real_
  ))

# Backwards-compat with 505: tolerate absence of the protocol-contamination flag
# and the clean F-stage column.
if (!"exclude_stage_analysis" %in% names(donor_meta)) {
  donor_meta$exclude_stage_analysis <- FALSE
}
if (!"F_stage_augmented_clean" %in% names(donor_meta) &&
    "F_stage_augmented" %in% names(donor_meta)) {
  donor_meta$F_stage_augmented_clean <- donor_meta$F_stage_augmented
}

# Per-cell-type donor_scores.tsv uses PLAIN INTEGER module ids (not the
# "<celltype>__<int>" namespaced form used by donor_scores_all.tsv).
scores_path <- file.path(NOCIRR_DIR, "donor_scores.tsv")
if (!file.exists(scores_path)) {
  stop(sprintf("Missing %s — run 501 with --exclude-stage Cirrhosis first.", scores_path))
}
scores <- read_tsv(scores_path, show_col_types = FALSE) |>
  mutate(cell_type = CELL_TYPE, module = as.integer(module))

joined <- scores |>
  inner_join(donor_meta, by = "sample") |>
  filter(!is.na(score))

# Defensive: the cirrhosis donors are already gone (excluded at the 501 cell-level
# filter), but re-assert here so the table can never silently include them, and
# also apply the protocol-contamination flag (GSE136103 / Liver_Atlas) for parity
# with 505. NA flags treated as not-excluded so old metadata still works.
n_pre <- dplyr::n_distinct(joined$sample)
joined <- joined |>
  filter(disease_stage_coarse != "Cirrhosis" | is.na(disease_stage_coarse)) |>
  filter(is.na(exclude_stage_analysis) | exclude_stage_analysis == FALSE)
n_post <- dplyr::n_distinct(joined$sample)
cat(sprintf("After cirrhosis + protocol-contamination filters: %d donors retained (was %d, dropped %d)\n",
            n_post, n_pre, n_pre - n_post))
cat(sprintf("Distinct donors actually used (any axis): %d\n", n_post))
cat("Donors by disease_stage_coarse:\n")
print(joined |> distinct(sample, disease_stage_coarse) |> count(disease_stage_coarse))

axes <- list(
  disease_stage = list(col = "disease_stage_ordinal", filter = function(d) d),
  F_stage       = list(col = "F_stage_inferred",      filter = function(d) d),
  F_stage_documented_only = list(
    col = "F_stage_inferred",
    filter = function(d) filter(d, F_stage_source == "documented")
  ),
  NAS           = list(col = "nas_score",             filter = function(d) d)
)

fit_one <- function(df, axis_col) {
  if (sum(!is.na(df[[axis_col]])) < 20 || dplyr::n_distinct(df$dataset) < 2) {
    return(tibble(beta = NA, SE = NA, t = NA, p = NA,
                  n = sum(!is.na(df[[axis_col]])),
                  n_donors = dplyr::n_distinct(df$sample[!is.na(df[[axis_col]])])))
  }
  fm <- as.formula(sprintf("score ~ %s + (1|dataset)", axis_col))
  m <- tryCatch(lmerTest::lmer(fm, data = df, REML = FALSE), error = function(e) NULL)
  if (is.null(m)) return(tibble(beta = NA, SE = NA, t = NA, p = NA,
                                n = nrow(df),
                                n_donors = dplyr::n_distinct(df$sample)))
  cf <- tryCatch(summary(m)$coefficients, error = function(e) NULL)
  if (is.null(cf) || !(axis_col %in% rownames(cf))) {
    return(tibble(beta = NA, SE = NA, t = NA, p = NA,
                  n = nrow(df), n_donors = dplyr::n_distinct(df$sample)))
  }
  row <- cf[axis_col, ]
  p_val <- if ("Pr(>|t|)" %in% colnames(cf)) row["Pr(>|t|)"] else NA_real_
  tibble(beta = unname(row["Estimate"]),
         SE   = unname(row["Std. Error"]),
         t    = unname(row["t value"]),
         p    = unname(p_val),
         n    = nrow(df),
         n_donors = dplyr::n_distinct(df$sample))
}

results <- tibble()
for (ax in names(axes)) {
  cfg <- axes[[ax]]
  df_ax <- cfg$filter(joined)
  per_module <- df_ax |>
    group_by(cell_type, module) |>
    group_modify(~ fit_one(.x, cfg$col)) |>
    ungroup() |>
    mutate(axis = ax,
           q = p.adjust(p, method = "BH"))
  results <- bind_rows(results, per_module)
}

dir.create(NOCIRR_DIR, recursive = TRUE, showWarnings = FALSE)
out_path <- file.path(NOCIRR_DIR, "phenotype_correlations_nocirrhosis.tsv")
write_tsv(results, out_path)
cat(sprintf("Wrote %d rows across %d axes to %s\n",
            nrow(results), length(axes), out_path))
