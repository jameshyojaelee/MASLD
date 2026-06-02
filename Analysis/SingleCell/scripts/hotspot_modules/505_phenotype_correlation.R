#!/usr/bin/env Rscript
# Mixed-effects regression of donor module score on disease_stage / F_stage / NAS.
#   lmer(score ~ axis + (1|dataset))
# Output: phenotype_correlations.tsv with beta, SE, t, p, q per (cell_type, module, axis).
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(lme4); library(lmerTest)
})

ROOT <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
RES <- file.path(ROOT, "Analysis/SingleCell/results_gpu_v2/hotspot_modules")

donor_meta <- read_tsv(file.path(ROOT, "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/donor_metadata_extended.tsv"),
                       show_col_types = FALSE) |>
  mutate(disease_stage_ordinal = case_when(
    disease_stage_coarse == "Healthy" ~ 0,
    disease_stage_coarse == "Steatosis" ~ 1,
    disease_stage_coarse == "Steatohepatitis" ~ 2,
    disease_stage_coarse == "Cirrhosis" ~ 3,
    TRUE ~ NA_real_
  ))

# Backwards-compat: handle absence of protocol-contamination flag and clean F-stage
# columns gracefully so the script still runs if 344 has not been re-run.
if (!"exclude_stage_analysis" %in% names(donor_meta)) {
  donor_meta$exclude_stage_analysis <- FALSE
}
if (!"F_stage_augmented_clean" %in% names(donor_meta) &&
    "F_stage_augmented" %in% names(donor_meta)) {
  donor_meta$F_stage_augmented_clean <- donor_meta$F_stage_augmented
}

scores <- read_tsv(file.path(RES, "donor_scores_all.tsv"), show_col_types = FALSE)

joined <- scores |>
  inner_join(donor_meta, by = "sample") |>
  filter(!is.na(score))

# Drop GSE136103 / Liver_Atlas (CD45+ enrichment + non-NAFLD protocol) before
# any stage / progression correlation. See donor_metadata_extended.tsv for the
# canonical flag; this is the protocol-contamination exclusion. NA flags are
# treated as not-excluded so old metadata still works.
n_pre <- dplyr::n_distinct(joined$sample)
joined <- joined |>
  filter(is.na(exclude_stage_analysis) | exclude_stage_analysis == FALSE)
n_post <- dplyr::n_distinct(joined$sample)
cat(sprintf("After protocol-contamination filter: %d donors retained (was %d, dropped %d)\n",
            n_post, n_pre, n_pre - n_post))

axes <- list(
  disease_stage = list(col = "disease_stage_ordinal", filter = function(d) d),
  F_stage       = list(col = "F_stage_inferred",     filter = function(d) d),
  F_stage_documented_only = list(
    col = "F_stage_inferred",
    filter = function(d) filter(d, F_stage_source == "documented")
  ),
  NAS           = list(col = "nas_score",             filter = function(d) d)
)

fit_one <- function(df, axis_col) {
  if (sum(!is.na(df[[axis_col]])) < 20 || dplyr::n_distinct(df$dataset) < 2) {
    return(tibble(beta = NA, SE = NA, t = NA, p = NA, n = sum(!is.na(df[[axis_col]]))))
  }
  fm <- as.formula(sprintf("score ~ %s + (1|dataset)", axis_col))
  m <- tryCatch(lmerTest::lmer(fm, data = df, REML = FALSE), error = function(e) NULL)
  if (is.null(m)) return(tibble(beta = NA, SE = NA, t = NA, p = NA, n = nrow(df)))
  # Extract fixed-effect row for axis_col directly from lmerTest summary
  # (avoids broom.mixed dep, which is not in rnaseq env).
  cf <- tryCatch(summary(m)$coefficients, error = function(e) NULL)
  if (is.null(cf) || !(axis_col %in% rownames(cf))) {
    return(tibble(beta = NA, SE = NA, t = NA, p = NA, n = nrow(df)))
  }
  row <- cf[axis_col, ]
  # lmerTest::lmer summary uses "Pr(>|t|)"; lme4::lmer uses no p column
  p_val <- if ("Pr(>|t|)" %in% colnames(cf)) row["Pr(>|t|)"] else NA_real_
  tibble(beta = unname(row["Estimate"]),
         SE   = unname(row["Std. Error"]),
         t    = unname(row["t value"]),
         p    = unname(p_val),
         n    = nrow(df))
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

write_tsv(results, file.path(RES, "phenotype_correlations.tsv"))
cat(sprintf("Wrote %d rows across %d axes\n", nrow(results), length(axes)))
