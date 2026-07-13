#!/usr/bin/env Rscript
# Mixed-effects regression of donor module score on disease_stage / F_stage / NAS.
#   lmer(score ~ axis + (1|dataset))
# Output: phenotype_correlations.tsv with beta, SE, t, p, q per (cell_type, module, axis).
#
# DONOR-LEVEL (canonical 2026-07-12): run-vs-donor pseudoreplication fix. The atlas
# `sample` in donor_scores_all.tsv is a SEQUENCING RUN / sort-fraction library, not a
# biological donor; 4 datasets carry multiple runs per donor (211 runs -> 73 donors in
# the disease-stage fit), so a run-level lmer inflates n and deflates SE. `(1|dataset)`
# does NOT absorb the within-donor correlation. FIX: collapse per-run module scores to
# one score per TRUE biological donor (mean over the donor's runs, via
# lib_donor_collapse::build_srr_to_donor_map) BEFORE the lmer, keeping (1|dataset).
# Module DEFINITIONS (Hotspot co-expression membership) are a co-expression structure
# unaffected by pseudoreplication and are NOT recomputed here. This reproduces the
# donor-level statistic that 511_donorlevel_disease_stage.R writes to
# donor_collapse/phenotype_correlations_donor.tsv. Set HOTSPOT_PHENO_OUT to redirect
# the output (e.g. for a temp validation run without touching the canonical file).
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(lme4); library(lmerTest)
})

ROOT <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
RES <- file.path(ROOT, "Analysis/SingleCell/results_gpu_v2/hotspot_modules")

# Shared run->true-donor map (named vector: SRR/GSM run id -> "{dataset}_{donor}").
source(file.path(ROOT, "Analysis/SingleCell/scripts/lib_donor_collapse.R"))
srr_to_donor <- build_srr_to_donor_map(ROOT)
to_donor <- function(s) ifelse(s %in% names(srr_to_donor), srr_to_donor[s], s)
mode_chr <- function(v) { t <- table(v); if (length(t) == 0L) NA_character_ else names(sort(t, decreasing = TRUE))[1] }

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
# Defensive: the axis columns must exist even if the metadata predates them, so the
# all-NA edge case below collapses to NA rather than erroring.
for (col in c("F_stage_inferred", "nas_score")) if (!col %in% names(donor_meta)) donor_meta[[col]] <- NA_real_
if (!"F_stage_source" %in% names(donor_meta)) donor_meta$F_stage_source <- NA_character_

# ---------------------------------------------------------------------------
# DONOR-LEVEL metadata: one row per biological donor. stage / dataset / exclude
# are invariant within a donor (verified 0 disagreements); F_stage_inferred and
# nas_score are averaged across the donor's runs (all-NA -> NA, not NaN).
# ---------------------------------------------------------------------------
meta_donor <- donor_meta |>
  distinct(sample, .keep_all = TRUE) |>
  mutate(donor = to_donor(sample)) |>
  group_by(donor) |>
  summarise(
    disease_stage_ordinal  = dplyr::first(disease_stage_ordinal),
    dataset                = dplyr::first(dataset),
    exclude_stage_analysis = any(exclude_stage_analysis %in% TRUE),
    F_stage_inferred = if (all(is.na(F_stage_inferred))) NA_real_ else mean(F_stage_inferred, na.rm = TRUE),
    F_stage_source   = mode_chr(F_stage_source),
    nas_score        = if (all(is.na(nas_score))) NA_real_ else mean(nas_score, na.rm = TRUE),
    .groups = "drop"
  )

# DONOR-LEVEL scores: collapse per-run module score to the unweighted mean per
# (cell_type, module, donor). Module DEFINITIONS are unchanged (scores reused verbatim).
scores <- read_tsv(file.path(RES, "donor_scores_all.tsv"), show_col_types = FALSE) |>
  mutate(donor = to_donor(sample)) |>
  group_by(cell_type, module, donor) |>
  summarise(score = mean(score, na.rm = TRUE), .groups = "drop")

joined <- scores |>
  inner_join(meta_donor, by = "donor") |>
  filter(!is.na(score))

# Drop GSE136103 / Liver_Atlas (CD45+ enrichment + non-NAFLD protocol) before
# any stage / progression correlation. See donor_metadata_extended.tsv for the
# canonical flag; this is the protocol-contamination exclusion. NA flags are
# treated as not-excluded so old metadata still works.
n_pre <- dplyr::n_distinct(joined$donor)
joined <- joined |>
  filter(is.na(exclude_stage_analysis) | exclude_stage_analysis == FALSE)
n_post <- dplyr::n_distinct(joined$donor)
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

out_pheno <- Sys.getenv("HOTSPOT_PHENO_OUT", unset = file.path(RES, "phenotype_correlations.tsv"))
write_tsv(results, out_pheno)
cat(sprintf("Wrote %d rows across %d axes to %s\n", nrow(results), length(axes), out_pheno))
