#!/usr/bin/env Rscript
# Arm A4: recorded histology as the axis. This is the reference denominator for
# Q1, not a competitor.
#
# Two variants:
#   A4        fibrosis stage itself, ties allowed. Every Q1a statistic is rank
#             based, so ties are handled by mid-ranks rather than broken.
#   A4_tie    fibrosis stage with WITHIN-STAGE ties broken by the A1 axis. The
#             difference between A4 and A4_tie isolates exactly one thing: does a
#             derived ordering add resolution INSIDE a histology stage, which is
#             the only place it can add anything the clinical label does not
#             already carry.

suppressPackageStartupMessages({ library(data.table) })
source(file.path(Sys.getenv("MASLD_PROJECT_ROOT"),
                 "scripts/analysis/continuous_axis_benchmark/lib_axis_common.R"))

OUT <- Sys.getenv("CAB_OUT_ROOT", "")
assert_true(nzchar(OUT), "CAB_OUT_ROOT is unset")
pre <- read_prespec()
meta <- load_manifest()

a1_path <- file.path(OUT, "arms", "A1_axis.tsv")
a1 <- if (file.exists(a1_path)) fread(a1_path)[, .(sample_id, a1 = axis_raw)] else NULL

axis_dt <- meta[, .(sample_id, dataset, fibrosis_stage, nas_score)]
axis_dt[, axis_raw := as.numeric(fibrosis_stage)]

if (!is.null(a1)) {
  axis_dt <- merge(axis_dt, a1, by = "sample_id", all.x = TRUE)
  # Rank within stage, then add as a fraction of one stage step, so the tie-break
  # can never reorder two donors across different stages.
  axis_dt[, within := rank_within(a1, fibrosis_stage)]
  axis_dt[, axis_raw_tiebreak := fifelse(is.na(fibrosis_stage), NA_real_,
                                         as.numeric(fibrosis_stage) + 0.99 * (within - 0.5))]
  axis_dt[, within := NULL]
} else {
  axis_dt[, axis_raw_tiebreak := NA_real_]
}

axis_dt[, axis_rank_within_cohort := rank_within(axis_raw, dataset)]
axis_dt[, `:=`(arm = "A4", matrix_id = "recorded_histology")]
write_tsv_once(axis_dt, file.path(OUT, "arms", "A4_axis.tsv"))

popC <- population(meta, "POP-C")
summary_dt <- data.table(
  arm = "A4", matrix_id = "recorded_histology",
  n_with_fibrosis = sum(!is.na(meta$fibrosis_stage)),
  n_pop_c = nrow(popC),
  n_distinct_levels = uniqueN(na.omit(meta$fibrosis_stage)),
  tiebreak_available = !is.null(a1),
  note = "reference denominator for Q1; A4 vs A4_tie isolates within-stage resolution"
)
write_tsv_once(summary_dt, file.path(OUT, "arms", "A4_axis_summary.tsv"))
print(summary_dt)
log_step("ARM_A4_COMPLETE")
