#!/usr/bin/env Rscript
# 03: one tidy table of every axis, plus how much they agree with each other.
#
# The cross-arm agreement matrix is not bookkeeping. If two arms that differ only
# in a preprocessing step produce uncorrelated orderings, then "the axis" is not a
# property of the data, and every downstream comparison is really a comparison of
# analyst choices. That has to be visible before Q1 is judged, not discovered
# afterwards.

suppressPackageStartupMessages({ library(data.table) })
source(file.path(Sys.getenv("MASLD_PROJECT_ROOT"),
                 "scripts/analysis/continuous_axis_benchmark/lib_axis_common.R"))

OUT <- Sys.getenv("CAB_OUT_ROOT", "")
assert_true(nzchar(OUT), "CAB_OUT_ROOT is unset")
meta <- load_manifest()
arms <- c("A1", "A2", "A3", "A4", "A5", "A6")
available <- arms[file.exists(file.path(OUT, "arms", paste0(arms, "_axis.tsv")))]
assert_true(length(available) >= 2L, "Need at least two arms to build a registry")
log_step("registry over arms: ", paste(available, collapse = ", "))

long <- rbindlist(lapply(available, function(a) {
  d <- fread(file.path(OUT, "arms", paste0(a, "_axis.tsv")))
  d[, .(arm = a, sample_id, axis_raw,
        axis_rank_within_cohort = if ("axis_rank_within_cohort" %in% names(d))
          axis_rank_within_cohort else NA_real_)]
}), fill = TRUE)
write_tsv_once(long, file.path(OUT, "arms", "axis_registry_long.tsv.gz"))

wide <- dcast(long, sample_id ~ arm, value.var = "axis_raw")
M <- as.matrix(wide[, -1, with = FALSE])
rownames(M) <- wide$sample_id
S <- cor(M, method = "spearman", use = "pairwise.complete.obs")
write_tsv_once(as.data.table(S, keep.rownames = "arm"),
               file.path(OUT, "arms", "axis_cross_arm_spearman.tsv"))

j <- merge(wide, meta[, .(sample_id, dataset, fibrosis_stage, nas_score)], by = "sample_id")
per_arm <- rbindlist(lapply(available, function(a) {
  v <- j[[a]]
  ok <- is.finite(v)
  data.table(
    arm = a, n_scored = sum(ok), coverage = sum(ok) / nrow(meta),
    spearman_vs_fibrosis = suppressWarnings(
      cor(v, j$fibrosis_stage, method = "spearman", use = "complete.obs")),
    spearman_vs_nas = suppressWarnings(
      cor(v, j$nas_score, method = "spearman", use = "complete.obs")),
    kruskal_eta2_on_cohort = {
      d <- j[ok]
      if (uniqueN(d$dataset) < 2L) NA_real_ else {
        kw <- kruskal.test(d[[a]] ~ droplevels(factor(d$dataset)))
        (unname(kw$statistic) - uniqueN(d$dataset) + 1) / (nrow(d) - uniqueN(d$dataset))
      }
    })
}))
write_tsv_once(per_arm, file.path(OUT, "arms", "axis_registry_summary.tsv"))
print(per_arm)
cat("\ncross-arm Spearman:\n"); print(round(S, 3))
log_step("AXIS_REGISTRY_COMPLETE")
