#!/usr/bin/env Rscript
# Arm A3: the sanctioned comparator. No new ordering is inferred here.
#
# The frozen prespecified NMF programme scores already exist and are described in
# the Resource as continuous interpretive axes. Their first principal component is
# the strongest continuous readout the project already owns, so it is the bar any
# newly derived axis has to clear. Taking PC1 of all six programmes rather than
# picking one programme avoids post-hoc selection.
#
# Coverage is 685 of 844, so every statistic computed against this arm is also
# recomputed for the other arms restricted to the same 685 donors, otherwise the
# comparison is confounded by sample size.

suppressPackageStartupMessages({ library(data.table) })
source(file.path(Sys.getenv("MASLD_PROJECT_ROOT"),
                 "scripts/analysis/continuous_axis_benchmark/lib_axis_common.R"))

OUT <- Sys.getenv("CAB_OUT_ROOT", "")
assert_true(nzchar(OUT), "CAB_OUT_ROOT is unset")
pre <- read_prespec()
set.seed(pre$seeds$master)
meta <- load_manifest()

build <- function(path, cols, label) {
  nmf <- fread(path)
  assert_true(all(c("sample_id", cols) %in% names(nmf)),
              paste0("Expected columns absent from ", path))
  # P1..P6 only. The P_*_score columns in this file duplicate them and are not
  # independent features.
  dt <- nmf[, c("sample_id", cols), with = FALSE]
  dt <- dt[sample_id %in% meta$sample_id]
  m <- as.matrix(dt[, ..cols])
  m <- scale(m, center = TRUE, scale = TRUE)
  keep <- complete.cases(m)
  p <- prcomp(m[keep, , drop = FALSE], center = FALSE, scale. = FALSE)
  out <- data.table(sample_id = dt$sample_id[keep], axis_raw = p$x[, 1])
  # Orient so that higher axis means more disease, using the frozen cohort-adjusted
  # association with fibrosis. Sign is a display convention, not a result.
  j <- merge(out, meta[, .(sample_id, fibrosis_stage)], by = "sample_id")
  r <- suppressWarnings(cor(j$axis_raw, j$fibrosis_stage,
                            method = "spearman", use = "complete.obs"))
  if (is.finite(r) && r < 0) out[, axis_raw := -axis_raw]
  list(axis = out, pct_var = 100 * p$sdev[1]^2 / sum(p$sdev^2), n = nrow(out), label = label)
}

k6 <- build(substrate_path("nmf_k6"), paste0("P", 1:6), "k6")
k4 <- build(substrate_path("nmf_k4"), paste0("P", 1:4), "k4")

axis_dt <- merge(meta[, .(sample_id, dataset, fibrosis_stage, nas_score)],
                 k6$axis, by = "sample_id", all.x = TRUE)
axis_dt <- merge(axis_dt, setnames(copy(k4$axis), "axis_raw", "axis_raw_k4"),
                 by = "sample_id", all.x = TRUE)
axis_dt[, axis_rank_within_cohort := rank_within(axis_raw, dataset)]
axis_dt[, `:=`(arm = "A3", matrix_id = "frozen_nmf")]
write_tsv_once(axis_dt, file.path(OUT, "arms", "A3_axis.tsv"))

summary_dt <- data.table(
  arm = "A3", matrix_id = "frozen_nmf",
  n_k6 = k6$n, n_k4 = k4$n, n_manifest = nrow(meta),
  coverage_k6 = k6$n / nrow(meta),
  pc1_pct_variance_k6 = k6$pct_var, pc1_pct_variance_k4 = k4$pct_var,
  spearman_k6_vs_k4 = cor(axis_dt$axis_raw, axis_dt$axis_raw_k4,
                          method = "spearman", use = "complete.obs"),
  note = "no new ordering inferred; frozen prespecified programme scores only"
)
write_tsv_once(summary_dt, file.path(OUT, "arms", "A3_axis_summary.tsv"))
print(summary_dt)
log_step("ARM_A3_COMPLETE")
