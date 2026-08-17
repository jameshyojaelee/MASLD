#!/usr/bin/env Rscript
# 05: the decisive Q1 test. Does a derived axis explain expression that recorded
# histology does not, and vice versa?
#
# Both directions are read off ONE fit containing both terms, so the comparison is
# symmetric and neither term gets the advantage of being fitted alone:
#
#   AH : ~ dataset + inferred_sex + factor(fibrosis_stage) + ns(axis_rank, df = 4)
#
#   axis | histology   moderated F over the 4 spline coefficients
#   histology | axis   moderated F over the 4 stage coefficients
#
# The spline df is matched to the 4 df that five fibrosis levels spend, so neither
# side wins on flexibility alone.
#
# Nominal BH is not trusted here. The axis is a continuous function of the same
# expression matrix, so some association is guaranteed by construction. The
# calibrated threshold comes from permuting the axis WITHIN cohort, which
# preserves both the axis marginal and the cohort structure and destroys only the
# donor-to-axis assignment.

suppressPackageStartupMessages({
  library(data.table); library(limma); library(edgeR); library(splines)
})
source(file.path(Sys.getenv("MASLD_PROJECT_ROOT"),
                 "scripts/analysis/continuous_axis_benchmark/lib_axis_common.R"))

OUT <- Sys.getenv("CAB_OUT_ROOT", "")
ARM <- Sys.getenv("CAB_ARM", "A1")
assert_true(nzchar(OUT), "CAB_OUT_ROOT is unset")
dir.create(file.path(OUT, "q1"), recursive = TRUE, showWarnings = FALSE)

pre <- read_prespec()
set.seed(pre$seeds$perm_base)
N_PERM <- pre$q1_permutation$n_perm
SPLINE_DF <- 4L

meta <- load_manifest()
dge <- readRDS(substrate_path("dge"))[, meta$sample_id]

axis <- fread(file.path(OUT, "arms", sprintf("%s_axis.tsv", ARM)))
assert_true("axis_raw" %in% names(axis), "Arm axis file has no axis_raw column")

# Anti-circularity restriction. A6 exists precisely so the tested genes are ones
# the axis never saw. Without this the test is close to tautological: the axis is a
# linear combination of the genes, so it explains them. The first run of this
# script omitted the restriction and reported 16,509 of 23,370 genes -- 71% of the
# transcriptome -- which is the signature of measuring circularity, not biology.
#
# The split is regenerated deterministically from the same seed and the same PCA
# gene list 02f used, so H1 here is byte-identical to the H1 that built the axis.
test_genes <- rownames(dge)
if (ARM == "A6") {
  sz <- fread(file.path(OUT, "arms", "A6_split_sizes.tsv"))
  medoid <- fread(file.path(OUT, "arms", "A6_axis_summary.tsv"))$medoid_split[1]
  seed <- sz[split == medoid, seed]
  assert_true(length(seed) == 1L, "Could not resolve the medoid split seed")
  pca_genes <- fread(file.path(OUT, "arms", "pca_input_genes.tsv"))$gene_id
  set.seed(seed)
  h1 <- sample(pca_genes, floor(length(pca_genes) / 2))
  assert_true(length(h1) == sz[split == medoid, n_axis_genes],
              "Regenerated H1 does not match the recorded axis-gene count; the split is not reproducible")
  test_genes <- setdiff(rownames(dge), h1)
  log_step("A6 anti-circularity: testing ", length(test_genes),
           " held-out genes; ", length(h1), " axis genes excluded")
}

popC <- population(meta, "POP-C")
info <- merge(popC, axis[, .(sample_id, axis_raw)], by = "sample_id")
info <- info[is.finite(axis_raw)]
info[, dataset := droplevels(factor(dataset))]
info[, stage := factor(fibrosis_stage)]
info[, axis_rank := rank_within(axis_raw, dataset)]
log_step("arm ", ARM, ": POP-C n=", nrow(info), " across ", uniqueN(info$dataset), " cohorts")
assert_true(nrow(info) >= 300L, "POP-C shrank unexpectedly after joining the axis")

sub <- dge[test_genes, info$sample_id]
BH_FAMILY <- length(test_genes)

build_design <- function(d) {
  X <- model.matrix(~ dataset + inferred_sex + stage + ns(axis_rank, df = SPLINE_DF), data = d)
  assert_true(qr(X)$rank == ncol(X), "The AH design is rank deficient")
  X
}
X <- build_design(info)
axis_cols <- grep("^ns\\(axis_rank", colnames(X))
stage_cols <- grep("^stage", colnames(X))
assert_true(length(axis_cols) == SPLINE_DF, "Unexpected number of spline coefficients")
assert_true(length(stage_cols) >= 1L, "No stage coefficients are estimable")

log_step("voomWithQualityWeights on the AH design")
v <- voomWithQualityWeights(sub, X, plot = FALSE)

fit_F <- function(design, cols) {
  fit <- eBayes(lmFit(v$E, design, weights = v$weights))
  tt <- topTable(fit, coef = cols, number = Inf, sort.by = "none")
  list(F = tt$F, P = tt$P.Value)
}
obs_axis <- fit_F(X, axis_cols)
obs_stage <- fit_F(X, stage_cols)
# The BH family is the tested set. For every arm but A6 that is the full 23,370;
# for A6 it is the held-out half, and correcting over 23,370 there would be wrong.
assert_true(length(obs_axis$P) == BH_FAMILY,
            sprintf("BH family is %d, expected %d", length(obs_axis$P), BH_FAMILY))

# --------------------------------------------------------------- permutation --
# voom weights are held at the observed fit: they depend on the mean-variance
# trend, which the axis permutation does not change. Only the linear model is
# refit per permutation. This is what makes 200 permutations affordable.
log_step("permuting the axis within cohort, ", N_PERM, " draws")
perm_max <- numeric(N_PERM)
perm_counts <- integer(N_PERM)
nominal_thresh <- 0.05
for (b in seq_len(N_PERM)) {
  d <- copy(info)
  d[, axis_rank := unlist(lapply(split(axis_rank, dataset), sample), use.names = FALSE)[
      order(order(dataset))]]
  Xp <- model.matrix(~ dataset + inferred_sex + stage + ns(axis_rank, df = SPLINE_DF), data = d)
  if (qr(Xp)$rank < ncol(Xp)) { perm_max[b] <- NA_real_; perm_counts[b] <- NA_integer_; next }
  fp <- fit_F(Xp, grep("^ns\\(axis_rank", colnames(Xp)))
  perm_max[b] <- min(fp$P, na.rm = TRUE)
  perm_counts[b] <- sum(p.adjust(fp$P, "BH") < nominal_thresh, na.rm = TRUE)
  if (b %% 25 == 0) log_step("  perm ", b, "/", N_PERM)
}

# Calibrated threshold: the p-value cut whose median permutation yield is <= 5% of
# the observed yield. Reported alongside the nominal BH count so the inflation is
# visible rather than silently corrected away.
obs_bh <- p.adjust(obs_axis$P, "BH")
n_axis_nominal <- sum(obs_bh < 0.05, na.rm = TRUE)
perm_median <- median(perm_counts, na.rm = TRUE)
emp_p <- (1 + sum(perm_counts >= n_axis_nominal, na.rm = TRUE)) / (1 + sum(!is.na(perm_counts)))
n_axis_calibrated <- max(0L, n_axis_nominal - as.integer(ceiling(perm_median)))

stage_bh <- p.adjust(obs_stage$P, "BH")
n_stage <- sum(stage_bh < 0.05, na.rm = TRUE)

res <- data.table(
  arm = ARM, population = "POP-C", n_donors = nrow(info), n_genes = length(obs_bh),
  bh_family_size = BH_FAMILY,
  genes_excluded_as_axis_inputs = length(rownames(dge)) - BH_FAMILY,
  anti_circularity_restricted = ARM == "A6",
  spline_df = SPLINE_DF,
  n_axis_given_histology_nominalBH = n_axis_nominal,
  n_axis_given_histology_calibrated = n_axis_calibrated,
  perm_median_false_yield = perm_median,
  perm_q95_false_yield = quantile(perm_counts, 0.95, na.rm = TRUE),
  empirical_p_axis = emp_p,
  n_histology_given_axis = n_stage,
  asymmetry_ratio_axis_over_histology = n_axis_calibrated / max(1L, n_stage),
  gate_g1_threshold = pre$q1_gate$an_arm_passes_only_if_all_four_hold$g1_de_increment$pop_c_min_genes,
  gate_g1_pass = n_axis_calibrated >=
    pre$q1_gate$an_arm_passes_only_if_all_four_hold$g1_de_increment$pop_c_min_genes
)
write_tsv_once(res, file.path(OUT, "q1", sprintf("%s_de_increment.tsv", ARM)))
write_tsv_once(data.table(perm = seq_len(N_PERM), n_bh05 = perm_counts, min_p = perm_max),
               file.path(OUT, "q1", sprintf("%s_de_increment_permutations.tsv", ARM)))
write_tsv_once(data.table(gene_id = rownames(v$E),
                          F_axis = obs_axis$F, P_axis = obs_axis$P, BH_axis = obs_bh,
                          F_stage = obs_stage$F, P_stage = obs_stage$P, BH_stage = stage_bh),
               file.path(OUT, "q1", sprintf("%s_de_increment_per_gene.tsv.gz", ARM)))
print(res)
log_step("Q1_DE_INCREMENT_COMPLETE arm=", ARM)
