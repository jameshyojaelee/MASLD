#!/usr/bin/env Rscript
# Arm A1: the published ordering, reimplemented as faithfully as the Methods allow.
#
# Kamzolas et al. state: quantile normalisation, ComBat with sex as covariate,
# prcomp, then "pseudo-temporal ordering was inferred using Slingshot v2.18.0",
# with the resulting order "aligned with PC1", low = early.
#
# They state NONE of: the clustering Slingshot was given, how many reduced
# dimensions, the start cluster, the distance metric, or any seed. Slingshot
# requires a clustering, so those are unavoidable choices. Rather than pick one
# silently, this arm runs a 27-setting grid and reports how much the axis moves.
# If the axis is not stable across reasonable settings then it is not identified
# from the published methods, and that partially answers Q1 before any
# downstream test runs.

suppressPackageStartupMessages({
  library(data.table)
  library(slingshot)
  library(SingleCellExperiment)
})

source(file.path(Sys.getenv("MASLD_PROJECT_ROOT"),
                 "scripts/analysis/continuous_axis_benchmark/lib_axis_common.R"))

OUT <- Sys.getenv("CAB_OUT_ROOT", "")
assert_true(nzchar(OUT), "CAB_OUT_ROOT is unset")
MATRIX_ID <- Sys.getenv("CAB_MATRIX", "E_cb")   # A1 = E_cb, A2 = E_qn, A2b = E_cbd
ARM_ID <- Sys.getenv("CAB_ARM", "A1")
# E_cbd was added 2026-08-17 for the post-hoc arm A2b. This widens an input
# whitelist only; A1 and A2 pass the same assertion before and after, so no
# completed result changes. The pre-edit copy is sealed under
# RNA-seq/results/continuous_axis_benchmark/20260814T184338Z/scripts_asrun_20260816/.
assert_true(MATRIX_ID %in% c("E_cb", "E_qn", "E_raw", "E_cbd"),
            "CAB_MATRIX must be E_cb, E_qn, E_raw or E_cbd")

pre <- read_prespec()
meta <- load_manifest()
E <- readRDS(file.path(OUT, "arms", paste0(MATRIX_ID, ".rds")))
pca_genes <- fread(file.path(OUT, "arms", "pca_input_genes.tsv"))$gene_id
assert_true(identical(colnames(E), meta$sample_id), "Matrix columns are not in manifest order")

# Label-blind alternative to a histology-chosen root: project each sample onto
# the frozen pooled disease logFC vector. The vector comes from a prior frozen
# analysis, so no histology label of these samples is consulted at fit time.
disease_projection <- function() {
  path <- file.path(project_root(),
                    "RNA-seq/results/manuscript_release/candidates",
                    "resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BG001-DECISION",
                    "arms/F_five/results/integration/deg_results.csv")
  if (!file.exists(path)) return(NULL)
  deg <- fread(path, select = c("gene", "logFC"))
  shared <- intersect(deg$gene, rownames(E))
  if (length(shared) < 1000L) return(NULL)
  w <- deg[match(shared, gene), logFC]
  as.numeric(scale(t(E[shared, ])) %*% w)
}
proj <- disease_projection()

# One axis for one (k, nPC, seed, root rule) setting.
fit_axis <- function(k, npc, seed, root_rule) {
  set.seed(seed)
  p <- prcomp(t(E[pca_genes, ]), center = TRUE, scale. = FALSE)
  rd <- p$x[, seq_len(npc), drop = FALSE]
  km <- kmeans(rd, centers = k, nstart = 25, iter.max = 100)
  cl <- km$cluster
  score <- if (root_rule == "histology") {
    st <- meta$fibrosis_stage
    vapply(split(st, cl), function(v) mean(v, na.rm = TRUE), numeric(1))
  } else {
    assert_true(!is.null(proj), "Disease projection unavailable for the label-blind root rule")
    vapply(split(proj, cl), mean, numeric(1))
  }
  score[!is.finite(score)] <- Inf
  start <- names(score)[which.min(score)]
  sce <- SingleCellExperiment(assays = list(logcounts = t(rd)),
                              reducedDims = list(PCA = rd))
  # Slingshot's default inter-cluster distance uses each cluster's covariance, so a
  # cluster with fewer members than dimensions yields NaN MST weights and the fit
  # dies. That is not a bug to route around silently: a setting that cannot be fit
  # at all is evidence about identifiability. Try the published default, fall back
  # to centroid distance, and if both fail record the setting as unfittable.
  sce <- tryCatch(
    slingshot(sce, clusterLabels = cl, reducedDim = "PCA", start.clus = start),
    error = function(e1) tryCatch(
      slingshot(sce, clusterLabels = cl, reducedDim = "PCA", start.clus = start,
                dist.method = "simple"),
      error = function(e2) NULL))
  if (is.null(sce)) {
    return(list(pt = rep(NA_real_, nrow(rd)), pt_l1 = rep(NA_real_, nrow(rd)),
                n_lineages = NA_integer_, n_scored = 0L, n_scored_l1 = 0L,
                pc1 = p$x[, 1], pc1_var = 100 * p$sdev[1]^2 / sum(p$sdev^2),
                start = start, fitted = FALSE,
                min_cluster_n = min(table(cl))))
  }

  # Slingshot fits a TREE, not necessarily one path. With k clusters it can return
  # several lineages, and slingPseudotime() is NA for any sample not assigned to a
  # given lineage. Taking lineage 1 alone silently dropped 385 of 844 donors in the
  # first run. The scalar axis is therefore the curve-weighted average across all
  # lineages, which is defined for every donor; lineage 1 is retained separately so
  # the choice stays visible.
  ptm <- slingPseudotime(sce)
  w <- slingCurveWeights(sce)
  ptm[!is.finite(ptm)] <- NA_real_
  w[!is.finite(w)] <- 0
  w[is.na(ptm)] <- 0
  denom <- rowSums(w)
  pt <- ifelse(denom > 0, rowSums(ptm * w, na.rm = TRUE) / denom, NA_real_)
  pt_l1 <- ptm[, 1]
  n_lineages <- ncol(ptm)

  # Sign convention: low = early. Anchor on PC1 and on the root cluster, never on
  # the histology of individual samples.
  if (cor(pt, p$x[, 1], use = "complete.obs") < 0) { pt <- -pt; pt_l1 <- -pt_l1 }
  if (mean(pt[cl == start], na.rm = TRUE) > mean(pt, na.rm = TRUE)) { pt <- -pt; pt_l1 <- -pt_l1 }
  list(pt = pt, pt_l1 = pt_l1, n_lineages = n_lineages,
       n_scored = sum(is.finite(pt)), n_scored_l1 = sum(is.finite(pt_l1)),
       pc1 = p$x[, 1], pc1_var = 100 * p$sdev[1]^2 / sum(p$sdev^2), start = start,
       fitted = TRUE, min_cluster_n = min(table(cl)))
}

log_step("arm ", ARM_ID, " on ", MATRIX_ID, ": primary axis (k=8, nPC=10, seed=42)")
primary <- fit_axis(8L, 10L, pre$seeds$master, "histology")
primary_blind <- fit_axis(8L, 10L, pre$seeds$master, "projection")

axis_dt <- data.table(
  sample_id = meta$sample_id, dataset = meta$dataset,
  fibrosis_stage = meta$fibrosis_stage, nas_score = meta$nas_score,
  axis_raw = primary$pt,
  axis_rank_within_cohort = rank_within(primary$pt, meta$dataset),
  axis_lineage1_only = primary$pt_l1,
  axis_raw_blind_root = primary_blind$pt,
  pc1 = primary$pc1, arm = ARM_ID, matrix_id = MATRIX_ID
)
write_tsv_once(axis_dt, file.path(OUT, "arms", sprintf("%s_axis.tsv", ARM_ID)))

# ------------------------------------------------------- identifiability grid --
log_step("identifiability grid: 27 settings")
grid <- CJ(k = c(5L, 8L, 12L), npc = c(5L, 10L, 20L), seed = c(42L, 43L, 44L))
axes <- vector("list", nrow(grid))
grid[, `:=`(n_lineages = NA_integer_, n_scored = NA_integer_,
            start_cluster = NA_character_, fitted = NA, min_cluster_n = NA_integer_)]
for (i in seq_len(nrow(grid))) {
  a <- fit_axis(grid$k[i], grid$npc[i], grid$seed[i], "histology")
  axes[[i]] <- a$pt
  set(grid, i, "n_lineages", a$n_lineages)
  set(grid, i, "n_scored", a$n_scored)
  set(grid, i, "start_cluster", as.character(a$start))
  set(grid, i, "fitted", a$fitted)
  set(grid, i, "min_cluster_n", as.integer(a$min_cluster_n))
  log_step(sprintf("  k=%2d nPC=%2d seed=%d start=%s fitted=%s lineages=%s scored=%d minclust=%d PC1var=%.2f%%",
                   grid$k[i], grid$npc[i], grid$seed[i], a$start, a$fitted,
                   ifelse(is.na(a$n_lineages), "NA", a$n_lineages),
                   a$n_scored, a$min_cluster_n, a$pc1_var))
}
write_tsv_once(grid, file.path(OUT, "arms", sprintf("%s_grid_settings.tsv", ARM_ID)))
A <- do.call(cbind, axes)
colnames(A) <- sprintf("k%d_npc%d_seed%d", grid$k, grid$npc, grid$seed)
A <- A[, vapply(seq_len(ncol(A)), function(j) any(is.finite(A[, j])), logical(1)), drop = FALSE]
S <- cor(A, method = "spearman", use = "pairwise.complete.obs")
offdiag <- S[upper.tri(S)]
offdiag <- offdiag[is.finite(offdiag)]

grid_dt <- as.data.table(S, keep.rownames = "setting")
write_tsv_once(grid_dt, file.path(OUT, "arms", sprintf("%s_axis_stability.tsv", ARM_ID)))

thr <- pre$q1_gate$an_arm_passes_only_if_all_four_hold$g4_axis_identifiability$min_median_abs_rho
# The prespecified g4 criterion is a MEDIAN, which is the wrong summary if the
# distribution is bimodal: a median of 0.95 can coexist with a quarter of settings
# producing an unrelated ordering. The criterion is applied exactly as written, and
# the quantiles below are reported alongside so the pass is never read as
# "the axis is stable" when it means "the axis is stable most of the time".
summary_dt <- data.table(
  arm = ARM_ID, matrix_id = MATRIX_ID, n_settings = ncol(A),
  median_abs_spearman = median(abs(offdiag)),
  min_abs_spearman = min(abs(offdiag)),
  p05_abs_spearman = quantile(abs(offdiag), 0.05),
  p10_abs_spearman = quantile(abs(offdiag), 0.10),
  q25_abs_spearman = quantile(abs(offdiag), 0.25),
  frac_pairs_below_0.3 = mean(abs(offdiag) < 0.3),
  frac_pairs_below_0.5 = mean(abs(offdiag) < 0.5),
  frac_pairs_above_0.8 = mean(abs(offdiag) >= 0.8),
  n_lineages_primary = primary$n_lineages,
  n_scored_primary = primary$n_scored,
  n_scored_lineage1_only = primary$n_scored_l1,
  n_settings_unfittable = sum(!grid$fitted),
  median_n_lineages_grid = median(grid$n_lineages, na.rm = TRUE),
  min_n_lineages_grid = min(grid$n_lineages, na.rm = TRUE),
  max_n_lineages_grid = max(grid$n_lineages, na.rm = TRUE),
  n_distinct_lineage_counts = uniqueN(na.omit(grid$n_lineages)),
  pc1_pct_variance = primary$pc1_var,
  kamzolas_reported_pc1_pct_variance = 9,
  spearman_primary_vs_blind_root = cor(primary$pt, primary_blind$pt, method = "spearman"),
  g4_threshold = thr,
  g4_pass = median(abs(offdiag)) >= thr
)
write_tsv_once(summary_dt, file.path(OUT, "arms", sprintf("%s_axis_summary.tsv", ARM_ID)))
print(summary_dt)
log_step("ARM_", ARM_ID, "_COMPLETE")
