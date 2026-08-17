#!/usr/bin/env Rscript
# Arm A6: the only arm immune to the circularity that makes every other arm's Q1b
# result suspect.
#
# The problem. Slingshot is aligned to PC1, and PC1 is a linear combination of the
# top-variance genes. Testing those same genes for association with the axis is
# guaranteed to succeed. A permuted-axis null does not fix this: permutation fixes
# the MAGNITUDE of association expected by chance, but the genes were still chosen
# because they vary along the axis.
#
# The fix. Split the genes at random into halves. Build the axis from H1 only.
# Test only H2, with the BH family restricted to H2. A gene in H2 never entered
# the axis, so any association it shows is not selection.
#
# Twenty splits, because a single split is one draw from a noisy procedure.

suppressPackageStartupMessages({
  library(data.table); library(slingshot); library(SingleCellExperiment)
})
source(file.path(Sys.getenv("MASLD_PROJECT_ROOT"),
                 "scripts/analysis/continuous_axis_benchmark/lib_axis_common.R"))

OUT <- Sys.getenv("CAB_OUT_ROOT", "")
assert_true(nzchar(OUT), "CAB_OUT_ROOT is unset")
MATRIX_ID <- Sys.getenv("CAB_MATRIX", "E_cb")
N_SPLIT <- as.integer(Sys.getenv("CAB_N_SPLIT", "20"))

pre <- read_prespec()
meta <- load_manifest()
E <- readRDS(file.path(OUT, "arms", paste0(MATRIX_ID, ".rds")))
pca_genes <- fread(file.path(OUT, "arms", "pca_input_genes.tsv"))$gene_id
assert_true(identical(colnames(E), meta$sample_id), "Matrix columns are not in manifest order")

axis_from <- function(genes, seed) {
  set.seed(seed)
  p <- prcomp(t(E[genes, , drop = FALSE]), center = TRUE, scale. = FALSE)
  rd <- p$x[, 1:10, drop = FALSE]
  km <- kmeans(rd, centers = 8L, nstart = 25, iter.max = 100)
  st <- vapply(split(meta$fibrosis_stage, km$cluster),
               function(v) mean(v, na.rm = TRUE), numeric(1))
  st[!is.finite(st)] <- Inf
  sce <- SingleCellExperiment(assays = list(logcounts = t(rd)),
                              reducedDims = list(PCA = rd))
  sce <- slingshot(sce, clusterLabels = km$cluster, reducedDim = "PCA",
                   start.clus = names(st)[which.min(st)])
  pt <- slingPseudotime(sce)[, 1]
  if (cor(pt, p$x[, 1], use = "complete.obs") < 0) pt <- -pt
  pt
}

log_step("arm A6: ", N_SPLIT, " gene splits on ", MATRIX_ID)
rows <- vector("list", N_SPLIT)
axes <- matrix(NA_real_, nrow = ncol(E), ncol = N_SPLIT,
               dimnames = list(colnames(E), paste0("split", seq_len(N_SPLIT))))
half_sets <- vector("list", N_SPLIT)

for (s in seq_len(N_SPLIT)) {
  seed <- pre$seeds$split_base + s
  set.seed(seed)
  h1 <- sample(pca_genes, floor(length(pca_genes) / 2))
  h2 <- setdiff(rownames(E), h1)          # everything the axis never saw
  pt <- axis_from(h1, seed)
  axes[, s] <- pt
  half_sets[[s]] <- data.table(split = s, seed = seed,
                               n_axis_genes = length(h1), n_test_genes = length(h2))
  log_step(sprintf("  split %2d/%d  axis genes=%d  test genes=%d", s, N_SPLIT,
                   length(h1), length(h2)))
  rows[[s]] <- data.table(split = s, sample_id = colnames(E), axis_raw = pt)
}

long <- rbindlist(rows)
write_tsv_once(long, file.path(OUT, "arms", "A6_axis_by_split.tsv.gz"))
write_tsv_once(rbindlist(half_sets), file.path(OUT, "arms", "A6_split_sizes.tsv"))

# The representative A6 axis is the split whose ordering is most typical of all
# twenty, i.e. the medoid under Spearman distance. Choosing the medoid rather than
# the best-performing split avoids selecting on the outcome.
S <- cor(axes, method = "spearman", use = "pairwise.complete.obs")
medoid <- which.max(colMeans(abs(S)))
rep_dt <- data.table(sample_id = colnames(E), dataset = meta$dataset,
                     fibrosis_stage = meta$fibrosis_stage, nas_score = meta$nas_score,
                     axis_raw = axes[, medoid])
rep_dt[, axis_rank_within_cohort := rank_within(axis_raw, dataset)]
rep_dt[, `:=`(arm = "A6", matrix_id = MATRIX_ID, split_used = medoid)]
write_tsv_once(rep_dt, file.path(OUT, "arms", "A6_axis.tsv"))

offdiag <- S[upper.tri(S)]
summary_dt <- data.table(
  arm = "A6", matrix_id = MATRIX_ID, n_splits = N_SPLIT,
  medoid_split = medoid,
  median_abs_spearman_between_splits = median(abs(offdiag)),
  min_abs_spearman_between_splits = min(abs(offdiag)),
  note = "axis built on gene half H1; all Q1b testing restricted to H2 with BH over H2 only"
)
write_tsv_once(summary_dt, file.path(OUT, "arms", "A6_axis_summary.tsv"))
write_tsv_once(as.data.table(S, keep.rownames = "split"),
               file.path(OUT, "arms", "A6_split_agreement.tsv"))
print(summary_dt)
log_step("ARM_A6_COMPLETE")
