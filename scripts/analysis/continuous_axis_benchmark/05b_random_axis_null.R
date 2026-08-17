#!/usr/bin/env Rscript
# 05b: is the inferred axis special, or does ANY expression-derived direction
# explain this many genes?
#
# Why this script exists. The prespecification claimed arm A6 -- axis built on one
# random half of genes, tested on the other -- was "the only arm immune to the
# PC1-axis / PC1-genes circularity". That claim was WRONG, and the data said so:
# restricting to 11,685 held-out genes still returned 8,303 hits (71%), barely
# below the un-restricted arms. Gene splitting does not create independence,
# because genes are co-expressed. Any random half of the transcriptome spans the
# same dominant expression components as the other half, so an axis built on H1
# predicts H2 nearly as well as it predicts H1.
#
# The right control is therefore not a held-out gene set but a held-out AXIS. This
# is the same logic that retired the staging classifier in this repo: the question
# was never "does the curated set beat zero", it was "does it beat 500 random
# genes". Here: does the published construction beat an axis built from a random
# subset of genes, scored on identical test genes with identical machinery?
#
# Design. Hold out a fixed test set T (25% of genes) from ALL axis construction.
# Build the real axis and N random axes from the remaining 75%. Score every axis
# on T with the same AH model. Where the real axis falls in the random
# distribution is the answer.

suppressPackageStartupMessages({
  library(data.table); library(limma); library(edgeR); library(splines)
  library(slingshot); library(SingleCellExperiment)
})
source(file.path(Sys.getenv("MASLD_PROJECT_ROOT"),
                 "scripts/analysis/continuous_axis_benchmark/lib_axis_common.R"))

OUT <- Sys.getenv("CAB_OUT_ROOT", "")
MATRIX_ID <- Sys.getenv("CAB_MATRIX", "E_cb")
N_RANDOM <- as.integer(Sys.getenv("CAB_N_RANDOM", "30"))
assert_true(nzchar(OUT), "CAB_OUT_ROOT is unset")
dir.create(file.path(OUT, "q1"), recursive = TRUE, showWarnings = FALSE)

pre <- read_prespec()
set.seed(pre$seeds$split_base)
SPLINE_DF <- 4L

meta <- load_manifest()
dge <- readRDS(substrate_path("dge"))[, meta$sample_id]
E <- readRDS(file.path(OUT, "arms", paste0(MATRIX_ID, ".rds")))
pca_genes <- fread(file.path(OUT, "arms", "pca_input_genes.tsv"))$gene_id

all_genes <- rownames(dge)
test_genes <- sample(all_genes, floor(length(all_genes) * 0.25))
build_pool <- setdiff(pca_genes, test_genes)
log_step("test genes (never used to build any axis): ", length(test_genes))
log_step("axis construction pool: ", length(build_pool))

popC <- population(meta, "POP-C")
info <- copy(popC)
info[, dataset := droplevels(factor(dataset))]
info[, stage := factor(fibrosis_stage)]

# One voom fit under the histology-only design. Weights come from the
# mean-variance trend, which does not depend on which axis is added, so holding
# them fixed lets 31 axes be scored without 31 voom runs. Identical treatment for
# the real and random axes, so the comparison is unaffected.
X0 <- model.matrix(~ dataset + inferred_sex + stage, data = info)
v <- voomWithQualityWeights(dge[test_genes, info$sample_id], X0, plot = FALSE)

axis_from <- function(genes, seed) {
  set.seed(seed)
  p <- prcomp(t(E[genes, info$sample_id, drop = FALSE]), center = TRUE, scale. = FALSE)
  rd <- p$x[, 1:10, drop = FALSE]
  km <- kmeans(rd, centers = 8L, nstart = 10, iter.max = 100)
  st <- vapply(split(info$fibrosis_stage, km$cluster),
               function(z) mean(z, na.rm = TRUE), numeric(1))
  st[!is.finite(st)] <- Inf
  sce <- SingleCellExperiment(assays = list(logcounts = t(rd)),
                              reducedDims = list(PCA = rd))
  sce <- tryCatch(slingshot(sce, clusterLabels = km$cluster, reducedDim = "PCA",
                            start.clus = names(st)[which.min(st)]),
                  error = function(e) NULL)
  if (is.null(sce)) return(NULL)
  ptm <- slingPseudotime(sce); w <- slingCurveWeights(sce)
  ptm[!is.finite(ptm)] <- NA_real_; w[!is.finite(w)] <- 0; w[is.na(ptm)] <- 0
  den <- rowSums(w)
  pt <- ifelse(den > 0, rowSums(ptm * w, na.rm = TRUE) / den, NA_real_)
  if (cor(pt, p$x[, 1], use = "complete.obs") < 0) pt <- -pt
  pt
}

score_axis <- function(pt) {
  d <- copy(info); d[, axis_rank := rank_within(pt, dataset)]
  d <- d[is.finite(axis_rank)]
  if (nrow(d) < 200L) return(NULL)
  X <- model.matrix(~ dataset + inferred_sex + stage + ns(axis_rank, df = SPLINE_DF), data = d)
  if (qr(X)$rank < ncol(X)) return(NULL)
  idx <- match(d$sample_id, info$sample_id)
  fit <- eBayes(lmFit(v$E[, idx], X, weights = v$weights[, idx]))
  tt <- topTable(fit, coef = grep("^ns\\(axis_rank", colnames(X)),
                 number = Inf, sort.by = "none")
  list(n_hits = sum(p.adjust(tt$P.Value, "BH") < 0.05, na.rm = TRUE),
       rho_fib = suppressWarnings(cor(pt, info$fibrosis_stage,
                                      method = "spearman", use = "complete.obs")))
}

log_step("real axis from the construction pool")
real_pt <- axis_from(build_pool, pre$seeds$master)
assert_true(!is.null(real_pt), "The real axis could not be fit")
real <- score_axis(real_pt)
log_step("real axis: ", real$n_hits, " of ", length(test_genes),
         " test genes; rho_fib=", round(real$rho_fib, 4))

log_step(N_RANDOM, " random axes from random subsets of the same pool")
res <- vector("list", N_RANDOM)
for (b in seq_len(N_RANDOM)) {
  set.seed(pre$seeds$sim_base + 1000L + b)
  g <- sample(build_pool, length(build_pool) %/% 2L)
  pt <- axis_from(g, pre$seeds$sim_base + 1000L + b)
  if (is.null(pt)) { res[[b]] <- NULL; next }
  s <- score_axis(pt)
  if (is.null(s)) { res[[b]] <- NULL; next }
  res[[b]] <- data.table(draw = b, n_hits = s$n_hits, rho_fib = s$rho_fib)
  if (b %% 5 == 0) log_step("  random axis ", b, "/", N_RANDOM, ": ", s$n_hits, " hits")
}
rnd <- rbindlist(res)
assert_true(nrow(rnd) >= 10L, "Too few random axes fit successfully")

emp_p <- (1 + sum(rnd$n_hits >= real$n_hits)) / (1 + nrow(rnd))
verdict <- data.table(
  matrix_id = MATRIX_ID, n_test_genes = length(test_genes),
  n_build_pool = length(build_pool), n_random_axes = nrow(rnd),
  real_n_hits = real$n_hits, real_frac = real$n_hits / length(test_genes),
  real_rho_fibrosis = real$rho_fib,
  random_median_hits = median(rnd$n_hits), random_q95_hits = quantile(rnd$n_hits, 0.95),
  random_min_hits = min(rnd$n_hits), random_max_hits = max(rnd$n_hits),
  random_median_rho_fibrosis = median(rnd$rho_fib, na.rm = TRUE),
  empirical_p_real_vs_random = emp_p,
  real_beats_random_axes = real$n_hits > quantile(rnd$n_hits, 0.95),
  interpretation = if (real$n_hits > quantile(rnd$n_hits, 0.95))
    "the inferred axis explains more held-out genes than a random expression direction"
  else
    "the inferred axis is indistinguishable from a random expression direction; the DE-increment counts measure that any global expression component explains most genes, not that this axis is meaningful")
write_tsv_once(verdict, file.path(OUT, "q1", "random_axis_null_verdict.tsv"))
write_tsv_once(rnd, file.path(OUT, "q1", "random_axis_null_draws.tsv"))
write_tsv_once(data.table(gene_id = test_genes),
               file.path(OUT, "q1", "random_axis_null_test_genes.tsv.gz"))
print(verdict)
log_step("RANDOM_AXIS_NULL_COMPLETE beats_random=", verdict$real_beats_random_axes)
