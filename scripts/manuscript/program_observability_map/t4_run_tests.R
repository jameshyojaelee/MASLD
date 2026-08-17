#!/usr/bin/env Rscript

# Unit tests for T4 on synthetic data with known structure.
#
# The two tests that matter are the ones the design exists to pass:
#   - a planted latent factor shared by two DISJOINT programs must be recovered;
#   - two programs that share half their genes and nothing else must NOT be
#     called, even though their raw residual correlation is large.
# A third test reproduces the sibling system's retraction on purpose: a U-shaped
# stage response shared by two programs is shown to manufacture an edge under a
# linear stage term and to be removed by the saturated one.
#
# Light enough to run on the login node before any job is submitted.

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_dir <- dirname(normalizePath(script_arg))
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})
source(file.path(script_dir, "analysis_lib.R"))
source(file.path(script_dir, "t4_interaction_lib.R"))

set.seed(4104)
passed <- 0L
check <- function(label, condition, detail = "") {
  if (!isTRUE(condition)) fail("FAIL: ", label, if (nzchar(detail)) paste0(" [", detail, "]") else "")
  passed <<- passed + 1L
  cat("  ok  ", label, if (nzchar(detail)) paste0("  (", detail, ")") else "", "\n", sep = "")
}

# ------------------------------------------------------------------- fixture

cohorts <- paste0("C", 1:4)
n_per <- 100L
n_genes <- 1000L
genes <- sprintf("G%03d", seq_len(n_genes))

meta <- rbindlist(lapply(cohorts, function(cc) {
  data.table(sample_id = paste0(cc, "_", sprintf("%03d", seq_len(n_per))),
             dataset = cc,
             sex_final = sample(c("M", "F"), n_per, TRUE),
             fibrosis_stage = sample(0:4, n_per, TRUE),
             nas_score = sample(0:8, n_per, TRUE))
}))
n_total <- nrow(meta)

# Twenty-four programs. P1/P2 disjoint and coordinated; P3/P4 share 20 of 40
# genes and are otherwise unrelated; P5/P6 disjoint and both respond to stage as
# a U. The eighteen fillers exist so that the program space is large enough for
# projecting out k components to be a small perturbation, as it is at p = 113 in
# the real data, rather than half the space.
sets <- list(
  P1 = genes[1:40], P2 = genes[41:80],
  P3 = genes[81:120], P4 = genes[c(101:120, 121:140)],
  P5 = genes[141:180], P6 = genes[181:220]
)
for (k in seq_len(18L)) {
  sets[[sprintf("F%02d", k)]] <- genes[(220 + (k - 1) * 40) + seq_len(40)]
}
features <- names(sets)

latent <- rnorm(n_total)                                   # P1-P2 coordination
u_shape <- with(meta, (fibrosis_stage - 2)^2)              # P5-P6 stage response
u_shape <- as.numeric(scale(u_shape))

expr <- matrix(rnorm(n_genes * n_total), n_genes, n_total,
               dimnames = list(genes, meta$sample_id))
add <- function(mat, gs, v, amp) {
  mat[gs, ] <- mat[gs, ] + rep(amp * v, each = length(gs))
  mat
}
expr <- add(expr, sets$P1, latent, 0.55)
expr <- add(expr, sets$P2, latent, 0.55)
expr <- add(expr, sets$P5, u_shape, 0.60)
expr <- add(expr, sets$P6, u_shape, 0.60)
# P3 and P4 each get their own private factor so that they are internally
# coherent, which is what makes their shared genes able to induce a correlation.
expr <- add(expr, sets$P3, rnorm(n_total), 0.55)
expr <- add(expr, sets$P4, rnorm(n_total), 0.55)
# Each filler gets its own private factor: internally coherent, mutually
# independent, and connected to nothing.
for (k in seq_len(18L)) {
  expr <- add(expr, sets[[sprintf("F%02d", k)]], rnorm(n_total), 0.55)
}
# The shared genes carry both private factors, exactly as a real shared gene
# would; nothing links P3 to P4 beyond the genes they hold in common.

membership <- rbindlist(lapply(features, function(f) data.table(
  program_uid = f, mapped_symbol = sets[[f]],
  original_l1_weight = 1 / length(sets[[f]]))))
registry <- data.table(program_uid = features)

weights <- t4_membership_weights(membership, features)
blocks <- lapply(cohorts, function(cc)
  t4_build_block(expr, meta[dataset == cc], weights, features, 0.80))
names(blocks) <- cohorts
pairs <- t4_pair_table(features, weights)
ops <- lapply(blocks, t4_disjoint_operators, pairs = pairs, features = features,
              min_retained_weight = 0.50)

pair_of <- function(a, b) pairs[(program_a == a & program_b == b) |
                                  (program_a == b & program_b == a), pair_index]

# ------------------------------------------------ 1. the residual shortcut

# U %*% w must be the studentized residual of the program score that
# score_programs() produces, on the real scoring path, in every cohort.
scored <- score_programs(expr, meta, membership, registry, 0.80)
max_diff <- 0
for (cc in cohorts) {
  block <- blocks[[cc]]
  Y <- scored$primary[features, block$samples, drop = FALSE]
  ref <- hc3_fit_matrix(Y, block$X)$studentized_residuals      # features x samples
  mine <- t4_residual_gene_matrix(block) %*% block$W           # samples x features
  d <- max(abs(stats::cor(t(ref)) - stats::cor(mine)))
  max_diff <- max(max_diff, d)
}
check("U %*% w reproduces hc3_fit_matrix studentized-residual correlations",
      max_diff < 1e-9, sprintf("max |delta r| = %.2e", max_diff))

# ------------------------------------------------ 2. pair geometry

check("Jaccard is 0 for disjoint programs and 1/3 for the 50 percent overlap",
      pairs[pair_index == pair_of("P1", "P2"), jaccard] == 0 &&
        abs(pairs[pair_index == pair_of("P3", "P4"), jaccard] - 20 / 60) < 1e-12)
check("weight cosine matches the analytic overlap",
      abs(pairs[pair_index == pair_of("P3", "P4"), weight_cosine] - 0.5) < 1e-12,
      sprintf("%.4f", pairs[pair_index == pair_of("P3", "P4"), weight_cosine]))

# ------------------------------------------------ 3. the null preserves genes

Zp <- t4_permute_genes(blocks[[1L]]$Z)
check("per-gene permutation preserves every gene's own values",
      all(vapply(seq_len(nrow(Zp)), function(k)
        isTRUE(all.equal(unname(sort(Zp[k, ])),
                         unname(sort(blocks[[1L]]$Z[k, ])))), logical(1))))
check("per-gene permutation actually reorders donors",
      mean(vapply(seq_len(nrow(Zp)), function(k)
        identical(Zp[k, ], blocks[[1L]]$Z[k, ]), logical(1))) < 0.01)

# ------------------------------------------------ 4. observed pass and null

observed <- t4_one_pass(blocks, pairs, ops)
n_null <- 300L
null_full <- matrix(NA_real_, n_null, nrow(pairs))
null_disj <- matrix(NA_real_, n_null, nrow(pairs))
for (b in seq_len(n_null)) {
  Zl <- lapply(blocks, function(bl) t4_permute_genes(bl$Z))
  p <- t4_one_pass(blocks, pairs, ops, Zlist = Zl)
  null_full[b, ] <- p$z_full
  null_disj[b, ] <- p$z_disjoint
}
std <- function(obs, null) (obs - colMeans(null, na.rm = TRUE)) /
  apply(null, 2L, stats::sd, na.rm = TRUE)
z_full <- std(observed$z_full, null_full)
z_disj <- std(observed$z_disjoint, null_disj)

i12 <- pair_of("P1", "P2"); i34 <- pair_of("P3", "P4"); i56 <- pair_of("P5", "P6")

check("planted coordination between disjoint programs is recovered",
      z_disj[i12] > 6, sprintf("standardized z = %.1f, r = %.2f",
                               z_disj[i12], tanh(observed$z_disjoint[i12])))

# THE TEST THIS DESIGN EXISTS FOR.
check("shared membership alone produces a large raw correlation",
      tanh(observed$z_full[i34]) > 0.30,
      sprintf("raw r = %.2f", tanh(observed$z_full[i34])))
# The per-gene null centres the edge on the correlation shared membership
# produces UNDER GENE INDEPENDENCE, which is the weight cosine. When the shared
# genes are themselves coordinated, as they are here and as they are in real
# data, the induced correlation exceeds that and the full-membership arm is left
# looking significant with nothing biological behind it. This is not a bug to be
# tuned away; it is the reason the leave-shared-genes-out arm is primary and the
# full-membership arm is reported only alongside its overlap.
check("the per-gene null recovers the analytic cosine as the edge centre",
      abs(tanh(mean(null_full[, i34])) - pairs[pair_index == i34, weight_cosine]) < 0.03,
      sprintf("null mean r = %.3f vs cosine %.3f",
              tanh(mean(null_full[, i34])), pairs[pair_index == i34, weight_cosine]))
check("...but it under-corrects when the shared genes are coordinated, so the full-membership arm is inflated",
      z_full[i34] > 3,
      sprintf("standardized z = %.2f on raw r = %.2f, null mean r = %.2f",
              z_full[i34], tanh(observed$z_full[i34]), tanh(mean(null_full[, i34]))))
check("...and the leave-shared-genes-out correlation is null",
      abs(z_disj[i34]) < 3 && abs(tanh(observed$z_disjoint[i34])) < 0.12,
      sprintf("disjoint r = %.3f, standardized z = %.2f",
              tanh(observed$z_disjoint[i34]), z_disj[i34]))

# ------------------------------------------------ 5. the retraction test

# Two disjoint programs that share only a U-shaped stage response. Under the
# saturated residual model they must not be called. Under a linear stage term
# they are, which is the artifact this design is built to avoid.
linear_blocks <- lapply(cohorts, function(cc) {
  d <- meta[dataset == cc]
  Z <- zscore_rows(expr[, d$sample_id, drop = FALSE])
  variable <- rownames(Z)[rowSums(is.finite(Z)) == ncol(Z)]
  Z <- Z[variable, , drop = FALSE]
  X <- stats::model.matrix(~ sex_final + fibrosis_stage + nas_score, data = d)
  xtx_inv <- solve(crossprod(X))
  lev <- rowSums((X %*% xtx_inv) * X)
  bl <- blocks[[cc]]
  list(cohort = cc, samples = d$sample_id, genes = variable, Z = Z, X = X,
       xtx_inv = xtx_inv, dinv = 1 / sqrt(pmax(1 - lev, 1e-8)),
       n = nrow(X), n_par = ncol(X), df = nrow(X) - ncol(X), W = bl$W,
       meta = d[, .(sample_id, fibrosis_stage, nas_score, sex_final)])
})
names(linear_blocks) <- cohorts
lin_obs <- t4_one_pass(linear_blocks, pairs, ops)
lin_null <- matrix(NA_real_, 200L, nrow(pairs))
for (b in seq_len(200L)) {
  Zl <- lapply(linear_blocks, function(bl) t4_permute_genes(bl$Z))
  lin_null[b, ] <- t4_one_pass(linear_blocks, pairs, ops, Zlist = Zl)$z_disjoint
}
z_lin <- std(lin_obs$z_disjoint, lin_null)
check("a linear stage term manufactures an edge from a U-shaped stage response",
      z_lin[i56] > 5, sprintf("standardized z = %.1f, r = %.2f",
                              z_lin[i56], tanh(lin_obs$z_disjoint[i56])))
check("the saturated residual model removes it",
      abs(z_disj[i56]) < 3, sprintf("standardized z = %.2f, r = %.3f",
                                    z_disj[i56], tanh(observed$z_disjoint[i56])))
check("the planted coordination survives the saturated model",
      z_disj[i12] > 6)

# ------------------------------------------------ 6. stage-leakage diagnostic

leak_sat <- t4_stage_leakage(blocks)
leak_lin <- t4_stage_leakage(linear_blocks)
check("stage-leakage diagnostic flags the linear model and clears the saturated one",
      max(abs(leak_sat$spearman_stage)) < max(abs(leak_lin$spearman_stage)) &&
        mean(leak_sat$stage_r2) < 0.02,
      sprintf("saturated mean stage R2 = %.4f, linear = %.4f",
              mean(leak_sat$stage_r2), mean(leak_lin$stage_r2)))

# ------------------------------------------------ 7. entanglement eligibility

# A pair that cannot give up its shared genes and keep 80 percent of its weight
# must be refused rather than tested.
ops80 <- lapply(blocks, t4_disjoint_operators, pairs = pairs, features = features,
                min_retained_weight = 0.80)
check("a 50 percent-overlap pair is refused at an 80 percent retention floor",
      !ops80[[1L]]$eligible[match(i34, ops80[[1L]]$pair_index)])
check("...and is admitted at a 50 percent floor",
      ops[[1L]]$eligible[match(i34, ops[[1L]]$pair_index)])

# ------------------------------------------------ 8. global statistics

null_globals <- vapply(seq_len(50L), function(b) {
  Zl <- lapply(blocks, function(bl) t4_permute_genes(bl$Z))
  t4_one_pass(blocks, pairs, ops, Zlist = Zl)$globals
}, numeric(length(observed$globals)))
check("planted structure raises the leading eigenvalue above the null",
      observed$globals[["lambda1"]] > max(null_globals["lambda1", ]),
      sprintf("observed %.3f vs null max %.3f", observed$globals[["lambda1"]],
              max(null_globals["lambda1", ])))
check("planted structure raises cross-cohort agreement above the null",
      observed$globals[["cross_cohort_agreement"]] >
        max(null_globals["cross_cohort_agreement", ]),
      sprintf("observed %.3f vs null max %.3f",
              observed$globals[["cross_cohort_agreement"]],
              max(null_globals["cross_cohort_agreement", ])))

# ------------------------------------------------ 9. the global-component arm

# Every program score built from one bulk library shares variance that belongs
# to no pair. Added here as a factor common to all genes, it makes essentially
# every pair correlated; the partial arm has to strip it while leaving a real
# pairwise factor alone.
global_factor <- rnorm(n_total)
expr_g <- expr + rep(1.4 * global_factor, each = n_genes)
blocks_g <- lapply(cohorts, function(cc)
  t4_build_block(expr_g, meta[dataset == cc], weights, features, 0.80))
names(blocks_g) <- cohorts
ops_g <- lapply(blocks_g, t4_disjoint_operators, pairs = pairs, features = features,
                min_retained_weight = 0.50)
obs_g <- t4_one_pass(blocks_g, pairs, ops_g)
null_g_d <- matrix(NA_real_, 200L, nrow(pairs))
null_g_p <- matrix(NA_real_, 200L, nrow(pairs))
null_g_p3 <- matrix(NA_real_, 200L, nrow(pairs))
for (b in seq_len(200L)) {
  Zl <- lapply(blocks_g, function(bl) t4_permute_genes(bl$Z))
  p <- t4_one_pass(blocks_g, pairs, ops_g, Zlist = Zl)
  null_g_d[b, ] <- p$z_disjoint
  null_g_p[b, ] <- p$z_partial_k1
  null_g_p3[b, ] <- p$z_partial_k3
}
zg_d <- std(obs_g$z_disjoint, null_g_d)
zg_p <- std(obs_g$z_partial_k1, null_g_p)
zg_p3 <- std(obs_g$z_partial_k3, null_g_p3)
noise_pairs <- setdiff(seq_len(nrow(pairs)), i12)

check("a shared global component makes almost every pair correlated",
      mean(abs(zg_d[noise_pairs]) > 3) > 0.8,
      sprintf("%.0f%% of non-planted pairs called in the disjoint arm",
              100 * mean(abs(zg_d[noise_pairs]) > 3)))
# ONE COMPONENT IS NOT ALWAYS THE WHOLE GLOBAL AXIS. The leading component is
# pulled toward the programs with the most extra shared variance (here P1-P2),
# so it under-removes the global factor from the rest and a residue can survive.
# That is why k = 3 is carried as a prespecified sensitivity and why the k = 1
# count is never reported on its own.
# Calls are made the way the pipeline makes them: a positive excess surviving BH
# across the whole pair family, never an uncorrected z.
called_positive <- function(z) {
  q <- p.adjust(2 * stats::pnorm(-abs(z)), method = "BH")
  which(q < 0.05 & z > 0)
}
leak_k1 <- intersect(called_positive(zg_p), noise_pairs)
leak_k3 <- intersect(called_positive(zg_p3), noise_pairs)
check("the k = 1 partial arm removes the global component from at least 99 percent of pairs",
      length(leak_k1) / length(noise_pairs) < 0.01,
      sprintf("%d of %d non-planted pairs still called at k = 1 (%s)",
              length(leak_k1), length(noise_pairs),
              paste0(pairs$program_a[leak_k1], "-", pairs$program_b[leak_k1], collapse = ", ")))
check("the k = 3 partial arm does the same",
      length(leak_k3) / length(noise_pairs) < 0.01,
      sprintf("%d of %d non-planted pairs still called at k = 3",
              length(leak_k3), length(noise_pairs)))
# THE STATED LIMIT OF THIS ARM. The component is estimated, not known. P5 and P6
# are the only two programs here with no private factor, so after the U-shaped
# stage response is removed they are almost pure global axis; the part of the
# global factor the estimated component misses is the same in both and they stay
# correlated. Programs that are little more than a readout of the global axis
# can therefore still pair up, and the projection cannot fix that. Every other
# non-planted pair, all of which carry a private factor, is removed.
check("the residue is confined to programs that carry no private structure",
      all(pairs$program_a[leak_k1] %in% c("P5", "P6")) &&
        all(pairs$program_b[leak_k1] %in% c("P5", "P6")),
      sprintf("k = 1 residue: %s",
              paste0(pairs$program_a[leak_k1], "-", pairs$program_b[leak_k1], collapse = ", ")))
check("the planted pairwise factor is called at both k = 1 and k = 3",
      i12 %in% called_positive(zg_p) && i12 %in% called_positive(zg_p3),
      sprintf("k1 z = %.1f r = %.2f | k3 z = %.1f r = %.2f",
              zg_p[i12], tanh(obs_g$z_partial_k1[i12]),
              zg_p3[i12], tanh(obs_g$z_partial_k3[i12])))
# The leftovers come back NEGATIVE, which is the projection artifact the support
# rule refuses to call: removing a dimension from p programs forces the leftover
# covariances to sum down.
check("...and the leftovers are negative, which is why only positive partial edges are called",
      mean(zg_p[noise_pairs] < 0) > 0.8,
      sprintf("%.0f%% of non-planted partial edges are negative, median r = %.2f",
              100 * mean(zg_p[noise_pairs] < 0),
              stats::median(tanh(obs_g$z_partial_k1[noise_pairs]))))
check("the global component is measured and reported",
      obs_g$globals[["global_component_share"]] >
        observed$globals[["global_component_share"]],
      sprintf("with global factor %.3f vs without %.3f",
              obs_g$globals[["global_component_share"]],
              observed$globals[["global_component_share"]]))

# ------------------------------------------------ 10. MDE is monotone and sane

mde <- t4_minimum_detectable_excess(0, apply(null_disj, 2L, stats::sd, na.rm = TRUE),
                                    alpha = 0.05)
check("minimum detectable excess is finite and below 1 everywhere",
      all(is.finite(mde)) && all(mde > 0 & mde < 1),
      sprintf("median %.3f", stats::median(mde)))

cat("\n", passed, " checks passed\n", sep = "")
