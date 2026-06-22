# ─────────────────────────────────────────────────────────────────────────────
# Phase 0b — Genome-wide, cohort-blocked permutation-null engine (SHARED)
#
# The single calibration backbone for the heterogeneity-aware program (T1 DV
# selection, T2 state-discovery, T3 program convergence). It builds the null of
# ANY per-gene decision statistic by permuting the group label (e.g. stage)
# WITHIN strata (e.g. dataset), so cohort structure / composition is preserved
# and only the group↔statistic association is destroyed (GATE-B).
#
# CRITICAL (GATE-B): the statistic function `stat_fn` must recompute EVERY
# data-dependent transform (voom floor, dispersion trend, loess detrend,
# cameraPR VIF, batch term) from the permuted labels — i.e. transforms are
# refit inside each replicate. The engine only supplies permuted label vectors;
# correctness of refit-per-replicate is the caller's contract.
#
# p-values use the add-one estimator (Phipson & Smyth 2010, "Permutation
# p-values should never be zero"): p = (1 + #{null exceeds obs}) / (n_perm + 1).
#
# Usage:
#   src <- new.env(); ... ; stat_fn <- function(group) { <recompute per-gene stat> }
#   res <- perm_null(stat_fn, group = meta$stage, strata = meta$dataset,
#                    n_perm = 1000, alternative = "greater", n_cores = 16)
#   res$p_perm   # per-gene permutation p-value (add-one)
#   res$q_perm   # BH across genes
# ─────────────────────────────────────────────────────────────────────────────

# Shuffle `group` independently within each level of `strata` (cohort block).
within_stratum_permute <- function(group, strata) {
  g <- group
  for (s in unique(strata)) {
    idx <- which(strata == s)
    if (length(idx) > 1L) g[idx] <- sample(group[idx])
  }
  g
}

#' Cohort-blocked permutation null for a per-gene statistic.
#'
#' @param stat_fn  function(group) -> numeric vector (length = n_genes). MUST
#'                 recompute all data-dependent transforms internally (GATE-B).
#' @param group    observed group labels (length n_samples).
#' @param strata   blocking variable (e.g. dataset); permutation is within-block.
#' @param n_perm   number of permutations (>= 200 for a usable tail; 1000 default).
#' @param seed     RNG seed (reproducible permutation set).
#' @param alternative "greater" (DV: larger = more), "less", or "two.sided".
#' @param n_cores  parallel workers over permutations (mclapply).
#' @param return_null keep the full n_genes x n_perm null matrix (memory-heavy).
#' @return list(observed, p_perm, q_perm, exceed_count, n_perm, alternative[, null])
perm_null <- function(stat_fn, group, strata, n_perm = 1000L, seed = 42L,
                      alternative = c("greater", "two.sided", "less"),
                      n_cores = 1L, return_null = FALSE) {
  alternative <- match.arg(alternative)
  stopifnot(length(group) == length(strata), n_perm >= 1L)

  set.seed(seed)
  obs <- as.numeric(stat_fn(group))
  G <- length(obs)
  if (G == 0L) stop("stat_fn returned length-0 vector")

  # pre-generate the (reproducible) permuted label sets before any parallelism
  perms <- lapply(seq_len(n_perm), function(i) within_stratum_permute(group, strata))

  # single pass: collect the full permuted statistic per replicate (enables both
  # the per-gene exceedance count AND a decile-pooled null without recomputation)
  run_one <- function(g) as.numeric(stat_fn(g))
  if (n_cores > 1L && requireNamespace("parallel", quietly = TRUE)) {
    res <- parallel::mclapply(perms, run_one, mc.cores = n_cores)
    if (any(vapply(res, function(x) inherits(x, "try-error") || length(x) != G, logical(1))))
      stop("a permutation worker failed or returned wrong length")
  } else {
    res <- lapply(perms, run_one)
  }
  null_mat <- do.call(cbind, res)                    # G x n_perm

  exceed_count <- switch(alternative,
    greater   = rowSums(null_mat >= obs,            na.rm = TRUE),
    less      = rowSums(null_mat <= obs,            na.rm = TRUE),
    two.sided = rowSums(abs(null_mat) >= abs(obs),  na.rm = TRUE))
  p_perm <- (exceed_count + 1) / (n_perm + 1)        # add-one (Phipson & Smyth 2010)
  q_perm <- p.adjust(p_perm, method = "BH")

  out <- list(observed = obs, p_perm = p_perm, q_perm = q_perm,
              exceed_count = exceed_count, n_perm = n_perm, alternative = alternative)
  if (return_null) out$null <- null_mat
  out
}
