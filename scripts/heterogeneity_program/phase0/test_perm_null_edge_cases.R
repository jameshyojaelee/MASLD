#!/usr/bin/env Rscript
root <- Sys.getenv('MASLD_PROJECT_ROOT', '/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
path <- Sys.getenv('MASLD_PERM_SOURCE', file.path(root, 'scripts/heterogeneity_program/lib/perm_null.R'))
source(path)
before <- new.env()
source(file.path(root, 'docs/technical/code_logic_remaining_20261022/before/scripts/heterogeneity_program/lib/perm_null.R'), local = before)
set.seed(91)
X <- matrix(rnorm(13 * 12), 13, 12)
group <- factor(rep(c('a', 'b'), 6))
strata <- rep(seq_len(3), each = 4)
stat <- function(g) rowMeans(X[, g == 'b', drop = FALSE]) - rowMeans(X[, g == 'a', drop = FALSE])

# Independently build the shuffled sample-index assignments and fixed tail counts.
set.seed(23)
null <- vapply(seq_len(19), function(i) {
  g <- group
  for (s in unique(strata)) {
    index <- which(strata == s)
    g[index] <- group[index[sample.int(length(index))]]
  }
  stat(g)
}, numeric(nrow(X)))
for (alternative in c('greater', 'less', 'two.sided')) {
  observed <- stat(group)
  exceeds <- switch(alternative, greater = null >= observed, less = null <= observed,
                    two.sided = abs(null) >= abs(observed))
  expected <- (1 + rowSums(exceeds)) / 20
  for (cores in c(1L, 2L)) for (keep_null in c(FALSE, TRUE)) {
    actual <- perm_null(stat, group, strata, 19L, 23L, alternative, cores, keep_null)
    original <- before$perm_null(stat, group, strata, 19L, 23L, alternative, cores, keep_null)
    stopifnot(identical(actual, original), identical(actual$p_perm, expected),
              identical(actual$q_perm, p.adjust(expected, 'BH')))
    if (keep_null) stopifnot(identical(actual$null, null))
  }
}
cat('FINITE_PERMUTATION_REPLAY_PASS: all alternatives, serial/parallel, retained/omitted null\n')

missing_observed <- function(g) c(NA_real_, NaN, Inf, -Inf, mean(g == 'b'))
answer <- perm_null(missing_observed, group, strata, n_perm = 19L)
stopifnot(all(is.na(answer$p_perm[1:4])), all(is.na(answer$q_perm[1:4])), answer$p_perm[5] == 1)
missing_draw <- function(g) c(if (identical(g, group)) 1 else NA_real_, mean(g == 'b'))
answer <- perm_null(missing_draw, group, strata, n_perm = 19L, return_null = TRUE)
stopifnot(is.na(answer$p_perm[1]), is.na(answer$q_perm[1]), answer$p_perm[2] == 1)

must_fail <- function(expression) {
  stopifnot(inherits(tryCatch({force(expression); NULL}, error = identity), 'error'))
}
for (cores in c(1L, 2L)) {
  malformed <- function(g) if (identical(g, group)) c(1, 2) else 1
  must_fail(suppressWarnings(perm_null(malformed, group, strata, n_perm = 19L, n_cores = cores)))
}
must_fail(perm_null(stat, group, strata, n_perm = 2.5))
must_fail(perm_null(stat, group, strata, n_perm = Inf))
must_fail(perm_null(stat, group, replace(strata, 1, NA), n_perm = 3L))

# Serial streaming must preserve RNG use even if the statistic itself draws numbers.
random_stat <- function(g) c(mean(g == 'b'), runif(3))
old <- before$perm_null(random_stat, group, strata, n_perm = 9L, seed = 81L)
old_state <- .Random.seed
new <- perm_null(random_stat, group, strata, n_perm = 9L, seed = 81L)
stopifnot(identical(old, new), identical(old_state, .Random.seed))
cat('UNDEFINED_PERMUTATION_REGRESSION_PASS: missing tests stay missing; malformed draws are rejected\n')
