# H49 shared functions (sourced by h49_score.R and test_h49_lib.R).

suppressPackageStartupMessages(library(Matrix))

cmp_key <- function(s) gsub("[^a-z0-9]", "", tolower(s))

# Signed z from a moderated t with df (vectors), keeping precision in the tails.
t_to_z <- function(t, df) {
  z <- numeric(length(t))
  pos <- !is.na(t) & t > 0
  neg <- !is.na(t) & t <= 0
  z[pos] <- qnorm(pt(t[pos], df[pos], lower.tail = FALSE), lower.tail = FALSE)
  z[neg] <- qnorm(pt(t[neg], df[neg]))
  z[is.na(t)] <- NA_real_
  pmax(pmin(z, 37), -37)
}

t_to_p <- function(t, df) 2 * pt(-abs(t), df)

# Random gene sets as a sparse averaging matrix: n_draw x n_universe, rows sum to 1.
draw_sets <- function(n_universe, size, n_draw) {
  idx <- vapply(seq_len(n_draw), function(i) sample.int(n_universe, size), integer(size))
  if (size == 1) idx <- matrix(idx, nrow = 1)
  sparseMatrix(i = rep(seq_len(n_draw), each = size), j = as.vector(idx),
               x = 1 / size, dims = c(n_draw, n_universe))
}

# Mean over set rows of a genes x conditions matrix.
set_mean <- function(mat, rows) colMeans(mat[rows, , drop = FALSE])

# Two-sided permutation p with +1 correction; null is n_draw x n_cond, obs length n_cond.
perm_p_two <- function(obs, null) {
  n <- nrow(null)
  ge <- colSums(sweep(null, 2, obs, ">="))
  le <- colSums(sweep(null, 2, obs, "<="))
  pmin(2 * pmin((1 + ge) / (n + 1), (1 + le) / (n + 1)), 1)  # first argument keeps the names
}

perm_p_greater <- function(obs, null) {
  n <- nrow(null)
  (1 + colSums(sweep(null, 2, obs, ">="))) / (n + 1)
}

col_median <- function(m) apply(m, 2, median)

# Exact two-sided permutation for difference in means between two labelled groups.
exact_perm_diff <- function(x, is_a) {
  n <- length(x); k <- sum(is_a)
  obs <- mean(x[is_a]) - mean(x[!is_a])
  cmb <- combn(n, k)
  stat <- apply(cmb, 2, function(ix) mean(x[ix]) - mean(x[-ix]))
  list(obs = obs, p = mean(abs(stat) >= abs(obs) - 1e-12), n_perm = ncol(cmb))
}
