# Unit checks for h49_lib.R. Run: Rscript --vanilla test_h49_lib.R (from the script directory).
source("h49_lib.R")
set.seed(1)

# t_to_z agrees with the direct formula where that is precise, and keeps sign
t <- c(-5, -1.2, 0, 0.7, 3, 12); df <- c(4, 10, 20, 8, 30, 15)
direct <- qnorm(pt(t, df))
stopifnot(all(abs(t_to_z(t, df)[1:5] - direct[1:5]) < 1e-8))
stopifnot(t_to_z(12, 15) > 5, sign(t_to_z(t, df))[c(1, 2, 4, 5, 6)] == sign(t)[c(1, 2, 4, 5, 6)])
stopifnot(abs(t_to_z(50, 1e6) - 37) < 1e-9 || t_to_z(50, 1e6) <= 37)

# draw_sets: rows average exactly `size` distinct genes
M <- draw_sets(100, 7, 50)
stopifnot(all(abs(rowSums(M) - 1) < 1e-12), all(rowSums(M > 0) == 7))
M1 <- draw_sets(10, 1, 5); stopifnot(all(rowSums(M1 > 0) == 1))

# sparse product equals a loop over draws
mat <- matrix(rnorm(100 * 3), 100, 3)
idx <- lapply(seq_len(nrow(M)), function(i) which(M[i, ] > 0))
loop <- t(sapply(idx, function(ix) colMeans(mat[ix, , drop = FALSE])))
stopifnot(max(abs(as.matrix(M %*% mat) - loop)) < 1e-12)

# perm p-values on a known null
null <- matrix(c(1:99, 1:99 + 0.5), ncol = 2)
p <- perm_p_two(c(99, 20), null)
stopifnot(abs(p[1] - 2 * (1 + 1) / 100) < 1e-12)               # one null value >= 99
stopifnot(abs(p[2] - 2 * (1 + 19) / 100) < 1e-12)              # 1.5..19.5 are <= 20
stopifnot(abs(perm_p_greater(c(99, 200), null)[1] - 2 / 100) < 1e-12, perm_p_greater(c(99, 200), null)[2] == 1 / 100)
# names of the conditions are kept (results are indexed by condition name)
nn <- null; colnames(nn) <- c("a", "b")
stopifnot(identical(names(perm_p_two(c(a = 99, b = 20), nn)), c("a", "b")),
          identical(names(perm_p_greater(c(a = 99, b = 20), nn)), c("a", "b")))

# exact permutation equals brute force on a tiny example
x <- c(3, 1, 4, 1, 5); a <- c(TRUE, TRUE, FALSE, FALSE, FALSE)
e <- exact_perm_diff(x, a)
bf <- combn(5, 2, function(ix) mean(x[ix]) - mean(x[-ix]))
stopifnot(e$n_perm == 10, abs(e$p - mean(abs(bf) >= abs(e$obs) - 1e-12)) < 1e-12)

cat("test_h49_lib.R: all checks passed\n")
