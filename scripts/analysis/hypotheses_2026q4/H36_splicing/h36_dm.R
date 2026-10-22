# H36: Dirichlet-multinomial regression in base R (new code; validated by
# h36_02_validate_dm.R before use). Modality: junction read counts within an
# intron cluster. Model: multinomial logit mean (reference = last column) and one
# concentration parameter phi per fit. Optimizer: L-BFGS-B, memory 100, factr 1e3, maxit 5000
# (the first settings, memory 5 / factr 1e5 / maxit 300, stopped at the iteration limit in 25% of fits). Log-likelihood without the multinomial
# coefficient (identical for nested models on the same data).

softmax_rows <- function(eta) {
  mx <- eta[cbind(seq_len(nrow(eta)), max.col(eta, ties.method = "first"))]
  e <- exp(eta - mx)
  e / rowSums(e)
}

dm_make_fn <- function(Y, X) {
  N <- nrow(Y); J <- ncol(Y); p <- ncol(X); n <- rowSums(Y)
  kB <- p * (J - 1L)
  last_par <- NULL; last <- NULL
  ev <- function(par) {
    if (!is.null(last_par) && identical(par, last_par)) return(last)
    B <- matrix(par[seq_len(kB)], p, J - 1L)
    phi <- exp(par[kB + 1L])
    pi <- pmax(softmax_rows(cbind(X %*% B, 0)), 1e-10)
    a <- phi * pi
    ll <- sum(lgamma(phi) - lgamma(phi + n)) + sum(lgamma(Y + a) - lgamma(a))
    g <- digamma(Y + a) - digamma(a) + (digamma(phi) - digamma(phi + n))
    h <- phi * g
    s <- rowSums(pi * h)
    deta <- pi * (h - s)
    gB <- crossprod(X, deta[, seq_len(J - 1L), drop = FALSE])
    gphi <- phi * sum(g * pi)
    last_par <<- par
    last <<- list(nll = -ll, grad = -c(as.vector(gB), gphi))
    last
  }
  list(fn = function(par) ev(par)$nll, gr = function(par) ev(par)$grad, p = p, J = J)
}

# X must have an intercept as its first column. Y: N x J integer matrix, last
# column is the reference junction.
dm_fit <- function(Y, X, start = NULL, maxit = 5000L) {
  fnl <- dm_make_fn(Y, X)
  J <- fnl$J; p <- fnl$p; kB <- p * (J - 1L)
  if (is.null(start)) {
    tot <- colSums(Y) + 0.5
    B0 <- matrix(0, p, J - 1L)
    B0[1, ] <- log(tot[-J] / tot[J])
    start <- c(as.vector(B0), log(20))
  }
  k <- length(start)
  o <- optim(start, fnl$fn, fnl$gr, method = "L-BFGS-B",
             lower = c(rep(-Inf, k - 1L), log(0.05)), upper = c(rep(Inf, k - 1L), log(1e5)),
             control = list(maxit = maxit, factr = 1e3, lmm = 100L))
  list(par = o$par, ll = -o$value, conv = o$convergence,
       B = matrix(o$par[seq_len(kB)], p, J - 1L), phi = exp(o$par[k]))
}

# Likelihood-ratio test for adding the columns x_new (N x q matrix or vector) to X0.
dm_lrt <- function(Y, X0, x_new, f0 = NULL) {
  x_new <- as.matrix(x_new); q <- ncol(x_new); J <- ncol(Y)
  if (is.null(f0)) f0 <- dm_fit(Y, X0)
  X1 <- cbind(X0, x_new)
  start1 <- c(as.vector(rbind(f0$B, matrix(0, q, J - 1L))), log(f0$phi))
  f1 <- dm_fit(Y, X1, start = start1)
  lr <- max(0, 2 * (f1$ll - f0$ll))
  df <- q * (J - 1L)
  list(lr = lr, df = df, p = pchisq(lr, df, lower.tail = FALSE), f0 = f0, f1 = f1, X1 = X1)
}

# Drop linearly dependent columns (keeps the intercept and earlier columns).
independent_cols <- function(X) {
  q <- qr(X)
  X[, q$pivot[seq_len(q$rank)], drop = FALSE][, order(q$pivot[seq_len(q$rank)]), drop = FALSE]
}

# Predicted junction proportions at covariate profile x (row vector matching X1).
dm_pred <- function(f1, x) softmax_rows(cbind(matrix(x, nrow = 1) %*% f1$B, 0))
