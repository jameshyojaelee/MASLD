# Beta-binomial regression by maximum likelihood with exact gradients (used by 02_analysis.R).
# Model: y_i | n_i ~ BetaBinomial(n_i, mu_i, rho), logit(mu_i) = X_i beta, rho constant (aod's phi:
# Var = n mu (1 - mu) [1 + (n - 1) rho]). Parameters: c(beta, g) with rho = plogis(g).
# Why this exists: aod::betabin (default BFGS, numerical gradient) stopped at different points
# for the same model depending on covariate centering, with convergence code 0 (H35 audit).
# Per-person read counts here are 1e6 to 1e8 and rho is about 3e-4, so the likelihood is very flat
# in some directions and a numerical-gradient optimizer stops early.

bb_negll <- function(par, X, y, n) {
  p <- ncol(X)
  mu <- plogis(drop(X %*% par[seq_len(p)]))
  s <- exp(-par[p + 1L])                      # s = a + b = (1 - rho) / rho
  a <- mu * s; b <- (1 - mu) * s
  -sum(lbeta(y + a, n - y + b) - lbeta(a, b))
}

bb_grad <- function(par, X, y, n) {
  p <- ncol(X)
  mu <- plogis(drop(X %*% par[seq_len(p)]))
  s <- exp(-par[p + 1L])
  a <- mu * s; b <- (1 - mu) * s
  d_mu <- s * (digamma(y + a) - digamma(n - y + b) - digamma(a) + digamma(b))
  d_s <- mu * digamma(y + a) + (1 - mu) * digamma(n - y + b) - digamma(n + s) -
    mu * digamma(a) - (1 - mu) * digamma(b) + digamma(s)
  -c(crossprod(X, d_mu * mu * (1 - mu)), sum(d_s * (-s)))
}

# Newton polishing of a parameter vector with the exact gradient and a numerical Jacobian of it.
bb_hess <- function(par, X, y, n) { H <- numDeriv::jacobian(function(q) bb_grad(q, X, y, n), par); (H + t(H)) / 2 }
bb_polish <- function(par, X, y, n, max_it = 25L) {
  for (it in seq_len(max_it)) {
    g <- bb_grad(par, X, y, n); H <- bb_hess(par, X, y, n)
    step <- tryCatch(solve(H, g), error = function(e) NULL)
    if (is.null(step) || anyNA(step)) break
    f0 <- bb_negll(par, X, y, n); tt <- 1; moved <- FALSE
    while (tt > 1e-6) {                                 # step halving
      new <- par - tt * step; f1 <- bb_negll(new, X, y, n)
      if (is.finite(f1) && f1 <= f0 + 1e-12) { moved <- TRUE; break }
      tt <- tt / 2
    }
    if (!moved) break
    par <- new
    if (abs(f0 - f1) < 1e-13) break
  }
  par
}

# Multi-start maximum likelihood, every start polished by Newton steps. Returns the best optimum and a verdict:
# converged = positive-definite observed information, Newton decrement g'H^-1 g < tol_nd at the best optimum,
# and at least two starts reaching the same objective (within 1e-6), which guards against a local optimum.
# nlminb's own convergence code is recorded but not used: PORT reports false convergence at tight tolerances.
bb_fit <- function(X, y, n, tol_nd = 1e-6) {
  p <- ncol(X)
  b_glm <- tryCatch(suppressWarnings(coef(glm(cbind(y, n - y) ~ X - 1, family = quasibinomial()))), error = function(e) NULL)
  b_ols <- tryCatch(coef(lm(qlogis((y + 0.5) / (n + 1)) ~ X - 1)), error = function(e) NULL)
  bases <- Filter(function(b) !is.null(b) && !anyNA(b), list(unname(b_glm), unname(b_ols)))
  starts <- list()
  for (b0 in bases) for (rho0 in c(1e-5, 1e-4, 1e-3, 1e-2)) starts[[length(starts) + 1L]] <- c(b0, qlogis(rho0))
  objs <- rep(NA_real_, length(starts)); pars <- vector("list", length(starts)); codes <- rep(NA_integer_, length(starts))
  for (i in seq_along(starts)) {
    o <- tryCatch(suppressWarnings(nlminb(starts[[i]], bb_negll, bb_grad, X = X, y = y, n = n,
                                           control = list(iter.max = 2000L, eval.max = 4000L))),
                  error = function(e) NULL)
    if (is.null(o) || !is.finite(o$objective)) next
    pp <- bb_polish(o$par, X, y, n); ob <- bb_negll(pp, X, y, n)
    if (!is.finite(ob) || ob > o$objective + 1e-8) { pp <- o$par; ob <- o$objective }
    pars[[i]] <- pp; objs[i] <- ob; codes[i] <- o$convergence
  }
  if (all(is.na(objs))) return(list(converged = FALSE, why = "no start produced a finite likelihood"))
  ib <- which.min(objs); par <- pars[[ib]]
  n_agree <- sum(abs(objs - objs[ib]) < 1e-6, na.rm = TRUE)
  g <- bb_grad(par, X, y, n); H <- bb_hess(par, X, y, n)
  V <- tryCatch(solve(H), error = function(e) NULL)
  pd <- !is.null(V) && all(is.finite(V)) && all(diag(V) > 0) && all(eigen(H, symmetric = TRUE, only.values = TRUE)$values > 0)
  nd <- if (pd) drop(t(g) %*% V %*% g) else NA_real_
  conv <- pd && is.finite(nd) && nd < tol_nd && n_agree >= 2L
  why <- if (conv) "" else paste0("pd ", pd, "; newton decrement ", signif(nd, 3), "; starts agreeing ", n_agree, "; nlminb code ", codes[ib])
  se <- if (pd) sqrt(diag(V)) else rep(NA_real_, p + 1L)
  list(converged = conv, why = why, par = par, se = se, nd = nd, n_agree = n_agree,
       logL = -objs[ib] + sum(lchoose(n, y)), rho = plogis(par[p + 1L]))
}
