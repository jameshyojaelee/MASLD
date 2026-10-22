#!/usr/bin/env Rscript
# Checks for lib_bb.R: exact gradient, recovery on simulated data at the real read depths, agreement with
# aod::betabin where aod converges, and invariance to covariate centering (the property aod failed).
suppressPackageStartupMessages({ library(aod); library(numDeriv) })
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "lib_bb.R"))
set.seed(20260923)

sim <- function(N, nmed, beta, rho) {
  x1 <- rnorm(N); x2 <- rnorm(N, 5, 0.5); x3 <- sample(0:4, N, TRUE)
  X <- cbind(1, x1, x2, x3)
  n <- pmax(100, round(nmed * exp(rnorm(N, 0, 0.5))))
  mu <- plogis(drop(X %*% beta)); s <- (1 - rho) / rho
  p <- rbeta(N, mu * s, (1 - mu) * s)
  y <- rbinom(N, n, p)
  list(X = X, y = y, n = n)
}
ok <- TRUE
chk <- function(name, cond) { cat(sprintf("%-70s %s\n", name, if (isTRUE(cond)) "PASS" else "FAIL")); if (!isTRUE(cond)) ok <<- FALSE }

# 1. gradient: exact check at small counts (no cancellation), looser relative-to-norm check at real depth
beta <- c(qlogis(0.03), 0.2, 0.0, -0.02)
d <- sim(60, 300, beta, 0.02)
for (par in list(c(beta, qlogis(0.02)), c(beta + 0.05, qlogis(0.05)), c(beta - 0.1, qlogis(0.005)))) {
  g_an <- bb_grad(par, d$X, d$y, d$n); g_nu <- numDeriv::grad(bb_negll, par, X = d$X, y = d$y, n = d$n)
  e <- max(abs(g_an - g_nu) / pmax(1, abs(g_nu)))
  chk(sprintf("gradient matches numDeriv at n~300 (max rel err %.2e)", e), e < 1e-6)
}
d <- sim(100, 4e6, beta, 3e-4)
for (par in list(c(beta + 0.05, qlogis(1e-3)), c(beta - 0.1, qlogis(1e-4)))) {
  g_an <- bb_grad(par, d$X, d$y, d$n); g_nu <- numDeriv::grad(bb_negll, par, X = d$X, y = d$y, n = d$n)
  e <- max(abs(g_an - g_nu)) / max(1, sqrt(sum(g_nu^2)))
  chk(sprintf("gradient matches numDeriv at n~4e6 (error / norm %.2e)", e), e < 1e-5)
}

# 2. recovery at real depth (4e6 reads/person, rho 3e-4), 40 replicates
cover <- matrix(NA, 40, 4); rho_hat <- numeric(40); conv <- logical(40)
for (r in 1:40) {
  d <- sim(100, 4e6, beta, 3e-4)
  f <- bb_fit(d$X, d$y, d$n); conv[r] <- f$converged
  if (f$converged) { cover[r, ] <- abs(f$par[1:4] - beta) < 1.96 * f$se[1:4]; rho_hat[r] <- f$rho } else { rho_hat[r] <- NA; cat("  replicate", r, "not converged:", f$why, "\n") }
}
chk(sprintf("converged in %d/40 replicates", sum(conv)), sum(conv) >= 39)
chk(sprintf("95%% CI coverage per coefficient %s", paste(sprintf("%.2f", colMeans(cover, na.rm = TRUE)), collapse = " ")),
    all(colMeans(cover, na.rm = TRUE) > 0.85))
chk(sprintf("median rho_hat %.2e vs truth 3e-4", median(rho_hat, na.rm = TRUE)), abs(log(median(rho_hat, na.rm = TRUE) / 3e-4)) < 0.3)

# 3. agreement with aod where aod converges (moderate counts)
d <- sim(150, 300, beta, 0.02)
df <- data.frame(y = d$y, f = d$n - d$y, x1 = d$X[, 2], x2 = d$X[, 3], x3 = d$X[, 4])
fa <- aod::betabin(cbind(y, f) ~ x1 + x2 + x3, ~ 1, data = df, control = list(maxit = 5000))
fm <- bb_fit(d$X, d$y, d$n)
if (!fm$converged) cat("  mine:", fm$why, "\n")
chk(sprintf("aod code %d, mine converged %s", fa@code, fm$converged), fa@code == 0 && fm$converged)
cat("  coef mine:", signif(fm$par[1:4], 4), " aod:", signif(fa@param[1:4], 4), " se mine:", signif(fm$se[1:4], 3), "\n")
chk(sprintf("mine is at least as good as aod (logL mine - aod = %.4f)", fm$logL - fa@logL), fm$logL - fa@logL > -1e-6)
chk(sprintf("coefficients agree with aod within 0.25 SE (max %.2f SE)", max(abs(fm$par[1:4] - fa@param[1:4]) / fm$se[1:4])), max(abs(fm$par[1:4] - fa@param[1:4]) / fm$se[1:4]) < 0.25)
chk(sprintf("rho agrees with aod (%.4f vs %.4f)", fm$rho, as.numeric(fa@random.param[1])), abs(fm$rho - as.numeric(fa@random.param[1])) < 1e-3)
# my likelihood evaluated at aod's own solution must equal aod's logL (same likelihood, same constants)
ll_at_aod <- -bb_negll(c(fa@param[1:4], qlogis(as.numeric(fa@random.param[1]))), d$X, d$y, d$n) + sum(lchoose(d$n, d$y))
chk(sprintf("likelihood at aod's solution equals aod logL (%.4f vs %.4f)", ll_at_aod, fa@logL), abs(ll_at_aod - fa@logL) < 1e-3)

# 4. centering invariance (aod failed this on the real cohorts)
d <- sim(100, 4e6, beta, 3e-4)
Xc <- d$X; Xc[, 2:4] <- sweep(d$X[, 2:4], 2, colMeans(d$X[, 2:4]))
f1 <- bb_fit(d$X, d$y, d$n); f2 <- bb_fit(Xc, d$y, d$n)
if (!f1$converged) cat("  uncentered:", f1$why, "\n"); if (!f2$converged) cat("  centered:", f2$why, "\n")
chk(sprintf("both converge (%s, %s)", f1$converged, f2$converged), f1$converged && f2$converged)
chk(sprintf("same maximum with and without centering (dlogL %.2e)", abs(f1$logL - f2$logL)), abs(f1$logL - f2$logL) < 1e-6)
chk(sprintf("same slopes with and without centering (max abs diff %.2e)", max(abs(f1$par[2:4] - f2$par[2:4]))), max(abs(f1$par[2:4] - f2$par[2:4])) < 1e-6)
cat(if (ok) "ALL PASS\n" else "SOME FAILED\n")
quit(status = if (ok) 0L else 1L)
