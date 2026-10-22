# Debug on synthetic AR(1) LD only (no real outcome): RNG save/restore around wrapped fits, purity-step draws,
# coloc_state paths, noise covariance, and the NC-3 statistic's null distribution.
H80 <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H80-20261015T130802Z"
Sys.setenv(H80_SCRIPTS = file.path(H80, "scripts"))
source(file.path(H80, "scripts/h80_lib.R"))
streams <- make_streams(20L)
ar1 <- function(p, rho) rho^abs(outer(seq_len(p), seq_len(p), "-"))
p1 <- 300; p2 <- 300
R <- matrix(0, p1 + p2, p1 + p2)
R[1:p1, 1:p1] <- ar1(p1, 0.97); R[p1 + 1:p2, p1 + 1:p2] <- ar1(p2, 0.6)
blk <- rep(1:2, c(p1, p2))
ids <- paste0("1:", 1000 + seq_len(p1 + p2))
be <- block_eigen(R, blk)
Rr <- reg_matrix(R, ids)
cat("block_eigen: blocks", be$n_blocks, "clipped", be$n_neg, "min eig", signif(be$min_eig, 3), "\n")

# 1. noise covariance against R + 1e-3 I
set_stream(streams, 1L)
Z <- vapply(1:4000, function(i) noise_from(be, rnorm(p1 + p2)), numeric(p1 + p2))
S_emp <- tcrossprod(Z) / ncol(Z)
S_true <- R; diag(S_true) <- diag(S_true) + LD_REGULARIZE
cat("noise covariance: max |emp - true| ", signif(max(abs(S_emp - S_true)), 3), " (sampling SD about ", signif(sqrt(2 / 4000), 3), ")\n")

# 2. RNG restore and purity draws
set_stream(streams, 2L)
z <- 6 * R[, 50] + noise_from(be, rnorm(p1 + p2)); names(z) <- ids
before <- get(".Random.seed", envir = .GlobalEnv)
fr <- fit_eqtl(z, Rr)
after <- get(".Random.seed", envir = .GlobalEnv)
cat("wrapped eQTL fit: restored ", fr$restored, ", stream identical ", identical(before, after), ", kind ", RNGkind()[1],
    ", drew ", fr$drew, ", n_cs ", fit_ncs(fr$res), ", cs sizes ", paste(lengths(fr$res$sets$cs), collapse = ";"), "\n")
x1 <- runif(1)
set_stream(streams, 2L); invisible(rnorm(p1 + p2)); x2 <- runif(1)
cat("stream continues as if no fit ran: ", identical(x1, x2), "\n")

# 3. coloc_state paths
mk <- function(lam, c, eps) { z <- lam * R[, c] + noise_from(be, eps); names(z) <- ids; z }
set_stream(streams, 3L)
eE <- rnorm(p1 + p2); eG <- rnorm(p1 + p2)
c <- 450; k <- 150   # different LD blocks, r2 = 0
fe15 <- fit_eqtl(mk(15, c, eE), Rr); fg15 <- fit_gwas(mk(15, c, eG), Rr, 50000)
fgk <- fit_gwas(mk(15, k, eG), Rr, 50000); fe0 <- fit_eqtl(mk(0, c, eE), Rr)
s1 <- coloc_state(fe15, fg15); s2 <- coloc_state(fe15, fgk); s3 <- coloc_state(fe0, fg15)
cat("shared (15,15): ", s1$state, " PP.H4 ", signif(s1$pp4, 3), " restored ", s1$coloc_restored, " drew ", s1$coloc_drew, "\n")
cat("distinct (15,15): ", s2$state, " PP.H4 ", signif(s2$pp4, 3), " PP.H3 ", signif(s2$pp3, 3), "\n")
cat("no eQTL (0,15): ", s3$state, "\n")
err <- structure(list(msg = "x"), class = "h80_error")
cat("error path: ", coloc_state(list(res = err), fg15)$state, "\n")
cat("RNG kind after all calls: ", RNGkind()[1], "\n")

# 4. NC-3 statistic under the simulator: 300 synthetic windows x 10 replicates
set_stream(streams, 4L)
T <- vapply(1:300, function(w) {
  s <- sum(vapply(1:10, function(r) sum(noise_from(be, rnorm(p1 + p2))^2), numeric(1)))
  (s - 10 * be$trS) / sqrt(2 * 10 * be$trS2)
}, numeric(1))
cat("NC-3 T_w over 300 synthetic windows: mean ", signif(mean(T), 3), " sd ", signif(sd(T), 3), " max |T| ", signif(max(abs(T)), 3),
    " eff df ", round(10 * be$trS^2 / be$trS2), "\n")

# 5. lambda-cross helper used by the assembler
src <- readLines(file.path(H80, "scripts/04_assemble.R"))
i0 <- grep("^lam_cross <- function", src); i1 <- i0 + 6
eval(parse(text = src[i0:i1]))
print(lam_cross(c(3, 4, 5, 6, 7, 8, 10, 15), c(0.1, 0.2, 0.4, 0.6, 0.7, 0.85, 0.9, 0.95), 0.5))
print(lam_cross(c(3, 4, 5), c(0.6, 0.7, 0.9), 0.5)); print(lam_cross(c(3, 4, 5), c(0.1, 0.2, 0.3), 0.5))
