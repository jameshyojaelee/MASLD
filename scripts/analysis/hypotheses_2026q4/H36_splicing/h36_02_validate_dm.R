#!/usr/bin/env Rscript
# H36 step 2: validation of the base-R Dirichlet-multinomial likelihood-ratio test.
# PRESPEC_H36.md section 6 and Addendum 2.
#   (i)   analytic vs numerical gradient on 20 real clusters (gate)
#   (ii)  simulated null with the real staged design, stage permuted within cohort (gate)
#   (iii) spike-in power and dPSI_0to4 recovery (reported)
#   (iv)  cross-check of the LR statistic against VGAM::dirmultinomial (reported, not a gate)
# Usage: Rscript h36_02_validate_dm.R <repo> <build_dir> <out_dir>
# Exit status 3 when a gate fails, so a dependent job does not start.
suppressPackageStartupMessages({ library(data.table); library(parallel); requireNamespace("VGAM") })
args <- commandArgs(trailingOnly = TRUE)
REPO <- args[1]; BUILD <- args[2]; OUT <- args[3]
NC <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8"))
HERE <- file.path(REPO, "scripts/analysis/hypotheses_2026q4/H36_splicing")
source(file.path(HERE, "h36_dm.R"))
source(file.path(HERE, "h36_prep.R"))
load_build(BUILD)
RNGkind("L'Ecuyer-CMRG")
t0 <- Sys.time()
elapsed <- function() round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
gate <- list()

# ---- (i) gradient check -------------------------------------------------------
num_grad <- function(f, x, h = 1e-3) {
  vapply(seq_along(x), function(k) {
    e <- numeric(length(x)); e[k] <- h
    (-f(x + 2 * e) + 8 * f(x + e) - 8 * f(x - e) + f(x - 2 * e)) / (12 * h)
  }, numeric(1))
}
set.seed(20260923)
ok_cl <- which(cl_tab$n_staged >= 150L)
gcl <- sort(sample(ok_cl, 20L))
grad_one <- function(k) {
  set.seed(20260923L + k)
  p <- prep(k, TRUE)
  X1 <- cbind(design0(p$pi), stage0[p$pi] - mean(stage0[p$pi]))
  fn <- dm_make_fn(p$Y, X1)
  J <- ncol(p$Y); kpar <- ncol(X1) * (J - 1L) + 1L
  fit <- dm_fit(p$Y, X1)
  start <- c(rnorm(kpar - 1L, 0, 0.3), log(20))
  rbindlist(lapply(list(perturbed = start, fitted = fit$par), function(par) {
    ga <- fn$gr(par); gn <- num_grad(fn$fn, par)
    data.table(cl = k, J = J, n_par = length(par), max_abs_grad = max(abs(ga)),
               max_err = max(abs(ga - gn) / pmax(1, abs(ga))))
  }), idcol = "point")
}
G <- rbindlist(mclapply(gcl, grad_one, mc.cores = NC))
fwrite(G, file.path(OUT, "validation_gradient.tsv"), sep = "\t")
gate$gradient <- max(G$max_err) < 1e-5
cat(sprintf("[%s min] (i) gradient: max error %.2e over %d clusters x 2 points -> %s\n", elapsed(),
            max(G$max_err), length(gcl), if (gate$gradient) "PASS" else "FAIL"))

# ---- (ii) null simulation and (iii) spike-in -----------------------------------
staged <- which(!is.na(stage0) & ds0 != "GSE126848")
N <- length(staged)
cohort_f <- factor(ds0[staged]); coh <- as.integer(cohort_f); nC <- nlevels(cohort_f)
sex_s <- sex0[staged]; st_real <- stage0[staged]
X0s <- independent_cols(model.matrix(~ d + sx, data.frame(d = cohort_f, sx = sex_s)))
cat("simulation design:", N, "staged persons;", paste(names(table(cohort_f)), table(cohort_f), collapse = ", "), "\n")

rdirmult <- function(P, phi, n) {
  G <- matrix(rgamma(length(P), shape = phi * as.vector(P)), nrow(P), ncol(P))
  Pr <- G / rowSums(G)
  t(vapply(seq_len(nrow(P)), function(i) as.integer(rmultinom(1L, n[i], Pr[i, ])), integer(ncol(P))))
}
make_setting <- function(J, phi, s) {
  set.seed(20260923L + 7919L * s)
  base <- c(0.3, rep(0.7 / (J - 1L), J - 1L))
  list(J = J, phi = phi, s = s, eta0 = log(base[-J] / base[J]),
       cs = matrix(rnorm(nC * (J - 1L), 0, 0.5), nC, J - 1L), sx = c(0.2, rep(0, J - 2L)))
}
eta_nostage <- function(S) {
  matrix(S$eta0, N, S$J - 1L, byrow = TRUE) + S$cs[coh, , drop = FALSE] + outer(sex_s, S$sx)
}
gen_counts <- function(S, beta, st) {
  eta <- eta_nostage(S) + outer(st, beta)
  n <- pmax(20L, as.integer(round(rlnorm(N, log(100), 0.7))))
  rdirmult(softmax_rows(cbind(eta, 0)), S$phi, n)
}
permute_within <- function() {
  s <- st_real
  for (d in seq_len(nC)) { i <- which(coh == d); s[i] <- st_real[i][sample.int(length(i))] }
  s
}
test_counts <- function(Y, st) {
  o <- order(colSums(Y)); Y <- Y[, o, drop = FALSE]
  r <- dm_lrt(Y, X0s, st - mean(st))
  xbar <- colMeans(X0s); m <- mean(st)
  dp <- as.numeric(dm_pred(r$f1, c(xbar, 4 - m)) - dm_pred(r$f1, c(xbar, 0 - m)))
  list(p = r$p, lr = r$lr, conv = max(r$f0$conv, r$f1$conv), dpsi_junction1 = dp[which(o == 1L)])
}

sets <- CJ(J = c(2L, 3L, 5L), phi = c(5, 20, 100))
sets[, s := .I]
run_null <- function(S, b, seed_base) {
  set.seed(seed_base + b)
  st <- permute_within()
  Y <- gen_counts(S, rep(0, S$J - 1L), st_real)    # counts carry no stage effect
  r <- test_counts(Y, st)
  c(p = r$p, conv = r$conv)
}
null_res <- rbindlist(lapply(seq_len(nrow(sets)), function(i) {
  S <- make_setting(sets$J[i], sets$phi[i], sets$s[i])
  m <- do.call(rbind, mclapply(seq_len(1000L), function(b) run_null(S, b, 20260923L), mc.cores = NC))
  data.table(J = S$J, phi = S$phi, n_sim = 1000L, pass1_rate05 = mean(m[, "p"] < 0.05),
             pass1_rate01 = mean(m[, "p"] < 0.01), conv_fail = sum(m[, "conv"] != 0),
             p = list(m[, "p"]))
}))
cat(sprintf("[%s min] (ii) pass 1 done\n", elapsed()))
# Addendum 2: settings outside 0.035-0.065 get 5,000 more simulations; gate on the pooled 6,000
null_res[, `:=`(n_total = n_sim, rate05 = pass1_rate05, rate01 = pass1_rate01, extra = FALSE)]
for (i in which(null_res$pass1_rate05 < 0.035 | null_res$pass1_rate05 > 0.065)) {
  S <- make_setting(sets$J[i], sets$phi[i], sets$s[i])
  m2 <- do.call(rbind, mclapply(seq_len(5000L), function(b) run_null(S, b, 20260924L + 10000L * sets$s[i]), mc.cores = NC))
  allp <- c(null_res$p[[i]], m2[, "p"])
  set(null_res, i, "n_total", 6000L); set(null_res, i, "rate05", mean(allp < 0.05))
  set(null_res, i, "rate01", mean(allp < 0.01)); set(null_res, i, "extra", TRUE)
  set(null_res, i, "conv_fail", null_res$conv_fail[i] + sum(m2[, "conv"] != 0))
}
null_res[, in_window := rate05 >= 0.035 & rate05 <= 0.065]
fwrite(null_res[, setdiff(names(null_res), "p"), with = FALSE], file.path(OUT, "validation_null_sim.tsv"), sep = "\t")
gate$null <- all(null_res$in_window)
print(null_res[, .(J, phi, n_total, rate05, rate01, conv_fail, extra, in_window)])
cat(sprintf("[%s min] (ii) null simulation -> %s\n", elapsed(), if (gate$null) "PASS" else "FAIL"))

# spike-in: stage effect on junction 1 (log-odds per stage step), phi = 20
spk <- CJ(J = c(2L, 3L, 5L), shift = c(0.05, 0.10, 0.20))
spk[, s := .I + 100L]
spike_res <- rbindlist(lapply(seq_len(nrow(spk)), function(i) {
  S <- make_setting(spk$J[i], 20, spk$s[i])
  beta <- c(spk$shift[i], rep(0, S$J - 2L))
  m <- do.call(rbind, mclapply(seq_len(1000L), function(b) {
    set.seed(20260923L + 500000L + b)
    r <- test_counts(gen_counts(S, beta, st_real), st_real)
    c(p = r$p, dpsi = r$dpsi_junction1, conv = r$conv)
  }, mc.cores = NC))
  eb <- colMeans(eta_nostage(S))
  true_dpsi <- softmax_rows(rbind(c(eb + beta * 4, 0)))[1] - softmax_rows(rbind(c(eb, 0)))[1]   # mean non-stage logit; stage 0 vs 4
  data.table(J = S$J, shift_per_step = spk$shift[i], phi = 20, n_sim = 1000L,
             true_dpsi_0to4 = true_dpsi, mean_est_dpsi = mean(m[, "dpsi"]),
             power_05 = mean(m[, "p"] < 0.05), power_1e4 = mean(m[, "p"] < 1e-4),
             frac_est_dpsi_ge_0.10 = mean(abs(m[, "dpsi"]) >= 0.10), conv_fail = sum(m[, "conv"] != 0))
}))
fwrite(spike_res, file.path(OUT, "validation_spike.tsv"), sep = "\t")
print(spike_res)
cat(sprintf("[%s min] (iii) spike-in done\n", elapsed()))

# ---- (iv) VGAM cross-check on 30 real clusters ---------------------------------
set.seed(20260923)
vcl <- sort(sample(ok_cl, 30L))
vg_one <- function(k) {
  p <- prep(k, TRUE); Y <- p$Y; X0 <- design0(p$pi); xs <- stage0[p$pi] - mean(stage0[p$pi])
  ours <- dm_lrt(Y, X0, xs)
  dat <- as.data.frame(X0[, -1, drop = FALSE]); names(dat) <- paste0("c", seq_len(ncol(dat))); dat$stage <- xs
  Yc <- Y
  fit <- function(rhs) {
    setTimeLimit(elapsed = 180, transient = TRUE)
    VGAM::vglm(as.formula(paste("Yc ~", rhs)), VGAM::dirmultinomial, data = dat, maxit = 200, trace = FALSE)
  }
  r0 <- if (ncol(dat) > 1L) paste(names(dat)[-ncol(dat)], collapse = " + ") else "1"
  v <- tryCatch({
    f0 <- fit(r0); f1 <- fit(paste(r0, "+ stage"))
    list(lr = 2 * (as.numeric(stats::logLik(f1)) - as.numeric(stats::logLik(f0))), status = "ok")
  }, error = function(e) list(lr = NA_real_, status = paste("fail:", substr(conditionMessage(e), 1, 80))))
  data.table(cl = k, J_used = ncol(Y), n_persons = nrow(Y), lr_ours = ours$lr, lr_vgam = v$lr, vgam_status = v$status)
}
V <- rbindlist(mclapply(vcl, vg_one, mc.cores = NC))
fwrite(V, file.path(OUT, "validation_vgam.tsv"), sep = "\t")
vok <- V[vgam_status == "ok"]
cat(sprintf("[%s min] (iv) VGAM: %d of %d fits ok; correlation of LR (ours vs VGAM) = %s; median |difference| = %s\n",
            elapsed(), nrow(vok), nrow(V),
            if (nrow(vok) > 2) signif(cor(vok$lr_ours, vok$lr_vgam), 4) else NA,
            if (nrow(vok) > 0) signif(median(abs(vok$lr_ours - vok$lr_vgam)), 3) else NA))

gate_pass <- isTRUE(gate$gradient) && isTRUE(gate$null)
writeLines(c(if (gate_pass) "PASS" else "FAIL",
             sprintf("gradient=%s null=%s", gate$gradient, gate$null)), file.path(OUT, "GATE.txt"))
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo_02.txt"))
cat("GATE:", if (gate_pass) "PASS" else "FAIL", "| total", elapsed(), "min\n")
if (!gate_pass) quit(status = 3L)
