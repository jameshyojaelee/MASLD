#!/usr/bin/env Rscript
# H36 step 2b: re-run of validation part (iii) only (spike-in), after correcting the true dPSI_0to4 formula
# in the first run's validation_spike.tsv (it was evaluated at the stage-centred baseline).
# Part (iii) of h36_02_validate_dm.R, unchanged except for that formula.
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

writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo_02b.txt"))
cat("total", elapsed(), "min\n")
