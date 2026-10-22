#!/usr/bin/env Rscript
# H36 optimizer check, part 2 (L-BFGS-B memory size and tolerance; nlminb as the reference optimum): the first stage-test run (L-BFGS-B, maxit 300) hit the iteration limit in about 25% of
# null fits (share rises with the number of junctions J). Compare optimizers on non-converged and converged
# clusters using only log-likelihoods, convergence codes and run time (no stage association is printed).
# Usage: Rscript optim_compare.R <repo> <build_dir> <chunk_dir> <out_file>
suppressPackageStartupMessages({ library(data.table); library(parallel) })
args <- commandArgs(trailingOnly = TRUE); REPO <- args[1]; BUILD <- args[2]; CH <- args[3]; OUTF <- args[4]
NC <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8"))
HERE <- file.path(REPO, "scripts/analysis/hypotheses_2026q4/H36_splicing")
source(file.path(HERE, "h36_dm.R")); source(file.path(HERE, "h36_prep.R")); load_build(BUILD)
R <- rbindlist(lapply(list.files(CH, "^chunk_.*rds$", full.names = TRUE), readRDS), fill = TRUE)[status == "ok", .(cl, J_used, conv0, conv1)]
set.seed(20260923)
VARIANTS <- list(lmm50_f1e5 = list(lmm = 50L, factr = 1e5), lmm50_f1e3 = list(lmm = 50L, factr = 1e3),
                 lmm50_f10 = list(lmm = 50L, factr = 10), lmm100_f1e3 = list(lmm = 100L, factr = 1e3))
nonc <- R[conv0 != 0 & J_used >= 4]; conv <- R[conv0 == 0 & J_used >= 4]
pick <- c(sample(nonc$cl, min(48L, nrow(nonc))), sample(conv$cl, min(16L, nrow(conv))))
fit_with <- function(method, Y, X, start = NULL) {
  fnl <- dm_make_fn(Y, X); J <- fnl$J; p <- fnl$p; kB <- p * (J - 1L)
  if (is.null(start)) { tot <- colSums(Y) + 0.5; B0 <- matrix(0, p, J - 1L); B0[1, ] <- log(tot[-J] / tot[J]); start <- c(as.vector(B0), log(20)) }
  k <- length(start); lo <- c(rep(-Inf, k - 1L), log(0.05)); hi <- c(rep(Inf, k - 1L), log(1e5))
  t0 <- proc.time()[["elapsed"]]
  r <- if (method == "nlminb") {
    o <- nlminb(start, fnl$fn, fnl$gr, lower = lo, upper = hi, control = list(iter.max = 3000L, eval.max = 6000L, rel.tol = 1e-12, x.tol = 1e-10))
    list(ll = -o$objective, code = o$convergence, it = o$evaluations[["function"]], par = o$par)
  } else {
    v <- VARIANTS[[method]]
    o <- optim(start, fnl$fn, fnl$gr, method = "L-BFGS-B", lower = lo, upper = hi,
               control = list(maxit = 5000L, factr = v$factr, lmm = v$lmm))
    list(ll = -o$value, code = o$convergence, it = o$counts[[1]], par = o$par)
  }
  r$sec <- proc.time()[["elapsed"]] - t0; r
}
one <- function(k) {
  p <- prep(k, TRUE); X0 <- design0(p$pi); xs <- stage0[p$pi] - mean(stage0[p$pi]); X1 <- cbind(X0, xs)
  rbindlist(lapply(c(names(VARIANTS), "nlminb"), function(m) {
    f0 <- fit_with(m, p$Y, X0)
    J <- ncol(p$Y); q <- ncol(X0)
    s1 <- c(as.vector(rbind(matrix(f0$par[seq_len(q * (J - 1L))], q, J - 1L), matrix(0, 1, J - 1L))), f0$par[length(f0$par)])
    f1 <- fit_with(m, p$Y, X1, s1)
    data.table(cl = k, J = J, n = nrow(p$Y), method = m, ll0 = f0$ll, code0 = f0$code, it0 = f0$it, sec0 = f0$sec,
               ll1 = f1$ll, code1 = f1$code, it1 = f1$it, sec1 = f1$sec)
  }))
}
res <- rbindlist(mclapply(pick, one, mc.cores = NC, mc.preschedule = FALSE), fill = TRUE)
res[, `:=`(ref0 = max(ll0), ref1 = max(ll1), lr_raw = 2 * pmax(0, ll1 - ll0)), by = cl]
fwrite(res, OUTF, sep = "\t")
res[, grp := ifelse(cl %in% nonc$cl, "nonconverged_in_run1", "converged_in_run1")]
print(res[, .(n = .N, code0_ok = mean(code0 == 0), code1_ok = mean(code1 == 0), med_it0 = median(it0), med_sec0 = round(median(sec0), 2),
              med_sec1 = round(median(sec1), 2), max_gap_ll0 = max(ref0 - ll0), med_gap_ll0 = median(ref0 - ll0),
              max_gap_ll1 = max(ref1 - ll1)), by = .(grp, method)])
print(res[, .(max_abs_lr_diff_vs_nlminb = max(abs(lr_raw - lr_raw[method == "nlminb"]))), by = .(cl, method)][, .(max_abs_lr_diff = max(max_abs_lr_diff_vs_nlminb)), by = method])
