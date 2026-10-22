#!/usr/bin/env Rscript
# H36 step 3: Dirichlet-multinomial stage test per cluster, within-cohort stage
# permutations, and the cohort-effect pipeline check. PRESPEC_H36.md sections 4, 5.
# Usage: Rscript h36_03_stage_test.R <repo> <build_dir> <out_dir> [n_perm=10] [limit=0]
suppressPackageStartupMessages({ library(data.table); library(parallel) })
args <- commandArgs(trailingOnly = TRUE)
REPO <- args[1]; BUILD <- args[2]; OUT <- args[3]
K_PERM <- if (length(args) >= 4) as.integer(args[4]) else 10L
LIMIT <- if (length(args) >= 5) as.integer(args[5]) else 0L
NC <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8"))
source(file.path(REPO, "scripts/analysis/hypotheses_2026q4/H36_splicing/h36_dm.R"))
set.seed(20260923)

source(file.path(REPO, "scripts/analysis/hypotheses_2026q4/H36_splicing/h36_prep.R"))
load_build(BUILD)

# within-cohort permutations of stage among staged persons (same labels for all clusters)
perm_stage <- lapply(seq_len(K_PERM), function(k) {
  set.seed(20260923L + k)
  s <- stage0
  for (d in unique(ds0)) {
    i <- which(ds0 == d & !is.na(stage0))
    if (length(i) > 1L) s[i] <- stage0[i][sample.int(length(i))]
  }
  s
})
saveRDS(perm_stage, file.path(OUT, "perm_stage_labels.rds"))

test_cluster <- function(k) {
  tryCatch({
    p <- prep(k, TRUE)
    Y <- p$Y; pi <- p$pi
    if (ncol(Y) < 2L || nrow(Y) < 150L) return(data.table(cl = k, status = "too_few_after_filters"))
    X0 <- design0(pi)
    st <- stage0[pi]; xs <- st - mean(st)
    f0 <- dm_fit(Y, X0)
    obs <- dm_lrt(Y, X0, xs, f0 = f0)
    xbar <- colMeans(X0); m <- mean(st)
    pi0 <- dm_pred(obs$f1, c(xbar, 0 - m)); pi4 <- dm_pred(obs$f1, c(xbar, 4 - m))
    dpsi <- as.numeric(pi4 - pi0)
    top <- which.max(abs(dpsi))
    pl <- vapply(perm_stage, function(ps) {
      xp <- ps[pi]; mk <- mean(xp)
      r <- dm_lrt(Y, X0, xp - mk, f0 = f0)
      d <- dm_pred(r$f1, c(xbar, 4 - mk)) - dm_pred(r$f1, c(xbar, 0 - mk))
      c(r$lr, max(abs(d)), max(r$f0$conv, r$f1$conv))
    }, numeric(3))
    dt <- data.table(cl = k, status = "ok", J_used = ncol(Y), n_persons = nrow(Y),
                     cohorts = paste(names(table(ds0[pi])), as.integer(table(ds0[pi])), sep = ":", collapse = ";"),
                     lr = obs$lr, df = obs$df, p = obs$p, conv0 = f0$conv, conv1 = obs$f1$conv,
                     phi0 = f0$phi, phi1 = obs$f1$phi,
                     dpsi_jr = paste(p$jr, signif(dpsi, 4), sep = ":", collapse = ";"),
                     top_jr = p$jr[top], top_dpsi = dpsi[top], max_abs_dpsi = max(abs(dpsi)))
    for (j in seq_len(ncol(pl))) {
      set(dt, j = paste0("lr_perm", j), value = pl[1, j])
      set(dt, j = paste0("dpsi_perm", j), value = pl[2, j])
      set(dt, j = paste0("conv_perm", j), value = pl[3, j])
    }
    dt
  }, error = function(e) data.table(cl = k, status = paste0("error: ", conditionMessage(e))))
}

todo <- seq_len(nrow(cl_tab)); if (LIMIT > 0L) todo <- todo[seq_len(min(LIMIT, length(todo)))]
chunks <- split(todo, ceiling(seq_along(todo) / 2000L))
t0 <- Sys.time(); res <- list()
for (ci in seq_along(chunks)) {
  f <- file.path(OUT, sprintf("chunk_%03d.rds", ci))
  if (file.exists(f)) { res[[ci]] <- readRDS(f); next }
  r <- rbindlist(mclapply(chunks[[ci]], test_cluster, mc.cores = NC, mc.preschedule = TRUE), fill = TRUE)
  saveRDS(r, f); res[[ci]] <- r
  cat(sprintf("chunk %d/%d done, %d clusters, elapsed %.1f min\n", ci, length(chunks), nrow(r),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
R <- rbindlist(res, fill = TRUE)
R <- merge(cl_tab[, .(cl, chr, str, start_min, end_max, J_tested_in_build = J)], R, by = "cl")
pcols <- grep("^lr_perm", names(R), value = TRUE)
allp <- grep("_perm[0-9]+$", names(R), value = TRUE)
fwrite(R[, setdiff(names(R), allp), with = FALSE], file.path(OUT, "stage_test_observed.tsv.gz"), sep = "\t")
fwrite(R[status == "ok", c("cl", "df", allp), with = FALSE], file.path(OUT, "stage_test_perm.tsv.gz"), sep = "\t")

# calibration summary: fraction of nominal p < 0.05 in observed and in each permutation
ok <- R[status == "ok"]
cal <- data.table(set = c("observed", paste0("perm", seq_along(pcols))),
                  n_clusters = nrow(ok),
                  frac_p_lt_0.05 = c(mean(ok$p < 0.05),
                                     vapply(pcols, function(cn) mean(pchisq(ok[[cn]], ok$df, lower.tail = FALSE) < 0.05), numeric(1))),
                  n_p_lt_1e4 = c(sum(ok$p < 1e-4),
                                 vapply(pcols, function(cn) sum(pchisq(ok[[cn]], ok$df, lower.tail = FALSE) < 1e-4), numeric(1))))
fwrite(cal, file.path(OUT, "calibration_by_set.tsv"), sep = "\t")
print(cal)
cat("status counts:\n"); print(table(R$status))

# ---- pipeline check: cohort effect on 2,000 random clusters (not biological) -----
set.seed(20260923)
samp <- sample(todo, min(2000L, length(todo)))
cohort_check <- function(k) {
  tryCatch({
    p <- prep(k, FALSE); Y <- p$Y; pi <- p$pi
    d <- factor(ds0[pi]); if (ncol(Y) < 2L || nlevels(d) < 2L) return(NULL)
    X0 <- independent_cols(model.matrix(~ sx, data.frame(sx = sex0[pi])))
    Xc <- model.matrix(~ d)[, -1, drop = FALSE]
    r <- dm_lrt(Y, X0, Xc)
    data.table(cl = k, n_persons = nrow(Y), n_cohorts = nlevels(d), lr = r$lr, df = r$df, p = r$p)
  }, error = function(e) NULL)
}
cc <- rbindlist(mclapply(samp, cohort_check, mc.cores = NC))
fwrite(cc, file.path(OUT, "cohort_effect_pipeline_check.tsv"), sep = "\t")
cat("cohort-effect check: clusters", nrow(cc), "| median p", signif(median(cc$p), 3),
    "| fraction p<1e-4", round(mean(cc$p < 1e-4), 3), "\n")
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo_03.txt"))
cat("total elapsed", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
