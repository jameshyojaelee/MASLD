# plan-H52: full-transport benchmark for S1 at GSE213621 df, from discovery-only edges.
# Reads only the frozen coarse rediscovery (job 3811, discovery cohorts only). No GSE213621 data.
suppressPackageStartupMessages(library(data.table))
rr <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z"
d <- fread(file.path(rr, "jobs/3811/out/approx_rediscovery_coarse_371.tsv.gz"))
f <- d[in_family == TRUE]
cat("family", nrow(f), " coordinated", f[coordinated == TRUE, .N], " noncoord", f[coordinated == FALSE, .N],
    " noncoord with z_part>0", f[coordinated == FALSE & z_part > 0, .N], " coord all positive", all(f[coordinated == TRUE]$z_part > 0), "\n")
cat("draw overlap", round(f[coordinated == TRUE, .N] / f[coordinated == FALSE, .N], 4), "\n")
cat("coordinated r quantiles:\n"); print(round(quantile(tanh(f[coordinated == TRUE]$z_part), c(.05,.25,.5,.75,.95)), 3))
pred <- function(z_part, zs_part, df, c) {
  pr <- pnorm(abs(atanh(c * tanh(z_part))) * sqrt(df - 3))
  pd <- pnorm(abs(zs_part))
  pd * pr + (1 - pd) * (1 - pr)
}
out <- rbindlist(lapply(c(254, 193, 105, 61), function(df) rbindlist(lapply(c(1, 0.75, 0.5, 0.25), function(c) {
  co <- f[coordinated == TRUE]; nc <- f[coordinated == FALSE]; ncp <- nc[z_part > 0]
  a <- mean(pred(co$z_part, co$zs_part, df, c)); b <- mean(pred(nc$z_part, nc$zs_part, df, c)); p <- mean(pred(ncp$z_part, ncp$zs_part, df, c))
  data.table(df = df, shrink = c, s1_coord = a, s1_noncoord_mixed = b, excess_mixed = a - b, s1_noncoord_pos = p, excess_pos = a - p)
}))))
print(out[, lapply(.SD, round, 4)])
boot_sd_197 <- 0.01038428
sd258 <- boot_sd_197 * sqrt(197 / 258)
cat(sprintf("bootstrap SD scaled to n=258: %.4f; MDE one-sided alpha 0.0125, power 0.8: %.4f; two-sided 95%% CI excludes 0 with power 0.8: %.4f\n",
            sd258, (qnorm(1 - 0.0125) + qnorm(0.8)) * sd258, (qnorm(0.975) + qnorm(0.8)) * sd258))
fwrite(out, file.path(rr, "scratch/plan-H52/benchmark_s1.tsv"), sep = "\t")
