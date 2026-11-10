# U8a: PC3 known-truth simulation with the prespec seed and draw order (55 normals per draw,
# beta 0 block then beta 0.30 block), fitted with the independent GLS.
here <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U8a/rederive"
source(file.path(here, "u8a_lib.R"))
ds <- rep(FIT_DS, c(37, 18)); st <- c(rep(0:2, c(4, 7, 26)), rep(0:2, c(5, 4, 9)))
sdv <- ifelse(ds == FIT_DS[1], 0.33, 0.95)
set.seed(20260923)
sims <- lapply(c(0, 0.30), function(b) lapply(1:2000, function(r) b * st + rnorm(55) * sdv))
res <- rbindlist(lapply(1:2, function(i) {
  b <- c(0, 0.30)[i]
  r <- t(sapply(sims[[i]], function(y) { g <- gls_fit(y, st, ds); c(g$beta, g$se) }))
  data.table(true_beta = b, coverage = mean(abs(r[, 1] - b) <= qt(0.975, 52) * r[, 2]),
             rejection = mean(2 * pt(abs(r[, 1] / r[, 2]), 52, lower.tail = FALSE) < 0.05), bias = mean(r[, 1] - b))
}))
print(res)
dir.create(file.path(here, "out/pc3_sameseed"), showWarnings = FALSE)
fwrite(res, file.path(here, "out/pc3_sameseed/pc3.tsv"), sep = "\t")
