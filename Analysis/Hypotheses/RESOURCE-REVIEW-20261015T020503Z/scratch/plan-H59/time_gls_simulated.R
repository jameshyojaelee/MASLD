# Timing only, on simulated data with the Ductular injury design (no real scores read).
suppressPackageStartupMessages(library(nlme))
set.seed(20260923)
ds <- factor(c(rep("A",37), rep("B",18)))
st <- c(rep(0,4), rep(1,7), rep(2,26), rep(0,5), rep(1,4), rep(2,9))
y <- rnorm(55, sd = ifelse(ds=="A", 0.33, 0.95))
t0 <- proc.time()[["elapsed"]]
for (i in 1:200) {
  f <- gls(y ~ st + ds, weights = varIdent(form = ~ 1 | ds), method = "REML")
}
el <- proc.time()[["elapsed"]] - t0
cat(sprintf("seconds per gls fit: %.4f\n", el/200))
