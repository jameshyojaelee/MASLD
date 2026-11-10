# DIG1 job 35b (post hoc sensitivity to job 35, labelled post hoc; not a test). The cytosol factor CF
# is a mean of quantile-normalised protein z-scores, so a fall in mitochondrial proteins raises it
# mechanically; adjusting for it may remove part of the effect being measured. Same K1-K3 and CCT
# (C1, C2) models as jobs 33/35 with MF + six composition scores, without CF; and with the six
# composition scores alone. Seed 20261024. Usage: Rscript 35b_ci_composition_noCF.R <new_out_dir>
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(TRUE); stopifnot(length(args) == 1L, !dir.exists(args[1]))
out <- args[1]; dir.create(out, recursive = TRUE)
set.seed(20261024)
RUN <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/DIG1-20261024T113951Z"
COMP <- c("immune", "macrophage", "stroma", "ductular", "endothelial", "hepatocyte")
rd <- function(k) { a <- fread(file.path(RUN, "out_cct_composition", paste0("scores_", k, ".tsv")))
  b <- fread(file.path(RUN, "out_membrane_topology", paste0("scores_", k, ".tsv")))
  stopifnot(identical(a$sample, b$sample)); cbind(a, b[, .(CI, CII, CI_minus_CII, OXPHOS_mitoribo)]) }
D1 <- rd("PXD051911"); D2 <- rd("PXD027722"); D3 <- rd("PXD065892")
D1[, sex := factor(sex)][, batch := factor(batch)]
pcor <- function(z, a, b, vars) {
  z <- z[complete.cases(z[, c(a, b, vars), with = FALSE])]
  xx <- model.matrix(reformulate(if (length(vars)) vars else "1"), data = z)
  rr <- qr.resid(qr(xx), cbind(rank(z[[a]]), rank(z[[b]]))); r <- cor(rr[, 1], rr[, 2]); df <- nrow(xx) - ncol(xx) - 1
  data.table(n = nrow(z), df = df, estimate = r, p = 2 * pt(-abs(r * sqrt(df / (1 - r^2))), df)) }
fl <- function(z, y, g, vars, nperm) {
  z <- z[complete.cases(z[, c(y, vars), with = FALSE])]; lab <- as.numeric(z$grp == g)
  X0 <- model.matrix(reformulate(if (length(vars)) vars else "1"), data = z); X1 <- cbind(lab, X0)
  b1 <- function(yy) qr.coef(qr(X1), yy)[1]
  f0 <- qr.fitted(qr(X0), z[[y]]); e0 <- z[[y]] - f0; obs <- b1(z[[y]])
  nd <- replicate(nperm, b1(f0 + e0[sample.int(length(e0))]))
  data.table(n = nrow(z), n1 = sum(lab), df = nrow(z) - ncol(X1), estimate = unname(obs), sd_y = sd(z[[y]]), p = (1 + sum(abs(nd) >= abs(obs) - 1e-12)) / (nperm + 1)) }
BASE1 <- c("batch", "age", "bmi", "sex")
M1 <- D1[grp == "disease"]
ADJ <- list(MF_comp = c("MF", COMP), comp = COMP)
res <- rbindlist(lapply(names(ADJ), function(a) rbindlist(lapply(c("CI_minus_CII", "CCT"), function(y) rbind(
  cbind(adjust = a, y = y, what = "PXD051911 MASLD vs ballooning (rho)", pcor(M1, y, "Ballooning", c(BASE1, ADJ[[a]]))),
  cbind(adjust = a, y = y, what = "PXD065892 NAFLD vs healthy (coef)", fl(D3[grp %in% c("disease", "control")], y, "disease", ADJ[[a]], 20000)),
  cbind(adjust = a, y = y, what = "PXD065892 AIH vs healthy (coef)", fl(D3[grp %in% c("AIH", "control")], y, "AIH", ADJ[[a]], 20000)),
  cbind(adjust = a, y = y, what = "PXD051911 MASLD vs control (coef)", fl(D1[grp %in% c("disease", "control")], y, "disease", c(BASE1, ADJ[[a]]), 20000)), fill = TRUE)))))
res[, sd_units := estimate / sd_y]
fwrite(res, file.path(out, "posthoc_noCF.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
options(width = 200); print(res[, .(adjust, y, what, n, df, estimate, sd_units, p)], digits = 3)
