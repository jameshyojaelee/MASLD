# V5 independent recompute of the primary statistic from saved tables (V4 verification pattern).
# Uses lm() + hatvalues for HC3 and refits Freedman-Lane with lm(), not the closed-form helpers
# of 04_outcome.R. Inputs: outputs/job3/donor_scores.tsv, outputs/donor_qc.tsv,
# outputs/job3/permutation_matrix.tsv.gz, outputs/job3/results_primary.tsv.
suppressPackageStartupMessages(library(data.table))
RUN <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/V5-20261021T122723Z"
J3 <- file.path(RUN, "outputs/job3")
out <- file.path(J3, "verification_recompute.tsv")
if (file.exists(out)) stop("refusing to overwrite")
ds <- fread(file.path(J3, "donor_scores.tsv"))
dq <- fread(file.path(RUN, "outputs/donor_qc.tsv"), colClasses = list(character = "id"))
d <- merge(ds[, .(gsm, ductular)], dq[, .(gsm, disease, sex_f, n_hep)], by = "gsm")
setorder(d, gsm)
P <- as.matrix(fread(cmd = paste("zcat", shQuote(file.path(J3, "permutation_matrix.tsv.gz"))), header = FALSE))
res <- fread(file.path(J3, "results_primary.tsv"))[1]

hc3_t <- function(fit, term) {
  X <- model.matrix(fit); e <- residuals(fit); h <- hatvalues(fit)
  B <- solve(crossprod(X)) %*% t(X)
  V <- B %*% diag(e^2 / (1 - h)^2) %*% t(B)
  unname(coef(fit)[term] / sqrt(V[term, term]))
}
d[, lnh := log10(n_hep)]
fit <- lm(ductular ~ disease + sex_f + lnh, data = d)
t_obs <- hc3_t(fit, "disease")
red <- lm(ductular ~ sex_f + lnh, data = d)
f0 <- fitted(red); r0 <- residuals(red)
tnull <- apply(P, 1, function(p) {
  d$yp <- f0 + r0[p]
  hc3_t(lm(yp ~ disease + sex_f + lnh, data = d), "disease")
})
perm_p <- (1 + sum(abs(tnull) >= abs(t_obs))) / (length(tnull) + 1)
v <- data.table(quantity = c("beta", "hc3_t", "perm_p"),
                original = c(res$beta, res$hc3_t, res$perm_p),
                recomputed = c(unname(coef(fit)["disease"]), t_obs, perm_p))
v[, abs_diff := abs(original - recomputed)][, match := abs_diff < c(1e-8, 1e-8, 1e-12)]
fwrite(v, out, sep = "\t")
print(v)
