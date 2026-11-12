#!/usr/bin/env Rscript
# ENZFAT1 S1b: same draws as 10_shares.R (seed 20261110, 10,000 draws of vech(S) from V of the
# 9-trait ldsc() call in out_S1), but a draw is dropped only when the matrix of the trait set that a
# quantity uses is not positive definite (trio: 3x3; trio plus one criterion: 4x4), as in the
# prespecification's trait sets. 10_shares.R required the full 9x9 matrix to be positive definite,
# which dropped 9,279 of 10,000 draws. Unit: GWAS participants (200-block jackknife V).
# Usage: 11_shares_subpd.R <S1 out dir with ldsc_all.rds> <S1b out dir>
suppressMessages({library(MASS); library(data.table)})
args <- commandArgs(trailingOnly = TRUE)
cov9 <- readRDS(file.path(args[1], "ldsc_all.rds")); out <- args[2]
SEED <- 20261110L; NDRAW <- 10000L
S <- cov9$S; V <- cov9$V; traits <- colnames(S); dimnames(S) <- list(traits, traits); k <- length(traits)
pd <- function(M) min(eigen(M, symmetric = TRUE, only.values = TRUE)$values) > 0
r2_two <- function(r1, r2, r12) (r1^2 + r2^2 - 2 * r1 * r2 * r12) / (1 - r12^2)
trio_shares <- function(Rm, alt, fat, cir) {
  raf <- Rm[alt, fat]; rac <- Rm[alt, cir]; rfc <- Rm[fat, cir]
  r2_alt <- r2_two(raf, rac, rfc); r2_cir <- r2_two(rfc, rac, raf); l <- (rac - raf * rfc) / sqrt(1 - raf^2)
  c(alt_shared_fat = raf^2, u_fat = 1 - raf^2, u_both = 1 - r2_alt, cir_fat_first = rfc^2,
    cir_alt_beyond_fat = l^2, cir_alt_first = rac^2, cir_fat_beyond_alt = r2_cir - rac^2,
    cir_unique_fat = r2_cir - rac^2, cir_unique_alt = r2_cir - rfc^2, cir_common = rfc^2 + rac^2 - r2_cir,
    cir_neither = 1 - r2_cir, load_alt_beyond_fat_on_cir = l)
}
ext_load <- function(Rm, alt, fat, x) {
  raf <- Rm[alt, fat]; l <- (Rm[alt, x] - raf * Rm[fat, x]) / sqrt(1 - raf^2)
  setNames(c(l, l^2), paste0(c("load_alt_beyond_fat_on_", "sq_load_alt_beyond_fat_on_"), x))
}
ext <- intersect(c("nafld_finngen", "nafld_mvp", "ggt_ukb", "ast_ukb"), traits)
trios <- list(primary = c("alt_ukb", "pdff", "cirr_ghouse"), replication = c("alt_mvp", "pdff", "cirr_mvp"))
# one function returns every quantity, with NA where that quantity's trait-set matrix is not PD
all_q <- function(Sm, check = TRUE) {
  Rm <- cov2cor(Sm); o <- c()
  for (nm in names(trios)) {
    t3 <- trios[[nm]]
    s <- trio_shares(Rm, t3[1], t3[2], t3[3])
    if (check && !pd(Sm[t3, t3])) s[] <- NA
    o <- c(o, setNames(s, paste(nm, names(s), sep = ":")))
    for (x in ext) {
      e <- ext_load(Rm, t3[1], t3[2], x); s4 <- c(t3, x)
      if (check && !pd(Sm[s4, s4])) e[] <- NA
      o <- c(o, setNames(e, paste(nm, names(e), sep = ":")))
    }
  }
  o
}
est <- all_q(S, check = FALSE)
set.seed(SEED)
low <- lower.tri(S, diag = TRUE)
D <- mvrnorm(NDRAW, mu = S[low], Sigma = V)
draws <- matrix(NA_real_, NDRAW, length(est), dimnames = list(NULL, names(est)))
pd9 <- logical(NDRAW)
for (i in seq_len(NDRAW)) {
  Si <- matrix(0, k, k, dimnames = dimnames(S)); Si[low] <- D[i, ]; Si <- Si + t(Si) - diag(diag(Si))
  pd9[i] <- pd(Si); draws[i, ] <- all_q(Si)
}
nvalid <- colSums(is.finite(draws))
cat("draws with full 9x9 matrix positive definite:", sum(pd9), "of", NDRAW, "\n")
q <- function(f, ...) apply(draws, 2, function(v) f(v[is.finite(v)], ...))
res <- data.table(quantity = names(est), estimate = est, se = q(sd), lo95 = q(quantile, 0.025),
                  hi95 = q(quantile, 0.975), valid_draws = nvalid)
res[, c("trio", "name") := tstrsplit(quantity, ":", fixed = TRUE)]
fwrite(res[, .(trio, name, estimate, se, lo95, hi95, valid_draws)], file.path(out, "shares.tsv"), sep = "\t")
pick <- function(x) res[quantity == x]
u <- draws[, "primary:u_both"]; u <- u[is.finite(u)]
p_primary <- (sum(u <= 0.5) + 1) / (length(u) + 1)
two_sided <- function(x) { r <- pick(x); 2 * pnorm(-abs(r$estimate / r$se)) }
prim <- data.table(family = "primary", test = c("u_both > 0.5 (one-sided, draws)", "ALT-beyond-fat loading on cirrhosis != 0"),
                   estimate = c(pick("primary:u_both")$estimate, pick("primary:load_alt_beyond_fat_on_cir")$estimate),
                   p = c(p_primary, two_sided("primary:load_alt_beyond_fat_on_cir")))
prim[, p_holm := p.adjust(p, "holm")]
f1q <- paste0("primary:load_alt_beyond_fat_on_", ext)
f1 <- data.table(family = "F1", test = f1q, estimate = sapply(f1q, function(x) pick(x)$estimate), p = sapply(f1q, two_sided))
f1[, p_holm := p.adjust(p, "holm")]
fwrite(rbind(prim, f1), file.path(out, "tests.tsv"), sep = "\t")
stops <- rbindlist(lapply(c("primary:u_both", "primary:cir_neither", "replication:u_both", "replication:cir_neither"), function(x) {
  r <- pick(x); data.table(check = paste(x, "more than 2 SE above zero"), value = r$estimate / r$se, pass = r$estimate / r$se > 2)
}))
a <- pick("primary:u_both"); b <- pick("replication:u_both"); dd <- a$estimate - b$estimate
both <- is.finite(draws[, "primary:u_both"]) & is.finite(draws[, "replication:u_both"])
se_joint <- sd(draws[both, "primary:u_both"] - draws[both, "replication:u_both"])
stops <- rbind(stops, data.table(check = c("trios agree on u_both within 2 SE (covariance ignored; prespecified)",
                                           "same, SE from joint draws (descriptive)"),
                                 value = c(dd / sqrt(a$se^2 + b$se^2), dd / se_joint),
                                 pass = abs(c(dd / sqrt(a$se^2 + b$se^2), dd / se_joint)) <= 2))
fwrite(stops, file.path(out, "stop_rules.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
print(res[, .(trio, name, estimate = round(estimate, 3), se = round(se, 3), lo95 = round(lo95, 3), hi95 = round(hi95, 3), valid_draws)])
print(fread(file.path(out, "tests.tsv"))); print(stops)
cat("S1b DONE\n")
