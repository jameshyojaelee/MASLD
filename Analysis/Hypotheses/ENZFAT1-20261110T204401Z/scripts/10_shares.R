#!/usr/bin/env Rscript
# ENZFAT1 S1: shares of ALT and cirrhosis heritability, from the Genomic SEM ldsc() covariance.
# Input modality: GWAS summary statistics (LD score regression). Units: GWAS participants;
# sampling covariance V comes from the 200-block genome jackknife inside ldsc().
# Shares depend on genetic correlations only, so the main call uses the observed scale.
# Usage: 10_shares.R <S0 out dir> <S1 out dir> <ld score dir>
suppressMessages({library(GenomicSEM); library(MASS); library(data.table)})
args <- commandArgs(trailingOnly = TRUE)
s0 <- args[1]; out <- args[2]; ld <- args[3]
SEED <- 20261110L; NDRAW <- 10000L
setwd(out)

traits <- c("alt_ukb", "pdff", "cirr_ghouse", "alt_mvp", "cirr_mvp",
            "nafld_finngen", "nafld_mvp", "ggt_ukb", "ast_ukb")
h2 <- fread(file.path(s0, "h2_summary.tsv"))
weak <- h2[h2_z < 4, trait]
if (length(weak)) cat("h2 z < 4, dropped:", weak, "\n")
if (any(c("alt_ukb", "pdff", "cirr_ghouse") %in% weak)) stop("a primary indicator has h2 z < 4")
traits <- setdiff(traits, weak)

run_ldsc <- function(tr, sprev, pprev, tag) {
  ldsc(traits = file.path(s0, "munged", paste0(tr, ".sumstats.gz")), sample.prev = sprev,
       population.prev = pprev, ld = ld, wld = ld, trait.names = tr,
       ldsc.log = file.path(out, paste0("ldsc_", tag)))
}
cov9 <- run_ldsc(traits, rep(NA, length(traits)), rep(NA, length(traits)), "all")
saveRDS(cov9, file.path(out, "ldsc_all.rds"))
S <- cov9$S; V <- cov9$V
dimnames(S) <- list(traits, traits)
k <- length(traits)
fwrite(data.table(trait = traits, S), file.path(out, "S_observed.tsv"), sep = "\t")
fwrite(data.table(trait = traits, cov2cor(S)), file.path(out, "rg_matrix.tsv"), sep = "\t")
fwrite(data.table(trait = traits, cov9$I), file.path(out, "ldsc_intercepts.tsv"), sep = "\t")

min_eig <- min(eigen(S, symmetric = TRUE, only.values = TRUE)$values)
cat("smallest eigenvalue of S:", min_eig, "\n")
if (min_eig <= 0) stop("S is not positive definite; stop rule: report pairwise values only")

# shares for one ALT, fat, cirrhosis trio and optional external criteria, from a correlation matrix
r2_two <- function(r1, r2, r12) (r1^2 + r2^2 - 2 * r1 * r2 * r12) / (1 - r12^2)
shares <- function(Rm, alt, fat, cir, ext) {
  raf <- Rm[alt, fat]; rac <- Rm[alt, cir]; rfc <- Rm[fat, cir]
  r2_alt <- r2_two(raf, rac, rfc)          # ALT on fat and cirrhosis
  r2_cir <- r2_two(rfc, rac, raf)          # cirrhosis on fat and ALT
  load_cir <- (rac - raf * rfc) / sqrt(1 - raf^2)   # ALT beyond fat -> cirrhosis
  o <- c(alt_shared_fat = raf^2, u_fat = 1 - raf^2, u_both = 1 - r2_alt,
         cir_fat_first = rfc^2, cir_alt_beyond_fat = load_cir^2,
         cir_alt_first = rac^2, cir_fat_beyond_alt = r2_cir - rac^2,
         cir_unique_fat = r2_cir - rac^2, cir_unique_alt = r2_cir - rfc^2,
         cir_common = rfc^2 + rac^2 - r2_cir, cir_neither = 1 - r2_cir,
         load_alt_beyond_fat_on_cir = load_cir)
  for (x in ext) {
    l <- (Rm[alt, x] - raf * Rm[fat, x]) / sqrt(1 - raf^2)
    o[paste0("load_alt_beyond_fat_on_", x)] <- l
    o[paste0("sq_load_alt_beyond_fat_on_", x)] <- l^2
  }
  o
}
ext_all <- intersect(c("nafld_finngen", "nafld_mvp", "ggt_ukb", "ast_ukb"), traits)
trios <- list(primary = c("alt_ukb", "pdff", "cirr_ghouse"), replication = c("alt_mvp", "pdff", "cirr_mvp"))
trios <- trios[sapply(trios, function(x) all(x %in% traits))]
all_shares <- function(Sm) {
  Rm <- cov2cor(Sm)
  unlist(lapply(names(trios), function(nm) {
    t3 <- trios[[nm]]
    s <- shares(Rm, t3[1], t3[2], t3[3], ext_all)
    setNames(s, paste(nm, names(s), sep = ":"))
  }))
}
est <- all_shares(S)

# draws of vech(S) from its sampling distribution (V is ordered as the lower triangle of S, by column)
set.seed(SEED)
low <- lower.tri(S, diag = TRUE)
stopifnot(nrow(V) == sum(low))
D <- mvrnorm(NDRAW, mu = S[low], Sigma = V)
draws <- matrix(NA_real_, NDRAW, length(est), dimnames = list(NULL, names(est)))
n_bad <- 0L
for (i in seq_len(NDRAW)) {
  Si <- matrix(0, k, k, dimnames = dimnames(S)); Si[low] <- D[i, ]
  Si <- Si + t(Si) - diag(diag(Si))
  if (min(eigen(Si, symmetric = TRUE, only.values = TRUE)$values) <= 0) { n_bad <- n_bad + 1L; next }
  draws[i, ] <- all_shares(Si)
}
ok <- stats::complete.cases(draws)
cat("draws dropped (matrix not positive definite):", n_bad, "of", NDRAW, "\n")
res <- data.table(quantity = names(est), estimate = est, se = apply(draws[ok, ], 2, sd),
                  lo95 = apply(draws[ok, ], 2, quantile, 0.025), hi95 = apply(draws[ok, ], 2, quantile, 0.975))
res[, c("trio", "name") := tstrsplit(quantity, ":", fixed = TRUE)]
fwrite(res[, .(trio, name, estimate, se, lo95, hi95)], file.path(out, "shares.tsv"), sep = "\t")

# tests
pick <- function(q) res[quantity == q]
u <- draws[ok, "primary:u_both"]
p_primary <- (sum(u <= 0.5) + 1) / (length(u) + 1)              # one-sided, H0: u_both <= 0.5
two_sided <- function(q) { r <- pick(q); 2 * pnorm(-abs(r$estimate / r$se)) }
p_co <- two_sided("primary:load_alt_beyond_fat_on_cir")
prim <- data.table(family = "primary", test = c("u_both > 0.5 (one-sided, draws)", "ALT-beyond-fat loading on cirrhosis != 0"),
                   estimate = c(pick("primary:u_both")$estimate, pick("primary:load_alt_beyond_fat_on_cir")$estimate),
                   p = c(p_primary, p_co))
prim[, p_holm := p.adjust(p, "holm")]
f1q <- paste0("primary:load_alt_beyond_fat_on_", ext_all)
f1 <- data.table(family = "F1", test = f1q, estimate = sapply(f1q, function(q) pick(q)$estimate),
                 p = sapply(f1q, two_sided))
f1[, p_holm := p.adjust(p, "holm")]
fwrite(rbind(prim, f1), file.path(out, "tests.tsv"), sep = "\t")

# stop rules and trio agreement
stops <- list()
for (q in c("primary:u_both", "primary:cir_neither")) {
  r <- pick(q); stops[[q]] <- data.table(check = paste(q, "more than 2 SE above zero"), value = r$estimate / r$se,
                                         pass = r$estimate / r$se > 2)
}
if ("replication" %in% names(trios)) {
  a <- pick("primary:u_both"); b <- pick("replication:u_both")
  d <- a$estimate - b$estimate
  se_indep <- sqrt(a$se^2 + b$se^2)                                   # prespecified: covariance ignored
  se_joint <- sd(draws[ok, "primary:u_both"] - draws[ok, "replication:u_both"])
  stops[["trio"]] <- data.table(check = c("trios agree on u_both within 2 SE (covariance ignored; prespecified)",
                                          "same, SE from joint draws (descriptive)"),
                                value = c(d / se_indep, d / se_joint), pass = abs(c(d / se_indep, d / se_joint)) <= 2)
}
fwrite(rbindlist(stops), file.path(out, "stop_rules.tsv"), sep = "\t")

# sanity: shares do not depend on the assumed population prevalence of cirrhosis
t3 <- trios$primary
sp <- c(NA, NA, 15225 / (15225 + 1564786))
prev <- rbindlist(lapply(c(0.005, 0.01, 0.02), function(pp) {
  cv <- run_ldsc(t3, sp, c(NA, NA, pp), paste0("prev_", pp))
  Sm <- cv$S; dimnames(Sm) <- list(t3, t3)
  s <- shares(cov2cor(Sm), t3[1], t3[2], t3[3], character(0))
  data.table(population_prevalence = pp, name = names(s), estimate = s)
}))
fwrite(dcast(prev, name ~ population_prevalence, value.var = "estimate"), file.path(out, "prevalence_check.tsv"), sep = "\t")

writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
print(res[, .(trio, name, estimate = round(estimate, 3), se = round(se, 3), lo95 = round(lo95, 3), hi95 = round(hi95, 3))])
print(fread(file.path(out, "tests.tsv"))); print(fread(file.path(out, "stop_rules.tsv")))
cat("S1 DONE\n")
