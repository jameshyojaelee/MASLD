#!/usr/bin/env Rscript
# ENZFAT1 S1b (2026-11-18 update): same draws as 10_shares.R (seed 20261110, 10,000 draws of
# vech(S) from V of the 9-trait ldsc() call in out_S1), fixing the s4 -> s3 bug on external
# criteria so (alt, fat, x) checks 3x3 PD rather than 4x4 with cirrhosis, and adding:
#   - analytical numDeriv::jacobian delta-method SEs for (u_alt_fat, u_alt_cir, u_alt_joint)
#     and all share quantities alongside MVN draws,
#   - explicit ldsc_intercepts summary table from cov9$I,
#   - clear labeling distinguishing 2-trait u_alt_fat (0.786) from 3-trait u_alt_joint (0.701).
# Usage: 11_shares_subpd_20261118.R <S1 out dir with ldsc_all.rds> <S1b out dir>
suppressMessages({library(MASS); library(data.table); library(numDeriv)})
args <- commandArgs(trailingOnly = TRUE)
s1_dir <- args[1]
out <- args[2]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
cov9 <- readRDS(file.path(s1_dir, "ldsc_all.rds"))
SEED <- 20261110L; NDRAW <- 10000L
S <- cov9$S; V <- cov9$V; I_mat <- cov9$I
traits <- colnames(S); dimnames(S) <- list(traits, traits); dimnames(I_mat) <- list(traits, traits)
k <- length(traits)

pd <- function(M) min(eigen(M, symmetric = TRUE, only.values = TRUE)$values) > 0
r2_two <- function(r1, r2, r12) (r1^2 + r2^2 - 2 * r1 * r2 * r12) / (1 - r12^2)

trio_shares <- function(Rm, alt, fat, cir) {
  raf <- Rm[alt, fat]; rac <- Rm[alt, cir]; rfc <- Rm[fat, cir]
  r2_alt <- r2_two(raf, rac, rfc); r2_cir <- r2_two(rfc, rac, raf); l <- (rac - raf * rfc) / sqrt(1 - raf^2)
  c(alt_shared_fat = raf^2, u_fat = 1 - raf^2, u_cir = 1 - rac^2, u_both = 1 - r2_alt,
    cir_fat_first = rfc^2, cir_alt_beyond_fat = l^2, cir_alt_first = rac^2,
    cir_fat_beyond_alt = r2_cir - rac^2, cir_unique_fat = r2_cir - rac^2,
    cir_unique_alt = r2_cir - rfc^2, cir_common = rfc^2 + rac^2 - r2_cir,
    cir_neither = 1 - r2_cir, load_alt_beyond_fat_on_cir = l)
}

ext_load <- function(Rm, alt, fat, x) {
  raf <- Rm[alt, fat]; l <- (Rm[alt, x] - raf * Rm[fat, x]) / sqrt(1 - raf^2)
  setNames(c(l, l^2), paste0(c("load_alt_beyond_fat_on_", "sq_load_alt_beyond_fat_on_"), x))
}

ext <- intersect(c("nafld_finngen", "nafld_mvp", "ggt_ukb", "ast_ukb"), traits)
trios <- list(primary = c("alt_ukb", "pdff", "cirr_ghouse"),
              replication = c("alt_mvp", "pdff", "cirr_mvp"))

# Returns every quantity, with NA where that quantity's trait-set matrix is not PD
all_q <- function(Sm, check = TRUE) {
  Rm <- cov2cor(Sm); o <- c()
  for (nm in names(trios)) {
    t3 <- trios[[nm]]
    s <- trio_shares(Rm, t3[1], t3[2], t3[3])
    if (check) {
      if (!pd(Sm[t3, t3])) s[] <- NA_real_
      # 2-trait submatrices for u_fat (alt, fat) and u_cir (alt, cir)
      if (!pd(Sm[t3[c(1, 2)], t3[c(1, 2)]])) s[c("alt_shared_fat", "u_fat")] <- NA_real_ else if (is.na(s["u_fat"])) {
        s_full <- trio_shares(Rm, t3[1], t3[2], t3[3])
        s[c("alt_shared_fat", "u_fat")] <- s_full[c("alt_shared_fat", "u_fat")]
      }
      if (!pd(Sm[t3[c(1, 3)], t3[c(1, 3)]])) s["u_cir"] <- NA_real_ else if (is.na(s["u_cir"])) {
        s_full <- trio_shares(Rm, t3[1], t3[2], t3[3])
        s["u_cir"] <- s_full["u_cir"]
      }
    }
    o <- c(o, setNames(s, paste(nm, names(s), sep = ":")))
    for (x in ext) {
      e <- ext_load(Rm, t3[1], t3[2], x); s3 <- c(t3[1], t3[2], x)
      if (check && !pd(Sm[s3, s3])) e[] <- NA_real_
      o <- c(o, setNames(e, paste(nm, names(e), sep = ":")))
    }
  }
  o
}

# Explicit 3-quantity extractor (u_alt_fat, u_alt_cir, u_alt_joint) for clean delta-method labeling
u_key_fn <- function(Sm) {
  Rm <- cov2cor(Sm); o <- c()
  for (nm in names(trios)) {
    t3 <- trios[[nm]]
    raf <- Rm[t3[1], t3[2]]; rac <- Rm[t3[1], t3[3]]; rfc <- Rm[t3[2], t3[3]]
    r2_alt <- r2_two(raf, rac, rfc)
    v <- c(u_alt_fat = 1 - raf^2, u_alt_cir = 1 - rac^2, u_alt_joint = 1 - r2_alt)
    o <- c(o, setNames(v, paste(nm, names(v), sep = ":")))
  }
  o
}

est <- all_q(S, check = FALSE)
low <- lower.tri(S, diag = TRUE)
stopifnot(nrow(V) == sum(low))
s_vec <- S[low]
vec_to_S <- function(v) {
  Sm <- matrix(0, k, k, dimnames = dimnames(S))
  Sm[low] <- v
  Sm + t(Sm) - diag(diag(Sm))
}

# Analytical delta-method SEs via numDeriv::jacobian on vech(S) with V
J_all <- numDeriv::jacobian(func = function(v) all_q(vec_to_S(v), check = FALSE), x = s_vec)
se_delta_all <- sqrt(pmax(0, diag(J_all %*% V %*% t(J_all))))
names(se_delta_all) <- names(est)

u_est <- u_key_fn(S)
J_u <- numDeriv::jacobian(func = function(v) u_key_fn(vec_to_S(v)), x = s_vec)
se_delta_u <- sqrt(pmax(0, diag(J_u %*% V %*% t(J_u))))
names(se_delta_u) <- names(u_est)

# MVN draws from sampling covariance V
set.seed(SEED)
D <- mvrnorm(NDRAW, mu = S[low], Sigma = V)
draws <- matrix(NA_real_, NDRAW, length(est), dimnames = list(NULL, names(est)))
pd9 <- logical(NDRAW)
for (i in seq_len(NDRAW)) {
  Si <- vec_to_S(D[i, ])
  pd9[i] <- pd(Si); draws[i, ] <- all_q(Si, check = TRUE)
}
nvalid <- colSums(is.finite(draws))
cat("draws with full 9x9 matrix positive definite:", sum(pd9), "of", NDRAW, "\n")

q <- function(f, ...) apply(draws, 2, function(v) f(v[is.finite(v)], ...))
res <- data.table(quantity = names(est), estimate = est, se = q(sd), lo95 = q(quantile, 0.025),
                  hi95 = q(quantile, 0.975), se_delta = se_delta_all, valid_draws = nvalid)
res[, c("trio", "name") := tstrsplit(quantity, ":", fixed = TRUE)]
res[, model_scope := fifelse(name %in% c("alt_shared_fat", "u_fat"), "2-trait (alt, fat)",
                     fifelse(name == "u_cir", "2-trait (alt, cir)",
                     fifelse(grepl("^load_alt_beyond_fat_on_(nafld|ggt|ast)|^sq_load_alt_beyond_fat_on_", name),
                             "3-trait (alt, fat, criterion)", "3-trait (alt, fat, cir)")))]
fwrite(res[, .(trio, name, model_scope, estimate, se, lo95, hi95, se_delta, valid_draws)],
       file.path(out, "shares.tsv"), sep = "\t")

# Explicit summary table for (u_alt_fat, u_alt_cir, u_alt_joint): 2-trait 0.786 vs 3-trait 0.701
alias_map <- c(u_alt_fat = "u_fat", u_alt_cir = "u_cir", u_alt_joint = "u_both")
scope_map <- c(u_alt_fat = "2-trait (ALT, liver fat)",
               u_alt_cir = "2-trait (ALT, cirrhosis)",
               u_alt_joint = "3-trait (ALT, liver fat, cirrhosis)")
desc_map <- c(u_alt_fat = "2-trait share: 1 - rg(ALT, liver fat)^2 [primary=0.786]",
              u_alt_cir = "2-trait share: 1 - rg(ALT, cirrhosis)^2",
              u_alt_joint = "3-trait share: 1 - R^2(ALT ~ liver fat + cirrhosis) [primary=0.701]")
delta_tab <- rbindlist(lapply(names(u_est), function(key) {
  parts <- strsplit(key, ":", fixed = TRUE)[[1]]
  tr <- parts[1]; met <- parts[2]; al <- alias_map[[met]]
  r <- res[trio == tr & name == al]
  sed <- se_delta_u[[key]]; estv <- u_est[[key]]
  data.table(
    trio = tr,
    metric = met,
    alias_in_shares = al,
    model_scope = scope_map[[met]],
    label = desc_map[[met]],
    estimate = estv,
    se_mvn = r$se,
    lo95_mvn = r$lo95,
    hi95_mvn = r$hi95,
    se_delta = sed,
    lo95_delta = estv - qnorm(0.975) * sed,
    hi95_delta = estv + qnorm(0.975) * sed,
    rel_se_diff = (sed - r$se) / r$se,
    valid_draws = r$valid_draws
  )
}))
fwrite(delta_tab, file.path(out, "delta_method_se.tsv"), sep = "\t")

# Explicit ldsc_intercepts matrix + pairwise summary table
fwrite(data.table(trait = traits, I_mat), file.path(out, "ldsc_intercepts.tsv"), sep = "\t")
Rg <- cov2cor(S)
int_rows <- rbindlist(lapply(seq_len(k), function(i) {
  rbindlist(lapply(seq(i, k), function(j) {
    data.table(
      trait1 = traits[i],
      trait2 = traits[j],
      intercept_type = if (i == j) "univariate_diagonal" else "cross_trait_overlap",
      ldsc_intercept = I_mat[i, j],
      genetic_cov_S = S[i, j],
      genetic_corr_rg = Rg[i, j]
    )
  }))
}))
int_rows <- int_rows[order(intercept_type != "univariate_diagonal", -abs(ldsc_intercept))]
fwrite(int_rows, file.path(out, "ldsc_intercepts_summary.tsv"), sep = "\t")

# Hypothesis tests (using MVN draws and delta-method z comparison)
pick <- function(x) res[quantity == x]
u <- draws[, "primary:u_both"]; u <- u[is.finite(u)]
p_primary <- (sum(u <= 0.5) + 1) / (length(u) + 1)
two_sided <- function(x) { r <- pick(x); 2 * pnorm(-abs(r$estimate / r$se)) }
two_sided_delta <- function(x) { r <- pick(x); 2 * pnorm(-abs(r$estimate / r$se_delta)) }
prim <- data.table(
  family = "primary",
  test = c("u_both (3-trait u_alt_joint) > 0.5 (one-sided, draws)",
           "ALT-beyond-fat loading on cirrhosis != 0"),
  estimate = c(pick("primary:u_both")$estimate, pick("primary:load_alt_beyond_fat_on_cir")$estimate),
  p = c(p_primary, two_sided("primary:load_alt_beyond_fat_on_cir")),
  p_delta = c(pnorm((pick("primary:u_both")$estimate - 0.5) / pick("primary:u_both")$se_delta, lower.tail = FALSE),
              two_sided_delta("primary:load_alt_beyond_fat_on_cir"))
)
prim[, p_holm := p.adjust(p, "holm")]
f1q <- paste0("primary:load_alt_beyond_fat_on_", ext)
f1 <- data.table(
  family = "F1",
  test = f1q,
  estimate = sapply(f1q, function(x) pick(x)$estimate),
  p = sapply(f1q, two_sided),
  p_delta = sapply(f1q, two_sided_delta)
)
f1[, p_holm := p.adjust(p, "holm")]
fwrite(rbind(prim, f1), file.path(out, "tests.tsv"), sep = "\t")

# Stop rules and trio agreement
stops <- rbindlist(lapply(c("primary:u_both", "primary:cir_neither", "replication:u_both", "replication:cir_neither"), function(x) {
  r <- pick(x); data.table(check = paste(x, "more than 2 SE above zero"), value = r$estimate / r$se, pass = r$estimate / r$se > 2)
}))
a <- pick("primary:u_both"); b <- pick("replication:u_both"); dd <- a$estimate - b$estimate
both <- is.finite(draws[, "primary:u_both"]) & is.finite(draws[, "replication:u_both"])
se_joint <- sd(draws[both, "primary:u_both"] - draws[both, "replication:u_both"])
stops <- rbind(stops, data.table(
  check = c("trios agree on u_both within 2 SE (covariance ignored; prespecified)",
            "same, SE from joint draws (descriptive)"),
  value = c(dd / sqrt(a$se^2 + b$se^2), dd / se_joint),
  pass = abs(c(dd / sqrt(a$se^2 + b$se^2), dd / se_joint)) <= 2
))
fwrite(stops, file.path(out, "stop_rules.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))

cat("\n--- 2-trait vs 3-trait Delta-Method vs MVN SEs ---\n")
print(delta_tab[, .(trio, metric, model_scope, estimate = round(estimate, 4),
                    se_mvn = round(se_mvn, 4), se_delta = round(se_delta, 4),
                    lo95_mvn = round(lo95_mvn, 4), hi95_mvn = round(hi95_mvn, 4),
                    valid_draws)])
cat("\n--- All Shares (including external criteria 3x3 PD valid_draws) ---\n")
print(res[, .(trio, name, model_scope, estimate = round(estimate, 4), se = round(se, 4),
              se_delta = round(se_delta, 4), lo95 = round(lo95, 4), hi95 = round(hi95, 4), valid_draws)])
cat("\n--- LDSC Intercepts (diagonal & top cross-trait overlaps) ---\n")
print(int_rows[intercept_type == "univariate_diagonal" | abs(ldsc_intercept) >= 0.05,
               .(trait1, trait2, intercept_type, ldsc_intercept = round(ldsc_intercept, 4),
                 genetic_corr_rg = round(genetic_corr_rg, 4))])
cat("\n--- Tests & Stop Rules ---\n")
print(fread(file.path(out, "tests.tsv"))); print(stops)
cat("S1b_20261118 DONE\n")
