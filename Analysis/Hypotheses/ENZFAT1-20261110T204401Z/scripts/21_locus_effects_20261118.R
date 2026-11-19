#!/usr/bin/env Rscript
# ENZFAT1 S2 step 2 (2026-11-18 remediation, N14):
# - Aligns ALT-defined loci to Pan-UKBB ALT-raising allele and AST-defined loci cleanly to enzyme_raising_allele
#   (and ctrl_GOT1_lead to Pan-UKBB AST-raising allele, other controls to Pan-UKBB ALT-raising allele).
# - Replaces ashr::ash(..., fixg = TRUE) control evaluation with exact log-space halfuniform Gaussian convolution
#   (eval_halfuniform_ash_logspace) so extreme-|z| controls never underflow to NaN or lfsr=1/pm=0.
# Usage: 21_locus_effects_20261118.R <locus_rows.tsv> <loci.tsv> <out dir>
suppressMessages({library(data.table); library(ashr)})
args <- commandArgs(trailingOnly = TRUE)
rows <- fread(args[1]); loci <- fread(args[2]); out <- args[3]
set.seed(20261110)
rows[, locus := as.character(locus)]; loci[, locus := as.character(locus)]

# Exact log-space halfuniform Gaussian-convolution posterior evaluator for fixed prior g (unimix).
eval_halfuniform_ash_logspace <- function(betahat, sebetahat, g) {
  stopifnot(length(betahat) == length(sebetahat))
  n <- length(betahat)
  pm <- rep(NA_real_, n); psd <- rep(NA_real_, n)
  lfsr <- rep(NA_real_, n); lfdr <- rep(NA_real_, n)

  gl_x <- c(
    -0.9894009349916499326, -0.9445750230732325761, -0.8656312023878317439, -0.7554044083550030339,
    -0.6178762444026437484, -0.4580167776572273863, -0.2816035507792589132, -0.0950125098376374402,
     0.0950125098376374402,  0.2816035507792589132,  0.4580167776572273863,  0.6178762444026437484,
     0.7554044083550030339,  0.8656312023878317439,  0.9445750230732325761,  0.9894009349916499326
  )
  gl_w <- c(
    0.0271524594117540949, 0.0622535239386478929, 0.0951585116824927848, 0.1246289712555338721,
    0.1495959888165767321, 0.1691565193950025382, 0.1826034150449235889, 0.1894506104550684963,
    0.1894506104550684963, 0.1826034150449235889, 0.1691565193950025382, 0.1495959888165767321,
    0.1246289712555338721, 0.0951585116824927848, 0.0622535239386478929, 0.0271524594117540949
  )

  log1mexp <- function(d) {
    if (d >= 0) return(-Inf)
    if (d > -log(2)) log(-expm1(d)) else log1p(-exp(d))
  }

  log_pnorm_diff <- function(l, u) {
    if (u <= 0) {
      lu <- pnorm(u, lower.tail = TRUE, log.p = TRUE)
      ll <- pnorm(l, lower.tail = TRUE, log.p = TRUE)
      return(lu + log1mexp(ll - lu))
    } else if (l >= 0) {
      sl <- pnorm(l, lower.tail = FALSE, log.p = TRUE)
      su <- pnorm(u, lower.tail = FALSE, log.p = TRUE)
      return(sl + log1mexp(su - sl))
    } else {
      if (abs(l) < 1e-4 && abs(u) < 1e-4) {
        du <- u - l; d3 <- u^3 - l^3; d5 <- u^5 - l^5; d7 <- u^7 - l^7
        return(log((du - d3 / 6 + d5 / 40 - d7 / 336) / sqrt(2 * pi)))
      }
      return(log(pnorm(u) - pnorm(l)))
    }
  }

  tail_shifted_moments <- function(z, w) {
    M <- min(z * w, 45)
    edges <- unique(c(0, pmin(M, c(0.5, 2, 6, 15, 45))))
    edges <- edges[edges <= M]
    if (tail(edges, 1) < M) edges <- c(edges, M)
    vs <- numeric(0); ws <- numeric(0)
    for (j in seq_len(length(edges) - 1L)) {
      e0 <- edges[j]; e1 <- edges[j + 1L]
      if (e1 <= e0) next
      half <- 0.5 * (e1 - e0); mid <- 0.5 * (e1 + e0)
      vs <- c(vs, mid + half * gl_x)
      ws <- c(ws, half * gl_w)
    }
    fv <- ws * exp(-vs - (vs^2) / (2 * z^2))
    I0 <- sum(fv)
    ev <- sum(vs * fv) / I0
    vv <- max(0, sum(((vs - ev)^2) * fv) / I0)
    list(mean = ev / z, var = vv / (z^2))
  }

  truncnorm_moments <- function(l, u) {
    w <- u - l
    if (w < 1e-3) {
      half <- 0.5 * w; mid <- 0.5 * (l + u)
      ts <- mid + half * gl_x
      ld <- -0.5 * ts^2
      fv <- (half * gl_w) * exp(ld - max(ld))
      I0 <- sum(fv)
      mu <- sum(ts * fv) / I0
      vr <- max(0, sum(((ts - mu)^2) * fv) / I0)
      return(list(mean = mu, var = vr))
    }
    if (u < -8) {
      tm <- tail_shifted_moments(-u, w)
      return(list(mean = u - tm$mean, var = tm$var))
    }
    if (l > 8) {
      tm <- tail_shifted_moments(l, w)
      return(list(mean = l + tm$mean, var = tm$var))
    }
    D <- log_pnorm_diff(l, u)
    rl <- exp(dnorm(l, log = TRUE) - D)
    ru <- exp(dnorm(u, log = TRUE) - D)
    mu <- rl - ru
    vr <- max(0, 1 + (if (rl == 0) 0 else l * rl) - (if (ru == 0) 0 else u * ru) - mu^2)
    list(mean = mu, var = vr)
  }

  idx <- which(g$pi > 0)
  pi_k <- g$pi[idx]; a_k <- g$a[idx]; b_k <- g$b[idx]
  log_pi <- log(pi_k)
  K <- length(idx)

  for (i in seq_len(n)) {
    x <- betahat[i]; s <- sebetahat[i]
    if (!is.finite(x) || !is.finite(s) || s <= 0) next
    log_lik <- numeric(K); mu_k <- numeric(K); var_k <- numeric(K)
    p0_k <- numeric(K); pp_k <- numeric(K); pn_k <- numeric(K)
    for (k in seq_len(K)) {
      ak <- a_k[k]; bk <- b_k[k]
      if (ak == bk) {
        log_lik[k] <- dnorm(x, mean = ak, sd = s, log = TRUE)
        mu_k[k] <- ak; var_k[k] <- 0
        p0_k[k] <- as.numeric(ak == 0)
        pp_k[k] <- as.numeric(ak > 0)
        pn_k[k] <- as.numeric(ak < 0)
      } else {
        l <- (ak - x) / s; u <- (bk - x) / s
        D <- log_pnorm_diff(l, u)
        log_lik[k] <- D - log(bk - ak)
        tm <- truncnorm_moments(l, u)
        mu_k[k] <- x + s * tm$mean
        var_k[k] <- (s^2) * tm$var
        p0_k[k] <- 0
        if (ak >= 0) {
          pp_k[k] <- 1; pn_k[k] <- 0
        } else if (bk <= 0) {
          pp_k[k] <- 0; pn_k[k] <- 1
        } else {
          z0 <- -x / s
          pp_k[k] <- min(1, max(0, exp(log_pnorm_diff(z0, u) - D)))
          pn_k[k] <- min(1, max(0, exp(log_pnorm_diff(l, z0) - D)))
        }
      }
    }
    lw <- log_pi + log_lik
    m <- max(lw)
    w_post <- exp(lw - m)
    w_post <- w_post / sum(w_post)
    pm_i <- sum(w_post * mu_k)
    psd_i <- sqrt(max(0, sum(w_post * (var_k + (mu_k - pm_i)^2))))
    p0 <- sum(w_post * p0_k); pp <- sum(w_post * pp_k); pn <- sum(w_post * pn_k)
    pm[i] <- pm_i; psd[i] <- psd_i; lfdr[i] <- p0; lfsr[i] <- min(p0 + pn, p0 + pp)
  }
  list(pm = pm, psd = psd, lfsr = lfsr, lfdr = lfdr)
}

# Alignment: ALT-defined loci -> Pan-UKBB ALT-raising allele; AST-defined loci -> enzyme_raising_allele
# (verified against Pan-UKBB AST); controls -> AST-raising for ctrl_GOT1_lead, ALT-raising otherwise.
ref_alt <- rows[trait == "alt_ukb", .(locus, alt_ukb_raise = ifelse(beta > 0, ea, oa))]
ref_ast <- rows[trait == "ast_ukb", .(locus, ast_ukb_raise = ifelse(beta > 0, ea, oa))]

all_loci <- unique(rows$locus)
ref <- data.table(locus = all_loci)
ref <- merge(ref, loci[, .(locus, defining_enzyme, enzyme_raising_allele)], by = "locus", all.x = TRUE)
ref <- merge(ref, ref_alt, by = "locus", all.x = TRUE)
ref <- merge(ref, ref_ast, by = "locus", all.x = TRUE)
ref[, ref_allele := fcase(
  defining_enzyme == "MVP_AST_EUR", enzyme_raising_allele,
  defining_enzyme == "MVP_ALT_EUR", alt_ukb_raise,
  grepl("GOT1", locus), ast_ukb_raise,
  default = alt_ukb_raise
)]
ref <- ref[!is.na(ref_allele)]

cat("AST-defined loci alignment check:\n")
print(ref[defining_enzyme == "MVP_AST_EUR", .(locus, defining_enzyme, enzyme_raising_allele, ast_ukb_raise, alt_ukb_raise, ref_allele)])

d <- merge(rows, ref[, .(locus, ref_allele, defining_enzyme)], by = "locus")
d[, beta_al := ifelse(ea == ref_allele, beta, -beta)]
stopifnot(d[trait == "alt_ukb" & defining_enzyme == "MVP_ALT_EUR", all(beta_al > 0)])
stopifnot(d[trait == "ast_ukb" & grepl("GOT1", locus), all(beta_al > 0)])

wide <- dcast(d, locus + ref_allele ~ trait, value.var = c("beta_al", "se"))
info <- loci[, .(locus, chr, pos, genes, cls_DIG1 = cls, defining_enzyme, enzyme_raising_allele)]
wide <- merge(info, wide, by = "locus", all.y = TRUE)
is_locus <- !grepl("^ctrl_", wide$locus)

# adaptive shrinkage per outcome over the lead variants of the locus table (controls excluded from the fit)
for (tr in c("pdff", "cirr_ghouse", "nafld_finngen")) {
  b <- wide[[paste0("beta_al_", tr)]]; s <- wide[[paste0("se_", tr)]]
  use <- is_locus & is.finite(b) & is.finite(s)
  fit <- ash(b[use], s[use], mixcompdist = "halfuniform")
  wide[use, paste0(c("postmean_", "postsd_", "lfsr_", "lfdr_"), tr) :=
         .(get_pm(fit), get_psd(fit), get_lfsr(fit), get_lfdr(fit))]

  # Parity check of exact log-space evaluator against ashr on moderate-|z| fitted loci
  mod <- use & (abs(b / s) < 5)
  if (any(mod)) {
    chk_log <- eval_halfuniform_ash_logspace(b[mod], s[mod], get_fitted_g(fit))
    chk_ash <- ash(b[mod], s[mod], g = get_fitted_g(fit), fixg = TRUE, mixcompdist = "halfuniform")
    d_pm <- max(abs(chk_log$pm - get_pm(chk_ash)))
    d_psd <- max(abs(chk_log$psd - get_psd(chk_ash)))
    d_lfsr <- max(abs(chk_log$lfsr - get_lfsr(chk_ash)))
    d_lfdr <- max(abs(chk_log$lfdr - get_lfdr(chk_ash)))
    cat(sprintf("%s log-space parity (n=%d, |z|<5): max|d_pm|=%.2e max|d_psd|=%.2e max|d_lfsr|=%.2e max|d_lfdr|=%.2e\n",
                tr, sum(mod), d_pm, d_psd, d_lfsr, d_lfdr))
    stopifnot(d_pm < 1e-8, d_psd < 1e-8, d_lfsr < 1e-8, d_lfdr < 1e-8)
  }

  # Control rows outside the locus table: exact log-space posterior under fixed prior g
  uc <- !is_locus & is.finite(b) & is.finite(s)
  if (any(uc)) {
    fc <- eval_halfuniform_ash_logspace(b[uc], s[uc], get_fitted_g(fit))
    stopifnot(all(is.finite(fc$pm)), all(is.finite(fc$psd)), all(is.finite(fc$lfsr)), all(is.finite(fc$lfdr)))
    wide[uc, paste0(c("postmean_", "postsd_", "lfsr_", "lfdr_"), tr) :=
           .(fc$pm, fc$psd, fc$lfsr, fc$lfdr)]
  }
  cat(tr, ": loci fitted", sum(use), "; pi0", round(get_pi0(fit), 3), "; lfsr < 0.05:",
      sum(get_lfsr(fit) < 0.05), "\n")
}

ratio <- function(f, sf, a, sa) list(f / a, sqrt(sf^2 / a^2 + f^2 * sa^2 / a^4))
wide[, c("fat_per_alt", "fat_per_alt_se") := ratio(beta_al_pdff, se_pdff, beta_al_alt_ukb, se_alt_ukb)]
wide[, c("cirr_per_alt", "cirr_per_alt_se") := ratio(beta_al_cirr_ghouse, se_cirr_ghouse, beta_al_alt_ukb, se_alt_ukb)]
wide[, c("fat_per_alt_mvp", "fat_per_alt_mvp_se") := ratio(beta_al_pdff, se_pdff, beta_al_alt_mvp, se_alt_mvp)]
wide[, enz_ukb_beta := fifelse(defining_enzyme == "MVP_AST_EUR", beta_al_ast_ukb, beta_al_alt_ukb)]
wide[, enz_ukb_se := fifelse(defining_enzyme == "MVP_AST_EUR", se_ast_ukb, se_alt_ukb)]
wide[, c("fat_per_enzyme", "fat_per_enzyme_se") := ratio(beta_al_pdff, se_pdff, enz_ukb_beta, enz_ukb_se)]
wide[, c("cirr_per_enzyme", "cirr_per_enzyme_se") := ratio(beta_al_cirr_ghouse, se_cirr_ghouse, enz_ukb_beta, enz_ukb_se)]
for (tr in c("pdff", "cirr_ghouse", "nafld_finngen", "nafld_mvp", "cirr_mvp", "albumin_mvp", "platelet_mvp", "alt_mvp", "ast_mvp")) {
  if (paste0("beta_al_", tr) %in% names(wide))
    wide[, paste0("z_", tr) := get(paste0("beta_al_", tr)) / get(paste0("se_", tr))]
}
fwrite(wide, file.path(out, "locus_effects.tsv"), sep = "\t")

# controls
has <- function(g) sapply(strsplit(ifelse(is.na(wide$genes), "", wide$genes), "[,;]"), function(x) g %in% x)
gpt <- wide[chr == 8 & abs(pos - 145733794) < 5e5]
ctl <- rbind(
  wide[locus == "ctrl_PNPLA3_rs738409" | (chr == 19 & pos == 19379549) | has("GCKR"),
       .(control = "liver-fat effect expected (lfsr < 0.05)", locus, genes, z_pdff, postmean_pdff, lfsr_pdff, z_cirr_ghouse, lfsr_cirr_ghouse, pass = lfsr_pdff < 0.05)],
  gpt[, .(control = "GPT locus: liver-fat |z| < 2", locus, genes, z_pdff, postmean_pdff, lfsr_pdff, z_cirr_ghouse, lfsr_cirr_ghouse, pass = abs(z_pdff) < 2)],
  wide[has("GOT1") | locus == "ctrl_GOT1_lead",
       .(control = "GOT1 locus: liver-fat |z| < 2", locus, genes, z_pdff, postmean_pdff, lfsr_pdff, z_cirr_ghouse, lfsr_cirr_ghouse, pass = abs(z_pdff) < 2)],
  wide[grepl("^ctrl_HSD17B13", locus),
       .(control = "HSD17B13: ALT-raising allele raises cirrhosis, liver-fat |z| < 2", locus, genes, z_pdff, postmean_pdff, lfsr_pdff, z_cirr_ghouse, lfsr_cirr_ghouse, pass = abs(z_pdff) < 2 & z_cirr_ghouse > 0)])
fwrite(ctl, file.path(out, "controls.tsv"), sep = "\t")
print(ctl)
# Verify zero NaN among evaluated controls
eval_ctl <- ctl[is.finite(z_pdff)]
stopifnot(all(is.finite(eval_ctl$postmean_pdff)), all(is.finite(eval_ctl$lfsr_pdff)))
eval_ctl_c <- ctl[is.finite(z_cirr_ghouse)]
stopifnot(all(is.finite(eval_ctl_c$lfsr_cirr_ghouse)))

if (nrow(gpt)) cat("GPT locus z(MVP 571.5) / z(MVP ALT):", round(gpt$z_nafld_mvp / gpt$z_alt_mvp, 3), "\n")
g68 <- loci[chr == 8 & abs(pos - 145733794) < 5e5]
if (nrow(g68)) cat("GPT locus, DIG1 table z(MVP 571.5) / z(MVP ALT):", round(g68$z_MVP_NAFLD_EUR / g68$z_MVP_ALT_EUR, 3), "\n")

cat("\n== Positive liver-fat loci (lfsr_pdff < 0.05) ==\n")
print(wide[is_locus & lfsr_pdff < 0.05, .(locus, chr, pos, genes, defining_enzyme, z_pdff = round(z_pdff, 2), postmean_pdff = signif(postmean_pdff, 3), lfsr_pdff = signif(lfsr_pdff, 3), z_cirr_ghouse = round(z_cirr_ghouse, 2), lfsr_cirr_ghouse = signif(lfsr_cirr_ghouse, 3))])

summ <- data.table(
  quantity = c("loci with liver-fat and enzyme rows", "liver fat lfsr < 0.05 (F3)",
               "liver fat lfsr < 0.05 & postmean > 0 (concordant)",
               "liver fat lfsr < 0.05 & postmean < 0 (reversed)",
               "cirrhosis lfsr < 0.05 (secondary)",
               "FinnGen NAFLD lfsr < 0.05 (out of sample)",
               "both liver fat and cirrhosis lfsr < 0.05",
               "median fat effect per unit ALT effect"),
  value = c(sum(is_locus & is.finite(wide$lfsr_pdff)),
            sum(is_locus & wide$lfsr_pdff < 0.05, na.rm = TRUE),
            sum(is_locus & wide$lfsr_pdff < 0.05 & wide$postmean_pdff > 0, na.rm = TRUE),
            sum(is_locus & wide$lfsr_pdff < 0.05 & wide$postmean_pdff < 0, na.rm = TRUE),
            sum(is_locus & wide$lfsr_cirr_ghouse < 0.05, na.rm = TRUE),
            sum(is_locus & wide$lfsr_nafld_finngen < 0.05, na.rm = TRUE),
            sum(is_locus & wide$lfsr_pdff < 0.05 & wide$lfsr_cirr_ghouse < 0.05, na.rm = TRUE),
            median(wide$fat_per_alt[is_locus], na.rm = TRUE)))
fwrite(summ, file.path(out, "summary.tsv"), sep = "\t"); print(summ)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
cat("S2 20261118 DONE\n")
