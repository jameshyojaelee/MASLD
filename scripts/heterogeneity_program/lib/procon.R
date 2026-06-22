# ─────────────────────────────────────────────────────────────────────────────
# T3 — PROCON statistical core: program-level multimodal convergence.
#
#  · acat()            — Cauchy combination test (Liu & Xie 2019); dependence-
#                        robust omnibus p, no covariance to estimate.
#  · program_camera()  — limma::cameraPR competitive enrichment (inter-gene-
#                        correlation corrected, rank-based) → signed z + p per
#                        (program, modality).
#  · procon_combine()  — between-modality combine: metafor REML (pooled signed
#                        effect + I² + p; PRIMARY/interpretable) AND ACAT
#                        (omnibus p; CO-PRIMARY). Output is a heuristic
#                        convergence score, NOT a posterior (only the genetic arm
#                        is a true Bayes factor). Calibration = permutation (perm_null).
# ─────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({ library(data.table); library(limma) })

# Cauchy combination test (ACAT). p: vector of p-values; w: optional weights.
acat <- function(p, w = NULL) {
  keep <- is.finite(p); pp <- p[keep]                       # M1 fix: filter p and w together
  if (length(pp) == 0L) return(NA_real_)
  if (length(pp) == 1L) return(pp)
  w <- if (is.null(w)) rep(1 / length(pp), length(pp)) else { ww <- w[keep]; ww / sum(ww) }
  eps <- 1e-15; pp <- pmin(pmax(pp, eps), 1 - eps)
  stat <- sum(w * tan((0.5 - pp) * pi))         # weighted sum of Cauchy ~ Cauchy
  pcauchy(stat, lower.tail = FALSE)             # numerically stable tail
}

# Competitive enrichment of each program in a per-gene ranking statistic.
# stat: named numeric (gene → signed modality statistic). sets: named list of gene vectors.
# Returns data.table(program_id, n, direction, p, z) — z = signed standardized effect.
program_camera <- function(stat, sets, min_set = 5L) {
  stat <- stat[is.finite(stat)]
  idx <- lapply(sets, function(g) which(names(stat) %in% g))
  idx <- idx[vapply(idx, length, 1L) >= min_set]
  if (length(idx) == 0L) return(NULL)
  cr <- limma::cameraPR(stat, idx, use.ranks = TRUE)     # rank-based (non-parametric), VIF-corrected
  data.table(program_id = rownames(cr), n = cr$NGenes, direction = as.character(cr$Direction),
             p = cr$PValue,
             z = ifelse(cr$Direction == "Up", 1, -1) *
                 qnorm(pmin(pmax(1 - cr$PValue / 2, 1e-15), 1 - 1e-15)))
}

# Between-modality combine per program. dt: long (program_id, modality, z, p[, sei]).
# `sei` (optional) lets thin/low-coverage modalities be down-weighted (review H1):
# metafor weights by 1/sei², so a 848-gene spatial arm can carry less than bulk.
procon_combine <- function(dt) {
  has_metafor <- requireNamespace("metafor", quietly = TRUE)
  if (!"sei" %in% names(dt)) dt[, sei := 1]
  dt[, {
    ok <- is.finite(z) & is.finite(sei); zz <- z[ok]; se <- sei[ok]; pp <- p[is.finite(p)]
    k <- length(zz)
    sign_frac <- if (k >= 1L) max(mean(zz > 0), mean(zz < 0)) else NA_real_   # dominant-direction agreement
    if (k >= 2L && has_metafor) {
      fit <- tryCatch(metafor::rma(yi = zz, sei = se, method = "REML"),
                error = function(e) tryCatch(metafor::rma(yi = zz, sei = se, method = "DL"),
                error = function(e) NULL))
      eff <- if (!is.null(fit)) as.numeric(fit$beta) else weighted.mean(zz, 1/se^2)
      i2  <- if (!is.null(fit)) fit$I2 else NA_real_
      pm  <- if (!is.null(fit)) fit$pval else 2 * pnorm(-abs(eff) * sqrt(sum(1/se^2)))
    } else if (k >= 1L) {
      eff <- weighted.mean(zz, 1/se^2); i2 <- NA_real_
      pm  <- 2 * pnorm(-abs(eff))                 # H3 fix: unit-variance single z, NOT *sqrt(k)
    } else { eff <- NA_real_; i2 <- NA_real_; pm <- NA_real_ }
    .(n_modalities = k, procon_effect = eff, procon_I2 = i2,
      procon_p_meta = pm, procon_p_acat = acat(pp), sign_frac = sign_frac)
  }, by = program_id]
}

# Concordance label. M2 fix: I² is unreliable for n_mod<4 → don't let it force
# "Conflicted". M5 fix: a genuinely sign-SPLIT program (modalities disagree in
# direction) gets its own "Split" state instead of collapsing to "ns"/"Conflicted".
procon_label <- function(effect, I2, q, n_mod, sign_frac, sig_thresh = 0.05) {
  fcase(
    is.na(q) | q >= sig_thresh,                          "ns",
    !is.na(sign_frac) & sign_frac < 0.6,                 "Split",
    effect > 0 & (n_mod < 4 | is.na(I2) | I2 < 50),      "Concordant_up",
    effect < 0 & (n_mod < 4 | is.na(I2) | I2 < 50),      "Concordant_down",
    default = "Conflicted")
}
