# Statistical and evidence rules for the continuum translation extension.
# Marginal directions do not identify conditional multi-signal effects.
suppressPackageStartupMessages(library(data.table))

tr_bh <- function(p, n = 117L) {
  stopifnot(length(p) <= n)
  ans <- rep(NA_real_, length(p))
  ok <- is.finite(p)
  ans[ok] <- p.adjust(p[ok], method = "BH", n = n)
  ans
}

tr_harmonize <- function(a1, a2, ea, nea, gb, eb) {
  a1 <- toupper(a1); a2 <- toupper(a2)
  ea <- toupper(ea); nea <- toupper(nea)
  valid <- !is.na(a1) & !is.na(a2) & !is.na(ea) & !is.na(nea)
  snv <- valid & a1 %in% c("A", "C", "G", "T") &
    a2 %in% c("A", "C", "G", "T") & a1 != a2
  direct <- valid & a1 == ea & a2 == nea
  reverse <- valid & a1 == nea & a2 == ea
  palindrome <- snv & paste0(a1, a2) %in% c("AT", "TA", "CG", "GC")
  state <- rep("allele_pair_unresolved", length(a1))
  state[snv & direct] <- "direct"
  state[snv & reverse] <- "swapped"
  state[palindrome] <- "palindromic_unresolved"
  state[!is.finite(gb) | !is.finite(eb)] <- "effect_missing_or_duplicate"
  aligned <- rep(NA_real_, length(a1))
  aligned[state == "direct"] <- eb[state == "direct"]
  aligned[state == "swapped"] <- -eb[state == "swapped"]
  direction <- sign(gb * aligned)
  direction[direction == 0] <- NA_real_
  data.table(allele_state = state, aligned_eqtl_beta = aligned,
             marginal_expression_direction = direction)
}

tr_direction <- function(pp, direction, cutoff = 0.95) {
  stopifnot(all(is.finite(pp)), all(pp >= 0), sum(pp) <= 1.00001)
  positive <- sum(pp[!is.na(direction) & direction > 0])
  negative <- sum(pp[!is.na(direction) & direction < 0])
  # Do not renormalize after removing unresolvable alleles or absent effects.
  data.table(positive_mass = positive, negative_mass = negative,
             unresolved_mass = max(0, 1 - positive - negative),
             marginal_direction = if (positive >= cutoff) 1L else
               if (negative >= cutoff) -1L else NA_integer_)
}

tr_trait <- function(trait, tier) {
  ifelse(grepl("PDFF|liver.?fat", trait, ignore.case = TRUE), "liver_fat",
         ifelse(tier == 1, "direct_MASLD", ifelse(tier == 2,
           "liver_enzyme", "other_trait")))
}

tr_fit <- function(d, composition = FALSE) {
  cols <- if (composition) grep("^ilr_", names(d), value = TRUE) else character()
  form <- reformulate(c("axis_z", "factor(fibrosis_stage)",
                        "factor(inferred_sex)", cols), "outcome_z")
  blank <- data.table(estimable = FALSE, n_participants = nrow(d),
    beta = NA_real_, se = NA_real_, p_value = NA_real_, hc3_se = NA_real_,
    hc3_p_value = NA_real_, residual_df = NA_integer_, max_leverage = NA_real_,
    failure_reason = "insufficient_or_singular_design")
  if (nrow(d) < 30 || uniqueN(d$fibrosis_stage) < 2 ||
      uniqueN(d$inferred_sex) < 2) return(blank)
  fit <- tryCatch(lm(form, data = d), error = function(e) NULL)
  if (is.null(fit)) return(blank)
  x <- model.matrix(fit)
  if (fit$rank < ncol(x) || df.residual(fit) < 10) return(blank)
  cf <- coef(summary(fit))["axis_z", ]
  h <- hatvalues(fit)
  inv <- chol2inv(qr.R(qr(x)))
  # qr pivoting is immaterial for a full-rank, unpivoted lm design checked here.
  if (!identical(qr(x)$pivot, seq_len(ncol(x)))) return(blank)
  meat <- crossprod(x * as.numeric(residuals(fit) / (1 - h)))
  vc <- inv %*% meat %*% inv
  k <- match("axis_z", colnames(x)); hc3 <- sqrt(vc[k, k])
  data.table(estimable = TRUE, n_participants = nrow(d), beta = cf[[1]],
    se = cf[[2]], p_value = cf[[4]], hc3_se = hc3,
    hc3_p_value = 2 * pt(-abs(cf[[1]] / hc3), df.residual(fit)),
    residual_df = df.residual(fit), max_leverage = max(h), failure_reason = NA_character_)
}

tr_meta <- function(d) {
  ok <- d[estimable == TRUE & is.finite(beta) & is.finite(se) & se > 0]
  if (nrow(ok) != 2L || uniqueN(ok$dataset) != 2L) return(data.table(
    beta = NA_real_, se = NA_real_, p_value = NA_real_,
    direction_consistent = FALSE, n_participants = sum(ok$n_participants),
    both_cohort_p_below_05 = FALSE))
  w <- 1 / ok$se^2; b <- sum(w * ok$beta) / sum(w); s <- sqrt(1 / sum(w))
  data.table(beta = b, se = s, p_value = 2 * pnorm(-abs(b / s)),
    direction_consistent = all(sign(ok$beta) == sign(b)),
    n_participants = sum(ok$n_participants),
    both_cohort_p_below_05 = all(ok$p_value < 0.05))
}
