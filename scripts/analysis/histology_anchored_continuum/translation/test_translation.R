#!/usr/bin/env Rscript
source("scripts/analysis/histology_anchored_continuum/translation/lib_translation.R")
h <- tr_harmonize(c("A", "A", "A", "A", "A", "A"),
  c("C", "C", "T", "C", "C", "C"),
  c("A", "C", "A", "T", "A", "A"),
  c("C", "A", "T", "G", "C", "C"),
  c(1, 1, 1, 1, NA, -1), rep(2, 6))
stopifnot(identical(h$marginal_expression_direction, c(1, -1, NA, NA, NA, -1)),
  tr_direction(c(.5, .5), c(1, -1))$unresolved_mass == 0,
  is.na(tr_direction(c(.5, .5), c(1, -1))$marginal_direction),
  is.na(tr_direction(c(.94, .06), c(1, NA))$marginal_direction),
  tr_direction(c(.96, .04), c(-1, NA))$marginal_direction == -1,
  identical(tr_bh(c(.001, NA), 117), c(.117, NA)),
  identical(tr_trait(c("PDFF", "NAFLD", "ALT", "BMI"), c(1, 1, 2, 3)),
            c("liver_fat", "direct_MASLD", "liver_enzyme", "other_trait")))
# Regression check against direct residualization; composition is correlated
# with both outcome and axis and must recover the conditional coefficient.
set.seed(20260906)
n <- 160L
d <- data.table(axis_z = rnorm(n), fibrosis_stage = rep(0:3, 40),
                inferred_sex = rep(c("F", "M"), each = 80))
d[, ilr_01 := .8 * axis_z + rnorm(n)]
d[, outcome_z := .5 * axis_z + 2 * ilr_01 + rnorm(n, sd = .3)]
a <- tr_fit(d, TRUE)
f <- lm(outcome_z ~ axis_z + factor(fibrosis_stage) + factor(inferred_sex) + ilr_01, d)
stopifnot(a$estimable, abs(a$beta - coef(f)[["axis_z"]]) < 1e-10,
          a$hc3_se > 0, tr_fit(d, FALSE)$beta > a$beta + 1)
cat("PASS: allele swaps, ambiguity, missing-mass direction, BH family, traits, composition regression\n")
