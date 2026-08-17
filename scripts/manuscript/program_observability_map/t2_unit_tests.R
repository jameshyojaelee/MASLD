# Unit tests for the T2 liver-proteome system, on synthetic data only.
#
# Each test corresponds to a way this system could be silently wrong: the
# membership schema trap returning zero genes, the matrix HC3 disagreeing with
# the reference per-feature fit, a batch label read off the wrong part of a
# filename, an MDE that does not match the interval it is supposed to describe,
# and a calibration harness that passes a view it should refuse.

suppressPackageStartupMessages({
  library(data.table)
})

script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file <- sub("^--file=", "", args[grepl("^--file=", args)])
  if (length(file)) return(dirname(normalizePath(file[[1]])))
  getwd()
}
here <- script_dir()
source(file.path(here, "analysis_lib.R"))
source(file.path(here, "calibration_lib.R"))
source(file.path(here, "systems_contract.R"))
source(file.path(here, "t2_lib.R"))

passed <- 0L
check <- function(label, expression) {
  ok <- tryCatch(isTRUE(expression), error = function(e) {
    cat("  ERROR in ", label, ": ", conditionMessage(e), "\n", sep = "")
    FALSE
  })
  cat(sprintf("  [%s] %s\n", if (ok) "ok" else "FAIL", label))
  if (!ok) fail("Unit test failed: ", label)
  passed <<- passed + 1L
}
expect_error <- function(label, expr) {
  ok <- inherits(tryCatch(force(expr), error = function(e) e), "error")
  cat(sprintf("  [%s] %s\n", if (ok) "ok" else "FAIL", label))
  if (!ok) fail("Unit test failed (expected an error): ", label)
  passed <<- passed + 1L
}

set.seed(20260812L)

cat("\n== batch derivation ==\n")
check("2019 and 2020 prefixes map to distinct levels",
      identical(as.character(t2_batch_from_filename(
        c("20191017_FP_Qexn2_1015_22cmr2um_Human_female_10.raw",
          "11122020_FP_ATLASLiver_vis1_ID1.raw"))), c("b2019", "b2020")))
check("b2020 is the reference level",
      levels(t2_batch_from_filename("11122020_FP_x.raw"))[[1]] == "b2020")
expect_error("an unrecognised acquisition prefix is refused",
             t2_batch_from_filename("20220101_FP_x.raw"))

cat("\n== membership schema adapter ==\n")
v2 <- data.table(program_uid = rep(c("P1", "P2"), each = 3L),
                 mapped_symbol = c("AAA", "BBB", "CCC", "DDD", "EEE", "FFF"),
                 original_l1_weight = rep(c(0.5, 0.3, 0.2), 2L))
legacy <- data.table(program_uid = rep(c("P1", "P2"), each = 3L),
                     mapped_symbol = c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE),
                     gene_symbol = c("AAA", "BBB", "CCC", "DDD", "EEE", "FFF"),
                     original_l1_weight = rep(c(0.5, 0.3, 0.2), 2L))
check("v2 schema yields the symbol strings, not a boolean filter",
      identical(sort(t2_membership_adapter(v2)$mapped_symbol),
                sort(c("AAA", "BBB", "CCC", "DDD", "EEE", "FFF"))))
check("legacy schema is filtered on the boolean and reads gene_symbol",
      identical(sort(t2_membership_adapter(legacy)$mapped_symbol),
                sort(c("AAA", "BBB", "DDD", "EEE", "FFF"))))
check("adapter lowercases nothing and uppercases input",
      t2_membership_adapter(data.table(program_uid = "P", mapped_symbol = " aaa ",
                                       original_l1_weight = 1))$mapped_symbol == "AAA")
expect_error("a membership with no usable symbols is refused",
             t2_membership_adapter(data.table(program_uid = "P", mapped_symbol = NA_character_,
                                              original_l1_weight = 1)))

cat("\n== coverage ==\n")
cov_all <- t2_program_coverage(t2_membership_adapter(v2), c("AAA", "BBB", "CCC",
                                                            "DDD", "EEE", "FFF"), "all")
check("full coverage is 1", all(abs(cov_all$retained_l1 - 1) < 1e-12))
cov_half <- t2_program_coverage(t2_membership_adapter(v2), c("AAA", "DDD"), "half")
check("retained L1 is weight-based, not gene-count-based",
      all(abs(cov_half$retained_l1 - 0.5) < 1e-12))
rule <- t2_apply_coverage_rule(cov_half, min_retained_l1 = 0.50, min_observed_genes = 2L)
check("the gene-count arm of the rule can veto a program that passes on weight",
      all(rule$testable == FALSE))

cat("\n== matrix HC3 agrees with the reference per-feature fit ==\n")
n <- 58L
meta <- data.table(
  sample_id = paste0("S", seq_len(n)),
  batch = factor(rep(c("b2020", "b2019"), c(25L, 33L)), levels = c("b2020", "b2019")),
  sex_final = factor(sample(c("Female", "Male"), n, TRUE)),
  age_z = rnorm(n), bmi_z = rnorm(n),
  steatosis = sample(0:3, n, TRUE), ballooning = sample(0:2, n, TRUE),
  inflammation = sample(0:3, n, TRUE))
meta[, nas := steatosis + ballooning + inflammation]
n_feature <- 20L
Y <- matrix(rnorm(n_feature * n), n_feature, n,
            dimnames = list(paste0("P", seq_len(n_feature)), meta$sample_id))
Y["P1", ] <- 0.6 * meta$ballooning + rnorm(n, sd = 0.7)
view <- t2_fit_view(Y, meta, T2_HISTOLOGY_COMPONENTS, "unit_joint")
ref_data <- as.data.frame(meta)
ref_data$score <- Y["P1", ]
ref <- hc3_fit(ref_data$score, ref_data,
               score ~ batch + sex_final + age_z + bmi_z + steatosis + ballooning + inflammation)
got <- view$rows[program_uid == "P1"]
check("matrix HC3 coefficients match hc3_fit()",
      max(abs(got[order(term), beta] -
                ref$coefficients[sort(T2_HISTOLOGY_COMPONENTS)])) < 1e-9)
check("matrix HC3 standard errors match hc3_fit()",
      max(abs(got[order(term), se_hc3] - ref$se[sort(T2_HISTOLOGY_COMPONENTS)])) < 1e-9)
check("matrix HC3 p-values match hc3_fit()",
      max(abs(got[order(term), p_value] - ref$p_value[sort(T2_HISTOLOGY_COMPONENTS)])) < 1e-9)
check("residual df is n minus the number of model columns",
      unique(view$rows$residual_df) == n - unique(view$rows$n_model_columns))
check("the planted ballooning effect is recovered",
      got[term == "ballooning", p_value] < 0.01)

cat("\n== MDE and support classification ==\n")
cls <- t2_classify(copy(view$rows)[, q_value := p.adjust(p_value, "BH")], 0.20)
check("every negative carries a finite MDE",
      all(is.finite(cls[support != "supported", mde])))
check("MDE is the two-sided 80%-power effect at the model's own df",
      max(abs(cls$mde - (stats::qt(0.975, df = cls$residual_df) + stats::qnorm(0.80)) *
                cls$se_hc3)) < 1e-12)
check("an informative null has a per-axis-SD interval inside the threshold",
      all(pmax(abs(cls[support == "informative_null", ci_lower_per_axis_sd]),
               abs(cls[support == "informative_null", ci_upper_per_axis_sd])) < 0.20))
check("an indeterminate call still admits a threshold-sized effect",
      all(pmax(abs(cls[support == "indeterminate", ci_lower_per_axis_sd]),
               abs(cls[support == "indeterminate", ci_upper_per_axis_sd])) >= 0.20))

cat("\n== the threshold is applied on a scale that is comparable across axes ==\n")
check("an ordinal term is rescaled by its own SD",
      all(abs(cls[term == "steatosis", threshold_scale_factor] -
                stats::sd(meta$steatosis)) < 1e-12))
meta_bin <- copy(meta)[, ballooning_binary := as.integer(ballooning >= 1)]
bin_view <- t2_fit_view(Y, meta_bin, "ballooning_binary", "unit_binary")
bin <- t2_classify(copy(bin_view$rows)[, q_value := p.adjust(p_value, "BH")], 0.20)
check("a binary term is recognised as binary and left in native units",
      all(bin$term_scale == "binary") && all(bin$threshold_scale_factor == 1))
check("rescaling an axis changes the native beta but not the verdict or the MDE",
      {
        meta_k <- copy(meta)[, steatosis := steatosis * 4]
        k <- t2_fit_view(Y, meta_k, T2_HISTOLOGY_COMPONENTS, "unit_scaled")
        kc <- t2_classify(copy(k$rows)[, q_value := p.adjust(p_value, "BH")], 0.20)
        a <- cls[term == "steatosis"][order(program_uid)]
        b <- kc[term == "steatosis"][order(program_uid)]
        max(abs(b$beta * 4 - a$beta)) < 1e-9 &&
          max(abs(b$mde_per_axis_sd - a$mde_per_axis_sd)) < 1e-9 &&
          identical(a$support, b$support)
      })

cat("\n== contrast Wald ==\n")
ref_fit <- hc3_fit(ref_data$score, ref_data,
                   score ~ batch + sex_final + age_z + bmi_z +
                     steatosis + ballooning + inflammation)
check("the identity contrast reduces exactly to wald_test()",
      {
        a <- contrast_wald(ref_fit, CONTRAST_OMNIBUS)
        b <- wald_test(ref_fit, T2_HISTOLOGY_COMPONENTS)
        abs(a$f_statistic - b$f_statistic) < 1e-9 && abs(a$p_value - b$p_value) < 1e-12
      })
# Equal coefficients on (s, b, l) is algebraically the same claim as zero
# coefficients on (s - l) and (b - l) once the composite s + b + l is in the
# model. Fitting that reparameterisation and testing it with the library's own
# wald_test() checks the contrast against something that was not written here.
check("the equal-weight contrast equals wald_test() on the reparameterised model",
      {
        rp <- ref_data
        rp$nas_sum <- rp$steatosis + rp$ballooning + rp$inflammation
        rp$d_stea <- rp$steatosis - rp$inflammation
        rp$d_ball <- rp$ballooning - rp$inflammation
        rp_fit <- hc3_fit(rp$score, rp,
                          score ~ batch + sex_final + age_z + bmi_z +
                            nas_sum + d_stea + d_ball)
        a <- contrast_wald(ref_fit, CONTRAST_EQUAL_WEIGHT)
        b <- wald_test(rp_fit, c("d_stea", "d_ball"))
        abs(a$f_statistic - b$f_statistic) < 1e-7 && abs(a$p_value - b$p_value) < 1e-9
      })
check("the equal-weight contrast has 2 df and the omnibus 3",
      contrast_wald(ref_fit, CONTRAST_EQUAL_WEIGHT)$df1 == 2L &&
        contrast_wald(ref_fit, CONTRAST_OMNIBUS)$df1 == 3L)
check("a program with three equal coefficients does not reject equal weighting",
      {
        eq <- ref_data
        eq$score <- 0.5 * (eq$steatosis + eq$ballooning + eq$inflammation) +
          rnorm(nrow(eq), sd = 0.5)
        f <- hc3_fit(eq$score, eq, score ~ batch + sex_final + age_z + bmi_z +
                       steatosis + ballooning + inflammation)
        contrast_wald(f, CONTRAST_EQUAL_WEIGHT)$p_value > 0.05 &&
          contrast_wald(f, CONTRAST_OMNIBUS)$p_value < 0.01
      })
check("a program driven by one component alone does reject equal weighting",
      {
        one <- ref_data
        one$score <- 1.2 * one$steatosis + rnorm(nrow(one), sd = 0.5)
        f <- hc3_fit(one$score, one, score ~ batch + sex_final + age_z + bmi_z +
                       steatosis + ballooning + inflammation)
        contrast_wald(f, CONTRAST_EQUAL_WEIGHT)$p_value < 0.01
      })
check("a contrast naming an absent term is refused, not silently dropped",
      {
        bad <- CONTRAST_EQUAL_WEIGHT
        colnames(bad) <- c("steatosis", "ballooning", "not_a_term")
        isFALSE(contrast_wald(ref_fit, bad)$estimable)
      })

cat("\n== score dimensionality ==\n")
dims <- t2_score_dimensionality(Y)
check("effective dimension of independent features is near the feature count",
      dims$effective_dimension_participation_ratio > 0.7 * nrow(Y))
Yc <- rbind(Y, Y[1, , drop = FALSE] + 1e-9)
rownames(Yc)[nrow(Yc)] <- "P_dup"
check("duplicating a feature does not raise the effective dimension",
      t2_score_dimensionality(Yc)$effective_dimension_participation_ratio <
        dims$effective_dimension_participation_ratio + 1)
check("a near-duplicate is caught by the maximum pairwise correlation",
      t2_score_dimensionality(Yc)$max_abs_pairwise_correlation > 0.999)
assert_negatives_have_mde(cls, "support", "mde")
check("contract accepts these negatives", TRUE)

cat("\n== the shared calibration harness on this view ==\n")
perm_p <- function(i) {
  pm <- permute_histology_within_cohort(copy(meta), T2_HISTOLOGY_COMPONENTS,
                                        cohort_column = "batch")
  t2_view_p(Y, pm, T2_HISTOLOGY_COMPONENTS)
}
cal <- calibrate_view(view$rows$p_value, perm_p, n_reps = 200L,
                      n_tests = nrow(view$rows), label = "unit_joint")
cat(sprintf("    null false calls per %d tests: mean %.2f, max %d\n",
            cal$n_tests, cal$null_call_mean, as.integer(cal$null_call_max)))
check("a parametric view carries no permutation resolution floor",
      cal$resolution_floor == 0 && isTRUE(cal$resolvable))
check("the joint HC3 view is calibrated on synthetic data", isTRUE(cal$calibrated))
check("gate_view passes a calibrated view", isTRUE(gate_view(cal)$passed))

cat("\n== within-batch permutation preserves what it must ==\n")
pm <- permute_histology_within_cohort(copy(meta), T2_PERMUTED_COLUMNS[1:4],
                                      cohort_column = "batch")
check("permutation preserves the batch-histology association exactly",
      identical(pm[, sort(ballooning), by = batch], meta[, sort(ballooning), by = batch]))
check("permutation preserves the mutual correlation of the axes",
      abs(cor(pm$steatosis, pm$ballooning) - cor(meta$steatosis, meta$ballooning)) < 1e-12)
check("permutation does move donors",
      !identical(pm$ballooning, meta$ballooning))

cat("\n== the gate actually refuses ==\n")
bad <- list(label = "deliberately_uncalibrated", calibrated = FALSE, resolvable = TRUE,
            null_call_mean = 9.3, n_tests = 60L, tolerance = 1.0,
            resolution_floor = 0, bh_threshold_needed = 0.05 / 60)
expect_error("gate_view stops on an uncalibrated view", gate_view(bad))
expect_error("assert_counts_calibrated stops when any view is unreportable",
             assert_counts_calibrated(1L, data.table(view = "x", reportable = FALSE)))

cat("\n== language ==\n")
assert_language("cross-sectional association of program score with recorded histologic grade")
check("permitted wording passes", TRUE)
expect_error("forbidden ordering vocabulary is caught",
             assert_language("the program follows a severity axis"))

cat(sprintf("\nAll %d T2 unit tests passed.\n", passed))
