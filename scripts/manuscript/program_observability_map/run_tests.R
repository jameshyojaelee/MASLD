#!/usr/bin/env Rscript

# Unit tests for the v9 additions. Deliberately small and synthetic so this runs
# as light debugging before any compute job is submitted.

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_dir <- dirname(normalizePath(script_arg))
suppressPackageStartupMessages(library(data.table))
source(file.path(script_dir, "analysis_lib.R"))

set.seed(1)
passed <- 0L
check <- function(label, condition, detail = "") {
  if (!isTRUE(condition)) fail("FAIL: ", label, if (nzchar(detail)) paste0(" [", detail, "]") else "")
  passed <<- passed + 1L
  cat("  ok  ", label, "\n", sep = "")
}

# ---- wald_test reduces to the squared t on a single term ------------------
n <- 120L
d <- data.table(x = rnorm(n), z = rnorm(n), sex_final = sample(c("M", "F"), n, TRUE))
d[, score := 0.4 * x + 0.2 * z + rnorm(n)]
fit <- hc3_fit(d$score, d, score ~ sex_final + x + z)
w <- wald_test(fit, "x")
check("wald_test on one term equals the squared HC3 t",
      abs(w$chisq - fit$statistic[["x"]]^2) < 1e-9,
      sprintf("%.10f vs %.10f", w$chisq, fit$statistic[["x"]]^2))
check("wald_test one-term p equals the two-sided t p",
      abs(w$p_value - fit$p_value[["x"]]) < 1e-9)
w2 <- wald_test(fit, c("x", "z"))
check("wald_test on two terms carries 2 numerator df", w2$df1 == 2L)

# ---- shape_basis is centered and near-orthogonal --------------------------
b <- shape_basis(c(0, 1, 1, 2, 2, 2, 3, 3, 4))
check("shape_basis linear column is mean-zero", abs(mean(b$linear)) < 1e-12)
check("shape_basis quadratic column is mean-zero", abs(mean(b$quadratic)) < 1e-12)
check("shape_basis scaling is independent of sample size",
      identical(shape_basis(c(0, 2, 4))$linear, c(-2, 0, 2)))

# ---- a quadratic axis is invisible to the linear term but not the 2-df test
m <- 6L
cohorts <- paste0("C", 1:4)
meta <- rbindlist(lapply(cohorts, function(cc) {
  nn <- 90L
  data.table(sample_id = paste0(cc, "_", seq_len(nn)), dataset = cc,
             sex_final = sample(c("M", "F"), nn, TRUE),
             fibrosis_stage = sample(0:4, nn, TRUE),
             nas_score = sample(0:8, nn, TRUE))
}))
scores <- matrix(rnorm(m * nrow(meta)), nrow = m,
                 dimnames = list(paste0("P", seq_len(m)), meta$sample_id))
# P1 responds to NAS as a symmetric U, which has zero linear component
u <- (meta$nas_score - mean(meta$nas_score))^2
scores["P1", ] <- as.numeric(scale(u)) * 1.2 + rnorm(nrow(meta), sd = 0.4)
# P2 responds linearly to fibrosis
scores["P2", ] <- 0.5 * meta$fibrosis_stage + rnorm(nrow(meta), sd = 0.6)

lin <- fit_feature_models(scores, meta, cohorts, "primary")
shp <- fit_axis_shape_models(scores, meta, cohorts)
# The shape model's linear term is conditional on curvature, so it is close to
# but not identical with the v8 linear estimand. This asserts the closeness and
# guards against the two arms diverging, without pretending they are the same
# quantity. The v8 estimand itself comes from the standalone linear arm.
cmp <- merge(lin[axis == "fibrosis", .(feature_id, cohort, beta_linear_arm = beta)],
             shp[axis == "fibrosis", .(feature_id, cohort, beta_shape_arm = beta_linear)],
             by = c("feature_id", "cohort"))
check("linear and shape arms track each other closely on fibrosis",
      cor(cmp$beta_linear_arm, cmp$beta_shape_arm) > 0.95,
      sprintf("pearson %.4f, max abs diff %.4f",
              cor(cmp$beta_linear_arm, cmp$beta_shape_arm),
              max(abs(cmp$beta_linear_arm - cmp$beta_shape_arm))))
check("they are not identical, because the shape linear term is conditional on curvature",
      max(abs(cmp$beta_linear_arm - cmp$beta_shape_arm)) > 0)
lin_nas_p1 <- lin[feature_id == "P1" & axis == "nas", median(abs(statistic_hc3))]
shape_nas_p1 <- shp[feature_id == "P1" & axis == "nas", median(joint_f)]
check("a U-shaped NAS response is weak under the linear term", lin_nas_p1 < 3,
      sprintf("median |t| = %.2f", lin_nas_p1))
check("the same response is strong under the 2-df shape test", shape_nas_p1 > 10,
      sprintf("median F = %.2f", shape_nas_p1))

sm <- meta_analyze_shape(shp)
check("meta_analyze_shape returns one row per program and axis", nrow(sm) == m * 2L)
check("meta_analyze_shape recovers the U-shaped NAS program",
      sm[feature_id == "P1" & axis == "nas", is.finite(p_joint_meta) & p_joint_meta < 0.01])

# Regression test for the reference bug that inflated the first v9 discovery.
# p_joint_meta must come from F(df1, k-1), never chi-square(df1). The chi-square
# reference is retained only as a diagnostic column, and it must be uniformly the
# more liberal of the two, which is exactly why it produced 13.0 false positives
# per 113 programs under a permutation null against 0.0 for F.
sm_ok <- sm[estimable == TRUE]
check("the joint test declares its reference", all(sm_ok$joint_reference == "F(df1, k-1)"))
check("joint denominator df is k-1", all(sm_ok$joint_df2 == sm_ok$n_cohorts - 1L))
check("p_joint_meta matches an F reference, not chi-square",
      all(abs(sm_ok$p_joint_meta -
              stats::pf(sm_ok$joint_chisq / sm_ok$joint_df1,
                        df1 = sm_ok$joint_df1, df2 = sm_ok$joint_df2,
                        lower.tail = FALSE)) < 1e-12))
check("the chi-square diagnostic is retained and is uniformly more liberal",
      all(sm_ok$p_joint_chisq_diagnostic <= sm_ok$p_joint_meta + 1e-12) &&
        any(sm_ok$p_joint_chisq_diagnostic < sm_ok$p_joint_meta),
      sprintf("median chisq p %.3g vs F p %.3g",
              median(sm_ok$p_joint_chisq_diagnostic), median(sm_ok$p_joint_meta)))

# ---- categorical omnibus and its combination ------------------------------
om <- categorical_omnibus(scores, meta, cohorts)
check("categorical_omnibus returns one row per program, axis and cohort",
      nrow(om) == m * 2L * length(cohorts))
oc <- combine_omnibus(om)
check("combine_omnibus sums degrees of freedom across cohorts",
      oc[feature_id == "P1" & axis == "nas",
         combined_df == om[feature_id == "P1" & axis == "nas" & estimable == TRUE, sum(df1)]])
check("combine_omnibus also recovers the U-shaped NAS program",
      oc[feature_id == "P1" & axis == "nas", p_combined < 0.01])

# ---- minimum detectable effect and informative nulls ----------------------
check("mde grows with the standard error",
      minimum_detectable_effect(0.1, 4L) > minimum_detectable_effect(0.05, 4L))
check("a tight interval around zero is an informative null",
      informative_null(-0.05, 0.05, 0.20))
check("a wide interval around zero is not", !informative_null(-0.5, 0.5, 0.20))

# ---- decomposition --------------------------------------------------------
dec <- decompose_effect(c(0.4, 0.4, 0), c(0.1, -0.2, 0.1))
check("retained fraction is the coefficient ratio", abs(dec$retained_fraction[[1]] - 0.25) < 1e-12)
check("a sign flip is reported as such", !dec$sign_preserved[[2]])
check("a zero unadjusted effect yields NA rather than a division", is.na(dec$retained_fraction[[3]]))

# ---- separability ---------------------------------------------------------
sep <- canonical_separability(scores, meta[, .(fibrosis_stage, nas_score)],
                              n_components = 4L, n_permutations = 200L, seed = 7L)
check("canonical_separability returns one correlation per histology column",
      length(sep$canonical_correlations) == 2L)
check("canonical correlations are ordered", sep$canonical_correlations[1] >= sep$canonical_correlations[2])
check("empirical p values are bounded by the permutation resolution",
      all(sep$empirical_p >= 1 / 201) && all(sep$empirical_p <= 1))

# ---- composition: the v8 rank-deficiency regression test -------------------
# Two structurally absent lineages must NOT silently drop a cohort. Under v8 the
# pseudocount mapped both onto the same log-ratio and the design lost rank.
lineages <- c("Hepatocytes", "Macrophages", "Cholangiocytes", "Fibroblasts", "T cells")
comp <- data.table(sample_id = meta$sample_id, dataset = meta$dataset)
comp[, Hepatocytes := runif(.N, 0.80, 0.95)]
comp[, Macrophages := runif(.N, 0.01, 0.06)]
comp[, Cholangiocytes := 0]
comp[, Fibroblasts := 0]
comp[, `T cells` := runif(.N, 0.005, 0.03)]
cs <- build_composition_scores(comp, meta, lineages, 1e-6,
                               eligible_lineages = "Macrophages", n_components = 3L)
check("only the eligible lineage produces a log-ratio column",
      identical(cs$logratio_columns, "logratio_macrophages_hepatocytes"))
check("CLR components are returned for every cohort",
      cs$pcs[, sum(!is.na(composition_pc1))] == nrow(meta))
cm <- merge(meta, cs$pcs, by = "sample_id", sort = FALSE)
cm[, score := as.numeric(scores["P2", match(sample_id, colnames(scores))])]
cm[, `:=`(fibrosis_centered = fibrosis_stage - mean(fibrosis_stage),
          nas_centered = nas_score - mean(nas_score))]
pcfit <- hc3_fit(cm$score, cm[dataset == "C1"],
                 score ~ sex_final + fibrosis_centered + nas_centered +
                   composition_pc1 + composition_pc2 + composition_pc3)
check("the CLR-PC design is estimable with two structurally absent lineages",
      isTRUE(pcfit$estimable), if (isTRUE(pcfit$estimable)) "" else pcfit$failure_reason)

# ---- matrix-form HC3 must equal the per-feature reference ------------------
# The fast path carries the 1,000-permutation discordance null, so it has to be
# provably the same estimator as hc3_fit(), not merely similar.
mm <- copy(meta[dataset == "C1"])
mm[, `:=`(fibrosis_centered = fibrosis_stage - mean(fibrosis_stage),
          nas_centered = nas_score - mean(nas_score), sex_final = factor(sex_final))]
Xm <- stats::model.matrix(~ sex_final + fibrosis_centered + nas_centered, data = mm)
Ym <- scores[, mm$sample_id, drop = FALSE]
fm <- hc3_fit_matrix(Ym, Xm)
check("hc3_fit_matrix is estimable on a full-rank design", isTRUE(fm$estimable))
ref_b <- ref_se <- numeric(nrow(Ym))
for (i in seq_len(nrow(Ym))) {
  mm[, score := as.numeric(Ym[i, ])]
  f1 <- hc3_fit(mm$score, mm, score ~ sex_final + fibrosis_centered + nas_centered)
  ref_b[i] <- f1$coefficients[["fibrosis_centered"]]
  ref_se[i] <- f1$se[["fibrosis_centered"]]
}
check("matrix-form coefficients equal hc3_fit to 1e-10",
      max(abs(fm$coefficients["fibrosis_centered", ] - ref_b)) < 1e-10,
      sprintf("max diff %.3g", max(abs(fm$coefficients["fibrosis_centered", ] - ref_b))))
check("matrix-form HC3 standard errors equal hc3_fit to 1e-10",
      max(abs(fm$se["fibrosis_centered", ] - ref_se)) < 1e-10,
      sprintf("max diff %.3g", max(abs(fm$se["fibrosis_centered", ] - ref_se))))
check("matrix-form studentized residuals match the reference to 1e-10", {
  mm[, score := as.numeric(Ym[1, ])]
  max(abs(fm$studentized_residuals[1, ] -
          studentized_residuals(mm$score, mm,
                                score ~ sex_final + fibrosis_centered + nas_centered))) < 1e-10
})
check("hc3_fit_matrix refuses a rank-deficient design",
      !isTRUE(hc3_fit_matrix(Ym, cbind(Xm, Xm[, 2]))$estimable))

# ---- shared calibration harness -------------------------------------------
source(file.path(script_dir, "calibration_lib.R"))
source(file.path(script_dir, "systems_contract.R"))

# A well-behaved view: p-values uniform under the null, so BH rejects almost never.
set.seed(11)
good <- calibrate_view(
  observed_p = c(runif(100, 0, 1e-6), runif(13)),
  permuted_p_fn = function(i) runif(113),
  n_reps = 300L, n_tests = 113L, label = "well_behaved")
check("a calibrated view passes the gate", good$calibrated && good$reportable,
      sprintf("null mean %.2f", good$null_call_mean))
check("gate_view does not raise on a calibrated view",
      isTRUE(gate_view(good)$passed))

# A miscalibrated view: p-values skewed to zero under the null, like the HC3
# categorical Wald on sparse NAS levels that gave 37.9 false calls per 113.
bad <- calibrate_view(
  observed_p = runif(113),
  permuted_p_fn = function(i) runif(113)^6,
  n_reps = 300L, n_tests = 113L, label = "anticonservative")
check("an anticonservative view is caught", !bad$calibrated,
      sprintf("null mean %.1f false calls per 113", bad$null_call_mean))
check("gate_view raises on it",
      inherits(try(gate_view(bad), silent = TRUE), "try-error"))

# A resolution failure: too few permutations for the family size, which is what
# made "0 discordant donors" arithmetic rather than biology.
thin <- calibrate_view(
  observed_p = runif(469), permuted_p_fn = function(i) runif(469),
  n_reps = 100L, n_tests = 469L, label = "under_resolved",
  observed_p_is_empirical = TRUE)
check("an under-resolved view is caught", !thin$resolvable,
      sprintf("floor %.3g vs threshold %.3g", thin$resolution_floor, thin$bh_threshold_needed))
check("required_permutations gives the fix", required_permutations(469L) == 9380L)
check("calibration_row carries reportability",
      isFALSE(calibration_row(thin)$reportable))
# The distinction the harness got wrong at first: a parametric p has no floor,
# so the same permutation count must NOT mark it unresolvable.
parametric <- calibrate_view(
  observed_p = runif(469), permuted_p_fn = function(i) runif(469),
  n_reps = 100L, n_tests = 469L, label = "parametric",
  observed_p_is_empirical = FALSE)
check("a parametric view is not penalised by the permutation floor",
      parametric$resolvable && !thin$resolvable)

# Joint permutation must preserve the correlation between the histologic axes.
mtest <- data.table(dataset = rep(c("A", "B"), each = 200))
mtest[, fibrosis_stage := sample(0:4, .N, TRUE)]
mtest[, nas_score := pmin(8, pmax(0, fibrosis_stage + sample(-1:3, .N, TRUE)))]
r_before <- mtest[, cor(fibrosis_stage, nas_score, method = "spearman")]
pm <- permute_histology_within_cohort(mtest, c("fibrosis_stage", "nas_score"))
check("joint permutation preserves the fibrosis-NAS correlation",
      abs(pm[, cor(fibrosis_stage, nas_score, method = "spearman")] - r_before) < 1e-9,
      sprintf("%.4f vs %.4f", pm[, cor(fibrosis_stage, nas_score, method = "spearman")], r_before))
check("joint permutation still breaks the sample linkage",
      !identical(pm$fibrosis_stage, mtest$fibrosis_stage))

# ---- contract assertions ---------------------------------------------------
vac <- data.table(sample_id = paste0("S", 1:10))
vac[, participant_id := sample_id]
check("a vacuous biological-unit check is rejected",
      inherits(try(assert_biological_unit(vac, "participant_id"), silent = TRUE), "try-error"))
check("the same check passes when external pairing evidence is supplied",
      isTRUE(assert_biological_unit(vac, "participant_id", cohorts = "NOPE",
                                    pairing_root = file.path(tempdir(), "absent"))))
check("forbidden ordering vocabulary is caught",
      inherits(try(assert_language("programs along the disease trajectory"), silent = TRUE),
               "try-error"))
check("permitted language passes",
      isTRUE(assert_language("cross-sectional stage-associated remodeling")))

# Scoped exemptions. The guard may only open for a live exemption id whose
# scope_paths contain the caller; every other combination must still refuse.
EXEMPT_ID <- "ORDERING-EXEMPTION-2026-08-14"
IN_SCOPE <- "scripts/analysis/continuous_axis_benchmark/02a_arm_kamzolas_exact.R"
OUT_SCOPE <- "scripts/manuscript/program_observability_map/t1_run_cellstate.R"
check("an exempt term is still rejected without an exemption id",
      inherits(try(assert_language("pseudotime"), silent = TRUE), "try-error"))
check("an exempt term passes with a live id and an in-scope caller",
      isTRUE(assert_language("pseudotime and sliding window",
                             exemption_id = EXEMPT_ID, caller_path = IN_SCOPE)))
check("an exemption does not travel outside its scope_paths",
      inherits(try(assert_language("pseudotime", exemption_id = EXEMPT_ID,
                                   caller_path = OUT_SCOPE), silent = TRUE), "try-error"))
check("an unknown exemption id licenses nothing",
      inherits(try(assert_language("pseudotime", exemption_id = "NO-SUCH-EXEMPTION",
                                   caller_path = IN_SCOPE), silent = TRUE), "try-error"))
check("a live in-scope exemption still cannot license 'progression'",
      inherits(try(assert_language("progression", exemption_id = EXEMPT_ID,
                                   caller_path = IN_SCOPE), silent = TRUE), "try-error"))
check("assert_no_exemption refuses an exempt term outright",
      inherits(try(assert_no_exemption("pseudotime"), silent = TRUE), "try-error"))
check("assert_no_exemption passes on Resource wording",
      isTRUE(assert_no_exemption("cross-sectional stage-associated remodeling")))

check("an unnamed inclusion criterion is rejected",
      inherits(try(assert_inclusion_criterion("novelty"), silent = TRUE), "try-error"))
check("a named inclusion criterion passes",
      isTRUE(assert_inclusion_criterion(c("observability", "evidence_interpretation"))))

cat("\nALL ", passed, " CHECKS PASSED\n", sep = "")
