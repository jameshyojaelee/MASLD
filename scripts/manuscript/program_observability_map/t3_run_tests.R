#!/usr/bin/env Rscript

# T3 unit tests, on synthetic data, run before the real substrate is touched.
#
# Two of these tests exist because the fast paths they check are the only
# reason the matched-random denominator is affordable at all. If either
# equivalence fails, the denominator is measuring a different estimator from the
# one the headline count came from, and the comparison is void.

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_dir <- dirname(normalizePath(script_arg))
source(file.path(script_dir, "analysis_lib.R"))
source(file.path(script_dir, "calibration_lib.R"))
source(file.path(script_dir, "systems_contract.R"))
source(file.path(script_dir, "t3_vocabulary_lib.R"))

passed <- 0L
check <- function(label, expr) {
  ok <- isTRUE(tryCatch(expr, error = function(e) {
    message("    error: ", conditionMessage(e)); FALSE
  }))
  cat(sprintf("[%s] %s\n", if (ok) "PASS" else "FAIL", label))
  if (!ok) fail("Test failed: ", label)
  passed <<- passed + 1L
}

set.seed(20260812L)

# --- synthetic substrate ---------------------------------------------------
n_genes <- 400L
cohorts <- c("C1", "C2", "C3", "C4")
n_per <- c(C1 = 60L, C2 = 80L, C3 = 55L, C4 = 70L)
meta <- rbindlist(lapply(cohorts, function(co) data.table(
  sample_id = paste0(co, "_S", seq_len(n_per[[co]])),
  dataset = co,
  sex_final = sample(c("F", "M"), n_per[[co]], replace = TRUE),
  fibrosis_stage = sample(0:4, n_per[[co]], replace = TRUE),
  nas_score = sample(0:8, n_per[[co]], replace = TRUE)
)))
X <- matrix(rnorm(n_genes * nrow(meta)), nrow = n_genes,
            dimnames = list(paste0("G", seq_len(n_genes)), meta$sample_id))
# inject a real fibrosis association into the first 40 genes
X[1:40, ] <- X[1:40, ] + rep(meta$fibrosis_stage * 0.35, each = 40L)
X <- X + 6  # keep on a logCPM-like scale so the decile matching has structure
X[1:100, ] <- X[1:100, ] + 3

membership <- rbindlist(lapply(1:6, function(p) data.table(
  feature_id = paste0("P", p),
  gene_symbol = paste0("G", sample.int(n_genes, 30L)),
  direction = "up"
)))
membership <- unique(membership, by = c("feature_id", "gene_symbol"))

# --- 1. hc3_fit_matrix reproduces the reference per-feature hc3_fit ---------
check("hc3_fit_matrix equals hc3_fit to 1e-10", {
  md <- copy(meta[dataset == "C2"])
  md[, `:=`(fibrosis_centered = fibrosis_stage - mean(fibrosis_stage),
            nas_centered = nas_score - mean(nas_score))]
  Xd <- model.matrix(T3_PRIMARY_FORMULA, data = md)
  Y <- X[1:20, md$sample_id, drop = FALSE]
  fm <- hc3_fit_matrix(Y, Xd)
  ref <- lapply(seq_len(nrow(Y)), function(i) {
    md[, score := Y[i, ]]
    hc3_fit(md$score, md, score ~ sex_final + fibrosis_centered + nas_centered)
  })
  b_ref <- vapply(ref, function(f) f$coefficients[["fibrosis_centered"]], numeric(1))
  s_ref <- vapply(ref, function(f) f$se[["fibrosis_centered"]], numeric(1))
  max(abs(fm$coefficients["fibrosis_centered", ] - b_ref)) < 1e-10 &&
    max(abs(fm$se["fibrosis_centered", ] - s_ref)) < 1e-10
})

# --- 2. fast scoring reproduces score_gene_sets -----------------------------
check("t3_score_sets_fast equals score_gene_sets to 1e-10", {
  ref <- score_gene_sets(X, meta, membership, min_genes = 10L,
                         coverage_threshold = 0.80, directional = FALSE)
  universe <- rownames(X)
  zlist <- t3_cohort_z(X, meta, cohorts, universe)
  sets <- lapply(split(membership$gene_symbol, membership$feature_id),
                 function(g) sort(match(g, universe)))
  fast <- t3_score_sets_fast(zlist, sets)
  common <- intersect(rownames(fast), rownames(ref$scores))
  ref_mat <- ref$scores[common, colnames(fast), drop = FALSE]
  length(common) == 6L && max(abs(fast[common, ] - ref_mat)) < 1e-10
})

# --- 3. fast fit reproduces fit_feature_models + meta_analyze ---------------
check("t3_fit_vocabulary equals fit_feature_models + meta_analyze to 1e-8", {
  universe <- rownames(X)
  zlist <- t3_cohort_z(X, meta, cohorts, universe)
  sets <- lapply(split(membership$gene_symbol, membership$feature_id),
                 function(g) sort(match(g, universe)))
  S <- t3_score_sets_fast(zlist, sets)
  designs <- t3_designs(meta, cohorts)
  fast <- t3_fit_vocabulary(S, designs, mc_cores = 1L)

  ref_cohort <- fit_feature_models(S, meta, cohorts, "primary")
  ref_meta <- meta_analyze(ref_cohort)
  ref_meta[, q_value := p.adjust(p_value_meta, method = "BH", n = .N), by = model_kind]

  m <- merge(fast, ref_meta, by = c("feature_id", "axis"), suffixes = c("", "_ref"))
  nrow(m) == nrow(fast) &&
    max(abs(m$beta_meta - m$beta_meta_ref)) < 1e-8 &&
    max(abs(m$se_meta - m$se_meta_ref)) < 1e-8 &&
    max(abs(m$p_value_meta - m$p_value_meta_ref)) < 1e-8 &&
    max(abs(m$q_value - m$q_value_ref)) < 1e-8 &&
    identical(m$n_cohorts, m$n_cohorts_ref)
})

# --- 4. matched-random draws match size and decile composition exactly ------
check("matched-random draws preserve size and expression-decile composition", {
  bins <- t3_expression_deciles(X)
  universe <- rownames(X)
  members <- split(membership$gene_symbol, membership$feature_id)
  draws <- t3_matched_random_sets(members, universe, bins, n_draws = 25L, seed = 7L)
  all(vapply(names(members), function(f) {
    obs <- match(members[[f]], universe)
    obs_tab <- tabulate(bins[universe[obs]], nbins = 10L)
    all(vapply(draws[[f]], function(d) {
      length(d) == length(obs) && identical(tabulate(bins[universe[d]], nbins = 10L), obs_tab)
    }, logical(1)))
  }, logical(1)))
})

check("matched-random draws are reproducible under the same seed", {
  bins <- t3_expression_deciles(X)
  universe <- rownames(X)
  members <- split(membership$gene_symbol, membership$feature_id)
  a <- t3_matched_random_sets(members, universe, bins, n_draws = 5L, seed = 11L)
  b <- t3_matched_random_sets(members, universe, bins, n_draws = 5L, seed = 11L)
  identical(a, b)
})

check("matched-random draws are not simply the observed set", {
  bins <- t3_expression_deciles(X)
  universe <- rownames(X)
  members <- split(membership$gene_symbol, membership$feature_id)
  d <- t3_matched_random_sets(members, universe, bins, n_draws = 20L, seed = 3L)
  obs <- match(members[["P1"]], universe)
  mean(vapply(d[["P1"]], function(x) length(intersect(x, obs)) / length(obs), numeric(1))) < 0.5
})

# --- 5. the calibration gate actually refuses an anticonservative view ------
check("gate_view refuses a view with an inflated null call rate", {
  bad <- calibrate_view(
    observed_p = runif(50, 0, 0.01),
    permuted_p_fn = function(i) runif(50, 0, 0.001),
    n_reps = 5L, label = "deliberately_anticonservative"
  )
  g <- gate_view(bad, strict = FALSE)
  !bad$reportable && !g$passed
})

check("gate_view accepts a calibrated parametric view", {
  good <- calibrate_view(
    observed_p = c(runif(5, 0, 1e-6), runif(45, 0, 1)),
    permuted_p_fn = function(i) runif(50, 0, 1),
    n_reps = 20L, label = "calibrated"
  )
  good$reportable && gate_view(good, strict = FALSE)$passed
})

# --- 6. Wilson interval agrees with prop.test -------------------------------
check("t3_wilson equals prop.test's Wilson interval", {
  all(vapply(list(c(27, 113), c(17, 47), c(0, 16), c(4, 4)), function(kn) {
    w <- t3_wilson(kn[[1]], kn[[2]])
    p <- suppressWarnings(stats::prop.test(kn[[1]], kn[[2]], correct = FALSE)$conf.int)
    abs(w$lower - p[[1]]) < 1e-8 && abs(w$upper - p[[2]]) < 1e-8
  }, logical(1)))
})

# --- 7. symbol resolution and language ---------------------------------------
check("t3_resolve_symbols maps ENSG ids and leaves symbols alone", {
  ens2sym <- c(ENSG00000141510 = "TP53", ENSG00000012048 = "BRCA1")
  out <- t3_resolve_symbols(c("ENSG00000141510.12", "col1a1", "ENSG00000012048"), ens2sym)
  identical(out, c("TP53", "COL1A1", "BRCA1"))
})

check("assert_language rejects retired ordering vocabulary", {
  isTRUE(tryCatch({ assert_language("stage-associated remodeling"); TRUE }, error = function(e) FALSE)) &&
    isTRUE(tryCatch({ assert_language("disease progression over time"); FALSE },
                    error = function(e) TRUE))
})

# --- 8. the biological-unit assertion is non-vacuous -------------------------
check("assert_biological_unit refuses a unit column copied from sample_id", {
  m <- data.table(sample_id = c("a", "b"), participant_id = c("a", "b"))
  isTRUE(tryCatch({ assert_biological_unit(m, "participant_id"); FALSE },
                  error = function(e) TRUE))
})

cat(sprintf("\nAll %d T3 tests passed.\n", passed))
