# Cross-assay modules: shared library.
#
# The question. Every gene set in this Resource is discovered in one assay and
# then projected into the others, where most of it is not measured: of the 117
# frozen Hotspot programs, 26 are testable in the liver proteome, 12 in CosMx,
# and none reach a supported state in snATAC. A coverage failure and a real
# negative then look the same. This family inverts the order: the gene universe
# is fixed first to genes measured in bulk RNA, snRNA, Visium and the liver
# proteome, modules are built inside it, and every assay therefore has a real
# test to fail.
#
# TWO RULES THAT SHAPE THE CODE.
#
# 1. Modules are built LABEL-BLIND on co-expression alone and frozen before any
#    outcome is read. Selecting genes on stage association first would make the
#    external bulk cohorts re-ask a question already answered in the discovery
#    cohorts, and would leave no module without discovery support to act as a
#    negative control.
# 2. A random gene set is the wrong competitor for a co-expressed module: on
#    this substrate, expression-matched random sets fail to reach a real
#    module's coherence in half the arms tested. The primary competitive null
#    is therefore random CONNECTED subgraphs of the same co-expression graph,
#    size-matched. The expression-decile-matched draw is retained as the weaker
#    reference, not as the deciding one.
#
# LANGUAGE. Scores order donors within a cohort. Nothing here is comparable in
# level between cohorts or between assays, and no quantity combines assays.

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

CAM_PROJECT_ROOT <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)

# Reuse the audited helpers rather than reimplementing them. Sourcing only
# defines functions; neither file reads a contract at load time.
source(file.path(CAM_PROJECT_ROOT,
  "scripts/analysis/histology_anchored_continuum/molecular_layers/lib_molecular_layers.R"))
source(file.path(CAM_PROJECT_ROOT,
  "scripts/analysis/histology_anchored_continuum/lib_continuum.R"))
# t3_expression_deciles / t3_matched_random_sets: the expression-matched draw
# already used for the program-vocabulary nulls, reused rather than rewritten.
source(file.path(CAM_PROJECT_ROOT,
  "scripts/manuscript/program_observability_map/t3_vocabulary_lib.R"))

# ---------------------------------------------------------------------------
# Contract, paths, write-once IO
# ---------------------------------------------------------------------------

cam_fail <- function(...) stop(paste0(...), call. = FALSE)

cam_assert <- function(condition, message) {
  if (length(condition) != 1L || is.na(condition) || !condition) cam_fail(message)
  invisible(TRUE)
}

cam_script_dir <- function() file.path(CAM_PROJECT_ROOT, "scripts/analysis/cross_assay_modules")

cam_contract <- function() {
  path <- file.path(cam_script_dir(), "00_contract.json")
  cam_assert(file.exists(path), paste0("Missing contract: ", path))
  jsonlite::fromJSON(path, simplifyVector = TRUE)
}

cam_out_root <- function(must_exist = TRUE) {
  value <- Sys.getenv("CAM_OUT_ROOT", unset = "")
  cam_assert(nzchar(value), "CAM_OUT_ROOT is unset")
  if (!dir.exists(value) && !must_exist) dir.create(value, recursive = TRUE, showWarnings = FALSE)
  normalizePath(value, mustWork = must_exist)
}

cam_input <- function(key, contract = cam_contract(), must_exist = TRUE) {
  rel <- contract$inputs[[key]]
  cam_assert(!is.null(rel) && nzchar(rel), paste0("Contract has no input named ", key))
  path <- if (startsWith(rel, "/")) rel else file.path(CAM_PROJECT_ROOT, rel)
  if (must_exist) cam_assert(file.exists(path), paste0("Missing input ", key, ": ", path))
  path
}

cam_dir <- function(...) {
  path <- file.path(cam_out_root(must_exist = FALSE), ...)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

# Write-once. A rerun that would change a recorded number must land in a new
# output root rather than silently overwrite the one a figure already cites.
cam_write_tsv <- function(x, path) {
  cam_assert(!file.exists(path), paste0("Refusing to overwrite existing file: ", path))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(x, path, sep = "\t", quote = FALSE, na = "NA")
  invisible(path)
}

cam_write_json <- function(x, path) {
  cam_assert(!file.exists(path), paste0("Refusing to overwrite existing file: ", path))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(jsonlite::toJSON(x, auto_unbox = TRUE, pretty = TRUE, digits = 12), path)
  invisible(path)
}

cam_write_rds <- function(x, path) {
  cam_assert(!file.exists(path), paste0("Refusing to overwrite existing file: ", path))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(x, path, compress = "xz")
  invisible(path)
}

# sha256 via the system utility so values match the sha256sum records used
# elsewhere in the repository.
cam_sha256 <- function(path) {
  cam_assert(file.exists(path), paste0("Cannot hash a missing file: ", path))
  out <- system2("sha256sum", shQuote(path), stdout = TRUE)
  sub("\\s.*$", "", out[1])
}

cam_say <- function(...) {
  message(format(Sys.time(), "[%H:%M:%S] "), paste0(...))
}

# ---------------------------------------------------------------------------
# Guardrails
# ---------------------------------------------------------------------------

# The repository's spatial contract raises on these names because a column
# called combined_score is a cross-assay ranking whatever its contents are.
# Applying the same list here keeps this family inside the same rule.
cam_assert_no_prohibited_columns <- function(x, contract = cam_contract()) {
  bad <- intersect(tolower(names(x)), tolower(contract$prohibited_columns))
  cam_assert(length(bad) == 0L,
             paste0("Prohibited cross-assay column(s): ", paste(bad, collapse = ", ")))
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------

# Gene-wise z within a single cohort/assay block. Rows with no variance carry
# no information about donor ordering and are dropped rather than set to zero,
# which would dilute the mean toward an arbitrary centre.
cam_zscore_rows <- function(m) {
  mu <- rowMeans(m, na.rm = TRUE)
  sdv <- apply(m, 1L, stats::sd, na.rm = TRUE)
  keep <- is.finite(sdv) & sdv > 0
  z <- (m[keep, , drop = FALSE] - mu[keep]) / sdv[keep]
  z
}

# Equal-weight module score inside one cohort block, then standardised across
# donors of that block. Equal weights are deliberate: these modules have no
# native weight, and weighting one family but not its null would confound
# "better genes" with "has weights".
cam_score_block <- function(z, members, direction = 1) {
  idx <- intersect(members, rownames(z))
  if (length(idx) == 0L) return(rep(NA_real_, ncol(z)))
  s <- colMeans(z[idx, , drop = FALSE], na.rm = TRUE) * direction
  standardize_vector(s)
}

# Coverage decision, recorded per module per assay so that "not measured" and
# "measured and unsupported" never collapse into one cell.
cam_testability <- function(members, measured, contract = cam_contract()) {
  n_meas <- length(intersect(members, measured))
  frac <- if (length(members)) n_meas / length(members) else 0
  list(
    n_members = length(members),
    n_measured = n_meas,
    fraction_measured = frac,
    testable = n_meas >= contract$scoring$testability$min_members_measured &&
      frac >= contract$scoring$testability$min_fraction_measured
  )
}

# ---------------------------------------------------------------------------
# Association models
# ---------------------------------------------------------------------------

# One linear fit of a standardised module score on an ordinal endpoint with
# covariates. Returns the endpoint coefficient only. The endpoint enters as a
# recorded cross-sectional grade, never as time.
cam_fit_score <- function(dt, score_col, endpoint_col, covariates = character(0)) {
  keep <- is.finite(dt[[score_col]]) & is.finite(dt[[endpoint_col]])
  for (cv in covariates) keep <- keep & !is.na(dt[[cv]])
  d <- dt[keep]
  if (nrow(d) < 10L) return(list(estimable = FALSE, beta = NA_real_, se = NA_real_,
                                 p_two_sided = NA_real_, n = nrow(d)))
  terms <- c(endpoint_col, covariates)
  usable <- terms[vapply(terms, function(tt) length(unique(d[[tt]])) > 1L, logical(1))]
  if (!(endpoint_col %in% usable)) {
    return(list(estimable = FALSE, beta = NA_real_, se = NA_real_,
                p_two_sided = NA_real_, n = nrow(d)))
  }
  rhs <- paste(vapply(usable, function(tt) {
    if (is.character(d[[tt]]) || is.factor(d[[tt]])) paste0("factor(", tt, ")") else tt
  }, character(1)), collapse = " + ")
  fit <- try(stats::lm(stats::as.formula(paste(score_col, "~", rhs)), data = d), silent = TRUE)
  if (inherits(fit, "try-error")) {
    return(list(estimable = FALSE, beta = NA_real_, se = NA_real_,
                p_two_sided = NA_real_, n = nrow(d)))
  }
  cf <- summary(fit)$coefficients
  if (!(endpoint_col %in% rownames(cf))) {
    return(list(estimable = FALSE, beta = NA_real_, se = NA_real_,
                p_two_sided = NA_real_, n = nrow(d)))
  }
  list(estimable = TRUE,
       beta = unname(cf[endpoint_col, 1]),
       se = unname(cf[endpoint_col, 2]),
       p_two_sided = unname(cf[endpoint_col, 4]),
       n = nrow(d), df = fit$df.residual)
}

cam_fit_empty <- function(n = 0L) list(estimable = FALSE, beta = NA_real_, se = NA_real_,
                                       p_two_sided = NA_real_, n = n, df = NA_integer_)

# Endpoint standardised within one cohort/assay block over the participants a
# fit would use (finite endpoint, complete covariates). Every reported effect
# is then a standardized slope: SD of score per SD of endpoint. This does not
# make a fibrosis grade and a NAS the same estimand; it makes their slopes
# comparable in scale, and the endpoint name travels with every cell.
cam_endpoint_z <- function(y, covar = NULL) {
  keep <- is.finite(y)
  if (!is.null(covar)) keep <- keep & stats::complete.cases(as.data.frame(covar))
  out <- rep(NA_real_, length(y))
  if (sum(keep) < 2L || stats::sd(y[keep]) == 0) return(out)
  out[keep] <- (y[keep] - mean(y[keep])) / stats::sd(y[keep])
  out
}

# Equivalence-type one-sided p for H0: oriented slope >= sesoi against
# H1: oriented slope < sesoi, on the fit's own t distribution. A small p says
# the assay excludes an effect of at least `sesoi` in the discovery direction.
# This replaces the minimum-detectable-effect rule of v1, which compared a power
# quantity to a discovery slope in different endpoint units and could not
# establish absence.
cam_equivalence_p <- function(beta_oriented, se, df, sesoi) {
  n <- length(beta_oriented); se <- rep_len(se, n); df <- rep_len(df, n)
  ok <- is.finite(beta_oriented) & is.finite(se) & se > 0 & is.finite(df) & df > 0
  out <- rep(NA_real_, length(beta_oriented))
  out[ok] <- stats::pt((beta_oriented[ok] - sesoi) / se[ok], df = df[ok], lower.tail = TRUE)
  out
}

# One-sided 95 percent upper bound of the oriented slope, printed for every
# non-significant cell so "not resolved" is readable as a bound, not a blank.
cam_upper_bound <- function(beta_oriented, se, df, level = 0.95) {
  n <- length(beta_oriented); se <- rep_len(se, n); df <- rep_len(df, n)
  ok <- is.finite(beta_oriented) & is.finite(se) & is.finite(df) & df > 0
  out <- rep(NA_real_, length(beta_oriented))
  out[ok] <- beta_oriented[ok] + stats::qt(level, df = df[ok]) * se[ok]
  out
}

cam_partial_spearman <- function(x, y, Z = NULL) {
  ok <- is.finite(x) & is.finite(y)
  if (!is.null(Z)) {
    Z <- as.matrix(Z)
    ok <- ok & stats::complete.cases(Z)
  }
  if (sum(ok) < 10L) return(NA_real_)
  rx <- rank(x[ok]); ry <- rank(y[ok])
  if (is.null(Z)) return(suppressWarnings(stats::cor(rx, ry)))
  Zr <- apply(Z[ok, , drop = FALSE], 2L, rank)
  Zr <- Zr[, apply(Zr, 2L, function(v) length(unique(v)) > 1L), drop = FALSE]
  if (ncol(Zr) == 0L) return(suppressWarnings(stats::cor(rx, ry)))
  ex <- stats::residuals(stats::lm(rx ~ Zr))
  ey <- stats::residuals(stats::lm(ry ~ Zr))
  suppressWarnings(stats::cor(ex, ey))
}

# ---------------------------------------------------------------------------
# Nulls
# ---------------------------------------------------------------------------

# Permute the endpoint inside strata so the null keeps the covariate structure
# the real test conditions on. Permuting across strata would inflate the
# observed statistic by reintroducing between-stratum contrast.
cam_permute_within <- function(values, strata, rng_seed) {
  set.seed(rng_seed)
  out <- values
  for (g in split(seq_along(values), strata)) {
    if (length(g) > 1L) out[g] <- values[sample(g)]
  }
  out
}

cam_empirical_p_ge <- function(observed, null_values) {
  nv <- null_values[is.finite(null_values)]
  if (!is.finite(observed) || length(nv) == 0L) return(NA_real_)
  (1 + sum(nv >= observed)) / (1 + length(nv))
}

# Size-matched connected subgraphs of the same co-expression graph. This is the
# competitive reference: it asks whether a coherent module of this size scores
# this well by construction on this substrate, which an unstructured random set
# cannot answer.
cam_draw_connected_sets <- function(adjacency, target_size, n_draws, seed) {
  set.seed(seed)
  out <- vector("list", n_draws)
  failures <- 0L
  for (i in seq_len(n_draws)) {
    drawn <- ml_draw_connected_subgraph(adjacency, target_size)
    if (is.null(drawn)) { failures <- failures + 1L; next }
    out[[i]] <- drawn
  }
  list(sets = Filter(Negate(is.null), out), n_failed = failures,
       failure_rate = failures / n_draws)
}

# ---------------------------------------------------------------------------
# Evidence states
# ---------------------------------------------------------------------------

# The only place states are assigned. Two separate families:
#
#   association_state answers "does the standardized association replicate in
#   the frozen discovery direction", from a two-sided BH q and the SIGN of the
#   oriented effect. Coverage is decided first, and an assay with no endpoint of
#   the module's family is not_applicable rather than silent.
#
#   specificity_state answers "does it exceed co-expressed sets of the same
#   size", from a BH-corrected competitive q. It is never folded into the first.
#
# v1 applied direction twice and tested one tail, which made `discordant`
# unreachable and mislabelled reversals; every branch here is exercised by the
# self-tests in 15 before any real row is written.
cam_association_state <- function(applicable, testable, q_two_sided, effect_oriented,
                                  equivalence_q, cap = NA_character_,
                                  q_threshold = 0.05) {
  state <- if (!isTRUE(applicable)) {
    "not_applicable"
  } else if (!isTRUE(testable)) {
    "untestable"
  } else if (!is.finite(q_two_sided) || !is.finite(effect_oriented)) {
    "indeterminate"
  } else if (q_two_sided < q_threshold && effect_oriented > 0) {
    "replicates"
  } else if (q_two_sided < q_threshold && effect_oriented < 0) {
    "discordant"
  } else if (is.finite(equivalence_q) && equivalence_q < q_threshold) {
    "tested_negative"
  } else {
    "indeterminate"
  }
  if (!is.na(cap) && state == "replicates") state <- cap
  state
}

cam_specificity_state <- function(testable, competitive_q, q_threshold = 0.05) {
  if (!isTRUE(testable) || !is.finite(competitive_q)) return("not_assessed")
  if (competitive_q < q_threshold) "exceeds_coexpressed_sets" else "does_not_exceed"
}

# Non-significance against matched genes is NOT absence of organisation: the
# test carries no bound on spatial excess, and most non-significant modules sit
# above their null mean. The state is therefore `coherence_unresolved`.
cam_coherence_state <- function(testable, q_value, source_dependent = FALSE,
                                q_threshold = 0.05) {
  if (!isTRUE(testable)) return("untestable")
  if (isTRUE(source_dependent)) return("source_dependent")
  if (!is.finite(q_value)) return("coherence_unresolved")
  if (q_value < q_threshold) "coherent" else "coherence_unresolved"
}

# Synthetic checks that the rule does what its documentation says. 15 refuses
# to run unless every one of these holds.
cam_state_self_test <- function() {
  # strong reverse effect must be discordant, whatever the discovery direction
  cam_assert(cam_association_state(TRUE, TRUE, 1e-6, -1 * 1, NA) == "discordant",
             "self-test: reverse effect did not print discordant")
  # a decreasing module (direction -1) whose raw slope is negative agrees
  raw <- -0.5; dir <- -1
  cam_assert(cam_association_state(TRUE, TRUE, 1e-6, raw * dir, NA) == "replicates",
             "self-test: agreement for a direction -1 module did not print replicates")
  # a decreasing module whose raw slope is positive is a reversal
  cam_assert(cam_association_state(TRUE, TRUE, 1e-6, 0.5 * dir, NA) == "discordant",
             "self-test: reversal for a direction -1 module did not print discordant")
  # tight null excludes the SESOI
  cam_assert(cam_association_state(TRUE, TRUE, 0.9, 0.01, 0.001) == "tested_negative",
             "self-test: bound below SESOI did not print tested_negative")
  # wide null cannot exclude it
  cam_assert(cam_association_state(TRUE, TRUE, 0.9, 0.01, 0.6) == "indeterminate",
             "self-test: wide null did not print indeterminate")
  cam_assert(cam_association_state(TRUE, FALSE, 1e-6, 1, 1e-6) == "untestable",
             "self-test: coverage failure did not print untestable")
  cam_assert(cam_association_state(FALSE, TRUE, 1e-6, 1, 1e-6) == "not_applicable",
             "self-test: missing endpoint did not print not_applicable")
  cam_assert(cam_association_state(TRUE, TRUE, 1e-6, 1, NA, cap = "indeterminate") == "indeterminate",
             "self-test: cap did not apply")
  cam_assert(cam_specificity_state(TRUE, 0.01) == "exceeds_coexpressed_sets" &&
             cam_specificity_state(TRUE, 0.2) == "does_not_exceed" &&
             cam_specificity_state(FALSE, 0.01) == "not_assessed",
             "self-test: specificity states wrong")
  cam_assert(cam_coherence_state(TRUE, 0.01, source_dependent = TRUE) == "source_dependent",
             "self-test: Vu cap did not apply to a significant cell")
  cam_assert(cam_coherence_state(TRUE, 0.4) == "coherence_unresolved",
             "self-test: spatial non-significance printed as absence")
  # equivalence p: a slope far below the SESOI with a tiny SE gives p ~ 0
  cam_assert(cam_equivalence_p(0, 0.01, 50, 0.2) < 1e-6,
             "self-test: equivalence p not small for an excluded effect")
  cam_assert(cam_equivalence_p(0.2, 0.01, 50, 0.2) > 0.49,
             "self-test: equivalence p at the SESOI is not ~0.5")
  invisible(TRUE)
}

isTRUE_vec <- function(x) !is.na(x) & as.logical(x)

cam_session_info <- function(path) {
  if (file.exists(path)) return(invisible(path))
  writeLines(utils::capture.output(utils::sessionInfo()), path)
  invisible(path)
}

# ---------------------------------------------------------------------------
# Vectorised association, for the nulls
# ---------------------------------------------------------------------------

# The competitive null needs the same fit repeated for thousands of gene sets.
# Fitting them one at a time with lm() is the difference between minutes and
# days, so the endpoint coefficient is obtained by the Frisch-Waugh-Lovell
# route: residualise the endpoint and the scores on the covariates once, then
# every set's coefficient is one inner product. The result is algebraically the
# coefficient lm() would report, and 05 checks that on the observed modules.
#
# S: sets-by-samples score matrix. y: endpoint. X: covariate design WITHOUT an
# intercept column (one is added here).
cam_fast_assoc <- function(S, y, X = NULL) {
  n <- length(y)
  design <- if (is.null(X)) matrix(1, n, 1) else cbind(1, as.matrix(X))
  qrd <- qr(design)
  resid_on <- function(v) v - design %*% qr.coef(qrd, v)
  e <- as.vector(resid_on(y))
  see <- sum(e^2)
  if (!is.finite(see) || see <= 0) {
    return(list(beta = rep(NA_real_, nrow(S)), se = rep(NA_real_, nrow(S)),
                p_two_sided = rep(NA_real_, nrow(S)), n = n))
  }
  Sr <- t(resid_on(t(S)))
  beta <- as.vector(Sr %*% e) / see
  fitted <- beta %o% e
  rss <- rowSums((Sr - fitted)^2)
  df <- n - qrd$rank - 1L
  sigma2 <- rss / df
  se <- sqrt(sigma2 / see)
  tstat <- beta / se
  list(beta = beta, se = se,
       p_two_sided = 2 * stats::pt(abs(tstat), df = df, lower.tail = FALSE),
       n = n, df = df)
}

# Isometric log-ratio coordinates for a composition, via normalised Helmert
# contrasts. Proportions are not independent coordinates; regressing on the raw
# 16 lineage fractions would be singular, and on 15 of them would make the
# answer depend on which lineage was dropped.
cam_ilr <- function(P, floor_value = 1e-6) {
  P <- as.matrix(P)
  P[!is.finite(P)] <- 0
  P <- P / rowSums(P)
  P[P < floor_value] <- floor_value
  P <- P / rowSums(P)
  L <- log(P)
  D <- ncol(P)
  helmert <- stats::contr.helmert(D)
  helmert <- t(t(helmert) / sqrt(colSums(helmert^2)))
  L %*% helmert
}

# Connected size-matched draws restricted to components that can actually hold
# a set of the requested size. Seeding uniformly over all nodes would score a
# 16 percent failure rate against the graph rather than against the estimator:
# an isolated gene cannot begin a connected set of eight, and counting those
# attempts as failures would say nothing about how special the real module is.
cam_connected_draws <- function(graph, target_size, n_draws, seed) {
  comp <- igraph::components(graph)
  eligible_comps <- which(comp$csize >= target_size)
  nodes <- which(comp$membership %in% eligible_comps)
  if (length(nodes) < target_size) {
    return(list(sets = list(), n_failed = n_draws, failure_rate = 1,
                n_eligible_nodes = length(nodes)))
  }
  sub <- igraph::induced_subgraph(graph, nodes)
  adjacency <- lapply(igraph::as_adj_list(sub, mode = "all"), as.integer)
  names_sub <- igraph::V(sub)$name
  set.seed(seed)
  sets <- vector("list", n_draws); failed <- 0L
  for (i in seq_len(n_draws)) {
    drawn <- ml_draw_connected_subgraph(adjacency, target_size)
    if (is.null(drawn)) { failed <- failed + 1L; next }
    sets[[i]] <- names_sub[drawn]
  }
  list(sets = Filter(Negate(is.null), sets), n_failed = failed,
       failure_rate = failed / n_draws, n_eligible_nodes = length(nodes))
}

# Mutual-kNN graph on positive Spearman, WITH the edge weights 03 attaches to
# the frozen graph. 04's leave-one-cohort-out rebuilds use this so the rebuilt
# graphs differ from the full graph only by the cohort removed, not by their
# weighting. 03 keeps its own inline construction so the frozen graph is
# reproduced byte for byte.
cam_build_mutual_graph <- function(zz, k) {
  rho <- stats::cor(t(zz), method = "spearman")
  diag(rho) <- 0; rho[rho < 0] <- 0
  nn <- t(apply(rho, 1L, function(r) order(r, decreasing = TRUE)[seq_len(k)]))
  flag <- matrix(FALSE, nrow(rho), nrow(rho))
  for (i in seq_len(nrow(rho))) flag[i, nn[i, ]] <- TRUE
  m <- flag & t(flag); m[rho <= 0] <- FALSE
  gg <- igraph::graph_from_adjacency_matrix(m, mode = "undirected", diag = FALSE)
  igraph::V(gg)$name <- rownames(zz)
  igraph::E(gg)$weight <- rho[igraph::as_edgelist(gg, names = FALSE)]
  gg
}

# Mean pairwise Spearman inside a gene set, from a precomputed rho matrix.
# Used to compare the coherence of modules with that of their competitors.
cam_set_coherence <- function(rho, genes) {
  idx <- match(genes, rownames(rho)); idx <- idx[!is.na(idx)]
  if (length(idx) < 2L) return(NA_real_)
  sub <- rho[idx, idx, drop = FALSE]
  mean(sub[upper.tri(sub)])
}

# Inverse-variance fixed-effect meta over cohorts, vectorised across sets.
# beta_mat and se_mat are sets-by-cohorts.
cam_meta_rows <- function(beta_mat, se_mat, min_cohorts = 2L) {
  w <- 1 / se_mat^2
  ok <- is.finite(beta_mat) & is.finite(se_mat) & se_mat > 0
  w[!ok] <- 0; b <- beta_mat; b[!ok] <- 0
  sw <- rowSums(w)
  est <- rowSums(w * b) / sw
  est_se <- sqrt(1 / sw)
  n_ok <- rowSums(ok)
  est[n_ok < min_cohorts] <- NA_real_
  est_se[n_ok < min_cohorts] <- NA_real_
  z <- est / est_se
  list(beta = est, se = est_se, z = z,
       p_two_sided = 2 * stats::pnorm(abs(z), lower.tail = FALSE),
       n_cohorts = n_ok)
}

# Fast partial Spearman against a FIXED covariate block. The rank-residual
# projection is built once and reused, which is what makes 2,000 Freedman-Lane
# draws per module affordable; the value is identical to cam_partial_spearman.
cam_partial_engine <- function(Z) {
  Zr <- apply(as.matrix(Z), 2L, rank)
  Zr <- Zr[, apply(Zr, 2L, function(v) length(unique(v)) > 1L), drop = FALSE]
  design <- cbind(1, Zr)
  qrd <- qr(design)
  list(
    resid = function(v) as.vector(v - design %*% qr.coef(qrd, v)),
    design = design
  )
}

cam_partial_from_engine <- function(engine, x_rank, y_resid) {
  ex <- engine$resid(x_rank)
  suppressWarnings(stats::cor(ex, y_resid))
}
