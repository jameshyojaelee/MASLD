#!/usr/bin/env Rscript
# Unit tests for the molecular-system stress-test helpers.

suppressPackageStartupMessages({
  library(data.table)
  library(igraph)
})
root <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
script_dir <- file.path(root, "scripts/analysis/histology_anchored_continuum/molecular_layers")
source(file.path(script_dir, "lib_molecular_layers.R"))

# --- collection-balanced scoring -------------------------------------------
# Median within collection, then unweighted mean across collections, so a
# collection contributing many features must not outweigh a sparse one.
scores <- data.table(
  sample_id = rep(c("s1", "s2"), each = 5L),
  dataset = "D",
  community_id = "system_01",
  collection = rep(c("go_bp", "go_bp", "go_bp", "hallmark", "hotspot"), 2L),
  feature_score = c(1, 2, 3, 10, 20, 2, 4, 6, 30, 50)
)
balanced <- ml_collection_balanced_scores(scores, "community_id")
setorder(balanced, sample_id)
stopifnot(
  nrow(balanced) == 2L,
  all(balanced$n_collections == 3L),
  abs(balanced$system_score_raw[[1L]] - (2 + 10 + 20) / 3) < 1e-12,
  abs(balanced$system_score_raw[[2L]] - (4 + 30 + 50) / 3) < 1e-12
)

# --- Li and Ji effective number of tests ------------------------------------
stopifnot(abs(ml_effective_tests(diag(10L)) - 10) < 1e-8)
stopifnot(ml_effective_tests(matrix(1, nrow = 10L, ncol = 10L)) < 1.5)
set.seed(20260817)
loading <- rnorm(20L, mean = 0.8, sd = 0.05)
correlated <- outer(loading, loading)
diag(correlated) <- 1
effective <- ml_effective_tests(correlated)
stopifnot(effective > 1, effective < 20)

# --- connected-subgraph sampler ---------------------------------------------
ring <- make_ring(50L)
adjacency <- lapply(as_adj_list(ring, mode = "all"), as.integer)
set.seed(20260817)
for (attempt in seq_len(20L)) {
  drawn <- ml_draw_connected_subgraph(adjacency, 12L)
  stopifnot(
    length(drawn) == 12L,
    !anyDuplicated(drawn),
    all(drawn >= 1L & drawn <= 50L),
    count_components(induced_subgraph(ring, drawn)) == 1L
  )
}
# A star graph exercises the boundary bookkeeping: every leaf is reachable only
# through the hub, so the sampler must not revisit a node or stall.
star <- make_star(30L, mode = "undirected")
star_adjacency <- lapply(as_adj_list(star, mode = "all"), as.integer)
set.seed(1L)
star_draws <- replicate(20L, {
  drawn <- ml_draw_connected_subgraph(star_adjacency, 5L)
  isTRUE(!is.null(drawn) && count_components(induced_subgraph(star, drawn)) == 1L)
})
stopifnot(all(star_draws))
# Whole-graph draws must return every node exactly once.
stopifnot(setequal(ml_draw_connected_subgraph(adjacency, 50L), seq_len(50L)))

# --- residualized fit equals lm ---------------------------------------------
# The stress tests replace lm() with Frisch-Waugh residualization for speed;
# the coefficient and standard error must be identical, not merely close.
set.seed(4242L)
n <- 140L
stage <- sample(0:4, n, replace = TRUE)
sex <- sample(c("F", "M"), n, replace = TRUE)
x <- rnorm(n)
y <- 0.4 * x + 0.3 * stage + rnorm(n)
nuisance <- model.matrix(~ factor(stage) + factor(sex))
nuisance_qr <- qr(nuisance)
x_residual <- as.numeric(qr.resid(nuisance_qr, x))
x_residual <- x_residual - mean(x_residual)
sxx <- sum(x_residual^2)
residual_df <- n - ncol(nuisance) - 1L
y_residual <- as.numeric(qr.resid(nuisance_qr, y))
y_residual <- y_residual - mean(y_residual)
beta <- sum(x_residual * y_residual) / sxx
rss <- sum(y_residual^2) - beta^2 * sxx
se <- sqrt(rss / residual_df / sxx)
reference <- coef(summary(lm(y ~ x + factor(stage) + factor(sex))))["x", ]
stopifnot(
  abs(beta - reference[["Estimate"]]) < 1e-10,
  abs(se - reference[["Std. Error"]]) < 1e-10
)

# --- script contract ---------------------------------------------------------
stress_text <- readLines(
  file.path(script_dir, "73_stress_test_molecular_systems.R"), warn = FALSE
)
stopifnot(
  any(grepl("Refusing to overwrite output candidate", stress_text, fixed = TRUE)),
  any(grepl("Promoted atlas source hash drift", stress_text, fixed = TRUE)),
  any(grepl("Frozen inference source hash drift", stress_text, fixed = TRUE)),
  any(grepl("connected_subgraph", stress_text, fixed = TRUE)),
  any(grepl("size_matched_uniform", stress_text, fixed = TRUE)),
  any(grepl("joint_maxT_fwer_p_value", stress_text, fixed = TRUE)),
  any(grepl("permutation_maxT_supported_all_four", stress_text, fixed = TRUE)),
  any(grepl("Production stress tests must use", stress_text, fixed = TRUE)),
  !any(grepl("longitudinal", stress_text, fixed = TRUE))
)

cat("SYSTEM_STRESS_TESTS_PASS\n")
