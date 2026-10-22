#!/usr/bin/env Rscript
# Exercise the actual assignment in each producer without running the pipeline.
root <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
base <- file.path(root, "scripts/analysis/histology_anchored_continuum")
specs <- list(
  list(path = file.path(base, "04_score_programs_and_validate.R"),
       lhs = quote(rx_perm[g]), input = "rx", output = "rx_perm", index = "g"),
  list(path = file.path(base, "molecular_layers/70_infer_molecular_systems.R"),
       lhs = quote(permuted[group]), input = "x_residual", output = "permuted", index = "group"),
  list(path = file.path(base, "molecular_layers/73_stress_test_molecular_systems.R"),
       lhs = quote(values[group]), input = "values", output = "values", index = "group")
)

assignments_for <- function(expressions, lhs) {
  found <- list()
  visit <- function(node) {
    if (missing(node)) return(invisible(NULL))
    if (is.call(node) && identical(node[[1L]], as.name("<-")) &&
        identical(node[[2L]], lhs)) found[[length(found) + 1L]] <<- node
    if (is.call(node) || is.expression(node) || is.list(node)) {
      parts <- as.list(node)
      for (index in seq_along(parts)) {
        if (!identical(parts[[index]], quote(expr = ))) visit(parts[[index]])
      }
    }
  }
  visit(expressions)
  stopifnot(length(found) == 1L)
  found[[1L]]
}

for (spec in specs) {
  assignment <- assignments_for(parse(spec$path), spec$lhs)
  run <- function(values, groups) {
    env <- new.env(parent = baseenv())
    assign(spec$input, values, env)
    assign(spec$output, values, env)
    for (group in groups) {
      assign(spec$index, group, env)
      eval(assignment, env)
    }
    get(spec$output, env)
  }
  # A singleton must keep its original value, including a positive noninteger
  # that sample(x, 1) would interpret as a population size.
  for (value in c(1.7, 5, -2, 0.25, NA_real_)) {
    set.seed(20261021L)
    stopifnot(identical(run(value, list(1L)), value))
  }
  # Existing nonsingleton paths must keep both the permutation and RNG state.
  values <- c(1.7, -4, 2.3, 0.25, 9, -3, 8.5, 0)
  groups <- list(1:3, 4:6, 7:8)
  for (seed in seq_len(100L)) {
    set.seed(seed)
    expected <- values
    for (group in groups) expected[group] <- sample(values[group], length(group), replace = FALSE)
    expected_rng <- .Random.seed
    set.seed(seed)
    observed <- run(values, groups)
    stopifnot(identical(observed, expected), identical(.Random.seed, expected_rng))
  }
  # Never move values across stage-by-sex strata, including empty groups.
  mixed_groups <- list(integer(), 1L, 2:4, 5L, 6:8)
  set.seed(42L)
  observed <- run(values, mixed_groups)
  for (group in mixed_groups) {
    stopifnot(identical(sort(observed[group]), sort(values[group])))
  }
}
cat("SINGLETON_PERMUTATION_REGRESSION_PASS: three producers; singleton values retained; nonsingleton outputs and RNG identical\n")
