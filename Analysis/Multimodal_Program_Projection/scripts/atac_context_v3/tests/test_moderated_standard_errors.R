#!/usr/bin/env Rscript
# Evaluate the actual producer expressions without running either pipeline.
suppressPackageStartupMessages({library(data.table); library(edgeR); library(limma)})
setDTthreads(1L)
set.seed(20261021)
package <- Sys.getenv("MASLD_ATAC_SOURCE_ROOT",
                      "Analysis/Multimodal_Program_Projection/scripts/atac_context_v3")

assignment <- function(path, target, rhs_name) {
  found <- list()
  visit <- function(node) {
    if (!is.call(node)) return(invisible(NULL))
    if (identical(node[[1L]], as.name("<-")) &&
        identical(node[[2L]], as.name(target)) && is.call(node[[3L]]) &&
        identical(node[[3L]][[1L]], as.name(rhs_name))) {
      found[[length(found) + 1L]] <<- node
    }
    parts <- as.list(node)
    for (i in seq_along(parts)[-1L]) {
      if (!identical(parts[[i]], quote(expr = ))) visit(parts[[i]])
    }
  }
  for (node in parse(path)) visit(node)
  stopifnot(length(found) == 1L)
  found[[1L]]
}

check_uncertainty <- function(table, fit, coefficient, effect_column) {
  expected <- fit$stdev.unscaled[, coefficient] * sqrt(fit$s2.post)
  # These are the standard errors used by limma's t and confidence intervals.
  stopifnot(isTRUE(all.equal(unname(table$standard_error), unname(expected),
                            tolerance = 0)))
  stopifnot(max(abs(table[[effect_column]] / table$standard_error -
                    fit$t[, coefficient])) < 1e-10)
  interval <- topTable(fit, coef = coefficient, number = Inf,
                       sort.by = "none", confint = 0.95)
  margin <- qt(0.975, df = fit$df.total) * table$standard_error
  stopifnot(max(abs(interval$CI.R - (table[[effect_column]] + margin))) < 1e-10,
            max(abs(interval$CI.L - (table[[effect_column]] - margin))) < 1e-10)
  ordinary <- fit$stdev.unscaled[, coefficient] * fit$sigma
  stopifnot(any(abs(expected - ordinary) > 1e-6))
}

donors <- data.table(donor_id = paste0("d", 1:12),
                     condition = rep(c("NORMAL", "MASH"), each = 6L))
mu <- matrix(exp(runif(300L * 12L, 3, 6)), nrow = 300L)
counts <- t(matrix(rnbinom(length(mu), mu = mu, size = 4), nrow = 300L))
colnames(counts) <- paste0("peak", seq_len(ncol(counts)))
design <- model.matrix(~ factor(donors$condition, levels = c("NORMAL", "MASH")))
y <- calcNormFactors(DGEList(counts = t(counts)), method = "TMM")
v <- voomWithQualityWeights(y, design, plot = FALSE)
reference <- eBayes(lmFit(v, design), robust = TRUE)
env <- new.env(parent = globalenv())
eval(assignment(file.path(package, "04_run_da.R"), "fit_voom", "function"), env)
actual <- env$fit_voom(counts, donors)
check_uncertainty(actual, reference, 2L, "logFC")
stopifnot(identical(actual$pvalue,
                    topTable(reference, coef = 2L, number = Inf,
                             sort.by = "none")$P.Value))

scores <- matrix(rnorm(64L * 12L), nrow = 64L) * seq(0.1, 3, length.out = 64L)
rownames(scores) <- paste0("program", seq_len(nrow(scores)))
scores <- scores - rowMeans(scores)
fit_expression <- assignment(file.path(package, "05_run_programs.R"), "fit", "eBayes")
env <- list2env(list(score_matrix = scores, design = design), parent = globalenv())
eval(fit_expression, env)
fit <- env$fit
shifted <- scores + seq_len(nrow(scores)) * 1e-12
env <- list2env(list(score_matrix = shifted, design = design), parent = globalenv())
eval(fit_expression, env)
stopifnot(max(abs(fit$p.value[, 2L] - env$fit$p.value[, 2L])) < 1e-10,
          max(abs(fit$coefficients[, 2L] - env$fit$coefficients[, 2L])) < 1e-10)
cat("PROGRAM_OFFSET_INVARIANCE_PASS: arbitrary program offsets do not change disease tests\n")
top <- topTable(fit, coef = 2L, number = Inf, sort.by = "none")
env <- list2env(list(fit = fit, top = top), parent = globalenv())
eval(assignment(file.path(package, "05_run_programs.R"), "fit_table", "data.table"), env)
check_uncertainty(env$fit_table, fit, 2L, "effect")
stopifnot(identical(env$fit_table$effect, top$logFC),
          identical(env$fit_table$pvalue, top$P.Value))
cat("MODERATED_STANDARD_ERROR_REGRESSION_PASS: both producers; t statistics and confidence intervals agree; effects and p-values match each fitted model\n")
