#!/usr/bin/env Rscript
# Run only the source's Monte Carlo p-value expressions, without fitting LMMs
# or writing scientific outputs.
suppressPackageStartupMessages(library(data.table))
source_file <- Sys.getenv("MASLD_CONCORDANCE_SOURCE",
  "Analysis/SingleCell/scripts/v2/353_three_way_concordance.R")
expressions <- as.list(parse(source_file))
text <- vapply(expressions, function(x) paste(deparse(x), collapse = " "), character(1))
pair_expression <- expressions[grepl("^pair_dt\\[, `:=`\\(perm_p_value", text)]
all3_expression <- expressions[vapply(expressions, function(x)
  is.call(x) && identical(x[[1]], as.name("<-")) && identical(x[[2]], as.name("all3_p")), logical(1))]
stopifnot(length(pair_expression) == 1L, length(all3_expression) == 1L)

N_PERMS <- 4L
pair_dt <- data.table(jaccard = c(1, 0, 0.5))
perm_pair <- matrix(rep(c(0, 0.25, 0.5, 0.75), 3), ncol = 3)
eval(pair_expression[[1]])
# Under exchangeability, rank the observed arrangement among all five draws.
oracle <- vapply(seq_len(3), function(k) {
  arrangements <- c(pair_dt$jaccard[k], perm_pair[, k])
  length(which(arrangements >= arrangements[1])) / length(arrangements)
}, numeric(1))
stopifnot(identical(pair_dt$perm_p_value, oracle),
          identical(pair_dt$perm_p_value, c(0.2, 1, 0.6)))
perm_all3 <- c(0, 0.25, 0.5, 0.75)
all3_observed <- 1
eval(all3_expression[[1]])
stopifnot(identical(all3_p, 0.2), all3_p >= 1 / (N_PERMS + 1))
cat("MONTE_CARLO_JACCARD_RESOLUTION_AND_TIES_PASS\n")
