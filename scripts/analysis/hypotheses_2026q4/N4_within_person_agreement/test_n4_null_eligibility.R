suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
source(file.path(args[1], "scripts/analysis/hypotheses_2026q4/N4_within_person_agreement/n4_null_eligibility.R"))
d <- data.table(cohort = "A", fibrosis_stage = c(0, 0, 1, 1))
pools <- list(A = data.table(sample_id = letters[1:6], fibrosis_stage = c(0, 0, 0, 1, 1, 1)))
stopifnot(all(n4_null_capacity(d, pools)$sufficient))
set.seed(7); x <- n4_draw_null(d, pools)
stopifnot(!anyDuplicated(x$sample_id), x$substitutions == 0L,
          identical(pools$A$fibrosis_stage[match(x$sample_id, pools$A$sample_id)], d$fibrosis_stage))
set.seed(7); stopifnot(identical(x, n4_draw_null(d, pools)))
bad <- copy(d); bad$fibrosis_stage[1] <- 4
stopifnot(!all(n4_null_capacity(bad, pools)$sufficient),
          inherits(try(n4_draw_null(bad, pools), silent = TRUE), "try-error"))
short <- list(A = pools$A[1:3]); stopifnot(!all(n4_null_capacity(d, short)$sufficient))
cat("N4 null eligibility invariants PASS\n")
