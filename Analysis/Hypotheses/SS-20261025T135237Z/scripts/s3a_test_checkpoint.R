# SS S3a checkpoint unit test (no data): run_blocks() fresh, stopped and resumed, all blocks done, and a changed draw.
# Usage: Rscript s3a_test_checkpoint.R <new dir>   (case "stop" is run in a child process because it calls quit())
source("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/scripts/s3a_common.R")
a <- commandArgs(trailingOnly = TRUE)
if (length(a) == 2L && a[1] == "--child-stop") {  # child: write block 1 of 3, then quit(status = 3)
  run_blocks("t", 1:10, function(k) c(v = k^2), function(ids) ids * 10, a[2], 4L, stop_after = 1L); quit(status = 0) }
d <- new_outdir(a[1]); fun <- function(k) c(v = k^2); key <- function(ids) ids * 10
ref <- run_blocks("t", 1:10, fun, key, new_outdir(file.path(d, "fresh")), 4L)
stopifnot(length(ref) == 10L, identical(unname(sapply(ref, `[[`, "v")), (1:10)^2))
ck <- new_outdir(file.path(d, "stopped"))
rc <- system2(file.path(R.home("bin"), "Rscript"), c(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "--child-stop", ck))
stopifnot(rc == 3L, identical(list.files(ck), "t_block001.rds"))
res <- run_blocks("t", 1:10, fun, key, ck, 4L)
stopifnot(identical(res, ref), length(list.files(ck)) == 3L)
res2 <- run_blocks("t", 1:10, function(k) stop("must not refit"), key, ck, 4L)  # all blocks done: nothing refit
stopifnot(identical(res2, ref))
bad <- tryCatch(run_blocks("t", 1:10, fun, function(ids) ids * 11, ck, 4L), error = function(e) conditionMessage(e))
stopifnot(is.character(bad), grepl("does not match", bad))
nl <- run_blocks("n", 1:7, function(k) if (k %in% c(5L, 7L)) NULL else k, function(ids) ids, new_outdir(file.path(d, "nulls")), 4L)
stopifnot(length(nl) == 7L, is.null(nl[[5]]), is.null(nl[[7]]), identical(nl[[6]], 6L))  # skipped fits keep their slot
cat("CHECKPOINT UNIT TEST PASSED: fresh, stop (exit 3) + resume, all-done reuse, changed-draw refusal, NULL results kept\n")
