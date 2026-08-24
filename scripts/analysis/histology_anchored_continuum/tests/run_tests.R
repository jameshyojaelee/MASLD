#!/usr/bin/env Rscript

# Execute every isolated R contract in a fresh process so constants and package
# state from one test cannot leak into another.

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
stopifnot(length(script_arg) == 1L)
test_dir <- dirname(normalizePath(script_arg))
analysis_dir <- normalizePath(file.path(test_dir, ".."))
project_root <- normalizePath(file.path(analysis_dir, "../../.."))

tests <- sort(list.files(
  test_dir,
  pattern = "^test_[A-Za-z0-9_]+\\.R$",
  full.names = TRUE
))
stopifnot(length(tests) >= 1L)

failed <- character()
for (test in tests) {
  cat(sprintf("RUN\t%s\n", basename(test)))
  output <- system2(
    file.path(R.home("bin"), "Rscript"),
    c("--vanilla", shQuote(test)),
    env = paste0("MASLD_PROJECT_ROOT=", shQuote(project_root)),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(output, "status")
  if (length(output)) cat(paste0(output, "\n"), sep = "")
  if (!is.null(status) && status != 0L) failed <- c(failed, basename(test))
}

if (length(failed)) {
  stop("Contract tests failed: ", paste(failed, collapse = ", "), call. = FALSE)
}
cat(sprintf("ALL_CONTRACT_TESTS_PASS\t%d\n", length(tests)))
