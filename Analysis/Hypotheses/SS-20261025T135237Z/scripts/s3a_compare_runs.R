# SS S3a restart test: are two S3a output directories identical in every table? (uninterrupted vs killed + resumed)
# Usage: Rscript s3a_compare_runs.R <dir_a> <dir_b> <out_tsv>
suppressMessages(library(data.table))
a <- commandArgs(trailingOnly = TRUE); stopifnot(length(a) == 3L)
f <- sort(list.files(a[1], pattern = "\\.tsv(\\.gz)?$"))
f <- setdiff(f, c("input_sha256.tsv", "frozen_sha256_check.tsv", "script_sha256.tsv"))
r <- rbindlist(lapply(f, function(x) data.table(file = x, identical = file.exists(file.path(a[2], x)) &&
  isTRUE(all.equal(fread(file.path(a[1], x)), fread(file.path(a[2], x)), tolerance = 0)))))
oa <- readRDS(file.path(a[1], "observed_fits.rds")); ob <- readRDS(file.path(a[2], "observed_fits.rds"))
r <- rbind(r, data.table(file = "observed_fits.rds", identical = identical(oa, ob)))
fwrite(r, a[3], sep = "\t"); print(r)
if (!all(r$identical)) stop("RESTART TEST FAILED") else cat("RESTART TEST PASSED:", nrow(r), "outputs identical\n")
