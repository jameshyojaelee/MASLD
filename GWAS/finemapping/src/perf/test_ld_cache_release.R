#!/usr/bin/env Rscript
# Run on a compute node in the rnaseq environment. Uses small synthetic LD blocks.
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})
setDTthreads(1L)

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
stopifnot(length(script_arg) == 1L)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg)))
fixture_dir <- tempfile("ld_cache_release_")
dir.create(fixture_dir)
paths <- file.path(fixture_dir, paste0("block", 1:3))

for (k in seq_along(paths)) {
  n <- 24L
  bim <- data.table(chr = 1L, rsid = paste0("rs", seq_len(n)), dk = 0L,
                    pos = k * 100L + seq_len(n), alt = "A", ref = "C")
  bim[24L, pos := bim$pos[23L]]  # duplicate identity must be removed
  ld <- matrix(0.03, n, n)
  diag(ld) <- 1
  if (k == 3L) {
    ld[2L, ] <- NA_real_
    ld[, 2L] <- NA_real_
  }
  fwrite(bim, paste0(paths[k], ".bim"), sep = "\t", col.names = FALSE)
  fwrite(as.data.table(ld), paste0(paths[k], ".ld"), sep = "\t", col.names = FALSE)
}

cache <- new.env(parent = globalenv())
sys.source(file.path(script_dir, "get_ld_per_locus_cached.R"), envir = cache)
flag <- new.env(parent = emptyenv())
flag$freed <- FALSE
flag$at_read <- logical()
cache$find_LD_block <- function(CHR, START, END, ancestry) paths[START]
cache$fread <- function(file, ...) {
  if (file %in% paste0(paths[2L], c(".bim", ".ld"))) {
    flag$at_read <- c(flag$at_read, flag$freed)
  }
  data.table::fread(file, ...)
}

summary_stats <- function(k, flip = FALSE) {
  data.frame(chromosome = 1L, position = k * 100L + seq_len(23L),
             allele1 = if (flip) "C" else "A", allele2 = if (flip) "A" else "C",
             beta = seq_len(23L) / 100, standard_error = 0.1)
}
invoke <- function(k, flip = FALSE) {
  suppressWarnings(cache$get_ld_per_locus(summary_stats(k, flip), "fixture",
                                         1L, k, k + 1L))
}
check_values <- function(result, k, flip = FALSE, keep = seq_len(23L)) {
  expected_pos <- k * 100L + keep
  expected_snp <- paste(1L, expected_pos, "A", "C", sep = ":")
  expected_ld <- matrix(0.03, length(keep), length(keep),
                        dimnames = list(expected_snp, expected_snp))
  diag(expected_ld) <- 1
  expected_beta <- keep / 100 * if (flip) -1 else 1
  stopifnot(identical(result[[1]]$position, expected_pos),
            identical(result[[1]]$SNP, expected_snp),
            identical(result[[1]]$beta, expected_beta),
            all(result[[1]]$allele1 == "A"), all(result[[1]]$allele2 == "C"),
            identical(result[[2]], expected_ld))
}

first <- invoke(1L)
check_values(first, 1L)
hit <- invoke(1L, flip = TRUE)
check_values(hit, 1L, flip = TRUE)
stopifnot(identical(first[[2]], hit[[2]]),
          !"ss_matched" %in% names(cache$.LD_CACHE$entry$bim))

# A finalizer makes the old cache's lifetime observable without allocating a
# large block. Collection must happen before reading its replacement.
cache$.LD_CACHE$entry$sentinel <- local({
  sentinel <- new.env(parent = emptyenv())
  reg.finalizer(sentinel, local({
    observed <- flag
    function(x) observed$freed <- TRUE
  }), onexit = TRUE)
  sentinel
})
invisible(NULL)  # do not retain the sentinel in R's top-level .Last.value
check_values(invoke(2L), 2L)
stopifnot(length(flag$at_read) == 2L, all(flag$at_read))
check_values(invoke(3L), 3L, keep = seq_len(23L)[-2L])
cat("PASS: cache replacement releases old block before reading; exact LD,",
    "allele direction, duplicate removal, cache hit and NA removal verified.\n")
