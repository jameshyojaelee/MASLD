# discriminate_cache_identity.R
#
# THE QUESTION: memoizing the LD block load changed lambda_s_locus by up to
# 4.1e-14. Two mutually exclusive explanations, with very different consequences:
#
#   (A) The cache returns a DIFFERENT MATRIX than a fresh read. That would be a
#       real defect -- the cache would be corrupting data, and must not ship,
#       regardless of how small today's effect looks.
#
#   (B) The cache returns the IDENTICAL matrix, and the 4.1e-14 arises downstream
#       inside estimate_s_rss (eigen + a Brent optimiser with finite tolerance),
#       where a different allocation pattern changes BLAS summation order. That is
#       benign floating-point noise, not a data change.
#
# This decides between them directly, in ONE process, by calling the original and
# the memoized function on the SAME inputs and comparing the returned LD matrix
# with identical(). No statistics involved -- if the matrices are identical, (A)
# is dead and (B) is the only survivor.
#
# The summary stats are synthesised from the block's own .bim so the comparison
# needs no GWAS file: the object under test is the LD matrix path, and both arms
# receive byte-identical inputs either way.

suppressMessages(library(data.table))
suppressMessages(library(Matrix))
suppressMessages(library(dplyr))

FM <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
setwd(FM)

source(file.path(FM, "src/finemapping_functions.R"))
orig_fn <- get_ld_per_locus                       # capture BEFORE the override

source(file.path(FM, "src/perf/get_ld_per_locus_cached_sparse.R"))
cached_fn <- get_ld_per_locus                     # memoized version

CHR <- 21
blocks <- fread(file.path(get_ld_base_dir("EUR"), "approx_LD_blocks.txt"))
b <- blocks[chr == CHR][order(start)][2]          # a mid-size chr21 block
cat(sprintf("block: chr%d %d-%d\n", CHR, b$start, b$stop))

stem <- file.path(get_ld_base_dir("EUR"), paste0("chr", CHR),
                  paste0(b$start, ".", b$stop), paste0(b$start, ".", b$stop))
bim <- fread(paste0(stem, ".bim"))
setnames(bim, c("chr", "rsid", "dk", "pos", "alt", "ref"))
cat(sprintf("block variants: %d\n\n", nrow(bim)))

set.seed(1)
# Several overlapping windows inside the SAME block, so the memoized arm gets a
# cache MISS on the first and HITS on the rest -- which is exactly the code path
# that produced the discrepancy.
n_win <- 6
res <- data.table()
for (w in seq_len(n_win)) {
  lo <- as.integer(quantile(bim$pos, 0.05 * w))
  hi <- as.integer(quantile(bim$pos, 0.95))
  sub <- bim[pos >= lo & pos <= hi]
  if (nrow(sub) < 200) next
  sub <- sub[sample(.N, min(.N, 3000))][order(pos)]

  ss <- data.frame(
    chromosome     = CHR,
    position       = sub$pos,
    allele1        = sub$alt,
    allele2        = sub$ref,
    beta           = rnorm(nrow(sub), 0, 0.02),
    standard_error = rep(0.01, nrow(sub))
  )

  a <- orig_fn  (ss, LOCUS = paste0("w", w), CHR = CHR, START = lo, END = hi, ancestry = "EUR")
  c_ <- cached_fn(ss, LOCUS = paste0("w", w), CHR = CHR, START = lo, END = hi, ancestry = "EUR")

  if (is.null(a) || is.null(c_)) { cat(sprintf("window %d: one arm returned NULL\n", w)); next }

  A <- as.matrix(a[[2]]); C <- as.matrix(c_[[2]])
  same_dim <- identical(dim(A), dim(C))
  bit_id   <- same_dim && identical(A, C)
  maxdiff  <- if (same_dim) max(abs(A - C)) else NA_real_
  ss_id    <- identical(as.data.frame(a[[1]]), as.data.frame(c_[[1]]))

  res <- rbind(res, data.table(window = w, n = nrow(A), same_dim = same_dim,
                               ld_identical = bit_id, max_abs_diff = maxdiff,
                               ss_identical = ss_id))
  cat(sprintf("window %d: n=%5d  dim_match=%s  LD identical=%s  max|diff|=%.3g  ss identical=%s\n",
              w, nrow(A), same_dim, bit_id, maxdiff, ss_id))
}

cat("\n=================== VERDICT ===================\n")
if (nrow(res) == 0) {
  cat("INCONCLUSIVE: no window produced output from both arms\n")
} else if (all(res$ld_identical) && all(res$ss_identical)) {
  cat("(A) REFUTED / (B) SUPPORTED:\n")
  cat("  The memoized function returns a BIT-IDENTICAL LD matrix and summary-stat\n")
  cat("  table on every window tested. The cache does not alter data. The 4.1e-14\n")
  cat("  shift in lambda_s_locus therefore arises downstream, inside estimate_s_rss.\n")
} else {
  cat("(A) SUPPORTED -- THE CACHE ALTERS DATA. DO NOT ADOPT.\n")
  print(res[ld_identical == FALSE | ss_identical == FALSE])
}
cat("===============================================\n")
print(res)
