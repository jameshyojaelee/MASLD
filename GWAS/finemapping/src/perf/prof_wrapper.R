# prof_wrapper.R -- statistical profile of a REAL 06_susie_coloc.R run.
#
# Answers empirically: how much of the ~88 s/gene is (a) reading and reshaping
# the LD matrix, versus (b) the actual statistics (susie_rss / estimate_s_rss /
# coloc). Everything so far on that split has been arithmetic on work-counts,
# not measurement.
#
# Method: Rprof samples the R call stack every 20 ms and attributes time to the
# function on top. Time spent inside compiled code (fread's C parser, LAPACK's
# eigen, susieR's internals) is attributed to the R function that called it,
# which is exactly the granularity the question needs.
#
# The run is capped by elapsed setTimeLimit and the partial profile summarised in
# the same session -- Rprof flushes continuously, so a capped run yields a valid
# sample. Nothing is written to any results directory that matters: the caller
# sets COLOC_OUT_SUFFIX=_profile and the run is killed long before the write step.

prof_file <- Sys.getenv("PROF_OUT", unset = "/tmp/coloc_prof.out")
limit_s   <- as.numeric(Sys.getenv("PROF_SECONDS", unset = "1200"))
base      <- Sys.getenv("MASLD_PROJECT_ROOT",
                        unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
script    <- file.path(base, "GWAS/finemapping/src/06_susie_coloc.R")

cat("profiling:", script, "\n")
cat("cap:", limit_s, "s | interval: 0.02 s\n\n")

Rprof(prof_file, interval = 0.02, line.profiling = FALSE)
setTimeLimit(elapsed = limit_s, transient = FALSE)
res <- try(source(script, echo = FALSE), silent = TRUE)
setTimeLimit(elapsed = Inf)
Rprof(NULL)

cat("\n\n=================== PROFILE ===================\n")
if (inherits(res, "try-error")) cat("(run capped as designed:", trimws(sub(".*:", "", res[1])), ")\n\n")

s <- try(summaryRprof(prof_file), silent = TRUE)
if (inherits(s, "try-error")) {
  cat("summaryRprof failed:", as.character(s), "\n")
} else {
  tot <- s$sampling.time
  cat(sprintf("total sampled: %.1f s\n\n", tot))

  cat("--- BY SELF (where the CPU actually sits) ---\n")
  bs <- head(s$by.self, 25)
  bs$pct <- round(100 * bs$self.time / tot, 1)
  print(bs[, c("self.time", "pct", "total.time")])

  cat("\n--- BY TOTAL (inclusive; find the stage owners) ---\n")
  bt <- s$by.total
  keep <- rownames(bt) %in% c("get_ld_per_locus", "fread", "bdiag", "as.matrix",
                              "susie_rss", "estimate_s_rss", "eigen", "coloc.susie",
                              "coloc.abf", "runsusie", "find_LD_block", "rbind",
                              "duplicated", "match", "readRDS", "%in%")
  bt2 <- bt[keep, , drop = FALSE]
  bt2$pct_of_total <- round(100 * bt2$total.time / tot, 1)
  print(bt2[order(-bt2$total.time), c("total.time", "pct_of_total", "self.time")])

  cat("\n--- HEADLINE SPLIT ---\n")
  gt <- function(fn) if (fn %in% rownames(bt)) bt[fn, "total.time"] else 0
  ld_stage   <- gt("get_ld_per_locus")
  stats_est  <- gt("estimate_s_rss")
  stats_sus  <- gt("susie_rss")
  stats_col  <- gt("coloc.susie") + gt("coloc.abf")
  cat(sprintf("  LD load + reshape (get_ld_per_locus) : %7.1f s  %5.1f%%\n", ld_stage, 100*ld_stage/tot))
  cat(sprintf("  estimate_s_rss (diagnostic, O(M^3))  : %7.1f s  %5.1f%%\n", stats_est, 100*stats_est/tot))
  cat(sprintf("  susie_rss (the fine-mapping)         : %7.1f s  %5.1f%%\n", stats_sus, 100*stats_sus/tot))
  cat(sprintf("  coloc.susie + coloc.abf              : %7.1f s  %5.1f%%\n", stats_col, 100*stats_col/tot))
  cat(sprintf("  --- statistics total                 : %7.1f s  %5.1f%%\n",
              stats_est+stats_sus+stats_col, 100*(stats_est+stats_sus+stats_col)/tot))
}
cat("===============================================\n")
