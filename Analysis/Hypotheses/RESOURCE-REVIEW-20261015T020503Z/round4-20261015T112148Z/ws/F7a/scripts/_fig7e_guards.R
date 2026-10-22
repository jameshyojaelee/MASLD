# Guards for the 7E verdict column. Both can fail, and both DO fail on the pre-fix
# table (see 99_demo_verdict_guard_fires.R).
#
# G1 catches a class string that is not a member of the figure's vocabulary - which is
#    what a conditional field looks like before anyone shortens it.
# G2 does not look at any field name at all. It asserts the invariant the vocabulary is
#    supposed to encode: a row is 'supported' exactly when its observed statistic
#    reached its own threshold. G2 still fires if someone launders a conditional label
#    down to a bare vocabulary word, which G1 alone would let through.

# The vocabulary is read from the framework table, not typed here, so the two panels
# cannot drift apart.
verdict_vocabulary <- function(framework_tab) {
  d <- framework_tab$detail[framework_tab$stage == "6. verdict vocabulary" &
                            framework_tab$item  == "states a question can receive"]
  stopifnot(length(d) == 1L)
  v <- trimws(strsplit(d, "/", fixed = TRUE)[[1]])
  stopifnot(length(v) == 4L, all(nzchar(v)))
  v
}

check_verdict_rows <- function(v, vocab) {
  stopifnot(all(c("verdict", "question", "observed", "threshold") %in% names(v)))
  obs <- suppressWarnings(as.numeric(v$observed))
  thr <- suppressWarnings(as.numeric(v$threshold))
  if (any(is.na(obs)) || any(is.na(thr)) || any(thr <= 0))
    stop("G0 FAILED: observed/threshold missing or threshold not positive")

  bad <- setdiff(unique(v$verdict), vocab)
  if (length(bad))
    stop("G1 FAILED: verdict value(s) outside the four-word vocabulary, which is what a ",
         "conditional 'if it had been negative' field looks like: ",
         paste(sprintf("%s", bad), collapse = " | "))

  reached <- obs >= thr
  claims  <- v$verdict == "supported"
  wrong   <- which(claims != reached)
  if (length(wrong))
    stop("G2 FAILED: verdict disagrees with the realised comparison of observed against ",
         "its own threshold for ", length(wrong), " row(s): ",
         paste(sprintf("'%s' labelled %s but observed %.4f vs threshold %.4f (ratio %.3f)",
                       v$question[wrong], v$verdict[wrong], obs[wrong], thr[wrong],
                       obs[wrong] / thr[wrong]), collapse = "; "))

  message("[GUARD] G1 pass: all ", nrow(v), " verdicts are vocabulary members.")
  message("[GUARD] G2 pass: 'supported' agrees with observed >= threshold on all ",
          nrow(v), " rows (", sum(claims), " supported, ", sum(!claims), " not).")
  invisible(TRUE)
}
