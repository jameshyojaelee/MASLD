# Exact-stage null for distinct source/genetic units. No stage substitutions.
n4_null_capacity <- function(d, pools) {
  demand <- d[, .(required = .N), by = .(cohort, fibrosis_stage)]
  capacity <- data.table::rbindlist(lapply(names(pools), function(co) {
    pools[[co]][, .(available = .N), by = fibrosis_stage][, cohort := co]
  }), fill = TRUE)
  ans <- merge(demand, capacity, by = c("cohort", "fibrosis_stage"), all.x = TRUE)
  ans[is.na(available), available := 0L]
  ans[, sufficient := available >= required]
  ans
}
n4_draw_null <- function(d, pools) {
  capacity <- n4_null_capacity(d, pools)
  if (!all(capacity$sufficient)) stop("insufficient exact-stage eligible null units")
  picks <- character(nrow(d))
  for (co in unique(d$cohort)) for (stage in unique(d[cohort == co, fibrosis_stage])) {
    ii <- which(d$cohort == co & d$fibrosis_stage == stage)
    candidates <- pools[[co]][fibrosis_stage == stage, sample_id]
    picks[ii] <- candidates[sample.int(length(candidates), length(ii), replace = FALSE)]
  }
  stopifnot(!anyDuplicated(picks))
  list(sample_id = picks, substitutions = 0L)
}
