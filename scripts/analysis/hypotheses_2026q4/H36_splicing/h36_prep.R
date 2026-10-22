# H36: shared loader. Reads the output of h36_01 and builds per-cluster count
# matrices exactly as the stage test uses them (h36_02 validates this code path).
# Defines globals: long, cl_tab, pers, idx, stage0, ds0, sex0.

load_build <- function(BUILD) {
  long <<- readRDS(file.path(BUILD, "long.rds"))
  cl_tab <<- fread(file.path(BUILD, "clusters.tsv.gz"))
  setorder(cl_tab, cl)
  stopifnot(identical(cl_tab$cl, seq_len(nrow(cl_tab))))
  pers <<- fread(file.path(BUILD, "persons.tsv"))
  idx <<- split(seq_len(nrow(long)), long$cl)
  stage0 <<- pers$fibrosis_stage
  ds0 <<- pers$dataset
  sex0 <<- as.integer(pers$inferred_sex == "M")
  stopifnot(!anyNA(ds0), !anyNA(sex0))
  invisible(NULL)
}

# Counts of cluster k: persons x junctions (column c is junction rank J + 1 - c).
get_Y <- function(k) {
  ii <- idx[[k]]; J <- cl_tab$J[k]
  up <- sort(unique(long$pi[ii]))
  Y <- matrix(0L, length(up), J)
  Y[cbind(match(long$pi[ii], up), J + 1L - long$jr[ii])] <- long$n[ii]
  list(Y = Y, pi = up, J = J)
}

# Persons used in the test: staged persons in cohorts with >= 10 included staged
# persons (need_stage = TRUE), or all included persons. Zero-count junctions are
# dropped and columns are ordered by total reads so the reference is the most-read
# junction. jr gives the junction rank of each returned column.
prep <- function(k, need_stage = TRUE) {
  g <- get_Y(k); pi <- g$pi
  ok <- rep(TRUE, length(pi))
  if (need_stage) {
    ok <- !is.na(stage0[pi])
    tab <- table(ds0[pi][ok]); keepc <- names(tab)[tab >= 10L]
    ok <- ok & ds0[pi] %in% keepc
  }
  Y <- g$Y[ok, , drop = FALSE]; pi <- pi[ok]
  jr_col <- g$J + 1L - seq_len(g$J)
  nz <- colSums(Y) > 0L
  Y <- Y[, nz, drop = FALSE]; jr_col <- jr_col[nz]
  o <- order(colSums(Y)); Y <- Y[, o, drop = FALSE]; jr_col <- jr_col[o]
  list(Y = Y, pi = pi, jr = jr_col)
}

# Covariates of the null model: cohort (when more than one) and sex.
design0 <- function(pi) {
  d <- factor(ds0[pi])
  X <- if (nlevels(d) > 1L) model.matrix(~ d + sx, data.frame(d = d, sx = sex0[pi])) else
    model.matrix(~ sx, data.frame(sx = sex0[pi]))
  independent_cols(X)
}
