# U8a independent re-derivation library (code review 2026-10-23).
# Written from METHODS / PRESPEC_H59 text, not from the H59 producer code.
suppressPackageStartupMessages(library(data.table))

ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
SCORE_FILE <- file.path(ROOT, "Analysis/Multimodal_Program_Projection/candidates",
  "program-context-v2-stage-corrected-candidate-2026-08-13-v3/hotspot/donor_program_scores_primary.tsv")
META_FILE <- file.path(ROOT, "Analysis/SingleCell/candidates/source-diagnosis-primary-correction-20261014",
  "metadata/donor_metadata_stage_corrected_donor.tsv")
REGISTRY_FILE <- file.path(ROOT, "figures/main/fig4_singlecell_programs/source_tables",
  "source_diagnosis_correction_20261014/program_registry_v2.tsv")
H59 <- file.path(ROOT, "Analysis/Hypotheses/H59-20261015T112432Z")
FIT_DS <- c("GSE202379", "GSE244832")
STAGES <- c(Healthy = 0, Steatosis = 1, Steatohepatitis = 2)

say <- function(...) cat(sprintf(...), "\n", sep = "")

# Donor table: score file joined to the 2026-10-14 source-diagnosis metadata.
# Stage, dataset and exclusion come only from the metadata.
load_donor_scores <- function(score_file = SCORE_FILE) {
  s <- fread(score_file)
  if (!"dataset" %in% names(s)) s[, dataset := NA_character_]
  s <- s[, .(cell_type, module, donor, score, dataset_scorefile = dataset)]
  m <- fread(META_FILE)
  m <- m[, .(donor = sample, dataset, stage_label = disease_stage_coarse,
             exclude = as.logical(exclude_stage_analysis))]
  say("score rows %d, score donors %d, metadata donors %d", nrow(s), uniqueN(s$donor), nrow(m))
  j <- merge(s, m, by = "donor", all.x = TRUE)
  say("after left join: %d rows; rows without metadata %d", nrow(j), sum(is.na(j$dataset)))
  stopifnot(!anyNA(j$dataset), all(is.na(j$dataset_scorefile) | j$dataset == j$dataset_scorefile))
  j[, stage := unname(STAGES[stage_label])]
  u62 <- j[exclude == FALSE & !is.na(stage)]
  say("62-universe rows %d donors %d", nrow(u62), uniqueN(u62$donor))
  u62[, program := paste0(cell_type, "_", module)]
  u62[]
}

# Profiled REML for score ~ stage + dataset with a separate residual variance per
# dataset (two datasets).  Closed form given the variance ratio; 1-D search on it.
gls_fit <- function(y, s, d) {
  lv <- sort(unique(d))
  stopifnot(length(lv) == 2L)
  st <- lapply(lv, function(l) {
    i <- d == l
    yc <- y[i] - mean(y[i]); sc <- s[i] - mean(s[i])
    c(n = sum(i), sxx = sum(sc^2), sxy = sum(sc * yc), syy = sum(yc^2))
  })
  st <- do.call(rbind, st)
  n <- sum(st[, "n"]); p <- 3
  core <- function(rho) {             # rho = log(sd_2 / sd_1)
    w <- c(1, exp(-2 * rho))
    sw <- sum(w * st[, "sxx"])
    b <- sum(w * st[, "sxy"]) / sw
    rss <- sum(w * (st[, "syy"] - 2 * b * st[, "sxy"] + b^2 * st[, "sxx"]))
    list(b = b, sw = sw, rss = rss, w = w)
  }
  negll <- function(rho) {
    k <- core(rho)
    0.5 * (n - p) * log(k$rss / (n - p)) + st[2, "n"] * rho - rho + 0.5 * log(k$sw)
  }
  o <- optimize(negll, c(-12, 12), tol = 1e-11)
  k <- core(o$minimum)
  s2 <- k$rss / (n - p)
  se <- sqrt(s2 / k$sw)
  list(beta = k$b, se = se, t = k$b / se, sigma1 = sqrt(s2), sigma2 = sqrt(s2) * exp(o$minimum),
       rho = o$minimum, n = n, levels = lv, negll = o$objective)
}

# Equal-variance OLS with dataset intercepts; HC3 SE.  Returns NA HC3 when a leverage is 1.
ols_fit <- function(y, s, d) {
  d <- factor(d)
  X <- if (nlevels(d) > 1) cbind(model.matrix(~ d), stage = s) else cbind(1, stage = s)
  q <- qr(X)
  if (q$rank < ncol(X)) return(list(beta = NA_real_, se = NA_real_, t = NA_real_, df = NA_integer_,
                                   se_hc3 = NA_real_, maxlev = NA_real_))
  XtXi <- chol2inv(qr.R(q))
  b <- drop(XtXi %*% crossprod(X, y))
  r <- y - drop(X %*% b)
  df <- nrow(X) - ncol(X)
  s2 <- sum(r^2) / df
  j <- ncol(X)
  h <- rowSums((X %*% XtXi) * X)
  hc3 <- if (any(h > 1 - 1e-10)) NA_real_ else {
    m <- XtXi %*% crossprod(X * (r / (1 - h))) %*% XtXi
    sqrt(m[j, j])
  }
  list(beta = b[j], se = sqrt(s2 * XtXi[j, j]), t = b[j] / sqrt(s2 * XtXi[j, j]), df = df,
       se_hc3 = hc3, maxlev = max(h))
}

# Within-dataset permutation per PRESPEC_H59: the donor with the j-th smallest key
# receives the j-th stage of the dataset's donors in sorted donor-ID order.
permute_stage <- function(donor, stage, dataset, keyrow) {
  out <- stage
  for (l in unique(dataset)) {
    i <- which(dataset == l)
    i <- i[order(donor[i], method = "radix")]
    s_sorted <- stage[i]
    k <- keyrow[donor[i]]
    out[i[order(k)]] <- s_sorted
  }
  out
}

# The prespec draws a 19,999 x 62 matrix (column-major); smaller B uses its first rows.
make_keys <- function(seed, B, donors62, B_full = 19999L) {
  set.seed(seed)
  K <- matrix(runif(B_full * length(donors62)), nrow = B_full)
  colnames(K) <- sort(donors62, method = "radix")
  K[seq_len(B), , drop = FALSE]
}

perm_p <- function(tobs, tnull) {
  ok <- is.finite(tnull)
  (1 + sum(abs(tnull[ok]) >= abs(tobs))) / (1 + sum(ok))
}
