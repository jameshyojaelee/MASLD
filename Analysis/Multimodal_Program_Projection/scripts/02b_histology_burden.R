#!/usr/bin/env Rscript

# Summarize the four nonredundant pathologist-scored features into a single
# outcome-only histologic-burden axis, then estimate its association with each
# testable frozen-program DIA-MS score in MASLD samples. NAS is deliberately
# excluded from the PCA because it is a composite of steatosis, ballooning,
# and lobular inflammation and would otherwise double-weight those features.
#
# The primary association is the conventional rank-residual partial Spearman:
# rank the program score and burden, residualize both ranks on the nuisance
# design, and correlate the residuals. P values use batch-stratified
# Freedman-Lane permutations. Confidence intervals resample patients within
# acquisition batch and repeat the PCA, rank residualization, and association.
# The resulting intervals therefore represent uncertainty in both the burden
# axis and the program association; they are not the min/max of outcomes.

suppressPackageStartupMessages({
  library(data.table)
})

BASE <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
OUT <- file.path(
  BASE, "Analysis/Multimodal_Program_Projection/results/proteomics"
)
REGISTRY_FILE <- file.path(
  BASE, "Analysis/Multimodal_Program_Projection/results/frozen_programs.tsv"
)
SCORE_FILE <- file.path(OUT, "module_protein_scores.tsv")
META_FILE <- file.path(OUT, "panel4c_metadata.tsv")

inputs <- c(REGISTRY_FILE, SCORE_FILE, META_FILE)
if (any(!file.exists(inputs))) {
  stop("Missing histology-burden input(s): ",
       paste(inputs[!file.exists(inputs)], collapse = ", "))
}

registry <- fread(REGISTRY_FILE)
scores <- fread(SCORE_FILE)
meta <- fread(META_FILE)
components <- c("Steatosis", "Ballooning", "Inflammation", "Fibrosis")
required_meta <- c(
  "sample_id", "group", "acquisition_batch", "sex", "age", "bmi",
  components
)
if (!all(required_meta %in% names(meta))) {
  stop("Metadata lacks required columns: ",
       paste(setdiff(required_meta, names(meta)), collapse = ", "))
}

meta_masld <- meta[
  group == "MASLD" & complete.cases(meta[, ..components])
]
if (nrow(meta_masld) < 20L) {
  stop("Too few complete MASLD samples for histology-burden analysis")
}

derive_burden <- function(d) {
  z <- scale(as.matrix(d[, ..components]))
  if (any(!is.finite(z))) return(NULL)
  fit <- prcomp(z, center = FALSE, scale. = FALSE)
  loading <- fit$rotation[, 1L]
  score <- fit$x[, 1L]
  if (sum(loading) < 0) {
    loading <- -loading
    score <- -score
  }
  list(
    score = as.numeric(score),
    loading = loading,
    variance_explained = summary(fit)$importance[2L, 1L]
  )
}

axis_fit <- derive_burden(meta_masld)
if (is.null(axis_fit) || any(axis_fit$loading <= 0)) {
  stop("Histology PC1 failed the prespecified coherent-positive-loading check")
}
meta_masld[, histology_burden := axis_fit$score]

axis_summary <- data.table(
  feature = components,
  loading = unname(axis_fit$loading[components]),
  variance_explained = axis_fit$variance_explained,
  n_masld = nrow(meta_masld),
  definition = paste(components, collapse = " + "),
  nas_excluded_reason = "NAS is composite of steatosis, ballooning, and inflammation"
)
fwrite(
  axis_summary,
  file.path(OUT, "histology_burden_axis.tsv"),
  sep = "\t",
  quote = FALSE
)

rank_residualize <- function(y, d, include_batch = TRUE) {
  rhs <- c("age", "bmi")
  if (include_batch && uniqueN(d$acquisition_batch[!is.na(d$acquisition_batch)]) > 1L) {
    rhs <- c("acquisition_batch", rhs)
  }
  if (uniqueN(d$sex[!is.na(d$sex)]) > 1L) rhs <- c(rhs, "sex")
  nuisance <- model.matrix(reformulate(rhs), data = d)
  ok <- is.finite(y) & apply(nuisance, 1L, function(z) all(is.finite(z)))
  out <- rep(NA_real_, length(y))
  if (sum(ok) > ncol(nuisance) + 2L) {
    out[ok] <- residuals(lm.fit(
      nuisance[ok, , drop = FALSE],
      rank(y[ok], ties.method = "average")
    ))
  }
  out
}

partial_spearman <- function(d, score_col, include_batch = TRUE) {
  score_resid <- rank_residualize(d[[score_col]], d, include_batch)
  burden_resid <- rank_residualize(d$histology_burden, d, include_batch)
  ok <- is.finite(score_resid) & is.finite(burden_resid)
  if (sum(ok) < 8L) {
    return(c(n = sum(ok), rho = NA_real_, pvalue_asymptotic = NA_real_))
  }
  tst <- suppressWarnings(cor.test(score_resid[ok], burden_resid[ok], method = "pearson"))
  c(n = sum(ok), rho = unname(tst$estimate), pvalue_asymptotic = tst$p.value)
}

freedman_lane_pvalue <- function(d, score_col, n_perm, seed) {
  rhs <- c("acquisition_batch", "age", "bmi")
  if (uniqueN(d$sex[!is.na(d$sex)]) > 1L) rhs <- c(rhs, "sex")
  nuisance <- model.matrix(reformulate(rhs), data = d)
  ok <- is.finite(d[[score_col]]) & is.finite(d$histology_burden) &
    apply(nuisance, 1L, function(z) all(is.finite(z)))
  if (sum(ok) < 8L) return(NA_real_)

  nuisance <- nuisance[ok, , drop = FALSE]
  x_rank <- rank(d[[score_col]][ok], ties.method = "average")
  y_rank <- rank(d$histology_burden[ok], ties.method = "average")
  batch <- droplevels(factor(d$acquisition_batch[ok]))
  nuisance_qr <- qr(nuisance)
  x_resid <- qr.resid(nuisance_qr, x_rank)
  y_resid <- qr.resid(nuisance_qr, y_rank)
  observed <- cor(x_resid, y_resid)
  fitted_reduced <- y_rank - y_resid
  by_batch <- split(seq_along(y_resid), batch)

  set.seed(seed)
  exceed <- 0L
  for (i in seq_len(n_perm)) {
    perm_index <- seq_along(y_resid)
    for (ii in by_batch) perm_index[ii] <- sample(ii, length(ii), replace = FALSE)
    y_star <- fitted_reduced + y_resid[perm_index]
    perm_stat <- cor(x_resid, qr.resid(nuisance_qr, y_star))
    if (is.finite(perm_stat) && abs(perm_stat) >= abs(observed)) exceed <- exceed + 1L
  }
  (exceed + 1) / (n_perm + 1)
}

bootstrap_rho <- function(d, score_col) {
  index_by_batch <- split(seq_len(nrow(d)), d$acquisition_batch)
  idx <- unlist(lapply(
    index_by_batch,
    function(ii) sample(ii, length(ii), replace = TRUE)
  ), use.names = FALSE)
  b <- copy(d[idx])
  burden <- derive_burden(b)
  if (is.null(burden) || any(burden$loading <= 0)) return(NA_real_)
  b[, histology_burden := burden$score]
  unname(partial_spearman(b, score_col)["rho"])
}

analysis_data <- merge(
  scores[group == "MASLD"],
  meta_masld[, .(
    sample_id, acquisition_batch, sex, age, bmi,
    Steatosis, Ballooning, Inflammation, Fibrosis, histology_burden
  )],
  by = c("sample_id", "acquisition_batch"),
  all = FALSE,
  sort = FALSE
)

n_boot <- as.integer(Sys.getenv("FIG4_HIST_BOOTSTRAPS", "4999"))
if (!is.finite(n_boot) || n_boot < 999L) {
  stop("FIG4_HIST_BOOTSTRAPS must be an integer >= 999")
}
n_perm <- as.integer(Sys.getenv("FIG4_HIST_PERMUTATIONS", "9999"))
if (!is.finite(n_perm) || n_perm < 999L) {
  stop("FIG4_HIST_PERMUTATIONS must be an integer >= 999")
}

program_ids <- unique(analysis_data$program_id)
rows <- rbindlist(lapply(seq_along(program_ids), function(program_index) {
  pid <- program_ids[[program_index]]
  d <- analysis_data[program_id == pid]
  primary <- partial_spearman(d, "score")
  equal <- partial_spearman(d, "equal_score")
  leave_top <- partial_spearman(d, "leave_top_score")
  batch_rho <- vapply(c("2019", "2020"), function(batch_id) {
    unname(partial_spearman(
      d[acquisition_batch == batch_id], "score", include_batch = FALSE
    )["rho"])
  }, numeric(1))
  permutation_p <- freedman_lane_pvalue(
    d, "score", n_perm, seed = 20260806L + program_index
  )
  set.seed(30260806L + program_index)
  boot <- replicate(n_boot, bootstrap_rho(d, "score"))
  ci <- quantile(boot, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
  data.table(
    program_id = pid,
    n = as.integer(primary["n"]),
    rho = unname(primary["rho"]),
    pvalue = permutation_p,
    pvalue_asymptotic = unname(primary["pvalue_asymptotic"]),
    ci_low = ci[1L],
    ci_high = ci[2L],
    n_boot = n_boot,
    n_boot_finite = sum(is.finite(boot)),
    n_permutations = n_perm,
    equal_rho = unname(equal["rho"]),
    leave_top_rho = unname(leave_top["rho"]),
    rho_2019 = batch_rho[1L],
    rho_2020 = batch_rho[2L]
  )
}))
rows[, padj := p.adjust(pvalue, method = "BH")]
rows[, sensitivity_sign_agree := is.finite(rho) & is.finite(equal_rho) &
  is.finite(leave_top_rho) & sign(rho) == sign(equal_rho) &
  sign(rho) == sign(leave_top_rho)]
rows[, batch_sign_agree := is.finite(rho) & is.finite(rho_2019) &
  is.finite(rho_2020) & sign(rho) == sign(rho_2019) &
  sign(rho) == sign(rho_2020)]
rows[, robust := is.finite(padj) & padj < 0.05 &
  sensitivity_sign_agree & batch_sign_agree]

out <- merge(
  registry[, .(
    program_id, display_order, cell_type, module, program_name
  )],
  rows,
  by = "program_id",
  all.x = TRUE,
  sort = FALSE
)
out[, testable := is.finite(rho)]
out[, inference_method := "rank-residual partial Spearman; batch-stratified Freedman-Lane permutation"]
setorder(out, display_order)
fwrite(
  out,
  file.path(OUT, "module_protein_histology_burden_masld.tsv"),
  sep = "\t",
  quote = FALSE
)

if (any(rows$n_boot_finite < 0.95 * n_boot)) {
  stop("Fewer than 95% finite bootstrap estimates for at least one program")
}
message(
  "[histology burden] MASLD n=", nrow(meta_masld),
  "; PC1 variance=", sprintf("%.3f", axis_fit$variance_explained),
  "; testable programs=", nrow(rows),
  "; robust programs=", sum(rows$robust),
  "; bootstraps/program=", n_boot,
  "; permutations/program=", n_perm
)
