suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

ml_fail <- function(...) stop(paste0(...), call. = FALSE)

ml_assert <- function(condition, message) {
  if (length(condition) != 1L || is.na(condition) || !condition) ml_fail(message)
  invisible(TRUE)
}

ml_project_root <- function() {
  value <- Sys.getenv(
    "MASLD_PROJECT_ROOT",
    unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
  )
  normalizePath(value, mustWork = TRUE)
}

ml_script_dir <- function() {
  file.path(
    ml_project_root(),
    "scripts/analysis/histology_anchored_continuum/molecular_layers"
  )
}

ml_out_root <- function(must_exist = TRUE) {
  value <- Sys.getenv("HAC_ML_OUT_ROOT", unset = "")
  ml_assert(nzchar(value), "HAC_ML_OUT_ROOT is unset")
  normalizePath(value, mustWork = must_exist)
}

ml_read_contract <- function() {
  path <- file.path(ml_script_dir(), "00_contract.json")
  ml_assert(file.exists(path), paste0("Missing molecular-layer contract: ", path))
  jsonlite::fromJSON(path, simplifyVector = TRUE)
}

ml_resolve <- function(path, must_exist = TRUE) {
  if (!grepl("^/", path)) path <- file.path(ml_project_root(), path)
  if (must_exist) ml_assert(file.exists(path), paste0("Required input missing: ", path))
  normalizePath(path, mustWork = must_exist)
}

ml_source_root <- function(contract = ml_read_contract()) {
  override <- Sys.getenv("HAC_SOURCE_ROOT", unset = "")
  ml_resolve(if (nzchar(override)) override else contract$source_hac_candidate)
}

ml_ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  ml_assert(dir.exists(path), paste0("Could not create directory: ", path))
  invisible(path)
}

ml_write_tsv_once <- function(x, path) {
  ml_assert(!file.exists(path), paste0("Refusing to overwrite: ", path))
  ml_ensure_dir(dirname(path))
  data.table::fwrite(x, path, sep = "\t", quote = FALSE, na = "NA")
  invisible(path)
}

ml_write_json_once <- function(x, path) {
  ml_assert(!file.exists(path), paste0("Refusing to overwrite: ", path))
  ml_ensure_dir(dirname(path))
  jsonlite::write_json(
    x, path, auto_unbox = TRUE, pretty = TRUE, null = "null", na = "string"
  )
  invisible(path)
}

ml_sha256 <- function(path) {
  ml_assert(file.exists(path), paste0("Cannot hash missing file: ", path))
  output <- system2("sha256sum", path, stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  ml_assert(is.null(status) || identical(status, 0L), paste0("Hash failed: ", path))
  sub("[[:space:]].*$", "", output[[1L]])
}

ml_base_gene_id <- function(x) sub("\\.[0-9]+$", "", as.character(x))

ml_standardize <- function(x) {
  answer <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  if (sum(ok) >= 2L && stats::sd(x[ok]) > 0) answer[ok] <- as.numeric(scale(x[ok]))
  answer
}

ml_percentile <- function(x) {
  answer <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  if (sum(ok) == 1L) answer[ok] <- 0.5
  if (sum(ok) > 1L) {
    answer[ok] <- (rank(x[ok], ties.method = "average") - 1) / (sum(ok) - 1)
  }
  answer
}

ml_complete_bh <- function(p, family_size) {
  answer <- rep(NA_real_, length(p))
  ok <- is.finite(p)
  if (any(ok)) answer[ok] <- p.adjust(p[ok], method = "BH", n = family_size)
  answer
}

ml_fixed_meta <- function(beta, se) {
  ok <- is.finite(beta) & is.finite(se) & se > 0
  if (sum(ok) < 2L) {
    return(data.table(
      estimable = FALSE, beta = NA_real_, se = NA_real_, z = NA_real_,
      p_value = NA_real_, ci_low = NA_real_, ci_high = NA_real_, n_cohorts = sum(ok)
    ))
  }
  weight <- 1 / se[ok]^2
  estimate <- sum(weight * beta[ok]) / sum(weight)
  estimate_se <- sqrt(1 / sum(weight))
  z <- estimate / estimate_se
  data.table(
    estimable = TRUE,
    beta = estimate,
    se = estimate_se,
    z = z,
    p_value = 2 * pnorm(abs(z), lower.tail = FALSE),
    ci_low = estimate - qnorm(0.975) * estimate_se,
    ci_high = estimate + qnorm(0.975) * estimate_se,
    n_cohorts = sum(ok)
  )
}

ml_random_meta <- function(beta, se) {
  ok <- is.finite(beta) & is.finite(se) & se > 0
  beta <- beta[ok]
  se <- se[ok]
  k <- length(beta)
  if (k < 2L) {
    return(data.table(
      estimable = FALSE, beta = NA_real_, se = NA_real_, p_value = NA_real_,
      ci_low = NA_real_, ci_high = NA_real_, tau2 = NA_real_, n_cohorts = k
    ))
  }
  fixed_weight <- 1 / se^2
  fixed_beta <- sum(fixed_weight * beta) / sum(fixed_weight)
  q <- sum(fixed_weight * (beta - fixed_beta)^2)
  c_value <- sum(fixed_weight) - sum(fixed_weight^2) / sum(fixed_weight)
  tau2 <- max(0, (q - (k - 1)) / c_value)
  weight <- 1 / (se^2 + tau2)
  estimate <- sum(weight * beta) / sum(weight)
  estimate_se <- sqrt(1 / sum(weight))
  z <- estimate / estimate_se
  data.table(
    estimable = TRUE, beta = estimate, se = estimate_se,
    p_value = 2 * pnorm(abs(z), lower.tail = FALSE),
    ci_low = estimate - qnorm(0.975) * estimate_se,
    ci_high = estimate + qnorm(0.975) * estimate_se,
    tau2 = tau2, n_cohorts = k
  )
}

ml_load_axes <- function(contract = ml_read_contract()) {
  source_root <- ml_source_root(contract)
  unsupervised <- fread(file.path(source_root, "unsupervised", "participant_scores.tsv"))
  projection <- fread(file.path(source_root, "projection", "participant_scores.tsv"))
  axes <- rbindlist(list(unsupervised, projection), use.names = TRUE, fill = TRUE)
  axes <- unique(
    axes[axis_id %in% contract$co_primary_axes,
         .(sample_id, dataset, axis_id, axis_raw, axis_percentile)],
    by = c("sample_id", "axis_id")
  )
  ml_assert(
    all(contract$co_primary_axes %in% axes$axis_id),
    "A frozen co-primary continuum score is missing"
  )
  axes
}

ml_fit_outcome <- function(data, outcome = "outcome_z", stage = "fibrosis_stage") {
  d <- copy(data)
  required <- c(outcome, "axis_raw", stage, "inferred_sex")
  ml_assert(all(required %in% names(d)), "Outcome model columns are incomplete")
  d <- d[
    is.finite(get(outcome)) & is.finite(axis_raw) &
      !is.na(get(stage)) & !is.na(inferred_sex)
  ]
  if (nrow(d) < 20L || uniqueN(d[[stage]]) < 2L || uniqueN(d$inferred_sex) < 2L) {
    return(data.table(
      estimable = FALSE, n = nrow(d), beta = NA_real_, se = NA_real_,
      p_value = NA_real_, ci_low = NA_real_, ci_high = NA_real_
    ))
  }
  d[, axis_z := ml_standardize(axis_raw)]
  d[, stage_factor := factor(get(stage))]
  d[, sex_factor := factor(inferred_sex)]
  formula <- as.formula(paste(outcome, "~ axis_z + stage_factor + sex_factor"))
  fit <- tryCatch(lm(formula, data = d), error = function(e) NULL)
  if (is.null(fit) || !"axis_z" %in% rownames(coef(summary(fit)))) {
    return(data.table(
      estimable = FALSE, n = nrow(d), beta = NA_real_, se = NA_real_,
      p_value = NA_real_, ci_low = NA_real_, ci_high = NA_real_
    ))
  }
  coefficient <- coef(summary(fit))["axis_z", ]
  data.table(
    estimable = TRUE, n = nrow(d), beta = coefficient[["Estimate"]],
    se = coefficient[["Std. Error"]], p_value = coefficient[["Pr(>|t|)"]],
    ci_low = coefficient[["Estimate"]] - qt(0.975, df.residual(fit)) * coefficient[["Std. Error"]],
    ci_high = coefficient[["Estimate"]] + qt(0.975, df.residual(fit)) * coefficient[["Std. Error"]]
  )
}

ml_signature_symbols <- function(contract = ml_read_contract()) {
  path <- file.path(ml_source_root(contract), "programs", "signature_gencode_v49_mapping.tsv")
  mapping <- fread(path)
  column <- if ("gencode_v49_gene_symbol" %in% names(mapping)) {
    "gencode_v49_gene_symbol"
  } else {
    "released_gene_symbol"
  }
  unique(toupper(mapping[[column]]))
}

ml_write_session_info <- function(path) {
  ml_assert(!file.exists(path), paste0("Refusing to overwrite: ", path))
  ml_ensure_dir(dirname(path))
  writeLines(capture.output(sessionInfo()), path)
  invisible(path)
}

ml_fixed_windows <- function(percentile, centers, width) {
  rbindlist(lapply(seq_along(centers), function(index) {
    lower <- centers[[index]] - width / 2
    upper <- centers[[index]] + width / 2
    included <- if (index == length(centers)) {
      is.finite(percentile) & percentile >= lower & percentile <= upper
    } else {
      is.finite(percentile) & percentile >= lower & percentile < upper
    }
    # An empty window would recycle the scalar columns against integer(0) and
    # vanish as a 0-row table, making "no participants in window" indistinguishable
    # from "window never computed" downstream.
    if (!any(included)) {
      stop(sprintf("Window %d (center %.4f, [%.4f, %.4f)) contains no participants",
                   index, centers[[index]], lower, upper))
    }
    data.table(
      window_id = index, center = centers[[index]], lower = lower,
      upper = upper, sample_index = which(included)
    )
  }))
}
