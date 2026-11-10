#!/usr/bin/env Rscript
# Candidate donor-level accessibility comparisons; named results are read-only.
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(edgeR)
  library(limma)
  library(jsonlite)
})

`%||%` <- function(x, y) if (is.null(x)) y else x
require_true <- function(x, message) if (!isTRUE(x)) stop(message, call. = FALSE)

parallel_da <- function(items, fun, workers = 1L) {
  if (workers == 1L) return(lapply(items, fun))
  require_true(.Platform$OS.type == "unix", "Parallel resampling requires Unix fork support")
  values <- parallel::mclapply(items, function(item) {
    Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
      BLIS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
    data.table::setDTthreads(1L)
    if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
      RhpcBLASctl::blas_set_num_threads(1L)
      RhpcBLASctl::omp_set_num_threads(1L)
    }
    fun(item)
  }, mc.cores = workers, mc.preschedule = TRUE, mc.set.seed = FALSE, mc.silent = TRUE)
  require_true(length(values) == length(items) && !any(vapply(values, is.null, logical(1))),
    "Parallel resampling lost a child result; inspect allocation memory/worker errors")
  failed <- vapply(values, inherits, logical(1), "try-error")
  require_true(!any(failed), paste("Parallel resampling failed:", paste(values[failed], collapse = "; ")))
  values
}

read_sparse <- function(path) {
  con <- gzfile(path, "rb")
  on.exit(close(con))
  as(readMM(con), "CsparseMatrix")
}

load_da_counts <- function(root, cohort, lineage) {
  prefix <- file.path(root, cohort, lineage)
  paths <- paste0(prefix, c(".mtx.gz", ".donors.tsv", ".peaks.tsv"))
  require_true(all(file.exists(paths)), paste("Missing count input:", prefix))
  counts <- read_sparse(paths[1L])
  donors <- fread(paths[2L])
  peaks <- fread(paths[3L])
  require_true(all(c("donor_id", "condition", "n_cells", "contrast_eligible") %in% names(donors)), "Missing donor eligibility columns")
  require_true("peak_coordinate" %in% names(peaks), "Missing peak_coordinate")
  require_true(nrow(counts) == nrow(donors) && ncol(counts) == nrow(peaks), "Count dimensions differ")
  require_true(!anyDuplicated(donors$donor_id) && !anyDuplicated(peaks$peak_coordinate), "Duplicated donor or peak")
  require_true(all(is.finite(counts@x) & counts@x >= 0 & counts@x == floor(counts@x)), "Counts must be nonnegative integers")
  donors[, condition := toupper(condition)]
  keep <- toupper(as.character(donors$contrast_eligible)) %in% c("TRUE", "1") & donors$n_cells >= 20L & donors$condition %chin% c("NORMAL", "MASH")
  donors <- donors[keep]
  counts <- counts[keep, , drop = FALSE]
  rownames(counts) <- donors$donor_id
  colnames(counts) <- peaks$peak_coordinate
  list(counts = counts, donors = donors, peaks = peaks)
}

filter_da_peaks <- function(counts, rule) {
  if (rule == "primary") {
    y <- DGEList(counts = t(counts))
    return(filterByExpr(y, design = matrix(1, ncol(y), 1L)))
  }
  require_true(rule == "prevalence", "Unknown filter")
  lib <- Matrix::rowSums(counts)
  require_true(all(lib > 0), "Zero donor library")
  raw_cpm <- Diagonal(x = 1e6 / lib) %*% counts
  Matrix::colSums(raw_cpm >= 1) >= 4L & Matrix::colSums(counts) >= 20
}

make_da_design <- function(donors, covariates = character()) {
  require_true(all(donors$condition %chin% c("NORMAL", "MASH")), "Invalid primary condition")
  d <- as.data.frame(donors)
  d$group <- factor(d$condition, levels = c("NORMAL", "MASH"))
  require_true(all(covariates %in% names(d)), "Missing covariate")
  require_true(!anyNA(d[, c("group", covariates), drop = FALSE]), "Missing covariate values")
  for (column in covariates) {
    if (is.character(d[[column]]) || grepl("batch|sex", column, ignore.case = TRUE)) d[[column]] <- factor(d[[column]])
    require_true(length(unique(d[[column]])) > 1L, paste("Constant covariate:", column))
  }
  design <- model.matrix(reformulate(c(covariates, "group")), d)
  require_true(qr(design)$rank == ncol(design), "DA design is not full rank")
  require_true(nrow(design) - ncol(design) >= 2L, "Fewer than two residual degrees of freedom")
  coefficient <- match("groupMASH", colnames(design))
  require_true(!is.na(coefficient), "MASH coefficient absent")
  list(design = design, coefficient = coefficient)
}

fit_da <- function(counts, donors, method, covariates = character()) {
  specification <- make_da_design(donors, covariates)
  y <- calcNormFactors(DGEList(counts = t(counts)), method = "TMM")
  require_true(all(y$samples$lib.size > 0), "Zero library after filtering")
  design <- specification$design
  coefficient <- specification$coefficient
  if (method == "voom_qw") {
    v <- voomWithQualityWeights(y, design, plot = FALSE)
    fit <- eBayes(lmFit(v, design), robust = TRUE)
    tab <- topTable(fit, coef = coefficient, number = Inf, sort.by = "none")
    se <- fit$stdev.unscaled[, coefficient] * sqrt(fit$s2.post)
    critical <- qt(0.975, df = fit$df.total)
    answer <- data.table(peak_coordinate = rownames(tab), logFC = tab$logFC,
      standard_error = as.numeric(se), ci_lower = tab$logFC - critical * se,
      ci_upper = tab$logFC + critical * se, pvalue = tab$P.Value,
      standard_error_method = "moderated_t")
    weight <- v$targets$sample.weights %||% rep(NA_real_, nrow(donors))
  } else {
    require_true(method == "edger_ql", "Unknown DA method")
    y <- estimateDisp(y, design, robust = TRUE)
    fit <- glmQLFit(y, design, robust = TRUE)
    tab <- glmQLFTest(fit, coef = coefficient)$table
    # A QL F statistic is not a Wald t statistic. Do not fabricate SE = logFC/sqrt(F).
    answer <- data.table(peak_coordinate = rownames(tab), logFC = tab$logFC,
      standard_error = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
      pvalue = tab$PValue, standard_error_method = "not_reported_by_qlf")
    weight <- rep(NA_real_, nrow(donors))
  }
  answer[, qvalue := p.adjust(pvalue, "BH")]
  qc <- data.table(donor_id = donors$donor_id, condition = donors$condition,
    raw_library_size = y$samples$lib.size, tmm_normalization_factor = y$samples$norm.factors,
    effective_library_size = y$samples$lib.size * y$samples$norm.factors,
    sample_quality_weight = weight, design_rank = qr(design)$rank,
    n_model_columns = ncol(design), residual_df = nrow(design) - ncol(design), n_tested_peaks = nrow(answer))
  list(results = answer, qc = qc)
}

append_metadata <- function(x, cohort, lineage, representation, rule, method) {
  x[, `:=`(cohort = cohort, lineage = lineage, representation = representation, filter = rule, method = method)]
  x
}

join_donor_metadata <- function(donors, path) {
  if (is.null(path) || !file.exists(path)) return(copy(donors))
  meta <- fread(path)
  require_true("donor_id" %in% names(meta) && !anyDuplicated(meta$donor_id), "Invalid donor metadata")
  out <- copy(donors)
  for (column in setdiff(names(meta), names(donors))) out[[column]] <- meta[[column]][match(out$donor_id, meta$donor_id)]
  if ("condition" %in% names(meta)) {
    expected <- toupper(meta$condition[match(out$donor_id, meta$donor_id)])
    require_true(!anyNA(expected) && identical(expected, out$condition), "Condition metadata conflict")
  }
  out
}

covariate_availability <- function(donors, column) {
  if (!column %in% names(donors)) return(data.table(covariate = column, status = "missing_column",
    reason = "not_recorded_in_supplied_donor_metadata", n_donors = nrow(donors), n_nonmissing = 0L))
  values <- donors[[column]]
  complete <- !is.na(values) & nzchar(trimws(as.character(values)))
  if (is.numeric(values)) complete <- complete & is.finite(values)
  status <- if (!any(complete)) "all_missing" else if (!all(complete)) "incomplete" else if (length(unique(values)) < 2L) "constant" else "estimable"
  reason <- if (status == "estimable") "" else "retain_all_primary_donors; no_complete_case_donor_deletion"
  if (status == "estimable") {
    design <- tryCatch(make_da_design(donors, column), error = function(e) e)
    if (inherits(design, "error")) { status <- "nonestimable"; reason <- conditionMessage(design) }
  }
  data.table(covariate = column, status = status, reason = reason, n_donors = nrow(donors), n_nonmissing = sum(complete))
}

arm_seed <- function(seed, cohort, lineage, representation, rule) {
  as.integer(seed + match(cohort, c("GSE244832", "GSE281367")) * 100000L +
    match(lineage, c("hepatocyte", "stellate", "macrophage")) * 1000L +
    match(representation, c("fragment", "binary")) * 100L + match(rule, c("primary", "prevalence")) * 10L)
}

baseline_reproduction <- function(answer, baseline_root, cohort, lineage, tolerance) {
  path <- file.path(baseline_root, "da", paste0("cohort_primary_", cohort, "_", lineage, ".tsv"))
  require_true(file.exists(path), paste("Corrected baseline missing:", path))
  reference <- fread(path)
  require_true(setequal(reference$peak_coordinate, answer$peak_coordinate), "Primary baseline peak universe drift")
  ref <- reference[match(answer$peak_coordinate, peak_coordinate)]
  rbindlist(lapply(c("logFC", "standard_error", "pvalue", "qvalue"), function(column) {
    require_true(column %in% names(ref), paste("Baseline lacks", column))
    delta <- max(abs(answer[[column]] - ref[[column]]), na.rm = TRUE)
    require_true(is.finite(delta) && delta <= tolerance, paste("Primary baseline mismatch:", column, delta))
    data.table(cohort = cohort, lineage = lineage, quantity = column, n_peaks = nrow(answer), max_absolute_difference = delta,
      tolerance = tolerance, status = "pass", baseline_path = path)
  }))
}

lodo_da <- function(counts, donors, method, full, workers = 1L) {
  rows <- parallel_da(seq_len(nrow(donors)), function(index) {
    reduced <- donors[-index]
    if (min(table(factor(reduced$condition, levels = c("NORMAL", "MASH")))) < 3L) {
      return(data.table(omitted_donor = donors$donor_id[index], status = "untestable", reason = "fewer_than_three_donors_per_group"))
    }
    value <- tryCatch(fit_da(counts[-index, , drop = FALSE], reduced, method)$results, error = function(e) e)
    if (inherits(value, "error")) return(data.table(omitted_donor = donors$donor_id[index], status = "untestable", reason = conditionMessage(value)))
    reference <- full[match(value$peak_coordinate, peak_coordinate)]
    data.table(omitted_donor = donors$donor_id[index], status = "testable", reason = "", n_peaks = nrow(value),
      effect_spearman = cor(value$logFC, reference$logFC, method = "spearman"),
      median_absolute_effect_change = median(abs(value$logFC - reference$logFC)),
      maximum_absolute_effect_change = max(abs(value$logFC - reference$logFC)),
      sign_agreement = mean(sign(value$logFC) == sign(reference$logFC)))
  }, workers)
  rbindlist(rows, fill = TRUE)
}

null_da <- function(counts, donors, method, reps, max_peaks, seed, workers = 1L) {
  if (reps < 1L) return(data.table(status = "not_run", reason = "null_reps_zero"))
  set.seed(seed)
  keep <- sort(sample(seq_len(ncol(counts)), min(ncol(counts), max_peaks)))
  counts <- counts[, keep, drop = FALSE]
  reduced <- matrix(1, nrow(donors), 1L)
  y <- calcNormFactors(DGEList(counts = t(counts)), method = "TMM")
  y <- estimateDisp(y, reduced, robust = TRUE)
  fit <- glmFit(y, reduced)
  mu <- fit$fitted.values
  dispersion <- y$tagwise.dispersion
  rows <- parallel_da(seq_len(reps), function(iteration) {
    # Each method receives exactly the same simulated counts for this arm.
    set.seed(seed + iteration)
    simulated <- matrix(rnbinom(length(mu), mu = as.vector(mu),
      size = rep(1 / pmax(dispersion, 1e-8), ncol(mu))), nrow(mu), ncol(mu))
    dimnames(simulated) <- dimnames(mu)
    result <- tryCatch(fit_da(t(simulated), donors, method)$results, error = function(e) e)
    if (inherits(result, "error")) return(data.table(replicate = iteration, status = "failed", reason = conditionMessage(result), n_tested_peaks = nrow(mu)))
    data.table(replicate = iteration, status = "testable", reason = "", n_tested_peaks = nrow(result),
      n_nominal_p05 = sum(result$pvalue < 0.05), n_bh_q05 = sum(result$qvalue < 0.05),
      simulated_count_sum = sum(simulated), simulated_count_sum_squares = sum(as.double(simulated)^2),
      null_model = "fitted_intercept_NB_global_null_fixed_features", seed = seed,
      interpretation = "model_based_diagnostic_not_biological_ground_truth")
  }, workers)
  rbindlist(rows, fill = TRUE)
}

restricted_label_design <- function(donors) {
  recorded <- names(donors)[grepl("^batch", names(donors), ignore.case = TRUE)]
  incomplete <- recorded[vapply(donors[, ..recorded], function(x) anyNA(x) || any(!nzchar(trimws(as.character(x)))), logical(1))]
  columns <- setdiff(recorded, incomplete)
  columns <- columns[vapply(donors[, ..columns], function(x) length(unique(x)) > 1L, logical(1))]
  blocks <- if (length(columns)) interaction(donors[, ..columns], drop = TRUE) else factor(rep("all", nrow(donors)))
  strata <- split(seq_len(nrow(donors)), blocks)
  n_mash <- vapply(strata, function(index) sum(donors$condition[index] == "MASH"), integer(1))
  ways <- mapply(function(index, n) choose(length(index), n), strata, n_mash)
  list(strata = strata, n_mash = n_mash, ways = ways, assignment_count = prod(ways),
    columns = columns, incomplete_columns = incomplete,
    restriction = if (length(columns)) paste(columns, collapse = "+") else "none_variable_recorded",
    assumption = "conditional_exchangeability_within_recorded_batches_assumed; unmeasured_confounding_unaddressed")
}

restricted_label_choices <- function(design) {
  mapply(function(index, n) {
    if (n == 0L) list(integer()) else if (n == length(index)) list(index) else combn(index, n, simplify = FALSE)
  }, design$strata, design$n_mash, SIMPLIFY = FALSE)
}

restricted_labels <- function(design, choices, iteration, n_donors) {
  labels <- rep("NORMAL", n_donors)
  rank <- iteration - 1
  for (block in seq_along(choices)) {
    chosen <- rank %% design$ways[block] + 1
    rank <- floor(rank / design$ways[block])
    labels[choices[[block]][[chosen]]] <- "MASH"
  }
  labels
}

directional_da <- function(counts, donors, method, discoveries, observed, reps, seed, workers = 1L) {
  design <- restricted_label_design(donors)
  unavailable <- function(status) list(p = NA_real_, n = 0L, status = status, mode = "not_evaluated",
    assignment_count = design$assignment_count, restriction = design$restriction, assumption = design$assumption)
  if (reps < 1L || !length(discoveries)) return(unavailable(if (length(discoveries)) "not_run" else "no_reference_discoveries"))
  if (length(design$incomplete_columns)) return(unavailable("untestable_incomplete_recorded_batch_metadata"))
  if (design$assignment_count <= 1) return(unavailable("untestable_no_within_batch_label_exchange"))
  exact <- is.finite(design$assignment_count) && design$assignment_count <= reps
  n_evaluated <- if (exact) as.integer(design$assignment_count) else reps
  choices <- if (exact) restricted_label_choices(design) else NULL
  exceed <- parallel_da(seq_len(n_evaluated), function(iteration) {
    set.seed(seed + iteration)
    shuffled <- copy(donors)
    if (exact) shuffled$condition <- restricted_labels(design, choices, iteration, nrow(donors))
    else for (index in design$strata) shuffled$condition[index] <- sample(donors$condition[index])
    fit <- fit_da(counts, shuffled, method)$results
    effect <- fit$logFC[match(names(discoveries), fit$peak_coordinate)]
    as.integer(sum(sign(effect) == sign(discoveries)) >= observed)
  }, workers)
  p <- if (exact) sum(unlist(exceed)) / n_evaluated else (1 + sum(unlist(exceed))) / (1 + n_evaluated)
  status <- if (length(design$columns)) {
    if (exact) "batch_restricted_exact" else "batch_restricted_monte_carlo"
  } else if (exact) "donor_relabeling_exact_no_recorded_batch" else "donor_relabeling_no_recorded_batch"
  list(p = p, n = n_evaluated, status = status, mode = if (exact) "exact" else "monte_carlo",
    assignment_count = design$assignment_count, restriction = design$restriction, assumption = design$assumption)
}

adjust_directional_rows <- function(summary) {
  summary <- copy(summary)
  tested_status <- c("batch_restricted_monte_carlo", "donor_relabeling_no_recorded_batch",
    "batch_restricted_exact", "donor_relabeling_exact_no_recorded_batch")
  summary[cohort == "cross_cohort", c("directional_pvalue", "directional_qvalue", "directional_tested", "directional_family_size") := {
    tested <- status %chin% tested_status & !is.na(directional_reps) & directional_reps > 0L &
      is.finite(directional_pvalue) & directional_pvalue >= 0 & directional_pvalue <= 1
    # Unavailable hypotheses remain in the fixed three-lineage BH family, but
    # their internal p=1 placeholders must never appear as measured test results.
    internal <- fifelse(tested, directional_pvalue, 1)
    adjusted <- p.adjust(internal, "BH", n = 3L)
    list(fifelse(tested, directional_pvalue, NA_real_), fifelse(tested, adjusted, NA_real_), tested, rep(3L, .N))
  }, by = .(representation, filter, method)]
  summary
}

parse_da_args <- function(args) {
  values <- list(null_reps = NULL, skip_lodo = FALSE, binary_root = NULL, workers = "1",
    no_baseline_check = FALSE, sensitivity_label = NULL, counts_root = NULL)
  i <- 1L
  while (i <= length(args)) {
    key <- args[i]
    if (key == "--skip-lodo") { values$skip_lodo <- TRUE; i <- i + 1L; next }
    if (key == "--no-baseline-check") { values$no_baseline_check <- TRUE; i <- i + 1L; next }
    require_true(key %in% c("--config", "--output", "--binary-root", "--null-reps", "--workers", "--counts-root", "--sensitivity-label") && i < length(args),
      "Usage: compare_da.R --config JSON --output NEWDIR [--binary-root ROOT] [--null-reps N] [--workers N] [--skip-lodo] [--counts-root ROOT --no-baseline-check --sensitivity-label cell_downsample]")
    values[[gsub("-", "_", substring(key, 3))]] <- args[i + 1L]
    i <- i + 2L
  }
  require_true(!is.null(values$config) && !is.null(values$output), "--config and --output are required")
  values
}

run_da <- function(args) {
  require_true(nzchar(Sys.getenv("SLURM_JOB_ID")), "Run analytical comparisons in a Slurm compute allocation")
  cfg <- fromJSON(args$config, simplifyVector = FALSE)
  resource_config <- cfg$execution
  if (!is.null(cfg$da)) cfg <- modifyList(cfg, cfg$da)
  workers <- as.integer(args$workers %||% 1L)
  allocated <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = as.character(resource_config$da_cpus %||% 1L)))
  allowed <- min(allocated, as.integer(resource_config$da_cpus %||% allocated))
  require_true(is.finite(workers) && workers >= 1L && as.numeric(args$workers %||% 1L) == workers && workers <= allowed,
    paste("--workers must be a positive integer within the compute allocation/configured limit:", allowed))
  if (workers > 1L) {
    require_true(Sys.getenv("OPENBLAS_NUM_THREADS") == "1" && Sys.getenv("OMP_NUM_THREADS") == "1",
      "Start parallel resampling with OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1")
    data.table::setDTthreads(1L)
  }
  root <- cfg$project_root %||% "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
  baseline <- cfg$baseline_root %||% file.path(root, "Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-candidate-2026-08-11-r1")
  corrected <- cfg$corrected_baseline_root %||% file.path(root, "Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-statistics-correction-20261022")
  counts_root <- args$counts_root %||% cfg$counts_root %||% file.path(baseline, "counts")
  if (isTRUE(args$no_baseline_check)) {
    require_true(!is.null(args$sensitivity_label) && grepl("downsample", args$sensitivity_label, fixed = TRUE),
      "Skipping baseline reproduction requires a clearly labeled downsample sensitivity")
    require_true(normalizePath(counts_root, mustWork = TRUE) != normalizePath(file.path(baseline, "counts"), mustWork = TRUE),
      "Baseline reproduction is mandatory for the original primary counts")
  }
  cohorts <- if (is.null(cfg$cohorts)) c("GSE244832", "GSE281367") else names(cfg$cohorts)
  if (is.null(cohorts)) cohorts <- unlist(cfg$cohorts)
  lineages <- unlist(cfg$lineages %||% list("hepatocyte", "stellate", "macrophage"))
  require_true(identical(cohorts, c("GSE244832", "GSE281367")), "Prescribed cohorts must remain in order")
  require_true(setequal(lineages, c("hepatocyte", "stellate", "macrophage")), "Preserve the three-lineage family")
  seed <- as.integer(cfg$seed %||% 42L)
  null_reps <- as.integer(args$null_reps %||% cfg$null_reps %||% 20L)
  directional_reps <- as.integer(cfg$directional_reps %||% null_reps)
  require_true(null_reps >= 0 && directional_reps >= 0, "Negative resampling count")
  out <- args$output
  require_true(!file.exists(out), "Refusing to overwrite output")
  dir.create(out, recursive = TRUE)
  set.seed(seed)
  all_results <- all_qc <- filter_qc <- summaries <- checks <- lodo <- sensitivities <- availability <- nulls <- replication <- list()
  representations <- c(fragment = counts_root)
  binary_root <- args$binary_root %||% cfg$binary_root
  if (!is.null(binary_root)) representations <- c(representations, binary = binary_root)
  for (representation in names(representations)) for (lineage in lineages) {
    message("[DA] ", representation, " / ", lineage)
    inputs <- lapply(cohorts, function(cohort) {
      x <- load_da_counts(representations[[representation]], cohort, lineage)
      meta <- cfg$cohort_metadata[[cohort]] %||% file.path(root, "data", cohort, "metadata/donor_pairing.csv")
      x$donors <- join_donor_metadata(x$donors, meta)
      x$donors <- join_donor_metadata(x$donors, cfg$cohort_covariate_metadata[[cohort]])
      x
    })
    names(inputs) <- cohorts
    for (cohort in cohorts) for (column in c("batch_mm", "batch", "sex", "age", "BMI", "bmi")) {
      availability[[length(availability) + 1L]] <- append_metadata(covariate_availability(inputs[[cohort]]$donors, column),
        cohort, lineage, representation, "metadata", "not_applicable")
    }
    require_true(all(vapply(inputs, function(x) min(table(factor(x$donors$condition,
      levels = c("NORMAL", "MASH")))) >= 4L, logical(1))), "Fewer than four eligible donors per group")
    filter_inputs <- if (representation == "fragment") inputs else setNames(lapply(cohorts,
      function(cohort) load_da_counts(counts_root, cohort, lineage)), cohorts)
    if (representation == "binary") for (cohort in cohorts) {
      require_true(identical(rownames(inputs[[cohort]]$counts), rownames(filter_inputs[[cohort]]$counts)) &&
        identical(colnames(inputs[[cohort]]$counts), colnames(filter_inputs[[cohort]]$counts)),
        "Binary and fragment rosters/features differ")
      require_true(identical(inputs[[cohort]]$donors$n_cells, filter_inputs[[cohort]]$donors$n_cells),
        "Binary and fragment nucleus denominators differ")
      difference <- filter_inputs[[cohort]]$counts - inputs[[cohort]]$counts
      require_true(all(difference@x >= 0), "Accessible nucleus counts exceed fragment counts")
      for (index in seq_len(nrow(inputs[[cohort]]$counts))) {
        require_true(max(inputs[[cohort]]$counts[index, ]) <= inputs[[cohort]]$donors$n_cells[index],
          "Binary count exceeds the number of nuclei in a donor")
      }
    }
    for (rule in c("primary", "prevalence")) {
      # Lock representation comparisons to fragment-defined condition-blind peaks.
      retained <- lapply(filter_inputs, function(x) colnames(x$counts)[filter_da_peaks(x$counts, rule)])
      common <- Reduce(intersect, retained)
      filter_qc[[length(filter_qc) + 1L]] <- data.table(representation = representation, lineage = lineage,
        filter = rule, feature_filter_representation = "fragment", n_retained_gse244832 = length(retained[[1L]]),
        n_retained_gse281367 = length(retained[[2L]]), n_jointly_testable = length(common))
      for (method in c("voom_qw", "edger_ql")) {
        fits <- list()
        for (cohort in cohorts) {
          x <- inputs[[cohort]]
          message("[DA arm] ", cohort, " / ", representation, " / ", lineage, " / ", rule, " / ", method,
            " / ", length(common), " peaks")
          if (!length(common)) {
            summaries[[length(summaries) + 1L]] <- data.table(cohort = cohort, lineage = lineage,
              representation = representation, filter = rule, method = method, status = "untestable",
              reason = "no_jointly_testable_peaks", n_tested_peaks = 0L,
              n_normal = sum(x$donors$condition == "NORMAL"), n_mash = sum(x$donors$condition == "MASH"))
            next
          }
          counts <- x$counts[, common, drop = FALSE]
          fit <- fit_da(counts, x$donors, method)
          fits[[cohort]] <- fit$results
          result <- copy(fit$results)
          result[, `:=`(effect_units = if (representation == "fragment") "log2_normalized_fragment_count" else "log2_normalized_accessible_nucleus_count",
            bh_family = "all_jointly_testable_peaks_within_cohort_lineage_representation_filter_method",
            evidence_status = "candidate_sensitivity_not_adopted")]
          all_results[[length(all_results) + 1L]] <- append_metadata(result, cohort, lineage, representation, rule, method)
          all_qc[[length(all_qc) + 1L]] <- append_metadata(fit$qc, cohort, lineage, representation, rule, method)
          summaries[[length(summaries) + 1L]] <- data.table(cohort = cohort, lineage = lineage,
            representation = representation, filter = rule, method = method, status = "testable", reason = "",
            n_tested_peaks = nrow(result), n_bh_q05 = sum(result$qvalue < 0.05),
            n_normal = sum(x$donors$condition == "NORMAL"), n_mash = sum(x$donors$condition == "MASH"))
          if (representation == "fragment" && rule == "primary" && method == "voom_qw") {
            checks[[length(checks) + 1L]] <- if (isTRUE(args$no_baseline_check)) {
              data.table(cohort = cohort, lineage = lineage, status = "not_applicable_downsample_sensitivity",
                reason = args$sensitivity_label, n_peaks = nrow(result))
            } else baseline_reproduction(result, corrected, cohort, lineage, as.numeric(cfg$baseline_tolerance %||% 1e-7))
            fwrite(rbindlist(checks, fill = TRUE), file.path(out, "baseline_reproduction.tsv"), sep = "\t")
            message("[DA baseline] ", cohort, " / ", lineage, " / ", checks[[length(checks)]]$status[1L])
          }
          if (!isTRUE(args$skip_lodo)) lodo[[length(lodo) + 1L]] <- append_metadata(
            lodo_da(counts, x$donors, method, fit$results, workers), cohort, lineage, representation, rule, method)
          for (column in c("batch_mm", "batch", "sex", "age", "BMI", "bmi")) {
            available <- covariate_availability(x$donors, column)
            if (available$status != "estimable") next
            adjusted <- tryCatch(fit_da(counts, x$donors, method, column)$results, error = function(e) e)
            if (inherits(adjusted, "error")) {
              value <- data.table(covariate = column, status = "untestable", reason = conditionMessage(adjusted))
            } else {
              value <- adjusted[, .(peak_coordinate, adjusted_logFC = logFC, adjusted_pvalue = pvalue, adjusted_qvalue = qvalue)]
              value[, `:=`(covariate = column, status = "testable", reason = "")]
            }
            sensitivities[[length(sensitivities) + 1L]] <- append_metadata(value, cohort, lineage, representation, rule, method)
          }
          nulls[[length(nulls) + 1L]] <- append_metadata(null_da(counts, x$donors, method,
            null_reps, as.integer(cfg$null_max_peaks %||% 2000L), arm_seed(seed, cohort, lineage, representation, rule), workers),
            cohort, lineage, representation, rule, method)
        }
        if (length(fits) == 2L) {
          joined <- merge(fits[[1L]], fits[[2L]], by = "peak_coordinate", suffixes = c("_gse244832", "_gse281367"))
          joined[, evidence_state := fifelse(qvalue_gse244832 < 0.05 & qvalue_gse281367 < 0.05,
            fifelse(logFC_gse244832 * logFC_gse281367 > 0, "supported", "discordant"),
            fifelse(xor(qvalue_gse244832 < 0.05, qvalue_gse281367 < 0.05), "source_dependent", "indeterminate"))]
          joined[, `:=`(lineage = lineage, representation = representation, filter = rule, method = method,
            evidence_status = "candidate_only_shared_outcome_blind_feature_ascertainment")]
          replication[[length(replication) + 1L]] <- joined
          reference <- fits[[2L]][qvalue < 0.05]
          discoveries <- setNames(reference$logFC, reference$peak_coordinate)
          effect <- fits[[1L]]$logFC[match(names(discoveries), fits[[1L]]$peak_coordinate)]
          observed <- sum(sign(effect) == sign(discoveries))
          d <- directional_da(inputs[[1L]]$counts[, common, drop = FALSE], inputs[[1L]]$donors,
            method, discoveries, observed, directional_reps, arm_seed(seed, cohorts[1L], lineage, representation, rule) + 10000L, workers)
          summaries[[length(summaries) + 1L]] <- data.table(cohort = "cross_cohort", lineage = lineage,
            representation = representation, filter = rule, method = method, status = d$status,
            n_tested_peaks = length(common), n_reference_discoveries = length(discoveries),
            n_direction_concordant = observed, n_supported = sum(joined$evidence_state == "supported"),
            directional_pvalue = d$p, directional_reps = d$n, directional_mode = d$mode,
            directional_assignment_count = d$assignment_count, directional_restriction = d$restriction,
            directional_assumption = d$assumption)
        } else summaries[[length(summaries) + 1L]] <- data.table(cohort = "cross_cohort", lineage = lineage,
          representation = representation, filter = rule, method = method, status = "untestable",
          reason = "no_jointly_testable_peaks", n_tested_peaks = 0L, directional_pvalue = 1, directional_reps = 0L)
      }
    }
  }
  summary <- adjust_directional_rows(rbindlist(summaries, fill = TRUE))
  write_rows <- function(rows, name) {
    tab <- if (length(rows)) rbindlist(rows, fill = TRUE) else data.table(status = "not_run")
    fwrite(tab, file.path(out, name), sep = "\t", compress = if (endsWith(name, ".gz")) "gzip" else "none")
  }
  write_rows(all_results, "peak_results.tsv.gz")
  write_rows(all_qc, "model_qc.tsv.gz")
  write_rows(filter_qc, "filter_qc.tsv")
  fwrite(summary, file.path(out, "cohort_summary.tsv"), sep = "\t")
  write_rows(replication, "replication.tsv.gz")
  write_rows(checks, "baseline_reproduction.tsv")
  write_rows(lodo, "lodo_summary.tsv")
  write_rows(sensitivities, "covariate_sensitivity.tsv.gz")
  write_rows(availability, "covariate_availability.tsv")
  write_rows(nulls, "null_diagnostics.tsv")
  write_json(list(status = "candidate_complete_not_adopted", seed = seed, workers = workers, null_reps = null_reps,
    directional_reps = directional_reps, binary_status = if (is.null(binary_root)) "not_run_missing_per_cell_binary_counts" else "completed",
    directional_resampling_policy = "exact_including_identity_once_if_assignment_space_le_requested_reps; plus_one_Monte_Carlo_otherwise",
    corrected_baseline_root = corrected, counts_root = counts_root, primary_method = "existing_v3_voom_quality_weights",
    biological_unit = "donor", directional_bh_family = "three_prespecified_lineages_per_representation_filter_method",
    feature_filter_representation = "fragment_for_both_count_representations",
    qlf_intervals = "not_reported; donor_influence_is_a_sensitivity_not_a_confidence_interval",
    null_max_peaks = as.integer(cfg$null_max_peaks %||% 2000L),
    sensitivity_label = args$sensitivity_label %||% "original_counts",
    baseline_reproduction_required = !isTRUE(args$no_baseline_check),
    null_interpretation = "NB_simulation_is_a_model_based_diagnostic_not_actual_biological_FDR"),
    file.path(out, "run_summary.json"), auto_unbox = TRUE, pretty = TRUE)
  writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
}

if (sys.nframe() == 0L) run_da(parse_da_args(commandArgs(trailingOnly = TRUE)))
