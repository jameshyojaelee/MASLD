#!/usr/bin/env Rscript

# Precontract composition acceptance for BULK-PROGRAM-MAP-v9.
#
# WHY THIS EXISTS
# The v8 composition arm used the deposited MuSiC estimates. Those estimates are
# degenerate on the real bulk substrate: median hepatocyte 0.942, cholangiocytes
# exactly zero in 97.8% of samples, fibroblasts 87%, macrophages 73%,
# endothelial 65%. Because cholangiocyte and fibroblast fractions are both
# identically zero in GSE135251, the 1e-6 pseudocount collapsed their two
# hepatocyte log-ratios onto the same vector, the design lost rank, and
# GSE135251 (214 of 469 discovery samples) silently left the composition model.
# The project's pseudobulk benchmark does not catch this: MuSiC is well
# calibrated there (hepatocyte 0.626 against a truth of 0.623). That benchmark
# is a cell-number simplex while bulk RNA is mRNA-mass weighted, so it does not
# transfer to the real substrate.
#
# WHAT THIS SCRIPT LEARNED THE HARD WAY
# Absence of exact zeros is not sufficient. BayesPrism produces no exact zeros
# but does produce solver underflow: in GSE135251 the 5th-percentile
# cholangiocyte proportion is 3.3e-13. log(3.3e-13 / hep) against
# log(1e-6 / hep) is a fifteen-unit swing driven entirely by solver noise, so a
# lineage can be non-zero everywhere and still carry no usable signal. The
# dynamic-range floor below, not the zero-fraction screen, is what decides
# whether a lineage may enter the model as a covariate.
#
# CRITERIA PROVENANCE, STATED PLAINLY
# The zero-fraction, coverage and stage-direction criteria were evaluated during
# planning on 2026-08-12 and are not blind. Thresholds sit far from every
# observed value so the decision does not turn on where a cut was placed: MuSiC
# records a zero fraction of 0.978 against a threshold of 0.05, BayesPrism
# records 0.000. The stage-direction criterion uses Kleiner stage, which is an
# outcome. It constrains only the composition estimator and the lineage set,
# never program selection, and its use must be repeated wherever a
# composition-adjusted program result is reported.
#
# THIS SCRIPT PUBLISHES ITS AUDIT BEFORE IT DECIDES. An earlier version raised
# an error when no estimator qualified, which destroyed the evidence needed to
# see why. A run in which nothing qualifies is a result and must be readable.

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_dir <- dirname(normalizePath(script_arg))
source(file.path(script_dir, "config.R"))
source(file.path(script_dir, "analysis_lib.R"))

options(digits = 17)

COVARIATE_LINEAGES <- c("Cholangiocytes", "Fibroblasts", "Macrophages", "T cells")
DENOMINATOR_LINEAGE <- "Hepatocytes"
DIRECTION_EXPECTED <- c(Hepatocytes = -1, Fibroblasts = 1, Cholangiocytes = 1)

# ------------------------------------------------------------------ estimators

# Read per-cohort proportions, restricted to the cohorts in the corrected
# census. This restriction is load-bearing, not tidiness: GSE156918, GSE205974
# and inhouse_MCD were deconvolved against a different reference with a
# different lineage set (17 columns including Kupffer cells and NK cells).
# Binding them in produced a column union in which no sample row was complete,
# so every estimator scored zero coverage. Only lineage sets that agree exactly
# across the census cohorts are usable, and that agreement is asserted here.
read_per_cohort <- function(suffix, cohorts) {
  rows <- list()
  lineage_sets <- list()
  for (cohort in cohorts) {
    path <- file.path(deconvolution_root, cohort, paste0(cohort, suffix))
    if (!file.exists(path)) next
    value <- suppressWarnings(fread(path, check.names = FALSE))
    setnames(value, 1L, "sample_id")
    value[, sample_id := as.character(sample_id)]
    value[, dataset := cohort]
    lineage_sets[[cohort]] <- sort(setdiff(
      names(value)[vapply(value, is.numeric, logical(1))], c("sample_id", "dataset")
    ))
    rows[[length(rows) + 1L]] <- value
  }
  if (!length(rows)) return(NULL)
  reference <- lineage_sets[[1L]]
  divergent <- names(lineage_sets)[
    !vapply(lineage_sets, function(x) identical(x, reference), logical(1))
  ]
  if (length(divergent)) {
    fail("Lineage sets disagree across census cohorts for ", suffix, ": ",
         paste(divergent, collapse = ", "))
  }
  unique(rbindlist(rows, use.names = TRUE), by = c("sample_id", "dataset"))
}

read_music_deposited <- function() {
  value <- fread(composition_candidate, check.names = FALSE)
  contract <- fread(composition_testability)
  keep <- contract[testability == "testable", celltype]
  unique(value[, c("sample_id", "dataset", keep), with = FALSE], by = c("sample_id", "dataset"))
}

lineage_columns <- function(value) {
  setdiff(names(value)[vapply(value, is.numeric, logical(1))], c("sample_id", "dataset"))
}

# --------------------------------------------------------------------- scoring

census_samples <- function() {
  dge <- readRDS(input_dge)
  meta <- as.data.table(readRDS(input_meta))
  samples <- colnames(dge)
  dataset <- as.character(dge$samples$dataset)
  if (length(dataset) != length(samples) || anyNA(dataset)) {
    dataset <- meta$dataset[match(samples, meta$sample_id)]
  }
  assert_true(!anyNA(dataset), "Every census sample must carry a dataset")
  rm(dge); gc()
  data.table(sample_id = samples, dataset = dataset)
}

audit_method <- function(label, proportions, census, phenotype) {
  empty <- list(
    coverage = data.table(), lineage = data.table(), dynamic = data.table(),
    direction = data.table(), collinearity = data.table(), matched = NULL,
    lineages = character(0)
  )
  if (is.null(proportions)) {
    empty$coverage <- data.table(method = label, dataset = NA_character_,
                                 n_census = NA_integer_, n_covered = 0L, coverage = 0)
    return(empty)
  }
  matched <- merge(census, proportions, by = c("sample_id", "dataset"),
                   all.x = TRUE, sort = FALSE)
  matched <- matched[match(census$sample_id, sample_id)]
  lineages <- lineage_columns(matched)
  matched[, .row_complete := rowSums(!is.na(.SD)) == length(lineages), .SDcols = lineages]

  coverage <- matched[, .(
    method = label, n_census = .N, n_covered = sum(.row_complete),
    coverage = mean(.row_complete)
  ), by = dataset]

  observed <- matched[.row_complete == TRUE]
  M <- as.matrix(observed[, ..lineages]); storage.mode(M) <- "double"
  lineage_audit <- data.table(
    method = label, celltype = lineages,
    zero_fraction = colMeans(M == 0),
    median_fraction = apply(M, 2L, stats::median),
    max_fraction = apply(M, 2L, max)
  )

  # Dynamic range decides covariate eligibility. A lineage whose lower tail sits
  # at solver underflow cannot carry a log-ratio no matter how few exact zeros
  # it has.
  dynamic <- rbindlist(lapply(lineages, function(lineage) {
    observed[dataset %in% composition_required_cohorts, .(
      method = label, celltype = lineage,
      p05 = stats::quantile(get(lineage), 0.05, names = FALSE),
      p50 = stats::median(get(lineage)),
      fraction_below_floor = mean(get(lineage) < composition_dynamic_range_floor)
    ), by = dataset]
  }))

  direction <- rbindlist(lapply(intersect(names(DIRECTION_EXPECTED), lineages), function(lineage) {
    rbindlist(lapply(discovery_cohorts, function(cohort) {
      d <- merge(observed[dataset == cohort, c("sample_id", lineage), with = FALSE],
                 phenotype[, .(sample_id, fibrosis_stage)], by = "sample_id")
      d <- d[is.finite(fibrosis_stage)]
      rho <- if (nrow(d) > 10L && stats::sd(d[[lineage]]) > 0) {
        suppressWarnings(stats::cor(d[[lineage]], d$fibrosis_stage, method = "spearman"))
      } else NA_real_
      data.table(
        method = label, celltype = lineage, cohort = cohort, n = nrow(d),
        spearman_fibrosis = rho, expected_sign = DIRECTION_EXPECTED[[lineage]],
        agrees = isTRUE(is.finite(rho) && sign(rho) == DIRECTION_EXPECTED[[lineage]])
      )
    }))
  }))

  list(coverage = coverage, lineage = lineage_audit, dynamic = dynamic,
       direction = direction, matched = matched, lineages = lineages)
}

logratio_collinearity <- function(label, matched, eligible) {
  if (length(eligible) < 2L) return(data.table())
  rbindlist(lapply(discovery_cohorts, function(cohort) {
    d <- matched[dataset == cohort & .row_complete == TRUE]
    if (!nrow(d) || !DENOMINATOR_LINEAGE %in% names(d)) return(data.table())
    den <- d[[DENOMINATOR_LINEAGE]]
    den[den <= 0] <- composition_pseudocount
    L <- vapply(eligible, function(l) {
      v <- d[[l]]; v[v <= 0] <- composition_pseudocount; log(v / den)
    }, numeric(nrow(d)))
    r <- suppressWarnings(stats::cor(L))
    data.table(
      method = label, cohort = cohort, n_eligible = length(eligible),
      max_abs_correlation = if (length(eligible) > 1L) max(abs(r[upper.tri(r)])) else NA_real_,
      design_rank = qr(cbind(1, L))$rank, design_columns = length(eligible) + 1L,
      full_rank = qr(cbind(1, L))$rank == (length(eligible) + 1L)
    )
  }))
}

# ------------------------------------------------------------------------ main

if (dir.exists(composition_acceptance_root)) {
  fail("Refusing to overwrite acceptance seal: ", composition_acceptance_root)
}

message("[1/6] Building corrected census")
census <- census_samples()
assert_true(nrow(census) == 1257L, "Corrected census must contain 1,257 QC-passing samples")
phenotype <- as.data.table(readRDS(input_meta))[, .(sample_id, fibrosis_stage, nas_score)]

message("[2/6] Reading candidate estimators")
census_cohorts <- sort(unique(census$dataset))
message("      census cohorts: ", paste(census_cohorts, collapse = ", "))
methods <- list(
  music_deposited = read_music_deposited(),
  bayesprism = read_per_cohort("_bayesprism_proportions.tsv", census_cohorts),
  rectangle = read_per_cohort("_rectangle_proportions.tsv", census_cohorts)
)

message("[3/6] Auditing every estimator")
audits <- lapply(names(methods), function(label) {
  message("      - ", label)
  audit_method(label, methods[[label]], census, phenotype)
})
names(audits) <- names(methods)

coverage_table <- rbindlist(lapply(audits, `[[`, "coverage"), use.names = TRUE, fill = TRUE)
lineage_table <- rbindlist(lapply(audits, `[[`, "lineage"), use.names = TRUE, fill = TRUE)
dynamic_table <- rbindlist(lapply(audits, `[[`, "dynamic"), use.names = TRUE, fill = TRUE)
direction_table <- rbindlist(lapply(audits, `[[`, "direction"), use.names = TRUE, fill = TRUE)

message("[4/6] Applying the frozen gate and criteria")

# Gate first. An estimator that cannot cover every discovery cohort and the
# holdout is not a candidate, however well it scores elsewhere. This is what
# removes Rectangle, which has no GSE174478 or GSE193066 output.
gate <- rbindlist(lapply(names(methods), function(label) {
  cov <- coverage_table[method == label & dataset %in% composition_required_cohorts]
  data.table(
    method = label,
    n_required_cohorts_present = nrow(cov),
    min_cohort_coverage = if (nrow(cov)) min(cov$coverage) else 0,
    passes_gate = nrow(cov) == length(composition_required_cohorts) &&
      all(cov$coverage >= 1 - 1e-12)
  )
}))

eligibility <- rbindlist(lapply(names(methods), function(label) {
  dyn <- dynamic_table[method == label]
  if (!nrow(dyn)) return(data.table())
  dyn[, .(
    method = label,
    n_cohorts = .N,
    min_p05 = min(p05),
    passes_dynamic_range = all(p05 >= composition_dynamic_range_floor)
  ), by = celltype]
}))

criteria <- rbindlist(lapply(names(methods), function(label) {
  lin <- lineage_table[method == label]
  dir <- direction_table[method == label]
  elig <- eligibility[method == label & celltype %in% COVARIATE_LINEAGES & passes_dynamic_range]
  hep <- lin[celltype == DENOMINATOR_LINEAGE, median_fraction]
  coll <- if (nrow(elig) >= 2L)
    logratio_collinearity(label, audits[[label]]$matched, elig$celltype) else data.table()
  data.table(
    method = label,
    passes_gate = gate[method == label, passes_gate],
    zero_inflation_max = if (nrow(lin)) max(lin$zero_fraction) else NA_real_,
    passes_zero_inflation = nrow(lin) > 0L &&
      max(lin$zero_fraction) <= composition_max_zero_fraction,
    median_hepatocyte = if (length(hep) == 1L) hep else NA_real_,
    passes_plausibility = length(hep) == 1L && is.finite(hep) &&
      hep >= composition_hepatocyte_median_range[[1L]] &&
      hep <= composition_hepatocyte_median_range[[2L]],
    n_direction_lineages_ok = if (nrow(dir))
      dir[, .(ok = sum(agrees) >= composition_min_direction_cohorts), by = celltype][, sum(ok)]
      else 0L,
    passes_direction = nrow(dir) > 0L &&
      dir[, .(ok = sum(agrees) >= composition_min_direction_cohorts), by = celltype][, all(ok)],
    n_eligible_covariate_lineages = nrow(elig),
    eligible_covariate_lineages = paste(elig$celltype, collapse = ";"),
    max_logratio_correlation = if (nrow(coll)) max(coll$max_abs_correlation, na.rm = TRUE) else NA_real_,
    all_cohorts_full_rank = if (nrow(coll)) all(coll$full_rank) else NA
  )
}))
criteria[, eligible := passes_gate & passes_zero_inflation & passes_plausibility &
           passes_direction]

collinearity_table <- rbindlist(lapply(names(methods), function(label) {
  elig <- eligibility[method == label & celltype %in% COVARIATE_LINEAGES & passes_dynamic_range]
  if (nrow(elig) < 2L) return(data.table())
  logratio_collinearity(label, audits[[label]]$matched, elig$celltype)
}), use.names = TRUE, fill = TRUE)

qualified <- criteria[eligible == TRUE]
accepted_method <- if (nrow(qualified)) {
  # Decision rule, fixed: among qualified estimators prefer the one resolving the
  # most covariate lineages, then the lowest zero inflation, then alphabetical.
  setorder(qualified, -n_eligible_covariate_lineages, zero_inflation_max, method)
  qualified$method[[1L]]
} else NA_character_

message("[5/6] Publishing acceptance audit",
        if (is.na(accepted_method)) " (no estimator qualified)" else
          paste0(" (accepted: ", accepted_method, ")"))

tmp <- atomic_dir(composition_acceptance_root)
on.exit(if (dir.exists(tmp)) unlink(tmp, recursive = TRUE), add = TRUE)

write_tsv(coverage_table, file.path(tmp, "cohort_coverage.tsv"))
write_tsv(lineage_table, file.path(tmp, "lineage_audit.tsv"))
write_tsv(dynamic_table, file.path(tmp, "dynamic_range_audit.tsv"))
write_tsv(direction_table, file.path(tmp, "stage_direction_audit.tsv"))
write_tsv(eligibility, file.path(tmp, "lineage_eligibility.tsv"))
write_tsv(collinearity_table, file.path(tmp, "logratio_collinearity.tsv"))
write_tsv(gate, file.path(tmp, "coverage_gate.tsv"))
write_tsv(criteria, file.path(tmp, "method_criteria.tsv"))

rule <- data.table(
  criterion = c("gate_coverage", "zero_inflation", "dynamic_range", "plausibility",
                "stage_direction", "collinearity", "decision"),
  rule = c(
    "complete proportions for every discovery cohort and the holdout; applied before any quality criterion",
    sprintf("max exact-zero fraction over all lineages <= %.2f", composition_max_zero_fraction),
    sprintf("a lineage may enter the model only if its 5th-percentile proportion is >= %.0e in every required cohort",
            composition_dynamic_range_floor),
    sprintf("median hepatocyte fraction within [%.2f, %.2f]",
            composition_hepatocyte_median_range[[1L]], composition_hepatocyte_median_range[[2L]]),
    sprintf("expected Kleiner-stage sign in at least %d of four discovery cohorts for hepatocytes, fibroblasts and cholangiocytes",
            composition_min_direction_cohorts),
    sprintf("max absolute correlation among eligible log-ratios < %.2f, design full rank in every discovery cohort",
            composition_max_logratio_correlation),
    "among qualified estimators: most eligible covariate lineages, then lowest zero inflation, then alphabetical"
  ),
  uses_outcome = c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE)
)
write_tsv(rule, file.path(tmp, "acceptance_rule.tsv"))

if (!is.na(accepted_method)) {
  accepted <- copy(audits[[accepted_method]]$matched)
  accepted[, .row_complete := NULL]
  fwrite(accepted, file.path(tmp, "accepted_composition.tsv.gz"),
         sep = "\t", quote = FALSE, na = "NA")
}

inputs <- c(input_dge, input_meta, composition_candidate, composition_testability)
write_tsv(data.table(
  input = c("corrected_dge", "corrected_meta", "music_deposited", "composition_testability"),
  path = inputs, sha256 = vapply(inputs, sha256_file, character(1))
), file.path(tmp, "input_manifest.tsv"))
writeLines(capture.output(sessionInfo()), file.path(tmp, "sessionInfo.txt"))

eligible_lineages <- if (!is.na(accepted_method)) {
  eligibility[method == accepted_method & celltype %in% COVARIATE_LINEAGES &
                passes_dynamic_range, celltype]
} else character(0)

payload <- list(
  state = if (is.na(accepted_method)) "NO_ESTIMATOR_QUALIFIED" else "COMPOSITION_ESTIMATOR_ACCEPTED",
  release_id = release_id,
  workstream_id = workstream_id,
  accepted_method = accepted_method,
  eligible_covariate_lineages = eligible_lineages,
  n_eligible_covariate_lineages = length(eligible_lineages),
  logratio_specification_available = length(eligible_lineages) >= 2L,
  clr_component_specification = paste(
    "The prespecified parallel adjustment is the first",
    composition_n_clr_components,
    "principal components of the within-cohort CLR matrix over all retained lineages.",
    "It is well conditioned regardless of which individual lineage underflows and is",
    "the sole composition adjustment when fewer than two lineages clear the dynamic-range floor."
  ),
  n_census_samples = nrow(census),
  criteria_are_blind = FALSE,
  criteria_provenance = paste(
    "Zero-fraction, coverage and stage-direction criteria were evaluated during planning",
    "on 2026-08-12 and are not blind. Thresholds sit far from every observed value.",
    "The stage-direction criterion uses Kleiner stage and constrains only the composition",
    "estimator and lineage set, never program selection."
  ),
  supersedes = "BULK-PROGRAM-MAP-v8 used the deposited MuSiC estimates with no real-bulk acceptance test"
)
jsonlite::write_json(payload, file.path(tmp, "COMPOSITION_ACCEPTED.json"),
                     pretty = TRUE, auto_unbox = TRUE)

artifacts <- setdiff(list.files(tmp), "artifact_manifest.tsv")
write_tsv(data.table(
  artifact = artifacts,
  sha256 = vapply(file.path(tmp, artifacts), sha256_file, character(1)),
  size_bytes = file.info(file.path(tmp, artifacts))$size
), file.path(tmp, "artifact_manifest.tsv"))

message("[6/6] Sealing")
publish_dir(tmp, composition_acceptance_root)
Sys.chmod(list.files(composition_acceptance_root, full.names = TRUE), mode = "0440")
cat("COMPOSITION_ACCEPTANCE_COMPLETE\t", accepted_method, "\t",
    length(eligible_lineages), " eligible covariate lineages\t",
    composition_acceptance_root, "\n", sep = "")
