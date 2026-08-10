#!/usr/bin/env Rscript
# Figure-independent, outcome-blind Hotspot v2 donor refit.
#
# Reuses the existing 117 Hotspot module definitions and scores. The primary
# score is the pooled-cell donor mean; the unweighted mean across a donor's runs
# is a scoring sensitivity. No external spatial, genetics, protein, or screen
# outcome is read. All generated files are isolated under the candidate root.

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
})

options(stringsAsFactors = FALSE)

BASE <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
RELEASE_ID <- "program-context-v2-candidate-2026-08-07"
release_override <- Sys.getenv("HOTSPOT_V2_RELEASE_ID", unset = "")
if (nzchar(release_override) && release_override != RELEASE_ID) {
  stop("HOTSPOT_V2_RELEASE_ID must equal the ratified ID: ", RELEASE_ID)
}

OUT <- file.path(
  BASE, "Analysis/Multimodal_Program_Projection/candidates",
  RELEASE_ID, "hotspot"
)
WORK <- file.path(OUT, "work")
CHECK_ONLY <- "--check-only" %in% commandArgs(trailingOnly = TRUE) ||
  tolower(Sys.getenv("HOTSPOT_V2_CHECK_ONLY", unset = "false")) %in%
    c("true", "t", "1", "yes", "y")

HS <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/hotspot_modules")
SCRIPT_DIR <- file.path(BASE, "Analysis/SingleCell/scripts/hotspot_modules")

PRIMARY_SCORE_FILE <- file.path(HS, "donor_collapse/donor_scores_all_weighted.tsv")
EQUAL_RUN_FILE <- file.path(HS, "donor_scores_all.tsv")
LOO_FILE <- file.path(HS, "loo_stability.tsv")
NAME_FILE <- file.path(HS, "module_names.tsv")
META_FILE <- file.path(
  BASE,
  "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/donor_metadata_extended.tsv"
)
FSTAGE_FILE <- file.path(
  BASE,
  "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/donor_fstage_documented.tsv"
)
DONOR_LIB <- file.path(BASE, "Analysis/SingleCell/scripts/lib_donor_collapse.R")
V1_ROOT <- file.path(BASE, "Analysis/Multimodal_Program_Projection/results")
GENCODE_FILE <- file.path(BASE, "data/gencode_v49_gene_metadata.tsv.gz")
V1_ANCHOR_SHA256 <- c(
  "freeze_manifest.tsv" = "fa32594484a1c7776e5923259876cf038e33a1d287ebd2474f3badc27a92eda7",
  "frozen_programs.tsv" = "79fe355dc36a764f67a04171947de9dc1df867b2260a66d6319fc66756ea2e69",
  "frozen_program_membership.tsv" = "900b6e7bd514994fef2dea2c976f5cd0eb715c715a859cb86ced4bdcaf7f5a4c"
)

CELL_TYPES <- c(
  "hepatocytes", "fibroblasts", "macrophages", "cholangiocytes", "tcells"
)
EXPECTED_CENSUS <- c(
  hepatocytes = 30L, fibroblasts = 29L, macrophages = 19L,
  cholangiocytes = 28L, tcells = 11L
)
EXPECTED_PROGRAMS <- sum(EXPECTED_CENSUS)
PRIMARY_STAGES <- c("Healthy", "Steatosis", "Steatohepatitis")

source(DONOR_LIB)

fail <- function(...) stop(..., call. = FALSE)

assert_true <- function(ok, message) {
  if (!isTRUE(ok)) fail(message)
}

sha256_file <- function(path) {
  digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
}

rel_path <- function(path) {
  root <- paste0(normalizePath(BASE, winslash = "/", mustWork = TRUE), "/")
  resolved <- normalizePath(path, winslash = "/", mustWork = TRUE)
  if (startsWith(resolved, root)) substring(resolved, nchar(root) + 1L) else resolved
}

manifest_for <- function(paths, role) {
  paths <- unique(paths)
  missing <- paths[!file.exists(paths)]
  if (length(missing)) fail("Missing required file(s): ", paste(missing, collapse = ", "))
  info <- file.info(paths)
  data.table(
    role = role,
    relative_path = vapply(paths, rel_path, character(1)),
    bytes = as.numeric(info$size),
    sha256 = vapply(paths, sha256_file, character(1)),
    modified_utc = format(info$mtime, tz = "UTC", usetz = TRUE),
    row_count = NA_integer_,
    gene_count = NA_integer_,
    count_status = "not_measured"
  )
}

add_tabular_counts <- function(manifest, paths) {
  metrics <- rbindlist(lapply(paths, function(path) {
    relative_path <- rel_path(path)
    if (!grepl("\\.(csv|tsv)(\\.gz)?$", path, ignore.case = TRUE)) {
      return(data.table(
        relative_path = relative_path, row_count = NA_integer_,
        gene_count = NA_integer_, count_status = "not_applicable"
      ))
    }
    x <- tryCatch(
      fread(path, showProgress = FALSE),
      error = function(e) fail("Could not count tabular input ", path, ": ", e$message)
    )
    gene_col <- intersect(
      c("gene", "source_gene", "gene_id", "ensembl_base", "gene_name", "symbol"),
      names(x)
    )
    gene_count <- NA_integer_
    if (length(gene_col)) {
      values <- as.character(x[[gene_col[1L]]])
      gene_count <- uniqueN(values[!is.na(values) & nzchar(values)])
    }
    data.table(
      relative_path = relative_path, row_count = nrow(x),
      gene_count = gene_count, count_status = "measured"
    )
  }))
  manifest[metrics, on = "relative_path", `:=`(
    row_count = i.row_count,
    gene_count = i.gene_count,
    count_status = i.count_status
  )]
  assert_true(!any(manifest$count_status == "not_measured"),
              "Input count audit left an unclassified file")
  manifest[]
}

write_new_tsv <- function(x, path) {
  if (file.exists(path)) fail("Refusing to overwrite candidate work product: ", path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp.", Sys.getpid())
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  fwrite(x, tmp, sep = "\t", quote = FALSE, na = "NA")
  if (!file.rename(tmp, path)) fail("Atomic rename failed for ", path)
}

parse_flag <- function(x) {
  if (is.logical(x)) return(fifelse(is.na(x), FALSE, x))
  y <- tolower(trimws(as.character(x)))
  y %in% c("true", "t", "1", "yes", "y")
}

first_nonmissing <- function(x) {
  keep <- !is.na(x) & trimws(as.character(x)) != ""
  if (!any(keep)) return(NA_character_)
  as.character(x[which(keep)[1L]])
}

direction_of <- function(x) {
  fifelse(is.na(x), NA_character_, fifelse(x > 0, "positive", fifelse(x < 0, "negative", "zero")))
}

normalize_score_modules <- function(x, source_name) {
  required <- c("cell_type", "module", "score")
  if (!all(required %in% names(x))) {
    fail(source_name, " lacks columns: ", paste(setdiff(required, names(x)), collapse = ", "))
  }
  x <- copy(x)
  x[, module_source := as.character(module)]
  prefixed <- grepl("__", x$module_source, fixed = TRUE)
  prefix <- sub("__.*$", "", x$module_source)
  if (any(prefixed & prefix != x$cell_type)) {
    fail(source_name, " has module prefixes inconsistent with cell_type")
  }
  x[, module := suppressWarnings(as.integer(sub("^.*__", "", module_source)))]
  if (anyNA(x$module) || any(!is.finite(x$score))) {
    fail(source_name, " contains nonnumeric module IDs or nonfinite scores")
  }
  x[, module_source := NULL]
  x[]
}

load_membership <- function() {
  rbindlist(lapply(CELL_TYPES, function(ct) {
    path <- file.path(HS, ct, "module_genes.tsv")
    if (!file.exists(path)) fail("Missing module membership: ", path)
    x <- fread(path)
    if (!all(c("gene", "module", "weight") %in% names(x))) {
      fail("Malformed module membership: ", path)
    }
    x[, `:=`(
      cell_type = ct,
      module = as.integer(module),
      gene = as.character(gene),
      weight = as.numeric(weight)
    )]
    if (anyNA(x$module) || any(!nzchar(x$gene)) ||
        any(!is.finite(x$weight)) || any(x$weight <= 0)) {
      fail("Invalid module/gene/positive weight in ", path)
    }
    x[, .(cell_type, module, gene, weight)]
  }), use.names = TRUE)
}

validate_program_universe <- function(membership, primary_scores, equal_scores,
                                      stability, names_dt) {
  programs <- unique(membership[, .(cell_type, module)])
  setorder(programs, cell_type, module)
  census <- programs[, .N, by = cell_type]
  observed <- setNames(census$N, census$cell_type)
  assert_true(setequal(names(observed), names(EXPECTED_CENSUS)),
              "Unexpected Hotspot lineage set")
  assert_true(all(observed[names(EXPECTED_CENSUS)] == EXPECTED_CENSUS),
              paste0("Hotspot lineage census drift: ",
                     paste(names(observed), observed, sep = "=", collapse = ", ")))
  assert_true(nrow(programs) == EXPECTED_PROGRAMS,
              "Hotspot program universe is not exactly 117")

  key <- function(z) unique(z[, .(cell_type, module)])
  assert_true(fsetequal(programs, key(primary_scores)),
              "Pooled-cell score program universe differs from memberships")
  assert_true(fsetequal(programs, key(equal_scores)),
              "Equal-run score program universe differs from memberships")
  assert_true(fsetequal(programs, key(stability)),
              "LOO stability program universe differs from memberships")
  assert_true(fsetequal(programs, key(names_dt)),
              "Module-name program universe differs from memberships")
  programs
}

build_donor_metadata <- function() {
  required <- c("sample", "dataset", "disease_stage_coarse", "exclude_stage_analysis")
  header <- names(fread(META_FILE, nrows = 0L))
  if (!all(required %in% header)) {
    fail("Extended donor metadata lacks allowed columns: ",
         paste(setdiff(required, header), collapse = ", "))
  }

  # Deliberately select only the four approved columns. Inferred/augmented
  # fibrosis fields present in the source file never enter this process.
  x <- fread(META_FILE, select = required)
  x[, exclude_stage_analysis := parse_flag(exclude_stage_analysis)]
  donor_map <- build_srr_to_donor_map(BASE)
  x[, donor := fifelse(sample %in% names(donor_map), donor_map[sample], sample)]

  audit <- x[, .(
    n_source_records = .N,
    n_dataset_values = uniqueN(dataset[!is.na(dataset) & dataset != ""]),
    n_stage_values = uniqueN(disease_stage_coarse[
      !is.na(disease_stage_coarse) & disease_stage_coarse != ""
    ]),
    any_excluded = any(exclude_stage_analysis)
  ), by = donor]
  bad <- audit[n_dataset_values > 1L | n_stage_values > 1L]
  if (nrow(bad)) {
    fail("Dataset/stage disagreement within biological donor(s): ",
         paste(bad$donor, collapse = ", "))
  }

  meta <- x[, .(
    dataset = first_nonmissing(dataset),
    disease_stage_coarse = first_nonmissing(disease_stage_coarse),
    exclude_stage_analysis = any(exclude_stage_analysis),
    n_metadata_source_records = .N
  ), by = donor]
  meta[, stage_ordinal := fcase(
    disease_stage_coarse == "Healthy", 0,
    disease_stage_coarse == "Steatosis", 1,
    disease_stage_coarse == "Steatohepatitis", 2,
    default = NA_real_
  )]
  list(meta = meta, audit = audit, donor_map = donor_map)
}

fit_fixed_ols <- function(d, stage_col, dataset_adjust = TRUE) {
  if (!nrow(d)) {
    return(data.table(
      estimable = FALSE, failure_reason = "no_complete_rows",
      beta = NA_real_, se = NA_real_, statistic = NA_real_, pvalue = NA_real_,
      hc3_se = NA_real_, hc3_statistic = NA_real_, hc3_pvalue = NA_real_,
      hc3_failure_reason = "no_complete_rows", max_leverage = NA_real_,
      residual_df = NA_integer_, n_donors = 0L, n_datasets = 0L,
      n_stage_levels = 0L, n_stage0 = 0L, n_stage1 = 0L, n_stage2 = 0L,
      n_stage3 = 0L, n_stage4 = 0L, design_rank = 0L,
      n_model_columns = 0L, condition_number = NA_real_,
      stage_residual_variance = NA_real_
    ))
  }

  stage <- suppressWarnings(as.numeric(d[[stage_col]]))
  score <- suppressWarnings(as.numeric(d$score))
  dataset <- as.character(d$dataset)
  keep <- is.finite(stage) & is.finite(score) & !is.na(dataset) & nzchar(dataset)
  d <- d[keep]
  stage <- stage[keep]
  score <- score[keep]
  dataset <- dataset[keep]

  stage_counts <- tabulate(match(stage, 0:4), nbins = 5L)
  n <- length(score)
  n_ds <- uniqueN(dataset)
  n_stage <- uniqueN(stage)

  base_result <- function(reason, rank = NA_integer_, pcols = NA_integer_,
                          cond = NA_real_, stage_var = NA_real_, rdf = NA_integer_) {
    data.table(
      estimable = FALSE, failure_reason = reason,
      beta = NA_real_, se = NA_real_, statistic = NA_real_, pvalue = NA_real_,
      hc3_se = NA_real_, hc3_statistic = NA_real_, hc3_pvalue = NA_real_,
      hc3_failure_reason = reason, max_leverage = NA_real_,
      residual_df = rdf, n_donors = n, n_datasets = n_ds,
      n_stage_levels = n_stage,
      n_stage0 = stage_counts[1L], n_stage1 = stage_counts[2L],
      n_stage2 = stage_counts[3L], n_stage3 = stage_counts[4L],
      n_stage4 = stage_counts[5L], design_rank = rank,
      n_model_columns = pcols, condition_number = cond,
      stage_residual_variance = stage_var
    )
  }

  if (n == 0L) return(base_result("no_complete_rows"))
  if (n_stage < 2L) return(base_result("fewer_than_two_stage_levels"))

  ds_levels <- sort(unique(dataset))
  dummies <- matrix(numeric(0), nrow = n, ncol = 0L)
  if (dataset_adjust && length(ds_levels) > 1L) {
    dummies <- do.call(cbind, lapply(ds_levels[-1L], function(level) {
      as.numeric(dataset == level)
    }))
    if (is.null(dim(dummies))) dummies <- matrix(dummies, ncol = 1L)
    colnames(dummies) <- paste0("dataset_", make.names(ds_levels[-1L], unique = TRUE))
  }
  z <- cbind(`(Intercept)` = 1, dummies)
  x <- cbind(`(Intercept)` = 1, stage_effect = stage, dummies)
  qrx <- qr(x)
  design_rank <- qrx$rank
  pcols <- ncol(x)
  cond <- tryCatch(kappa(x, exact = TRUE), error = function(e) Inf)
  stage_resid <- tryCatch(lm.fit(z, stage)$residuals, error = function(e) rep(NA_real_, n))
  stage_var <- if (sum(is.finite(stage_resid)) > 1L) var(stage_resid, na.rm = TRUE) else NA_real_
  rdf <- n - design_rank

  if (!is.finite(stage_var) || stage_var <= sqrt(.Machine$double.eps)) {
    return(base_result("zero_stage_variance_after_dataset", design_rank, pcols,
                       cond, stage_var, rdf))
  }
  if (design_rank < pcols) {
    return(base_result("rank_deficient_design", design_rank, pcols, cond, stage_var, rdf))
  }
  if (rdf <= 0L) {
    return(base_result("nonpositive_residual_df", design_rank, pcols, cond, stage_var, rdf))
  }

  fit <- lm.fit(x, score)
  sigma2 <- sum(fit$residuals^2) / rdf
  inv_xtx <- tryCatch({
    r <- qr.R(fit$qr)[seq_len(pcols), seq_len(pcols), drop = FALSE]
    inv_pivoted <- chol2inv(r)
    ans <- matrix(0, nrow = pcols, ncol = pcols,
                  dimnames = list(colnames(x), colnames(x)))
    pivot <- fit$qr$pivot[seq_len(pcols)]
    ans[pivot, pivot] <- inv_pivoted
    ans
  }, error = function(e) NULL)
  if (is.null(inv_xtx) || !is.finite(sigma2) || sigma2 < 0) {
    return(base_result("nonfinite_covariance", design_rank, pcols, cond, stage_var, rdf))
  }
  beta <- unname(fit$coefficients["stage_effect"])
  se <- sqrt(unname(sigma2 * inv_xtx["stage_effect", "stage_effect"]))
  if (!is.finite(beta) || !is.finite(se) || se <= 0) {
    return(base_result("nonfinite_stage_coefficient", design_rank, pcols, cond,
                       stage_var, rdf))
  }
  statistic <- beta / se
  pvalue <- 2 * pt(abs(statistic), df = rdf, lower.tail = FALSE)

  leverage <- rowSums((x %*% inv_xtx) * x)
  max_leverage <- max(leverage)
  hc3_failure_reason <- NA_character_
  hc3_se <- hc3_statistic <- hc3_pvalue <- NA_real_
  if (any(!is.finite(leverage)) || any(1 - leverage <= sqrt(.Machine$double.eps))) {
    hc3_failure_reason <- "nonfinite_or_unit_leverage"
  } else {
    adjusted_sq_residual <- (fit$residuals / (1 - leverage))^2
    meat <- crossprod(x, sweep(x, 1L, adjusted_sq_residual, `*`))
    hc3_vcov <- inv_xtx %*% meat %*% inv_xtx
    hc3_se <- sqrt(unname(hc3_vcov["stage_effect", "stage_effect"]))
    if (!is.finite(hc3_se) || hc3_se <= 0) {
      hc3_failure_reason <- "nonfinite_hc3_stage_uncertainty"
      hc3_se <- NA_real_
    } else {
      hc3_statistic <- beta / hc3_se
      hc3_pvalue <- 2 * pt(abs(hc3_statistic), df = rdf, lower.tail = FALSE)
    }
  }

  data.table(
    estimable = TRUE, failure_reason = NA_character_, beta = beta, se = se,
    statistic = statistic, pvalue = pvalue,
    hc3_se = hc3_se, hc3_statistic = hc3_statistic, hc3_pvalue = hc3_pvalue,
    hc3_failure_reason = hc3_failure_reason, max_leverage = max_leverage,
    residual_df = rdf,
    n_donors = n, n_datasets = n_ds, n_stage_levels = n_stage,
    n_stage0 = stage_counts[1L], n_stage1 = stage_counts[2L],
    n_stage2 = stage_counts[3L], n_stage3 = stage_counts[4L],
    n_stage4 = stage_counts[5L], design_rank = design_rank,
    n_model_columns = pcols, condition_number = cond,
    stage_residual_variance = stage_var
  )
}

fit_program_family <- function(scores, programs, stage_col,
                               dataset_adjust = TRUE, adjust_bh = FALSE) {
  ans <- rbindlist(lapply(seq_len(nrow(programs)), function(i) {
    ct <- programs$cell_type[i]
    mod <- programs$module[i]
    d <- scores[cell_type == ct & module == mod]
    cbind(
      data.table(cell_type = ct, module = mod),
      fit_fixed_ols(d, stage_col = stage_col, dataset_adjust = dataset_adjust)
    )
  }), use.names = TRUE, fill = TRUE)
  setorder(ans, cell_type, module)
  ans[, direction := direction_of(beta)]
  ans[, `:=`(qvalue = NA_real_, hc3_qvalue = NA_real_)]
  if (adjust_bh) {
    idx <- which(is.finite(ans$pvalue))
    ans$qvalue[idx] <- p.adjust(ans$pvalue[idx], method = "BH", n = EXPECTED_PROGRAMS)
    hc3_idx <- which(is.finite(ans$hc3_pvalue))
    ans$hc3_qvalue[hc3_idx] <- p.adjust(
      ans$hc3_pvalue[hc3_idx], method = "BH", n = EXPECTED_PROGRAMS
    )
  }
  ans[]
}

fit_context_effects <- function(primary_scores, programs, datasets) {
  rbindlist(lapply(seq_len(nrow(programs)), function(i) {
    ct <- programs$cell_type[i]
    mod <- programs$module[i]
    d <- primary_scores[cell_type == ct & module == mod]
    rbindlist(lapply(datasets, function(ds) {
      cohort <- cbind(
        data.table(
          cell_type = ct, module = mod, analysis_type = "cohort",
          dataset_context = ds, held_out_dataset = NA_character_
        ),
        fit_fixed_ols(d[dataset == ds], "stage_ordinal", dataset_adjust = FALSE)
      )
      lodo <- cbind(
        data.table(
          cell_type = ct, module = mod, analysis_type = "lodo",
          dataset_context = NA_character_, held_out_dataset = ds
        ),
        fit_fixed_ols(d[dataset != ds], "stage_ordinal", dataset_adjust = TRUE)
      )
      rbind(cohort, lodo, use.names = TRUE, fill = TRUE)
    }), use.names = TRUE, fill = TRUE)
  }), use.names = TRUE, fill = TRUE)
}

classify_source_stability <- function(primary_effects, context_effects) {
  rbindlist(lapply(seq_len(nrow(primary_effects)), function(i) {
    p <- primary_effects[i]
    cfx <- context_effects[cell_type == p$cell_type & module == p$module]
    reason <- NA_character_
    label <- "not_assessable"
    n_cohort_same <- 0L
    n_lodo_same <- 0L
    if (isTRUE(p$estimable) && is.finite(p$beta) && p$beta != 0) {
      expected_sign <- sign(p$beta)
      gse <- cfx[analysis_type == "lodo" & held_out_dataset == "GSE244832"]
      cohort_ok <- cfx[analysis_type == "cohort" & estimable == TRUE & is.finite(beta)]
      lodo_ok <- cfx[analysis_type == "lodo" & estimable == TRUE & is.finite(beta)]
      n_cohort_same <- cohort_ok[sign(beta) == expected_sign, .N]
      n_lodo_same <- lodo_ok[sign(beta) == expected_sign, .N]
      other_opposite <- rbind(cohort_ok, lodo_ok, use.names = TRUE, fill = TRUE)[
        beta != 0 & sign(beta) != expected_sign, .N
      ]
      if (nrow(gse) != 1L || !isTRUE(gse$estimable[1L])) {
        label <- "GSE244832_dependent"
        reason <- "GSE244832-excluded model non-estimable"
      } else if (!is.finite(gse$beta[1L]) || sign(gse$beta[1L]) != expected_sign) {
        label <- "GSE244832_dependent"
        reason <- "GSE244832-excluded direction differs from primary"
      } else if (other_opposite > 0L) {
        label <- "direction_heterogeneous"
        reason <- "at least one estimable cohort/LODO direction differs"
      } else if (n_cohort_same >= 2L) {
        label <- "multi_dataset_supported"
        reason <- "direction retained in at least two independently estimable cohort slopes"
      } else {
        reason <- paste0(
          "fewer than two independently estimable cohort slopes; LODO fits are ",
          "reported but not counted as independent support"
        )
      }
    } else {
      reason <- "primary model non-estimable or zero"
    }
    data.table(
      cell_type = p$cell_type, module = p$module,
      source_stability = label, source_stability_reason = reason,
      n_cohort_same_direction = n_cohort_same,
      n_lodo_same_direction = n_lodo_same
    )
  }))
}

# -------------------------------------------------------------------------
# Load only prespecified, internal Hotspot inputs.
# -------------------------------------------------------------------------
membership <- load_membership()
stability <- fread(LOO_FILE)
names_dt <- fread(NAME_FILE)
stability[, module := as.integer(module)]
names_dt[, module := as.integer(module)]
if (anyDuplicated(stability[, .(cell_type, module)]) ||
    anyDuplicated(names_dt[, .(cell_type, module)])) {
  fail("Duplicate program key in stability or module-name input")
}

primary_raw <- normalize_score_modules(fread(PRIMARY_SCORE_FILE), PRIMARY_SCORE_FILE)
if (!"sample" %in% names(primary_raw)) fail("Primary pooled-cell file lacks sample")
setnames(primary_raw, "sample", "donor")
primary_raw <- primary_raw[, .(cell_type, module, donor = as.character(donor), score)]
if (anyDuplicated(primary_raw[, .(cell_type, module, donor)])) {
  fail("Duplicate pooled-cell donor/program score")
}

equal_raw <- normalize_score_modules(fread(EQUAL_RUN_FILE), EQUAL_RUN_FILE)
if (!"sample" %in% names(equal_raw)) fail("Equal-run source lacks sample")
meta_obj <- build_donor_metadata()
donor_map <- meta_obj$donor_map
equal_raw[, donor := fifelse(sample %in% names(donor_map), donor_map[sample], sample)]
equal_donor <- equal_raw[, .(
  score = if (all(is.na(score))) NA_real_ else mean(score, na.rm = TRUE)
), by = .(cell_type, module, donor)]

programs <- validate_program_universe(
  membership, primary_raw, equal_donor, stability, names_dt
)

meta_donor <- meta_obj$meta
assert_joinable <- function(scores, label) {
  missing <- setdiff(unique(scores$donor), meta_donor$donor)
  if (length(missing)) {
    fail(label, " donor(s) missing metadata: ", paste(head(missing, 20L), collapse = ", "))
  }
}
assert_joinable(primary_raw, "Pooled-cell")
assert_joinable(equal_donor, "Equal-run")

eligible_meta <- meta_donor[
  !exclude_stage_analysis & disease_stage_coarse %in% PRIMARY_STAGES &
    is.finite(stage_ordinal)
]
primary_scores <- merge(primary_raw, eligible_meta, by = "donor", all = FALSE)
equal_scores <- merge(equal_donor, eligible_meta, by = "donor", all = FALSE)
setcolorder(primary_scores, c(
  "cell_type", "module", "donor", "score", "dataset",
  "disease_stage_coarse", "stage_ordinal", "exclude_stage_analysis",
  "n_metadata_source_records"
))
setcolorder(equal_scores, names(primary_scores))

if (anyDuplicated(primary_scores[, .(cell_type, module, donor)]) ||
    anyDuplicated(equal_scores[, .(cell_type, module, donor)])) {
  fail("More than one score per biological donor/program after collapse")
}
assert_true(
  fsetequal(
    primary_scores[, .(cell_type, module, donor)],
    equal_scores[, .(cell_type, module, donor)]
  ),
  "Pooled-cell and equal-run sensitivities do not contain identical donor/program keys"
)

primary_effects <- fit_program_family(
  primary_scores, programs, "stage_ordinal", dataset_adjust = TRUE, adjust_bh = TRUE
)
equal_effects <- fit_program_family(
  equal_scores, programs, "stage_ordinal", dataset_adjust = TRUE, adjust_bh = TRUE
)

loo_columns <- grep("^loo_", names(stability), value = TRUE)
assert_true(length(loo_columns) >= 5L,
            "LOO input lacks per-dataset Jaccard columns")
stability[, stability_median := apply(.SD, 1L, function(z) {
  z <- z[is.finite(z)]
  if (length(z)) median(z) else NA_real_
}), .SDcols = loo_columns]
preflight_selection <- Reduce(
  function(x, y) merge(x, y, by = c("cell_type", "module"), all = TRUE),
  list(
    primary_effects[, .(
      cell_type, module, primary_estimable = estimable,
      primary_qvalue = qvalue, primary_hc3_qvalue = hc3_qvalue,
      primary_beta = beta
    )],
    equal_effects[, .(
      cell_type, module, equal_estimable = estimable, equal_beta = beta
    )],
    stability[, .(cell_type, module, stability_median)]
  )
)
preflight_selection[, primary_selected :=
  primary_estimable & is.finite(primary_qvalue) & primary_qvalue < 0.05 &
    is.finite(stability_median) & stability_median >= 0.5]
preflight_selection[, direction_agree :=
  primary_estimable & equal_estimable & is.finite(primary_beta) &
    is.finite(equal_beta) & primary_beta != 0 & equal_beta != 0 &
    sign(primary_beta) == sign(equal_beta)]
preflight_selection[, hc3_supported :=
  primary_estimable & is.finite(primary_hc3_qvalue) & primary_hc3_qvalue < 0.05]
preflight_selection[, robust_display :=
  primary_selected & direction_agree & hc3_supported]

datasets <- sort(unique(primary_scores$dataset))
context_effects <- fit_context_effects(primary_scores, programs, datasets)
context_effects[, direction := direction_of(beta)]
source_stability <- classify_source_stability(primary_effects, context_effects)

stage_counts <- primary_scores[, .(stage_dataset_n_donors = uniqueN(donor)), by = .(
  cell_type, module, dataset, disease_stage_coarse, stage_ordinal
)]
program_grid <- copy(programs)[, join_key := 1L]
design_grid <- CJ(
  dataset = datasets,
  disease_stage_coarse = PRIMARY_STAGES,
  unique = TRUE
)
design_grid[, `:=`(
  stage_ordinal = match(disease_stage_coarse, PRIMARY_STAGES) - 1L,
  join_key = 1L
)]
stage_design <- merge(program_grid, design_grid, by = "join_key", allow.cartesian = TRUE)
stage_design[, join_key := NULL]
stage_design <- merge(
  stage_design, stage_counts,
  by = c("cell_type", "module", "dataset", "disease_stage_coarse", "stage_ordinal"),
  all.x = TRUE
)
stage_design[is.na(stage_dataset_n_donors), stage_dataset_n_donors := 0L]
primary_design_diagnostics <- primary_effects[, .(
  cell_type, module,
  primary_estimable = estimable,
  primary_failure_reason = failure_reason,
  primary_residual_df = residual_df,
  primary_n_donors = n_donors,
  primary_n_datasets = n_datasets,
  primary_n_stage_levels = n_stage_levels,
  primary_design_rank = design_rank,
  primary_n_model_columns = n_model_columns,
  primary_condition_number = condition_number,
  primary_stage_residual_variance = stage_residual_variance
)]
stage_design <- merge(
  stage_design, primary_design_diagnostics,
  by = c("cell_type", "module"), all.x = TRUE
)
setorder(stage_design, cell_type, module, dataset, stage_ordinal)

# Documented F-stage sensitivity: 58 source records collapse to 46 true donors.
fstage <- fread(FSTAGE_FILE, select = c("sample", "dataset", "F_stage_documented", "source"))
fstage <- fstage[is.finite(F_stage_documented)]
n_fstage_records <- nrow(fstage)
fstage[, donor := fifelse(sample %in% names(donor_map), donor_map[sample], sample)]
fstage_conflict <- fstage[, .(n_stage_values = uniqueN(F_stage_documented)), by = donor][
  n_stage_values > 1L
]
if (nrow(fstage_conflict)) {
  fail("Documented F-stage conflict within donor(s): ",
       paste(fstage_conflict$donor, collapse = ", "))
}
fstage_donor <- fstage[, .(
  dataset = first_nonmissing(dataset),
  F_stage_documented = as.numeric(F_stage_documented[1L]),
  fstage_source = paste(sort(unique(source)), collapse = ";"),
  n_documented_source_records = .N
), by = donor]
n_fstage_donors <- uniqueN(fstage_donor$donor)
if (n_fstage_records != 58L || n_fstage_donors != 46L) {
  fail("Documented F-stage contract drift: expected 58 records -> 46 donors; observed ",
       n_fstage_records, " -> ", n_fstage_donors)
}
fstage_scores <- merge(primary_raw, fstage_donor, by = "donor", all = FALSE)
fstage_effects <- fit_program_family(
  fstage_scores, programs, "F_stage_documented",
  dataset_adjust = FALSE, adjust_bh = TRUE
)
fstage_effects[, `:=`(
  n_documented_source_records_total = n_fstage_records,
  n_documented_donors_total = n_fstage_donors
)]

pairing_paths <- file.path(
  BASE,
  c(
    "data/GSE244832/metadata/donor_pairing.csv",
    "data/GSE185477/metadata/donor_pairing.csv",
    "data/GSE202379/metadata/donor_pairing.csv",
    "data/GSE136103/metadata/donor_pairing.csv"
  )
)
membership_paths <- file.path(HS, CELL_TYPES, "module_genes.tsv")
producer_path <- file.path(SCRIPT_DIR, "512_hotspot_v2_donor_refit.R")
input_paths <- c(
  PRIMARY_SCORE_FILE, EQUAL_RUN_FILE, LOO_FILE, NAME_FILE, META_FILE,
  FSTAGE_FILE, GENCODE_FILE, DONOR_LIB, pairing_paths, membership_paths,
  producer_path
)
input_manifest <- add_tabular_counts(
  manifest_for(input_paths, "input"), input_paths
)

v1_anchor_paths <- file.path(V1_ROOT, names(V1_ANCHOR_SHA256))
if (!all(file.exists(v1_anchor_paths))) {
  fail("One or more immutable v1 anchor files are missing")
}
observed_v1_anchor <- vapply(v1_anchor_paths, sha256_file, character(1))
assert_true(
  identical(unname(observed_v1_anchor), unname(V1_ANCHOR_SHA256)),
  paste0(
    "Immutable v1 anchor drift detected: ",
    paste(names(V1_ANCHOR_SHA256)[observed_v1_anchor != V1_ANCHOR_SHA256],
          collapse = ", ")
  )
)

v1_paths <- list.files(V1_ROOT, recursive = TRUE, full.names = TRUE, all.files = TRUE)
v1_paths <- v1_paths[file.exists(v1_paths) & !dir.exists(v1_paths)]
if (!length(v1_paths)) fail("No immutable v1 result files found under ", V1_ROOT)
v1_baseline <- manifest_for(v1_paths, "immutable_v1")

analysis_spec <- data.table(
  release_id = RELEASE_ID,
  primary_score = "pooled-cell donor mean",
  sensitivity_score = "unweighted mean of run scores within donor",
  included_stages = paste(PRIMARY_STAGES, collapse = ";"),
  primary_model = "score ~ stage_ordinal + factor(dataset)",
  stage_encoding = "Healthy=0;Steatosis=1;Steatohepatitis=2",
  dataset_effect = "fixed",
  primary_multiplicity = "BH over complete 117-program family",
  primary_uncertainty = "homoskedastic OLS SE with two-sided residual-df t test",
  robustness_uncertainty = "HC3 sandwich SE with two-sided residual-df t test; BH over 117",
  stability_statistic = "median best-match top-50 Jaccard across finite LOO datasets",
  stability_threshold = 0.5,
  robust_display_rule = paste0(
    "primary OLS q<0.05 & LOO median>=0.5 & equal-run sign agreement & ",
    "primary HC3 q<0.05"
  ),
  external_test_rule = "robust_display & cell_type==hepatocytes",
  fstage_model = "score ~ F_stage_documented; Andrews documented records only; donor collapsed",
  prohibited_fstage_inputs = "F_stage_inferred;F_stage_augmented;F_stage_augmented_clean",
  external_outcomes_read = FALSE
)

summary_dt <- data.table(
  release_id = RELEASE_ID,
  n_programs = nrow(programs),
  n_primary_score_rows = nrow(primary_scores),
  n_equal_run_score_rows = nrow(equal_scores),
  n_primary_estimable = sum(primary_effects$estimable),
  n_primary_q_lt_0_05 = sum(primary_effects$qvalue < 0.05, na.rm = TRUE),
  n_primary_hc3_q_lt_0_05 = sum(primary_effects$hc3_qvalue < 0.05, na.rm = TRUE),
  n_primary_selected = sum(preflight_selection$primary_selected),
  n_robust_display = sum(preflight_selection$robust_display),
  n_external_hepatocyte_programs = sum(
    preflight_selection$robust_display &
      preflight_selection$cell_type == "hepatocytes"
  ),
  n_documented_fstage_records = n_fstage_records,
  n_documented_fstage_donors = n_fstage_donors,
  n_v1_files_protected = nrow(v1_baseline),
  n_input_files = nrow(input_manifest),
  n_tabular_inputs_counted = sum(input_manifest$count_status == "measured"),
  n_gene_counted_inputs = sum(is.finite(input_manifest$gene_count)),
  n_live_plan_inputs = sum(grepl("^docs/plans/", input_manifest$relative_path))
)
assert_true(summary_dt$n_live_plan_inputs == 0L,
            "Mutable plan documents must not enter the frozen input manifest")

if (CHECK_ONLY) {
  cat("[512 check-only] input/model preflight PASS\n")
  print(summary_dt)
  quit(save = "no", status = 0L)
}

if (file.exists(file.path(OUT, "gate_status.tsv"))) {
  fail("Candidate is already sealed; refusing to refit in place: ", OUT)
}
dir.create(WORK, recursive = TRUE, showWarnings = FALSE)

targets <- list(
  input_manifest = input_manifest,
  v1_baseline_manifest = v1_baseline,
  analysis_specification = analysis_spec,
  preflight_summary = summary_dt,
  donor_metadata_audit = meta_obj$audit,
  donor_program_scores_primary = primary_scores,
  donor_program_scores_equal_run = equal_scores,
  primary_effects = primary_effects,
  equal_run_effects = equal_effects,
  stage_dataset_design_audit = stage_design,
  cohort_and_lodo_effects = context_effects,
  source_stability = source_stability,
  documented_fstage_donor_map = fstage_donor,
  documented_fstage_sensitivity = fstage_effects
)

for (nm in names(targets)) {
  write_new_tsv(targets[[nm]], file.path(WORK, paste0(nm, ".tsv")))
}

cat("[512] candidate refit work products written under: ", WORK, "\n", sep = "")
print(summary_dt)
