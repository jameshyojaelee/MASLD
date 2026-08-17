#!/usr/bin/env Rscript
# Freeze the figure-independent Hotspot v2 registry produced by Script 512.
# Refuses to overwrite a sealed candidate and verifies immutable Figure 4 v1.

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
})

options(stringsAsFactors = FALSE)

BASE <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
RELEASE_ID <- Sys.getenv(
  "HOTSPOT_V2_RELEASE_ID",
  unset = "program-context-v2-candidate-2026-08-07"
)
if (!grepl("^[A-Za-z0-9][A-Za-z0-9._-]+$", RELEASE_ID)) {
  stop("HOTSPOT_V2_RELEASE_ID contains unsafe path characters", call. = FALSE)
}

OUT <- file.path(
  BASE, "Analysis/Multimodal_Program_Projection/candidates",
  RELEASE_ID, "hotspot"
)
WORK <- file.path(OUT, "work")
HS <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/hotspot_modules")
SCRIPT_DIR <- file.path(BASE, "Analysis/SingleCell/scripts/hotspot_modules")
V1_ROOT <- file.path(BASE, "Analysis/Multimodal_Program_Projection/results")
GENCODE_FILE <- file.path(BASE, "data/gencode_v49_gene_metadata.tsv.gz")
CHECK_ONLY <- "--check-only" %in% commandArgs(trailingOnly = TRUE) ||
  tolower(Sys.getenv("HOTSPOT_V2_CHECK_ONLY", unset = "false")) %in%
    c("true", "t", "1", "yes", "y")

CELL_TYPES <- c(
  "hepatocytes", "fibroblasts", "macrophages", "cholangiocytes", "tcells"
)
EXPECTED_CENSUS <- c(
  hepatocytes = 30L, fibroblasts = 29L, macrophages = 19L,
  cholangiocytes = 28L, tcells = 11L
)
EXPECTED_PROGRAMS <- sum(EXPECTED_CENSUS)

fail <- function(...) stop(..., call. = FALSE)
assert_true <- function(ok, message) if (!isTRUE(ok)) fail(message)
sha256_file <- function(path) digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
sha256_text <- function(x) digest(x, algo = "sha256", serialize = FALSE)

rel_path <- function(path) {
  root <- paste0(normalizePath(BASE, winslash = "/", mustWork = TRUE), "/")
  resolved <- normalizePath(path, winslash = "/", mustWork = TRUE)
  if (startsWith(resolved, root)) substring(resolved, nchar(root) + 1L) else resolved
}

write_new_tsv <- function(x, path) {
  if (file.exists(path)) fail("Refusing to overwrite frozen candidate artifact: ", path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp.", Sys.getpid())
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  fwrite(x, tmp, sep = "\t", quote = FALSE, na = "NA")
  if (!file.rename(tmp, path)) fail("Atomic rename failed for ", path)
}

manifest_for <- function(paths, role) {
  paths <- unique(paths)
  missing <- paths[!file.exists(paths)]
  if (length(missing)) fail("Missing manifest file(s): ", paste(missing, collapse = ", "))
  info <- file.info(paths)
  data.table(
    role = role,
    relative_path = vapply(paths, rel_path, character(1)),
    bytes = as.numeric(info$size),
    sha256 = vapply(paths, sha256_file, character(1)),
    modified_utc = format(info$mtime, tz = "UTC", usetz = TRUE)
  )
}

build_environment_record <- function() {
  packages <- c("data.table", "digest")
  pkg_rows <- data.table(
    record_type = "r_package",
    name = packages,
    version = vapply(packages, function(pkg) {
      as.character(packageVersion(pkg))
    }, character(1)),
    value = vapply(packages, function(pkg) {
      normalizePath(find.package(pkg), winslash = "/", mustWork = TRUE)
    }, character(1))
  )
  sys <- Sys.info()
  software <- extSoftVersion()
  rbindlist(list(
    data.table(
      record_type = "runtime", name = "R", version = as.character(getRversion()),
      value = paste(R.version$platform, R.version$arch, sep = ";")
    ),
    data.table(
      record_type = "runtime", name = "Rscript", version = as.character(getRversion()),
      value = normalizePath(
        file.path(R.home("bin"), "Rscript"), winslash = "/", mustWork = TRUE
      )
    ),
    pkg_rows,
    data.table(
      record_type = "operating_system", name = as.character(sys[["sysname"]]),
      version = as.character(sys[["release"]]), value = as.character(sys[["machine"]])
    ),
    data.table(
      record_type = "environment", name = "CONDA_PREFIX", version = NA_character_,
      value = fifelse(nzchar(Sys.getenv("CONDA_PREFIX")),
                      Sys.getenv("CONDA_PREFIX"), NA_character_)
    ),
    data.table(
      record_type = "locale", name = "LC_ALL", version = NA_character_,
      value = Sys.getlocale()
    ),
    data.table(
      record_type = "external_library", name = names(software),
      version = as.character(software), value = NA_character_
    )
  ), use.names = TRUE)
}

load_membership <- function() {
  gencode <- fread(
    GENCODE_FILE,
    select = c("gene_id", "gene_name", "ensembl_base"),
    showProgress = FALSE
  )
  gencode[, `:=`(
    gene_name = as.character(gene_name),
    ensembl_base = sub("\\.[0-9]+$", "", as.character(ensembl_base))
  )]
  gencode <- unique(gencode[
    !is.na(gene_name) & nzchar(gene_name) &
      !is.na(ensembl_base) & nzchar(ensembl_base),
    .(ensembl_base, gene_name)
  ])
  base_cardinality <- gencode[, .(n_symbols_for_base = uniqueN(gene_name)),
                              by = ensembl_base]
  symbol_cardinality <- gencode[, .(n_bases_for_symbol = uniqueN(ensembl_base)),
                                by = gene_name]
  unambiguous_map <- merge(
    merge(gencode, base_cardinality, by = "ensembl_base"),
    symbol_cardinality, by = "gene_name"
  )[n_symbols_for_base == 1L & n_bases_for_symbol == 1L,
    .(ensembl_base, mapped_symbol = gene_name)]
  if (anyDuplicated(unambiguous_map$ensembl_base) ||
      anyDuplicated(unambiguous_map$mapped_symbol)) {
    fail("GENCODE v49 unambiguous symbol map is not one-to-one")
  }

  x <- rbindlist(lapply(CELL_TYPES, function(ct) {
    path <- file.path(HS, ct, "module_genes.tsv")
    if (!file.exists(path)) fail("Missing membership: ", path)
    z <- fread(path)
    if (!all(c("gene", "module", "weight") %in% names(z))) fail("Malformed: ", path)
    z[, `:=`(
      cell_type = ct,
      module = as.integer(module),
      source_gene = as.character(gene),
      canonical_gene = fifelse(
        grepl("^ENSG[0-9]+\\.[0-9]+$", as.character(gene)),
        sub("\\.[0-9]+$", "", as.character(gene)),
        as.character(gene)
      ),
      source_weight = as.numeric(weight)
    )]
    z[, .(cell_type, module, source_gene, canonical_gene, source_weight)]
  }))
  if (anyNA(x$module) || any(!nzchar(x$canonical_gene)) ||
      any(!is.finite(x$source_weight)) || any(x$source_weight <= 0)) {
    fail("Invalid module membership or positive weight")
  }
  if (anyDuplicated(x[, .(cell_type, module, canonical_gene)])) {
    fail("Ensembl-version stripping creates duplicate genes within a program")
  }

  x[, mapped_symbol := NA_character_]
  x[unambiguous_map, on = .(canonical_gene = ensembl_base),
    mapped_symbol := i.mapped_symbol]
  x[unambiguous_map, on = .(canonical_gene = mapped_symbol),
    mapped_symbol := i.mapped_symbol]
  x[, mapped_symbol_status := fcase(
    !is.na(mapped_symbol) & grepl("^ENSG[0-9]+$", canonical_gene),
      "gencode_v49_unambiguous_ensembl_to_symbol",
    !is.na(mapped_symbol), "gencode_v49_unique_symbol_confirmed",
    default = "unmapped_or_ambiguous"
  )]

  programs <- unique(x[, .(cell_type, module)])
  census <- programs[, .N, by = cell_type]
  observed <- setNames(census$N, census$cell_type)
  assert_true(setequal(names(observed), names(EXPECTED_CENSUS)), "Lineage set drift")
  assert_true(all(observed[names(EXPECTED_CENSUS)] == EXPECTED_CENSUS),
              "Lineage census drift")
  assert_true(nrow(programs) == EXPECTED_PROGRAMS, "Expected exactly 117 programs")

  x[, original_l1_weight := source_weight / sum(source_weight),
    by = .(cell_type, module)]
  program_summary <- x[, .(
    n_source_genes = uniqueN(source_gene),
    n_canonical_genes = uniqueN(canonical_gene),
    n_mapped_symbols = sum(!is.na(mapped_symbol)),
    n_positive_source_weights = sum(source_weight > 0),
    source_weight_sum = sum(source_weight),
    source_weight_min = min(source_weight),
    source_weight_median = median(source_weight),
    source_weight_max = max(source_weight),
    source_weights_all_positive = all(source_weight > 0),
    weight_transform = "original_l1_weight=source_weight/sum(source_weight)"
  ), by = .(cell_type, module)]
  x[, canonical_weight_text := sprintf("%.17e", source_weight)]
  hash_dt <- x[order(cell_type, module, canonical_gene, canonical_weight_text), .(
    canonical_serialization = paste(
      paste(cell_type, canonical_gene, canonical_weight_text, sep = "\t"),
      collapse = "\n"
    )
  ), by = .(cell_type, module)]
  hash_dt[, membership_sha256 := vapply(canonical_serialization, sha256_text, character(1))]
  hash_dt[, program_uid := paste0(
    "hotspot_", cell_type, "_", substr(membership_sha256, 1L, 16L)
  )]
  if (anyDuplicated(hash_dt$program_uid)) fail("program_uid collision")

  x <- merge(
    x, hash_dt[, .(cell_type, module, membership_sha256, program_uid)],
    by = c("cell_type", "module"), all.x = TRUE
  )
  setorder(x, cell_type, module, canonical_gene, source_gene)
  programs <- merge(
    hash_dt[, .(cell_type, module, program_uid, membership_sha256)],
    program_summary, by = c("cell_type", "module"), all.x = TRUE
  )
  list(membership = x, programs = programs)
}

membership_obj <- load_membership()
membership <- membership_obj$membership
program_ids <- membership_obj$programs
environment_record <- build_environment_record()
environment_record[, release_id := RELEASE_ID]
setcolorder(environment_record, "release_id")
assert_true(nrow(environment_record[record_type == "r_package"]) == 2L,
            "Environment record lacks required R packages")
assert_true(all(program_ids$n_source_genes > 0L) &&
              all(program_ids$n_positive_source_weights ==
                    program_ids$n_source_genes) &&
              all(program_ids$source_weights_all_positive),
            "Membership gene/weight preflight failed")

loo <- fread(file.path(HS, "loo_stability.tsv"))
module_names <- fread(file.path(HS, "module_names.tsv"))
loo[, module := as.integer(module)]
module_names[, module := as.integer(module)]
assert_true(nrow(loo) == EXPECTED_PROGRAMS &&
              !anyDuplicated(loo[, .(cell_type, module)]),
            "LOO stability must contain 117 unique programs")
assert_true(nrow(module_names) == EXPECTED_PROGRAMS &&
              !anyDuplicated(module_names[, .(cell_type, module)]),
            "Module names must contain 117 unique programs")

loo_columns <- grep("^loo_", names(loo), value = TRUE)
assert_true(length(loo_columns) >= 5L &&
              all(vapply(loo[, ..loo_columns], is.numeric, logical(1))),
            "LOO table lacks the numeric per-dataset Jaccard family")
setnames(loo, "stability_score", "source_stability_mean")
if ("stability_fail" %in% names(loo)) {
  setnames(loo, "stability_fail", "source_stability_fail_mean")
}
loo[, stability_mean := rowMeans(.SD, na.rm = TRUE), .SDcols = loo_columns]
loo[!is.finite(stability_mean), stability_mean := NA_real_]
loo[, stability_median := apply(.SD, 1L, function(z) {
  z <- z[is.finite(z)]
  if (length(z)) median(z) else NA_real_
}), .SDcols = loo_columns]
assert_true(
  max(abs(loo$source_stability_mean - loo$stability_mean), na.rm = TRUE) <= 5.1e-5,
  "Published stability_score does not reproduce the producer's rounded LOO mean"
)
loo[, `:=`(
  stability_score = stability_median,
  stability_statistic = "median_best_match_top50_jaccard"
)]

if (CHECK_ONLY && !dir.exists(WORK)) {
  cat("[513 check-only] membership/name/LOO/environment preflight PASS\n")
  print(data.table(
    n_programs = nrow(program_ids),
    n_membership_rows = nrow(membership),
    n_mapped_symbols = sum(!is.na(membership$mapped_symbol)),
    n_environment_records = nrow(environment_record),
    n_stable_by_mean = sum(loo$stability_mean >= 0.5, na.rm = TRUE),
    n_stable_by_median = sum(loo$stability_median >= 0.5, na.rm = TRUE),
    n_mean_median_threshold_flips = sum(
      (loo$stability_mean >= 0.5) != (loo$stability_median >= 0.5),
      na.rm = TRUE
    )
  ))
  quit(save = "no", status = 0L)
}

required_work <- c(
  "input_manifest.tsv", "v1_baseline_manifest.tsv", "analysis_specification.tsv",
  "preflight_summary.tsv", "donor_program_scores_primary.tsv",
  "donor_program_scores_equal_run.tsv", "primary_effects.tsv",
  "equal_run_effects.tsv", "stage_dataset_design_audit.tsv",
  "cohort_and_lodo_effects.tsv", "source_stability.tsv",
  "documented_fstage_donor_map.tsv", "documented_fstage_sensitivity.tsv"
)
work_paths <- file.path(WORK, required_work)
missing_work <- work_paths[!file.exists(work_paths)]
if (length(missing_work)) {
  fail("Run Script 512 first; missing: ", paste(basename(missing_work), collapse = ", "))
}

primary <- fread(file.path(WORK, "primary_effects.tsv"))
equal <- fread(file.path(WORK, "equal_run_effects.tsv"))
source_stability <- fread(file.path(WORK, "source_stability.tsv"))
assert_true(nrow(primary) == EXPECTED_PROGRAMS &&
              !anyDuplicated(primary[, .(cell_type, module)]),
            "Primary refit must contain 117 unique programs")
assert_true(nrow(equal) == EXPECTED_PROGRAMS &&
              !anyDuplicated(equal[, .(cell_type, module)]),
            "Equal-run refit must contain 117 unique programs")

q_rederived <- rep(NA_real_, nrow(primary))
q_idx <- which(is.finite(primary$pvalue))
q_rederived[q_idx] <- p.adjust(primary$pvalue[q_idx], method = "BH", n = EXPECTED_PROGRAMS)
assert_true(isTRUE(all.equal(primary$qvalue, q_rederived, tolerance = 1e-13,
                            check.attributes = FALSE)),
            "Primary q-values do not rederive over the 117-program family")
hc3_q_rederived <- rep(NA_real_, nrow(primary))
hc3_q_idx <- which(is.finite(primary$hc3_pvalue))
hc3_q_rederived[hc3_q_idx] <- p.adjust(
  primary$hc3_pvalue[hc3_q_idx], method = "BH", n = EXPECTED_PROGRAMS
)
assert_true(isTRUE(all.equal(primary$hc3_qvalue, hc3_q_rederived,
                            tolerance = 1e-13, check.attributes = FALSE)),
            "Primary HC3 q-values do not rederive over the 117-program family")

setnames(
  primary,
  setdiff(names(primary), c("cell_type", "module")),
  paste0("primary_", setdiff(names(primary), c("cell_type", "module")))
)
setnames(
  equal,
  setdiff(names(equal), c("cell_type", "module")),
  paste0("equal_run_", setdiff(names(equal), c("cell_type", "module")))
)

registry <- Reduce(
  function(x, y) merge(x, y, by = c("cell_type", "module"), all.x = TRUE),
  list(program_ids, module_names, loo, primary, equal, source_stability)
)
setorder(registry, cell_type, module)
assert_true(nrow(registry) == EXPECTED_PROGRAMS, "Registry join changed 117-row universe")

registry[, `:=`(
  release_id = RELEASE_ID,
  model = "score ~ stage_ordinal + factor(dataset)",
  included_stages = "Healthy;Steatosis;Steatohepatitis",
  score_definition = "pooled-cell donor mean",
  sensitivity_score_definition = "unweighted mean of run scores within donor",
  primary_selected = primary_estimable & is.finite(primary_qvalue) &
    primary_qvalue < 0.05 & is.finite(stability_median) & stability_median >= 0.5,
  scoring_direction_agree = primary_estimable & equal_run_estimable &
    is.finite(primary_beta) & is.finite(equal_run_beta) & primary_beta != 0 &
    equal_run_beta != 0 & sign(primary_beta) == sign(equal_run_beta),
  selected_unstable = primary_estimable & is.finite(primary_qvalue) &
    primary_qvalue < 0.05 & (!is.finite(stability_median) | stability_median < 0.5),
  hc3_supported = primary_estimable & is.finite(primary_hc3_qvalue) &
    primary_hc3_qvalue < 0.05
)]
registry[, robust_display := primary_selected & scoring_direction_agree & hc3_supported]
registry[, selected_hc3_fragile :=
           primary_selected & scoring_direction_agree & !hc3_supported]
registry[, external_test_eligible := robust_display & cell_type == "hepatocytes"]
registry[, tested_negative := primary_estimable & is.finite(primary_qvalue) &
           primary_qvalue >= 0.05]
registry[, registry_state := fcase(
  !primary_estimable, "untestable",
  selected_unstable, "selected_unstable",
  primary_selected & !scoring_direction_agree, "selected_direction_discordant",
  selected_hc3_fragile, "selected_hc3_fragile",
  robust_display, "robust_display",
  tested_negative, "tested_nonsignificant",
  default = "testable_unclassified"
)]
assert_true(all(registry$n_source_genes > 0L) &&
              all(registry$n_source_genes == registry$n_canonical_genes) &&
              all(registry$n_positive_source_weights == registry$n_source_genes) &&
              all(registry$source_weights_all_positive),
            "Registry gene/positive-weight summary failed")

producer_hash <- sha256_file(file.path(SCRIPT_DIR, "512_hotspot_v2_donor_refit.R"))
freezer_hash <- sha256_file(file.path(SCRIPT_DIR, "513_freeze_hotspot_v2_registry.R"))
registry[, `:=`(
  producer_sha256 = producer_hash,
  freezer_sha256 = freezer_hash,
  gencode_v49_sha256 = sha256_file(GENCODE_FILE),
  external_outcomes_read = FALSE
)]

join_uid <- function(x) {
  merge(program_ids, x, by = c("cell_type", "module"), all.y = TRUE)
}

donor_primary <- join_uid(fread(file.path(WORK, "donor_program_scores_primary.tsv")))
donor_equal <- join_uid(fread(file.path(WORK, "donor_program_scores_equal_run.tsv")))
design_audit <- join_uid(fread(file.path(WORK, "stage_dataset_design_audit.tsv")))
context_effects <- join_uid(fread(file.path(WORK, "cohort_and_lodo_effects.tsv")))
fstage_sensitivity <- join_uid(fread(file.path(WORK, "documented_fstage_sensitivity.tsv")))
fstage_map <- fread(file.path(WORK, "documented_fstage_donor_map.tsv"))
input_manifest <- fread(file.path(WORK, "input_manifest.tsv"))
analysis_spec <- fread(file.path(WORK, "analysis_specification.tsv"))
required_manifest_fields <- c("row_count", "gene_count", "count_status")
assert_true(all(required_manifest_fields %in% names(input_manifest)),
            "Input manifest lacks row/gene count audit fields")
assert_true(!any(input_manifest$count_status == "not_measured"),
            "Input manifest contains unmeasured count status")
gencode_rel <- rel_path(GENCODE_FILE)
gencode_manifest <- input_manifest[relative_path == gencode_rel]
assert_true(nrow(gencode_manifest) == 1L &&
              is.finite(gencode_manifest$row_count) &&
              is.finite(gencode_manifest$gene_count),
            "GENCODE mapping source lacks measured row/gene counts")

fig2_source <- registry[, .(
  release_id, program_uid, membership_sha256, cell_type, module, module_name,
  module_hallmark, module_top_pathway, beta = primary_beta, se = primary_se,
  statistic = primary_statistic, pvalue = primary_pvalue, qvalue = primary_qvalue,
  hc3_se = primary_hc3_se, hc3_statistic = primary_hc3_statistic,
  hc3_pvalue = primary_hc3_pvalue, hc3_qvalue = primary_hc3_qvalue,
  max_leverage = primary_max_leverage, direction = primary_direction,
  stability_score, stability_mean, stability_median, primary_selected,
  hc3_supported, scoring_direction_agree, selected_hc3_fragile,
  robust_display, source_stability,
  source_stability_reason, n_donors = primary_n_donors,
  n_datasets = primary_n_datasets, n_healthy = primary_n_stage0,
  n_steatosis = primary_n_stage1, n_steatohepatitis = primary_n_stage2
)]

tested_universe <- registry[, .(
  release_id, program_uid, membership_sha256, cell_type, module, module_name,
  estimable = primary_estimable, failure_reason = primary_failure_reason,
  pvalue = primary_pvalue, qvalue = primary_qvalue,
  hc3_pvalue = primary_hc3_pvalue, hc3_qvalue = primary_hc3_qvalue,
  primary_selected, selected_unstable, hc3_supported, selected_hc3_fragile,
  scoring_direction_agree, robust_display, tested_negative,
  external_test_eligible, registry_state
)]

if (CHECK_ONLY) {
  cat("[513 check-only] work-product and registry validation PASS\n")
  print(registry[, .(
    n_programs = .N,
    n_estimable = sum(primary_estimable),
    n_primary_selected = sum(primary_selected),
    n_hc3_supported = sum(hc3_supported),
    n_selected_hc3_fragile = sum(selected_hc3_fragile),
    n_robust_display = sum(robust_display),
    n_external = sum(external_test_eligible)
  )])
  quit(save = "no", status = 0L)
}

if (file.exists(file.path(OUT, "gate_status.tsv"))) {
  fail("Candidate registry is already sealed; refusing in-place rewrite")
}

baseline <- fread(file.path(WORK, "v1_baseline_manifest.tsv"))
current_v1_paths <- file.path(BASE, baseline$relative_path)
assert_true(all(file.exists(current_v1_paths)), "An immutable v1 file disappeared")
v1_now <- manifest_for(current_v1_paths, "immutable_v1_final")
v1_compare <- merge(
  baseline[, .(
    relative_path, baseline_bytes = bytes, baseline_sha256 = sha256
  )],
  v1_now[, .(
    relative_path, final_bytes = bytes, final_sha256 = sha256
  )],
  by = "relative_path", all = TRUE
)
v1_compare[, unchanged := baseline_bytes == final_bytes & baseline_sha256 == final_sha256]
assert_true(nrow(v1_compare) == nrow(baseline) && all(v1_compare$unchanged),
            "Immutable Figure 4 v1 changed during v2 production")

final_targets <- c(
  "program_membership_v2.tsv", "program_registry_v2.tsv",
  "donor_program_scores_primary.tsv", "donor_program_scores_equal_run.tsv",
  "stage_dataset_design_audit.tsv", "cohort_and_lodo_effects.tsv",
  "documented_fstage_donor_map.tsv", "documented_fstage_sensitivity.tsv",
  "fig2_program_source.tsv", "tested_universe.tsv", "input_manifest.tsv",
  "analysis_specification.tsv", "environment_record.tsv", "v1_preservation.tsv",
  "external_test_programs.tsv", "gate_status.tsv", "release_manifest.tsv"
)
existing <- file.path(OUT, final_targets)
if (any(file.exists(existing))) {
  fail("Refusing to overwrite existing final artifact(s): ",
       paste(basename(existing[file.exists(existing)]), collapse = ", "))
}

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
write_new_tsv(membership, file.path(OUT, "program_membership_v2.tsv"))
write_new_tsv(registry, file.path(OUT, "program_registry_v2.tsv"))
registry_sha <- sha256_file(file.path(OUT, "program_registry_v2.tsv"))
membership_sha <- sha256_file(file.path(OUT, "program_membership_v2.tsv"))

external <- registry[external_test_eligible == TRUE, .(
  release_id, registry_sha256 = registry_sha, program_uid, membership_sha256,
  cell_type, module, module_name, primary_beta, primary_se, primary_pvalue,
  primary_qvalue, primary_hc3_se, primary_hc3_pvalue, primary_hc3_qvalue,
  stability_score, stability_mean, stability_median, source_stability,
  expected_direction = primary_direction
)]

write_new_tsv(donor_primary, file.path(OUT, "donor_program_scores_primary.tsv"))
write_new_tsv(donor_equal, file.path(OUT, "donor_program_scores_equal_run.tsv"))
write_new_tsv(design_audit, file.path(OUT, "stage_dataset_design_audit.tsv"))
write_new_tsv(context_effects, file.path(OUT, "cohort_and_lodo_effects.tsv"))
write_new_tsv(fstage_map, file.path(OUT, "documented_fstage_donor_map.tsv"))
write_new_tsv(fstage_sensitivity, file.path(OUT, "documented_fstage_sensitivity.tsv"))
write_new_tsv(fig2_source, file.path(OUT, "fig2_program_source.tsv"))
write_new_tsv(tested_universe, file.path(OUT, "tested_universe.tsv"))
write_new_tsv(input_manifest, file.path(OUT, "input_manifest.tsv"))
write_new_tsv(analysis_spec, file.path(OUT, "analysis_specification.tsv"))
write_new_tsv(environment_record, file.path(OUT, "environment_record.tsv"))
write_new_tsv(v1_compare, file.path(OUT, "v1_preservation.tsv"))
write_new_tsv(external, file.path(OUT, "external_test_programs.tsv"))

environment_sha <- sha256_file(file.path(OUT, "environment_record.tsv"))

producer_paths <- file.path(
  SCRIPT_DIR,
  c(
    "511_pooled_cell_donor_scores.py",
    "512_hotspot_v2_donor_refit.R",
    "513_freeze_hotspot_v2_registry.R",
    "514_validate_hotspot_v2_registry.R",
    "515_correct_stage_metadata.R",
    Sys.getenv("HOTSPOT_V2_SBATCH_NAME", unset = "run_hotspot_v2_candidate.sbatch")
  )
)
environment_path <- file.path(OUT, "environment_record.tsv")
output_paths <- file.path(
  OUT,
  setdiff(
    final_targets,
    c("release_manifest.tsv", "environment_record.tsv", "gate_status.tsv")
  )
)
release_manifest <- rbindlist(list(
  input_manifest,
  manifest_for(work_paths, "intermediate"),
  manifest_for(producer_paths, "producer"),
  manifest_for(environment_path, "environment"),
  manifest_for(output_paths, "output")
), use.names = TRUE, fill = TRUE)
release_manifest[, release_id := RELEASE_ID]
write_new_tsv(release_manifest, file.path(OUT, "release_manifest.tsv"))
release_manifest_sha <- sha256_file(file.path(OUT, "release_manifest.tsv"))

# Final internal consistency before the freezer gate is written last.
assert_true(nrow(registry) == 117L && nrow(tested_universe) == 117L,
            "Frozen row count changed")
assert_true(fsetequal(
  external[, .(program_uid)],
  registry[external_test_eligible == TRUE, .(program_uid)]
), "External-test subset is not deterministic")
assert_true(all(v1_compare$unchanged), "v1 preservation check changed unexpectedly")

sealed_at <- format(Sys.time(), tz = "UTC", usetz = TRUE)
gate <- data.table(
  release_id = RELEASE_ID,
  status = "frozen_pending_independent_validation",
  registry_sha256 = registry_sha,
  membership_table_sha256 = membership_sha,
  environment_record_sha256 = environment_sha,
  release_manifest_sha256 = release_manifest_sha,
  n_programs = nrow(registry),
  n_primary_estimable = sum(registry$primary_estimable),
  n_primary_selected = sum(registry$primary_selected),
  n_hc3_supported = sum(registry$hc3_supported),
  n_selected_hc3_fragile = sum(registry$selected_hc3_fragile),
  n_robust_display = sum(registry$robust_display),
  n_external_test_programs = nrow(external),
  n_documented_fstage_records = 58L,
  n_documented_fstage_donors = 46L,
  external_outcomes_read = FALSE,
  independent_validation_passed = FALSE,
  v1_unchanged = TRUE,
  sealed_at_utc = sealed_at
)
write_new_tsv(gate, file.path(OUT, "gate_status.tsv"))

cat("[513] frozen Hotspot v2 candidate: ", OUT, "\n", sep = "")
print(gate)
