# Shared helpers for CONTINUOUS-AXIS-BENCHMARK-v1.
#
# Style and guarantees follow scripts/figures/build_pi_stage_extensions.R:
# refuse-to-overwrite writes, explicit failure, checksummed inputs and outputs,
# and every seed recorded. Nothing here computes a result; it only makes results
# reproducible and hard to overwrite by accident.

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

options(digits = 17, scipen = 999)

EXEMPTION_ID <- "ORDERING-EXEMPTION-2026-08-14"
WORKSTREAM_ID <- "CONTINUOUS-AXIS-BENCHMARK-v1"

project_root <- function() {
  root <- Sys.getenv("MASLD_PROJECT_ROOT", "")
  if (!nzchar(root)) stop("MASLD_PROJECT_ROOT is unset", call. = FALSE)
  root
}

fail <- function(...) stop(paste0(...), call. = FALSE)

assert_true <- function(value, message) {
  if (!isTRUE(value)) fail(message)
  invisible(TRUE)
}

sha256_file <- function(path) {
  assert_true(file.exists(path), paste0("Cannot hash a file that does not exist: ", path))
  sub("[[:space:]].*$", "", system2("sha256sum", shQuote(path), stdout = TRUE)[[1]])
}

log_step <- function(...) {
  cat(format(Sys.time(), "%H:%M:%S"), "|", paste0(...), "\n")
  flush.console()
}

# Writes refuse to clobber. A rerun into an existing directory is a mistake, not
# a convenience, because it silently mixes two executions in one output root.
write_tsv_once <- function(x, path) {
  assert_true(!file.exists(path), paste0("Refusing to overwrite an existing output: ", path))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  fwrite(x, path, sep = "\t")
  invisible(path)
}

write_json_once <- function(x, path) {
  assert_true(!file.exists(path), paste0("Refusing to overwrite an existing output: ", path))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write_json(x, path, auto_unbox = TRUE, pretty = TRUE, digits = 17, null = "null")
  invisible(path)
}

# ---------------------------------------------------------------- provenance --

# Records which exemption licensed this run and the exact content of the config
# at run time. The config is gitignored (.gitignore line 8 is a blanket '*'), so
# this hash is the only durable evidence of what was licensed.
stamp_exemption <- function(out_root) {
  cfg_path <- file.path(project_root(), "config/method_exemptions.json")
  assert_true(file.exists(cfg_path),
              "config/method_exemptions.json is absent, so no exemption is live and this workstream must not run")
  cfg <- fromJSON(cfg_path, simplifyVector = FALSE)
  live <- Filter(function(e) identical(e$exemption_id, EXEMPTION_ID) &&
                   identical(e$status, "active"), cfg$exemptions)
  assert_true(length(live) == 1L,
              paste0("Exemption ", EXEMPTION_ID, " is not active in config/method_exemptions.json"))
  stamp <- list(
    exemption_id = EXEMPTION_ID,
    workstream = WORKSTREAM_ID,
    config_path = "config/method_exemptions.json",
    config_sha256 = sha256_file(cfg_path),
    licensed_methods = live[[1]]$methods,
    licensed_terms = live[[1]]$terms,
    scope_paths = live[[1]]$scope_paths,
    claim_firewall = live[[1]]$claim_firewall,
    stamped_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  )
  write_json_once(stamp, file.path(out_root, "EXEMPTION.json"))
  invisible(stamp)
}

write_input_manifest <- function(paths, out_root) {
  rows <- rbindlist(lapply(names(paths), function(role) {
    p <- paths[[role]]
    data.table(input_role = role, path = p,
               exists = file.exists(p),
               size_bytes = if (file.exists(p)) file.size(p) else NA_real_,
               sha256 = if (file.exists(p)) sha256_file(p) else NA_character_)
  }))
  write_tsv_once(rows, file.path(out_root, "input_checksums.tsv"))
  rows
}

write_output_manifest <- function(out_root) {
  files <- list.files(out_root, recursive = TRUE, full.names = TRUE)
  files <- files[basename(files) != "output_checksums.tsv"]
  rows <- rbindlist(lapply(files, function(p) {
    data.table(relative_path = sub(paste0("^", out_root, "/"), "", p),
               size_bytes = file.size(p), sha256 = sha256_file(p))
  }))
  setorder(rows, relative_path)
  write_tsv_once(rows, file.path(out_root, "output_checksums.tsv"))
  rows
}

write_run_parameters <- function(out_root, extra = list()) {
  pre <- read_prespec()
  params <- c(list(
    workstream = WORKSTREAM_ID,
    exemption_id = EXEMPTION_ID,
    seed_master = pre$seeds$master,
    seed_perm_base = pre$seeds$perm_base,
    seed_sim_base = pre$seeds$sim_base,
    seed_split_base = pre$seeds$split_base,
    R_version = R.version.string,
    run_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  ), extra)
  writeLines(paste0(names(params), "\t", unlist(lapply(params, as.character))),
             file.path(out_root, "run_parameters.txt"))
  writeLines(capture.output(sessionInfo()), file.path(out_root, "sessionInfo.txt"))
  invisible(params)
}

# ------------------------------------------------------------- prespec + data --

read_prespec <- function() {
  path <- file.path(project_root(),
                    "scripts/analysis/continuous_axis_benchmark/00_prespecification.json")
  assert_true(file.exists(path), "00_prespecification.json is absent; nothing may run before it exists")
  fromJSON(path, simplifyVector = FALSE)
}

SUBSTRATE <- list(
  dge = paste0("RNA-seq/results/manuscript_release/candidates/",
               "resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BG001-DECISION/",
               "arms/F_five/results/integration/merged_dge.rds"),
  manifest = paste0("figures/candidates/pi-figure-redesign-2026-08-13-v3/analysis/",
                    "stage_extensions/five_cohort_sample_manifest.tsv"),
  stage_results = paste0("figures/candidates/pi-figure-redesign-2026-08-13-v3/analysis/",
                         "stage_extensions/stage_extension_all_gene_results.tsv"),
  nmf_k6 = "RNA-seq/results/subtypes/nmf_assignments.csv",
  nmf_k4 = "RNA-seq/results/subtypes/nmf_k4_programs.csv",
  composition = "RNA-seq/results/celltype_attribution/persample_celltype_proportions.csv"
)

substrate_path <- function(key) file.path(project_root(), SUBSTRATE[[key]])

# The 844 are already one cross-sectional biopsy per participant; this asserts
# that rather than aggregating, matching build_pi_stage_extensions.R.
load_manifest <- function() {
  meta <- fread(substrate_path("manifest"))
  assert_true(nrow(meta) == 844L,
              sprintf("Expected 844 participants in the manifest; found %d", nrow(meta)))
  assert_true(!anyDuplicated(meta$sample_id), "sample_id is not unique: the biological unit has drifted")
  assert_true(!anyDuplicated(meta$analysis_unit_id), "analysis_unit_id is not unique")
  assert_true(!anyNA(meta$dataset) && !anyNA(meta$inferred_sex),
              "dataset or inferred_sex is missing for at least one participant")
  meta[, dataset := factor(dataset)]
  meta[, inferred_sex := factor(inferred_sex)]
  meta[]
}

POP_C_COHORTS <- c("GSE130970", "GSE135251", "GSE162694")

population <- function(meta, id) {
  switch(id,
    "POP-A" = meta,
    "POP-B" = meta[!is.na(fibrosis_stage)],
    "POP-C" = meta[dataset %in% POP_C_COHORTS & fibrosis_stage %in% 0:4],
    "POP-D" = meta[!is.na(nas_score)],
    "POP-E" = meta[!is.na(fibrosis_stage) & !is.na(nas_score)],
    "POP-F" = meta[dataset == "GSE130970"],
    fail("Unknown population id: ", id))
}

# Cohort mean-shifts dominate a raw axis, so every evaluation ranks within cohort.
rank_within <- function(value, group) {
  out <- rep(NA_real_, length(value))
  for (g in unique(group)) {
    idx <- which(group == g & !is.na(value))
    if (!length(idx)) next
    out[idx] <- (rank(value[idx], ties.method = "average") - 0.5) / length(idx)
  }
  out
}
