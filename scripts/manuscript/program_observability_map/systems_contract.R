# The shared contract, as callable checks.
#
# Every system in the program framework calls these. They exist because each one
# corresponds to a mistake already made in this workstream, and a comment in a
# README does not stop a mistake from recurring.

suppressPackageStartupMessages({
  library(data.table)
})

# 1. Frozen programs. The registry and membership are SHA-pinned; a system that
# alters them is not projecting the frozen programs, it is making new ones.
assert_frozen_programs <- function(registry_path, membership_path,
                                   expected_programs = 117L) {
  registry <- fread(registry_path)
  membership <- fread(membership_path)
  assert_true(nrow(registry) == expected_programs,
              sprintf("Program registry must contain %d programs; found %d",
                      expected_programs, nrow(registry)))
  assert_true(uniqueN(registry$program_uid) == expected_programs,
              "Program UIDs must be unique")
  assert_true(all(registry$program_uid %in% membership$program_uid),
              "A registry program has no membership rows")
  invisible(list(registry = registry, membership = membership))
}

# 2. Declared biological unit. The anti-pattern this exists to stop: asserting
# uniqueN(participant_id) == nrow(meta) when participant_id was silently filled
# from sample_id, which compares a column to itself and always passes.
assert_biological_unit <- function(meta, unit_column, cohorts = NULL,
                                   pairing_root = NULL) {
  assert_true(unit_column %in% names(meta),
              sprintf("Biological unit column '%s' is absent; do not fall back to the sample id",
                      unit_column))
  assert_true(!identical(as.character(meta[[unit_column]]), as.character(meta$sample_id)) ||
                !is.null(pairing_root),
              paste("The unit column is identical to sample_id, so the uniqueness check is",
                    "vacuous. Supply pairing_root so absence of repeat sampling can be",
                    "verified externally."))
  if (!is.null(pairing_root) && !is.null(cohorts)) {
    pairing <- file.path(pairing_root, cohorts, "metadata", "donor_pairing.csv")
    present <- pairing[file.exists(pairing)]
    assert_true(!length(present),
                paste0("A cohort has a donor-pairing table, so repeat sampling is possible ",
                       "and the model needs a donor term: ", paste(present, collapse = ", ")))
  }
  invisible(TRUE)
}

# 3. No count without its calibration row.
assert_counts_calibrated <- function(counts, calibration_rows) {
  assert_true(nrow(calibration_rows) > 0L,
              "No calibration rows: every reported count needs a measured null false-call rate")
  assert_true(all(calibration_rows$reportable),
              paste0("A view failed the calibration gate and its counts must be withheld: ",
                     paste(calibration_rows[reportable == FALSE, view], collapse = ", ")))
  invisible(TRUE)
}

# 4. Negatives carry a minimum detectable effect.
assert_negatives_have_mde <- function(result, support_column, mde_column) {
  negatives <- result[get(support_column) %in%
                        c("unsupported", "informative_null", "indeterminate")]
  if (!nrow(negatives)) return(invisible(TRUE))
  assert_true(all(is.finite(negatives[[mde_column]])),
              "A negative call has no minimum detectable effect; it is not interpretable")
  invisible(TRUE)
}

# 6. Isolation. The v8 script root is hashed into a sealed contract.
assert_v8_untouched <- function(project_root, v8_contract_path) {
  if (!file.exists(v8_contract_path)) return(invisible(TRUE))
  contract <- jsonlite::read_json(v8_contract_path, simplifyVector = TRUE)
  manifest <- contract$code_manifest
  drift <- vapply(names(manifest), function(relative) {
    path <- file.path(project_root, relative)
    !file.exists(path) ||
      sub("[[:space:]].*$", "", system2("sha256sum", path, stdout = TRUE)[[1]]) != manifest[[relative]]
  }, logical(1))
  assert_true(!any(drift),
              paste0("The sealed v8 script root has drifted: ",
                     paste(names(manifest)[drift], collapse = ", ")))
  invisible(TRUE)
}

assert_holdout_sealed <- function(meta, holdout_cohort, crosswalk_path = NULL) {
  assert_true(!holdout_cohort %in% as.character(meta$dataset),
              paste0("HOLDOUT LEAK: ", holdout_cohort, " is present"))
  if (!is.null(crosswalk_path) && file.exists(crosswalk_path)) {
    tokens <- fread(crosswalk_path)[, unique(participant_token)]
    assert_true(!any(meta$sample_id %in% tokens),
                "HOLDOUT LEAK: a sample matches a holdout participant token")
  }
  invisible(TRUE)
}

# 7. Language. Cross-sectional data cannot carry an ordering claim, and the
# terms below have each leaked back into this project after being retired.
FORBIDDEN_TERMS <- c("progression", "trajectory", "continuum", "pseudotime",
                     "sliding window", "sliding-window", "severity axis",
                     "molecular clock", "disease ordering")

# 7a. Scoped exemptions. A workstream whose OBJECT OF STUDY is an ordering method
# cannot describe itself without the retired vocabulary. Rather than delete the
# guard -- which is why these terms leaked back before -- an exemption licenses
# named terms inside named paths only. Callers that pass neither exemption_id nor
# caller_path get the original unconditional behaviour, so every pre-existing call
# site is unchanged.
EXEMPTION_CONFIG_PATH <- "config/method_exemptions.json"

`%||%` <- function(a, b) if (is.null(a)) b else a

.read_exemptions <- function(project_root = Sys.getenv("MASLD_PROJECT_ROOT", ".")) {
  path <- file.path(project_root, EXEMPTION_CONFIG_PATH)
  # A missing config means zero exemptions. Fail closed, never open.
  if (!file.exists(path)) {
    return(list(exemptions = list(), never = character(0), sha256 = NA_character_))
  }
  cfg <- jsonlite::read_json(path, simplifyVector = FALSE)
  sha <- sub("[[:space:]].*$", "", system2("sha256sum", path, stdout = TRUE)[[1]])
  active <- Filter(function(e) identical(e$status, "active"), cfg$exemptions)
  list(exemptions = active,
       never = tolower(unlist(cfg$never_exemptible_terms %||% list())),
       sha256 = sha)
}

# Returns the terms an exemption licenses for this caller, or character(0).
.licensed_terms <- function(exemption_id, caller_path,
                            project_root = Sys.getenv("MASLD_PROJECT_ROOT", ".")) {
  if (is.null(exemption_id)) return(character(0))
  cfg <- .read_exemptions(project_root)
  match <- Filter(function(e) identical(e$exemption_id, exemption_id), cfg$exemptions)
  # An unknown, revoked or expired id licenses nothing.
  if (!length(match)) return(character(0))
  entry <- match[[1]]
  # An id without an in-scope caller licenses nothing. This is what stops an
  # exemption granted for one workstream from being borrowed by another.
  if (is.null(caller_path)) return(character(0))
  normalised <- gsub("\\\\", "/", caller_path)
  scopes <- unlist(entry$scope_paths)
  in_scope <- any(vapply(scopes, function(s) grepl(s, normalised, fixed = TRUE),
                         logical(1)))
  if (!in_scope) return(character(0))
  setdiff(tolower(unlist(entry$terms)), cfg$never)
}

assert_language <- function(text, allow = character(0),
                            exemption_id = NULL, caller_path = NULL) {
  lowered <- tolower(paste(text, collapse = " "))
  licensed <- .licensed_terms(exemption_id, caller_path)
  permitted <- union(allow, licensed)
  hits <- setdiff(FORBIDDEN_TERMS[vapply(FORBIDDEN_TERMS,
                                         function(t) grepl(t, lowered, fixed = TRUE),
                                         logical(1))], permitted)
  assert_true(!length(hits),
              paste0("Forbidden ordering vocabulary in output text: ",
                     paste(hits, collapse = ", ")))
  invisible(TRUE)
}

# 7b. The affirmative counterpart. A Resource script calls this to prove it did
# not rely on an exemption, so "no exemption was used" is a checked fact rather
# than an absence of evidence.
assert_no_exemption <- function(text, allow = character(0)) {
  assert_language(text, allow = allow, exemption_id = NULL, caller_path = NULL)
}

# 8. Inclusion test. A system that improves none of these does not belong in the
# Resource paper (PAPER.md:211-214).
INCLUSION_CRITERIA <- c("coverage", "observability", "evidence_interpretation",
                        "cross_assay_transportability", "experiment_routing")

assert_inclusion_criterion <- function(criteria) {
  assert_true(length(criteria) > 0L && all(criteria %in% INCLUSION_CRITERIA),
              paste0("Each system must name at least one inclusion criterion from: ",
                     paste(INCLUSION_CRITERIA, collapse = ", ")))
  invisible(TRUE)
}
