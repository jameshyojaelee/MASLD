#!/usr/bin/env Rscript
# 01: record what this run read, before it reads anything for science.
#
# Two things are pinned here. The sha256 of every input and of every script, so
# a later disagreement about a number can be traced to a changed file rather
# than argued about; and the accession overlap between the discovery cohorts
# and every evaluation cohort, because the transportability claim is only worth
# making if the evaluation sets are not the discovery set under another name.

source(file.path(Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"),
  "scripts/analysis/cross_assay_modules/lib_cross_assay_modules.R"))

contract <- cam_contract()
out <- cam_dir("inputs")

cam_say("hashing contract inputs")
keys <- names(contract$inputs)
rows <- rbindlist(lapply(keys, function(k) {
  path <- cam_input(k, contract, must_exist = FALSE)
  exists <- file.exists(path)
  is_dir <- exists && dir.exists(path)
  # A directory input is pinned by the hashes of the files inside it, so that
  # swapping one member file is as visible as swapping a named input.
  sha <- if (!exists) {
    NA_character_
  } else if (is_dir) {
    members <- sort(list.files(path, recursive = TRUE, full.names = TRUE))
    members <- members[!dir.exists(members)]
    if (!length(members)) NA_character_ else
      substr(digest::digest(paste(vapply(members, cam_sha256, character(1)),
                                  collapse = "\n"), algo = "sha256", serialize = FALSE), 1, 64)
  } else if (file.info(path)$size < 3e9) {
    cam_sha256(path)
  } else {
    NA_character_
  }
  data.table(
    input_key = k,
    path = path,
    exists = exists,
    is_directory = is_dir,
    size_bytes = if (exists && !is_dir) file.info(path)$size else NA_real_,
    sha256 = sha
  )
}))
cam_assert(all(rows$exists), paste0("Missing inputs: ",
  paste(rows[exists == FALSE, input_key], collapse = ", ")))
cam_write_tsv(rows, file.path(out, "input_hashes.tsv"))

cam_say("hashing this analysis's own code")
code <- list.files(cam_script_dir(), pattern = "\\.(R|py|sh|json)$", full.names = TRUE)
cam_write_tsv(
  data.table(script = basename(code), sha256 = vapply(code, cam_sha256, character(1))),
  file.path(out, "code_hashes.tsv")
)

# Cohort independence. Discovery is the five bulk cohorts; every evaluation
# cohort must be a different accession. This does not establish donor-level
# disjointness, which no public key supports, and the record says so.
cam_say("recording accession separation")
discovery <- contract$discovery_label$cohorts
evaluation <- c("GSE268273", "GSE276114", "PXD051911", "GSE267145", "GSE296875",
                "GSE192741", "Vu2025", "GSE312698")
atlas <- c("GSE244832", "GSE202379", "GSE185477", "GSE136103", "GSE189600",
           "GSE174748", "GSE192740")
overlap <- data.table(
  evaluation_accession = c(evaluation, atlas),
  role = c(rep("evaluation_assay", length(evaluation)), rep("snrna_atlas", length(atlas))),
  shares_accession_with_discovery = c(evaluation, atlas) %in% discovery,
  donor_level_disjointness = "unverifiable_no_public_donor_key"
)
cam_write_tsv(overlap, file.path(out, "accession_separation.tsv"))
cam_assert(!any(overlap$shares_accession_with_discovery),
           "An evaluation or atlas accession is also a discovery cohort")

cam_write_json(list(
  contract_id = contract$contract_id,
  seed = contract$seed,
  n_inputs = nrow(rows),
  n_scripts = length(code),
  discovery_cohorts = discovery,
  frozen_at_utc = format(Sys.time(), tz = "UTC", "%Y-%m-%dT%H:%M:%SZ")
), file.path(out, "freeze_summary.json"))

cam_session_info(file.path(cam_out_root(), "sessionInfo.txt"))
writeLines("inputs frozen", file.path(out, "READY"))
cam_say("01 complete: ", nrow(rows), " inputs, ", length(code), " scripts hashed")
