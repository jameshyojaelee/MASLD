#!/usr/bin/env Rscript
# sex_v3/v6_manifest_emit.R
# ---------------------------------------------------------------------------
# Final manifest mutation + v6_summary.md generator.
# Runs afterany on all v6 pillars. Reads the scaffold pipeline_manifest_v6.json,
# hashes all expected output csvs, counts n_reps_used where applicable, and
# writes outputs/v6_summary.md for the reviewer-response appendix.
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(digest)
})

BASE  <- Sys.getenv("MASLD_PROJECT_ROOT",
                    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
SEXV3 <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration",
                   "results/integration/sex_v3")
MANIFEST <- file.path(SEXV3, "pipeline_manifest_v6.json")
SUMMARY  <- file.path(SEXV3, "outputs/v6_summary.md")
dir.create(dirname(SUMMARY), showWarnings = FALSE, recursive = TRUE)

if (!file.exists(MANIFEST)) {
  cat("WARNING: manifest scaffold missing — emitting bare summary\n")
  base_manifest <- list(created_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ"),
                        job_ids = list(), status = "no_scaffold")
} else {
  base_manifest <- fromJSON(MANIFEST, simplifyVector = FALSE)
}

expected_outputs <- list(
  P1 = "interaction_classifier_v5_random.csv",
  P2 = "perm_fdr_v6.csv",
  P3 = "selection_bias_simulation.csv",
  P4 = "bootstrap_stability_v6.csv",
  P5 = "calibration_metrics_v6.csv",
  P6 = "per_cohort_forest_data_v6.csv",
  P7 = "sensitivity_sweep_v6.csv",
  P8A = "male_biased_characterization.csv",
  P8B = "loco_concordance_per_class_v6.csv",
  P9_FINAL = "sex_deg_classification_v6.csv"
)

hash_file <- function(p) {
  if (!file.exists(p)) return(NA_character_)
  digest::digest(p, file = TRUE, algo = "md5")
}

audit <- lapply(names(expected_outputs), function(k) {
  fn <- file.path(SEXV3, expected_outputs[[k]])
  exists_ok <- file.exists(fn)
  size <- if (exists_ok) file.info(fn)$size else NA_integer_
  md5  <- hash_file(fn)
  nrows <- if (exists_ok && grepl("\\.csv$", fn)) {
    tryCatch(nrow(fread(fn, nrows = -1L, showProgress = FALSE)),
             error = function(e) NA_integer_)
  } else NA_integer_
  list(pillar = k,
       expected_file = expected_outputs[[k]],
       exists = exists_ok,
       size_bytes = size,
       md5 = md5,
       nrows = nrows)
})

base_manifest$audit <- audit
base_manifest$completed_utc <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ")
all_pillars_ok <- all(vapply(audit, function(x) isTRUE(x$exists), logical(1)))
base_manifest$all_pillars_status <- if (all_pillars_ok) "ok" else "incomplete"

writeLines(toJSON(base_manifest, pretty = TRUE, auto_unbox = TRUE, null = "null"),
           MANIFEST)
cat("Updated manifest:", MANIFEST, "\n")

# ---------------------------------------------------------------------------
# v6_summary.md
# ---------------------------------------------------------------------------
md <- c(
  "# sex_v6 pipeline summary",
  "",
  paste0("Completed: ", base_manifest$completed_utc),
  paste0("Status: **", base_manifest$all_pillars_status, "**"),
  "",
  "## Output audit",
  "",
  "| Pillar | File | Exists | Size (B) | nRows | md5 |",
  "|---|---|---|---|---|---|"
)
for (a in audit) {
  md <- c(md, sprintf("| %s | `%s` | %s | %s | %s | `%s` |",
                      a$pillar,
                      a$expected_file,
                      if (isTRUE(a$exists)) "yes" else "MISSING",
                      ifelse(is.na(a$size_bytes), "—", format(a$size_bytes, big.mark = ",")),
                      ifelse(is.na(a$nrows), "—", a$nrows),
                      ifelse(is.na(a$md5), "—", substr(a$md5, 1, 10))))
}
md <- c(md, "",
        "## Job IDs",
        "",
        "| Pillar | JID |",
        "|---|---|")
for (k in names(expected_outputs)) {
  jid <- base_manifest$job_ids[[k]]
  if (is.null(jid)) jid <- "—"
  md <- c(md, sprintf("| %s | %s |", k, jid))
}

# Pull headline numbers if consensus exists
canon <- file.path(SEXV3, "sex_deg_classification_v6.csv")
if (file.exists(canon)) {
  cd <- tryCatch(fread(canon), error = function(e) NULL)
  if (!is.null(cd) && "class_v6_consensus" %in% names(cd)) {
    md <- c(md, "",
            "## class_v6_consensus tally",
            "",
            "| class | N |",
            "|---|---|")
    tt <- cd[, .N, by = class_v6_consensus][order(-N)]
    for (i in seq_len(nrow(tt))) {
      md <- c(md, sprintf("| %s | %d |", tt$class_v6_consensus[i], tt$N[i]))
    }
  }
}

writeLines(md, SUMMARY)
cat("Wrote v6_summary.md:", SUMMARY, "\n")
cat("\nDone:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
