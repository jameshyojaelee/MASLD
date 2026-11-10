#!/usr/bin/env Rscript
# Four primary voom fits plus explicit empty macrophage-universe checks.
if (!nzchar(Sys.getenv("SLURM_JOB_ID"))) stop("Requires a Slurm compute allocation")
arguments <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", arguments[grepl("^--file=", arguments)])
source(file.path(dirname(normalizePath(script)), "..", "scripts", "compare_da.R"))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L || length(args) > 3L) stop("Usage: check_real_baseline.R CONFIG NEW_OUTPUT [--profile]")
cfg <- fromJSON(args[1L], simplifyVector = FALSE)
out <- args[2L]
if (file.exists(out)) stop("Refusing to overwrite output")
dir.create(out, recursive = TRUE)
cohorts <- names(cfg$cohorts)
checks <- timings <- list()
first_fit <- TRUE
time_step <- function(label, expression) {
  message(format(Sys.time(), "%FT%T%z"), " START ", label)
  started <- proc.time()[["elapsed"]]
  value <- force(expression)
  elapsed <- proc.time()[["elapsed"]] - started
  message(format(Sys.time(), "%FT%T%z"), " DONE ", label, " elapsed_seconds=", elapsed)
  timings[[length(timings) + 1L]] <<- data.table(step = label, elapsed_seconds = elapsed)
  fwrite(rbindlist(timings), file.path(out, "timings.tsv"), sep = "\t")
  value
}
for (lineage in unlist(cfg$lineages)) {
  inputs <- setNames(lapply(cohorts, function(cohort) time_step(paste("load", cohort, lineage),
    load_da_counts(cfg$counts_root, cohort, lineage))), cohorts)
  retained <- lapply(cohorts, function(cohort) time_step(paste("filter", cohort, lineage),
    colnames(inputs[[cohort]]$counts)[filter_da_peaks(inputs[[cohort]]$counts, "primary")]))
  common <- Reduce(intersect, retained)
  message("COMMON ", lineage, " n=", length(common))
  for (cohort in cohorts) {
    if (!length(common)) {
      reference <- fread(file.path(cfg$corrected_baseline_root, "da",
        paste0("cohort_primary_", cohort, "_", lineage, ".tsv")))
      require_true(nrow(reference) == 0L, "Unexpected nonempty baseline")
      checks[[length(checks) + 1L]] <- data.table(cohort = cohort, lineage = lineage,
        status = "pass_empty_untestable_universe", n_peaks = 0L)
    } else {
      if (first_fit && "--profile" %in% args) Rprof(file.path(out, "first_voom.Rprof"), interval = 0.02)
      fit <- time_step(paste("fit", cohort, lineage), fit_da(inputs[[cohort]]$counts[, common, drop = FALSE],
        inputs[[cohort]]$donors, "voom_qw"))
      if (first_fit && "--profile" %in% args) {
        Rprof(NULL)
        capture.output(summaryRprof(file.path(out, "first_voom.Rprof")), file = file.path(out, "first_voom_profile.txt"))
      }
      first_fit <- FALSE
      checks[[length(checks) + 1L]] <- time_step(paste("compare", cohort, lineage),
        baseline_reproduction(fit$results, cfg$corrected_baseline_root, cohort, lineage,
          as.numeric(cfg$da$baseline_tolerance)))
    }
    fwrite(rbindlist(checks, fill = TRUE), file.path(out, "baseline_reproduction.tsv"), sep = "\t")
  }
}
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
message("PASS: four primary voom fits and both empty macrophage baseline universes")
