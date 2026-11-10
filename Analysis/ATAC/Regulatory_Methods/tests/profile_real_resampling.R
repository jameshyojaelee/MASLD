#!/usr/bin/env Rscript
# Bounded timing receipt only; two draws are not a calibration result.
if (!nzchar(Sys.getenv("SLURM_JOB_ID"))) stop("Requires a Slurm compute allocation")
arguments <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", arguments[grepl("^--file=", arguments)])
source(file.path(dirname(normalizePath(script)), "..", "scripts", "compare_da.R"))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: profile_real_resampling.R CONFIG INITIAL_DA NEW_OUTPUT")
cfg <- fromJSON(args[1L], simplifyVector = FALSE)
out <- args[3L]
if (file.exists(out)) stop("Refusing to overwrite output")
dir.create(out, recursive = TRUE)
data.table::setDTthreads(1L)
lineage <- "stellate"
cohorts <- names(cfg$cohorts)
timings <- nulls <- directions <- list()
time_step <- function(label, expression) {
  message(format(Sys.time(), "%FT%T%z"), " START ", label)
  started <- proc.time()
  value <- force(expression)
  elapsed <- proc.time() - started
  row <- data.table(step = label, elapsed_seconds = unname(elapsed["elapsed"]),
    user_seconds = unname(elapsed["user.self"]), system_seconds = unname(elapsed["sys.self"]))
  timings[[length(timings) + 1L]] <<- row
  fwrite(rbindlist(timings), file.path(out, "timings.tsv"), sep = "\t")
  message(format(Sys.time(), "%FT%T%z"), " DONE ", label, " elapsed_seconds=", row$elapsed_seconds)
  value
}
inputs <- setNames(lapply(cohorts, function(cohort) time_step(paste("load", cohort),
  load_da_counts(cfg$counts_root, cohort, lineage))), cohorts)
retained <- lapply(inputs, function(x) colnames(x$counts)[filter_da_peaks(x$counts, "prevalence")])
common <- Reduce(intersect, retained)
x <- inputs[["GSE244832"]]
x$donors <- join_donor_metadata(x$donors, cfg$cohort_metadata$GSE244832)
x$donors <- join_donor_metadata(x$donors, cfg$cohort_covariate_metadata$GSE244832)
counts <- x$counts[, common, drop = FALSE]
reference <- time_step("read fixed initial reference", fread(file.path(args[2L], "replication.tsv.gz"),
  select = c("peak_coordinate", "logFC_gse244832", "logFC_gse281367", "qvalue_gse281367", "lineage", "filter", "method")))
reference <- reference[lineage == "stellate" & filter == "prevalence"]
blocks <- x$donors[, .(n_normal = sum(condition == "NORMAL"), n_mash = sum(condition == "MASH")), by = .(batch_jb, batch_mm)]
blocks[, distinct_label_assignments := choose(n_normal + n_mash, n_normal)]
fwrite(blocks, file.path(out, "batch_strata.tsv"), sep = "\t")
space <- prod(blocks$distinct_label_assignments)
seed <- arm_seed(as.integer(cfg$seed), "GSE244832", lineage, "fragment", "prevalence")
for (method in c("voom_qw", "edger_ql")) {
  chosen_method <- method
  ref <- reference[method == chosen_method]
  require_true(setequal(ref$peak_coordinate, common), "Fixed reference and profiled peak universes differ")
  discoveries <- setNames(ref[qvalue_gse281367 < .05]$logFC_gse281367, ref[qvalue_gse281367 < .05]$peak_coordinate)
  observed <- sum(sign(ref$logFC_gse244832[match(names(discoveries), ref$peak_coordinate)]) == sign(discoveries))
  null <- time_step(paste("NB null two draws", method), null_da(counts, x$donors, method,
    reps = 2L, max_peaks = as.integer(cfg$da$null_max_peaks), seed = seed, workers = 1L))
  null[, `:=`(method = method, timing_only = TRUE)]
  nulls[[length(nulls) + 1L]] <- null
  fwrite(rbindlist(nulls, fill = TRUE), file.path(out, "null_profile.tsv"), sep = "\t")
  direction <- time_step(paste("directional two draws", method), directional_da(counts, x$donors, method,
    discoveries, observed, reps = 2L, seed = seed + 10000L, workers = 1L))
  directions[[length(directions) + 1L]] <- data.table(method = method, n_tested_peaks = length(common),
    n_reference_discoveries = length(discoveries), n_observed_concordant = observed,
    profile_directional_pvalue = direction$p, profile_draws = direction$n, profile_status = direction$status,
    timing_only = TRUE)
  fwrite(rbindlist(directions), file.path(out, "directional_profile.tsv"), sep = "\t")
}
write_json(list(status = "timing_only_not_a_scientific_calibration_result", cohort = "GSE244832", lineage = lineage,
  representation = "fragment", filter = "prevalence", n_donors = nrow(x$donors), n_tested_peaks = ncol(counts),
  workers = 1L, null_draws_per_method = 2L, directional_draws_per_method = 2L,
  seed = seed, directional_seed = seed + 10000L, distinct_batch_restricted_assignments = space,
  batch_policy = "joint_recorded_batch_jb_batch_mm; one_to_one_in_current_metadata",
  allocation = Sys.getenv(c("SLURM_JOB_ID", "SLURM_CPUS_PER_TASK", "SLURM_JOB_NODELIST")),
  hostname = Sys.info()[["nodename"]], threads = Sys.getenv(c("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS"))),
  file.path(out, "receipt.json"), auto_unbox = TRUE, pretty = TRUE)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
message("COMPLETE timing receipt; do not interpret two-draw p-values")
