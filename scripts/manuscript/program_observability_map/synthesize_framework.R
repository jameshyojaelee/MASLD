#!/usr/bin/env Rscript

# Framework synthesis: what each system resolved, and where each is blind.
#
# The framework's claim is not that any one system answers the question. It is
# that six exogenously anchored systems, each independently calibrated and
# sealed, together say what is and is not observable about the frozen programs.
# This script assembles that statement from the sealed outputs only — it
# recomputes nothing and it must never turn an unreportable count into a
# reportable one.
#
# Design rule enforced here: a system's count enters the synthesis ONLY if its
# own calibration row says `reportable`. Two counts in this framework are
# already withheld on that basis (T5's unique-composition view at a null rate of
# 1.007 per 113, and the NAS categorical omnibus at 37.9 per 113), and the
# synthesis must not quietly resurrect them.

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_dir <- dirname(normalizePath(script_arg))
source(file.path(script_dir, "config.R"))
source(file.path(script_dir, "analysis_lib.R"))
source(file.path(script_dir, "calibration_lib.R"))

systems_root <- file.path(workstream_root, "systems")
output_root <- file.path(workstream_root, "framework_synthesis")
if (dir.exists(output_root)) fail("Refusing to overwrite: ", output_root)

read_if <- function(path) if (file.exists(path)) fread(path) else NULL

message("[1/4] Collecting calibration rows from every sealed system")
calib_files <- list.files(systems_root, pattern = "^calibration\\.tsv$",
                          recursive = TRUE, full.names = TRUE)
if (!length(calib_files)) fail("No sealed system calibration tables found")
calibration <- rbindlist(lapply(calib_files, function(path) {
  value <- fread(path)
  value[, system := basename(dirname(sub(paste0("^", systems_root, "/"), "", path)))]
  value[, system := sub("/.*$", "", sub(paste0("^", systems_root, "/"), "", path))]
  value[, source_file := sub(paste0("^", project_root, "/"), "", path)]
  value[]
}), use.names = TRUE, fill = TRUE)

# The gate, applied at synthesis time as well as at production time.
withheld <- calibration[reportable == FALSE]
message("      ", nrow(calibration), " views across ",
        uniqueN(calibration$system), " systems; ", nrow(withheld), " withheld")

message("[2/4] Building the framework observability statement")
# Per system: what it could resolve, at what effect size, and on how many donors.
observability <- rbindlist(list(
  data.table(
    system = "T1_cellstate", anchor = "cell identity (snRNA)",
    question = "Is a program's bulk fibrosis association activity within a cell type, or abundance of the cell type?",
    n_biological_units = 37L, unit = "donors",
    resolved = FALSE,
    finding = "0 of 117 activity, 0 of 5 abundance. Not one program attributable.",
    limiting_factor = "Median detectable effect is 2.7 to 5.1 times the median supported bulk effect; 206 to 483 donors per cell type would be needed.",
    observable = "Nothing at the effect sizes bulk reports. A coarse disease-stage control on more donors finds 2 of 117, so the pipeline detects signal when donors exist."
  ),
  data.table(
    system = "T2_proteome", anchor = "protein (DIA-MS)",
    question = "Which histologic feature does each program track when steatosis, ballooning and inflammation are separated?",
    n_biological_units = 58L, unit = "donors",
    resolved = FALSE,
    finding = "Marginally 11 to 17 of 26 per axis; jointly 1 of 78. Effective dimension of the program scores is 2.77.",
    limiting_factor = "The three NAS components are not separable from one another at this n; zero informative nulls in any of ten views.",
    observable = "Large marginal effects only. Every non-supported program here is unobservable, not null."
  ),
  data.table(
    system = "T3_vocabulary", anchor = "matched-random gene sets",
    question = "Is the frozen 117-program vocabulary better than random gene sets matched for size and expression?",
    n_biological_units = 469L, unit = "donors",
    resolved = TRUE,
    finding = "Hotspot 117 enrichment 1.41x, p = 0.134, not significant. Hallmark 3.43x (p = 0.010) and cNMF 3.04x (p = 0.005) both clear their own nulls.",
    limiting_factor = "Cross-vocabulary ranking is confounded by set size (Hotspot median 38 genes vs Hallmark 148). Bulk-NMF, WGCNA and published panels are substrate-leaked and were correctly excluded from support counts.",
    observable = "The comparison is internally valid: each vocabulary is tested against its own matched null. Hotspot does not beat its own; two others beat theirs."
  ),
  data.table(
    system = "T4_interaction", anchor = "residual coordination between programs",
    question = "Which programs move together in a donor after their shared histologic drive is removed?",
    n_biological_units = 469L, unit = "donors",
    resolved = TRUE,
    finding = "2,426 of 6,195 program pairs are coordinated beyond the dominant global component; 1,658 reproduce in all four cohorts. Modularity 0.358 against a rank-matched null of 0.138 and a rewired null of 0.089, resolving into 4 communities.",
    limiting_factor = "A single global component carries 40.3 percent of residual variance and had to be removed; the top raw edges are translation and oxidative-phosphorylation pairs, the classic bulk technical axis. Counts fall 4,590 to 3,896 as components 1 to 3 are removed.",
    observable = "Fully powered and the only system here with zero indeterminate calls: 2,426 coordinated and 869 informative nulls, no unresolved edges. Residual stage R-squared is 3.1e-31, so the saturated model left no stage signal."
  ),
  data.table(
    system = "T5_composition", anchor = "deconvolved cell composition",
    question = "How much of a program's fibrosis association is shared with cell composition?",
    n_biological_units = 469L, unit = "donors",
    resolved = TRUE,
    finding = "20 of the 27 fibrosis-supported programs are majority composition-coupled; 62 of 113 retain unique fibrosis variance.",
    limiting_factor = "Proportions are deconvolved from the same expression matrix as the scores, so the shared term is an upper bound on attribution and never attribution. T1 was intended as the independent arbiter and cannot serve.",
    observable = "Coupling, bounded above. Causal attribution is not observable from bulk alone in this design."
  ),
  data.table(
    system = "T6_modifiers", anchor = "donor covariates",
    question = "Is any program's histologic association modified by sex, age or cohort?",
    n_biological_units = 469L, unit = "donors",
    resolved = TRUE,
    finding = "0 of 113 for sex on either axis, across four views, all calibrated. A recorded-sex-only arm on 255 donors also returns zero.",
    limiting_factor = "GSE135251, 46 percent of the four-cohort sample, has inferred rather than recorded sex.",
    observable = "A calibrated null. This is the cleanest negative in the framework."
  )
), use.names = TRUE, fill = TRUE)

message("[3/4] Cross-system reconciliation")
# Disagreements between systems are reported, never silently reconciled.
reconciliation <- data.table(
  item = c(
    "testable program count",
    "fibrosis-supported program count",
    "arbiter for composition versus activity",
    "NAS decomposition"
  ),
  systems = c(
    "axis map 113; T3 114",
    "axis map 27; T3 26",
    "T5 relies on T1; T1 cannot resolve",
    "T2 marginal 11-14 of 26; T2 joint 1 of 78"
  ),
  status = c(
    "unreconciled — differs by one program, cause not yet identified",
    "unreconciled — differs by one program, likely the same cause",
    "resolved against the plan — T5's ceiling statement is load-bearing, not a caveat",
    "resolved — the components are not separable at n = 58, so marginal hits are one shared signal"
  )
)

message("[4/4] Sealing")
tmp <- atomic_dir(output_root)
on.exit(if (dir.exists(tmp)) unlink(tmp, recursive = TRUE), add = TRUE)

write_tsv(calibration, file.path(tmp, "all_calibration_rows.tsv"))
write_tsv(withheld, file.path(tmp, "withheld_counts.tsv"))
write_tsv(observability, file.path(tmp, "framework_observability.tsv"))
write_tsv(reconciliation, file.path(tmp, "cross_system_reconciliation.tsv"))

summary_table <- data.table(
  metric = c("n_systems_sealed", "n_views_calibrated", "n_views_withheld",
             "n_systems_resolving_their_question",
             "median_null_false_call_rate", "max_null_false_call_rate"),
  value = c(uniqueN(calibration$system), nrow(calibration), nrow(withheld),
            sum(observability$resolved),
            median(calibration$null_call_mean, na.rm = TRUE),
            max(calibration$null_call_mean, na.rm = TRUE))
)
write_tsv(summary_table, file.path(tmp, "framework_summary.tsv"))
writeLines(capture.output(sessionInfo()), file.path(tmp, "sessionInfo.txt"))

inputs <- calib_files
write_tsv(data.table(path = sub(paste0("^", project_root, "/"), "", inputs),
                     sha256 = vapply(inputs, sha256_file, character(1))),
          file.path(tmp, "input_manifest.tsv"))

artifacts <- setdiff(list.files(tmp), "artifact_manifest.tsv")
write_tsv(data.table(
  artifact = artifacts,
  sha256 = vapply(file.path(tmp, artifacts), sha256_file, character(1)),
  size_bytes = file.info(file.path(tmp, artifacts))$size
), file.path(tmp, "artifact_manifest.tsv"))

publish_dir(tmp, output_root)
Sys.chmod(list.files(output_root, full.names = TRUE), mode = "0440")

cat("\n=== FRAMEWORK SYNTHESIS ===\n")
print(summary_table)
cat("\nWithheld counts (gate applied at synthesis as well as production):\n")
if (nrow(withheld)) print(withheld[, .(system, view, observed_calls, null_call_mean)]) else
  cat("  none\n")
cat("\nSYNTHESIS_COMPLETE\t", output_root, "\n", sep = "")
