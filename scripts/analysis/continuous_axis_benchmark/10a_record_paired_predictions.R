#!/usr/bin/env Rscript
# 10a: freeze the paired-biopsy predictions BEFORE the paired counts are loaded.
#
# GSE193066 was a sealed holdout and is now unsealed by user decision, which costs
# it as a confirmatory holdout. The substitute is this: write the predicted
# directions, hash them, make the file immutable, and have 10b refuse to run unless
# the hash still matches. git cannot serve this purpose here because .gitignore
# line 8 is a blanket '*', so the repo's own idiom -- content hashes in mode-440
# files, as freeze_contract.py uses -- is the durable mechanism.
#
# This script must not read a single GSE193066 count.

suppressPackageStartupMessages({ library(data.table); library(jsonlite) })
source(file.path(Sys.getenv("MASLD_PROJECT_ROOT"),
                 "scripts/analysis/continuous_axis_benchmark/lib_axis_common.R"))

OUT <- Sys.getenv("CAB_OUT_ROOT", "")
assert_true(nzchar(OUT), "CAB_OUT_ROOT is unset")
dir.create(file.path(OUT, "paired"), recursive = TRUE, showWarnings = FALSE)

# Refuse if any downstream output already exists: predictions recorded after the
# fact are not predictions.
existing <- list.files(file.path(OUT, "paired"), pattern = "^10b_", full.names = TRUE)
assert_true(!length(existing),
            paste0("10b outputs already exist; refusing to write predictions after the fact: ",
                   paste(basename(existing), collapse = ", ")))

pre <- read_prespec()
p <- pre$paired_gse193066$predictions

# Cross-sectional anchors these predictions are derived from, recorded so the
# prediction cannot later be described as vague.
q2 <- file.path(OUT, "q2", "M_gene_histology_rank1_verdict.tsv")
anchors <- if (file.exists(q2)) {
  v <- fread(q2)
  list(cross_sectional_pve1 = v$observed_pve1[1],
       cross_sectional_rank1_sufficient = v$rank1_sufficient[1],
       worst_window = v$worst_window_by_resid[1],
       frac_features_breaking_rank1 = v$frac_features_breaking_rank1[1])
} else list(note = "Q2 verdict not available at prediction time")

predictions <- list(
  workstream = "CONTINUOUS-AXIS-BENCHMARK-v1",
  recorded_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  dataset = "GSE193066",
  unit = "biological donor with two biopsies",
  expected_pairs = pre$paired_gse193066$expected_pairs,
  counts_loaded_at_prediction_time = FALSE,
  cross_sectional_anchors = anchors,
  predictions = list(
    P1 = list(statement = p$P1$statement, test = "one-sided Spearman",
              predicted_direction = "positive", powered = TRUE,
              falsified_if = "rho <= 0 or the one-sided p exceeds 0.05"),
    P2 = list(statement = p$P2$statement,
              threshold_disattenuated_cosine = p$P2$threshold_disattenuated_cosine,
              predicted_direction = "cosine > 0.5 between the within-person delta vector and the cross-sectional F0->F4 direction",
              powered = TRUE,
              falsified_if = "cosine <= 0.5, or the 54 delta vectors are not approximately rank 1 against their noise-model null"),
    P3 = list(statement = p$P3$statement, powered = FALSE,
              predicted_direction = "NULL interaction under amplification; non-null under ordered sequence",
              falsified_if = "a significant baseline-axis-by-programme interaction survives donor-level correction",
              caveat = "54 donors is underpowered for this interaction; reported descriptively with a simulated minimum detectable effect"),
    P4 = list(statement = p$P4$statement, powered = FALSE,
              predicted_direction = "F0->F1 transitioners show the largest orthogonal component of delta",
              falsified_if = "orthogonal component is not elevated in F0->F1 transitioners")
  ),
  interpretation_rule = paste(
    "P1 and P2 are the powered tests. If P2 holds, within-person change is parallel to",
    "the cross-sectional axis, which is the strongest available evidence that the",
    "cross-sectional amplification structure reflects change within people rather than",
    "differences between them. If P2 fails, the cross-sectional structure does not",
    "transport to within-person change and must be described as a between-donor",
    "property only.")
)

path <- file.path(OUT, "paired", "predictions.json")
write_json_once(predictions, path)
Sys.chmod(path, "0440")
sha <- sha256_file(path)
writeLines(c(sha, basename(path)), file.path(OUT, "paired", "predictions.sha256"))
Sys.chmod(file.path(OUT, "paired", "predictions.sha256"), "0440")

cat("\nPAIRED PREDICTIONS RECORDED AND SEALED\n")
cat("path   : ", path, "\n", sep = "")
cat("sha256 : ", sha, "\n", sep = "")
cat("\n10b will refuse to run unless this hash still matches.\n")
log_step("PAIRED_PREDICTIONS_COMPLETE")
