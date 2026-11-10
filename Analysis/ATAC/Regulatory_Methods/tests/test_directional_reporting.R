#!/usr/bin/env Rscript
if (!nzchar(Sys.getenv("SLURM_JOB_ID"))) stop("Requires a Slurm compute allocation")
arguments <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", arguments[grepl("^--file=", arguments)])
source(file.path(dirname(normalizePath(script)), "..", "scripts", "compare_da.R"))
fixture <- data.table(cohort = "cross_cohort", representation = "fragment", filter = "primary", method = "voom_qw",
  lineage = c("hepatocyte", "stellate", "macrophage"), status = c("batch_restricted_monte_carlo", "not_run", "untestable"),
  directional_pvalue = c(.01, NA_real_, 1), directional_reps = c(1000L, 0L, 0L))
actual <- adjust_directional_rows(fixture)
stopifnot(identical(actual$directional_qvalue, c(.03, NA_real_, NA_real_)),
  identical(actual$directional_pvalue, c(.01, NA_real_, NA_real_)),
  identical(actual$directional_tested, c(TRUE, FALSE, FALSE)), all(actual$directional_family_size == 3L))
fixture[, `:=`(status = c("not_run", "no_reference_discoveries", "untestable_no_within_batch_label_exchange"),
  directional_pvalue = c(NA_real_, 1, NA_real_), directional_reps = 0L)]
actual <- adjust_directional_rows(fixture)
stopifnot(all(is.na(actual$directional_qvalue)), all(is.na(actual$directional_pvalue)), !any(actual$directional_tested))
cat("PASS: fixed three-lineage BH with unavailable hypotheses internally retained and externally missing\n")

source(file.path(dirname(normalizePath(script)), "..", "scripts", "compare_da_sensitivity.R"))
full <- data.table(peak_coordinate = c("a", "b", "c"), logFC = c(1, 2, 3), qvalue = c(.01, .2, .03))
down <- data.table(peak_coordinate = c("b", "c", "d"), logFC = c(2.1, 3.1, 4.1), qvalue = c(.04, .5, .01))
mapping <- data.table(peak_coordinate = c("a", "b", "c", "d"), gene_id = c("ENSG1", "ENSG1", "ENSG2", "ENSG3"))
comparison <- compare_effect_arm(full, down, mapping)
stopifnot(comparison$n_common_peaks == 2L, comparison$n_full_only_peaks == 1L,
  comparison$n_downsample_only_peaks == 1L, comparison$bh_universe_changed,
  abs(comparison$median_absolute_log2_effect_change - .1) < 1e-12,
  comparison$n_bh_q05_full_original_family == 2L, comparison$n_bh_q05_full_on_common_peaks == 1L,
  comparison$n_mapped_genes_full == 2L, comparison$n_mapped_genes_downsample == 3L,
  comparison$n_mapped_genes_common_peaks == 2L, comparison$n_genes_shared_across_arms == 2L)
missing_map <- compare_effect_arm(full, down)
stopifnot(is.na(missing_map$n_mapped_genes_full), missing_map$gene_coverage_status == "unavailable_without_explicit_peak_gene_map")
cat("PASS: common-peak effects, changed equal-sized BH universes, original q families and explicit gene-map coverage\n")

# Exercise the actual compact-file interface, including donor roster checks.
temporary <- tempfile("da_sensitivity_fixture_")
dir.create(temporary)
for (name in c("full", "down")) {
  directory <- file.path(temporary, name)
  dir.create(directory)
  peaks <- copy(if (name == "full") full else down)
  peaks[, `:=`(cohort = "GSE244832", lineage = "hepatocyte", representation = "fragment", filter = "primary",
    method = "voom_qw", effect_units = "log2_normalized_fragment_count")]
  fwrite(peaks, file.path(directory, "peak_results.tsv.gz"), sep = "\t", compress = "gzip")
  arm <- unique(peaks[, ..SENSITIVITY_KEYS])
  qc <- rbindlist(list(cbind(arm, data.table(donor_id = "D01", condition = "NORMAL")),
    cbind(arm, data.table(donor_id = "D02", condition = "MASH"))))
  fwrite(qc, file.path(directory, "model_qc.tsv.gz"), sep = "\t", compress = "gzip")
  fwrite(cbind(arm, data.table(status = "testable", n_tested_peaks = 3L)), file.path(directory, "cohort_summary.tsv"), sep = "\t")
}
run_sensitivity_comparison(list(full = file.path(temporary, "full"), downsample = file.path(temporary, "down"),
  output = file.path(temporary, "comparison"), peak_gene_map = NULL, gene_map_label = NULL))
exported <- fread(file.path(temporary, "comparison", "sensitivity_comparison.tsv"))
stopifnot(nrow(exported) == 1L, exported$n_common_peaks == 2L,
  exported$donor_roster_status == "identical_donor_conditions", exported$bh_universe_changed)
cat("PASS: sensitivity comparison compact-file interface and unchanged donor roster\n")
