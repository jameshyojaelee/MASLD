#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
script <- commandArgs(trailingOnly = TRUE)[1L]
source(script)
fixture <- tempfile("atac_render_")
dir.create(fixture)
da <- file.path(fixture, "da")
dir.create(da)
writeLines("{}", file.path(da, "run_summary.json"))
summary <- CJ(cohort = c("GSE244832", "GSE281367"), lineage = "hepatocyte",
  representation = "fragment", filter = "primary", method = c("voom_qw", "edger_ql"))
summary[, n_tested_peaks := 5L]
fwrite(summary, file.path(da, "cohort_summary.tsv"), sep = "\t")
null <- copy(summary)
null[, `:=`(n_nominal_p05 = 0L, n_bh_q05 = 0L, null_model = "fitted_intercept_NB_global_null_fixed_features")]
fwrite(null, file.path(da, "null_diagnostics.tsv"), sep = "\t")
replication <- CJ(peak_coordinate = paste0("chr1:", 1:5, "-100"), method = c("voom_qw", "edger_ql"))
replication[, `:=`(lineage = "hepatocyte", representation = "fragment", filter = "primary",
  logFC_gse244832 = rep(c(-1, 0, 1, 2, 3), each = 2), logFC_gse281367 = rep(c(-2, 0, 1, 3, 4), each = 2),
  qvalue_gse244832 = rep(c(.01, .1, .2, .04, .9), each = 2), qvalue_gse281367 = .5)]
fwrite(replication, file.path(da, "replication.tsv.gz"), sep = "\t", compress = "gzip")
batch <- replication[, .(peak_coordinate, method, lineage, representation, filter, adjusted_logFC = -logFC_gse244832)]
batch[, `:=`(cohort = "GSE244832", covariate = "batch_mm", status = "testable", adjusted_qvalue = 1)]
fwrite(batch, file.path(da, "covariate_sensitivity.tsv.gz"), sep = "\t", compress = "gzip")
out <- file.path(fixture, "rendered")
render_results(list(da = da, output = out))
manifest <- fread(file.path(out, "render_manifest.tsv"))
stopifnot(sum(manifest$status == "rendered") == 9L,
  any(manifest$reason == "missing_table_or_not_run"),
  any(manifest$reason == "no_network_directory_supplied"))
pdfs <- list.files(out, pattern = "\\.pdf$", full.names = TRUE)
stopifnot(length(pdfs) == 9L, all(file.info(pdfs)$size > 1000),
  sum(grepl("da_cohort_effects", manifest$panel)) == 2L,
  sum(grepl("da_batch_effects", manifest$panel)) == 2L)
stopifnot(inherits(try(render_results(list(da = da, output = out)), silent = TRUE), "try-error"))
network <- file.path(fixture, "network")
dir.create(network)
writeLines('{"status":"completed"}', file.path(network, "completion.json"))
writeLines('{"pilot":true,"max_genes":5,"analysis_scope":"technical_pilot"}', file.path(network, "run_settings.json"))
fwrite(data.table(cohort = "test", lineage = "hepatocyte", endpoint = "promoter_accessibility",
  method = "ridge", baseline = "distance", improvement = 0.1, ci_low = -0.1, ci_high = 0.2,
  n_donors = 6L, status = "evaluated"), file.path(network, "comparisons.tsv"), sep = "\t")
dir.create(file.path(network, "test_hepatocyte"))
fwrite(data.table(outer_fold = 0L, gene_id = "TEST1", status = "evaluated", n_candidates = 5L),
  file.path(network, "test_hepatocyte", "coverage.tsv"), sep = "\t")
out_network <- file.path(fixture, "rendered_network")
render_results(list(da = da, network = network, output = out_network))
stopifnot(length(list.files(out_network, pattern = "\\.pdf$")) == 12L)
if (nzchar(Sys.which("pdftotext"))) {
  for (panel in list.files(out_network, pattern = "^network_.*\\.pdf$", full.names = TRUE)) {
    text <- paste(system2("pdftotext", c(shQuote(panel), "-"), stdout = TRUE), collapse = " ")
    stopifnot(grepl("Technical pilot; maximum 5 genes", text, fixed = TRUE))
  }
}
sample_fixture <- data.table(peak_coordinate = sprintf("p%08d", 1:5001), x = 1:5001, y = 1:5001)
omitted <- setdiff(1:5001, unique(round(seq(1, 5001, length.out = 5000))))
sample_fixture$y[omitted] <- -1000000
stats <- prepare_effect_scatter(sample_fixture, "x", "y")
stopifnot(stats$n == 5001L, nrow(stats$shown) == 5000L,
  identical(stats$rho, cor(sample_fixture$x, sample_fixture$y, method = "spearman")),
  stats$rho != cor(stats$shown$x, stats$shown$y, method = "spearman"))
message("PASS: source-only rendering, missing analyses skipped, DA/network PDFs, overwrite refused")
message("PASS: cross-cohort and batch panels; all-peak statistics independent of fixed display sample")
message("Fixture: ", fixture)
