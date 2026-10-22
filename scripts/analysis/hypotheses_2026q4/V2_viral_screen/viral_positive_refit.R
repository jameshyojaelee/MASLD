#!/usr/bin/env Rscript
# V2 sensitivity refit: do persons whose liver RNA-seq is positive for target class B, target class C or target class E
# change the disease calls of the 774-unit primary_source_controls model?
# (candidate bulk-corrected-controls-20261008T162308Z, unadopted).
#
# Model, filter, normalization, quality weights and tests are those of
# scripts/manuscript/resource_f_five/refit_corrected_controls.R, as copied in
# scripts/analysis/bulk_refit_followup/loco_refit.R (integrate_counts, fit_disease_model and
# deg_counts below are copies of that file). The unit is the person (774 units, one library
# each). filterByExpr(~ 0 + group_binary), RLE, voomWithQualityWeights,
# ~ dataset + inferred_sex (+ screen_positive) + group_binary, eBayes, TREAT lfc 0.25, ashr
# normal mixture; canonical call = BH padj < 0.05 and |logFC| > 0.5.
#
# Environment:
#   run mode            none | exclude | covariate
#                         none       all 774 units, model unchanged (must reproduce the candidate)
#                         exclude    target-positive persons removed before filtering and RLE
#                         covariate  all 774 units, model gains a 0/1 screen_positive term
#   VIRAL_PERSON_CALLS  person_calls.tsv written by v2_person_rates.py; a person is screen
#                       positive when basis == person_roster774, target in target class B/target class C/target class E and
#                       call == positive (flagged and indeterminate calls are not used)
#   VIRAL_OUTPUT_ROOT   must not exist

suppressPackageStartupMessages({
  library(ashr)
  library(data.table)
  library(edgeR)
  library(jsonlite)
  library(limma)
  library(yaml)
})

options(digits = 17, scipen = 999)
set.seed(20260923)
fail <- function(...) stop(..., call. = FALSE)

PROJECT_ROOT <- normalizePath(Sys.getenv(
  "MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"), mustWork = TRUE)
MODE <- Sys.getenv("VIRAL_MODE", "")
PERSON_CALLS <- Sys.getenv("VIRAL_PERSON_CALLS", "")
OUT_ROOT <- Sys.getenv("VIRAL_OUTPUT_ROOT", "")
if (!(MODE %in% c("none", "exclude", "covariate"))) fail("VIRAL_MODE must be none, exclude or covariate")
if (!nzchar(PERSON_CALLS) || !file.exists(PERSON_CALLS)) fail("VIRAL_PERSON_CALLS is required and must exist")
if (!nzchar(OUT_ROOT)) fail("VIRAL_OUTPUT_ROOT is required")
if (file.exists(OUT_ROOT)) fail("VIRAL_OUTPUT_ROOT must not exist: ", OUT_ROOT)
if (!dir.exists(dirname(OUT_ROOT))) fail("Output parent is absent: ", dirname(OUT_ROOT))
TREAT_LFC <- 0.25
FULL_ROOT <- file.path(PROJECT_ROOT,
  "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z")
ARM_DIR <- file.path(FULL_ROOT, "primary_source_controls")
V6 <- file.path(PROJECT_ROOT,
  "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10")
INPUTS <- c(
  counts_raw = file.path(PROJECT_ROOT,
    "figures/candidates/resource-corrections-20261005-v1/verified_counts_r8_20261007/assembled_774_primary/counts_raw.rds"),
  meta_matched = file.path(ARM_DIR, "meta_matched.rds"),
  model_design = file.path(ARM_DIR, "model_design.tsv"),
  full_deg = file.path(ARM_DIR, "deg_results.csv"),
  config = file.path(V6, "inputs/BULK-F-FIVE/frozen_model_inputs/human_datasets.yaml"),
  gene_metadata = file.path(V6, "inputs/BULK-F-FIVE/frozen_model_inputs/gencode_v49_gene_metadata.tsv.gz"),
  qualified_input_helper = file.path(PROJECT_ROOT, "scripts/manuscript/resource_f_five/qualified_bulk_inputs.R"),
  person_calls = PERSON_CALLS
)
PINNED <- c(
  counts_raw = "f931e340546d679dc69dd0d6f712f4e73b1d805e80c741b22d3ed25498c68675",
  config = "dc64f6f051022404689abcfeb4c71001f135123c7af65c5bc8d88dd663541f6a",
  gene_metadata = "49dcaae77323ffb20f5fd8399bb2899cd553e3634755a37b60a133d62549c189"
)
sha256_file <- function(path) {
  value <- system2("sha256sum", shQuote(normalizePath(path, mustWork = TRUE)), stdout = TRUE)
  strsplit(value[[1L]], "[[:space:]]+")[[1L]][[1L]]
}
for (id in names(INPUTS)) if (!file.exists(INPUTS[[id]])) fail("Missing input: ", INPUTS[[id]])
input_manifest <- data.table(input_id = names(INPUTS), path = normalizePath(INPUTS, mustWork = TRUE))
input_manifest[, sha256 := vapply(path, sha256_file, character(1))]
input_manifest[, pinned_sha256 := PINNED[input_id]]
if (input_manifest[!is.na(pinned_sha256), any(sha256 != pinned_sha256)]) fail("Input identity mismatch")
# The candidate's own manifest pins its fit artifacts.
full_manifest <- fread(file.path(FULL_ROOT, "output_manifest.tsv"))
for (f in c("meta_matched.rds", "model_design.tsv", "deg_results.csv")) {
  recorded <- full_manifest[relative_path == file.path("primary_source_controls", f), sha256]
  if (length(recorded) != 1L || sha256_file(file.path(ARM_DIR, f)) != recorded) fail("Candidate artifact changed: ", f)
}

source(file.path(PROJECT_ROOT, "scripts/manuscript/resource_f_five/qualified_bulk_inputs.R"))
config <- yaml::read_yaml(INPUTS[["config"]])$datasets
mega <- names(Filter(function(d) isTRUE(d$de$include_in_mega), config))
if (length(mega) != 5L) fail("Frozen config does not define exactly five cohorts")
cohorts_in <- mega
gene_metadata <- fread(INPUTS[["gene_metadata"]])
gene_metadata[, eb := sub("[.][0-9]+$", "", gene_id)]

# -------------------------------------------------- target-positive persons from the screen
person_calls <- fread(PERSON_CALLS, colClasses = "character")
need <- c("basis", "cohort", "unit", "virus", "call", "n_libraries", "libraries")
if (!all(need %in% names(person_calls))) fail("person_calls.tsv lacks columns: ",
                                              paste(setdiff(need, names(person_calls)), collapse = ", "))
positive_rows <- person_calls[basis == "person_roster774" & virus %in% c("HBV", "HCV", "HEV") &
                                call == "positive"]
viral_ids <- sort(unique(unlist(strsplit(positive_rows$libraries, ",", fixed = TRUE))))
cat("screen-positive libraries in the 774 roster:", length(viral_ids), "\n")
print(positive_rows[, .(cohort, unit, virus, call, n_libraries, libraries)])
excluded_ids <- if (MODE == "exclude") viral_ids else character(0)
covariate_ids <- if (MODE == "covariate") viral_ids else character(0)

# ------------------------------------------------ copies of the refit model steps
integrate_counts <- function(counts_raw, meta, pass_ids) {
  merged <- counts_raw[, pass_ids]
  meta <- meta[match(pass_ids, sample_id)]
  stopifnot(identical(meta$sample_id, colnames(merged)))
  keep_samples <- meta$dataset %in% cohorts_in & !(meta$sample_id %in% excluded_ids)
  merged <- merged[, keep_samples, drop = FALSE]
  meta <- meta[keep_samples]
  stopifnot(identical(meta$sample_id, colnames(merged)), uniqueN(meta$dataset) == length(cohorts_in))
  dge <- DGEList(counts = merged)
  dge$samples$dataset <- meta$dataset
  dge$samples$condition <- factor(meta$condition)
  dge$samples$group_binary <- factor(meta$group_binary, levels = c("Control", "Disease"))
  dge$samples$sex <- meta$sex_final
  design_filter <- model.matrix(~ 0 + group_binary, data = dge$samples)
  keep <- filterByExpr(dge, design = design_filter)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  calcNormFactors(dge, method = "RLE")
}

fit_disease_model <- function(dge, meta) {
  dge_mega <- dge[, dge$samples$dataset %in% cohorts_in]
  assert_bulk_unit_alignment(meta, colnames(dge_mega))
  sex <- meta$inferred_sex[match(colnames(dge_mega), meta$sample_id)]
  if (anyNA(sex)) fail("Missing inferred_sex for model samples")
  info <- data.frame(
    group_binary = factor(dge_mega$samples$group_binary, levels = c("Control", "Disease")),
    dataset = factor(dge_mega$samples$dataset),
    inferred_sex = factor(sex),
    viral_positive = factor(ifelse(colnames(dge_mega) %in% covariate_ids, "yes", "no"),
                            levels = c("no", "yes"))
  )
  rownames(info) <- colnames(dge_mega)
  design <- if (MODE == "covariate") {
    model.matrix(~ dataset + inferred_sex + viral_positive + group_binary, data = info)
  } else {
    model.matrix(~ dataset + inferred_sex + group_binary, data = info)
  }
  coef_name <- "group_binaryDisease"
  design_rank <- qr(design)$rank
  if (!coef_name %in% colnames(design)) fail("Disease coefficient is not estimable")
  if (MODE == "covariate" && !"viral_positiveyes" %in% colnames(design)) fail("viral_positive term is not estimable")
  if (design_rank != ncol(design)) fail("Design matrix is rank deficient")

  v <- voomWithQualityWeights(dge_mega, design, save.plot = TRUE)
  sw <- v$targets$sample.weights
  if (is.null(sw) || length(sw) != ncol(dge_mega) || !identical(rownames(v$targets), colnames(dge_mega)) ||
      any(!is.finite(sw)) || any(sw <= 0)) {
    fail("voomWithQualityWeights did not return ordered positive finite sample weights")
  }
  fit0 <- lmFit(v, design)
  fit <- eBayes(fit0)
  res <- topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
  res$gene <- rownames(res)
  dt <- as.data.table(res)
  setnames(dt, "adj.P.Val", "padj")
  coef_index <- match(coef_name, colnames(fit$coefficients))
  moderated_se <- fit$stdev.unscaled[, coef_index] * sqrt(fit$s2.post)
  names(moderated_se) <- rownames(fit$coefficients)
  dt[, SE := moderated_se[match(gene, names(moderated_se))]]

  ttm <- topTreat(treat(fit0, lfc = TREAT_LFC), coef = coef_name, number = Inf, sort.by = "none")
  dt[, treat_lfc := TREAT_LFC]
  dt[, treat_p := ttm[gene, "P.Value"]]
  dt[, treat_fdr := ttm[gene, "adj.P.Val"]]

  ok <- is.finite(dt$logFC) & is.finite(dt$SE) & dt$SE > 0
  ash <- ashr::ash(dt$logFC[ok], dt$SE[ok], mixcompdist = "normal")
  dt[, shrunk_logFC := NA_real_][ok, shrunk_logFC := ash$result$PosteriorMean]
  dt[, lfsr := NA_real_][ok, lfsr := ash$result$lfsr]
  out <- dt[, .(gene, logFC, SE, t, P.Value, padj, shrunk_logFC, lfsr,
                treat_lfc, treat_p, treat_fdr, AveExpr)]
  out[, symbol := gene_metadata[match(sub("[.][0-9]+$", "", gene), eb), gene_name]]

  numeric_fields <- c("logFC", "SE", "t", "P.Value", "padj", "treat_lfc", "treat_p", "treat_fdr", "AveExpr")
  if (!identical(out$gene, rownames(dge_mega)) || anyDuplicated(out$gene) ||
      any(vapply(numeric_fields, function(f) any(!is.finite(out[[f]])), logical(1))) ||
      any(out$SE <= 0)) {
    fail("Coefficient/TREAT output failed the numerical contract")
  }
  viral_coef <- NULL
  if (MODE == "covariate") {
    vt <- topTable(fit, coef = "viral_positiveyes", number = Inf, sort.by = "none")
    vt$gene <- rownames(vt)
    viral_coef <- as.data.table(vt)
    setnames(viral_coef, "adj.P.Val", "padj")
    viral_coef[, symbol := gene_metadata[match(sub("[.][0-9]+$", "", gene), eb), gene_name]]
    setcolorder(viral_coef, c("gene", "symbol"))
  }
  list(out = out, design = design, design_rank = design_rank, info = info,
       dge = dge_mega, weights = sw, viral_coef = viral_coef)
}

deg_counts <- function(out) {
  canonical <- out[padj < 0.05 & abs(logFC) > 0.5]
  treat_calls <- out[treat_fdr < 0.05]
  list(
    canonical = list(n = nrow(canonical), up = sum(canonical$logFC > 0), down = sum(canonical$logFC < 0)),
    treat = list(n = nrow(treat_calls), up = sum(treat_calls$logFC > 0), down = sum(treat_calls$logFC < 0))
  )
}

write_table_once <- function(x, path, sep = "\t") {
  if (file.exists(path)) fail("Refusing overwrite: ", path)
  connection <- if (grepl("[.]gz$", path)) gzfile(path, open = "wb") else file(path, open = "wb")
  on.exit(try(close(connection), silent = TRUE), add = TRUE)
  write.table(as.data.frame(x), connection, sep = sep, row.names = FALSE,
              col.names = TRUE, quote = FALSE, na = "NA", eol = "\n")
}

# ---------------------------------------------------------------------- fit
dir.create(OUT_ROOT, recursive = FALSE, mode = "0750")
OUT_ROOT <- normalizePath(OUT_ROOT, mustWork = TRUE)
cat("Mode:", MODE, "| cohorts in model:", paste(cohorts_in, collapse = ", "), "\n")

counts_raw <- readRDS(INPUTS[["counts_raw"]])
meta <- as.data.table(readRDS(INPUTS[["meta_matched"]]))
design_ids <- fread(INPUTS[["model_design"]], select = "sample_id")$sample_id
if (length(design_ids) != 774L || !setequal(design_ids, meta$sample_id)) fail("Candidate sample set is not the 774 units")
if (!all(viral_ids %in% design_ids)) fail("A screen-positive library is not one of the 774 model samples: ",
                                         paste(setdiff(viral_ids, design_ids), collapse = ", "))
if (MODE != "none" && length(viral_ids) == 0L) fail("No screen-positive persons: nothing to exclude or adjust for")
pass_ids <- design_ids
cat("screen-positive persons in the model (group, dataset, condition):\n")
print(meta[sample_id %in% viral_ids, .(sample_id, dataset, group_binary, condition, inferred_sex)])
dge <- integrate_counts(counts_raw, meta, pass_ids)
rm(counts_raw); invisible(gc())
meta_in <- meta[sample_id %in% colnames(dge)]
if (ncol(dge) != 774L - length(excluded_ids)) fail("Unexpected sample count after exclusion")
group_table <- table(dge$samples$group_binary)
cat(sprintf("samples %d (control %d, disease %d), genes after filter %d\n",
            ncol(dge), group_table[["Control"]], group_table[["Disease"]], nrow(dge)))
fitted <- fit_disease_model(dge, meta_in)
out <- fitted$out
counts <- deg_counts(out)
cat(sprintf("canonical %d (up %d / down %d) | TREAT %d (up %d / down %d)\n",
            counts$canonical$n, counts$canonical$up, counts$canonical$down,
            counts$treat$n, counts$treat$up, counts$treat$down))

write_table_once(out, file.path(OUT_ROOT, "deg_results.csv"), sep = ",")
write_table_once(data.table(
  sample_id = colnames(fitted$dge), dataset = fitted$dge$samples$dataset,
  group_binary = as.character(fitted$dge$samples$group_binary),
  inferred_sex = as.character(fitted$info$inferred_sex),
  viral_positive_screen = colnames(fitted$dge) %in% viral_ids,
  viral_term_in_model = MODE == "covariate" & colnames(fitted$dge) %in% viral_ids,
  lib_size = fitted$dge$samples$lib.size, norm_factor = fitted$dge$samples$norm.factors,
  sample_weight = as.numeric(fitted$weights)), file.path(OUT_ROOT, "viral_samples.tsv"))
write_table_once(data.table(sample_id = viral_ids,
                            used_as = fcase(MODE == "exclude", "excluded",
                                            MODE == "covariate", "covariate_term",
                                            default = "not_used"),
                            positive_rows[match(viral_ids, positive_rows$libraries), .(cohort, virus)]),
                 file.path(OUT_ROOT, "viral_positive_persons.tsv"))
if (!is.null(fitted$viral_coef)) {
  write_table_once(fitted$viral_coef, file.path(OUT_ROOT, "viral_positive_coefficient.tsv.gz"))
  cat(sprintf("viral_positive term: %d genes with BH padj < 0.05 of %d tested (descriptive; n = %d persons)\n",
              sum(fitted$viral_coef$padj < 0.05), nrow(fitted$viral_coef), length(viral_ids)))
}

# ------------------------------------------------------- comparison with the full fit
full <- fread(INPUTS[["full_deg"]])
full_calls <- full[, .(gene, symbol, full_logFC = logFC, full_padj = padj,
                       full_call = padj < 0.05 & abs(logFC) > 0.5)]
refit <- out[, .(gene, refit_logFC = logFC, refit_padj = padj, refit_call = padj < 0.05 & abs(logFC) > 0.5)]
cmp <- merge(full_calls, refit, by = "gene", all = TRUE)
cmp[, in_full_universe := !is.na(full_call)]
cmp[, in_refit_universe := !is.na(refit_call)]
cmp[is.na(full_call), full_call := FALSE]
cmp[is.na(refit_call), refit_call := FALSE]
cmp[, change := fcase(
  full_call & refit_call & sign(full_logFC) == sign(refit_logFC), "retained",
  full_call & refit_call, "retained_sign_flip",
  full_call & !in_refit_universe, "lost_untested",
  full_call, "lost",
  refit_call & !in_full_universe, "gained_new_gene",
  refit_call, "gained",
  default = "not_called")]
cmp[, full_call_opposite_sign_in_refit := full_call & in_refit_universe & sign(full_logFC) != sign(refit_logFC)]
setcolorder(cmp, c("gene", "symbol"))
write_table_once(cmp, file.path(OUT_ROOT, "refit_vs_full_gene_table.tsv.gz"))
changed <- cmp[!(change %in% c("retained", "not_called"))][order(change, -abs(fcoalesce(full_logFC, refit_logFC)))]
write_table_once(changed, file.path(OUT_ROOT, "changed_calls.tsv"))
cat("changed calls by class:\n"); print(changed[, .N, by = change])
common <- cmp[in_full_universe & in_refit_universe]
metrics <- data.table(
  metric = c("n_common_genes", "logFC_pearson_common", "logFC_spearman_common",
             "median_abs_delta_logFC", "sign_agreement_common",
             "logFC_pearson_full_calls", "n_full_calls", "n_retained", "n_lost", "n_lost_untested",
             "n_retained_sign_flip", "n_full_call_opposite_sign_in_refit", "n_gained", "n_gained_new_gene"),
  value = c(nrow(common), cor(common$full_logFC, common$refit_logFC),
            cor(common$full_logFC, common$refit_logFC, method = "spearman"),
            median(abs(common$refit_logFC - common$full_logFC)),
            mean(sign(common$full_logFC) == sign(common$refit_logFC)),
            cor(common[full_call == TRUE, full_logFC], common[full_call == TRUE, refit_logFC]),
            sum(cmp$full_call), sum(cmp$change == "retained"), sum(cmp$change == "lost"),
            sum(cmp$change == "lost_untested"), sum(cmp$change == "retained_sign_flip"),
            sum(cmp$full_call_opposite_sign_in_refit), sum(cmp$change == "gained"),
            sum(cmp$change == "gained_new_gene")))
write_table_once(metrics, file.path(OUT_ROOT, "refit_vs_full_metrics.tsv"))
print(metrics)

if (MODE == "none") {
  aligned <- identical(out$gene, full$gene)
  delta <- function(a, b) max(abs(a - b))
  check <- data.table(
    check = c("gene_order_identical", "n_genes", "logFC_max_abs_delta", "SE_max_abs_delta",
              "padj_max_abs_delta", "canonical_membership_identical", "canonical_n"),
    observed = c(aligned, nrow(out),
                 if (aligned) delta(out$logFC, full$logFC) else NA,
                 if (aligned) delta(out$SE, full$SE) else NA,
                 if (aligned) delta(out$padj, full$padj) else NA,
                 if (aligned) identical(out$padj < 0.05 & abs(out$logFC) > 0.5,
                                        full$padj < 0.05 & abs(full$logFC) > 0.5) else FALSE,
                 counts$canonical$n),
    expected = c(1, nrow(full), 1e-8, 1e-8, 1e-8, 1, 1854),
    rule = c("equal", "equal", "max", "max", "max", "equal", "equal"))
  check[, pass := !is.na(observed) & fifelse(rule == "equal", observed == expected, observed <= expected)]
  write_table_once(check, file.path(OUT_ROOT, "reproduction_check.tsv"))
  print(check)
  if (!all(check$pass)) fail("VIRAL_MODE=none did not reproduce the candidate fit")
}

summary <- c(
  list(mode = MODE, viral_positive_samples = viral_ids, n_excluded = length(excluded_ids),
       n_samples = ncol(fitted$dge), n_control = sum(fitted$info$group_binary == "Control"),
       n_disease = sum(fitted$info$group_binary == "Disease"), n_genes = nrow(out),
       design_columns = ncol(fitted$design), design_rank = fitted$design_rank,
       design_terms = colnames(fitted$design),
       viral_term_genes_padj_lt_0.05 = if (is.null(fitted$viral_coef)) NA_integer_ else sum(fitted$viral_coef$padj < 0.05),
       seed = 20260923L, treat_lfc = TREAT_LFC,
       per_cohort = lapply(split(fitted$info, fitted$info$dataset),
                           function(d) list(control = sum(d$group_binary == "Control"),
                                            disease = sum(d$group_binary == "Disease")))),
  counts)
writeLines(toJSON(summary, auto_unbox = TRUE, pretty = TRUE, digits = NA),
           file.path(OUT_ROOT, "viral_refit_summary.json"))
write_table_once(input_manifest, file.path(OUT_ROOT, "input_manifest.tsv"))
writeLines(capture.output(sessionInfo()), file.path(OUT_ROOT, "sessionInfo.txt"))
files <- setdiff(list.files(OUT_ROOT, recursive = TRUE), "output_manifest.tsv")
write_table_once(data.table(relative_path = files,
                            size_bytes = file.size(file.path(OUT_ROOT, files)),
                            sha256 = vapply(file.path(OUT_ROOT, files), sha256_file, character(1))),
                 file.path(OUT_ROOT, "output_manifest.tsv"))
cat("VIRAL_REFIT_COMPLETE", MODE, OUT_ROOT, "\n")
