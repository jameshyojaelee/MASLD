#!/usr/bin/env Rscript
# Refit the selected five-cohort F_five disease-versus-control model with the
# corrected GSE130970 control labels (review item 1, 2026-09-23).
#
# GSE130970 controls were set by fibrosis_stage == 0, which put 17 source NAFLD
# biopsies (16 of them steatotic) in the control arm. The source labels now come from
# RNA-seq/Human/Patient_Cohorts/metadata/source/GSE130970/GSE130970_source_diagnosis.tsv
# (Hoang et al. 2019 Supplementary Table 1).
#
# REFIT_MODE=fit (default), into REFIT_OUTPUT_ROOT:
#   reproduction_guard/       old labels, exact 03_integrate_counts.R filter + RLE
#                             and 05h_limma_voom_qw_canonical.R model. Must
#                             reproduce the stored F_five fit or the job stops.
#   primary_source_controls/  6 source controls (140 control / 704 disease).
#   sensitivity_strict_nas0/  4 controls with steatosis 0 and NAS 0 (138 / 706).
#   sensitivity_drop_sex_discordant/
#                             primary labels without the 3 GSE130970 biopsies whose
#                             GEO sex disagrees with inferred sex (841; 139 / 702).
# REFIT_MODE=reproduce, into REFIT_OUTPUT_ROOT, reading REFIT_SOURCE_ROOT:
#   refits each arm from its saved DGE and metadata, requires agreement with the
#   fit-mode deg_results.csv, and writes the model_input_manifest.tsv /
#   VALIDATED.json contract that the stage producers read as STAGE_RELEASE_ROOT.
#
# Design, filter and tests are unchanged: filterByExpr(~ 0 + group_binary), RLE,
# voomWithQualityWeights, ~ dataset + inferred_sex + group_binary, eBayes,
# TREAT lfc 0.25, ashr normal mixture. Unit of inference: one biopsy per
# participant (844 participants, 5 cohorts).

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
  "MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
), mustWork = TRUE)
MODE <- Sys.getenv("REFIT_MODE", "fit")
if (!MODE %in% c("fit", "reproduce")) fail("REFIT_MODE must be fit or reproduce")
OUT_ROOT <- Sys.getenv("REFIT_OUTPUT_ROOT", "")
if (!nzchar(OUT_ROOT)) fail("REFIT_OUTPUT_ROOT is required")
root_pattern <- if (MODE == "fit") {
  "^bulk-corrected-controls-[0-9]{8}T[0-9]{6}Z$"
} else {
  "^bulk-corrected-controls-[0-9]{8}T[0-9]{6}Z-reproduction$"
}
if (!grepl(root_pattern, basename(OUT_ROOT))) fail("Unexpected output root name: ", basename(OUT_ROOT))
if (file.exists(OUT_ROOT)) fail("REFIT_OUTPUT_ROOT must not exist: ", OUT_ROOT)
if (!dir.exists(dirname(OUT_ROOT))) fail("Output parent is absent: ", dirname(OUT_ROOT))
if (!dir.create(OUT_ROOT, recursive = FALSE, mode = "0750")) fail("Cannot create ", OUT_ROOT)
OUT_ROOT <- normalizePath(OUT_ROOT, mustWork = TRUE)

CANDIDATE_ROOT <- file.path(
  PROJECT_ROOT,
  "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10"
)
F5 <- file.path(CANDIDATE_ROOT, "inputs/BG001-DECISION/arms/F_five")
MODEL_INPUTS <- file.path(CANDIDATE_ROOT, "inputs/BULK-F-FIVE/frozen_model_inputs")
INPUTS <- c(
  counts_raw = file.path(F5, "results/integration/merged_counts_raw.rds"),
  meta_matched = file.path(F5, "results/integration/meta_matched.rds"),
  qc_report = file.path(F5, "qc/sample_qc_report.csv"),
  reference_dge = file.path(F5, "results/integration/merged_dge.rds"),
  reference_deg = file.path(F5, "results/integration/deg_results.csv"),
  config = file.path(MODEL_INPUTS, "human_datasets.yaml"),
  gene_metadata = file.path(MODEL_INPUTS, "gencode_v49_gene_metadata.tsv.gz"),
  source_diagnosis = file.path(
    PROJECT_ROOT,
    "RNA-seq/Human/Patient_Cohorts/metadata/source/GSE130970/GSE130970_source_diagnosis.tsv"
  )
)
# SHA-256 of the adopted F_five artifacts (arms/F_five/provenance/artifact_manifest.tsv).
PINNED <- c(
  counts_raw = "16afc00bf1db706731d07df1dfeabdad0bc423b1caa1b822631685a55ec7f225",
  meta_matched = "90ba6aca68643c4680c08fc982a6c72b76edb740ab726874cb50603711756a83",
  qc_report = "1c7fef991c599554369ecc79de2b2f441597cf9b6c1a7ff3ecbe05f1c5969560",
  reference_dge = "bc3ea19c8b064861339eba862abf58e82ddb5460c622b550dc45aa3b5fcd1f76",
  reference_deg = "691b09517f42ff24e472fdb9db460b1c339ec734876c40115f8b146293b5b633",
  config = "dc64f6f051022404689abcfeb4c71001f135123c7af65c5bc8d88dd663541f6a",
  gene_metadata = "49dcaae77323ffb20f5fd8399bb2899cd553e3634755a37b60a133d62549c189"
)
ARMS <- c("primary_source_controls", "sensitivity_strict_nas0", "sensitivity_drop_sex_discordant")
# GSE130970 biopsies whose GEO sex disagrees with the XIST/Y k-means sex
# (01_sample_qc.R:231-238). The model uses inferred_sex; the source-control
# reconstruction used GEO sex, and SRR9036376 is one of the 6 source controls.
SEX_DISCORDANT_GSE130970 <- c("SRR9036338", "SRR9036344", "SRR9036376")
# Every reported-versus-inferred sex mismatch in the 844 (checked 2026-09-23);
# the three GSE162694 samples stay in every arm.
SEX_DISCORDANT_ALL <- c(SEX_DISCORDANT_GSE130970, "SRR13199861", "SRR13199865", "SRR13199917")
ARM_EXCLUDED <- list(
  primary_source_controls = character(),
  sensitivity_strict_nas0 = character(),
  sensitivity_drop_sex_discordant = SEX_DISCORDANT_GSE130970
)
ARM_RULE <- c(
  primary_source_controls = "GSE130970 source_control_status (Hoang 2019)",
  sensitivity_strict_nas0 = "GSE130970 strict_control_nas0 (steatosis 0 and NAS 0)",
  sensitivity_drop_sex_discordant = paste0(
    "GSE130970 source_control_status; reported/inferred sex-discordant GSE130970 biopsies removed: ",
    paste(SEX_DISCORDANT_GSE130970, collapse = ", "))
)
TREAT_LFC <- 0.25
LOGFC_TOLERANCE <- 1e-8

sha256_file <- function(path) {
  value <- system2("sha256sum", shQuote(normalizePath(path, mustWork = TRUE)), stdout = TRUE)
  strsplit(value[[1L]], "[[:space:]]+")[[1L]][[1L]]
}
write_tsv_once <- function(x, path, sep = "\t") {
  if (file.exists(path)) fail("Refusing overwrite: ", path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  fwrite(x, path, sep = sep, quote = FALSE, na = "NA")
}
# deg_results.csv keeps the byte conventions of the reference producer
# (write.table, 15 significant digits, unquoted).
write_table_once <- function(x, path, sep = "\t") {
  if (file.exists(path)) fail("Refusing overwrite: ", path)
  connection <- if (grepl("[.]gz$", path)) gzfile(path, open = "wb") else file(path, open = "wb")
  on.exit(try(close(connection), silent = TRUE), add = TRUE)
  write.table(as.data.frame(x), connection, sep = sep, row.names = FALSE,
              col.names = TRUE, quote = FALSE, na = "NA", eol = "\n")
}
write_json_once <- function(x, path) {
  if (file.exists(path)) fail("Refusing overwrite: ", path)
  writeLines(toJSON(x, auto_unbox = TRUE, pretty = TRUE, digits = NA), path)
}
save_rds_once <- function(x, path) {
  if (file.exists(path)) fail("Refusing overwrite: ", path)
  saveRDS(x, path)
}

for (id in names(INPUTS)) if (!file.exists(INPUTS[[id]])) fail("Missing input: ", INPUTS[[id]])
input_manifest <- data.table(input_id = names(INPUTS), path = normalizePath(INPUTS, mustWork = TRUE))
input_manifest[, sha256 := vapply(path, sha256_file, character(1))]
input_manifest[, pinned_sha256 := PINNED[input_id]]
if (input_manifest[!is.na(pinned_sha256), any(sha256 != pinned_sha256)]) {
  print(input_manifest[!is.na(pinned_sha256) & sha256 != pinned_sha256])
  fail("Adopted F_five input identity mismatch")
}

config <- yaml::read_yaml(INPUTS[["config"]])$datasets
mega <- names(Filter(function(d) isTRUE(d$de$include_in_mega), config))
if (length(mega) != 5L) fail("Frozen config does not define exactly five cohorts")
gene_metadata <- fread(INPUTS[["gene_metadata"]])
gene_metadata[, eb := sub("[.][0-9]+$", "", gene_id)]

# ---------------------------------------------------------------- model steps
# 03_integrate_counts.R, remediation mode, canonical_five scope, native genes.
integrate_counts <- function(counts_raw, meta, pass_ids) {
  merged <- counts_raw[, pass_ids]
  meta <- meta[match(pass_ids, sample_id)]
  stopifnot(identical(meta$sample_id, colnames(merged)))
  keep_samples <- meta$dataset %in% mega
  merged <- merged[, keep_samples, drop = FALSE]
  meta <- meta[keep_samples]
  stopifnot(identical(meta$sample_id, colnames(merged)), uniqueN(meta$dataset) == 5L)
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

# 05h_limma_voom_qw_canonical.R, remediation mode.
fit_disease_model <- function(dge, meta) {
  dge_mega <- dge[, dge$samples$dataset %in% mega]
  sex <- meta$inferred_sex[match(colnames(dge_mega), meta$sample_id)]
  if (anyNA(sex)) fail("Missing inferred_sex for model samples")
  info <- data.frame(
    group_binary = factor(dge_mega$samples$group_binary, levels = c("Control", "Disease")),
    dataset = factor(dge_mega$samples$dataset),
    inferred_sex = factor(sex)
  )
  rownames(info) <- colnames(dge_mega)
  design <- model.matrix(~ dataset + inferred_sex + group_binary, data = info)
  coef_name <- "group_binaryDisease"
  design_rank <- qr(design)$rank
  if (!coef_name %in% colnames(design)) fail("Disease coefficient is not estimable")
  if (design_rank != ncol(design)) fail("C2 design matrix is rank deficient")

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
  list(out = out, design = design, design_rank = design_rank, info = info,
       dge = dge_mega, weights = sw, voom = v)
}

deg_counts <- function(out) {
  canonical <- out[padj < 0.05 & abs(logFC) > 0.5]
  treat_calls <- out[treat_fdr < 0.05]
  list(
    canonical = list(n = nrow(canonical), up = sum(canonical$logFC > 0), down = sum(canonical$logFC < 0)),
    treat = list(n = nrow(treat_calls), up = sum(treat_calls$logFC > 0), down = sum(treat_calls$logFC < 0))
  )
}

# Control status for control-dependent downstream rules: the source label where
# the cohort has one (GSE130970), group_binary elsewhere.
source_control <- function(meta) {
  if ("source_control_status" %in% names(meta)) {
    fifelse(!is.na(meta$source_control_status), meta$source_control_status == "Control",
            meta$group_binary == "Control")
  } else {
    meta$group_binary == "Control"
  }
}

# Figure 3D endpoint roster (fig3d_endpoint_validation.R) and the F0 reference
# composition of the F0-arm sensitivity (build_pi_stage_f0_arms.R), from labels only.
downstream_label_counts <- function(meta) {
  m <- copy(meta)
  m[, control := source_control(m)]
  m[, advanced := (!is.na(nas_score) & nas_score >= 5) | (!is.na(fibrosis_stage) & fibrosis_stage >= 3)]
  m[, strict_control := control & (is.na(nas_score) | nas_score <= 1) &
                        (is.na(fibrosis_stage) | fibrosis_stage <= 0)]
  m[, endpoint_code := fifelse(advanced, "Advanced", fifelse(strict_control, "StrictControl", NA_character_))]
  per_cohort <- dcast(m[!is.na(endpoint_code), .N, by = .(dataset, endpoint_code)],
                      dataset ~ endpoint_code, value.var = "N", fill = 0L)
  eligible <- sort(per_cohort[Advanced >= 5L & StrictControl >= 5L, dataset])
  endpoint <- m[dataset %in% eligible & !is.na(endpoint_code)]
  f0 <- m[dataset %in% c("GSE130970", "GSE135251", "GSE162694") & fibrosis_stage == 0L]
  list(
    fig3d_endpoint = list(
      n_total = nrow(endpoint), n_cohorts = length(eligible),
      n_advanced = sum(endpoint$endpoint_code == "Advanced"),
      n_strict_control = sum(endpoint$endpoint_code == "StrictControl"),
      eligible_cohorts = eligible,
      strict_control_by_cohort = as.list(setNames(per_cohort$StrictControl, per_cohort$dataset))
    ),
    f0_reference = list(
      A_all_F0 = list(n_ref_control = sum(f0$control), n_ref_disease = sum(!f0$control)),
      B_disease_only_F0 = list(n_ref_control = 0L, n_ref_disease = sum(!f0$control)),
      C_all_F0_plus_term = list(n_ref_control = sum(f0$control), n_ref_disease = sum(!f0$control))
    )
  )
}

# accepted_reference_path names the corrected fit a consumer should read: the
# arm's own deg_results.csv in fit mode, the fit-mode table it reproduced in
# reproduce mode. The pre-correction F_five table is recorded separately.
write_fit_outputs <- function(fitted, arm_dir, dge_path, meta_path,
                              reference_path = file.path(arm_dir, "deg_results.csv")) {
  out <- fitted$out
  write_table_once(out, file.path(arm_dir, "deg_results.csv"), sep = ",")
  write_table_once(out[, .(gene, symbol, logFC, SE, t, P.Value, padj, AveExpr)],
                   file.path(arm_dir, "coefficient_table.tsv.gz"))
  write_table_once(out[, .(gene, symbol, logFC, SE, AveExpr, treat_lfc, treat_p, treat_fdr)],
                   file.path(arm_dir, "treat_table.tsv.gz"))
  design_table <- as.data.table(fitted$design)
  design_table[, sample_id := rownames(fitted$info)]
  setcolorder(design_table, c("sample_id", setdiff(names(design_table), "sample_id")))
  write_table_once(design_table, file.path(arm_dir, "model_design.tsv"))
  write_table_once(data.table(sample_id = colnames(fitted$dge), sample_weight = as.numeric(fitted$weights)),
                   file.path(arm_dir, "lvqw_sample_weights.tsv"))
  save_rds_once(list(voom_xy = fitted$voom$voom.xy, voom_line = fitted$voom$voom.line),
                file.path(arm_dir, "voom_curve.rds"))
  write_table_once(data.table(
    dge_path = normalizePath(dge_path, mustWork = TRUE),
    meta_path = normalizePath(meta_path, mustWork = TRUE),
    config_path = normalizePath(INPUTS[["config"]], mustWork = TRUE),
    gene_metadata_path = normalizePath(INPUTS[["gene_metadata"]], mustWork = TRUE),
    accepted_reference_path = normalizePath(reference_path, mustWork = TRUE),
    pre_correction_fit_path = normalizePath(INPUTS[["reference_deg"]], mustWork = TRUE),
    n_genes = nrow(fitted$dge), n_samples = ncol(fitted$dge),
    n_cohorts = nlevels(fitted$info$dataset),
    n_control = sum(fitted$info$group_binary == "Control"),
    n_disease = sum(fitted$info$group_binary == "Disease"),
    treat_lfc = TREAT_LFC
  ), file.path(arm_dir, "model_input_manifest.tsv"))
  data.table(arm = basename(arm_dir), n_columns = ncol(fitted$design), rank = fitted$design_rank,
             full_rank = fitted$design_rank == ncol(fitted$design), coefficient = "group_binaryDisease",
             covariates = "dataset+inferred_sex")
}

max_abs_delta <- function(a, b) max(abs(a - b))

# ======================================================================= fit
if (MODE == "fit") {
  counts_raw <- readRDS(INPUTS[["counts_raw"]])
  meta_old <- as.data.table(readRDS(INPUTS[["meta_matched"]]))
  for (col in names(meta_old)) if (is.character(meta_old[[col]])) meta_old[get(col) == "", (col) := NA]
  qc <- fread(INPUTS[["qc_report"]])
  pass_ids <- qc[pass_technical == TRUE, sample_id]
  reference_dge <- readRDS(INPUTS[["reference_dge"]])
  reference_deg <- fread(INPUTS[["reference_deg"]])
  model_ids <- colnames(reference_dge)
  if (length(model_ids) != 844L) fail("Reference F_five DGE is not 844 participants")

  # ---- corrected GSE130970 labels, joined by run accession
  src <- fread(INPUTS[["source_diagnosis"]], colClasses = "character")
  if (nrow(src) != 78L || anyDuplicated(src$run) ||
      !setequal(src$source_control_status, c("Control", "NAFLD")) ||
      !setequal(src$strict_control_nas0, c("Control", "NAFLD"))) {
    fail("Source diagnosis table is not 78 unique runs with Control/NAFLD labels")
  }
  if (src[source_control_status == "Control" & steatosis_grade != "0", .N] ||
      src[strict_control_nas0 == "Control" & source_control_status != "Control", .N]) {
    fail("A steatotic biopsy or a non-source control is in a control arm")
  }
  gse <- meta_old[dataset == "GSE130970"]
  retained <- intersect(gse$sample_id, model_ids)
  if (length(retained) != 76L) fail("Expected 76 retained GSE130970 biopsies, found ", length(retained))
  unmatched <- setdiff(gse$sample_id, src$run)
  if (length(unmatched)) fail("GSE130970 samples without a source label: ", paste(unmatched, collapse = ", "))
  if (meta_old[sample_id %in% retained & group_binary == "Control", .N] != 23L) {
    fail("Old labels do not give 23 retained GSE130970 controls")
  }
  src_join <- src[match(gse$sample_id, run)]
  if (any(as.integer(src_join$fibrosis_stage) != gse$fibrosis_stage) ||
      any(as.integer(src_join$nafld_activity_score) != gse$nas_score)) {
    fail("Source table fibrosis/NAS disagree with the F_five metadata")
  }

  relabel <- function(meta, control_column) {
    m <- copy(meta)
    s <- src[match(m$sample_id, run)]
    is_gse <- m$dataset == "GSE130970"
    m[, `:=`(
      source_control_status = fifelse(is_gse, s$source_control_status, NA_character_),
      strict_control_nas0 = fifelse(is_gse, s$strict_control_nas0, NA_character_),
      steatosis_grade = fifelse(is_gse, as.integer(s$steatosis_grade), NA_integer_),
      lobular_grade = fifelse(is_gse, as.integer(s$lobular_inflammation_grade), NA_integer_),
      ballooning_grade = fifelse(is_gse, as.integer(s$cytological_ballooning_grade), NA_integer_)
    )]
    arm_control <- s[[control_column]] == "Control"
    # Mirrors the patched 00_harmonize_metadata.R GSE130970 branch.
    m[is_gse, group_binary := fifelse(arm_control[is_gse], "Control", "Disease")]
    m[is_gse, condition := fcase(
      group_binary == "Control", "Control",
      fibrosis_stage == 0L, "Fibrosis_F0",
      fibrosis_stage == 1L, "Fibrosis_F1",
      fibrosis_stage %in% 2:3, "Fibrosis_F2F3",
      fibrosis_stage == 4L, "Fibrosis_F4"
    )]
    m[is_gse & source_control_status == "Control", diagnosis_harmonized := "Control"]
    m
  }
  metas <- list(
    reproduction_guard = meta_old,
    primary_source_controls = relabel(meta_old, "source_control_status"),
    sensitivity_strict_nas0 = relabel(meta_old, "strict_control_nas0")
  )
  metas$sensitivity_drop_sex_discordant <- metas$primary_source_controls
  stopifnot(metas$primary_source_controls[sample_id == "SRR9036353", group_binary] == "Disease",
            meta_old[sample_id == "SRR9036353", group_binary] == "Control")

  # ---- expected control/disease totals, from the source table joined to the 844
  in_model <- function(m) m[sample_id %in% model_ids]
  old_counts <- in_model(meta_old)[, .(old_control = sum(group_binary == "Control"),
                                       old_disease = sum(group_binary == "Disease")), by = dataset]
  retained_src <- src[run %in% model_ids]
  old_non_gse <- in_model(meta_old)[dataset != "GSE130970", .(control = sum(group_binary == "Control"),
                                                             disease = sum(group_binary == "Disease"))]
  totals_from_source <- function(control_column) {
    n_gse_control <- retained_src[get(control_column) == "Control", .N]
    c(control = old_non_gse$control + n_gse_control,
      disease = old_non_gse$disease + nrow(retained_src) - n_gse_control)
  }
  dropped <- src[run %in% SEX_DISCORDANT_GSE130970]
  n_dropped_control <- dropped[source_control_status == "Control", .N]
  expected_totals <- list(
    primary_source_controls = totals_from_source("source_control_status"),
    sensitivity_strict_nas0 = totals_from_source("strict_control_nas0"),
    sensitivity_drop_sex_discordant = totals_from_source("source_control_status") -
      c(n_dropped_control, nrow(dropped) - n_dropped_control)
  )
  stopifnot(identical(unname(expected_totals$primary_source_controls), c(140L, 704L)),
            identical(unname(expected_totals$sensitivity_strict_nas0), c(138L, 706L)),
            identical(unname(expected_totals$sensitivity_drop_sex_discordant), c(139L, 702L)))

  # ---- reported-versus-inferred sex discordance among the 844
  sex_discordance <- in_model(meta_old)[!is.na(sex) & sex != inferred_sex,
    .(sample_id, dataset, reported_sex = sex, inferred_sex, sex_final, old_group_binary = group_binary,
      fibrosis_stage, nas_score)]
  sex_discordance[, primary_group_binary := metas$primary_source_controls[match(sex_discordance$sample_id, sample_id), group_binary]]
  sex_discordance[, excluded_in_sex_sensitivity_arm := sample_id %in% SEX_DISCORDANT_GSE130970]
  write_tsv_once(sex_discordance, file.path(OUT_ROOT, "sex_discordance_844.tsv"))
  print(sex_discordance)
  if (!setequal(sex_discordance[dataset == "GSE130970", sample_id], SEX_DISCORDANT_GSE130970)) {
    fail("GSE130970 sex-discordant set differs from ", paste(SEX_DISCORDANT_GSE130970, collapse = ", "))
  }
  if (!setequal(sex_discordance$sample_id, SEX_DISCORDANT_ALL)) {
    fail("Sex-discordant set in the 844 differs from the six recorded samples")
  }

  label_table <- merge(
    meta_old[dataset == "GSE130970", .(sample_id, old_group_binary = group_binary, old_condition = condition,
                                       fibrosis_stage, nas_score, retained_in_model = sample_id %in% model_ids)],
    metas$primary_source_controls[dataset == "GSE130970", .(sample_id, source_control_status, strict_control_nas0,
      steatosis_grade, lobular_grade, ballooning_grade,
      primary_group_binary = group_binary, primary_condition = condition)],
    by = "sample_id"
  )
  label_table[, strict_group_binary := metas$sensitivity_strict_nas0[match(label_table$sample_id, sample_id), group_binary]]
  label_table[, relabelled_primary := old_group_binary != primary_group_binary]
  write_tsv_once(label_table, file.path(OUT_ROOT, "gse130970_label_changes.tsv"))

  # ---- 1. reproduction guard (old labels)
  cat("\n=== reproduction guard: old labels ===\n")
  guard_dir <- file.path(OUT_ROOT, "reproduction_guard")
  guard_dge <- integrate_counts(counts_raw, meta_old, pass_ids)
  guard_fit <- fit_disease_model(guard_dge, meta_old)
  guard_out <- guard_fit$out
  guard_counts <- deg_counts(guard_out)
  ref_counts <- deg_counts(reference_deg)
  aligned <- identical(guard_out$gene, reference_deg$gene)
  guard_labels <- downstream_label_counts(in_model(meta_old))
  checks <- data.table(
    check = c("dge_gene_ids_identical", "dge_sample_ids_identical", "dge_counts_identical",
              "lib_size_max_abs_delta", "norm_factors_max_abs_delta", "n_genes",
              "canonical_n", "canonical_up", "canonical_down", "treat_n", "treat_up", "treat_down",
              "deg_gene_order_identical", "logFC_max_abs_delta", "SE_max_abs_delta",
              "padj_max_abs_delta", "treat_fdr_max_abs_delta", "lfsr_max_abs_delta",
              "canonical_membership_identical", "treat_membership_identical",
              "fig3d_endpoint_total", "fig3d_strict_controls", "fig3d_advanced", "fig3d_cohorts",
              "f0_ref_control", "f0_ref_disease"),
    observed = c(
      identical(rownames(guard_dge), rownames(reference_dge)),
      identical(colnames(guard_dge), colnames(reference_dge)),
      identical(unname(guard_dge$counts), unname(reference_dge$counts)),
      max_abs_delta(guard_dge$samples$lib.size, reference_dge$samples$lib.size),
      max_abs_delta(guard_dge$samples$norm.factors, reference_dge$samples$norm.factors),
      nrow(guard_out),
      guard_counts$canonical$n, guard_counts$canonical$up, guard_counts$canonical$down,
      guard_counts$treat$n, guard_counts$treat$up, guard_counts$treat$down,
      aligned,
      if (aligned) max_abs_delta(guard_out$logFC, reference_deg$logFC) else NA,
      if (aligned) max_abs_delta(guard_out$SE, reference_deg$SE) else NA,
      if (aligned) max_abs_delta(guard_out$padj, reference_deg$padj) else NA,
      if (aligned) max_abs_delta(guard_out$treat_fdr, reference_deg$treat_fdr) else NA,
      if (aligned) max_abs_delta(guard_out$lfsr, reference_deg$lfsr) else NA,
      if (aligned) identical(guard_out$padj < 0.05 & abs(guard_out$logFC) > 0.5,
                             reference_deg$padj < 0.05 & abs(reference_deg$logFC) > 0.5) else FALSE,
      if (aligned) identical(guard_out$treat_fdr < 0.05, reference_deg$treat_fdr < 0.05) else FALSE,
      guard_labels$fig3d_endpoint$n_total, guard_labels$fig3d_endpoint$n_strict_control,
      guard_labels$fig3d_endpoint$n_advanced, guard_labels$fig3d_endpoint$n_cohorts,
      guard_labels$f0_reference$A_all_F0$n_ref_control, guard_labels$f0_reference$A_all_F0$n_ref_disease
    ),
    expected = c(1, 1, 1, 0, 1e-12, 23370, 1347, 1003, 344, 1616, 1144, 472, 1,
                 LOGFC_TOLERANCE, 1e-8, 1e-8, 1e-8, 1e-6, 1, 1, 418, 115, 303, 4, 31, 73),
    rule = c(rep("equal", 3), "max", "max", rep("equal", 8), rep("max", 5), rep("equal", 8))
  )
  checks[, pass := !is.na(observed) & fifelse(rule == "equal", observed == expected, observed <= expected)]
  write_tsv_once(checks, file.path(guard_dir, "guard_checks.tsv"))
  print(checks)
  if (!all(checks$pass)) fail("Reproduction guard FAILED; no corrected arm was fitted")
  cat("reproduction guard PASS\n")
  rm(guard_fit, guard_dge); invisible(gc())

  # ---- 2. corrected arms
  reference_calls <- reference_deg[, .(gene, old_symbol = symbol, old_logFC = logFC, old_SE = SE,
                                       old_padj = padj, old_treat_fdr = treat_fdr,
                                       old_canonical = padj < 0.05 & abs(logFC) > 0.5,
                                       old_treat = treat_fdr < 0.05)]
  design_rows <- list()
  cohort_rows <- list()
  arm_count_rows <- list()
  for (arm in ARMS) {
    cat("\n=== arm ", arm, " ===\n", sep = "")
    arm_dir <- file.path(OUT_ROOT, arm)
    dir.create(arm_dir, recursive = FALSE, showWarnings = FALSE)
    meta_arm <- metas[[arm]]
    arm_ids <- setdiff(model_ids, ARM_EXCLUDED[[arm]])
    dge <- integrate_counts(counts_raw, meta_arm, setdiff(pass_ids, ARM_EXCLUDED[[arm]]))
    if (!setequal(colnames(dge), arm_ids)) fail("Arm ", arm, " sample set is not the 844 minus its exclusions")
    group_table <- table(dge$samples$group_binary)
    if (!identical(as.integer(group_table[c("Control", "Disease")]), unname(expected_totals[[arm]]))) {
      fail("Arm ", arm, " control/disease totals disagree with the source table")
    }
    dge_path <- file.path(arm_dir, "merged_dge.rds")
    meta_path <- file.path(arm_dir, "meta_matched.rds")
    save_rds_once(dge, dge_path)
    save_rds_once(meta_arm, meta_path)
    fitted <- fit_disease_model(dge, meta_arm)
    design_rows[[arm]] <- write_fit_outputs(fitted, arm_dir, dge_path, meta_path)
    out <- fitted$out
    counts <- deg_counts(out)
    cat(sprintf("  genes %d | canonical %d (up %d / down %d) | TREAT %d (up %d / down %d)\n",
                nrow(out), counts$canonical$n, counts$canonical$up, counts$canonical$down,
                counts$treat$n, counts$treat$up, counts$treat$down))

    universe <- merge(data.table(gene = rownames(reference_dge), in_reference_universe = TRUE),
                      data.table(gene = rownames(dge), in_arm_universe = TRUE), by = "gene", all = TRUE)
    universe[is.na(in_reference_universe), in_reference_universe := FALSE]
    universe[is.na(in_arm_universe), in_arm_universe := FALSE]
    universe[, status := fcase(in_reference_universe & in_arm_universe, "retained",
                               in_arm_universe, "gained", default = "lost")]
    universe[, symbol := gene_metadata[match(sub("[.][0-9]+$", "", gene), eb), gene_name]]
    write_tsv_once(universe, file.path(arm_dir, "gene_universe_diff.tsv.gz"))

    cohort <- merge(old_counts, meta_arm[sample_id %in% arm_ids][, .(new_control = sum(group_binary == "Control"),
                                                        new_disease = sum(group_binary == "Disease")), by = dataset],
                    by = "dataset")
    cohort[, `:=`(arm = arm, both_groups_present = new_control > 0L & new_disease > 0L)]
    write_tsv_once(cohort, file.path(arm_dir, "cohort_group_counts.tsv"))
    cohort_rows[[arm]] <- cohort

    comparison <- merge(reference_calls,
                        out[, .(gene, new_symbol = symbol, new_logFC = logFC, new_SE = SE, new_padj = padj,
                                new_treat_fdr = treat_fdr, new_canonical = padj < 0.05 & abs(logFC) > 0.5,
                                new_treat = treat_fdr < 0.05)],
                        by = "gene", all = TRUE)
    comparison[, symbol := fcoalesce(new_symbol, old_symbol)][, c("old_symbol", "new_symbol") := NULL]
    comparison[, delta_logFC := new_logFC - old_logFC]
    comparison[, canonical_change := fcase(
      is.na(old_canonical), fifelse(new_canonical, "gained_new_gene", "not_called_new_gene"),
      old_canonical & new_canonical & sign(old_logFC) == sign(new_logFC), "retained",
      old_canonical & new_canonical, "retained_sign_flip",
      old_canonical, "lost",
      new_canonical, "gained",
      default = "not_called")]
    comparison[, treat_change := fcase(
      is.na(old_treat), fifelse(new_treat, "gained_new_gene", "not_called_new_gene"),
      old_treat & new_treat & sign(old_logFC) == sign(new_logFC), "retained",
      old_treat & new_treat, "retained_sign_flip",
      old_treat, "lost",
      new_treat, "gained",
      default = "not_called")]
    setcolorder(comparison, c("gene", "symbol"))
    write_tsv_once(comparison, file.path(arm_dir, "old_vs_new_gene_table.tsv.gz"))
    common <- comparison[!is.na(old_logFC) & !is.na(new_logFC)]
    summary <- rbind(
      comparison[, .N, by = .(change = canonical_change)][, family := "canonical_padj0.05_absLFC0.5"],
      comparison[, .N, by = .(change = treat_change)][, family := "TREAT_lfc0.25"]
    )
    write_tsv_once(summary, file.path(arm_dir, "old_vs_new_summary.tsv"))
    metrics <- data.table(
      metric = c("n_common_genes", "logFC_spearman_common", "logFC_pearson_common",
                 "median_abs_delta_logFC", "max_abs_delta_logFC", "sign_agreement_common"),
      value = c(nrow(common), cor(common$old_logFC, common$new_logFC, method = "spearman"),
                cor(common$old_logFC, common$new_logFC), median(abs(common$delta_logFC)),
                max(abs(common$delta_logFC)), mean(sign(common$old_logFC) == sign(common$new_logFC)))
    )
    write_tsv_once(metrics, file.path(arm_dir, "old_vs_new_metrics.tsv"))
    print(summary); print(metrics)

    per_cohort <- lapply(split(cohort, by = "dataset"), function(r) list(control = r$new_control, disease = r$new_disease))
    expected <- c(
      list(arm = arm,
           control_rule = ARM_RULE[[arm]],
           excluded_samples = ARM_EXCLUDED[[arm]],
           n_samples = ncol(fitted$dge), n_cohorts = nlevels(fitted$info$dataset), n_genes = nrow(out),
           n_control = sum(fitted$info$group_binary == "Control"),
           n_disease = sum(fitted$info$group_binary == "Disease"),
           per_cohort = per_cohort),
      counts,
      downstream_label_counts(meta_arm[sample_id %in% arm_ids])
    )
    write_json_once(expected, file.path(arm_dir, "expected_counts.json"))
    arm_count_rows[[arm]] <- data.table(
      arm = arm, n_samples = ncol(fitted$dge), n_control = expected$n_control, n_disease = expected$n_disease,
      gse130970_control = cohort[dataset == "GSE130970", new_control], n_genes = nrow(out),
      canonical_n = counts$canonical$n, canonical_up = counts$canonical$up, canonical_down = counts$canonical$down,
      treat_n = counts$treat$n, treat_up = counts$treat$up, treat_down = counts$treat$down)
    rm(dge, fitted); invisible(gc())
  }
  write_tsv_once(rbindlist(design_rows), file.path(OUT_ROOT, "design_rank.tsv"))
  write_tsv_once(rbindlist(cohort_rows), file.path(OUT_ROOT, "cohort_group_counts_old_vs_new.tsv"))
  arm_counts <- rbind(
    data.table(arm = "adopted_2026-08-10", n_samples = 844L, n_control = 157L, n_disease = 687L,
               gse130970_control = 23L, n_genes = nrow(reference_deg),
               canonical_n = ref_counts$canonical$n, canonical_up = ref_counts$canonical$up,
               canonical_down = ref_counts$canonical$down, treat_n = ref_counts$treat$n,
               treat_up = ref_counts$treat$up, treat_down = ref_counts$treat$down),
    rbindlist(arm_count_rows))
  write_tsv_once(arm_counts, file.path(OUT_ROOT, "arm_deg_counts.tsv"))
  print(arm_counts)
}

# ================================================================= reproduce
if (MODE == "reproduce") {
  SOURCE_ROOT <- normalizePath(Sys.getenv("REFIT_SOURCE_ROOT", ""), mustWork = TRUE)
  if (!file.exists(file.path(SOURCE_ROOT, "REFIT_COMPLETE.json"))) fail("Fit-mode root is incomplete: ", SOURCE_ROOT)
  if (basename(OUT_ROOT) != paste0(basename(SOURCE_ROOT), "-reproduction")) fail("Reproduction root must pair with its fit root")
  source_manifest <- fread(file.path(SOURCE_ROOT, "output_manifest.tsv"))
  for (arm in ARMS) {
    cat("\n=== reproduce ", arm, " ===\n", sep = "")
    src_dir <- file.path(SOURCE_ROOT, arm)
    arm_dir <- file.path(OUT_ROOT, arm)
    dir.create(arm_dir, recursive = FALSE, showWarnings = FALSE)
    for (f in c("merged_dge.rds", "meta_matched.rds", "deg_results.csv", "expected_counts.json")) {
      recorded <- source_manifest[relative_path == file.path(arm, f), sha256]
      if (length(recorded) != 1L || sha256_file(file.path(src_dir, f)) != recorded) fail("Fit artifact changed: ", arm, "/", f)
    }
    dge <- readRDS(file.path(src_dir, "merged_dge.rds"))
    meta_arm <- as.data.table(readRDS(file.path(src_dir, "meta_matched.rds")))
    fitted <- fit_disease_model(dge, meta_arm)
    write_fit_outputs(fitted, arm_dir, file.path(src_dir, "merged_dge.rds"), file.path(src_dir, "meta_matched.rds"),
                      reference_path = file.path(src_dir, "deg_results.csv"))
    out <- fitted$out
    fit_deg <- fread(file.path(src_dir, "deg_results.csv"))
    expected <- fromJSON(file.path(src_dir, "expected_counts.json"), simplifyVector = TRUE)
    counts <- deg_counts(out)
    aligned <- identical(out$gene, fit_deg$gene)
    checks <- data.table(
      check = c("gene_order_identical", "n_genes", "n_samples", "n_control", "n_disease",
                "canonical_n", "treat_n", "treat_up", "treat_down",
                "logFC_max_abs_delta", "SE_max_abs_delta", "padj_max_abs_delta",
                "treat_fdr_max_abs_delta", "lfsr_max_abs_delta",
                "canonical_membership_identical", "treat_membership_identical", "design_full_rank"),
      observed = c(aligned, nrow(out), ncol(fitted$dge),
                   sum(fitted$info$group_binary == "Control"), sum(fitted$info$group_binary == "Disease"),
                   counts$canonical$n, counts$treat$n, counts$treat$up, counts$treat$down,
                   if (aligned) max_abs_delta(out$logFC, fit_deg$logFC) else NA,
                   if (aligned) max_abs_delta(out$SE, fit_deg$SE) else NA,
                   if (aligned) max_abs_delta(out$padj, fit_deg$padj) else NA,
                   if (aligned) max_abs_delta(out$treat_fdr, fit_deg$treat_fdr) else NA,
                   if (aligned) max_abs_delta(out$lfsr, fit_deg$lfsr) else NA,
                   if (aligned) identical(out$padj < 0.05 & abs(out$logFC) > 0.5,
                                          fit_deg$padj < 0.05 & abs(fit_deg$logFC) > 0.5) else FALSE,
                   if (aligned) identical(out$treat_fdr < 0.05, fit_deg$treat_fdr < 0.05) else FALSE,
                   fitted$design_rank == ncol(fitted$design)),
      expected = c(1, expected$n_genes, expected$n_samples, expected$n_control, expected$n_disease,
                   expected$canonical$n, expected$treat$n, expected$treat$up, expected$treat$down,
                   LOGFC_TOLERANCE, 1e-8, 1e-8, 1e-8, 1e-6, 1, 1, 1),
      rule = c(rep("equal", 9), rep("max", 5), rep("equal", 3))
    )
    checks[, pass := !is.na(observed) & fifelse(rule == "equal", observed == expected, observed <= expected)]
    print(checks)
    write_tsv_once(checks, file.path(arm_dir, "validation/validation_checks.tsv"))
    if (!all(checks$pass)) fail("Reproduction FAILED for ", arm)
    write_json_once(list(status = "PASS", arm = arm, fit_root = SOURCE_ROOT,
                         checks = nrow(checks), tolerance_logFC = LOGFC_TOLERANCE),
                    file.path(arm_dir, "validation/validation_report.json"))
    if (!file.copy(file.path(src_dir, "expected_counts.json"), file.path(arm_dir, "expected_counts.json"))) {
      fail("Could not copy expected_counts.json for ", arm)
    }
    writeLines(capture.output(sessionInfo()), file.path(arm_dir, "sessionInfo.txt"))
    write_json_once(list(status = "VALIDATED", arm = arm, fit_root = SOURCE_ROOT,
                         created_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)),
                    file.path(arm_dir, "VALIDATED.json"))
    # build_bulk_release.py validates every file listed here.
    files <- setdiff(list.files(arm_dir, recursive = TRUE), "output_manifest.tsv")
    write_tsv_once(data.table(relative_path = files,
                              size_bytes = file.info(file.path(arm_dir, files))$size,
                              sha256 = vapply(file.path(arm_dir, files), sha256_file, character(1))),
                   file.path(arm_dir, "output_manifest.tsv"))
    rm(dge, fitted); invisible(gc())
  }
}

# ================================================================= provenance
writeLines(capture.output(sessionInfo()), file.path(OUT_ROOT, "sessionInfo.txt"))
if (MODE == "fit") {
  files <- setdiff(list.files(OUT_ROOT, recursive = TRUE), c("output_manifest.tsv", "README.md"))
  output_manifest <- data.table(relative_path = files, size_bytes = file.info(file.path(OUT_ROOT, files))$size)
  output_manifest[, sha256 := vapply(file.path(OUT_ROOT, relative_path), sha256_file, character(1))]
  write_tsv_once(output_manifest, file.path(OUT_ROOT, "output_manifest.tsv"))
  readme <- c(
    paste0("# ", basename(OUT_ROOT)),
    "",
    "Five-cohort F_five disease-versus-control refit with corrected GSE130970 control labels",
    "(review item 1, docs/technical/SCIENTIFIC_REVIEW_CORRECTIONS_2026-09-23.md).",
    "Candidate only. Nothing here is adopted.",
    "",
    paste0("Producer: scripts/manuscript/resource_f_five/refit_corrected_controls.R (sha256 ",
           sha256_file(file.path(PROJECT_ROOT, "scripts/manuscript/resource_f_five/refit_corrected_controls.R")), ")"),
    paste0("SLURM job: ", Sys.getenv("SLURM_JOB_ID", "not_slurm"), " on ", Sys.info()[["nodename"]]),
    paste0("Created (UTC): ", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    "",
    "Model: filterByExpr(~ 0 + group_binary), RLE, voomWithQualityWeights,",
    "~ dataset + inferred_sex + group_binary, eBayes; canonical DEG = padj < 0.05 and |logFC| > 0.5;",
    "TREAT lfc 0.25 (BH) as sensitivity; ashr normal mixture. Unit: one biopsy per participant, n = 844.",
    "",
    "Arms:",
    "- reproduction_guard: old labels; must reproduce the adopted F_five fit (guard_checks.tsv).",
    "- primary_source_controls: 6 GSE130970 source controls (Hoang 2019); 140 control / 704 disease.",
    "- sensitivity_strict_nas0: 4 GSE130970 controls with steatosis 0 and NAS 0; 138 / 706.",
    "- sensitivity_drop_sex_discordant: primary labels without SRR9036338, SRR9036344 and SRR9036376",
    "  (GEO sex disagrees with XIST/Y inferred sex); 841 participants, 139 / 702.",
    "  sex_discordance_844.tsv lists every reported/inferred sex mismatch in the 844.",
    "Only primary_source_controls feeds the F0 arms, bulk figures and continuum.",
    "arm_deg_counts.tsv puts the DEG counts of every arm side by side.",
    "",
    "## Inputs",
    "",
    "| input | sha256 | path |",
    "|---|---|---|",
    input_manifest[, sprintf("| %s | %s | %s |", input_id, sha256, path)],
    "",
    "## Outputs",
    "",
    "Every file with its SHA-256 is listed in output_manifest.tsv."
  )
  writeLines(readme, file.path(OUT_ROOT, "README.md"))
  write_json_once(list(status = "REFIT_COMPLETE", root = OUT_ROOT, arms = ARMS,
                       created_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)),
                  file.path(OUT_ROOT, "REFIT_COMPLETE.json"))
}
cat("\nDONE ", MODE, ": ", OUT_ROOT, "\n", sep = "")
