#!/usr/bin/env Rscript
# Reproduce the accepted F_five five-cohort Disease-vs-Control fit in an
# isolated Resource-candidate workstream. This script has no canonical path.

suppressPackageStartupMessages({
  library(ashr)
  library(data.table)
  library(edgeR)
  library(jsonlite)
  library(limma)
  library(yaml)
})

options(digits = 17, scipen = 999)

PROJECT_ROOT <- normalizePath(
  Sys.getenv(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
  ),
  mustWork = TRUE
)
CANDIDATE_ID <- "resource-f-five-coloc-v6-candidate-2026-08-10"
CANDIDATE_ROOT <- normalizePath(
  Sys.getenv(
    "MASLD_RESOURCE_CANDIDATE_ROOT",
    file.path(
      PROJECT_ROOT,
      "RNA-seq/results/manuscript_release/candidates",
      CANDIDATE_ID
    )
  ),
  mustWork = TRUE
)
if (basename(CANDIDATE_ROOT) != CANDIDATE_ID) {
  stop("Candidate ID/path disagreement")
}
sentinel <- file.path(CANDIDATE_ROOT, ".resource_candidate_root")
if (!file.exists(sentinel) ||
    !identical(trimws(readLines(sentinel, warn = FALSE)), CANDIDATE_ID)) {
  stop("Missing or invalid Resource-candidate sentinel")
}
if (!file.exists(file.path(CANDIDATE_ROOT, "BASE_SNAPSHOT_COMPLETE.json"))) {
  stop("BASE snapshot has not completed")
}

MODEL_INPUT_ROOT <- file.path(CANDIDATE_ROOT, "inputs/BULK-F-FIVE/frozen_model_inputs")
ARM_ROOT <- file.path(CANDIDATE_ROOT, "inputs/BG001-DECISION/arms/F_five")
DGE_PATH <- file.path(ARM_ROOT, "results/integration/merged_dge.rds")
META_PATH <- file.path(ARM_ROOT, "results/integration/meta_matched.rds")
REFERENCE_PATH <- file.path(ARM_ROOT, "results/integration/deg_results.csv")
CONFIG_PATH <- file.path(MODEL_INPUT_ROOT, "human_datasets.yaml")
GENE_METADATA_PATH <- file.path(
  MODEL_INPUT_ROOT, "gencode_v49_gene_metadata.tsv.gz"
)
PUBLISHER_PATH <- file.path(CANDIDATE_ROOT, "code/publish_noreplace.py")
SNAPSHOT_IO_PATH <- file.path(CANDIDATE_ROOT, "code/snapshot_io.py")
PYTHON_PATH <- normalizePath(Sys.getenv("MASLD_PYTHON_PATH"), mustWork = TRUE)
FINAL_ROOT <- file.path(CANDIDATE_ROOT, "workstreams/BULK-POOLED-REPRO")
job_tag <- Sys.getenv("SLURM_JOB_ID", as.character(Sys.getpid()))
TMP_ROOT <- paste0(FINAL_ROOT, ".tmp.", job_tag)

is_symlink <- function(path) {
  target <- Sys.readlink(path)
  !is.na(target) && nzchar(target)
}
for (path in c(
  DGE_PATH, META_PATH, REFERENCE_PATH, CONFIG_PATH, GENE_METADATA_PATH,
  PUBLISHER_PATH, SNAPSHOT_IO_PATH
)) {
  if (!file.exists(path) || is_symlink(path)) {
    stop("Missing or symlinked frozen input: ", path)
  }
}
if (file.exists(FINAL_ROOT) || is_symlink(FINAL_ROOT) ||
    file.exists(TMP_ROOT) || is_symlink(TMP_ROOT)) {
  stop("Refusing existing pooled-reproduction output")
}

sha256_file <- function(path) {
  value <- system2(
    "sha256sum", shQuote(normalizePath(path, mustWork = TRUE)),
    stdout = TRUE, stderr = TRUE
  )
  status <- attr(value, "status")
  if ((!is.null(status) && status != 0L) || length(value) != 1L) {
    stop("sha256sum failed for ", path)
  }
  digest <- strsplit(value[[1L]], "[[:space:]]+")[[1L]][[1L]]
  if (!grepl("^[0-9a-f]{64}$", digest)) stop("Invalid SHA256 for ", path)
  digest
}

expected_hashes <- c(
  merged_dge = "bc3ea19c8b064861339eba862abf58e82ddb5460c622b550dc45aa3b5fcd1f76",
  meta_matched = "90ba6aca68643c4680c08fc982a6c72b76edb740ab726874cb50603711756a83",
  accepted_deg = "691b09517f42ff24e472fdb9db460b1c339ec734876c40115f8b146293b5b633"
)
observed_hashes <- c(
  merged_dge = sha256_file(DGE_PATH),
  meta_matched = sha256_file(META_PATH),
  accepted_deg = sha256_file(REFERENCE_PATH)
)
if (!identical(observed_hashes, expected_hashes)) {
  stop("Frozen F_five input identity mismatch")
}
RUNTIME_CONTRACT_PATH <- file.path(
  CANDIDATE_ROOT,
  "inputs/BG001-DECISION/contract/analysis_runtime_contract.json"
)
if (!file.exists(RUNTIME_CONTRACT_PATH) || is_symlink(RUNTIME_CONTRACT_PATH)) {
  stop("Missing or symlinked accepted runtime contract")
}
runtime_contract <- jsonlite::fromJSON(RUNTIME_CONTRACT_PATH, simplifyVector = FALSE)
if (!identical(runtime_contract$schema, "bg001-analysis-runtime-v1")) {
  stop("Accepted runtime contract schema drift")
}
runtime_packages <- setNames(
  vapply(runtime_contract$R_runtime$packages, `[[`, character(1), "version"),
  vapply(runtime_contract$R_runtime$packages, `[[`, character(1), "name")
)
required_packages <- c("ashr", "data.table", "edgeR", "jsonlite", "limma", "yaml")
observed_versions <- vapply(
  required_packages,
  function(package) as.character(packageVersion(package)),
  character(1)
)
if (!identical(observed_versions, runtime_packages[required_packages]) ||
    !identical(as.character(getRversion()), runtime_packages[["base"]])) {
  stop("R or package-version drift from accepted runtime")
}
rscript_path <- normalizePath(Sys.getenv("MASLD_RSCRIPT_PATH"), mustWork = TRUE)
if (!identical(rscript_path, normalizePath(runtime_contract$Rscript$path, mustWork = TRUE)) ||
    !identical(sha256_file(rscript_path), runtime_contract$Rscript$sha256)) {
  stop("Rscript executable drift from accepted runtime")
}
thread_variables <- c(
  "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
  "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"
)
if (!identical(Sys.getenv("SLURM_CPUS_PER_TASK"), "8") ||
    any(vapply(
      thread_variables,
      function(variable) !identical(Sys.getenv(variable), "8"),
      logical(1)
    ))) {
  stop("Thread-count contract drift; exact reproduction requires eight threads")
}

if (!isTRUE(dir.create(
  TMP_ROOT, recursive = FALSE, showWarnings = FALSE, mode = "0750"
))) {
  stop("Failed to create unique pooled temporary root")
}

write_table <- function(object, path, sep = "\t", compress = FALSE) {
  connection_mode <- if (compress) "wbx" else "wx"
  connection <- file(path, open = connection_mode)
  if (compress) {
    connection <- gzcon(connection, level = 6, text = TRUE)
  }
  on.exit(try(close(connection), silent = TRUE), add = TRUE)
  write.table(
    as.data.frame(object), connection, sep = sep, row.names = FALSE,
    col.names = TRUE, quote = FALSE, na = "NA", eol = "\n"
  )
  close(connection)
}

write_lines_once <- function(value, path) {
  connection <- file(path, open = "wx")
  on.exit(try(close(connection), silent = TRUE), add = TRUE)
  writeLines(value, connection, useBytes = TRUE)
  close(connection)
}

save_rds_once <- function(value, path) {
  connection <- file(path, open = "wxb")
  on.exit(try(close(connection), silent = TRUE), add = TRUE)
  saveRDS(value, connection)
  close(connection)
}

dge <- readRDS(DGE_PATH)
meta <- as.data.table(readRDS(META_PATH))
config <- yaml::read_yaml(CONFIG_PATH)$datasets
mega <- names(Filter(function(item) isTRUE(item$de$include_in_mega), config))
if (length(mega) != 5L) stop("Frozen config does not define exactly five mega cohorts")
keep <- dge$samples$dataset %in% mega
dge_mega <- dge[, keep]
sex <- meta$inferred_sex[match(colnames(dge_mega), meta$sample_id)]
if (anyNA(sex)) stop("Missing inferred sex in pooled model")
info <- data.frame(
  group_binary = factor(
    dge_mega$samples$group_binary,
    levels = c("Control", "Disease")
  ),
  dataset = factor(dge_mega$samples$dataset),
  inferred_sex = factor(sex),
  row.names = colnames(dge_mega)
)
if (nrow(dge_mega) != 23370L ||
    ncol(dge_mega) != 844L ||
    nlevels(info$dataset) != 5L ||
    sum(info$group_binary == "Control") != 157L ||
    sum(info$group_binary == "Disease") != 687L) {
  stop("F_five pooled model census drift")
}

design <- model.matrix(~ dataset + inferred_sex + group_binary, data = info)
coefficient <- "group_binaryDisease"
if (!coefficient %in% colnames(design) ||
    qr(design)$rank != ncol(design)) {
  stop("Pooled C2 design is not full-rank and estimable")
}

v <- voomWithQualityWeights(dge_mega, design, save.plot = TRUE)
sample_weights <- v$targets$sample.weights
if (is.null(sample_weights) ||
    length(sample_weights) != ncol(dge_mega) ||
    !identical(rownames(v$targets), colnames(dge_mega)) ||
    any(!is.finite(sample_weights)) ||
    any(sample_weights <= 0)) {
  stop("Invalid voom sample weights")
}
validate_voom_curve <- function(curve, label) {
  if (is.null(curve) || is.null(curve$x) || is.null(curve$y) ||
      !is.numeric(curve$x) || !is.numeric(curve$y) ||
      length(curve$x) == 0L || length(curve$x) != length(curve$y) ||
      any(!is.finite(curve$x)) || any(!is.finite(curve$y))) {
    stop("voomWithQualityWeights returned an invalid ", label, " curve")
  }
}
validate_voom_curve(v$voom.xy, "observed mean-variance")
validate_voom_curve(v$voom.line, "fitted mean-variance")

fit0 <- lmFit(v, design)
fit <- eBayes(fit0)
res <- topTable(fit, coef = coefficient, number = Inf, sort.by = "none")
res$gene <- rownames(res)
out <- as.data.table(res)
setnames(out, "adj.P.Val", "padj")
coef_index <- match(coefficient, colnames(fit$coefficients))
moderated_se <- fit$stdev.unscaled[, coef_index] * sqrt(fit$s2.post)
names(moderated_se) <- rownames(fit$coefficients)
out[, SE := moderated_se[match(gene, names(moderated_se))]]

treat_lfc <- 0.25
treat_table <- topTreat(
  treat(fit0, lfc = treat_lfc),
  coef = coefficient,
  number = Inf,
  sort.by = "none"
)
out[, `:=`(
  treat_lfc = treat_lfc,
  treat_p = treat_table[gene, "P.Value"],
  treat_fdr = treat_table[gene, "adj.P.Val"]
)]

ok <- is.finite(out$logFC) & is.finite(out$SE) & out$SE > 0
ash <- ashr::ash(out$logFC[ok], out$SE[ok], mixcompdist = "normal")
out[, shrunk_logFC := NA_real_][ok, shrunk_logFC := ash$result$PosteriorMean]
out[, lfsr := NA_real_][ok, lfsr := ash$result$lfsr]
out <- out[, .(
  gene, logFC, SE, t, P.Value, padj, shrunk_logFC, lfsr,
  treat_lfc, treat_p, treat_fdr, AveExpr
)]

genes <- fread(GENE_METADATA_PATH)
genes[, ensembl_base := sub("[.][0-9]+$", "", gene_id)]
out[, symbol := genes[
  match(sub("[.][0-9]+$", "", gene), ensembl_base),
  gene_name
]]

numeric_fields <- c(
  "logFC", "SE", "t", "P.Value", "padj", "treat_lfc",
  "treat_p", "treat_fdr", "AveExpr"
)
probability_fields <- c("P.Value", "padj", "treat_p", "treat_fdr")
if (!identical(out$gene, rownames(dge_mega)) ||
    anyDuplicated(out$gene) ||
    any(vapply(numeric_fields, function(field) {
      any(!is.finite(out[[field]]))
    }, logical(1))) ||
    any(out$SE <= 0) ||
    any(out$treat_lfc != 0.25) ||
    any(vapply(probability_fields, function(field) {
      any(out[[field]] < 0 | out[[field]] > 1)
    }, logical(1)))) {
  stop("Pooled coefficient table failed numerical validation")
}
selected <- out[treat_fdr < 0.05]
if (nrow(selected) != 1616L ||
    sum(selected$logFC > 0) != 1144L ||
    sum(selected$logFC < 0) != 472L) {
  stop("Pooled TREAT membership does not reproduce accepted headline counts")
}

write_table(out, file.path(TMP_ROOT, "deg_results.csv"), sep = ",")
write_table(
  out[, .(gene, symbol, logFC, SE, t, P.Value, padj, AveExpr)],
  file.path(TMP_ROOT, "coefficient_table.tsv.gz"),
  compress = TRUE
)
write_table(
  out[, .(
    gene, symbol, logFC, SE, AveExpr,
    treat_lfc, treat_p, treat_fdr
  )],
  file.path(TMP_ROOT, "treat_table.tsv.gz"),
  compress = TRUE
)
design_table <- as.data.table(design)
design_table[, sample_id := colnames(dge_mega)]
setcolorder(design_table, c("sample_id", setdiff(names(design_table), "sample_id")))
write_table(design_table, file.path(TMP_ROOT, "model_design.tsv"))
write_table(
  data.table(
    sample_id = colnames(dge_mega),
    sample_weight = as.numeric(sample_weights)
  ),
  file.path(TMP_ROOT, "lvqw_sample_weights.tsv")
)
save_rds_once(
  list(voom_xy = v$voom.xy, voom_line = v$voom.line),
  file.path(TMP_ROOT, "voom_curve.rds")
)
write_table(
  data.table(
    mean_log_count = v$voom.xy$x,
    sqrt_residual_sd = v$voom.xy$y
  ),
  file.path(TMP_ROOT, "voom_mean_variance.tsv")
)
write_table(
  data.table(
    dge_path = normalizePath(DGE_PATH, mustWork = TRUE),
    meta_path = normalizePath(META_PATH, mustWork = TRUE),
    config_path = normalizePath(CONFIG_PATH, mustWork = TRUE),
    gene_metadata_path = normalizePath(GENE_METADATA_PATH, mustWork = TRUE),
    accepted_reference_path = normalizePath(REFERENCE_PATH, mustWork = TRUE),
    n_genes = nrow(dge_mega),
    n_samples = ncol(dge_mega),
    n_cohorts = nlevels(info$dataset),
    n_control = sum(info$group_binary == "Control"),
    n_disease = sum(info$group_binary == "Disease"),
    treat_lfc = treat_lfc
  ),
  file.path(TMP_ROOT, "model_input_manifest.tsv")
)
write_table(
  data.table(
    component = c("R", "data.table", "edgeR", "limma", "ashr", "yaml", "jsonlite"),
    version = c(
      as.character(getRversion()),
      vapply(
        c("data.table", "edgeR", "limma", "ashr", "yaml", "jsonlite"),
        function(package) as.character(packageVersion(package)),
        character(1)
      )
    ),
    hostname = Sys.info()[["nodename"]],
    locale = Sys.getlocale(),
    slurm_job_id = Sys.getenv("SLURM_JOB_ID", "not_slurm"),
    slurm_partition = Sys.getenv("SLURM_JOB_PARTITION", "not_slurm"),
    slurm_qos = Sys.getenv("SLURM_JOB_QOS", "not_slurm"),
    slurm_nodelist = Sys.getenv("SLURM_JOB_NODELIST", "not_slurm"),
    slurm_cpus_per_task = Sys.getenv("SLURM_CPUS_PER_TASK", "not_slurm"),
    omp_num_threads = Sys.getenv("OMP_NUM_THREADS", "not_set"),
    openblas_num_threads = Sys.getenv("OPENBLAS_NUM_THREADS", "not_set"),
    mkl_num_threads = Sys.getenv("MKL_NUM_THREADS", "not_set"),
    veclib_maximum_threads = Sys.getenv("VECLIB_MAXIMUM_THREADS", "not_set"),
    numexpr_num_threads = Sys.getenv("NUMEXPR_NUM_THREADS", "not_set")
  ),
  file.path(TMP_ROOT, "environment_manifest.tsv")
)
write_lines_once(
  capture.output(sessionInfo()),
  file.path(TMP_ROOT, "sessionInfo.txt")
)

producer_path <- grep(
  "^--file=", commandArgs(trailingOnly = FALSE), value = TRUE
)
if (length(producer_path) != 1L) stop("Cannot resolve executed producer path")
producer_path <- sub("^--file=", "", producer_path)
write_table(
  data.table(
    artifact_id = c(
      "producer", "merged_dge", "meta_matched", "config",
      "gene_metadata", "accepted_reference", "runtime_contract", "Rscript",
      "atomic_publisher", "atomic_publisher_io", "Python"
    ),
    role = c(
      "resource_candidate_code",
      "bg001_f_five_artifact",
      "bg001_f_five_artifact",
      "pooled_config",
      "gene_annotation",
      "bg001_f_five_artifact",
      "bg001_analysis_bound_input;bg001_baseline_transaction",
      "analysis_runtime_rscript",
      "resource_candidate_code",
      "resource_candidate_code",
      "resource_candidate_python"
    ),
    path = c(
      normalizePath(producer_path, mustWork = TRUE),
      normalizePath(DGE_PATH, mustWork = TRUE),
      normalizePath(META_PATH, mustWork = TRUE),
      normalizePath(CONFIG_PATH, mustWork = TRUE),
      normalizePath(GENE_METADATA_PATH, mustWork = TRUE),
      normalizePath(REFERENCE_PATH, mustWork = TRUE),
      normalizePath(RUNTIME_CONTRACT_PATH, mustWork = TRUE),
      rscript_path,
      normalizePath(PUBLISHER_PATH, mustWork = TRUE),
      normalizePath(SNAPSHOT_IO_PATH, mustWork = TRUE),
      PYTHON_PATH
    ),
    sha256 = c(
      sha256_file(producer_path),
      sha256_file(DGE_PATH),
      sha256_file(META_PATH),
      sha256_file(CONFIG_PATH),
      sha256_file(GENE_METADATA_PATH),
      sha256_file(REFERENCE_PATH),
      sha256_file(RUNTIME_CONTRACT_PATH),
      sha256_file(rscript_path),
      sha256_file(PUBLISHER_PATH),
      sha256_file(SNAPSHOT_IO_PATH),
      sha256_file(PYTHON_PATH)
    )
  )[, .(artifact_id, path, role, sha256)],
  file.path(TMP_ROOT, "input_code_manifest.tsv")
)

artifact_paths <- sort(list.files(
  TMP_ROOT, recursive = TRUE, full.names = TRUE, all.files = TRUE,
  no.. = TRUE
))
artifact_paths <- artifact_paths[file.info(artifact_paths)$isdir == FALSE]
artifact_manifest <- data.table(
  artifact_path = substring(
    artifact_paths,
    nchar(normalizePath(TMP_ROOT, mustWork = TRUE)) + 2L
  ),
  bytes = file.info(artifact_paths)$size,
  sha256 = vapply(artifact_paths, sha256_file, character(1))
)
write_table(
  artifact_manifest,
  file.path(TMP_ROOT, "producer_artifact_manifest.tsv")
)
completion <- list(
  status = "PRODUCER_COMPLETE_PENDING_VALIDATION",
  candidate_id = CANDIDATE_ID,
  workstream = "BULK-POOLED-REPRO",
  n_genes = nrow(out),
  n_samples = ncol(dge_mega),
  n_treat = nrow(selected),
  n_up = sum(selected$logFC > 0),
  n_down = sum(selected$logFC < 0),
  producer_artifact_manifest_sha256 = sha256_file(
    file.path(TMP_ROOT, "producer_artifact_manifest.tsv")
  ),
  canonical_write = FALSE,
  created_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
)
write_lines_once(
  jsonlite::toJSON(completion, auto_unbox = TRUE, pretty = TRUE),
  file.path(TMP_ROOT, "PRODUCER_COMPLETE.json")
)
if (file.exists(FINAL_ROOT) || is_symlink(FINAL_ROOT)) {
  stop("Final pooled root appeared during execution")
}
publication <- system2(
  PYTHON_PATH,
  c(
    "-S",
    shQuote(normalizePath(PUBLISHER_PATH, mustWork = TRUE)),
    "--marker",
    "PRODUCER_COMPLETE.json",
    shQuote(normalizePath(TMP_ROOT, mustWork = TRUE)),
    shQuote(FINAL_ROOT)
  ),
  stdout = TRUE,
  stderr = TRUE
)
publication_status <- attr(publication, "status")
source_marker <- file.path(TMP_ROOT, "PRODUCER_COMPLETE.json")
official_marker <- file.path(FINAL_ROOT, "PRODUCER_COMPLETE.json")
if ((!is.null(publication_status) && publication_status != 0L) ||
    !dir.exists(FINAL_ROOT) || is_symlink(FINAL_ROOT) ||
    !dir.exists(TMP_ROOT) || is_symlink(TMP_ROOT) ||
    !file_test("-f", source_marker) || is_symlink(source_marker) ||
    !file_test("-f", official_marker) || is_symlink(official_marker) ||
    !identical(sha256_file(source_marker), sha256_file(official_marker))) {
  stop("Staged no-replace publication failed: ", paste(publication, collapse = " | "))
}
cat("PRODUCER_COMPLETE:", FINAL_ROOT, "\n")
