#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
})

options(digits = 17, scipen = 999)
set.seed(20260811)

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_argument) != 1L) stop("Cannot resolve producer path", call. = FALSE)
PRODUCER <- normalizePath(sub("^--file=", "", script_argument), mustWork = TRUE)
SCRIPT_DIR <- dirname(PRODUCER)
source(file.path(SCRIPT_DIR, "bulk_lncrna_common.R"), local = TRUE)
source(file.path(SCRIPT_DIR, "lncrna_biotype_association_common.R"), local = TRUE)

PROJECT_ROOT <- normalizePath(
  Sys.getenv(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
  ),
  mustWork = TRUE
)
CANDIDATE_ROOT <- file.path(
  PROJECT_ROOT,
  "RNA-seq/results/manuscript_release/candidates/",
  "resource-f-five-coloc-v6-candidate-2026-08-10"
)
ARM_ROOT <- file.path(CANDIDATE_ROOT, "inputs/BG001-DECISION/arms/F_five")
DGE_PATH <- file.path(ARM_ROOT, "results/integration/merged_dge.rds")
PRIMARY_PATH <- file.path(
  CANDIDATE_ROOT, "workstreams/BULK-POOLED-REPRO/deg_results.csv"
)

parse_arguments <- function(values) {
  expected <- c("source-gate", "gene-identity", "output")
  parsed <- list()
  for (value in values) {
    pieces <- strsplit(value, "=", fixed = TRUE)[[1L]]
    if (length(pieces) != 2L || !startsWith(pieces[[1L]], "--")) {
      fail("Arguments must use --name=value syntax: ", value)
    }
    name <- sub("^--", "", pieces[[1L]])
    if (!name %in% expected || name %in% names(parsed)) {
      fail("Unknown or duplicate argument: ", name)
    }
    parsed[[name]] <- pieces[[2L]]
  }
  missing <- setdiff(expected, names(parsed))
  if (length(missing)) fail("Missing arguments: ", paste(missing, collapse = ","))
  parsed
}

sha256_file <- function(path) {
  output <- system2("sha256sum", normalizePath(path, mustWork = TRUE), stdout = TRUE)
  strsplit(output[[1L]], "[[:space:]]+")[[1L]][[1L]]
}

write_tsv_once <- function(value, path) {
  connection <- file(path, open = "wx")
  on.exit(try(close(connection), silent = TRUE), add = TRUE)
  write.table(
    as.data.frame(value), connection, sep = "\t", row.names = FALSE,
    col.names = TRUE, quote = FALSE, na = "NA", eol = "\n"
  )
  close(connection)
}

arguments <- parse_arguments(commandArgs(trailingOnly = TRUE))
SOURCE_GATE <- normalizePath(arguments[["source-gate"]], mustWork = TRUE)
IDENTITY_PATH <- normalizePath(arguments[["gene-identity"]], mustWork = TRUE)
OUTPUT <- file.path(
  normalizePath(dirname(arguments[["output"]]), mustWork = TRUE),
  basename(arguments[["output"]])
)
if (file.exists(OUTPUT) || is_symlink(OUTPUT)) fail("Refusing output overwrite")
read_source_gate(SOURCE_GATE)

expected_hashes <- c(
  dge = "bc3ea19c8b064861339eba862abf58e82ddb5460c622b550dc45aa3b5fcd1f76",
  primary = "b2d2122916401077055e9c6c8019637fd49f2cddd10bd9819be2cc098ea40b84"
)
if (sha256_file(DGE_PATH) != expected_hashes[["dge"]] ||
    sha256_file(PRIMARY_PATH) != expected_hashes[["primary"]]) {
  fail("Frozen bulk input hash drift")
}

dge_all <- readRDS(DGE_PATH)
keep <- dge_all$samples$dataset %in% LNCRNA_COHORTS
dge <- dge_all[, keep, keep.lib.sizes = TRUE]
if (!identical(dim(dge), c(23370L, 844L))) fail("Five-cohort DGE census drift")

identity <- fread(IDENTITY_PATH, colClasses = "character")
required_identity <- c(
  "gene_id_versioned", "gene_type", "start_1based", "end_1based"
)
if (!all(required_identity %in% names(identity)) ||
    anyDuplicated(identity$gene_id_versioned)) {
  fail("GENCODE identity schema drift")
}
identity <- identity[match(rownames(dge), gene_id_versioned)]
if (anyNA(identity$gene_id_versioned)) fail("Incomplete stable-ID identity join")

primary <- fread(PRIMARY_PATH)
if (!identical(primary$gene, rownames(dge)) || nrow(primary) != 23370L) {
  fail("Primary result universe drift")
}
log_cpm <- cpm(dge, log = TRUE, prior.count = 0.5)
expression_variability <- apply(log_cpm, 1L, var)

cohort_detection <- integer(nrow(dge))
for (cohort in LNCRNA_COHORTS) {
  selected <- dge$samples$dataset == cohort
  cohort_group <- droplevels(factor(dge$samples$group_binary[selected]))
  if (nlevels(cohort_group) != 2L) fail("Cohort disease design is not estimable")
  cohort_detection <- cohort_detection + as.integer(
    filterByExpr(dge[, selected, keep.lib.sizes = TRUE], group = cohort_group)
  )
}

gene_table <- data.table(
  gene_id_versioned = rownames(dge),
  gene_type = identity$gene_type,
  treat_positive = primary$treat_fdr < LNCRNA_FDR,
  direction = fifelse(primary$logFC > 0, "up", fifelse(primary$logFC < 0, "down", "zero")),
  AveExpr = primary$AveExpr,
  expression_variability = expression_variability,
  gene_length_bp = as.integer(identity$end_1based) - as.integer(identity$start_1based) + 1L,
  cohort_detection_count = cohort_detection
)
gene_table[, display_biotype := fcase(
  gene_type == "protein_coding", "protein_coding",
  gene_type == "lncRNA", "lncRNA",
  default = "other"
)]
counts <- gene_table[, .(
  n_tested = .N,
  n_treat_positive = sum(treat_positive),
  n_treat_up = sum(treat_positive & direction == "up"),
  n_treat_down = sum(treat_positive & direction == "down")
), by = display_biotype]
counts[, display_order := match(
  display_biotype, c("protein_coding", "lncRNA", "other")
)]
setorder(counts, display_order)
counts[, display_order := NULL]
if (sum(counts$n_tested) != 23370L || sum(counts$n_treat_positive) != 1616L ||
    counts[display_biotype == "lncRNA", n_tested] != 6249L ||
    counts[display_biotype == "lncRNA", n_treat_positive] != 434L) {
  fail("Biotype count rederivation failed")
}

association <- fit_biotype_association(gene_table)
if (!dir.create(OUTPUT, recursive = FALSE, mode = "0750")) fail("Output creation failed")
write_tsv_once(counts, file.path(OUTPUT, "bulk_deg_counts_by_biotype.tsv"))
write_tsv_once(association$result, file.path(OUTPUT, "biotype_deg_association.tsv"))
write_tsv_once(gene_table, file.path(OUTPUT, "biotype_association_gene_covariates.tsv"))
write_tsv_once(
  data.table(
    role = c("source_gate", "gene_identity", "frozen_dge", "validated_primary", "producer", "association_common"),
    path = c(SOURCE_GATE, IDENTITY_PATH, DGE_PATH, PRIMARY_PATH, PRODUCER,
             file.path(SCRIPT_DIR, "lncrna_biotype_association_common.R")),
    sha256 = vapply(c(SOURCE_GATE, IDENTITY_PATH, DGE_PATH, PRIMARY_PATH, PRODUCER,
                      file.path(SCRIPT_DIR, "lncrna_biotype_association_common.R")),
                    sha256_file, character(1))
  ),
  file.path(OUTPUT, "source_manifest.tsv")
)
write_tsv_once(
  data.table(
    component = c("R", "data.table", "edgeR"),
    version = c(
      as.character(getRversion()),
      as.character(packageVersion("data.table")),
      as.character(packageVersion("edgeR"))
    ),
    seed = 20260811L,
    slurm_job_id = Sys.getenv("SLURM_JOB_ID", "not_slurm"),
    biological_unit = "gene_for_secondary_descriptive_association",
    primary_inferential_unit = "human_sample_in_the_upstream_TREAT_model"
  ),
  file.path(OUTPUT, "execution_manifest.tsv")
)
writeLines(capture.output(sessionInfo()), file.path(OUTPUT, "sessionInfo.txt"), useBytes = TRUE)
output_names <- c(
  "bulk_deg_counts_by_biotype.tsv", "biotype_deg_association.tsv",
  "biotype_association_gene_covariates.tsv", "source_manifest.tsv",
  "execution_manifest.tsv", "sessionInfo.txt"
)
write_tsv_once(
  data.table(
    relative_path = output_names,
    size_bytes = file.info(file.path(OUTPUT, output_names))$size,
    sha256 = vapply(file.path(OUTPUT, output_names), sha256_file, character(1))
  ),
  file.path(OUTPUT, "output_manifest.tsv")
)
