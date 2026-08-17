#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
})

options(digits = 17, scipen = 999)
set.seed(20260811)

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_argument) != 1L) stop("Cannot resolve producer path", call. = FALSE)
PRODUCER <- normalizePath(sub("^--file=", "", script_argument), mustWork = TRUE)
source(file.path(dirname(PRODUCER), "bulk_lncrna_common.R"), local = TRUE)

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
META_PATH <- file.path(ARM_ROOT, "results/integration/meta_matched.rds")
PRIMARY_PATH <- file.path(
  CANDIDATE_ROOT, "workstreams/BULK-POOLED-REPRO/deg_results.csv"
)
PRIMARY_VALIDATED <- file.path(
  CANDIDATE_ROOT, "workstreams/BULK-POOLED-REPRO/VALIDATED.json"
)
MODEL_DESIGN_PATH <- file.path(
  CANDIDATE_ROOT, "workstreams/BULK-POOLED-REPRO/model_design.tsv"
)

parse_arguments <- function(values) {
  expected <- c("source-gate", "gene-identity", "lncrna-class", "output")
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
  output <- system2("sha256sum", shQuote(normalizePath(path, mustWork = TRUE)),
                    stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if ((!is.null(status) && status != 0L) || length(output) != 1L) {
    fail("sha256sum failed for ", path)
  }
  digest <- strsplit(output[[1L]], "[[:space:]]+")[[1L]][[1L]]
  if (!grepl("^[0-9a-f]{64}$", digest)) fail("Invalid SHA256 for ", path)
  digest
}

write_tsv_once <- function(value, path) {
  connection <- file(path, open = "wx")
  on.exit(try(close(connection), silent = TRUE), add = TRUE)
  write.table(as.data.frame(value), connection, sep = "\t", row.names = FALSE,
              col.names = TRUE, quote = FALSE, na = "NA", eol = "\n")
  close(connection)
}

arguments <- parse_arguments(commandArgs(trailingOnly = TRUE))
SOURCE_GATE <- normalizePath(arguments[["source-gate"]], mustWork = TRUE)
IDENTITY_PATH <- normalizePath(arguments[["gene-identity"]], mustWork = TRUE)
CLASS_PATH <- normalizePath(arguments[["lncrna-class"]], mustWork = TRUE)
OUTPUT <- normalizePath(dirname(arguments[["output"]]), mustWork = TRUE)
OUTPUT <- file.path(OUTPUT, basename(arguments[["output"]]))
if (file.exists(OUTPUT) || is_symlink(OUTPUT)) {
  fail("Refusing to overwrite output: ", OUTPUT)
}

# This check occurs before any disease outcome or expression object is read.
read_source_gate(SOURCE_GATE)

required_inputs <- c(
  DGE_PATH, META_PATH, PRIMARY_PATH, PRIMARY_VALIDATED, MODEL_DESIGN_PATH,
  IDENTITY_PATH, CLASS_PATH, SOURCE_GATE,
  PRODUCER, file.path(dirname(PRODUCER), "bulk_lncrna_common.R")
)
for (path in required_inputs) {
  if (!file.exists(path) || is_symlink(path)) {
    fail("Missing or symlinked frozen input: ", path)
  }
}
expected_hashes <- c(
  DGE_PATH = "bc3ea19c8b064861339eba862abf58e82ddb5460c622b550dc45aa3b5fcd1f76",
  META_PATH = "90ba6aca68643c4680c08fc982a6c72b76edb740ab726874cb50603711756a83",
  PRIMARY_PATH = "b2d2122916401077055e9c6c8019637fd49f2cddd10bd9819be2cc098ea40b84",
  PRIMARY_VALIDATED = "51467ae61a57b6789e119263aae1fe1b40efd1ea886356187bb89f7ada619095",
  MODEL_DESIGN_PATH = "0e31fc1e0329e0951ffe21b140b74b0de88c8d583b39182adcd420e1a4051a4c"
)
names(expected_hashes) <- c(
  DGE_PATH, META_PATH, PRIMARY_PATH, PRIMARY_VALIDATED, MODEL_DESIGN_PATH
)
for (path in names(expected_hashes)) {
  assert_true(
    identical(sha256_file(path), unname(expected_hashes[[path]])),
    paste0("Frozen input SHA256 drift: ", path)
  )
}
assert_true(
  identical(as.character(getRversion()), "4.4.3") &&
    identical(as.character(packageVersion("edgeR")), "4.4.2") &&
    identical(as.character(packageVersion("limma")), "3.62.2") &&
    identical(as.character(packageVersion("data.table")), "1.18.2.1"),
  "R analysis runtime drift"
)

identity <- fread(IDENTITY_PATH, colClasses = "character")
identity_required <- c(
  "annotation_release", "gene_id_versioned", "gene_id_base", "gene_name",
  "gene_type", "chromosome", "is_canonical_chromosome"
)
assert_true(all(identity_required %in% names(identity)), "Identity schema drift")
assert_true(!anyDuplicated(identity$gene_id_versioned), "Duplicate versioned identity")
assert_true(all(identity$annotation_release == "GENCODE v49"), "Identity release drift")

classes <- fread(CLASS_PATH, colClasses = "character")
class_required <- c(
  "annotation_release", "gene_id_versioned", "lncrna_genomic_class",
  "mapping_status", "main_text_eligible", "main_text_exclusion_reason"
)
assert_true(all(class_required %in% names(classes)), "lncRNA class schema drift")
assert_true(!anyDuplicated(classes$gene_id_versioned), "Duplicate lncRNA class identity")
assert_true(all(classes$annotation_release == "GENCODE v49"), "lncRNA class release drift")

dge <- readRDS(DGE_PATH)
meta <- as.data.table(readRDS(META_PATH))
keep <- dge$samples$dataset %in% LNCRNA_COHORTS
dge_five <- dge[, keep, keep.lib.sizes = TRUE]
assert_true(
  nrow(dge_five) == 23370L && ncol(dge_five) == 844L,
  "Frozen five-cohort DGE census drift"
)
assert_true(!anyDuplicated(rownames(dge_five)), "Duplicate DGE versioned IDs")
sample_info <- data.table(
  sample_id = colnames(dge_five),
  dataset = as.character(dge_five$samples$dataset),
  inferred_sex = meta$inferred_sex[match(colnames(dge_five), meta$sample_id)],
  group_binary = as.character(dge_five$samples$group_binary)
)
assert_true(!anyNA(sample_info), "Missing five-cohort model metadata")
assert_true(
  identical(sort(unique(sample_info$dataset)), sort(LNCRNA_COHORTS)),
  "Five-cohort identity drift"
)
validated_design <- fread(MODEL_DESIGN_PATH, check.names = FALSE)
reconstructed_design <- make_disease_design(sample_info, include_dataset = TRUE)
assert_true(
  identical(validated_design$sample_id, sample_info$sample_id) &&
    identical(
      setdiff(names(validated_design), "sample_id"),
      colnames(reconstructed_design)
    ) &&
    isTRUE(all.equal(
      as.matrix(validated_design[, setdiff(names(validated_design), "sample_id"), with = FALSE]),
      unname(reconstructed_design),
      tolerance = 0,
      check.attributes = FALSE
    )),
  "Reconstructed pooled design differs from the validated primary model"
)

primary_all <- fread(PRIMARY_PATH)
assert_true(
  nrow(primary_all) == 23370L && !anyDuplicated(primary_all$gene),
  "Validated primary result universe drift"
)
assert_true(
  identical(primary_all$gene, rownames(dge_five)),
  "Validated primary/DGE order drift"
)
# TREAT check retained: it validates that the SEALED input table is unchanged.
assert_true(
  sum(primary_all$treat_fdr < LNCRNA_FDR) == 1616L &&
    all(primary_all$treat_lfc == LNCRNA_TREAT_LFC),
  "Validated primary TREAT family drift"
)
# CANONICAL 2026-08-12 family on the same table.
assert_true(
  sum(lncrna_canonical_positive(primary_all)) == 1347L,
  "Validated primary canonical family drift"
)

identity_tested <- identity[match(rownames(dge_five), gene_id_versioned)]
assert_true(
  nrow(identity_tested) == 23370L && !anyNA(identity_tested$gene_id_versioned) &&
    identical(identity_tested$gene_id_versioned, rownames(dge_five)),
  "Versioned-ID identity join is incomplete"
)
lncrna_ids <- identity_tested[gene_type == "lncRNA", gene_id_versioned]
assert_true(length(lncrna_ids) == 6249L, "Tested lncRNA census drift")
class_tested <- classes[match(lncrna_ids, gene_id_versioned)]
assert_true(
  nrow(class_tested) == length(lncrna_ids) &&
    !anyNA(class_tested$gene_id_versioned) &&
    identical(class_tested$gene_id_versioned, lncrna_ids),
  "Versioned-ID lncRNA class join is incomplete"
)

cohort_effects <- rbindlist(lapply(LNCRNA_COHORTS, function(cohort) {
  selected <- sample_info$dataset == cohort
  fit <- fit_disease_contrast(
    dge_five[, selected, keep.lib.sizes = TRUE],
    sample_info[selected],
    include_dataset = FALSE,
    quality_weights = TRUE
  )
  fit[gene_id_versioned %in% lncrna_ids, `:=`(
    cohort = cohort,
    n_samples = sum(selected),
    n_control = sum(sample_info[selected]$group_binary == "Control"),
    n_disease = sum(sample_info[selected]$group_binary == "Disease"),
    model = "voomWithQualityWeights:~inferred_sex+group_binary",
    bh_family_n = 23370L
  )]
  fit[gene_id_versioned %in% lncrna_ids]
}), use.names = TRUE)

loco_effects <- rbindlist(lapply(LNCRNA_COHORTS, function(excluded) {
  selected <- sample_info$dataset != excluded
  fit <- fit_disease_contrast(
    dge_five[, selected, keep.lib.sizes = TRUE],
    sample_info[selected],
    include_dataset = TRUE,
    quality_weights = TRUE
  )
  fit[gene_id_versioned %in% lncrna_ids, `:=`(
    excluded_cohort = excluded,
    n_samples = sum(selected),
    n_control = sum(sample_info[selected]$group_binary == "Control"),
    n_disease = sum(sample_info[selected]$group_binary == "Disease"),
    model = "voomWithQualityWeights:~dataset+inferred_sex+group_binary",
    bh_family_n = 23370L
  )]
  fit[gene_id_versioned %in% lncrna_ids]
}), use.names = TRUE)

equal_weight_all <- fit_disease_contrast(
  dge_five,
  sample_info,
  include_dataset = TRUE,
  quality_weights = FALSE
)
equal_weight <- equal_weight_all[gene_id_versioned %in% lncrna_ids]

primary <- primary_all[match(lncrna_ids, gene)]
primary[, gene_id_versioned := gene]
primary <- merge(
  primary,
  identity_tested[match(lncrna_ids, gene_id_versioned), .(
    gene_id_versioned, gene_id_base, gene_name, gene_type, chromosome,
    is_canonical_chromosome
  )],
  by = "gene_id_versioned",
  sort = FALSE
)
primary <- merge(
  primary,
  class_tested[, .(
    gene_id_versioned, lncrna_genomic_class, mapping_status,
    main_text_eligible, main_text_exclusion_reason
  )],
  by = "gene_id_versioned",
  sort = FALSE
)
assert_true(nrow(primary) == 6249L, "Primary lncRNA join drift")
assert_true(sum(primary$treat_fdr < LNCRNA_FDR) == 434L, "Primary lncRNA TREAT count drift")
assert_true(sum(lncrna_canonical_positive(primary)) == 402L, "Primary lncRNA canonical count drift")

bulk_results <- derive_high_confidence(
  primary,
  cohort_effects,
  loco_effects,
  equal_weight
)
setorder(bulk_results, gene_id_versioned)
setorder(cohort_effects, gene_id_versioned, cohort)
setorder(loco_effects, gene_id_versioned, excluded_cohort)

stage_results <- data.table(
  analysis_status = "blocked_no_fragment_native_stage_model",
  result_rows = 0L,
  biological_unit = "human_sample",
  multiple_testing_family = "not_run",
  detail = paste(
    "No fragment-native true-Kleiner stage model exists; no stage lncRNA",
    "effect, P value, FDR, or count is authorized."
  )
)
# Every count is computed live and carries its gate in the metric name. The
# hardcoded 1616/434 that used to sit here were TREAT values printed next to a
# live high_confidence count, so a single table mixed the two arms.
verdict <- data.table(
  metric = c(
    "source_gate", "tested_all_genes", "tested_lncrna",
    "primary_canonical_all_genes", "primary_canonical_lncrna",
    "high_confidence_lncrna_canonical",
    "primary_treat_all_genes", "primary_treat_lncrna",
    "complete_case_sensitivity", "stage_analysis"
  ),
  value = c(
    "pass", "23370", "6249",
    as.character(sum(lncrna_canonical_positive(primary_all))),
    as.character(sum(bulk_results$primary_canonical_positive)),
    as.character(sum(bulk_results$high_confidence)),
    as.character(sum(primary_all$treat_fdr < LNCRNA_FDR)),
    as.character(sum(bulk_results$primary_treat_positive)),
    "identical_to_primary_844_of_844_complete",
    "blocked_no_fragment_native_stage_model"
  ),
  interpretation = c(
    "Claim-bearing lncRNA analysis was permitted by the source gate.",
    "Frozen all-gene BH family.",
    "GENCODE v49 lncRNAs in the frozen tested universe.",
    "CANONICAL 2026-08-12 pooled all-gene family at padj<0.05 and |log2FC|>0.50.",
    "lncRNA subdivision of the same all-gene family; no lncRNA-only BH.",
    "Deterministic cohort, LOCO, mapping, and equal-library-weight rules, all on the canonical gate.",
    "COMPARATOR ARM: interval-null TREAT family at lfc=0.25 and BH FDR<0.05. Not the canonical gate.",
    "COMPARATOR ARM: lncRNA subdivision of the TREAT family. Do not report beside canonical counts without the label.",
    "All 844 primary-model samples have complete dataset, inferred-sex, and disease fields; a separate complete-case fit would be identical.",
    "No fragment-native stage inference has been run."
  )
)

if (!dir.create(OUTPUT, recursive = FALSE, mode = "0750")) {
  fail("Failed to create output directory: ", OUTPUT)
}
write_tsv_once(bulk_results, file.path(OUTPUT, "bulk_lncrna_results.tsv"))
write_tsv_once(cohort_effects, file.path(OUTPUT, "bulk_lncrna_cohort_effects.tsv"))
write_tsv_once(loco_effects, file.path(OUTPUT, "bulk_lncrna_leave_one_cohort_out.tsv"))
write_tsv_once(stage_results, file.path(OUTPUT, "stage_lncrna_results.tsv"))
write_tsv_once(verdict, file.path(OUTPUT, "bulk_lncrna_verdict.tsv"))

source_manifest <- data.table(
  role = c(
    "frozen_dge", "frozen_metadata", "validated_primary_all_gene_result",
    "validated_primary_marker", "validated_primary_design", "source_gate",
    "gene_identity", "lncrna_class", "producer", "analysis_library"
  ),
  path = required_inputs,
  size_bytes = file.info(required_inputs)$size,
  sha256 = vapply(required_inputs, sha256_file, character(1))
)
write_tsv_once(source_manifest, file.path(OUTPUT, "source_manifest.tsv"))
write_tsv_once(
  data.table(
    component = c("R", "data.table", "edgeR", "limma"),
    version = c(
      as.character(getRversion()),
      as.character(packageVersion("data.table")),
      as.character(packageVersion("edgeR")),
      as.character(packageVersion("limma"))
    ),
    seed = 20260811L,
    slurm_job_id = Sys.getenv("SLURM_JOB_ID", "not_slurm"),
    biological_unit = "human_sample",
    primary_model = "validated_voomWithQualityWeights:~dataset+inferred_sex+group_binary",
    primary_multiple_testing = "TREAT_|logFC|>0.25_BH_across_23370_genes"
  ),
  file.path(OUTPUT, "execution_manifest.tsv")
)
writeLines(capture.output(sessionInfo()), file.path(OUTPUT, "sessionInfo.txt"), useBytes = TRUE)

output_names <- c(
  "bulk_lncrna_results.tsv", "bulk_lncrna_cohort_effects.tsv",
  "bulk_lncrna_leave_one_cohort_out.tsv", "stage_lncrna_results.tsv",
  "bulk_lncrna_verdict.tsv", "source_manifest.tsv", "execution_manifest.tsv",
  "sessionInfo.txt"
)
output_manifest <- data.table(
  relative_path = output_names,
  size_bytes = file.info(file.path(OUTPUT, output_names))$size,
  sha256 = vapply(file.path(OUTPUT, output_names), sha256_file, character(1))
)
write_tsv_once(output_manifest, file.path(OUTPUT, "output_manifest.tsv"))
