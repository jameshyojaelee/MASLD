#!/usr/bin/env Rscript

# Candidate-only exact-Kleiner adjacent-stage models on the corrected
# fragment-count substrate. Significance is defined over the complete tested
# gene family for each contrast; lncRNAs are annotated only after fitting.

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
})

options(digits = 17, scipen = 999)
set.seed(20260811)

fail <- function(...) stop(..., call. = FALSE)
is_symlink <- function(path) {
  target <- Sys.readlink(path)
  !is.na(target) && nzchar(target)
}
sha256_file <- function(path) {
  path <- normalizePath(path, mustWork = TRUE)
  value <- system2("sha256sum", path, stdout = TRUE, stderr = TRUE)
  status <- attr(value, "status")
  if ((!is.null(status) && status != 0L) || length(value) != 1L) {
    fail("sha256sum failed for ", path)
  }
  digest <- strsplit(value[[1L]], "[[:space:]]+")[[1L]][[1L]]
  if (!grepl("^[0-9a-f]{64}$", digest)) fail("Invalid SHA256 for ", path)
  digest
}
write_tsv_once <- function(value, path) {
  if (file.exists(path) || is_symlink(path)) fail("Refusing overwrite: ", path)
  connection <- file(path, open = "wx")
  on.exit(try(close(connection), silent = TRUE), add = TRUE)
  write.table(
    as.data.frame(value), connection, sep = "\t", row.names = FALSE,
    col.names = TRUE, quote = FALSE, na = "NA", eol = "\n"
  )
  close(connection)
}

parse_arguments <- function(values) {
  expected <- c("source-gate", "gene-identity", "lncrna-class", "stage-preflight", "output")
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

arguments <- parse_arguments(commandArgs(trailingOnly = TRUE))
project_root <- normalizePath(
  Sys.getenv(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
  ),
  mustWork = TRUE
)
producer_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(producer_argument) != 1L) fail("Cannot resolve producer path")
producer <- normalizePath(sub("^--file=", "", producer_argument), mustWork = TRUE)
common_library <- normalizePath(file.path(dirname(producer), "bulk_lncrna_common.R"), mustWork = TRUE)
source(common_library, local = TRUE)

source_gate <- normalizePath(arguments[["source-gate"]], mustWork = TRUE)
identity_path <- normalizePath(arguments[["gene-identity"]], mustWork = TRUE)
lncrna_class_path <- normalizePath(arguments[["lncrna-class"]], mustWork = TRUE)
preflight_root <- normalizePath(arguments[["stage-preflight"]], mustWork = TRUE)
output_parent <- normalizePath(dirname(arguments[["output"]]), mustWork = TRUE)
output <- file.path(output_parent, basename(arguments[["output"]]))
if (file.exists(output) || is_symlink(output)) fail("Refusing existing output: ", output)

# The source gate is evaluated before expression or stage outcomes are read.
read_source_gate(source_gate)

candidate_root <- file.path(
  project_root,
  "RNA-seq/results/manuscript_release/candidates/",
  "resource-f-five-coloc-v6-candidate-2026-08-10"
)
dge_path <- file.path(
  candidate_root,
  "inputs/BG001-DECISION/arms/F_legacy/results/integration/merged_dge.rds"
)
ready_path <- file.path(preflight_root, "STAGE_PREFLIGHT_READY.json")
sample_path <- file.path(preflight_root, "audits/transition_sample_manifest.tsv")
transition_path <- file.path(preflight_root, "audits/transition_census.tsv")
validation_path <- file.path(preflight_root, "validation/validation_checks.tsv")
required_inputs <- c(
  dge_path, ready_path, sample_path, transition_path, validation_path,
  source_gate, identity_path, lncrna_class_path, producer, common_library
)
for (path in required_inputs) {
  if (!file.exists(path) || is_symlink(path)) fail("Missing or symlinked input: ", path)
}

expected_hashes <- c(
  dge = "56d1cf97f791a81404ddc6c5ea4c5df9e6a38e7cf011b39d0d7c3279c87afa4d",
  ready = "e0b34731685bc7b139ca4db15fd23dca96b45334d104f1d988bf441e1bbd38fa",
  samples = "41b8e27f9d9b3b2c85bb6c803e2c685292f2ff6eed15679dfbc18cebd01bb72b",
  transitions = "9487a21d17e34587b8b6d6eacb402eb2a714159b3c3ffef94fbb6a8fa945a1c0",
  validation = "54fcdeddc362d8291186d8ca3d6489b25db3b85ec3c316d022eb4bf727d5df68"
)
observed_hashes <- vapply(
  c(dge_path, ready_path, sample_path, transition_path, validation_path),
  sha256_file,
  character(1)
)
if (!identical(unname(observed_hashes), unname(expected_hashes))) {
  fail("Frozen stage input or validated preflight identity drift")
}

validation <- fread(validation_path, colClasses = "character")
if (!identical(names(validation), c("check_id", "status")) ||
    nrow(validation) != 7L || any(validation$status != "PASS")) {
  fail("Stage preflight validation is not complete")
}

expected_transitions <- data.table(
  transition = c("F0_to_F1", "F1_to_F2", "F2_to_F3", "F3_to_F4"),
  contrast = c("F1_vs_F0", "F2_vs_F1", "F3_vs_F2", "F4_vs_F3"),
  low_stage = 0:3,
  high_stage = 1:4,
  n_low = c(126L, 187L, 174L, 107L),
  n_high = c(187L, 174L, 132L, 42L),
  n_samples = c(313L, 361L, 306L, 149L),
  n_cohorts = c(6L, 6L, 6L, 5L)
)
transition_manifest <- fread(transition_path)
if (!identical(transition_manifest, expected_transitions)) {
  fail("Validated transition census drift")
}
samples <- fread(sample_path)
required_sample_fields <- c(
  "sample_id", "analysis_unit_id", "dataset", "fibrosis_stage",
  "inferred_sex", "transition", "contrast", "arm"
)
if (!identical(names(samples), required_sample_fields) ||
    anyDuplicated(samples[, .(transition, sample_id)]) ||
    anyDuplicated(samples[, .(transition, analysis_unit_id)])) {
  fail("Stage sample manifest schema or biological-unit uniqueness drift")
}

dge_all <- readRDS(dge_path)
if (!inherits(dge_all, "DGEList") || !identical(dim(dge_all), c(24196L, 1257L)) ||
    anyDuplicated(rownames(dge_all)) || anyDuplicated(colnames(dge_all))) {
  fail("Corrected fragment-count DGE identity drift")
}
identity <- fread(identity_path, colClasses = "character")
lncrna_class <- fread(lncrna_class_path, colClasses = "character")
if (anyDuplicated(identity$gene_id_versioned) ||
    anyDuplicated(lncrna_class$gene_id_versioned) ||
    !all(rownames(dge_all) %in% identity$gene_id_versioned)) {
  fail("GENCODE identity join is not one-to-one and complete")
}
annotation <- identity[
  match(rownames(dge_all), gene_id_versioned),
  .(
    gene_id_versioned, gene_id_base, gene_name, gene_type,
    chromosome, is_canonical_chromosome
  )
]
annotation[, display_biotype := fifelse(
  gene_type == "protein_coding", "protein_coding",
  fifelse(gene_type == "lncRNA", "lncRNA", "other")
)]
annotation <- merge(
  annotation,
  lncrna_class[, .(
    gene_id_versioned, lncrna_genomic_class, mapping_status,
    main_text_eligible, main_text_exclusion_reason
  )],
  by = "gene_id_versioned", all.x = TRUE, sort = FALSE
)
annotation <- annotation[match(rownames(dge_all), gene_id_versioned)]

tmp <- file.path(
  output_parent,
  paste0(".", basename(output), ".tmp.", Sys.getenv("SLURM_JOB_ID", Sys.getpid()))
)
if (file.exists(tmp) || is_symlink(tmp)) fail("Refusing existing temporary output: ", tmp)
dir.create(tmp, recursive = FALSE)
if (!dir.exists(tmp)) fail("Failed to create temporary output")

result_list <- list()
sample_list <- list()
design_list <- list()

for (index in seq_len(nrow(expected_transitions))) {
  expected <- expected_transitions[index]
  selected <- samples[transition == expected$transition]
  if (nrow(selected) != expected$n_samples ||
      sum(selected$arm == "low") != expected$n_low ||
      sum(selected$arm == "high") != expected$n_high ||
      uniqueN(selected$dataset) != expected$n_cohorts ||
      any(!selected$fibrosis_stage %in% c(expected$low_stage, expected$high_stage))) {
    fail("Stage manifest census drift for ", expected$transition)
  }
  if (!all(selected$sample_id %in% colnames(dge_all))) {
    fail("Stage sample missing from corrected DGE for ", expected$transition)
  }
  dge <- dge_all[, match(selected$sample_id, colnames(dge_all))]
  if (!identical(colnames(dge), selected$sample_id)) fail("DGE/sample order drift")

  info <- data.frame(
    dataset = droplevels(factor(selected$dataset)),
    inferred_sex = droplevels(factor(selected$inferred_sex)),
    fib_group = factor(selected$arm, levels = c("low", "high")),
    row.names = selected$sample_id
  )
  design <- model.matrix(~ dataset + inferred_sex + fib_group, data = info)
  coefficient <- "fib_grouphigh"
  if (!coefficient %in% colnames(design) || qr(design)$rank != ncol(design)) {
    fail("Stage design is not estimable for ", expected$transition)
  }
  transformed <- voomWithQualityWeights(dge, design, plot = FALSE)
  fitted <- eBayes(lmFit(transformed, design))
  coefficient_index <- match(coefficient, colnames(design))
  standard_error <- fitted$stdev.unscaled[, coefficient_index] * sqrt(fitted$s2.post)
  critical <- qt(0.975, df = fitted$df.total)
  names(standard_error) <- rownames(fitted$coefficients)
  names(critical) <- rownames(fitted$coefficients)
  table <- topTable(
    fitted, coef = coefficient_index, number = Inf, sort.by = "none"
  )
  genes <- rownames(dge)
  result <- data.table(
    gene_id_versioned = genes,
    logFC = table[genes, "logFC"],
    SE = standard_error[genes],
    CI_low = table[genes, "logFC"] - critical[genes] * standard_error[genes],
    CI_high = table[genes, "logFC"] + critical[genes] * standard_error[genes],
    t = table[genes, "t"],
    P.Value = table[genes, "P.Value"],
    FDR = table[genes, "adj.P.Val"],
    AveExpr = table[genes, "AveExpr"],
    transition = expected$transition,
    contrast = expected$contrast,
    low_stage = expected$low_stage,
    high_stage = expected$high_stage,
    n_low = expected$n_low,
    n_high = expected$n_high,
    n_samples = expected$n_samples,
    n_cohorts = expected$n_cohorts
  )
  if (max(abs(p.adjust(result$P.Value, method = "BH") - result$FDR)) > 1e-12 ||
      any(!is.finite(unlist(result[, .(logFC, SE, CI_low, CI_high, t, P.Value, FDR, AveExpr)])))) {
    fail("Stage statistical output failed internal checks for ", expected$transition)
  }
  result <- merge(result, annotation, by = "gene_id_versioned", sort = FALSE)
  result <- result[match(genes, gene_id_versioned)]
  result_list[[expected$transition]] <- result
  selected[, `:=`(
    low_stage = expected$low_stage,
    high_stage = expected$high_stage
  )]
  sample_list[[expected$transition]] <- selected
  design_list[[expected$transition]] <- data.table(
    transition = expected$transition,
    contrast = expected$contrast,
    formula = "~ dataset + inferred_sex + fib_group",
    coefficient = coefficient,
    design_rank = qr(design)$rank,
    n_design_columns = ncol(design),
    design_columns = paste(colnames(design), collapse = ";"),
    n_samples = nrow(design),
    n_biological_units = uniqueN(selected$analysis_unit_id),
    n_cohorts = uniqueN(selected$dataset),
    tested_family = nrow(result),
    multiple_testing = "BH_across_complete_contrast_family",
    significance_rule = "FDR<0.05_no_effect_size_threshold"
  )
}

all_results <- rbindlist(result_list, use.names = TRUE)
if (nrow(all_results) != 4L * 24196L ||
    any(all_results[, .N, by = transition]$N != 24196L)) {
  fail("Complete stage result family drift")
}
lncrna_results <- all_results[gene_type == "lncRNA"]
counts <- all_results[, .(
  n_tested = .N,
  n_fdr_positive = sum(FDR < 0.05),
  n_up = sum(FDR < 0.05 & logFC > 0),
  n_down = sum(FDR < 0.05 & logFC < 0)
), by = .(transition, contrast, display_biotype)]
counts[, `:=`(
  transition_order = match(transition, expected_transitions$transition),
  biotype_order = match(display_biotype, c("protein_coding", "lncRNA", "other"))
)]
setorder(counts, transition_order, biotype_order)
counts[, c("transition_order", "biotype_order") := NULL]
verdict <- all_results[, .(
  n_genes_tested = .N,
  n_lncrna_tested = sum(gene_type == "lncRNA"),
  n_genes_fdr_positive = sum(FDR < 0.05),
  n_lncrna_fdr_positive = sum(gene_type == "lncRNA" & FDR < 0.05),
  n_lncrna_main_text_eligible_fdr_positive = sum(
    gene_type == "lncRNA" & FDR < 0.05 & main_text_eligible == "true",
    na.rm = TRUE
  ),
  claim = "cross-sectional_stage-associated_remodeling_not_longitudinal_progression"
), by = .(transition, contrast, n_samples, n_cohorts)]

write_tsv_once(all_results, file.path(tmp, "stage_all_gene_results.tsv"))
write_tsv_once(lncrna_results, file.path(tmp, "stage_lncrna_results.tsv"))
write_tsv_once(counts, file.path(tmp, "stage_counts_by_biotype.tsv"))
write_tsv_once(rbindlist(sample_list), file.path(tmp, "stage_sample_manifest.tsv"))
write_tsv_once(rbindlist(design_list), file.path(tmp, "stage_design_audit.tsv"))
write_tsv_once(verdict, file.path(tmp, "stage_lncrna_verdict.tsv"))

source_manifest <- data.table(
  source_id = c(
    "corrected_fragment_dge", "validated_stage_preflight", "transition_samples",
    "transition_census", "preflight_validation", "lncrna_source_gate",
    "gene_identity", "lncrna_class", "producer", "analysis_library"
  ),
  path = required_inputs,
  size_bytes = file.info(required_inputs)$size,
  sha256 = vapply(required_inputs, sha256_file, character(1))
)
write_tsv_once(source_manifest, file.path(tmp, "source_manifest.tsv"))
write_tsv_once(
  data.table(
    component = c("R", "data.table", "edgeR", "limma"),
    version = c(
      as.character(getRversion()), as.character(packageVersion("data.table")),
      as.character(packageVersion("edgeR")), as.character(packageVersion("limma"))
    ),
    seed = 20260811L,
    slurm_job_id = Sys.getenv("SLURM_JOB_ID", "not_slurm"),
    biological_unit = "one_cross-sectional_human_biopsy_per_participant",
    model = "voomWithQualityWeights:~dataset+inferred_sex+adjacent_stage_group",
    multiple_testing = "BH_separately_across_24196_genes_per_adjacent_contrast"
  ),
  file.path(tmp, "execution_manifest.tsv")
)
writeLines(capture.output(sessionInfo()), file.path(tmp, "sessionInfo.txt"), useBytes = TRUE)

artifacts <- setdiff(list.files(tmp), "output_manifest.tsv")
output_manifest <- data.table(
  relative_path = artifacts,
  size_bytes = file.info(file.path(tmp, artifacts))$size,
  sha256 = vapply(file.path(tmp, artifacts), sha256_file, character(1))
)
write_tsv_once(output_manifest, file.path(tmp, "output_manifest.tsv"))
if (!file.rename(tmp, output)) fail("Atomic candidate publication failed; temporary output retained: ", tmp)
cat("STAGE_LNCRNA_CANDIDATE_COMPLETE\t", output, "\n", sep = "")
