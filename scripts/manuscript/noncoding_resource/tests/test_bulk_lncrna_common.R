#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
})

TEST_DIR <- dirname(normalizePath(
  sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
))
source(file.path(dirname(TEST_DIR), "bulk_lncrna_common.R"), local = TRUE)

expect_error <- function(expression, pattern) {
  observed <- tryCatch(
    {
      force(expression)
      ""
    },
    error = function(error) conditionMessage(error)
  )
  if (!grepl(pattern, observed, fixed = TRUE)) {
    stop("Expected error containing '", pattern, "', observed: ", observed)
  }
}

test_source_gate <- function() {
  root <- tempfile("lncrna-gate-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  pass_path <- file.path(root, "pass.tsv")
  fail_path <- file.path(root, "fail.tsv")
  fwrite(
    data.table(
      gate = c("count_contract", "strandedness", "lncrna_source_gate"),
      status = c("pass", "pass", "pass"),
      detail = c("ok", "ok", "eligible")
    ),
    pass_path,
    sep = "\t"
  )
  fwrite(
    data.table(
      gate = c("count_contract", "strandedness", "lncrna_source_gate"),
      status = c("pass", "fail", "fail"),
      detail = c("ok", "failed", "stop")
    ),
    fail_path,
    sep = "\t"
  )
  read_source_gate(pass_path)
  expect_error(read_source_gate(fail_path), "lncrna_source_gate is not pass")
}

test_producer_stops_before_outcomes <- function() {
  root <- tempfile("lncrna-producer-gate-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  gate_path <- file.path(root, "source_gate_status.tsv")
  identity_path <- file.path(root, "gene_identity.tsv")
  class_path <- file.path(root, "lncrna_genomic_class.tsv")
  output_path <- file.path(root, "forbidden-output")
  fwrite(
    data.table(
      gate = c("count_contract", "strandedness", "lncrna_source_gate"),
      status = c("pass", "fail", "fail"),
      detail = c("ok", "failed", "stop")
    ),
    gate_path,
    sep = "\t"
  )
  writeLines("not read because the gate failed", identity_path)
  writeLines("not read because the gate failed", class_path)
  producer <- file.path(dirname(TEST_DIR), "build_bulk_lncrna_candidate.R")
  output <- suppressWarnings(
    system2(
      file.path(R.home("bin"), "Rscript"),
      c(
        "--vanilla", producer,
        paste0("--source-gate=", gate_path),
        paste0("--gene-identity=", identity_path),
        paste0("--lncrna-class=", class_path),
        paste0("--output=", output_path)
      ),
      stdout = TRUE,
      stderr = TRUE
    )
  )
  status <- attr(output, "status")
  if (is.null(status) || status == 0L || file.exists(output_path) ||
      !any(grepl("lncrna_source_gate is not pass", output, fixed = TRUE))) {
    stop("Producer fail-gate contract mismatch:\n", paste(output, collapse = "\n"))
  }
}

test_all_gene_fit <- function() {
  set.seed(17)
  n_genes <- 100L
  n_samples <- 40L
  samples <- CJ(
    replicate = seq_len(5L),
    dataset = c("A", "B"),
    inferred_sex = c("F", "M"),
    group_binary = c("Control", "Disease")
  )
  samples[, sample_id := sprintf("S%03d", seq_len(.N))]
  setcolorder(
    samples,
    c("sample_id", "dataset", "inferred_sex", "group_binary", "replicate")
  )
  counts <- matrix(
    rnbinom(n_genes * n_samples, mu = 40, size = 8),
    nrow = n_genes,
    dimnames = list(
      sprintf("ENSG%011d.1", seq_len(n_genes)),
      samples$sample_id
    )
  )
  group <- samples$group_binary
  counts[1:5, group == "Disease"] <- counts[1:5, group == "Disease"] + 80L
  dge <- calcNormFactors(DGEList(counts = counts))
  samples[, replicate := NULL]
  quality <- fit_disease_contrast(
    dge, samples, include_dataset = TRUE, quality_weights = TRUE
  )
  equal <- fit_disease_contrast(
    dge, samples, include_dataset = TRUE, quality_weights = FALSE
  )
  stopifnot(
    nrow(quality) == n_genes,
    nrow(equal) == n_genes,
    identical(quality$gene_id_versioned, rownames(counts)),
    all(quality$treat_lfc == 0.25),
    all(quality$treat_fdr >= 0 & quality$treat_fdr <= 1),
    all(equal$treat_fdr >= 0 & equal$treat_fdr <= 1)
  )
}

test_high_confidence_rules <- function() {
  genes <- c("ENSG00000000001.1", "ENSG00000000002.1")
  primary <- data.table(
    gene_id_versioned = genes,
    logFC = c(0.7, 0.7),
    treat_fdr = c(0.001, 0.001),
    is_canonical_chromosome = c("true", "true"),
    mapping_status = c("mapping_unambiguous", "mapping_unambiguous"),
    main_text_eligible = c("true", "true")
  )
  cohort <- CJ(gene_id_versioned = genes, cohort = paste0("C", 1:5))
  cohort[, `:=`(logFC = 0.5, CI_low = 0.2, CI_high = 0.8)]
  cohort[gene_id_versioned == genes[[1L]] & cohort == "C5",
         `:=`(logFC = -0.05, CI_low = -0.3, CI_high = 0.2)]
  cohort[gene_id_versioned == genes[[2L]] & cohort == "C5",
         `:=`(logFC = -0.5, CI_low = -0.8, CI_high = -0.2)]
  loco <- CJ(gene_id_versioned = genes, excluded_cohort = paste0("C", 1:5))
  loco[, `:=`(logFC = 0.5, treat_fdr = 0.01)]
  loco[excluded_cohort == "C5", treat_fdr := 0.2]
  equal <- data.table(
    gene_id_versioned = genes,
    logFC = c(0.6, 0.6),
    treat_p = c(0.001, 0.001),
    treat_fdr = c(0.002, 0.002)
  )
  result <- derive_high_confidence(primary, cohort, loco, equal)
  stopifnot(
    result[gene_id_versioned == genes[[1L]], high_confidence],
    !result[gene_id_versioned == genes[[2L]], high_confidence],
    result[gene_id_versioned == genes[[1L]], cohort_direction_agreement_n] == 4L,
    result[gene_id_versioned == genes[[1L]], cohort_materially_opposite_n] == 0L,
    result[gene_id_versioned == genes[[2L]], cohort_materially_opposite_n] == 1L
  )
}

test_source_gate()
test_producer_stops_before_outcomes()
test_all_gene_fit()
test_high_confidence_rules()
cat("PASS: source gate, all-gene model, and high-confidence rules\n")
