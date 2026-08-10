#!/usr/bin/env Rscript

# Candidate-only v2 snATAC projection for SP-INT-05.
#
# This is a bounded copy of the validated dynamic promoter-accessibility logic
# in 03_atac_projection.R.  It never writes to the v1 result root.  Before
# emitting v2 results it re-runs the two corresponding v1 programs and requires
# their raw effects, exact permutation P values, sensitivities, coverage, and
# donor scores to reproduce the accepted v1 rows.  BH is then recomputed over
# the complete frozen two-program v2 family within each cohort.

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(edgeR)
  library(GenomicRanges)
  library(rtracklayer)
})

set.seed(42)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Cannot resolve candidate ATAC producer path")
SCRIPT_FILE <- normalizePath(sub("^--file=", "", script_arg))

BASE <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
RELEASE_ID <- "program-context-v2-candidate-2026-08-07"
CANDIDATE_ROOT <- Sys.getenv(
  "SPATIAL_CONTEXT_CANDIDATE_ROOT",
  file.path(
    BASE, "Analysis/Multimodal_Program_Projection/candidates",
    RELEASE_ID, "spatial_context"
  )
)
TARGET <- file.path(CANDIDATE_ROOT, "protein_atac_native/atac_v2")
if (dir.exists(TARGET)) stop("Refusing to overwrite candidate ATAC output: ", TARGET)
STAGE <- file.path(
  CANDIDATE_ROOT,
  paste0(".atac_v2.incomplete.", Sys.getpid())
)
if (dir.exists(STAGE)) stop("ATAC staging path already exists: ", STAGE)
dir.create(STAGE, recursive = TRUE, showWarnings = FALSE)

HOTSPOT <- file.path(
  BASE, "Analysis/Multimodal_Program_Projection/candidates",
  RELEASE_ID, "hotspot"
)
READY_FILE <- file.path(HOTSPOT, "READY")
REGISTRY_FILE <- file.path(HOTSPOT, "program_registry_v2.tsv")
MEMBERSHIP_FILE <- file.path(HOTSPOT, "program_membership_v2.tsv")
V1_ROOT <- file.path(BASE, "Analysis/Multimodal_Program_Projection")
V1_REGISTRY_FILE <- file.path(V1_ROOT, "results/frozen_programs.tsv")
V1_MEMBERSHIP_FILE <- file.path(V1_ROOT, "results/frozen_program_membership.tsv")
V1_EFFECT_FILE <- file.path(V1_ROOT, "results/atac/dynamic_program_accessibility.tsv")
V1_SCORE_FILE <- file.path(V1_ROOT, "results/atac/dynamic_program_scores.tsv")
GTF_FILE <- "/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"
G244 <- file.path(BASE, "Analysis/ATAC/Human_Multiome/results/snapatac2")
G244_META <- file.path(BASE, "Analysis/ATAC/Human_Multiome/metadata/donor_metadata_curated.tsv")
G281 <- file.path(BASE, "Analysis/ATAC/Human_External/pseudobulk")

required <- c(
  READY_FILE, REGISTRY_FILE, MEMBERSHIP_FILE, V1_REGISTRY_FILE,
  V1_MEMBERSHIP_FILE, V1_EFFECT_FILE, V1_SCORE_FILE, GTF_FILE, G244_META,
  file.path(G281, "hep_pseudobulk_counts_GSE281367.tsv.gz"),
  file.path(G281, "hep_pseudobulk_coldata_GSE281367.tsv"),
  file.path(G244, "hep_pseudobulk_counts.tsv.gz")
)
if (any(!file.exists(required))) {
  stop("Missing candidate ATAC input(s): ", paste(required[!file.exists(required)], collapse = ", "))
}

sha256 <- function(path) digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
ready <- fread(READY_FILE)
if (
  nrow(ready) != 1L || ready$release_id != RELEASE_ID ||
    ready$status != "ready_for_external_testing" ||
    ready$external_outcomes_read != FALSE
) stop("Plan 20 READY seal is not outcome-blind and ready for external testing")
if (sha256(REGISTRY_FILE) != ready$registry_sha256) stop("Plan 20 registry hash drift")
if (sha256(MEMBERSHIP_FILE) != ready$membership_table_sha256) stop("Plan 20 membership hash drift")

registry <- fread(REGISTRY_FILE)[external_test_eligible == TRUE]
if (nrow(registry) != 2L || uniqueN(registry$cell_type) != 1L || registry$cell_type[1] != "hepatocytes") {
  stop("Expected the sealed two-program hepatocyte external-test family")
}
setorder(registry, program_uid)
registry[, legacy_program_id := paste0(cell_type, "::", module)]

membership_raw <- fread(MEMBERSHIP_FILE)
membership <- membership_raw[
  program_uid %chin% registry$program_uid &
    !is.na(mapped_symbol) & mapped_symbol != "" &
    mapped_symbol_status != "unmapped_or_ambiguous",
  .(original_l1_weight = sum(original_l1_weight)),
  by = .(program_uid, membership_sha256, gene_symbol = mapped_symbol)
]
if (!setequal(unique(membership$program_uid), registry$program_uid)) {
  stop("One or more v2 programs has no mapped ATAC membership")
}

v1_registry <- fread(V1_REGISTRY_FILE)[program_id %chin% registry$legacy_program_id]
v1_membership <- fread(V1_MEMBERSHIP_FILE)[
  program_id %chin% registry$legacy_program_id & mapped_symbol == TRUE &
    !is.na(gene_symbol) & gene_symbol != "",
  .(original_l1_weight = sum(original_l1_weight)),
  by = .(program_id, gene_symbol)
]
if (nrow(v1_registry) != 2L || !setequal(v1_registry$program_id, registry$legacy_program_id)) {
  stop("The two v2 programs do not map uniquely to accepted v1 programs")
}

message("[ATAC-v2] importing GENCODE v49 gene coordinates")
genes <- import(GTF_FILE, format = "gtf", feature.type = "gene")
standard_chr <- paste0("chr", c(1:22, "X", "Y"))
genes <- genes[as.character(seqnames(genes)) %in% standard_chr]
names(genes) <- as.character(mcols(genes)$gene_name)
genes <- genes[!is.na(names(genes)) & names(genes) != ""]
genes <- genes[!duplicated(names(genes))]
promoter_gr <- promoters(genes, upstream = 2000L, downstream = 2001L)

fast_read_counts <- function(path) {
  con <- gzfile(path, "r")
  on.exit(close(con), add = TRUE)
  lines <- readLines(con)
  peaks <- strsplit(lines[1], "\t", fixed = TRUE)[[1]][-1]
  body <- lines[-1]
  mat <- matrix(0, nrow = length(body), ncol = length(peaks))
  donors <- character(length(body))
  for (i in seq_along(body)) {
    fields <- strsplit(body[i], "\t", fixed = TRUE)[[1]]
    donors[i] <- fields[1]
    mat[i, ] <- as.numeric(fields[-1])
  }
  rownames(mat) <- donors
  colnames(mat) <- peaks
  mat
}

parse_peaks <- function(x) {
  mt <- regexec("^(chr[^:]+):([0-9]+)-([0-9]+)$", x)
  bits <- regmatches(x, mt)
  good <- lengths(bits) == 4L
  out <- data.table(
    peak_index = which(good),
    peak_name = x[good],
    chr = vapply(bits[good], `[`, character(1), 2),
    start0 = as.integer(vapply(bits[good], `[`, character(1), 3)),
    end = as.integer(vapply(bits[good], `[`, character(1), 4))
  )
  out[, start := start0 + 1L]
  out
}

exact_two_group_p <- function(score, case) {
  score <- as.numeric(score)
  case <- as.logical(case)
  ok <- is.finite(score) & !is.na(case)
  score <- score[ok]
  case <- case[ok]
  n <- length(score)
  n_case <- sum(case)
  if (n < 6L || n_case == 0L || n_case == n) {
    return(c(effect = NA_real_, pvalue = NA_real_, std_error = NA_real_))
  }
  observed <- mean(score[case]) - mean(score[!case])
  cmb <- combn(n, n_case)
  total <- sum(score)
  case_means <- colSums(matrix(score[cmb], nrow = n_case)) / n_case
  perm <- case_means - (total - case_means * n_case) / (n - n_case)
  se <- sqrt(stats::var(score[case]) / n_case + stats::var(score[!case]) / (n - n_case))
  c(
    effect = observed,
    pvalue = mean(abs(perm) >= abs(observed) - 1e-12),
    std_error = se
  )
}

score_module <- function(z, m, samples, case, id_column, id_value) {
  measured <- m[get(id_column) == id_value & gene_symbol %in% rownames(z)]
  n_measured <- nrow(measured)
  retained <- sum(measured$original_l1_weight)
  testable <- n_measured >= 8L && retained >= 0.20
  if (!testable) {
    return(list(
      summary = data.table(
        n_measured = n_measured, retained_l1_weight = retained,
        testable = FALSE, effect = NA_real_, std_error = NA_real_,
        pvalue = NA_real_, equal_effect = NA_real_, leave_top_effect = NA_real_
      ),
      scores = NULL
    ))
  }
  measured[, weight := original_l1_weight / sum(original_l1_weight)]
  zz <- z[measured$gene_symbol, samples, drop = FALSE]
  observed_weight <- colSums(is.finite(zz) * measured$weight)
  weighted <- colSums(zz * measured$weight, na.rm = TRUE) / observed_weight
  weighted[!is.finite(weighted)] <- NA_real_
  equal <- colMeans(zz, na.rm = TRUE)
  equal[!is.finite(equal)] <- NA_real_
  top_gene <- measured$gene_symbol[which.max(measured$weight)]
  leave <- measured[gene_symbol != top_gene]
  leave[, weight := original_l1_weight / sum(original_l1_weight)]
  leave_zz <- z[leave$gene_symbol, samples, drop = FALSE]
  leave_observed_weight <- colSums(is.finite(leave_zz) * leave$weight)
  leave_score <- colSums(leave_zz * leave$weight, na.rm = TRUE) / leave_observed_weight
  leave_score[!is.finite(leave_score)] <- NA_real_
  primary <- exact_two_group_p(weighted, case)
  eq <- exact_two_group_p(equal, case)
  lo <- exact_two_group_p(leave_score, case)
  list(
    summary = data.table(
      n_measured = n_measured,
      retained_l1_weight = retained,
      testable = TRUE,
      effect = unname(primary["effect"]),
      std_error = unname(primary["std_error"]),
      pvalue = unname(primary["pvalue"]),
      equal_effect = unname(eq["effect"]),
      leave_top_effect = unname(lo["effect"])
    ),
    scores = data.table(
      sample_id = samples,
      score = weighted,
      equal_score = equal,
      leave_top_score = leave_score
    )
  )
}

numeric_equal <- function(a, b, tol = 1e-10) {
  if (length(a) != length(b)) return(FALSE)
  both_na <- is.na(a) & is.na(b)
  finite <- is.finite(a) & is.finite(b)
  all(both_na | (finite & abs(a - b) <= tol))
}

v2_results <- list()
v2_scores <- list()
v1_audit <- list()
v1_score_audit <- list()
coverage <- list()
v1_effects <- fread(V1_EFFECT_FILE)
v1_scores <- fread(V1_SCORE_FILE)

for (cohort in c("GSE281367", "GSE244832")) {
  current_cohort <- cohort
  message("[ATAC-v2] donor promoter accessibility: ", cohort)
  if (cohort == "GSE281367") {
    counts_file <- file.path(G281, "hep_pseudobulk_counts_GSE281367.tsv.gz")
    coldata <- fread(file.path(G281, "hep_pseudobulk_coldata_GSE281367.tsv"))
    coldata[, condition_curated := toupper(condition)]
  } else {
    counts_file <- file.path(G244, "hep_pseudobulk_counts.tsv.gz")
    coldata <- fread(G244_META, select = c("donor_id", "condition"))
    coldata[, condition_curated := toupper(condition)]
  }
  coldata <- coldata[condition_curated %in% c("NORMAL", "MASH")]
  counts <- fast_read_counts(counts_file)
  counts <- counts[, !duplicated(colnames(counts)), drop = FALSE]
  samples <- intersect(coldata$donor_id, rownames(counts))
  coldata <- coldata[match(samples, donor_id)]
  counts <- counts[samples, , drop = FALSE]
  expected_groups <- if (cohort == "GSE281367") c(6L, 6L) else c(5L, 9L)
  observed_groups <- c(
    sum(coldata$condition_curated == "NORMAL"),
    sum(coldata$condition_curated == "MASH")
  )
  if (!identical(observed_groups, expected_groups)) {
    stop(cohort, " donor contrast drift: ", paste(observed_groups, collapse = "/"))
  }

  peaks <- parse_peaks(colnames(counts))
  counts <- counts[, peaks$peak_index, drop = FALSE]
  peak_gr <- GRanges(
    seqnames = peaks$chr,
    ranges = IRanges(start = peaks$start, end = peaks$end)
  )
  all_program_genes <- unique(c(v1_membership$gene_symbol, membership$gene_symbol))
  all_program_genes <- intersect(all_program_genes, names(promoter_gr))
  prom <- promoter_gr[all_program_genes]
  hits <- findOverlaps(prom, peak_gr, ignore.strand = TRUE)
  by_gene <- split(subjectHits(hits), queryHits(hits))
  gene_counts <- matrix(
    0,
    nrow = length(all_program_genes),
    ncol = length(samples),
    dimnames = list(all_program_genes, samples)
  )
  for (idx in names(by_gene)) {
    peak_idx <- unique(by_gene[[idx]])
    gene_counts[as.integer(idx), ] <- rowSums(counts[, peak_idx, drop = FALSE])
  }
  measured_gene <- rowSums(gene_counts) >= 10 & rowSums(gene_counts > 0) >= 3
  gene_counts <- gene_counts[measured_gene, , drop = FALSE]
  dge <- DGEList(counts = t(counts))
  dge <- calcNormFactors(dge, method = "TMM")
  eff_lib <- dge$samples$lib.size * dge$samples$norm.factors
  logcpm <- cpm(gene_counts, lib.size = eff_lib, log = TRUE, prior.count = 0.5)
  z <- t(scale(t(logcpm)))
  z[!is.finite(z)] <- NA_real_
  case <- coldata$condition_curated == "MASH"

  for (i in seq_len(nrow(registry))) {
    uid <- registry$program_uid[i]
    legacy <- registry$legacy_program_id[i]
    old <- score_module(z, v1_membership, samples, case, "program_id", legacy)
    accepted <- v1_effects[cohort == current_cohort & program_id == legacy]
    if (nrow(accepted) != 1L) stop("Accepted v1 ATAC row is not unique: ", cohort, " / ", legacy)
    raw_fields <- c("retained_l1_weight", "effect", "pvalue", "equal_effect", "leave_top_effect")
    raw_match <- all(vapply(raw_fields, function(field) {
      numeric_equal(old$summary[[field]], accepted[[field]])
    }, logical(1)))
    discrete_match <- old$summary$n_measured == accepted$n_measured &&
      identical(as.logical(old$summary$testable), as.logical(accepted$testable))
    if (!raw_match || !discrete_match) {
      stop("Candidate engine fails v1 numeric regression for ", cohort, " / ", legacy)
    }
    accepted_scores <- v1_scores[cohort == current_cohort & program_id == legacy]
    old_scores <- copy(old$scores)
    if (is.null(old_scores)) stop("Accepted v1 selected program unexpectedly untestable")
    old_scores[, `:=`(cohort = cohort, program_id = legacy)]
    score_join <- merge(
      old_scores[, .(cohort, program_id, sample_id, score_rederived = score)],
      accepted_scores[, .(cohort, program_id, sample_id, score_accepted = score)],
      by = c("cohort", "program_id", "sample_id"),
      all = TRUE
    )
    score_pass <- nrow(score_join) == length(samples) &&
      !anyNA(score_join$score_rederived) && !anyNA(score_join$score_accepted) &&
      max(abs(score_join$score_rederived - score_join$score_accepted)) <= 1e-10
    if (!score_pass) stop("Candidate engine fails v1 donor-score regression for ", cohort, " / ", legacy)
    v1_audit[[paste(cohort, legacy)]] <- data.table(
      cohort = cohort,
      legacy_program_id = legacy,
      n_measured = old$summary$n_measured,
      retained_l1_weight = old$summary$retained_l1_weight,
      effect_abs_diff = abs(old$summary$effect - accepted$effect),
      pvalue_abs_diff = abs(old$summary$pvalue - accepted$pvalue),
      equal_effect_abs_diff = abs(old$summary$equal_effect - accepted$equal_effect),
      leave_top_effect_abs_diff = abs(old$summary$leave_top_effect - accepted$leave_top_effect),
      max_score_abs_diff = max(abs(score_join$score_rederived - score_join$score_accepted)),
      regression_pass = TRUE
    )

    ans <- score_module(z, membership, samples, case, "program_uid", uid)
    ans$summary[, `:=`(
      release_id = RELEASE_ID,
      registry_sha256 = ready$registry_sha256,
      ready_sha256 = sha256(READY_FILE),
      program_uid = uid,
      membership_sha256 = registry$membership_sha256[i],
      legacy_program_id = legacy,
      program_label = registry$module_name[i],
      cell_type = registry$cell_type[i],
      cohort = cohort,
      contrast = "MASH_vs_NORMAL",
      n_normal = observed_groups[1],
      n_mash = observed_groups[2],
      direction_expected = registry$primary_direction[i]
    )]
    v2_results[[paste(cohort, uid)]] <- ans$summary
    if (!is.null(ans$scores)) {
      ans$scores[, `:=`(
        release_id = RELEASE_ID,
        registry_sha256 = ready$registry_sha256,
        program_uid = uid,
        membership_sha256 = registry$membership_sha256[i],
        cohort = cohort,
        condition = coldata$condition_curated
      )]
      v2_scores[[paste(cohort, uid)]] <- ans$scores
    }
  }
  coverage[[cohort]] <- data.table(
    cohort = cohort,
    cell_type = "hepatocytes",
    n_donors = length(samples),
    n_normal = observed_groups[1],
    n_mash = observed_groups[2],
    n_unique_peaks = ncol(counts),
    n_candidate_program_genes_with_promoter_activity = nrow(gene_counts)
  )
}

dynamic <- rbindlist(v2_results, fill = TRUE)
dynamic[testable == TRUE, padj := p.adjust(pvalue, method = "BH"), by = cohort]
dynamic[, sensitivity_sign_agree := testable & is.finite(effect) &
  sign(effect) == sign(equal_effect) & sign(effect) == sign(leave_top_effect)]
dynamic[, robust := testable & is.finite(padj) & padj < 0.05 & sensitivity_sign_agree]
dynamic[, `:=`(
  pvalue_method = "exact_complete_donor_label_permutation_two_sided",
  multiplicity_family = "two_external_test_eligible_programs_within_cohort"
)]
setcolorder(dynamic, c(
  "release_id", "registry_sha256", "ready_sha256", "program_uid",
  "membership_sha256", "legacy_program_id", "program_label", "cell_type",
  "cohort", "contrast", "n_measured", "retained_l1_weight", "testable",
  "effect", "std_error", "pvalue", "equal_effect", "leave_top_effect",
  "padj", "sensitivity_sign_agree", "robust", "direction_expected",
  "n_normal", "n_mash", "pvalue_method", "multiplicity_family"
))
setorder(dynamic, cohort, program_uid)
scores <- rbindlist(v2_scores, fill = TRUE)
setorder(scores, cohort, program_uid, sample_id)
audit <- rbindlist(v1_audit, fill = TRUE)
setorder(audit, cohort, legacy_program_id)
coverage_table <- rbindlist(coverage, fill = TRUE)
setorder(coverage_table, cohort)

fwrite(dynamic, file.path(STAGE, "dynamic_program_accessibility.tsv"), sep = "\t", quote = FALSE)
fwrite(scores, file.path(STAGE, "dynamic_program_scores.tsv"), sep = "\t", quote = FALSE)
fwrite(coverage_table, file.path(STAGE, "dynamic_coverage.tsv"), sep = "\t", quote = FALSE)
fwrite(audit, file.path(STAGE, "v1_selected_regression.tsv"), sep = "\t", quote = FALSE)

input_paths <- unique(c(
  READY_FILE, REGISTRY_FILE, MEMBERSHIP_FILE, V1_REGISTRY_FILE,
  V1_MEMBERSHIP_FILE, V1_EFFECT_FILE, V1_SCORE_FILE, GTF_FILE, G244_META,
  file.path(G281, "hep_pseudobulk_counts_GSE281367.tsv.gz"),
  file.path(G281, "hep_pseudobulk_coldata_GSE281367.tsv"),
  file.path(G244, "hep_pseudobulk_counts.tsv.gz"), SCRIPT_FILE
))
input_manifest <- rbindlist(lapply(input_paths, function(path) {
  info <- file.info(path)
  data.table(
    absolute_path = normalizePath(path),
    bytes = as.numeric(info$size),
    sha256 = sha256(path),
    input_role = fifelse(
      normalizePath(path) == SCRIPT_FILE,
      "bounded_candidate_producer",
      "candidate_projection_input"
    )
  )
}))
fwrite(input_manifest, file.path(STAGE, "input_manifest.tsv"), sep = "\t", quote = FALSE)

status <- data.table(
  release_id = RELEASE_ID,
  status = "passed_v1_selected_numeric_regression_and_v2_projection",
  registry_sha256 = ready$registry_sha256,
  membership_table_sha256 = ready$membership_table_sha256,
  ready_sha256 = sha256(READY_FILE),
  producer_sha256 = sha256(SCRIPT_FILE),
  input_manifest_sha256 = sha256(file.path(STAGE, "input_manifest.tsv")),
  n_programs = uniqueN(dynamic$program_uid),
  n_cohorts = uniqueN(dynamic$cohort),
  n_v1_regression_rows = nrow(audit),
  all_v1_regression_pass = all(audit$regression_pass),
  generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
)
fwrite(status, file.path(STAGE, "producer_status.tsv"), sep = "\t", quote = FALSE)

dir.create(dirname(TARGET), recursive = TRUE, showWarnings = FALSE)
if (!file.rename(STAGE, TARGET)) stop("Could not atomically promote candidate ATAC staging directory")
cat("[ATAC-v2] PASS: v1 selected rows reproduced; v2 two-program families written\n")
print(dynamic[, .(cohort, program_uid, n_measured, retained_l1_weight, effect, pvalue, padj, robust)])
