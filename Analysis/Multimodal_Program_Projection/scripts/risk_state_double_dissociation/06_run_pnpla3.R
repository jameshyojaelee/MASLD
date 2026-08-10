#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
})

candidate_id <- "risk-state-double-dissociation-2026-08-09"
all_args <- commandArgs(trailingOnly = FALSE)
script_path <- normalizePath(sub("^--file=", "", all_args[grep("^--file=", all_args)][1]), mustWork = TRUE)
project_root <- normalizePath(file.path(dirname(script_path), "../../../.."), mustWork = TRUE)
candidate_root <- file.path(project_root, "Analysis/Multimodal_Program_Projection/candidates", candidate_id)

read_tsv <- function(path) read.delim(path, sep = "\t", quote = "", check.names = FALSE, stringsAsFactors = FALSE)
write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path))
  write.table(x, temporary, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  if (!file.rename(temporary, path)) stop("Atomic rename failed: ", path)
}
bind_rows_fill <- function(frames) {
  frames <- frames[vapply(frames, function(x) !is.null(x) && nrow(x) > 0L, logical(1))]
  if (!length(frames)) return(data.frame())
  fields <- unique(unlist(lapply(frames, names)))
  frames <- lapply(frames, function(frame) {
    for (field in setdiff(fields, names(frame))) frame[[field]] <- NA
    frame[, fields, drop = FALSE]
  })
  do.call(rbind, frames)
}

if (!file.exists(file.path(candidate_root, "SEALED.json")) || !file.exists(file.path(candidate_root, "SEAL_VALIDATED.json"))) stop("Validated seal required")
gate <- read_tsv(file.path(candidate_root, "source_gate_status.tsv"))
if (!any(gate$source_id == "GSE313544" & gate$status == "pass")) stop("GSE313544 source gate failed")
if (!any(gate$source_id == "GSE158182" & gate$status == "pass_source_dependent")) stop("GSE158182 source gate failed")

classes <- read_tsv(file.path(candidate_root, "frozen_inputs/evidence_classes__frozen_evidence_classes.tsv"))
membership <- read_tsv(file.path(candidate_root, "frozen_inputs/program_membership__program_membership_v2.tsv"))
registry <- read_tsv(file.path(candidate_root, "frozen_inputs/program_registry__program_registry_v2.tsv"))
external <- read_tsv(file.path(candidate_root, "frozen_inputs/external_programs__external_test_programs.tsv"))
sample_manifest <- read_tsv(file.path(candidate_root, "sample_manifest.tsv"))
if (nrow(registry) != 117L || nrow(external) != 2L) stop("Frozen program family drift")

gencode <- read.delim(gzfile(file.path(project_root, "data/gencode_v49_gene_metadata.tsv.gz")), sep = "\t", quote = "", stringsAsFactors = FALSE)
gencode <- gencode[!duplicated(gencode$ensembl_base), , drop = FALSE]

collapse_symbols <- function(expression, symbols) {
  keep <- !is.na(symbols) & nzchar(symbols)
  expression <- expression[keep, , drop = FALSE]
  symbols <- symbols[keep]
  medians <- apply(expression, 1L, median, na.rm = TRUE)
  ordering <- order(symbols, -medians, seq_along(symbols))
  selected <- ordering[!duplicated(symbols[ordering])]
  output <- expression[selected, , drop = FALSE]
  rownames(output) <- symbols[selected]
  output
}

load_gse313544 <- function() {
  path <- file.path(candidate_root, "source/GSE313544_COUNTS/GSE313544_normalized_counts.txt.gz")
  raw <- read.delim(gzfile(path), sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  expression <- as.matrix(raw[, -1, drop = FALSE])
  storage.mode(expression) <- "double"
  rownames(expression) <- raw[[1]]
  expression <- log2(expression + 0.5)
  metadata <- sample_manifest[sample_manifest$source_id == "GSE313544", , drop = FALSE]
  metadata$sample_id <- gsub("-", "_", metadata$title, fixed = TRUE)
  metadata <- metadata[match(colnames(expression), metadata$sample_id), , drop = FALSE]
  if (any(is.na(metadata$gsm))) stop("GSE313544 count columns do not map one-to-one to source samples")
  metadata$risk_allele_dosage <- as.numeric(metadata$risk_allele_dosage)
  metadata$treatment_num <- as.numeric(metadata$treatment == "palmitate")
  metadata$risk_treatment <- metadata$risk_allele_dosage * metadata$treatment_num
  design <- model.matrix(~ 0 + factor(clone) + treatment_num + risk_treatment, metadata)
  if (qr(design)$rank != ncol(design)) stop("GSE313544 design is not full rank")
  list(expression = collapse_symbols(expression, rownames(expression)), metadata = metadata, design = design, coefficient = "risk_treatment")
}

load_gse158182 <- function() {
  path <- file.path(candidate_root, "source/GSE158182_COUNTS/GSE158182_RawCounts_Subread_Genes.txt.gz")
  raw <- read.delim(gzfile(path), sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  counts <- as.matrix(raw[, -1, drop = FALSE])
  storage.mode(counts) <- "double"
  ensembl <- sub("\\..*$", "", raw[[1]])
  metadata <- sample_manifest[sample_manifest$source_id == "GSE158182", , drop = FALSE]
  metadata$sample_id <- metadata$title
  metadata <- metadata[match(colnames(counts), metadata$sample_id), , drop = FALSE]
  if (any(is.na(metadata$gsm))) stop("GSE158182 count columns do not map one-to-one to source samples")
  keep_sample <- metadata$treatment %in% c("control", "palmitate")
  metadata <- metadata[keep_sample, , drop = FALSE]
  counts <- counts[, keep_sample, drop = FALSE]
  group <- interaction(metadata$genotype, metadata$treatment, drop = TRUE)
  dge <- DGEList(counts)
  keep_gene <- filterByExpr(dge, group = group)
  dge <- calcNormFactors(dge[keep_gene, , keep.lib.sizes = FALSE])
  expression <- cpm(dge, log = TRUE, prior.count = 0.5)
  symbols <- gencode$gene_name[match(ensembl[keep_gene], gencode$ensembl_base)]
  expression <- collapse_symbols(expression, symbols)
  metadata$treatment_num <- as.numeric(metadata$treatment == "palmitate")
  metadata$i148m_treatment <- as.numeric(metadata$genotype == "I148M") * metadata$treatment_num
  metadata$ko_treatment <- as.numeric(metadata$genotype == "KO") * metadata$treatment_num
  design <- model.matrix(~ 0 + factor(line) + factor(differentiation) + treatment_num + i148m_treatment + ko_treatment, metadata)
  if (qr(design)$rank != ncol(design)) stop("GSE158182 design is not full rank")
  list(expression = expression, metadata = metadata, design = design, coefficient = "i148m_treatment")
}

gene_z <- function(expression) {
  keep <- apply(expression, 1L, function(x) all(is.finite(x)) && sd(x) > 0)
  t(scale(t(expression[keep, , drop = FALSE])))
}

state_alignment_scores <- function(expression) {
  z <- gene_z(expression)
  state <- classes[classes$static_class == "disease_state_only" & is.finite(classes$bulk_logFC) & classes$bulk_logFC != 0, , drop = FALSE]
  state <- state[!duplicated(state$gene_symbol), , drop = FALSE]
  common <- intersect(rownames(z), state$gene_symbol)
  state <- state[match(common, state$gene_symbol), , drop = FALSE]
  aligned <- z[common, , drop = FALSE] * sign(state$bulk_logFC)
  equal <- colMeans(aligned)
  weights <- abs(state$bulk_t)
  weights[!is.finite(weights)] <- 0
  weights <- weights / sum(weights)
  weighted <- as.numeric(crossprod(weights, aligned))
  names(weighted) <- colnames(aligned)
  list(equal = equal, bulk_t_weighted = weighted, n_mapped = length(common), mapped_weight = sum(abs(state$bulk_t), na.rm = TRUE))
}

score_programs <- function(expression, scheme = c("weighted", "equal", "leave_top")) {
  scheme <- match.arg(scheme)
  z <- gene_z(expression)
  uids <- registry$program_uid
  scores <- matrix(NA_real_, nrow = length(uids), ncol = ncol(z), dimnames = list(uids, colnames(z)))
  audit <- vector("list", length(uids))
  for (i in seq_along(uids)) {
    members <- membership[membership$program_uid == uids[i], , drop = FALSE]
    members <- members[!duplicated(members$mapped_symbol), , drop = FALSE]
    retained <- members[members$mapped_symbol %in% rownames(z), , drop = FALSE]
    retained_weight <- sum(retained$original_l1_weight)
    testable <- nrow(retained) >= 8L && retained_weight >= 0.20
    if (testable) {
      if (scheme == "leave_top") retained <- retained[-which.max(retained$original_l1_weight), , drop = FALSE]
      weights <- if (scheme == "equal") rep(1, nrow(retained)) else retained$original_l1_weight
      weights <- weights / sum(weights)
      scores[i, ] <- as.numeric(crossprod(weights, z[retained$mapped_symbol, , drop = FALSE]))
    }
    audit[[i]] <- data.frame(
      program_uid = uids[i], scoring_scheme = scheme, n_mapped_genes = nrow(retained),
      retained_l1_weight = retained_weight, testable = testable,
      failure_reason = if (testable) "" else if (nrow(retained) < 8L) "fewer_than_8_mapped_genes" else "less_than_20pct_original_l1_weight",
      stringsAsFactors = FALSE
    )
  }
  list(scores = scores, audit = do.call(rbind, audit))
}

fit_scalar <- function(score, design, coefficient) {
  fit <- eBayes(lmFit(matrix(score, nrow = 1L, dimnames = list("state_alignment", names(score))), design))
  index <- match(coefficient, colnames(design))
  if (is.na(index)) stop("Coefficient absent: ", coefficient)
  data.frame(
    estimate = fit$coefficients[1, index], se = fit$stdev.unscaled[1, index] * sqrt(fit$s2.post[1]),
    statistic = fit$t[1, index], p = fit$p.value[1, index], stringsAsFactors = FALSE
  )
}

fit_program_family <- function(scored, design, coefficient, dataset_id, scheme) {
  testable <- rowSums(is.finite(scored$scores)) == ncol(scored$scores)
  output <- registry[, c("program_uid", "cell_type", "module", "module_name", "primary_beta", "primary_qvalue"), drop = FALSE]
  output$estimate <- output$se <- output$statistic <- output$p <- NA_real_
  if (any(testable)) {
    fit <- eBayes(lmFit(scored$scores[testable, , drop = FALSE], design))
    index <- match(coefficient, colnames(design))
    output$estimate[testable] <- fit$coefficients[, index]
    output$se[testable] <- fit$stdev.unscaled[, index] * sqrt(fit$s2.post)
    output$statistic[testable] <- fit$t[, index]
    output$p[testable] <- fit$p.value[, index]
  }
  output$landscape_117_q <- p.adjust(output$p, method = "BH", n = 117L)
  primary <- output$program_uid %in% external$program_uid
  output$external_two_program_q <- NA_real_
  output$external_two_program_q[primary] <- p.adjust(output$p[primary], method = "BH", n = 2L)
  output$dataset_id <- dataset_id
  output$object_type <- "frozen_program"
  output$scoring_scheme <- scheme
  output$coefficient <- coefficient
  output$expected_direction <- sign(output$primary_beta)
  output$n_samples <- ncol(scored$scores)
  output$release_id <- candidate_id
  output
}

sample_score_frame <- function(dataset_id, metadata, scores) {
  frames <- list()
  for (scheme in names(scores)) {
    frame <- data.frame(
      dataset_id = dataset_id, sample_id = metadata$sample_id, gsm = metadata$gsm,
      clone = metadata$clone, line = metadata$line, genotype = metadata$genotype,
      differentiation = metadata$differentiation, treatment = metadata$treatment,
      score_type = "disease_state_alignment", scoring_scheme = scheme,
      score = as.numeric(scores[[scheme]]), release_id = candidate_id,
      stringsAsFactors = FALSE
    )
    frames[[scheme]] <- frame
  }
  do.call(rbind, frames)
}

objects <- list(GSE313544 = load_gse313544(), GSE158182 = load_gse158182())
primary_effects <- list()
program_effects <- list()
program_testability <- list()
sample_scores <- list()

for (dataset_id in names(objects)) {
  object <- objects[[dataset_id]]
  state_scores <- state_alignment_scores(object$expression)
  sample_scores[[dataset_id]] <- sample_score_frame(dataset_id, object$metadata, state_scores[c("equal", "bulk_t_weighted")])
  for (scheme in c("equal", "bulk_t_weighted")) {
    fit <- fit_scalar(state_scores[[scheme]], object$design, object$coefficient)
    fit$dataset_id <- dataset_id
    fit$object_type <- "disease_state_alignment"
    fit$scoring_scheme <- scheme
    fit$coefficient <- object$coefficient
    fit$n_mapped_genes <- state_scores$n_mapped
    fit$n_samples <- ncol(object$expression)
    fit$q <- fit$p
    fit$release_id <- candidate_id
    primary_effects[[paste(dataset_id, scheme, sep = ":")]] <- fit
  }
  for (scheme in c("weighted", "equal", "leave_top")) {
    scored <- score_programs(object$expression, scheme)
    program_effects[[paste(dataset_id, scheme, sep = ":")]] <- fit_program_family(scored, object$design, object$coefficient, dataset_id, scheme)
    audit <- scored$audit
    audit$dataset_id <- dataset_id
    audit$release_id <- candidate_id
    program_testability[[paste(dataset_id, scheme, sep = ":")]] <- audit
  }
}

primary_effects <- do.call(rbind, primary_effects)
rownames(primary_effects) <- NULL
program_effects <- do.call(rbind, program_effects)
program_testability <- do.call(rbind, program_testability)
sample_scores <- do.call(rbind, sample_scores)

g313 <- objects$GSE313544
equal313 <- state_alignment_scores(g313$expression)$equal
clone_values <- unique(g313$metadata$clone)
clone_delta <- do.call(rbind, lapply(clone_values, function(clone) {
  index <- g313$metadata$clone == clone
  metadata <- g313$metadata[index, , drop = FALSE]
  scores <- equal313[index]
  data.frame(
    clone = clone, genotype = unique(metadata$genotype), risk_allele_dosage = unique(metadata$risk_allele_dosage),
    palmitate_minus_control = scores[metadata$treatment == "palmitate"] - scores[metadata$treatment == "control"],
    stringsAsFactors = FALSE
  )
}))
genotype_medians <- aggregate(palmitate_minus_control ~ risk_allele_dosage + genotype, clone_delta, median)
genotype_medians <- genotype_medians[order(genotype_medians$risk_allele_dosage), , drop = FALSE]
monotonic <- all(diff(genotype_medians$palmitate_minus_control) >= 0)

# Exact clone-label permutation, preserving the observed 2/3/3 CC/GC/GG clone
# allocation.  The response is one paired palmitate-minus-control delta per
# clone, so cultures/samples are not treated as independent clone labels.
observed_clone_slope <- unname(coef(lm(palmitate_minus_control ~ risk_allele_dosage, clone_delta))[2])
clone_indices <- seq_len(nrow(clone_delta))
permuted_slopes <- numeric(0)
for (zero_index in combn(clone_indices, 2L, simplify = FALSE)) {
  remaining <- setdiff(clone_indices, zero_index)
  for (one_index in combn(remaining, 3L, simplify = FALSE)) {
    permuted_dosage <- rep(2, length(clone_indices))
    permuted_dosage[zero_index] <- 0
    permuted_dosage[one_index] <- 1
    permuted_slopes <- c(
      permuted_slopes,
      unname(coef(lm(clone_delta$palmitate_minus_control ~ permuted_dosage))[2])
    )
  }
}
clone_permutation <- data.frame(
  dataset_id = "GSE313544",
  statistic = "paired_clone_delta_dosage_slope",
  observed = observed_clone_slope,
  exact_permutation_p_greater = mean(permuted_slopes >= observed_clone_slope),
  null_mean = mean(permuted_slopes),
  null_sd = sd(permuted_slopes),
  n_exact_allocations = length(permuted_slopes),
  label_structure = "2_CC_3_GC_3_GG_clones",
  release_id = candidate_id,
  stringsAsFactors = FALSE
)

loo_clone <- do.call(rbind, lapply(clone_values, function(clone) {
  keep <- g313$metadata$clone != clone
  design <- model.matrix(~ 0 + factor(clone) + treatment_num + risk_treatment, g313$metadata[keep, , drop = FALSE])
  fit <- fit_scalar(equal313[keep], design, "risk_treatment")
  transform(fit, omitted_clone = clone)
}))
loo_stable <- all(loo_clone$estimate > 0)

g158 <- primary_effects[primary_effects$dataset_id == "GSE158182" & primary_effects$scoring_scheme == "equal", , drop = FALSE]
g313_primary <- primary_effects[primary_effects$dataset_id == "GSE313544" & primary_effects$scoring_scheme == "equal", , drop = FALSE]
stress_pass <- nrow(g313_primary) == 1L && g313_primary$q < 0.05 && g313_primary$estimate > 0 && monotonic && loo_stable && nrow(g158) == 1L && g158$estimate > 0

verdict <- data.frame(
  arm = "stress_contingent_PNPLA3_boundary",
  status = if (stress_pass) "pass_source_dependent_boundary" else "valid_null_or_nonconfirmatory",
  gse313544_estimate = g313_primary$estimate,
  gse313544_q = g313_primary$q,
  gse313544_exact_clone_permutation_p = clone_permutation$exact_permutation_p_greater,
  gse313544_clone_permutation_null_mean = clone_permutation$null_mean,
  monotonic_CC_GC_GG = monotonic,
  leave_one_clone_all_positive = loo_stable,
  gse158182_i148m_direction_agrees = nrow(g158) == 1L && g158$estimate > 0,
  gse158182_estimate = g158$estimate,
  promotion_gate_pass = stress_pass,
  claim_boundary = "coding-risk stress observability only; PNPLA3 is not in the frozen regulatory genetic class",
  stringsAsFactors = FALSE
)

source_reproduction <- data.frame(
  dataset_id = c("GSE313544", "GSE158182"),
  observed_samples = c(nrow(objects$GSE313544$metadata), nrow(objects$GSE158182$metadata)),
  expected_samples_in_primary_model = c(16L, 20L),
  biological_units = c(length(unique(objects$GSE313544$metadata$clone)), length(unique(objects$GSE158182$metadata$line))),
  design_rank = c(qr(objects$GSE313544$design)$rank, qr(objects$GSE158182$design)$rank),
  design_columns = c(ncol(objects$GSE313544$design), ncol(objects$GSE158182$design)),
  passed = c(nrow(objects$GSE313544$metadata) == 16L, nrow(objects$GSE158182$metadata) == 20L),
  stringsAsFactors = FALSE
)

write_tsv(primary_effects, file.path(candidate_root, "stress_response_effects.tsv"))
write_tsv(program_effects, file.path(candidate_root, "stress_response_program_effects.tsv"))
write_tsv(program_testability, file.path(candidate_root, "stress_response_program_testability.tsv"))
write_tsv(sample_scores, file.path(candidate_root, "stress_response_sample_scores.tsv"))
write_tsv(clone_delta, file.path(candidate_root, "stress_response_clone_deltas.tsv"))
write_tsv(genotype_medians, file.path(candidate_root, "stress_response_genotype_medians.tsv"))
write_tsv(clone_permutation, file.path(candidate_root, "stress_response_permutations.tsv"))
write_tsv(loo_clone, file.path(candidate_root, "stress_response_leave_one_clone.tsv"))
write_tsv(verdict, file.path(candidate_root, "stress_response_verdict.tsv"))
write_tsv(source_reproduction, file.path(candidate_root, "stress_response_source_reproduction.tsv"))
cat("PNPLA3_COMPLETE\t", verdict$status, "\n", sep = "")
