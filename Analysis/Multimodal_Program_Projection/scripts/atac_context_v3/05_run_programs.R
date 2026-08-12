#!/usr/bin/env Rscript

# Complete 117-program donor-level promoter-accessibility projection.

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(Matrix)
  library(GenomicRanges)
  library(rtracklayer)
  library(digest)
})

set.seed(42)

BASE <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
RELEASE_ID <- "atac-context-v3-candidate-2026-08-11-r1"
PROGRAM_RELEASE_ID <- "program-context-v2-candidate-2026-08-07"
CANDIDATE <- Sys.getenv(
  "ATAC_V3_CANDIDATE_ROOT",
  file.path(BASE, "Analysis/Multimodal_Program_Projection/candidates", RELEASE_ID)
)
EXPECTED <- normalizePath(
  file.path(BASE, "Analysis/Multimodal_Program_Projection/candidates", RELEASE_ID),
  mustWork = FALSE
)
if (!identical(normalizePath(CANDIDATE, mustWork = FALSE), EXPECTED)) stop("Unsafe candidate root")
OUT <- file.path(CANDIDATE, "programs")
if (dir.exists(OUT)) stop("Refusing to overwrite program output: ", OUT)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

HOTSPOT <- file.path(
  BASE, "Analysis/Multimodal_Program_Projection/candidates",
  PROGRAM_RELEASE_ID, "hotspot"
)
REGISTRY_FILE <- file.path(HOTSPOT, "program_registry_v2.tsv")
MEMBERSHIP_FILE <- file.path(HOTSPOT, "program_membership_v2.tsv")
READY_FILE <- file.path(HOTSPOT, "READY")
GTF <- "/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"
required <- c(REGISTRY_FILE, MEMBERSHIP_FILE, READY_FILE, GTF)
if (any(!file.exists(required))) stop("Missing program input(s)")

ready <- fread(READY_FILE)
if (nrow(ready) != 1L || ready$release_id != PROGRAM_RELEASE_ID ||
    ready$status != "ready_for_external_testing" || ready$external_outcomes_read != FALSE) {
  stop("Frozen program READY seal is invalid")
}
if (digest(REGISTRY_FILE, algo = "sha256", file = TRUE, serialize = FALSE) != ready$registry_sha256) {
  stop("Program registry hash drift")
}
if (digest(MEMBERSHIP_FILE, algo = "sha256", file = TRUE, serialize = FALSE) != ready$membership_table_sha256) {
  stop("Program membership hash drift")
}

registry <- fread(REGISTRY_FILE)
if (nrow(registry) != 117L || uniqueN(registry$program_uid) != 117L) stop("Expected 117 unique programs")
registry[, lineage := fcase(
  cell_type == "hepatocytes", "hepatocyte",
  cell_type == "fibroblasts", "stellate",
  cell_type == "macrophages", "macrophage",
  cell_type == "cholangiocytes", "cholangiocyte",
  cell_type == "tcells", "t_nk",
  default = NA_character_
)]
if (anyNA(registry$lineage)) stop("Unmapped program lineage")

membership_raw <- fread(MEMBERSHIP_FILE)
membership <- membership_raw[
  !is.na(mapped_symbol) & mapped_symbol != "" &
    mapped_symbol_status != "unmapped_or_ambiguous",
  .(original_l1_weight = sum(original_l1_weight)),
  by = .(program_uid, membership_sha256, gene_symbol = mapped_symbol)
]
if (!setequal(unique(membership$program_uid), registry$program_uid)) stop("Program membership incomplete")

read_sparse_gz <- function(path) {
  con <- gzfile(path, "rb")
  on.exit(close(con), add = TRUE)
  as(readMM(con), "dgCMatrix")
}

load_counts <- function(cohort, lineage) {
  root <- file.path(CANDIDATE, "counts", cohort)
  counts <- read_sparse_gz(file.path(root, paste0(lineage, ".mtx.gz")))
  donors <- fread(file.path(root, paste0(lineage, ".donors.tsv")))
  peaks <- fread(file.path(root, paste0(lineage, ".peaks.tsv")))
  if (nrow(counts) != nrow(donors) || ncol(counts) != nrow(peaks)) stop("Count dimension drift")
  rownames(counts) <- donors$donor_id
  colnames(counts) <- peaks$peak_coordinate
  list(counts = counts, donors = donors, peaks = peaks)
}

parse_peaks <- function(values) {
  match <- regexec("^(chr[^:]+):([0-9]+)-([0-9]+)$", values)
  pieces <- regmatches(values, match)
  if (any(lengths(pieces) != 4L)) stop("Malformed peak coordinate")
  GRanges(
    seqnames = vapply(pieces, `[`, character(1), 2L),
    ranges = IRanges(
      start = as.integer(vapply(pieces, `[`, character(1), 3L)) + 1L,
      end = as.integer(vapply(pieces, `[`, character(1), 4L))
    )
  )
}

message("[PROGRAM] Importing GENCODE v49 promoters")
genes <- import(GTF, format = "gtf", feature.type = "gene")
genes <- genes[as.character(seqnames(genes)) %chin% paste0("chr", c(1:22, "X", "Y"))]
names(genes) <- as.character(mcols(genes)$gene_name)
genes <- genes[!is.na(names(genes)) & names(genes) != ""]
promoters_v49 <- promoters(genes, upstream = 2000L, downstream = 2001L)

score_program <- function(z, program_membership) {
  measured <- program_membership[gene_symbol %chin% rownames(z)]
  n_measured <- nrow(measured)
  retained <- sum(measured$original_l1_weight)
  testable <- n_measured >= 8L && retained >= 0.20
  if (!testable) {
    return(list(
      n_measured = n_measured,
      retained = retained,
      testable = FALSE,
      weighted = NULL,
      equal = NULL,
      leave = NULL,
      top_gene = ""
    ))
  }
  measured[, weight := original_l1_weight / sum(original_l1_weight)]
  values <- z[measured$gene_symbol, , drop = FALSE]
  weighted <- colSums(values * measured$weight, na.rm = TRUE) /
    colSums(is.finite(values) * measured$weight)
  equal <- colMeans(values, na.rm = TRUE)
  top_gene <- measured$gene_symbol[which.max(measured$weight)]
  leave_membership <- measured[gene_symbol != top_gene]
  leave_membership[, weight := original_l1_weight / sum(original_l1_weight)]
  leave_values <- z[leave_membership$gene_symbol, , drop = FALSE]
  leave <- colSums(leave_values * leave_membership$weight, na.rm = TRUE) /
    colSums(is.finite(leave_values) * leave_membership$weight)
  list(
    n_measured = n_measured,
    retained = retained,
    testable = TRUE,
    weighted = weighted,
    equal = equal,
    leave = leave,
    top_gene = top_gene
  )
}

coverage_rows <- list()
score_rows <- list()
result_rows <- list()
model_qc_rows <- list()

for (cohort in c("GSE244832", "GSE281367")) {
  for (lineage in c("hepatocyte", "stellate", "macrophage", "cholangiocyte", "t_nk")) {
    message("[PROGRAM] ", cohort, " / ", lineage)
    input <- load_counts(cohort, lineage)
    target_lineage <- lineage
    lineage_programs <- registry[lineage == target_lineage]
    program_genes <- unique(membership[program_uid %chin% lineage_programs$program_uid, gene_symbol])
    program_genes <- intersect(program_genes, names(promoters_v49))
    peak_ranges <- parse_peaks(colnames(input$counts))
    promoter_ranges <- promoters_v49[names(promoters_v49) %chin% program_genes]
    hits <- findOverlaps(promoter_ranges, peak_ranges, ignore.strand = TRUE)
    by_gene <- split(subjectHits(hits), names(promoter_ranges)[queryHits(hits)])
    gene_counts <- matrix(
      0,
      nrow = length(program_genes),
      ncol = nrow(input$counts),
      dimnames = list(program_genes, input$donors$donor_id)
    )
    for (gene_symbol in names(by_gene)) {
      peak_index <- unique(by_gene[[gene_symbol]])
      gene_counts[gene_symbol, ] <- Matrix::rowSums(input$counts[, peak_index, drop = FALSE])
    }
    measured <- rowSums(gene_counts) >= 10 & rowSums(gene_counts > 0) >= 3
    gene_counts <- gene_counts[measured, , drop = FALSE]
    full_library <- Matrix::rowSums(input$counts)
    if (any(full_library <= 0)) stop("Zero donor peak library: ", cohort, "/", lineage)
    library_dge <- DGEList(counts = t(input$counts))
    library_dge <- calcNormFactors(library_dge, method = "TMM")
    effective_library <- library_dge$samples$lib.size * library_dge$samples$norm.factors
    logcpm <- cpm(gene_counts, lib.size = effective_library, log = TRUE, prior.count = 0.5)
    z <- t(scale(t(logcpm)))
    # A measured promoter with zero donor variance has a neutral standardized
    # contribution; it is not silently removed from the frozen weight mass.
    z[!is.finite(z)] <- 0

    scored <- list()
    for (uid in lineage_programs$program_uid) {
      answer <- score_program(z, membership[program_uid == uid])
      scored[[uid]] <- answer
      coverage_rows[[paste(cohort, uid, sep = "::")]] <- data.table(
        release_id = RELEASE_ID,
        program_release_id = PROGRAM_RELEASE_ID,
        cohort = cohort,
        program_uid = uid,
        membership_sha256 = lineage_programs[program_uid == uid, membership_sha256],
        lineage = lineage,
        n_program_genes = membership[program_uid == uid, uniqueN(gene_symbol)],
        n_measured_genes = answer$n_measured,
        retained_l1_weight = answer$retained,
        program_score_testable = answer$testable,
        program_score_reason = ifelse(answer$testable, "", "fewer_than_8_genes_or_20pct_weight"),
        top_weighted_gene = answer$top_gene
      )
      if (answer$testable) {
        score_rows[[paste(cohort, uid, sep = "::")]] <- data.table(
          release_id = RELEASE_ID,
          cohort = cohort,
          program_uid = uid,
          donor_id = names(answer$weighted),
          condition = input$donors$condition[match(names(answer$weighted), input$donors$donor_id)],
          lineage = lineage,
          score = as.numeric(answer$weighted),
          equal_weight_score = as.numeric(answer$equal),
          leave_top_gene_score = as.numeric(answer$leave)
        )
      }
    }

    eligible <- input$donors$contrast_eligible == TRUE & input$donors$condition %chin% c("NORMAL", "MASH")
    contrast_donors <- input$donors[eligible]
    counts_by_group <- contrast_donors[, .N, by = condition]
    n_normal <- counts_by_group[condition == "NORMAL", N]
    n_mash <- counts_by_group[condition == "MASH", N]
    if (length(n_normal) == 0L) n_normal <- 0L
    if (length(n_mash) == 0L) n_mash <- 0L
    contrast_testable <- n_normal >= 4L && n_mash >= 4L
    testable_uids <- names(scored)[vapply(scored, function(value) value$testable, logical(1))]

    fit_table <- NULL
    if (contrast_testable && length(testable_uids) > 0L) {
      sample_ids <- contrast_donors$donor_id
      score_matrix <- do.call(rbind, lapply(testable_uids, function(uid) scored[[uid]]$weighted[sample_ids]))
      rownames(score_matrix) <- testable_uids
      design <- model.matrix(~ factor(contrast_donors$condition, levels = c("NORMAL", "MASH")))
      if (qr(design)$rank != ncol(design)) stop("Program design is not full rank")
      fit <- eBayes(lmFit(score_matrix, design), robust = TRUE, trend = TRUE)
      top <- topTable(fit, coef = 2L, number = Inf, sort.by = "none")
      fit_table <- data.table(
        program_uid = rownames(top),
        effect = top$logFC,
        standard_error = fit$stdev.unscaled[, 2L] * fit$sigma,
        pvalue = top$P.Value
      )
    }
    model_fit <- contrast_testable && length(testable_uids) > 0L
    model_qc_rows[[paste(cohort, lineage, sep = "::")]] <- data.table(
      release_id = RELEASE_ID,
      cohort = cohort,
      lineage = lineage,
      model = "program_score_condition_only",
      contrast_testable = contrast_testable,
      n_normal = n_normal,
      n_mash = n_mash,
      design_rank = if (model_fit) qr(design)$rank else NA_integer_,
      n_model_columns = if (model_fit) ncol(design) else NA_integer_,
      n_tested_programs = if (is.null(fit_table)) 0L else nrow(fit_table)
    )

    for (uid in lineage_programs$program_uid) {
      answer <- scored[[uid]]
      fit_row <- if (!is.null(fit_table)) fit_table[program_uid == uid] else data.table()
      effect <- if (nrow(fit_row) == 1L) fit_row$effect else NA_real_
      equal_effect <- NA_real_
      leave_effect <- NA_real_
      if (answer$testable && contrast_testable) {
        cases <- input$donors$donor_id[eligible & input$donors$condition == "MASH"]
        controls <- input$donors$donor_id[eligible & input$donors$condition == "NORMAL"]
        equal_effect <- mean(answer$equal[cases]) - mean(answer$equal[controls])
        leave_effect <- mean(answer$leave[cases]) - mean(answer$leave[controls])
      }
      result_rows[[paste(cohort, uid, sep = "::")]] <- data.table(
        release_id = RELEASE_ID,
        program_release_id = PROGRAM_RELEASE_ID,
        cohort = cohort,
        program_uid = uid,
        membership_sha256 = lineage_programs[program_uid == uid, membership_sha256],
        program_name = lineage_programs[program_uid == uid, module_name],
        lineage = lineage,
        n_normal = n_normal,
        n_mash = n_mash,
        n_measured_genes = answer$n_measured,
        retained_l1_weight = answer$retained,
        program_score_testable = answer$testable,
        contrast_testable = answer$testable && contrast_testable,
        lineage_observability_state = ifelse(
          cohort == "GSE281367" && lineage == "t_nk",
          "source_dependent_combined_t_nk_label",
          "observed"
        ),
        testability_reason = fcase(
          !answer$testable, "fewer_than_8_genes_or_20pct_weight",
          !contrast_testable, "fewer_than_4_donors_per_condition",
          default = ""
        ),
        effect = effect,
        standard_error = if (nrow(fit_row) == 1L) fit_row$standard_error else NA_real_,
        pvalue = if (nrow(fit_row) == 1L) fit_row$pvalue else NA_real_,
        equal_weight_effect = equal_effect,
        leave_top_gene_effect = leave_effect,
        sensitivity_sign_agree = is.finite(effect) && is.finite(equal_effect) && is.finite(leave_effect) &&
          sign(effect) == sign(equal_effect) && sign(effect) == sign(leave_effect)
      )
    }
  }
}

coverage <- rbindlist(coverage_rows, fill = TRUE)
scores <- rbindlist(score_rows, fill = TRUE)
results <- rbindlist(result_rows, fill = TRUE)
model_qc <- rbindlist(model_qc_rows, fill = TRUE)
if (nrow(results) != 234L || uniqueN(results$program_uid) != 117L) stop("Expected 234 program-cohort rows")
results[, qvalue := NA_real_]
results[is.finite(pvalue), qvalue := p.adjust(pvalue, method = "BH"), by = cohort]
results[, within_source_state := fcase(
  !contrast_testable, "untestable",
  qvalue < 0.05 & sensitivity_sign_agree, "supported",
  default = "indeterminate"
)]

wide <- dcast(
  results,
  program_uid ~ cohort,
  value.var = c("effect", "qvalue", "contrast_testable", "sensitivity_sign_agree")
)
wide[, replicated :=
  contrast_testable_GSE244832 & contrast_testable_GSE281367 &
  qvalue_GSE244832 < 0.05 & qvalue_GSE281367 < 0.05 &
  effect_GSE244832 * effect_GSE281367 > 0 &
  sensitivity_sign_agree_GSE244832 & sensitivity_sign_agree_GSE281367]
wide[, cross_cohort_state := fcase(
  replicated, "supported",
  contrast_testable_GSE244832 & contrast_testable_GSE281367 &
    qvalue_GSE244832 < 0.05 & qvalue_GSE281367 < 0.05 &
    effect_GSE244832 * effect_GSE281367 < 0, "discordant",
  xor(qvalue_GSE244832 < 0.05, qvalue_GSE281367 < 0.05), "source_dependent",
  !contrast_testable_GSE244832 | !contrast_testable_GSE281367, "untestable",
  default = "indeterminate"
)]
results[wide, on = "program_uid", `:=`(
  replicated = i.replicated,
  cross_cohort_state = i.cross_cohort_state
)]

setorder(coverage, cohort, lineage, program_uid)
setorder(scores, cohort, lineage, program_uid, donor_id)
setorder(results, cohort, lineage, program_uid)
fwrite(coverage, file.path(OUT, "program_measurement_coverage.tsv"), sep = "\t")
fwrite(scores, file.path(OUT, "program_donor_scores.tsv.gz"), sep = "\t")
fwrite(results, file.path(OUT, "program_atac_results.tsv"), sep = "\t")
fwrite(model_qc, file.path(OUT, "program_model_qc.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
message("[PROGRAM] Wrote ", OUT)
