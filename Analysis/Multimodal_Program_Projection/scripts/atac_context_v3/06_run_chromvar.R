#!/usr/bin/env Rscript

# Official chromVAR donor-level motif deviations and separate-cohort inference.

suppressPackageStartupMessages({
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(chromVAR)
  library(data.table)
  library(digest)
  library(JASPAR2024)
  library(limma)
  library(Matrix)
  library(motifmatchr)
  library(SummarizedExperiment)
  library(TFBSTools)
})

set.seed(42)

BASE <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
RELEASE_ID <- "atac-context-v3-candidate-2026-08-11-r1"
CANDIDATE <- Sys.getenv(
  "ATAC_V3_CANDIDATE_ROOT",
  file.path(BASE, "Analysis/Multimodal_Program_Projection/candidates", RELEASE_ID)
)
EXPECTED <- normalizePath(
  file.path(BASE, "Analysis/Multimodal_Program_Projection/candidates", RELEASE_ID),
  mustWork = FALSE
)
if (!identical(normalizePath(CANDIDATE, mustWork = FALSE), EXPECTED)) stop("Unsafe candidate root")
OUT <- file.path(CANDIDATE, "chromvar")
if (dir.exists(OUT)) stop("Refusing to overwrite chromVAR output: ", OUT)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

read_sparse_gz <- function(path) {
  con <- gzfile(path, "rb")
  on.exit(close(con), add = TRUE)
  as(readMM(con), "dgCMatrix")
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

jaspar_object <- JASPAR2024()
jaspar_sqlite <- db(jaspar_object)
if (!file.exists(jaspar_sqlite)) stop("JASPAR2024 SQLite database is unavailable")
motifs <- getMatrixSet(
  jaspar_sqlite,
  opts = list(collection = "CORE", tax_group = "vertebrates", all_versions = FALSE)
)
if (length(motifs) == 0L) stop("JASPAR2024 CORE vertebrate motif universe is empty")
motif_ids <- names(motifs)
if (anyDuplicated(motif_ids)) stop("Duplicated JASPAR2024 motif IDs")
fwrite(data.table(
  release_id = RELEASE_ID,
  jaspar_package = "JASPAR2024",
  database_path = jaspar_sqlite,
  database_sha256 = digest(jaspar_sqlite, algo = "sha256", file = TRUE, serialize = FALSE),
  collection = "CORE",
  tax_group = "vertebrates",
  all_versions = FALSE,
  n_motifs = length(motifs)
), file.path(OUT, "jaspar_database.tsv"), sep = "\t")

coverage_rows <- list()
deviation_rows <- list()
result_rows <- list()
model_qc_rows <- list()
peak_filter_rows <- list()

for (cohort in c("GSE244832", "GSE281367")) {
  for (lineage in c("hepatocyte", "stellate", "macrophage", "cholangiocyte", "t_nk")) {
    message("[chromVAR] ", cohort, " / ", lineage)
    input <- load_counts(cohort, lineage)
    n_input_peaks <- ncol(input$counts)
    nonzero_library <- Matrix::colSums(input$counts) > 0
    if (!any(nonzero_library)) stop("No nonzero-library peaks for ", cohort, "/", lineage)
    input$counts <- input$counts[, nonzero_library, drop = FALSE]
    input$peaks <- input$peaks[nonzero_library]
    peak_filter_rows[[paste(cohort, lineage, sep = "::")]] <- data.table(
      release_id = RELEASE_ID,
      cohort = cohort,
      lineage = lineage,
      filter_method = "condition_blind_nonzero_total_donor_count",
      n_consensus_peaks = n_input_peaks,
      n_zero_library_peaks_excluded = n_input_peaks - ncol(input$counts),
      n_chromvar_input_peaks = ncol(input$counts)
    )
    row_ranges <- parse_peaks(colnames(input$counts))
    names(row_ranges) <- colnames(input$counts)
    se <- SummarizedExperiment(
      assays = list(counts = t(input$counts)),
      rowRanges = row_ranges,
      colData = DataFrame(input$donors)
    )
    se <- addGCBias(se, genome = BSgenome.Hsapiens.UCSC.hg38)
    motif_matches <- matchMotifs(
      motifs,
      se,
      genome = BSgenome.Hsapiens.UCSC.hg38,
      out = "matches"
    )
    match_matrix <- as(assay(motif_matches), "lgCMatrix")
    matched_peaks <- Matrix::colSums(match_matrix)
    background <- getBackgroundPeaks(se, niterations = 50)
    deviations <- computeDeviations(
      object = se,
      annotations = motif_matches,
      background_peaks = background
    )
    scores <- deviationScores(deviations)
    if (!identical(rownames(scores), motif_ids)) {
      scores <- scores[motif_ids, , drop = FALSE]
    }
    variance <- apply(scores, 1L, var, na.rm = TRUE)
    motif_testable <- matched_peaks[motif_ids] >= 10L & is.finite(variance) & variance > 0

    coverage_rows[[paste(cohort, lineage, sep = "::")]] <- data.table(
      release_id = RELEASE_ID,
      cohort = cohort,
      lineage = lineage,
      motif_id = motif_ids,
      n_matched_peaks = as.integer(matched_peaks[motif_ids]),
      donor_variance = as.numeric(variance[motif_ids]),
      motif_testable = as.logical(motif_testable),
      motif_testability_reason = ifelse(
        motif_testable,
        "",
        ifelse(matched_peaks[motif_ids] < 10L, "fewer_than_10_matched_peaks", "zero_donor_variance")
      )
    )
    deviation_rows[[paste(cohort, lineage, sep = "::")]] <- data.table(
      motif_id = rep(rownames(scores), times = ncol(scores)),
      donor_id = rep(colnames(scores), each = nrow(scores)),
      deviation_z = as.numeric(scores),
      release_id = RELEASE_ID,
      cohort = cohort,
      lineage = lineage
    )

    match_dir <- file.path(OUT, "motif_matches", cohort, lineage)
    dir.create(match_dir, recursive = TRUE, showWarnings = FALSE)
    Matrix::writeMM(match_matrix, file.path(match_dir, "matches.mtx"))
    fwrite(data.table(peak_coordinate = names(row_ranges)), file.path(match_dir, "peaks.tsv"), sep = "\t")
    fwrite(data.table(motif_id = motif_ids), file.path(match_dir, "motifs.tsv"), sep = "\t")

    eligible <- input$donors$contrast_eligible == TRUE & input$donors$condition %chin% c("NORMAL", "MASH")
    fit_donors <- input$donors[eligible]
    n_normal <- sum(fit_donors$condition == "NORMAL")
    n_mash <- sum(fit_donors$condition == "MASH")
    contrast_testable <- n_normal >= 4L && n_mash >= 4L
    tested_ids <- motif_ids[motif_testable]
    fit_table <- data.table(
      motif_id = character(),
      effect = numeric(),
      standard_error = numeric(),
      pvalue = numeric(),
      qvalue = numeric()
    )
    if (contrast_testable && length(tested_ids) > 0L) {
      design <- model.matrix(~ factor(fit_donors$condition, levels = c("NORMAL", "MASH")))
      if (qr(design)$rank != ncol(design)) stop("chromVAR design is not full rank")
      fit <- eBayes(
        lmFit(scores[tested_ids, fit_donors$donor_id, drop = FALSE], design),
        robust = TRUE,
        trend = TRUE
      )
      top <- topTable(fit, coef = 2L, number = Inf, sort.by = "none")
      fit_table <- data.table(
        motif_id = rownames(top),
        effect = top$logFC,
        standard_error = fit$stdev.unscaled[, 2L] * fit$sigma,
        pvalue = top$P.Value,
        qvalue = p.adjust(top$P.Value, method = "BH")
      )
    }
    model_fit <- contrast_testable && length(tested_ids) > 0L
    model_qc_rows[[paste(cohort, lineage, sep = "::")]] <- data.table(
      release_id = RELEASE_ID,
      cohort = cohort,
      lineage = lineage,
      model = "chromvar_deviation_condition_only",
      contrast_testable = contrast_testable,
      n_normal = n_normal,
      n_mash = n_mash,
      design_rank = if (model_fit) qr(design)$rank else NA_integer_,
      n_model_columns = if (model_fit) ncol(design) else NA_integer_,
      n_tested_motifs = if (nrow(fit_table) == 0L) 0L else nrow(fit_table)
    )
    family <- data.table(motif_id = motif_ids)
    family <- merge(family, fit_table, by = "motif_id", all.x = TRUE, sort = FALSE)
    family[, `:=`(
      release_id = RELEASE_ID,
      cohort = cohort,
      lineage = lineage,
      n_normal = n_normal,
      n_mash = n_mash,
      motif_testable = motif_testable[match(motif_id, motif_ids)],
      contrast_testable = motif_testable[match(motif_id, motif_ids)] & contrast_testable,
      within_source_state = fcase(
        !(motif_testable[match(motif_id, motif_ids)] & contrast_testable), "untestable",
        qvalue < 0.05, "supported",
        default = "indeterminate"
      )
    )]
    result_rows[[paste(cohort, lineage, sep = "::")]] <- family
  }
}

coverage <- rbindlist(coverage_rows, fill = TRUE)
deviation_table <- rbindlist(deviation_rows, fill = TRUE)
results <- rbindlist(result_rows, fill = TRUE)
model_qc <- rbindlist(model_qc_rows, fill = TRUE)
peak_filter <- rbindlist(peak_filter_rows, fill = TRUE)
wide <- dcast(results, lineage + motif_id ~ cohort, value.var = c("effect", "qvalue", "contrast_testable"))
wide[, replicated :=
  contrast_testable_GSE244832 & contrast_testable_GSE281367 &
  qvalue_GSE244832 < 0.05 & qvalue_GSE281367 < 0.05 &
  effect_GSE244832 * effect_GSE281367 > 0]
wide[, cross_cohort_state := fcase(
  replicated, "supported",
  contrast_testable_GSE244832 & contrast_testable_GSE281367 &
    qvalue_GSE244832 < 0.05 & qvalue_GSE281367 < 0.05 &
    effect_GSE244832 * effect_GSE281367 < 0, "discordant",
  (qvalue_GSE244832 < 0.05) != (qvalue_GSE281367 < 0.05), "source_dependent",
  !contrast_testable_GSE244832 | !contrast_testable_GSE281367, "untestable",
  default = "indeterminate"
)]
results[wide, on = c("lineage", "motif_id"), `:=`(
  replicated = i.replicated,
  cross_cohort_state = i.cross_cohort_state
)]

setorder(coverage, cohort, lineage, motif_id)
setorder(deviation_table, cohort, lineage, motif_id, donor_id)
setorder(results, cohort, lineage, motif_id)
fwrite(coverage, file.path(OUT, "motif_measurement_coverage.tsv.gz"), sep = "\t")
fwrite(deviation_table, file.path(OUT, "donor_deviations.tsv.gz"), sep = "\t")
fwrite(results, file.path(OUT, "chromvar_results.tsv.gz"), sep = "\t")
fwrite(model_qc, file.path(OUT, "chromvar_model_qc.tsv"), sep = "\t")
fwrite(peak_filter, file.path(OUT, "chromvar_peak_filter_qc.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
message("[chromVAR] Wrote ", OUT)
