#!/usr/bin/env Rscript

# Reproduce the released Kamzolas et al. UCAM/VCU PC1 and freeze discovery
# loadings for transport into the MASLD Resource cohorts.
#
# Required environment variables:
#   KAMZOLAS_V1_ROOT            clean checkout of the authors' v1.0 release
#   HAC_OUT_ROOT                new, write-once analysis output root
#   HAC_RESOURCE_GENE_IDS_PATH  TSV with a unique gene_id_base column
# Optional:
#   KAMZOLAS_V1_MANIFEST_PATH   TSV or sha256sum-style manifest for the checkout

options(digits = 17, scipen = 999)

EXPECTED_COMMIT <- "cbcf7764b625d851ca4e2e58b400b921cfb8f39a"
EXPECTED_TAG <- "v1.0"
EXPECTED_N_INPUT <- 136L
EXPECTED_OUTLIER <- "Sample 5"
EXPECTED_N_DISCOVERY <- 135L
EXPECTED_N_SIGNATURE <- 145L
EXPECTED_N_RESOURCE_SIGNATURE <- 139L
MIN_SIGNATURE_COVERAGE <- 0.90
MIN_RELEASED_PC1_SPEARMAN <- 0.99
SEED <- 20260817L
OPERATIONAL_PC1_MULTIPLIER <- -1
OPERATIONAL_SIGN_RELATIVE <- paste0(
  "src/01.preprocessing_and_trajectory_analysis/",
  "a.Linking_disease_state_to_phenotypes/",
  "UCAM_VCU - correlation between difference in positions on the trajectory and histological differences.R"
)
OPERATIONAL_SIGN_PATTERN <- paste0(
  "merged_data\\$PC1[[:space:]]*<-[[:space:]]*",
  "-[[:space:]]*merged_data\\$PC1"
)

fail <- function(...) stop(paste0(...), call. = FALSE)

assert_true <- function(value, message) {
  if (!isTRUE(value)) fail(message)
  invisible(TRUE)
}

require_env_path <- function(name, type = c("dir", "file")) {
  type <- match.arg(type)
  value <- Sys.getenv(name, "")
  assert_true(nzchar(value), paste0(name, " is unset"))
  value <- normalizePath(value, mustWork = FALSE)
  exists <- if (type == "dir") dir.exists(value) else file.exists(value)
  assert_true(exists, paste0(name, " does not identify an existing ", type, ": ", value))
  value
}

canonical_name <- function(x) tolower(gsub("[^[:alnum:]]", "", x))

find_column <- function(x, candidates, label) {
  idx <- match(canonical_name(candidates), canonical_name(names(x)), nomatch = 0L)
  idx <- idx[idx > 0L]
  assert_true(length(idx) >= 1L, paste0("Could not identify ", label, " column"))
  names(x)[idx[[1L]]]
}

strip_ensembl_version <- function(x) sub("[.][0-9]+$", "", trimws(as.character(x)))

sha256_file <- function(path) {
  assert_true(file.exists(path), paste0("Cannot hash absent file: ", path))
  ans <- system2("sha256sum", shQuote(path), stdout = TRUE, stderr = TRUE)
  assert_true(is.null(attr(ans, "status")) || attr(ans, "status") == 0L,
              paste0("sha256sum failed for ", path))
  sub("[[:space:]].*$", "", ans[[1L]])
}

git_output <- function(root, args) {
  ans <- system2("git", c("-C", shQuote(root), args), stdout = TRUE, stderr = TRUE)
  assert_true(is.null(attr(ans, "status")) || attr(ans, "status") == 0L,
              paste0("git ", paste(args, collapse = " "), " failed: ", paste(ans, collapse = "\n")))
  ans
}

write_tsv_once <- function(x, path) {
  assert_true(!file.exists(path), paste0("Refusing to overwrite output: ", path))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  invisible(path)
}

save_rds_once <- function(x, path) {
  assert_true(!file.exists(path), paste0("Refusing to overwrite output: ", path))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(x, path, compress = TRUE)
  invisible(path)
}

write_lines_once <- function(x, path) {
  assert_true(!file.exists(path), paste0("Refusing to overwrite output: ", path))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(x, path, useBytes = TRUE)
  invisible(path)
}

log_step <- function(...) {
  cat(format(Sys.time(), "%H:%M:%S"), "|", paste0(...), "\n")
  flush.console()
}

manifest_rows <- function(path, vendor_root) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines)) & !grepl("^[[:space:]]*#", lines)]
  assert_true(length(lines) > 0L, paste0("Checksum manifest is empty: ", path))

  sha_style <- all(grepl("^[[:xdigit:]]{64}[[:space:]]+[*]?.+", lines))
  if (sha_style) {
    sha <- tolower(sub("^([[:xdigit:]]{64}).*$", "\\1", lines))
    rel <- sub("^[[:xdigit:]]{64}[[:space:]]+[*]?", "", lines)
    tab <- data.frame(relative_path = rel, sha256 = sha, stringsAsFactors = FALSE)
  } else {
    tab <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE,
                      comment.char = "", quote = "\"")
    path_col <- find_column(tab, c("relative_path", "path", "file", "filename"),
                            "manifest path")
    sha_col <- find_column(tab, c("sha256", "sha256sum", "checksum"),
                           "manifest SHA-256")
    tab <- data.frame(relative_path = tab[[path_col]], sha256 = tolower(tab[[sha_col]]),
                      stringsAsFactors = FALSE)
  }

  root_prefix <- paste0(normalizePath(vendor_root, mustWork = TRUE), .Platform$file.sep)
  tab$relative_path <- trimws(as.character(tab$relative_path))
  is_absolute <- startsWith(tab$relative_path, .Platform$file.sep)
  if (any(is_absolute)) {
    absolute <- normalizePath(tab$relative_path[is_absolute], mustWork = FALSE)
    assert_true(all(startsWith(absolute, root_prefix)),
                "Manifest contains an absolute path outside KAMZOLAS_V1_ROOT")
    tab$relative_path[is_absolute] <- substring(absolute, nchar(root_prefix) + 1L)
  }
  tab$relative_path <- sub("^[.][/\\\\]", "", tab$relative_path)
  tab$relative_path <- gsub("\\\\", "/", tab$relative_path)
  assert_true(!any(grepl("(^|/)[.][.](/|$)", tab$relative_path)),
              "Manifest contains a parent-directory path")
  assert_true(all(grepl("^[[:xdigit:]]{64}$", tab$sha256)),
              "Manifest contains a malformed SHA-256 value")
  assert_true(!anyDuplicated(tab$relative_path), "Manifest paths are not unique")
  tab
}

verify_manifest <- function(manifest_path, vendor_root, required_relative_paths) {
  tab <- manifest_rows(manifest_path, vendor_root)
  idx <- match(required_relative_paths, tab$relative_path)
  assert_true(!anyNA(idx),
              paste0("Manifest omits required inputs: ",
                     paste(required_relative_paths[is.na(idx)], collapse = ", ")))
  checked <- tab[idx, , drop = FALSE]
  checked$observed_sha256 <- vapply(
    file.path(vendor_root, checked$relative_path), sha256_file, character(1L)
  )
  checked$pass <- checked$sha256 == checked$observed_sha256
  assert_true(all(checked$pass),
              paste0("Manifest checksum mismatch: ",
                     paste(checked$relative_path[!checked$pass], collapse = ", ")))
  checked
}

read_matrix_csv <- function(path) {
  tab <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  assert_true(ncol(tab) >= 2L, paste0("Matrix has fewer than two columns: ", path))
  ids <- strip_ensembl_version(tab[[1L]])
  assert_true(all(nzchar(ids)) && !anyDuplicated(ids),
              paste0("Gene IDs are empty or duplicated after version stripping: ", path))
  mat <- as.matrix(tab[-1L])
  suppressWarnings(storage.mode(mat) <- "double")
  assert_true(!anyNA(mat) && all(is.finite(mat)), paste0("Non-finite matrix value in ", path))
  rownames(mat) <- ids
  mat
}

safe_correlation <- function(score, histology, method) {
  keep <- is.finite(score) & is.finite(histology)
  assert_true(sum(keep) >= 3L, "Fewer than three complete participants for a correlation")
  c(n = sum(keep), estimate = unname(cor(score[keep], histology[keep], method = method)))
}

correlation_table <- function(scores, metadata) {
  score_columns <- c(
    "pc1_reproduced_raw", "pc1_released_raw",
    "pc1_reproduced_operational", "pc1_released_operational"
  )
  histology_columns <- c("NAS", "fibrosis", "steatosis", "ballooning", "inflammation")
  rows <- list()
  k <- 0L
  for (score_name in score_columns) {
    for (histology_name in histology_columns) {
      h <- metadata[[histology_name]]
      s <- scores[[score_name]]
      for (method in c("pearson", "spearman")) {
        donor <- safe_correlation(s, h, method)
        k <- k + 1L
        rows[[k]] <- data.frame(
          score = score_name, aggregation = "donor", histology = histology_name,
          method = method, n_units = unname(donor[["n"]]),
          estimate = unname(donor[["estimate"]]), stringsAsFactors = FALSE
        )

        grouped <- aggregate(s, list(histology_value = h), mean)
        names(grouped)[[2L]] <- "mean_score"
        group_result <- safe_correlation(grouped$mean_score, grouped$histology_value, method)
        k <- k + 1L
        rows[[k]] <- data.frame(
          score = score_name, aggregation = "histology_group_mean",
          histology = histology_name, method = method,
          n_units = unname(group_result[["n"]]),
          estimate = unname(group_result[["estimate"]]), stringsAsFactors = FALSE
        )
      }
    }
  }
  do.call(rbind, rows)
}

audit_code_paths <- function(vendor_root) {
  src_root <- file.path(vendor_root, "src")
  code_files <- list.files(src_root, pattern = "[.](R|r|py)$", recursive = TRUE,
                           full.names = TRUE)
  assert_true(length(code_files) > 0L, "No released source files found for code-path audit")
  patterns <- list(
    pc1_sorted_file = "PC1_sorted_samples",
    operational_pc1_negation = OPERATIONAL_SIGN_PATTERN,
    operational_rf_inverse_pc1 = "1[[:space:]]*-[[:space:]]*metadata\\$PC1",
    pseudotime_extraction = "slingPseudotime[[:space:]]*\\(|slingCurveWeights[[:space:]]*\\(",
    slingshot_fit = "(^|[^[:alnum:]_.])slingshot[[:space:]]*\\(",
    slingshot_plot = "SlingshotDataSet[[:space:]]*\\("
  )
  rows <- list()
  k <- 0L
  for (path in code_files) {
    text <- suppressWarnings(readLines(path, warn = FALSE))
    for (evidence_type in names(patterns)) {
      hits <- grep(patterns[[evidence_type]], text, perl = TRUE)
      if (!length(hits)) next
      for (line_number in hits) {
        k <- k + 1L
        relative <- substring(normalizePath(path),
                              nchar(normalizePath(vendor_root)) + 2L)
        rows[[k]] <- data.frame(
          code_path = relative,
          line_number = line_number,
          evidence_type = evidence_type,
          downstream = !grepl("^src/01[.]preprocessing_and_trajectory_analysis/[^/]+$", relative),
          line_text = gsub("[\t\r\n]", " ", trimws(text[[line_number]])),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (!length(rows)) {
    audit <- data.frame(code_path = NA_character_, line_number = NA_integer_,
                        evidence_type = "no_pattern_matches", downstream = NA,
                        line_text = NA_character_, stringsAsFactors = FALSE)
  } else {
    audit <- do.call(rbind, rows)
  }
  if (!any(audit$evidence_type == "pseudotime_extraction")) {
    audit <- rbind(
      audit,
      data.frame(code_path = NA_character_, line_number = NA_integer_,
                 evidence_type = "pseudotime_extraction_absent", downstream = NA,
                 line_text = paste0("No slingPseudotime() or slingCurveWeights() call in ",
                                    "released src/ code"), stringsAsFactors = FALSE)
    )
  }
  audit
}

summary_row <- function(metric, value, criterion = NA_character_, pass = NA,
                        detail = NA_character_) {
  data.frame(metric = metric, value = as.character(value), criterion = criterion,
             pass = if (is.na(pass)) NA else as.logical(pass), detail = detail,
             stringsAsFactors = FALSE)
}

# ------------------------------------------------------------------------- I/O

vendor_root <- require_env_path("KAMZOLAS_V1_ROOT", "dir")
out_root_raw <- Sys.getenv("HAC_OUT_ROOT", "")
assert_true(nzchar(out_root_raw), "HAC_OUT_ROOT is unset")
out_root <- normalizePath(out_root_raw, mustWork = FALSE)
resource_gene_path <- require_env_path("HAC_RESOURCE_GENE_IDS_PATH", "file")
out_dir <- file.path(out_root, "reproduction")

required_relative <- c(
  "data/merged_counts.csv",
  "data/metadata.csv",
  "data/PC1_sorted_samples.csv",
  "data/batch_corrected_counts_(dataset+gender).csv",
  "data/145 (from200) most important genes (variable importance with MeanDecreaseAccuracy>0).csv",
  OPERATIONAL_SIGN_RELATIVE
)
required_paths <- setNames(file.path(vendor_root, required_relative), required_relative)
missing_inputs <- names(required_paths)[!file.exists(required_paths)]
assert_true(!length(missing_inputs),
            paste0("Kamzolas v1 checkout is missing: ", paste(missing_inputs, collapse = ", ")))

output_paths <- file.path(out_dir, c(
  "participant_scores.tsv",
  "histology_correlations.tsv",
  "reproduction_summary.tsv",
  "code_path_audit.tsv",
  "pca_variance_explained.tsv",
  "signature_genes.tsv",
  "missing_signature_genes.tsv",
  "discovery_signature_model.rds",
  "discovery_normalized_expression.rds",
  "discovery_metadata.tsv",
  "input_checksums.tsv",
  "sessionInfo.txt"
))
existing_outputs <- output_paths[file.exists(output_paths)]
assert_true(!length(existing_outputs),
            paste0("Refusing to mix with existing reproduction outputs: ",
                   paste(existing_outputs, collapse = ", ")))

for (pkg in c("limma", "sva")) {
  assert_true(requireNamespace(pkg, quietly = TRUE),
              paste0("Required installed R package is unavailable: ", pkg))
}
set.seed(SEED)

# -------------------------------------------------------------- provenance

log_step("verifying frozen Kamzolas v1 checkout")
observed_commit <- trimws(git_output(vendor_root, c("rev-parse", "HEAD"))[[1L]])
assert_true(identical(observed_commit, EXPECTED_COMMIT),
            paste0("Expected Kamzolas commit ", EXPECTED_COMMIT, "; observed ", observed_commit))
observed_tags <- trimws(git_output(vendor_root, c("tag", "--points-at", "HEAD")))
assert_true(EXPECTED_TAG %in% observed_tags, paste0("Commit is not tagged ", EXPECTED_TAG))
tracked_status <- git_output(vendor_root, c("status", "--porcelain", "--untracked-files=no"))
tracked_status <- tracked_status[nzchar(tracked_status)]
assert_true(length(tracked_status) == 0L, "Kamzolas v1 checkout has modified tracked files")

operational_sign_source <- suppressWarnings(readLines(
  required_paths[[OPERATIONAL_SIGN_RELATIVE]], warn = FALSE
))
operational_sign_hits <- grep(
  OPERATIONAL_SIGN_PATTERN, operational_sign_source, perl = TRUE
)
assert_true(length(operational_sign_hits) == 1L,
            "Released source no longer contains exactly one operational PC1 negation")

manifest_path <- Sys.getenv("KAMZOLAS_V1_MANIFEST_PATH", "")
manifest_candidates <- c(
  manifest_path,
  file.path(out_root, "inputs", "kamzolas_v1_manifest.tsv"),
  file.path(vendor_root, "SHA256SUMS"),
  file.path(vendor_root, "MANIFEST.sha256")
)
manifest_candidates <- unique(manifest_candidates[nzchar(manifest_candidates)])
manifest_candidates <- manifest_candidates[file.exists(manifest_candidates)]
manifest_used <- if (length(manifest_candidates)) normalizePath(manifest_candidates[[1L]]) else NA_character_
manifest_check <- NULL
if (!is.na(manifest_used)) {
  log_step("verifying supplied checksum manifest: ", manifest_used)
  manifest_check <- verify_manifest(manifest_used, vendor_root, required_relative)
}

# ----------------------------------------------------------- exact upstream

log_step("reading released counts and metadata")
metadata_raw <- read.csv(required_paths[["data/metadata.csv"]], check.names = FALSE,
                         stringsAsFactors = FALSE)
sample_col <- find_column(metadata_raw, c("Sample name", "Sample.name"), "sample")
dataset_col <- find_column(metadata_raw, "Dataset", "dataset")
sex_col <- find_column(metadata_raw, "Sex", "sex")
nas_col <- find_column(metadata_raw, "NAS", "NAS")
fibrosis_col <- find_column(metadata_raw, "Fibrosis", "fibrosis")
steatosis_col <- find_column(metadata_raw, "Steatosis", "steatosis")
ballooning_col <- find_column(metadata_raw, "Ballooning", "ballooning")
inflammation_col <- find_column(metadata_raw, "Inflammation", "inflammation")

assert_true(nrow(metadata_raw) == EXPECTED_N_INPUT,
            paste0("Expected ", EXPECTED_N_INPUT, " metadata rows; observed ", nrow(metadata_raw)))
outlier_index <- which(metadata_raw[[sample_col]] == EXPECTED_OUTLIER)
assert_true(identical(outlier_index, 5L), "The released fifth-sample outlier is not Sample 5")
metadata <- metadata_raw[-outlier_index, , drop = FALSE]
sample_ids <- as.character(metadata[[sample_col]])
assert_true(length(sample_ids) == EXPECTED_N_DISCOVERY && !anyDuplicated(sample_ids),
            "Discovery participant count or identity uniqueness changed")

counts <- read_matrix_csv(required_paths[["data/merged_counts.csv"]])
assert_true(all(sample_ids %in% colnames(counts)), "Counts omit a discovery participant")
counts <- counts[, sample_ids, drop = FALSE]
assert_true(all(counts >= 0) && all(counts == floor(counts)),
            "Released count matrix is not non-negative integer counts")
keep <- rowSums(counts) > ncol(counts)
assert_true(any(keep), "Low-expression filter retained no genes")
counts_filtered <- counts[keep, , drop = FALSE]

log_step("log2(1+count), quantile normalization, and NAS-protected ComBat")
# limma's pure-R implementation is numerically equivalent to the released
# preprocessCore call. It avoids preprocessCore's cgroup-unsafe thread probe,
# which can abort before returning a value on the project SLURM nodes.
normalized <- limma::normalizeQuantiles(log2(1 + counts_filtered))
dimnames(normalized) <- dimnames(counts_filtered)
batch <- paste(metadata[[dataset_col]], metadata[[sex_col]])
mod <- model.matrix(~ as.factor(metadata[[nas_col]]))
corrected <- sva::ComBat(dat = normalized, batch = batch, mod = mod,
                         par.prior = TRUE, prior.plots = FALSE)
dimnames(corrected) <- dimnames(normalized)
assert_true(identical(colnames(corrected), sample_ids), "ComBat changed participant order")
assert_true(!anyNA(corrected) && all(is.finite(corrected)), "ComBat produced non-finite values")

log_step("PCA, released raw-PC1 reproduction, and operational -PC1 orientation")
full_pca <- prcomp(t(corrected), center = TRUE, scale. = FALSE)
variance <- full_pca$sdev^2
variance_table <- data.frame(
  component = paste0("PC", seq_along(variance)),
  eigenvalue = variance,
  variance_fraction = variance / sum(variance),
  variance_percent = 100 * variance / sum(variance),
  stringsAsFactors = FALSE
)

released_scores_raw <- read.csv(required_paths[["data/PC1_sorted_samples.csv"]],
                                check.names = FALSE, stringsAsFactors = FALSE)
released_sample_ids <- as.character(released_scores_raw[[1L]])
released_pc1_column <- find_column(released_scores_raw, "PC1", "PC1")
released_pc1_file_order <- as.numeric(released_scores_raw[[released_pc1_column]])
assert_true(all(is.finite(released_pc1_file_order)),
            "Released PC1 contains a non-finite value")
released_pc1_file_strictly_descending <- all(diff(released_pc1_file_order) < 0)
assert_true(released_pc1_file_strictly_descending,
            "Released PC1 file is not in strictly descending raw-PC1 order")
released_file_operational_percentile <-
  (rank(OPERATIONAL_PC1_MULTIPLIER * released_pc1_file_order,
        ties.method = "average") - 1) / (length(released_pc1_file_order) - 1)
assert_true(isTRUE(all.equal(
  released_file_operational_percentile,
  (seq_along(released_pc1_file_order) - 1) / (length(released_pc1_file_order) - 1),
  tolerance = 0
)), "Released file order is not the exact rank of operational -raw PC1")
assert_true(all(sample_ids %in% released_sample_ids),
            "Released PC1 file omits a discovery participant")
released_file_position <- match(sample_ids, released_sample_ids)
released_scores_raw <- released_scores_raw[
  released_file_position, , drop = FALSE
]
released_pc1_raw <- as.numeric(
  released_scores_raw[[released_pc1_column]]
)
released_pc2 <- as.numeric(released_scores_raw[[find_column(released_scores_raw, "PC2", "PC2")]])
assert_true(all(is.finite(released_pc1_raw)), "Released PC1 contains a non-finite value")
orientation_full_raw <- sign(cor(full_pca$x[, 1L], released_pc1_raw, method = "pearson"))
assert_true(orientation_full_raw != 0,
            "Reproduced PC1 has zero correlation with released raw PC1")
pc1_reproduced_raw <- unname(full_pca$x[, 1L] * orientation_full_raw)
released_pc1_operational <- released_pc1_raw * OPERATIONAL_PC1_MULTIPLIER
pc1_reproduced_operational <- pc1_reproduced_raw * OPERATIONAL_PC1_MULTIPLIER
orientation_full_operational <- orientation_full_raw * OPERATIONAL_PC1_MULTIPLIER
released_operational_rank_percentile <-
  (rank(released_pc1_operational, ties.method = "average") - 1) /
  (length(released_pc1_operational) - 1)
operational_percentile_file_order_max_abs_difference <- max(abs(
  released_operational_rank_percentile -
    (released_file_position - 1) / (length(released_file_position) - 1)
))
assert_true(isTRUE(all.equal(
  released_operational_rank_percentile,
  (released_file_position - 1) / (length(released_file_position) - 1),
  tolerance = 0
)), "Operational -raw-PC1 percentile disagrees with released descending file order")
pc1_spearman <- cor(pc1_reproduced_raw, released_pc1_raw, method = "spearman")
pc1_pearson <- cor(pc1_reproduced_raw, released_pc1_raw, method = "pearson")
operational_pc1_spearman <- cor(
  pc1_reproduced_operational, released_pc1_operational, method = "spearman"
)
assert_true(abs(pc1_spearman) >= MIN_RELEASED_PC1_SPEARMAN,
            paste0("Released-PC1 reproduction failed: |Spearman rho| = ", pc1_spearman))
assert_true(isTRUE(all.equal(
  pc1_reproduced_operational,
  OPERATIONAL_PC1_MULTIPLIER * pc1_reproduced_raw,
  tolerance = 0
)), "Operational PC1 is not the exact released-source sign transformation")

# The released corrected matrix is not substituted for the recomputation. It is
# used only to quantify numerical agreement with the independently recomputed matrix.
released_corrected <- read_matrix_csv(
  required_paths[["data/batch_corrected_counts_(dataset+gender).csv"]]
)
assert_true(setequal(rownames(corrected), rownames(released_corrected)) &&
              setequal(colnames(corrected), colnames(released_corrected)),
            "Released and recomputed corrected matrices have different dimensions or identifiers")
released_corrected <- released_corrected[rownames(corrected), colnames(corrected), drop = FALSE]
corrected_matrix_pearson <- cor(as.vector(corrected), as.vector(released_corrected))
corrected_matrix_max_abs_diff <- max(abs(corrected - released_corrected))

# ---------------------------------------------------------- frozen signatures

log_step("freezing 145-gene and full-transcriptome discovery models")
signature_raw <- read.csv(
  required_paths[["data/145 (from200) most important genes (variable importance with MeanDecreaseAccuracy>0).csv"]],
  check.names = FALSE, stringsAsFactors = FALSE
)
signature_id_col <- find_column(signature_raw, c("Sample.1", "ensembl", "gene_id"),
                                "signature Ensembl ID")
signature_symbol_col <- find_column(signature_raw, c("SYMBOL", "gene_symbol"),
                                    "signature gene symbol")
signature_ids <- strip_ensembl_version(signature_raw[[signature_id_col]])
signature_symbols <- as.character(signature_raw[[signature_symbol_col]])
assert_true(length(signature_ids) == EXPECTED_N_SIGNATURE && !anyDuplicated(signature_ids),
            "Released 145-gene signature cardinality or uniqueness changed")
assert_true(all(signature_ids %in% rownames(corrected)),
            "A released signature gene is absent from the discovery expression matrix")

resource_genes <- read.delim(resource_gene_path, check.names = FALSE,
                             stringsAsFactors = FALSE)
assert_true("gene_id_base" %in% names(resource_genes),
            "HAC_RESOURCE_GENE_IDS_PATH lacks required gene_id_base header")
resource_ids <- strip_ensembl_version(resource_genes$gene_id_base)
assert_true(all(nzchar(resource_ids)) && !anyDuplicated(resource_ids),
            "Resource gene IDs are empty or duplicated after version stripping")
in_resource <- signature_ids %in% resource_ids
common_ids <- signature_ids[in_resource & signature_ids %in% rownames(corrected)]
signature_coverage <- length(common_ids) / length(signature_ids)
assert_true(length(common_ids) == EXPECTED_N_RESOURCE_SIGNATURE,
            paste0("Expected ", EXPECTED_N_RESOURCE_SIGNATURE,
                   " Resource-common signature genes; observed ", length(common_ids)))
assert_true(signature_coverage >= MIN_SIGNATURE_COVERAGE,
            paste0("Signature coverage below ", MIN_SIGNATURE_COVERAGE))

signature_pca <- prcomp(t(corrected[common_ids, , drop = FALSE]),
                        center = TRUE, scale. = FALSE)
orientation_signature <- sign(cor(
  signature_pca$x[, 1L], released_pc1_operational, method = "pearson"
))
assert_true(orientation_signature != 0,
            "Discovery 145-gene PC1 has zero correlation with the operational released PC1")
signature_pc1 <- unname(signature_pca$x[, 1L] * orientation_signature)

discovery_model <- list(
  signature_gene_ids_all = signature_ids,
  common_gene_ids = common_ids,
  discovery_center = setNames(unname(signature_pca$center), names(signature_pca$center)),
  discovery_loading_oriented = setNames(
    unname(signature_pca$rotation[, 1L] * orientation_signature),
    rownames(signature_pca$rotation)
  ),
  discovery_full_center = setNames(unname(full_pca$center), names(full_pca$center)),
  discovery_full_loading_oriented = setNames(
    unname(full_pca$rotation[, 1L] * orientation_full_operational),
    rownames(full_pca$rotation)
  ),
  provenance = list(
    method = "log2(1+counts); quantile normalization; ComBat(Dataset x Sex, mod=NAS); PCA",
    expected_commit = EXPECTED_COMMIT,
    observed_commit = observed_commit,
    raw_full_pc1_orientation_reference = "released data/PC1_sorted_samples.csv raw PC1",
    full_pc1_orientation_reference = paste0(
      "released operational continuum: -raw PC1, asserted from ",
      OPERATIONAL_SIGN_RELATIVE, "; no target histology used"
    ),
    signature_pc1_orientation_reference = paste0(
      "released operational full-transcriptome continuum; no target histology used"
    ),
    operational_pc1_multiplier = OPERATIONAL_PC1_MULTIPLIER,
    scale = FALSE,
    n_discovery_participants = ncol(corrected),
    n_filtered_discovery_genes = nrow(corrected),
    n_resource_common_signature_genes = length(common_ids)
  )
)

discovery_metadata <- data.frame(
  sample_id = sample_ids,
  dataset = as.character(metadata[[dataset_col]]),
  sex = as.character(metadata[[sex_col]]),
  NAS = as.numeric(metadata[[nas_col]]),
  fibrosis = as.numeric(metadata[[fibrosis_col]]),
  stringsAsFactors = FALSE
)
metadata_for_cor <- cbind(
  discovery_metadata,
  steatosis = as.numeric(metadata[[steatosis_col]]),
  ballooning = as.numeric(metadata[[ballooning_col]]),
  inflammation = as.numeric(metadata[[inflammation_col]])
)
participant_scores <- data.frame(
  sample_id = sample_ids,
  dataset = discovery_metadata$dataset,
  sex = discovery_metadata$sex,
  NAS = discovery_metadata$NAS,
  fibrosis = discovery_metadata$fibrosis,
  pc1_reproduced_raw = pc1_reproduced_raw,
  pc1_released_raw = released_pc1_raw,
  pc1_reproduced_operational = pc1_reproduced_operational,
  pc1_released_operational = released_pc1_operational,
  released_file_position = released_file_position,
  released_operational_rank_percentile = released_operational_rank_percentile,
  pc2_released = released_pc2,
  signature_pc1 = signature_pc1,
  pc1_rank_percentile = (rank(pc1_reproduced_operational, ties.method = "average") - 1) /
    (length(pc1_reproduced_operational) - 1),
  stringsAsFactors = FALSE
)
histology_correlations <- correlation_table(participant_scores, metadata_for_cor)

# ------------------------------------------------------------- code-path audit

log_step("auditing released downstream ordering code")
code_audit <- audit_code_paths(vendor_root)
n_downstream_pc1_lines <- sum(code_audit$evidence_type == "pc1_sorted_file" &
                                code_audit$downstream %in% TRUE, na.rm = TRUE)
n_operational_pc1_negations <- sum(
  code_audit$evidence_type == "operational_pc1_negation", na.rm = TRUE
)
n_operational_rf_inverse_pc1 <- sum(
  code_audit$evidence_type == "operational_rf_inverse_pc1", na.rm = TRUE
)
n_pseudotime_extractions <- sum(code_audit$evidence_type == "pseudotime_extraction")
assert_true(n_downstream_pc1_lines > 0L,
            "No downstream consumer of PC1_sorted_samples.csv was found")
assert_true(n_operational_pc1_negations >= 1L,
            "Released source does not negate raw PC1 for the operational continuum")
assert_true(n_operational_rf_inverse_pc1 >= 1L,
            "Released RF source no longer identifies inverse PC1 as the real trajectory")
assert_true(n_pseudotime_extractions == 0L,
            "Released downstream source extracts Slingshot pseudotime")

# ----------------------------------------------------------- summary + writes

summary <- do.call(rbind, list(
  summary_row("vendor_commit", observed_commit, paste0("equals ", EXPECTED_COMMIT),
              observed_commit == EXPECTED_COMMIT),
  summary_row("vendor_v1_tag_present", EXPECTED_TAG, paste0(EXPECTED_TAG, " at HEAD"),
              EXPECTED_TAG %in% observed_tags),
  summary_row("vendor_tracked_worktree_clean", length(tracked_status), "equals 0",
              length(tracked_status) == 0L),
  summary_row("seed", SEED, "prespecified", TRUE),
  summary_row("supplied_manifest_used", ifelse(is.na(manifest_used), "none", manifest_used),
              "verified when present", if (is.na(manifest_used)) NA else TRUE),
  summary_row("n_input_participants", nrow(metadata_raw), paste0("equals ", EXPECTED_N_INPUT),
              nrow(metadata_raw) == EXPECTED_N_INPUT),
  summary_row("removed_outlier", EXPECTED_OUTLIER, "released fifth participant",
              identical(outlier_index, 5L)),
  summary_row("n_discovery_participants", ncol(corrected),
              paste0("equals ", EXPECTED_N_DISCOVERY), ncol(corrected) == EXPECTED_N_DISCOVERY),
  summary_row("n_input_genes", nrow(counts), NA_character_, NA),
  summary_row("n_filtered_genes", nrow(corrected), "rowSums(counts) > 135", TRUE),
  summary_row("pc1_variance_percent", variance_table$variance_percent[[1L]], NA_character_, NA),
  summary_row("pc2_variance_percent", variance_table$variance_percent[[2L]], NA_character_, NA),
  summary_row("released_pc1_spearman", pc1_spearman,
              paste0("absolute value >= ", MIN_RELEASED_PC1_SPEARMAN),
              abs(pc1_spearman) >= MIN_RELEASED_PC1_SPEARMAN),
  summary_row("released_pc1_pearson", pc1_pearson, NA_character_, NA),
  summary_row("operational_pc1_multiplier", OPERATIONAL_PC1_MULTIPLIER,
              "equals -1 from frozen released source", OPERATIONAL_PC1_MULTIPLIER == -1,
              paste0(OPERATIONAL_SIGN_RELATIVE, ": ",
                     trimws(operational_sign_source[[operational_sign_hits[[1L]]]]))),
  summary_row("operational_pc1_released_spearman", operational_pc1_spearman,
              paste0("absolute value >= ", MIN_RELEASED_PC1_SPEARMAN),
              abs(operational_pc1_spearman) >= MIN_RELEASED_PC1_SPEARMAN),
  summary_row("released_pc1_file_strictly_descending",
              released_pc1_file_strictly_descending, "equals TRUE",
              released_pc1_file_strictly_descending,
              "Operational percentile increases with file position as raw PC1 decreases"),
  summary_row("operational_percentile_file_order_max_abs_difference",
              operational_percentile_file_order_max_abs_difference,
              "equals 0", operational_percentile_file_order_max_abs_difference == 0,
              "Progression percentile is rank(-raw released PC1)"),
  summary_row("released_corrected_matrix_pearson", corrected_matrix_pearson,
              NA_character_, NA),
  summary_row("released_corrected_matrix_max_abs_difference", corrected_matrix_max_abs_diff,
              NA_character_, NA),
  summary_row("n_signature_genes", length(signature_ids),
              paste0("equals ", EXPECTED_N_SIGNATURE), length(signature_ids) == EXPECTED_N_SIGNATURE),
  summary_row("n_resource_common_signature_genes", length(common_ids),
              paste0("equals ", EXPECTED_N_RESOURCE_SIGNATURE),
              length(common_ids) == EXPECTED_N_RESOURCE_SIGNATURE),
  summary_row("signature_resource_coverage", signature_coverage,
              paste0(">= ", MIN_SIGNATURE_COVERAGE), signature_coverage >= MIN_SIGNATURE_COVERAGE),
  summary_row("signature_pc1_operational_spearman",
              cor(signature_pc1, released_pc1_operational, method = "spearman"),
              NA_character_, NA),
  summary_row("downstream_pc1_consumer_lines", n_downstream_pc1_lines, "> 0",
              n_downstream_pc1_lines > 0L,
              "Released downstream scripts read PC1_sorted_samples.csv"),
  summary_row("operational_pc1_negation_lines", n_operational_pc1_negations, ">= 1",
              n_operational_pc1_negations >= 1L,
              "Released source explicitly applies PC1 <- -PC1 for progression-facing analysis"),
  summary_row("operational_rf_inverse_pc1_lines", n_operational_rf_inverse_pc1, ">= 1",
              n_operational_rf_inverse_pc1 >= 1L,
              "Released RF source identifies inverse PC1 as the real trajectory"),
  summary_row("downstream_pseudotime_extractions", n_pseudotime_extractions, "equals 0",
              n_pseudotime_extractions == 0L,
              "No slingPseudotime() or slingCurveWeights() call in released source")
))

signature_genes <- data.frame(
  gene_id_base = signature_ids,
  gene_symbol = signature_symbols,
  in_resource = in_resource,
  mapping_reference = "GENCODE_v49_Ensembl_base_ID",
  missing_reason = ifelse(in_resource, "present_in_resource_substrate",
                          "absent_from_resource_substrate"),
  stringsAsFactors = FALSE
)
missing_signature_genes <- signature_genes[!signature_genes$in_resource, , drop = FALSE]
assert_true(nrow(missing_signature_genes) == EXPECTED_N_SIGNATURE -
              EXPECTED_N_RESOURCE_SIGNATURE,
            "Missing-signature audit must contain exactly six genes")

input_manifest <- data.frame(
  input_role = c(names(required_paths), "resource_gene_ids", "supplied_vendor_manifest"),
  path = c(unname(required_paths), resource_gene_path,
           ifelse(is.na(manifest_used), NA_character_, manifest_used)),
  size_bytes = c(unname(file.size(required_paths)), file.size(resource_gene_path),
                 ifelse(is.na(manifest_used), NA_real_, file.size(manifest_used))),
  sha256 = c(vapply(unname(required_paths), sha256_file, character(1L)),
             sha256_file(resource_gene_path),
             ifelse(is.na(manifest_used), NA_character_, sha256_file(manifest_used))),
  stringsAsFactors = FALSE
)

assert_true(all(summary$pass[!is.na(summary$pass)]),
            "At least one reproduction assertion failed before output was written")
assert_true(identical(names(discovery_model)[1:6], c(
  "signature_gene_ids_all", "common_gene_ids", "discovery_center",
  "discovery_loading_oriented", "discovery_full_center",
  "discovery_full_loading_oriented"
)), "Discovery model interface drifted")

log_step("writing immutable reproduction artifacts")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write_tsv_once(participant_scores, file.path(out_dir, "participant_scores.tsv"))
write_tsv_once(histology_correlations, file.path(out_dir, "histology_correlations.tsv"))
write_tsv_once(summary, file.path(out_dir, "reproduction_summary.tsv"))
write_tsv_once(code_audit, file.path(out_dir, "code_path_audit.tsv"))
write_tsv_once(variance_table, file.path(out_dir, "pca_variance_explained.tsv"))
write_tsv_once(signature_genes, file.path(out_dir, "signature_genes.tsv"))
write_tsv_once(missing_signature_genes,
               file.path(out_dir, "missing_signature_genes.tsv"))
save_rds_once(discovery_model, file.path(out_dir, "discovery_signature_model.rds"))
save_rds_once(corrected, file.path(out_dir, "discovery_normalized_expression.rds"))
write_tsv_once(discovery_metadata, file.path(out_dir, "discovery_metadata.tsv"))
write_tsv_once(input_manifest, file.path(out_dir, "input_checksums.tsv"))
write_lines_once(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

log_step("REPRODUCTION_COMPLETE: ", out_dir)
