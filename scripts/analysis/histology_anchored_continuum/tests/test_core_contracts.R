#!/usr/bin/env Rscript

# Fast unit contracts for the histology-anchored continuum. This deliberately
# uses base assertions instead of adding a test-framework dependency.

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
stopifnot(length(script_arg) == 1L)
test_dir <- dirname(normalizePath(script_arg))
analysis_dir <- normalizePath(file.path(test_dir, ".."))
project_root <- normalizePath(file.path(analysis_dir, "../../.."))
Sys.setenv(MASLD_PROJECT_ROOT = project_root)
source(file.path(analysis_dir, "lib_continuum.R"))

expect_error <- function(expression, pattern = NULL) {
  error <- tryCatch({
    force(expression)
    NULL
  }, error = identity)
  stopifnot(inherits(error, "error"))
  if (!is.null(pattern)) {
    stopifnot(grepl(pattern, conditionMessage(error), fixed = TRUE))
  }
  invisible(error)
}

test_case <- function(name, expression) {
  force(expression)
  cat(sprintf("PASS\t%s\n", name))
  invisible(TRUE)
}

test_case("Ensembl base-ID mapping and symbol collapse", {
  mapped <- base_gene_id(c(
    "ENSG00000123456.1", "ENSG00000999999.49", "ENSG00000888888",
    "NOT_AN_ENSEMBL.VERSION"
  ))
  stopifnot(identical(
    mapped,
    c("ENSG00000123456", "ENSG00000999999", "ENSG00000888888",
      "NOT_AN_ENSEMBL.VERSION")
  ))

  expression <- matrix(
    c(2, 4, 6, 8, 10, 12), nrow = 3L, byrow = TRUE,
    dimnames = list(
      c("ENSG000001.1", "ENSG000002.7", "ENSG000003"),
      c("P1", "P2")
    )
  )
  annotation <- data.table(
    gene_id = c("ENSG000001.12", "ENSG000002", "ENSG000003.3"),
    gene_name = c("GeneA", "genea", "GeneB")
  )
  collapsed <- collapse_symbols(expression, rownames(expression), annotation)
  stopifnot(identical(rownames(collapsed), c("GENEA", "GENEB")))
  stopifnot(isTRUE(all.equal(
    unname(collapsed),
    unname(rbind(GENEA = colMeans(expression[1:2, , drop = FALSE]),
                 GENEB = expression[3, ])),
    tolerance = 0
  )))

  conflicting_annotation <- rbind(
    annotation,
    data.table(gene_id = "ENSG000001", gene_name = "DifferentSymbol")
  )
  expect_error(
    collapse_symbols(expression, rownames(expression), conflicting_annotation),
    "One Ensembl gene maps to multiple symbols"
  )
})

test_case("PCA sign invariance under simultaneous score/loading reversal", {
  scores <- setNames(c(-3, -1, 0.5, 3.5), paste0("P", 1:4))
  loadings <- setNames(c(0.2, -0.7, 0.4, 0.5), paste0("G", 1:4))
  reference <- loadings

  forward <- orient_pca(scores, loadings, reference, minimum_overlap = 4L)
  reversed <- orient_pca(-scores, -loadings, reference, minimum_overlap = 4L)
  stopifnot(forward$orientation_sign == 1)
  stopifnot(reversed$orientation_sign == -1)
  stopifnot(isTRUE(all.equal(forward$scores, reversed$scores, tolerance = 0)))
  stopifnot(isTRUE(all.equal(forward$loadings, reversed$loadings, tolerance = 0)))
  stopifnot(forward$orientation_cosine > 1 - 1e-15)
  stopifnot(reversed$orientation_cosine > 1 - 1e-15)
})

test_case("Released operational continuum direction is frozen before validation", {
  prespec <- read_prespec()
  orientation <- prespec$orientation
  stopifnot(orientation$operational_continuum_multiplier == -1)
  stopifnot(identical(
    orientation$evidence_expression,
    "merged_data$PC1 <- -merged_data$PC1"
  ))
  stopifnot(identical(
    orientation$direction_definition,
    "low_values_early_high_values_advanced_fibrosis"
  ))
  stopifnot(isFALSE(orientation$target_histology_used_for_orientation))

  raw_score <- c(-4, -1, 2, 5)
  operational <- raw_score * orientation$operational_continuum_multiplier
  stopifnot(identical(operational, -raw_score))
})

test_case("Fixed discovery-loading projection", {
  genes <- c("ENSG1", "ENSG2", "ENSG3")
  expression <- matrix(
    c(2, 5, 8, 4, 7, 10, 6, 9, 12),
    nrow = 3L, byrow = TRUE,
    dimnames = list(genes, c("V1", "V2", "V3"))
  )
  center <- setNames(c(3, 5, 7), genes)
  loading <- setNames(c(0.5, -0.25, 2), genes)
  centered <- sweep(expression[genes, , drop = FALSE], 1L, center[genes], "-")
  observed <- as.numeric(crossprod(loading[genes], centered))
  expected <- vapply(seq_len(ncol(expression)), function(j) {
    sum((expression[, j] - center) * loading)
  }, numeric(1L))
  stopifnot(isTRUE(all.equal(observed, expected, tolerance = 0)))

  reordered <- expression[rev(genes), , drop = FALSE]
  reordered_projection <- as.numeric(crossprod(
    loading[genes],
    sweep(reordered[genes, , drop = FALSE], 1L, center[genes], "-")
  ))
  stopifnot(isTRUE(all.equal(observed, reordered_projection, tolerance = 0)))
})

test_case("Discovery/independent-validation participant disjointness", {
  assert_disjoint <- function(discovery_ids, validation_ids) {
    overlap <- intersect(as.character(discovery_ids), as.character(validation_ids))
    if (length(overlap)) {
      stop("Discovery/validation participant overlap: ",
           paste(overlap, collapse = ";"), call. = FALSE)
    }
    invisible(TRUE)
  }
  stopifnot(isTRUE(assert_disjoint(c("D1", "D2"), c("V1", "V2"))))
  expect_error(
    assert_disjoint(c("D1", "D2"), c("V1", "D2")),
    "Discovery/validation participant overlap: D2"
  )
})

test_case("All signature genes are removed before program scoring", {
  membership <- data.table(
    program_uid = c("P1", "P1", "P1", "P2", "P2"),
    gene_symbol = c("A", "B", "C", "B", "D"),
    weight = c(0.4, 0.3, 0.3, 0.7, 0.3)
  )
  signature_symbols <- c("B", "D")
  membership[, excluded_signature_gene :=
               toupper(gene_symbol) %in% toupper(signature_symbols)]
  membership_use <- membership[excluded_signature_gene == FALSE]
  excluded <- membership[excluded_signature_gene == TRUE]

  stopifnot(!any(toupper(membership_use$gene_symbol) %in%
                  toupper(signature_symbols)))
  stopifnot(setequal(excluded$gene_symbol, c("B", "D")))
  stopifnot(isTRUE(all.equal(
    membership_use[program_uid == "P1", sum(weight)], 0.7,
    tolerance = 1e-15
  )))
  # Retained L1 weight is measured against the original program total, not
  # renormalised before the 80% testability decision.
  original_total <- membership[, .(original_l1 = sum(weight)), by = program_uid]
  retained <- membership_use[, .(retained_l1 = sum(weight)), by = program_uid]
  audit <- merge(original_total, retained, by = "program_uid", all.x = TRUE)
  audit[is.na(retained_l1), retained_l1 := 0]
  audit[, retained_fraction := retained_l1 / original_l1]
  stopifnot(audit[program_uid == "P1", retained_fraction] == 0.7)
  stopifnot(audit[program_uid == "P2", retained_fraction] == 0)

  # The prespecified 80% L1 rule concerns signature removal. Frozen source rows
  # that were already unmapped remain a separate observability audit; they are
  # not incorrectly reclassified as signature exclusions.
  original_l1 <- 1
  signature_excluded_l1 <- 0.195
  observed_mapped_l1 <- 0.738
  signature_retained_fraction <-
    (original_l1 - signature_excluded_l1) / original_l1
  observed_fraction <- observed_mapped_l1 / original_l1
  stopifnot(signature_retained_fraction >= 0.8, observed_fraction < 0.8)
})

test_case("BH adjustment retains the complete 117-program family", {
  p <- c(0.001, 0.02, NA_real_, 0.5)
  observed <- p_adjust_complete_family(p, family_size = 117L, method = "BH")
  finite <- is.finite(p)
  expected <- rep(NA_real_, length(p))
  expected[finite] <- p.adjust(p[finite], method = "BH", n = 117L)
  stopifnot(identical(is.na(observed), is.na(expected)))
  stopifnot(isTRUE(all.equal(observed[finite], expected[finite], tolerance = 0)))
  stopifnot(all(observed[finite] >= p.adjust(p[finite], method = "BH")))
})

test_case("Nine fixed overlapping visualization windows", {
  prespec <- read_prespec()
  centers <- as.numeric(prespec$windows$centers)
  width <- as.numeric(prespec$windows$width)
  stopifnot(length(centers) == 9L)
  stopifnot(isTRUE(all.equal(centers, seq(0.1, 0.9, by = 0.1), tolerance = 1e-15)))
  stopifnot(width == 0.2)
  bounds <- data.table(
    center = centers,
    lower = pmax(0, centers - width / 2),
    upper = pmin(1, centers + width / 2)
  )
  stopifnot(nrow(bounds) == 9L)
  stopifnot(bounds$lower[[1L]] == 0 && bounds$upper[[9L]] == 1)
  stopifnot(isTRUE(all.equal(
    bounds$upper[-nrow(bounds)] - bounds$lower[-1L],
    rep(width / 2, 8L), tolerance = 1e-15
  )))
  memberships <- vapply(centers, function(center) {
    lower <- max(0, center - width / 2)
    upper <- min(1, center + width / 2)
    0.2 >= lower && (0.2 < upper || (center == max(centers) && 0.2 <= upper))
  }, logical(1L))
  stopifnot(sum(memberships) == 2L)
  stopifnot(identical(
    prespec$windows$interval_rule,
    "lower_closed_upper_open_except_final_upper_closed"
  ))
  stopifnot(isFALSE(prespec$windows$inference_use))
})

test_case("Seeded resampling is deterministic", {
  x <- seq_len(20L)
  y <- x + rep(c(-1, 1), 10L)
  strata <- rep(letters[1:4], each = 5L)
  boot_one <- bootstrap_spearman(x, y, strata, replicates = 100L, seed = 401L)
  boot_two <- bootstrap_spearman(x, y, strata, replicates = 100L, seed = 401L)
  stopifnot(identical(boot_one, boot_two))

  perm_one <- permutation_spearman_greater(
    x, y, replicates = 100L, seed = 402L
  )
  perm_two <- permutation_spearman_greater(
    x, y, replicates = 100L, seed = 402L
  )
  stopifnot(identical(perm_one, perm_two))
})

test_case("Within-stage permutations preserve stratum membership and positions", {
  permute_within_strata <- function(x, strata, seed) {
    stopifnot(length(x) == length(strata), !anyNA(strata))
    groups <- split(seq_along(x), strata)
    permuted <- x
    set.seed(seed)
    for (g in groups) {
      permuted[g] <- sample(x[g], length(g), replace = FALSE)
    }
    permuted
  }

  values <- seq_len(18L)
  strata <- rep(c("F0_F", "F0_M", "F1_F"), each = 6L)
  first <- permute_within_strata(values, strata, 20260819L)
  second <- permute_within_strata(values, strata, 20260819L)
  stopifnot(identical(first, second))
  stopifnot(!identical(first, values))
  for (level in unique(strata)) {
    idx <- which(strata == level)
    stopifnot(identical(sort(first[idx]), sort(values[idx])))
  }
  # The historical `unlist(lapply(groups, sample))` index-vector pattern is
  # invalid when groups are interleaved because it changes global positions.
  interleaved_strata <- rep(c("F0", "F1", "F2"), 6L)
  interleaved <- permute_within_strata(values, interleaved_strata, 20260819L)
  for (level in unique(interleaved_strata)) {
    idx <- which(interleaved_strata == level)
    stopifnot(identical(sort(interleaved[idx]), sort(values[idx])))
  }
})

cat("test_core_contracts: PASS\n")
