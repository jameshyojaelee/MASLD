# Assay-independent chronic-state scoring for the sealed Plan 43 candidate.

candidate_id <- "chronic-state-risk-bridge-2026-08-09"
`%||%` <- function(x, y) if (is.null(x) || !length(x) || is.na(x[1])) y else x

script_path <- tryCatch(normalizePath(sys.frame(1)$ofile, mustWork = TRUE), error = function(e) NA_character_)
if (is.na(script_path)) {
  all_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- all_args[grep("^--file=", all_args)][1]
  script_path <- normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
}
project_root <- normalizePath(file.path(dirname(script_path), "../../../.."), mustWork = TRUE)
candidate_root <- file.path(project_root, "Analysis/Multimodal_Program_Projection/candidates", candidate_id)

read_tsv <- function(path) {
  read.delim(path, sep = "\t", quote = "", check.names = FALSE,
             stringsAsFactors = FALSE, comment.char = "")
}

write_tsv <- function(x, path) {
  expected_root <- paste0(normalizePath(candidate_root, mustWork = TRUE), .Platform$file.sep)
  parent <- normalizePath(dirname(path), mustWork = FALSE)
  if (!startsWith(paste0(parent, .Platform$file.sep), expected_root)) {
    stop("Refusing Plan 43 write outside candidate root: ", path)
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path))
  write.table(x, temporary, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  if (!file.rename(temporary, path)) stop("Atomic rename failed: ", path)
}

require_plan43_seal <- function() {
  required <- file.path(candidate_root, c("SEALED.json", "SEAL_VALIDATED.json"))
  if (!all(file.exists(required))) stop("Plan 43 outcome access requires a validated seal")
  invisible(TRUE)
}

load_plan43_frozen <- function() {
  require_plan43_seal()
  membership <- read_tsv(file.path(candidate_root, "frozen_inputs/program_membership__program_membership_v2.tsv"))
  registry <- read_tsv(file.path(candidate_root, "frozen_inputs/program_registry__program_registry_v2.tsv"))
  axis <- read_tsv(file.path(candidate_root, "frozen_state_axis.tsv"))
  classes <- read_tsv(file.path(candidate_root, "integrity_correction_v2/corrected_evidence_classes_v2.tsv"))
  if (nrow(registry) != 117L || nrow(axis) != 117L) stop("Frozen 117-program geometry drift")
  if (!setequal(registry$program_uid, axis$program_uid)) stop("Registry/axis UID mismatch")
  list(membership = membership, registry = registry, axis = axis, classes = classes)
}

collapse_expression_to_symbol <- function(expression, symbols, method = c("highest_reference_median", "median"), reference_columns = seq_len(ncol(expression))) {
  method <- match.arg(method)
  stopifnot(nrow(expression) == length(symbols))
  valid <- !is.na(symbols) & nzchar(symbols)
  expression <- expression[valid, , drop = FALSE]
  symbols <- symbols[valid]
  if (method == "median") {
    output <- t(vapply(unique(symbols), function(symbol) {
      apply(expression[symbols == symbol, , drop = FALSE], 2L, median, na.rm = TRUE)
    }, numeric(ncol(expression))))
    rownames(output) <- unique(symbols)
    colnames(output) <- colnames(expression)
    return(output)
  }
  medians <- apply(expression[, reference_columns, drop = FALSE], 1L, median, na.rm = TRUE)
  ordering <- order(symbols, -medians, seq_along(symbols))
  selected <- ordering[!duplicated(symbols[ordering])]
  output <- expression[selected, , drop = FALSE]
  rownames(output) <- symbols[selected]
  output
}

baseline_z <- function(expression, reference_samples, minimum_sd = 1e-8) {
  reference_index <- match(reference_samples, colnames(expression))
  if (anyNA(reference_index) || anyDuplicated(reference_index)) {
    stop("Reference samples must be unique expression columns")
  }
  reference <- expression[, reference_index, drop = FALSE]
  means <- rowMeans(reference, na.rm = TRUE)
  sds <- apply(reference, 1L, sd, na.rm = TRUE)
  keep <- is.finite(means) & is.finite(sds) & sds > minimum_sd
  z <- sweep(expression[keep, , drop = FALSE], 1L, means[keep], "-")
  z <- sweep(z, 1L, sds[keep], "/")
  z[!is.finite(z)] <- NA_real_
  list(z = z, mean = means[keep], sd = sds[keep], n_reference = length(reference_index))
}

weighted_jaccard <- function(left, right) {
  genes <- union(names(left), names(right))
  l <- left[match(genes, names(left))]; l[is.na(l)] <- 0
  r <- right[match(genes, names(right))]; r[is.na(r)] <- 0
  denominator <- sum(pmax(l, r))
  if (denominator == 0) 0 else sum(pmin(l, r)) / denominator
}

membership_families <- function(membership, threshold = 0.20) {
  uids <- unique(membership$program_uid)
  vectors <- lapply(uids, function(uid) {
    rows <- membership[membership$program_uid == uid & nzchar(membership$mapped_symbol), , drop = FALSE]
    values <- tapply(rows$original_l1_weight, rows$mapped_symbol, sum)
    values / sum(values)
  })
  names(vectors) <- uids
  adjacency <- setNames(lapply(uids, function(x) character()), uids)
  if (length(uids) > 1L) {
    for (i in seq_len(length(uids) - 1L)) for (j in (i + 1L):length(uids)) {
      if (weighted_jaccard(vectors[[i]], vectors[[j]]) >= threshold) {
        adjacency[[i]] <- c(adjacency[[i]], uids[j])
        adjacency[[j]] <- c(adjacency[[j]], uids[i])
      }
    }
  }
  family <- setNames(rep(NA_integer_, length(uids)), uids)
  next_family <- 0L
  for (uid in uids) {
    if (!is.na(family[uid])) next
    next_family <- next_family + 1L
    queue <- uid
    family[uid] <- next_family
    while (length(queue)) {
      current <- queue[1]; queue <- queue[-1]
      unseen <- adjacency[[current]][is.na(family[adjacency[[current]]])]
      if (length(unseen)) {
        family[unseen] <- next_family
        queue <- c(queue, unseen)
      }
    }
  }
  data.frame(program_uid = names(family), overlap_family = unname(family), stringsAsFactors = FALSE)
}

program_scores <- function(expression, reference_samples, membership,
                           gene_scheme = c("weighted", "equal", "leave_top"),
                           minimum_genes = 8L, minimum_weight = 0.20) {
  gene_scheme <- match.arg(gene_scheme)
  standardization <- baseline_z(expression, reference_samples)
  z <- standardization$z
  uids <- unique(membership$program_uid)
  scores <- matrix(NA_real_, nrow = length(uids), ncol = ncol(z),
                   dimnames = list(uids, colnames(z)))
  audit <- vector("list", length(uids))
  for (i in seq_along(uids)) {
    uid <- uids[i]
    members <- membership[membership$program_uid == uid, , drop = FALSE]
    members <- members[nzchar(members$mapped_symbol), , drop = FALSE]
    members <- aggregate(original_l1_weight ~ mapped_symbol, members, sum)
    retained <- members[members$mapped_symbol %in% rownames(z), , drop = FALSE]
    retained_mass <- sum(retained$original_l1_weight)
    testable <- nrow(retained) >= minimum_genes && retained_mass >= minimum_weight
    if (testable) {
      if (gene_scheme == "leave_top") retained <- retained[-which.max(retained$original_l1_weight), , drop = FALSE]
      weights <- if (gene_scheme == "equal") rep(1, nrow(retained)) else retained$original_l1_weight
      weights <- weights / sum(weights)
      scores[i, ] <- as.numeric(crossprod(weights, z[retained$mapped_symbol, , drop = FALSE]))
    }
    audit[[i]] <- data.frame(
      program_uid = uid, gene_scheme = gene_scheme,
      n_mapped_genes = nrow(retained), retained_l1_weight = retained_mass,
      testable = testable,
      failure_reason = if (testable) "" else if (nrow(retained) < minimum_genes) "fewer_than_8_mapped_genes" else "less_than_20pct_original_l1_weight",
      stringsAsFactors = FALSE
    )
  }
  list(scores = scores, audit = do.call(rbind, audit), standardization = standardization)
}

state_score <- function(program_score_matrix, axis, loading_scheme = c("primary", "sign", "lineage_balanced", "overlap_balanced"),
                        membership = NULL, drop_lineage = NULL, drop_program = NULL,
                        minimum_programs = 90L, minimum_axis_mass = 0.80) {
  loading_scheme <- match.arg(loading_scheme)
  axis <- axis[match(rownames(program_score_matrix), axis$program_uid), , drop = FALSE]
  if (anyNA(axis$program_uid)) stop("Program score/axis mismatch")
  raw <- switch(
    loading_scheme,
    primary = axis$primary_loading,
    sign = axis$sign_loading,
    lineage_balanced = axis$lineage_balanced_loading,
    overlap_balanced = {
      if (is.null(membership)) stop("overlap_balanced requires membership")
      families <- membership_families(membership)
      family <- families$overlap_family[match(axis$program_uid, families$program_uid)]
      base <- axis$primary_loading
      family_mass <- tapply(abs(base), family, sum)
      base / family_mass[as.character(family)] / length(family_mass)
    }
  )
  eligible <- apply(program_score_matrix, 1L, function(x) all(is.finite(x)))
  if (!is.null(drop_lineage)) eligible <- eligible & axis$cell_type != drop_lineage
  if (!is.null(drop_program)) eligible <- eligible & axis$program_uid != drop_program
  n_testable <- sum(eligible)
  primary_mass <- sum(abs(axis$primary_loading[eligible]))
  passes <- n_testable >= minimum_programs && primary_mass >= minimum_axis_mass
  scores <- rep(NA_real_, ncol(program_score_matrix)); names(scores) <- colnames(program_score_matrix)
  if (passes) {
    loadings <- raw[eligible]
    loadings <- loadings / sum(abs(loadings))
    scores <- as.numeric(crossprod(loadings, program_score_matrix[eligible, , drop = FALSE]))
    names(scores) <- colnames(program_score_matrix)
  }
  list(
    scores = scores,
    audit = data.frame(
      loading_scheme = loading_scheme, dropped_lineage = drop_lineage %||% "",
      dropped_program = drop_program %||% "", n_testable_programs = n_testable,
      retained_primary_axis_mass = primary_mass, dataset_gate_pass = passes,
      stringsAsFactors = FALSE
    )
  )
}

score_state_schemes <- function(expression, reference_samples, frozen) {
  output <- list(); audits <- list(); standardization <- NULL
  for (gene_scheme in c("weighted", "equal", "leave_top")) {
    programs <- program_scores(expression, reference_samples, frozen$membership, gene_scheme)
    if (is.null(standardization)) standardization <- programs$standardization
    loadings <- if (gene_scheme == "weighted") c("primary", "sign", "lineage_balanced", "overlap_balanced") else "primary"
    for (loading in loadings) {
      name <- paste(gene_scheme, loading, sep = "__")
      state <- state_score(programs$scores, frozen$axis, loading, frozen$membership)
      output[[name]] <- state$scores
      audits[[name]] <- cbind(gene_scheme = gene_scheme, state$audit)
    }
  }
  primary_programs <- program_scores(expression, reference_samples, frozen$membership, "weighted")
  highest <- frozen$axis$program_uid[which.max(abs(frozen$axis$primary_loading))]
  dropped <- state_score(primary_programs$scores, frozen$axis, "primary", frozen$membership, drop_program = highest)
  output[["weighted__leave_highest_loading_program"]] <- dropped$scores
  audits[["weighted__leave_highest_loading_program"]] <- cbind(gene_scheme = "weighted", dropped$audit)
  for (lineage in sort(unique(frozen$axis$cell_type))) {
    name <- paste0("weighted__leave_lineage_", lineage)
    state <- state_score(primary_programs$scores, frozen$axis, "primary", frozen$membership, drop_lineage = lineage)
    output[[name]] <- state$scores
    audits[[name]] <- cbind(gene_scheme = "weighted", state$audit)
  }
  list(
    scores = do.call(rbind, output),
    audit = do.call(rbind, audits),
    program_scores = primary_programs$scores,
    program_testability = primary_programs$audit,
    standardization = standardization
  )
}
