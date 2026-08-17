#!/usr/bin/env Rscript
# Shared donor-collapse utility for the single-cell pseudobulk pipelines.
#
# FIXES a pseudoreplication bug: results_gpu_v2/pseudobulk/{CellType}_pseudobulk.csv
# matrices and the sample_meta tables built from the integrated atlas h5ad are keyed
# on `obs/sample`, which is a SEQUENCING RUN (or a sort-fraction library), not a
# biological donor, for 4 datasets:
#   GSE244832  117 runs -> 18 donors  (data/GSE244832/metadata/donor_pairing.csv)  [SRR run multiplicity]
#   GSE202379   67 runs -> 46 donors  (data/GSE202379/metadata/donor_pairing.csv)  [SRR anatomical-zone multiplicity]
#   GSE185477   21 runs ->  3 donors  (data/GSE185477/metadata/donor_pairing.csv)  [SRR run multiplicity]
#   GSE136103   20 GSMs -> 10 donors  (data/GSE136103/metadata/donor_pairing.csv)  [GSM CD45-sort-fraction multiplicity]
# NOTE GSE136103 is keyed on GSM (not SRR) in the atlas, and its multiplicity is by
# CD45+/CD45- sort fraction + A/B replicate (e.g. GSM4041150/151/152 = healthy1
# cd45+/cd45-A/cd45-B). The `rna_srrs` column therefore holds GSM IDs for this
# dataset; the collapse is by string match on the atlas sample ID, so SRR vs GSM is
# transparent. Summing a patient's per-cell-type counts across sort fractions is
# correct: sort gate changes WHICH cells are captured, not per-cell expression.
# GSE174748 and the two analyzed GSE189600 libraries pass through 1:1. Liver Atlas
# is resolved separately. Its source authority contains 48 human runs from 38
# assay samples and 19 biological donors; the canonical atlas contains 42 of
# those runs. Two low-cell GSE189600 libraries remain descriptive-only. The resulting
# rosters are 273 libraries -> 102 analyzed donors and 275 -> 104 descriptive donors.
# ⚠ LATENT RISK (GSE189600): the raw ENA has genuine run-PAIRS per replicate (e.g.
# SRR22677586 & SRR22677587 are both "human_healthy_rep1"); the current atlas ingests
# only ONE run per replicate (587/589/593/595), so 1 run/donor pass-through is benign
# NOW. A future re-integration that pulled BOTH runs of a replicate WOULD silently
# pseudoreplicate GSE189600 — add a data/GSE189600/metadata/donor_pairing.csv mapping
# the run-pairs to donor before any such re-integration. (Whether healthy_rep1/2/3 are
# biological vs technical replicates is unconfirmed from labels alone; the census
# treats the 4 ingested = 4 distinct donors, which is a ±2 edge case on a tiny dataset
# with no headline impact.)
# Donor IDs are dataset-prefixed (e.g. "GSE244832_D01", "GSE136103_healthy1") to
# avoid collisions across datasets.
#
# Usage (see patched pseudobulk_de.R for the canonical integration pattern):
#   source(file.path(BASE, "Analysis/SingleCell/scripts/lib_donor_collapse.R"))
#   sample_meta <- collapse_sample_meta_to_donor(sample_meta, BASE)   # once, after loading
#   counts_mat  <- collapse_counts_to_donor(counts_mat, BASE)          # per cell type, after loading

suppressPackageStartupMessages(library(data.table))

.DONOR_PAIRING_FILES <- c(
  GSE244832 = "data/GSE244832/metadata/donor_pairing.csv",
  GSE185477 = "data/GSE185477/metadata/donor_pairing.csv",
  GSE202379 = "data/GSE202379/metadata/donor_pairing.csv",
  GSE136103 = "data/GSE136103/metadata/donor_pairing.csv"
)

# Build the srr -> dataset-prefixed-donor-id map once. Returns a named character
# vector: names = raw SRR run IDs, values = "{dataset}_{donor_id}".
build_srr_to_donor_map <- function(base_dir) {
  pairs <- list()
  for (ds in names(.DONOR_PAIRING_FILES)) {
    fp <- file.path(base_dir, .DONOR_PAIRING_FILES[[ds]])
    if (!file.exists(fp)) {
      warning("donor_pairing.csv not found for ", ds, " at ", fp, " -- skipping collapse for this dataset")
      next
    }
    dp <- fread(fp, colClasses = "character")
    for (i in seq_len(nrow(dp))) {
      srrs <- strsplit(dp$rna_srrs[i], ";", fixed = TRUE)[[1]]
      donor_id <- paste0(ds, "_", dp$donor_id[i])
      pairs[[length(pairs) + 1]] <- setNames(rep(donor_id, length(srrs)), srrs)
    }
  }
  liver_runs_path <- file.path(base_dir, "Liver_Atlas/metadata/SraRunTable.csv")
  liver_samples_path <- file.path(
    base_dir, "Liver_Atlas/metadata/GSE192740_sampleInfo_scRNAseq.tsv")
  if (!file.exists(liver_runs_path) || !file.exists(liver_samples_path)) {
    warning("Liver Atlas run/sample authorities are missing -- refusing to infer run-level donors")
  } else {
    liver_runs <- fread(liver_runs_path, colClasses = "character")
    liver_samples <- fread(liver_samples_path, sep = "\t", colClasses = "character")
    short_col <- "characteristics: shortFileName"
    if (!all(c("Run", "shortfilename") %in% names(liver_runs)) ||
        !all(c(short_col, "title") %in% names(liver_samples))) {
      stop("Liver Atlas donor authorities lack required columns", call. = FALSE)
    }
    liver_join <- merge(
      liver_runs[, .(Run, shortfilename)],
      liver_samples[, .(shortfilename = get(short_col), title)],
      by = "shortfilename", all = FALSE)
    liver_join[, donor_short := ifelse(
      grepl("H[0-9]+", title),
      sub(".*?(H[0-9]+).*", "\\1", title, perl = TRUE),
      NA_character_)]
    liver_join <- liver_join[!is.na(donor_short) & nzchar(donor_short)]
    if (uniqueN(liver_join$Run) != 48L || uniqueN(liver_join$donor_short) != 19L) {
      stop("Liver Atlas source run-to-donor cardinality drift", call. = FALSE)
    }
    pairs[[length(pairs) + 1]] <- setNames(
      paste0("Liver_Atlas_", liver_join$donor_short), liver_join$Run)
  }
  if (length(pairs) == 0) return(character(0))
  result <- do.call(c, pairs)
  if (anyDuplicated(names(result))) {
    stop("A sequencing library maps to more than one biological donor", call. = FALSE)
  }
  result
}

#' Collapse a run-level sample_meta table (one row per SRR run) to donor-level
#' (one row per biological donor). Samples not covered by any donor_pairing.csv
#' pass through unchanged (their `sample` value becomes their own donor id).
#'
#' @param sample_meta data.table with at least a `sample` column (SRR run ID or
#'   already-a-donor ID for 1:1 datasets), plus any other per-sample columns
#'   (e.g. condition, dataset, preparation_method) that should be carried over.
#' @param base_dir project root
#' @return data.table, one row per donor. `sample` is renamed to the donor ID.
#'   Non-`sample` columns are taken by MAJORITY VOTE across the donor's runs
#'   (matches the convention already used to build sample_meta from cell_meta).
#'   A warning is raised if a donor's runs disagree on `condition` (data problem).
collapse_sample_meta_to_donor <- function(sample_meta, base_dir) {
  srr_to_donor <- build_srr_to_donor_map(base_dir)
  sample_meta <- copy(sample_meta)
  sample_meta[, .donor_id := ifelse(sample %in% names(srr_to_donor),
                                     srr_to_donor[sample], sample)]

  n_runs_before <- nrow(sample_meta)
  n_donors_after <- uniqueN(sample_meta$.donor_id)
  if (n_runs_before != n_donors_after) {
    message("  [donor-collapse] sample_meta: ", n_runs_before, " runs -> ",
            n_donors_after, " donors")
  }

  # Check for within-donor condition disagreement (would indicate a data issue,
  # e.g. a donor's runs spanning two different disease labels).
  other_cols <- setdiff(names(sample_meta), c("sample", ".donor_id"))
  if ("condition" %in% other_cols) {
    n_cond_per_donor <- sample_meta[, uniqueN(condition), by = .donor_id]
    bad <- n_cond_per_donor[V1 > 1]
    if (nrow(bad) > 0) {
      warning("collapse_sample_meta_to_donor: ", nrow(bad), " donor(s) have ",
              "inconsistent `condition` across runs: ",
              paste(bad$.donor_id, collapse = ", "),
              " -- using majority vote, but this should be investigated")
    }
  }

  # Per-donor aggregation: `dataset` is invariant within a donor (take first);
  # every other column is resolved by MAJORITY VOTE across the donor's runs
  # (matches the convention already used to build sample_meta from cell_meta).
  collapsed <- sample_meta[, {
    res <- lapply(other_cols, function(cn) {
      v <- .SD[[cn]]
      if (cn == "dataset") {
        v[1]
      } else {
        # Majority vote over NON-NA values. `table()` drops NA by default, so an
        # ALL-NA column yields a length-0 table and `names(sort(...))[1]` returns
        # NULL — which data.table rejects ("Column N of j's result for the first
        # group is NULL"). Return NA_character_ explicitly in that case.
        tb <- sort(table(v), decreasing = TRUE)
        if (length(tb) == 0L) NA_character_ else names(tb)[1]
      }
    })
    setNames(res, other_cols)
  }, by = .donor_id, .SDcols = other_cols]
  setnames(collapsed, ".donor_id", "sample")
  collapsed[]
}

#' Collapse a genes x samples raw-count matrix (columns = SRR run IDs) to
#' genes x donors, summing raw counts across a donor's runs. Columns not
#' covered by any donor_pairing.csv pass through unchanged.
#'
#' @param counts_mat numeric matrix, rownames = genes, colnames = sample IDs
#' @param base_dir project root
#' @return numeric matrix, genes x donors (colnames = dataset-prefixed donor IDs
#'   for collapsed columns; original sample ID otherwise)
collapse_counts_to_donor <- function(counts_mat, base_dir) {
  srr_to_donor <- build_srr_to_donor_map(base_dir)
  samples <- colnames(counts_mat)
  donor_ids <- ifelse(samples %in% names(srr_to_donor),
                       srr_to_donor[samples], samples)
  if (length(unique(donor_ids)) == length(samples)) {
    colnames(counts_mat) <- donor_ids  # no collapsing needed, just relabel
    return(counts_mat)
  }
  # Sum columns within each donor group. t() so rowsum() groups the (former)
  # columns; t() back to genes x donors.
  collapsed_t <- rowsum(t(counts_mat), group = donor_ids, reorder = FALSE)
  out <- t(collapsed_t)
  storage.mode(out) <- storage.mode(counts_mat)
  out
}
