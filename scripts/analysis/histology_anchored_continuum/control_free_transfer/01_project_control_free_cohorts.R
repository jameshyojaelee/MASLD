#!/usr/bin/env Rscript
# CFT-v1: project the FROZEN 139-gene discovery loading into never-scored,
# control-free cohorts and anchor it to their own recorded histology.
# No refit. No per-cohort PCA. No re-derivation of the frozen signature.
suppressPackageStartupMessages({
  library(digest); library(data.table); library(limma); library(jsonlite)
})

OUT <- Sys.getenv("CFT_OUT_DIR"); stopifnot(nzchar(OUT))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

ROOT       <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
RELEASE    <- file.path(ROOT, "RNA-seq/results/histology_anchored_continuum/candidates/hac-continuum-20260818T024923Z")
PRESPEC    <- file.path(ROOT, "RNA-seq/results/histology_anchored_continuum/prespec/mnc-fib-v1-20260828T133158Z")
DGE        <- file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BG001-DECISION/arms/F_five/results/integration/merged_dge.rds")
COHORT_RES <- file.path(ROOT, "RNA-seq/Human/Patient_Cohorts/results")
META_ROOT  <- file.path(ROOT, "RNA-seq/Human/Patient_Cohorts/pipelines/custom")
PROP_ROOT  <- file.path(ROOT, "Analysis/Deconvolution/results")

B_DRAWS   <- 1000L
BOOT_B    <- 10000L
BOOT_SEED <- 20260828L
CV_SEED   <- 20260828L
LIN6 <- c("Hepatocytes", "Macrophages", "Mono+mono derived cells", "cDC2s", "T cells", "pDCs")
# Floor (a): the measured drop-stage-0 fixed_projection range from the scored cohorts.
FLOOR_A <- c(lo = 0.465325, hi = 0.616087)

say <- function(...) { cat(sprintf(...), "\n", sep = ""); flush.console() }
abort <- function(...) stop(paste0("CFT-ABORT: ", sprintf(...)), call. = FALSE)
chk <- function(cond, ...) if (!isTRUE(cond)) abort(...)
base_id <- function(x) sub("\\.[0-9]+$", "", as.character(x))
sha_file <- function(p) digest(file = p, algo = "sha256")
sp_rho <- function(x, y) suppressWarnings(cor(x, y, method = "spearman"))

guard <- new.env(); guard$rows <- list()
record <- function(check, expected, observed, ok, fatal = TRUE) {
  guard$rows[[length(guard$rows) + 1L]] <- data.table(
    check = check, expected = as.character(expected),
    observed = as.character(observed), pass = ok, fatal = fatal)
  if (!ok && fatal) abort("guard failed: %s (expected %s, observed %s)",
                          check, as.character(expected), as.character(observed))
  invisible(TRUE)
}
close_to <- function(a, b, tol) is.finite(a) && is.finite(b) && abs(a - b) <= tol

## ============================ TIER 1: frozen axis provenance ================
say("[guard] frozen axis provenance")
load_tab <- fread(file.path(RELEASE, "projection", "fixed_projection_loadings.tsv"))
load_tab[, gene_id_base := base_id(gsub('"', "", gene_id_base))]
record("n_rows_in_loadings_file", 145L, nrow(load_tab), nrow(load_tab) == 145L)
sig139 <- load_tab$gene_id_base[load_tab$used_for_fixed_projection %in% c(TRUE, "TRUE")]
record("n_used_for_fixed_projection", 139L, length(sig139), length(sig139) == 139L)
sig_sha <- digest(paste0(paste(sort(sig139), collapse = "\n"), "\n"), algo = "sha256", serialize = FALSE)
record("signature_sorted_139_sha256", "165b9168648651b16d6ffe59b4d1f1c18515308e2fd63701cc21d1c6a43ba8ff",
       sig_sha, identical(sig_sha, "165b9168648651b16d6ffe59b4d1f1c18515308e2fd63701cc21d1c6a43ba8ff"))

model <- readRDS(file.path(RELEASE, "reproduction", "discovery_signature_model.rds"))
loading <- model$discovery_loading_oriented; names(loading) <- base_id(names(loading))
centre  <- model$discovery_center;            names(centre)  <- base_id(names(centre))
record("frozen_loading_sum", 9.774235120123, sprintf("%.12f", sum(loading)),
       close_to(sum(loading), 9.774235120123, 1e-9))
record("frozen_loading_sumsq", 1.0, sprintf("%.12f", sum(loading^2)), close_to(sum(loading^2), 1, 1e-9))
record("frozen_loading_max_abs", 0.220510763499, sprintf("%.12f", max(abs(loading))),
       close_to(max(abs(loading)), 0.220510763499, 1e-9))
record("discovery_centre_sum", 771.590375354, sprintf("%.9f", sum(centre)),
       close_to(sum(centre), 771.590375354, 1e-8))
chk(setequal(names(loading), sig139), "frozen loading names are not the 139-gene set")

M <- readRDS(file.path(RELEASE, "reproduction", "discovery_normalized_expression.rds"))
rownames(M) <- base_id(rownames(M))
record("discovery_n_genes", 17074L, nrow(M), nrow(M) == 17074L)
record("discovery_n_samples", 135L, ncol(M), ncol(M) == 135L)
rep_scores <- fread(file.path(RELEASE, "reproduction", "participant_scores.tsv"))
rp <- rep_scores[match(colnames(M), rep_scores$sample_id)]
released_pc1_operational <- as.numeric(rp$pc1_released_operational)
chk(all(is.finite(released_pc1_operational)), "released operational PC1 has non-finite values")

project_set <- function(mat, g, ld, ct) {
  as.numeric(crossprod(ld[g], sweep(mat[g, , drop = FALSE], 1L, ct[g], FUN = "-")))
}

## ============================ TIER 2: counting-convention audit =============
say("[guard] counting-convention audit: raw featureCounts vs merged_dge")
dge <- readRDS(DGE)
counts_dge <- dge$counts
chk(ncol(counts_dge) == 844L, "merged_dge does not carry 844 samples")
record("merged_dge_n_genes", 23370L, nrow(counts_dge), nrow(counts_dge) == 23370L)
dge_universe_versioned <- rownames(counts_dge)
dge_universe_base <- base_id(dge_universe_versioned)
chk(!anyDuplicated(dge_universe_base), "merged_dge base gene ids are not unique")

read_featurecounts <- function(cohort) {
  f <- file.path(COHORT_RES, cohort, "counts", "featurecounts", "gene_counts.txt")
  chk(file.exists(f), "featureCounts file missing for %s", cohort)
  cmdline <- readLines(f, n = 1L)
  tb <- fread(f, skip = 1L, sep = "\t", header = TRUE)
  gid <- tb[[1]]
  mat <- as.matrix(tb[, 7:ncol(tb), with = FALSE])
  colnames(mat) <- sub("\\.Aligned.*$", "", basename(colnames(mat)))
  rownames(mat) <- gid
  storage.mode(mat) <- "integer"
  list(mat = mat, cmdline = cmdline,
       paired_flag = grepl('"-p"', cmdline, fixed = TRUE),
       both_ends_flag = grepl('"-B"', cmdline, fixed = TRUE))
}

audit_rows <- list()
for (cid in c("GSE162694", "GSE213621", "GSE130970", "GSE135251")) {
  fc <- read_featurecounts(cid)
  smp <- intersect(colnames(fc$mat), rownames(dge$samples)[dge$samples$dataset == cid])
  chk(length(smp) > 0, "%s: zero-sample join between featureCounts and merged_dge -- RAISING", cid)
  g <- intersect(rownames(fc$mat), dge_universe_versioned)
  chk(length(g) > 1000L, "%s: fewer than 1000 shared genes with merged_dge", cid)
  A <- fc$mat[g, smp, drop = FALSE]; B <- counts_dge[g, smp, drop = FALSE]
  ident <- sum(A == B); tot <- length(A)
  ratio <- sum(as.numeric(B)) / sum(as.numeric(A))
  audit_rows[[length(audit_rows) + 1L]] <- data.table(
    cohort = cid, featurecounts_p_flag = fc$paired_flag, featurecounts_B_flag = fc$both_ends_flag,
    n_samples_compared = length(smp), n_genes_compared = length(g),
    fraction_cells_identical = ident / tot,
    total_count_ratio_dge_over_raw = ratio)
  record(paste0("counting_convention_identical:", cid), "1.0", sprintf("%.6f", ident / tot),
         ident == tot, fatal = FALSE)
  rm(fc, A, B); gc()
}
for (cid in c("GSE174478", "GSE240729", "GSE167523", "PRJNA512027")) {
  fc <- read_featurecounts(cid)
  audit_rows[[length(audit_rows) + 1L]] <- data.table(
    cohort = cid, featurecounts_p_flag = fc$paired_flag, featurecounts_B_flag = fc$both_ends_flag,
    n_samples_compared = NA_integer_, n_genes_compared = NA_integer_,
    fraction_cells_identical = NA_real_, total_count_ratio_dge_over_raw = NA_real_)
  rm(fc); gc()
}
fwrite(rbindlist(audit_rows), file.path(OUT, "counting_convention_audit.tsv"), sep = "\t")

## ============================ frozen normalisation recipe ==================
normalise_cohort <- function(counts, cohort_id) {
  chk(is.matrix(counts), "%s counts are not a matrix", cohort_id)
  chk(!anyNA(counts) && all(is.finite(counts)), "%s counts contain non-finite values", cohort_id)
  chk(all(counts >= 0) && all(counts == floor(counts)),
      "%s input is not a raw non-negative integer count matrix", cohort_id)
  keep <- rowSums(counts) > ncol(counts)
  chk(sum(keep) >= 1000L, "%s retains fewer than 1,000 genes after the row-sum filter", cohort_id)
  log_expression <- log2(counts[keep, , drop = FALSE] + 1)
  normalised <- normalizeQuantiles(log_expression)   # preprocessCore aborts under cgroups
  dimnames(normalised) <- dimnames(log_expression)
  chk(all(is.finite(normalised)), "%s quantile normalisation returned non-finite values", cohort_id)
  normalised
}

## reproduce a released cohort matrix straight from raw featureCounts, restricted
## to the merged_dge substrate: proves the whole raw -> score chain, not just parts.
say("[guard] end-to-end reproduction of a released cohort score from raw featureCounts")
repro_rows <- list()
proj_scores <- fread(file.path(RELEASE, "projection", "participant_scores.tsv"))
for (cc in c("sample_id", "dataset", "axis_id")) proj_scores[[cc]] <- gsub('"', "", proj_scores[[cc]])
for (cid in c("GSE162694", "GSE213621")) {
  fc <- read_featurecounts(cid)
  smp <- rownames(dge$samples)[dge$samples$dataset == cid]
  chk(all(smp %in% colnames(fc$mat)), "%s: released samples absent from raw featureCounts", cid)
  g <- intersect(dge_universe_versioned, rownames(fc$mat))
  cm <- fc$mat[g, smp, drop = FALSE]
  rownames(cm) <- base_id(rownames(cm))
  Xc <- normalise_cohort(cm, cid)
  gg <- sig139[sig139 %in% rownames(Xc)]
  sc <- setNames(project_set(Xc, gg, loading, centre), colnames(Xc))
  ref <- proj_scores[dataset == cid & axis_id == "fixed_projection"]
  ref <- setNames(as.numeric(ref$axis_raw), ref$sample_id)
  common <- intersect(names(sc), names(ref))
  chk(length(common) > 0, "%s: zero-row join to released scores -- RAISING", cid)
  repro_rows[[length(repro_rows) + 1L]] <- data.table(
    cohort = cid, n_samples = length(common), n_signature_present = length(gg),
    max_abs_difference = max(abs(sc[common] - ref[common])),
    spearman_vs_released = sp_rho(sc[common], ref[common]),
    pearson_vs_released = cor(sc[common], ref[common]))
  rm(fc, cm, Xc); gc()
}
rep_dt <- rbindlist(repro_rows)
fwrite(rep_dt, file.path(OUT, "released_cohort_reproduction_from_raw_counts.tsv"), sep = "\t")
record("raw_to_score_chain_min_spearman", ">=0.999", sprintf("%.6f", min(rep_dt$spearman_vs_released)),
       min(rep_dt$spearman_vs_released) >= 0.999, fatal = FALSE)

## ============================ sealed null gene sets =========================
say("[null] loading sealed 1000-draw gene sets")
read_sets <- function(gzf, digf, coll_sha) {
  tb <- fread(gzf)
  chk(nrow(tb) == B_DRAWS * 139L, "%s: expected %d rows, got %d", gzf, B_DRAWS * 139L, nrow(tb))
  sets <- split(tb$gene_id_base, tb$draw_index)
  sets <- sets[as.character(seq_len(B_DRAWS))]
  per <- vapply(sets, function(s) digest(paste(sort(s), collapse = "\n"), algo = "sha256", serialize = FALSE), character(1))
  dg <- fread(digf)
  chk(identical(unname(per), dg$set_sha256), "%s: per-draw digests do not match row for row", digf)
  coll <- digest(paste(unname(per), collapse = "\n"), algo = "sha256", serialize = FALSE)
  chk(identical(coll, coll_sha), "%s: collection digest mismatch", gzf)
  lapply(sets, as.character)
}
setsA <- read_sets(file.path(PRESPEC, "null_gene_sets.tsv.gz"),
                   file.path(PRESPEC, "null_gene_set_digests.tsv"),
                   "1b6e2f040329e25faf8e9be85e0a8e5ed2920fdb2c232b14553143733d3f7e36")
setsB <- read_sets(file.path(PRESPEC, "null_gene_sets_caliper.tsv.gz"),
                   file.path(PRESPEC, "null_gene_set_digests_caliper.tsv"),
                   "275454cce7a32de013b5efdb23a21c8ae3e1bba1d734908cc1518d5fd89c1fdb")
record("sealed_null_sets_loaded", "A+B verified", "A+B verified", TRUE)

precompute <- function(sets, tag) {
  LD <- matrix(NA_real_, nrow = 139L, ncol = B_DRAWS)
  CT <- matrix(NA_real_, nrow = 139L, ncol = B_DRAWS)
  GN <- matrix(NA_character_, nrow = 139L, ncol = B_DRAWS)
  bad <- integer(0)
  for (b in seq_len(B_DRAWS)) {
    gs <- sets[[b]]
    pca <- prcomp(t(M[gs, , drop = FALSE]), center = TRUE, scale. = FALSE)
    sgn <- sign(cor(pca$x[, 1L], released_pc1_operational, method = "pearson"))
    if (sgn == 0) { bad <- c(bad, b); next }
    ld <- pca$rotation[, 1L] * sgn
    GN[, b] <- names(ld); LD[, b] <- unname(ld); CT[, b] <- unname(pca$center[names(ld)])
  }
  if (length(bad)) say("[null %s] %d draws discarded for zero orientation sign", tag, length(bad))
  list(LD = LD, CT = CT, GN = GN, bad = bad)
}
NULL_UNIVERSE <- fread(file.path(PRESPEC, "sampling_universe.tsv"))$gene_id_base
say("[null] precomputing discovery PC1 loadings for the sealed draws")
PRE <- list(A = precompute(setsA, "A"), B = precompute(setsB, "B"))

## ============================ helpers ======================================
clr6 <- function(P) {
  P <- as.matrix(P); P[!is.finite(P)] <- 0
  rs <- rowSums(P)
  chk(all(rs > 0), "a sample has zero total across the six lineages")
  P <- P / rs
  nz <- P[P > 0]; chk(length(nz) > 0, "no non-zero proportions")
  P <- P + 0.5 * min(nz); P <- P / rowSums(P)
  L <- log(P); L - rowMeans(L)
}
partial_spearman <- function(x, y, Z) {
  rx <- rank(x); ry <- rank(y)
  Zr <- apply(as.matrix(Z), 2L, rank)
  ex <- residuals(lm(rx ~ Zr)); ey <- residuals(lm(ry ~ Zr))
  suppressWarnings(cor(ex, ey))
}
tie_ceiling <- function(outcome, tiebreak) {
  s <- split(seq_along(outcome), outcome)
  star <- numeric(length(outcome))
  for (g in s) star[g] <- outcome[g[1]] + rank(tiebreak[g], ties.method = "first") / (length(g) + 1)
  sp_rho(star, outcome)
}
mde80 <- function(q95, n) tanh(atanh(q95) + qnorm(0.80) * 1.06 / sqrt(n - 3))

## ============================ cohort metadata ==============================
say("[meta] harmonising recorded histology for the control-free cohorts")
meta_174478 <- {
  d <- fread(file.path(META_ROOT, "GSE174478/metadata/SraRunTable.csv"))
  chk(nrow(d) == 94L, "GSE174478 metadata has %d rows not 94", nrow(d))
  chk(all(c("Run", "fibrosis_stage", "nas_score", "age", "sex", "disease_state") %in% names(d)),
      "GSE174478 metadata is missing a required column")
  chk(all(d$disease_state == "NAFLD"), "GSE174478 carries a non-NAFLD disease_state value")
  data.table(sample_id = d$Run, dataset = "GSE174478",
             fibrosis_stage = as.numeric(d$fibrosis_stage), nas_score = as.numeric(d$nas_score),
             age = as.numeric(d$age), sex = as.character(d$sex),
             diagnosis = as.character(d$disease_state), is_control = FALSE)
}
meta_240729 <- {
  d <- fread(file.path(META_ROOT, "GSE240729/metadata/SraRunTable.csv"))
  chk(nrow(d) == 67L, "GSE240729 metadata has %d rows not 67", nrow(d))
  chk("fibrosisscore" %in% names(d), "GSE240729 metadata lacks fibrosisscore")
  fs <- as.numeric(sub("^F", "", d$fibrosisscore))
  chk(!anyNA(fs), "GSE240729 fibrosisscore did not parse to an ordinal stage")
  chk(all(grepl("FFPE", d$tissue, fixed = TRUE)), "GSE240729 tissue is not FFPE as recorded")
  data.table(sample_id = d$Run, dataset = "GSE240729",
             fibrosis_stage = fs, nas_score = NA_real_,
             age = NA_real_, sex = NA_character_,
             diagnosis = "MASLD_FFPE", is_control = FALSE)
}
meta_167523 <- {
  d <- fread(file.path(META_ROOT, "GSE167523/metadata/SraRunTable.csv"))
  chk(nrow(d) == 98L, "GSE167523 metadata has %d rows not 98", nrow(d))
  data.table(sample_id = d$Run, dataset = "GSE167523",
             fibrosis_stage = NA_real_, nas_score = NA_real_,
             age = suppressWarnings(as.numeric(d$AGE)), sex = as.character(d$gender),
             diagnosis = as.character(d$disease_subtype),
             is_control = FALSE,
             binary_nash = as.integer(d$disease_subtype == "NASH"))
}
meta_512027 <- {
  d <- fread(file.path(META_ROOT, "PRJNA512027/metadata/SraRunTable.csv"))
  chk("disease_stage" %in% names(d), "PRJNA512027 metadata lacks disease_stage")
  data.table(sample_id = d$Run, dataset = "PRJNA512027",
             fibrosis_stage = NA_real_, nas_score = NA_real_,
             age = suppressWarnings(as.numeric(d$AGE)), sex = as.character(d$sex),
             diagnosis = as.character(d$disease_stage),
             is_control = d$disease_stage == "NORMAL",
             binary_disease = as.integer(d$disease_stage != "NORMAL"))
}
meta_all <- list(GSE174478 = meta_174478, GSE240729 = meta_240729,
                 GSE167523 = meta_167523, PRJNA512027 = meta_512027)

read_props <- function(cohort) {
  f <- file.path(PROP_ROOT, cohort, paste0(cohort, "_bayesprism_proportions.tsv"))
  if (!file.exists(f)) return(NULL)
  # ragged: 16 header names, 17 data fields (unnamed row-name index at column 0)
  hdr <- strsplit(readLines(f, n = 1L), "\t", fixed = TRUE)[[1]]
  tb <- fread(f, skip = 1L, header = FALSE, sep = "\t")
  chk(ncol(tb) == length(hdr) + 1L,
      "%s proportions are not ragged by exactly one column (hdr %d, data %d)",
      cohort, length(hdr), ncol(tb))
  setnames(tb, c("sample_id", hdr))
  chk(all(LIN6 %in% names(tb)), "%s proportions lack a required lineage column", cohort)
  tb[, c("sample_id", LIN6), with = FALSE]
}

## ============================ scoring the new cohorts ======================
score_cohort <- function(cohort, universe = c("merged_dge_substrate", "full_featurecounts")) {
  universe <- match.arg(universe)
  fc <- read_featurecounts(cohort)
  cm <- fc$mat
  if (universe == "merged_dge_substrate") {
    g <- intersect(dge_universe_versioned, rownames(cm))
    chk(length(g) > 10000L, "%s: only %d genes shared with the merged_dge substrate", cohort, length(g))
    cm <- cm[g, , drop = FALSE]
  }
  rownames(cm) <- base_id(rownames(cm))
  cm <- cm[!duplicated(rownames(cm)), , drop = FALSE]
  lib <- colSums(cm)
  Xc <- normalise_cohort(cm, cohort)
  gg <- sig139[sig139 %in% rownames(Xc)]
  sc <- setNames(project_set(Xc, gg, loading, centre), colnames(Xc))
  list(X = Xc, score = sc, n_sig_present = length(gg), lib = lib,
       n_genes_input = nrow(cm), n_genes_kept = nrow(Xc))
}

PRIMARY   <- c("GSE174478", "GSE240729")
SECONDARY <- c("GSE167523", "PRJNA512027")

per_sample <- list(); anchors <- list(); nullres <- list(); XMAT <- list()
covrows <- list(); comp_rows <- list(); partrows <- list(); sens_rows <- list()

for (cohort in c(PRIMARY, SECONDARY)) {
  say("[score] %s", cohort)
  S <- score_cohort(cohort, "merged_dge_substrate")
  md <- meta_all[[cohort]]
  common <- intersect(names(S$score), md$sample_id)
  chk(length(common) > 0, "%s: zero-row join between scores and metadata -- RAISING", cohort)
  chk(length(common) >= 0.9 * nrow(md), "%s: only %d of %d metadata rows joined to counts",
      cohort, length(common), nrow(md))
  d <- md[match(common, md$sample_id)]
  d[, fixed_projection := as.numeric(S$score[common])]
  d[, log10_libsize := log10(as.numeric(S$lib[common]))]

  pr <- read_props(cohort)
  if (!is.null(pr)) {
    pj <- pr[match(d$sample_id, pr$sample_id)]
    chk(sum(!is.na(pj$sample_id)) > 0, "%s: zero-row join to BayesPrism proportions -- RAISING", cohort)
    for (L in LIN6) d[[L]] <- as.numeric(pj[[L]])
  } else {
    for (L in LIN6) d[[L]] <- NA_real_
  }
  d[, has_composition := rowSums(is.finite(as.matrix(.SD))) == length(LIN6), .SDcols = LIN6]

  covrows[[cohort]] <- data.table(
    cohort = cohort, universe = "merged_dge_substrate",
    n_samples_scored = length(common), n_metadata_rows = nrow(md),
    n_genes_in_universe = S$n_genes_input, n_genes_after_rowsum_filter = S$n_genes_kept,
    n_signature_present = S$n_sig_present, n_signature_imputed = 0L,
    n_with_composition = sum(d$has_composition),
    n_null_universe_genes_present = sum(NULL_UNIVERSE %in% rownames(S$X)),
    n_null_universe_total = length(NULL_UNIVERSE))
  record(paste0("signature_coverage:", cohort), 139L, S$n_sig_present, S$n_sig_present == 139L,
         fatal = FALSE)

  # sensitivity: full featureCounts universe
  S2 <- score_cohort(cohort, "full_featurecounts")
  d[, fixed_projection_full_universe := as.numeric(S2$score[d$sample_id])]
  sens_rows[[cohort]] <- data.table(
    cohort = cohort, n = nrow(d), n_sig_present_full = S2$n_sig_present,
    n_genes_kept_full = S2$n_genes_kept,
    spearman_substrate_vs_full = sp_rho(d$fixed_projection, d$fixed_projection_full_universe))

  per_sample[[cohort]] <- d
  XMAT[[cohort]] <- S$X
  rm(S, S2); gc()
}

## deposit the scoring layer before any inference runs
fwrite(rbindlist(per_sample, fill = TRUE), file.path(OUT, "per_sample_scores.tsv"), sep = "\t")
fwrite(rbindlist(covrows), file.path(OUT, "coverage_by_cohort.tsv"), sep = "\t")
fwrite(rbindlist(sens_rows), file.path(OUT, "gene_universe_sensitivity.tsv"), sep = "\t")
fwrite(rbindlist(guard$rows), file.path(OUT, "reproduction_guard.tsv"), sep = "\t")

## ---- strata: with zero controls, S2 == S1 by construction; record it -------
strata_for <- function(cohort, d) {
  out <- list()
  if (cohort %in% PRIMARY) {
    dd <- d[is.finite(fibrosis_stage)]
    out[["S1_all_sample"]] <- dd
    out[["S2_within_disease"]] <- dd[is_control == FALSE]
    out[["S3_drop_stage_zero"]] <- dd[fibrosis_stage != 0]
  } else if (cohort == "GSE167523") {
    out[["SEC_binary_NAFL_vs_NASH"]] <- d
  } else {
    out[["SEC_S1_all_sample"]] <- d
    out[["SEC_S2_within_disease_controls_dropped"]] <- d[is_control == FALSE]
  }
  out
}

outcome_for <- function(cohort, dd) {
  if (cohort == "GSE174478") return(list(fibrosis_stage = dd$fibrosis_stage,
                                         nas_score = dd$nas_score))
  if (cohort %in% PRIMARY) return(list(fibrosis_stage = dd$fibrosis_stage))
  if (cohort == "GSE167523") return(list(binary_nash = as.numeric(dd$binary_nash)))
  list(binary_disease = as.numeric(dd$binary_disease))
}

## ============================ anchors + matched null =======================
for (cohort in c(PRIMARY, SECONDARY)) {
  d <- per_sample[[cohort]]
  STR <- strata_for(cohort, d)
  Xc <- XMAT[[cohort]]
  for (st in names(STR)) {
    dd <- STR[[st]]
    if (nrow(dd) < 10L) next
    outs <- outcome_for(cohort, dd)
    for (onm in names(outs)) {
      y <- outs[[onm]]
      keep <- is.finite(y)
      ddk <- dd[keep]; yk <- y[keep]
      if (length(unique(yk)) < 2L) next
      obs <- sp_rho(ddk$fixed_projection, yk)
      obs_partial_lib <- partial_spearman(ddk$fixed_projection, yk, ddk$log10_libsize)
      ceil <- tie_ceiling(yk, ddk$fixed_projection)
      pos <- match(ddk$sample_id, colnames(Xc))
      chk(!anyNA(pos), "%s/%s: sample not found in the cohort matrix -- RAISING", cohort, st)
      for (nl in c("A", "B")) {
        P <- PRE[[nl]]
        rho <- rep(NA_real_, B_DRAWS); rhop <- rep(NA_real_, B_DRAWS)
        n_incomplete <- 0L
        for (b in seq_len(B_DRAWS)) {
          if (b %in% P$bad) next
          g <- P$GN[, b]
          if (!all(g %in% rownames(Xc))) { n_incomplete <- n_incomplete + 1L; next }
          ld <- setNames(P$LD[, b], g); ct <- setNames(P$CT[, b], g)
          s <- project_set(Xc, g, ld, ct)
          sv <- s[pos]
          rho[b]  <- sp_rho(sv, yk)
          rhop[b] <- partial_spearman(sv, yk, ddk$log10_libsize)
        }
        rv <- rho[is.finite(rho)]; rpv <- rhop[is.finite(rhop)]
        if (cohort %in% PRIMARY) {
          chk(length(rv) >= 800L, "%s/%s/%s null%s: only %d usable draws in a PRIMARY cohort -- RAISING",
              cohort, st, onm, nl, length(rv))
        } else if (length(rv) < 200L) {
          # secondary arm: the cohort's own row-sum filter removes null-draw genes.
          # Record the shortfall explicitly; never silently degrade to a skip.
          nullres[[length(nullres) + 1L]] <- data.table(
            cohort = cohort, stratum = st, outcome = onm, null = nl, claim = "SECONDARY",
            n = nrow(ddk), n_draws_retained = length(rv), n_draws_incomplete_coverage = n_incomplete,
            observed_rho = obs, null_median = NA_real_, null_q05 = NA_real_, null_q95 = NA_real_,
            null_q99 = NA_real_, null_mean = NA_real_, null_sd = NA_real_,
            margin_over_median = NA_real_, p_emp_signed = NA_real_,
            observed_abs_rho = abs(obs), null_abs_q95 = NA_real_, p_emp_abs = NA_real_,
            observed_partial_rho = obs_partial_lib, null_partial_q95 = NA_real_,
            p_emp_partial = NA_real_, tie_ceiling = ceil, obs_over_ceiling = obs / ceil,
            mde80_vs_null_q95 = NA_real_,
            inside_within_disease_floor_a = NA, above_floor_a_lower_bound = NA,
            null_not_evaluable_reason = sprintf("only %d of %d draws had full gene coverage in this cohort", length(rv), B_DRAWS))
          next
        }
        q <- quantile(rv, c(0.05, 0.5, 0.95, 0.99), type = 7, names = FALSE)
        nullres[[length(nullres) + 1L]] <- data.table(
          cohort = cohort, stratum = st, outcome = onm, null = nl,
          claim = ifelse(cohort == "GSE174478", "PRIMARY",
                  ifelse(cohort == "GSE240729", "BOUNDARY_FFPE", "SECONDARY")),
          n = nrow(ddk), n_draws_retained = length(rv), n_draws_incomplete_coverage = n_incomplete,
          observed_rho = obs, null_median = q[2], null_q05 = q[1], null_q95 = q[3], null_q99 = q[4],
          null_mean = mean(rv), null_sd = sd(rv),
          margin_over_median = obs - q[2],
          p_emp_signed = (1 + sum(rv >= obs)) / (1 + length(rv)),
          observed_abs_rho = abs(obs), null_abs_q95 = quantile(abs(rv), 0.95, type = 7, names = FALSE),
          p_emp_abs = (1 + sum(abs(rv) >= abs(obs))) / (1 + length(rv)),
          observed_partial_rho = obs_partial_lib,
          null_partial_q95 = quantile(rpv, 0.95, type = 7, names = FALSE),
          p_emp_partial = (1 + sum(rpv >= obs_partial_lib)) / (1 + length(rpv)),
          tie_ceiling = ceil, obs_over_ceiling = obs / ceil,
          mde80_vs_null_q95 = mde80(q[3], nrow(ddk)),
          inside_within_disease_floor_a = (obs >= FLOOR_A[["lo"]] & obs <= FLOOR_A[["hi"]]),
          above_floor_a_lower_bound = obs >= FLOOR_A[["lo"]],
          null_not_evaluable_reason = NA_character_)
        fwrite(data.table(draw_index = seq_len(B_DRAWS), rho_signed = rho,
                          rho_partial_libsize = rhop),
               file.path(OUT, sprintf("nullvec_%s_%s_%s_null%s.tsv", cohort, st, onm, nl)), sep = "\t")
      }
      anchors[[length(anchors) + 1L]] <- data.table(
        cohort = cohort, stratum = st, outcome = onm, n = nrow(ddk),
        observed_rho = obs, tie_ceiling = ceil,
        partial_adj_log_libsize = obs_partial_lib)
    }
  }
  gc()
}

## ============================ competitors (floor c) ========================
say("[C4] competitors: log library size, hepatocyte fraction, 6-lineage CLR")
for (cohort in c(PRIMARY, SECONDARY)) {
  d <- per_sample[[cohort]]
  STR <- strata_for(cohort, d)
  for (st in names(STR)) {
    dd <- STR[[st]]
    outs <- outcome_for(cohort, dd)
    for (onm in names(outs)) {
      y <- outs[[onm]]; keep <- is.finite(y) & dd$has_composition
      ddk <- dd[keep]; yk <- y[keep]
      if (nrow(ddk) < 20L || length(unique(yk)) < 2L) next
      Lc <- clr6(ddk[, LIN6, with = FALSE])
      colnames(Lc) <- paste0("CLR_", make.names(LIN6))
      nn <- nrow(ddk)
      set.seed(CV_SEED)
      fk <- sample(rep_len(1:10, nn))
      dfcv <- data.frame(y = yk, Lc)
      oof <- numeric(nn)
      for (k in 1:10) {
        tr <- fk != k
        fit <- lm(y ~ ., data = dfcv[tr, , drop = FALSE])
        oof[!tr] <- predict(fit, newdata = dfcv[!tr, , drop = FALSE])
      }
      cands <- list(continuum = ddk$fixed_projection,
                    C_libsize = ddk$log10_libsize,
                    C_hepatocyte = ddk$Hepatocytes,
                    C_lineage6_oof = oof)
      set.seed(BOOT_SEED)
      IDX <- matrix(sample.int(nn, nn * BOOT_B, replace = TRUE), nrow = nn, ncol = BOOT_B)
      saveRDS(IDX, file.path(OUT, sprintf("bootstrap_indices_%s_%s_%s.rds", cohort, st, onm)))
      absr <- lapply(cands, function(v) vapply(seq_len(BOOT_B), function(j) {
        i <- IDX[, j]; s <- yk[i]
        if (length(unique(s)) < 2L) return(NA_real_)
        abs(sp_rho(v[i], s))
      }, numeric(1)))
      rows <- list()
      for (nm in c("C_libsize", "C_hepatocyte", "C_lineage6_oof")) {
        del <- absr$continuum - absr[[nm]]; del <- del[is.finite(del)]
        ci <- quantile(del, c(0.025, 0.975), type = 7, names = FALSE)
        p <- max(2 * min(mean(del <= 0), mean(del >= 0)), 1 / (length(del) + 1))
        rows[[length(rows) + 1L]] <- data.table(
          cohort = cohort, stratum = st, outcome = onm, competitor = nm, n = nn,
          continuum_abs_rho = abs(sp_rho(cands$continuum, yk)),
          competitor_rho = sp_rho(cands[[nm]], yk),
          competitor_abs_rho = abs(sp_rho(cands[[nm]], yk)),
          diff_abs_rho = abs(sp_rho(cands$continuum, yk)) - abs(sp_rho(cands[[nm]], yk)),
          boot_diff_median = median(del), boot_ci_lo = ci[1], boot_ci_hi = ci[2],
          win = ci[1] > 0, boot_p = p, n_boot_used = length(del))
      }
      rr <- rbindlist(rows); rr[, boot_p_holm := p.adjust(boot_p, method = "holm")]
      comp_rows[[length(comp_rows) + 1L]] <- rr
      Lm <- as.matrix(Lc)
      partrows[[length(partrows) + 1L]] <- data.table(
        cohort = cohort, stratum = st, outcome = onm, n = nn,
        raw_spearman = sp_rho(ddk$fixed_projection, yk),
        partial_adj_log_libsize = partial_spearman(ddk$fixed_projection, yk, ddk$log10_libsize),
        partial_adj_hepatocyte = partial_spearman(ddk$fixed_projection, yk, ddk$Hepatocytes),
        partial_adj_6clr_joint = partial_spearman(ddk$fixed_projection, yk, Lm),
        partial_adj_age_sex = if (all(is.finite(ddk$age)) && length(unique(ddk$sex)) > 1L)
          partial_spearman(ddk$fixed_projection, yk,
                           cbind(ddk$age, as.numeric(factor(ddk$sex)))) else NA_real_)
      rm(IDX, absr); gc()
    }
  }
}

## ============================ deposits =====================================
fwrite(rbindlist(per_sample, fill = TRUE), file.path(OUT, "per_sample_scores.tsv"), sep = "\t")
fwrite(rbindlist(covrows), file.path(OUT, "coverage_by_cohort.tsv"), sep = "\t")
fwrite(rbindlist(sens_rows), file.path(OUT, "gene_universe_sensitivity.tsv"), sep = "\t")
fwrite(rbindlist(anchors), file.path(OUT, "anchors.tsv"), sep = "\t")
NR <- rbindlist(nullres, fill = TRUE)
fwrite(NR, file.path(OUT, "matched_null_results.tsv"), sep = "\t")
CMP <- rbindlist(comp_rows); fwrite(CMP, file.path(OUT, "competitor_comparisons.tsv"), sep = "\t")
fwrite(rbindlist(partrows), file.path(OUT, "partial_spearman_adjustments.tsv"), sep = "\t")
fwrite(rbindlist(guard$rows), file.path(OUT, "reproduction_guard.tsv"), sep = "\t")

## primary decision family: fixed_projection x fibrosis x Null A, S2 within-disease
prim <- NR[cohort %in% PRIMARY & outcome == "fibrosis_stage" &
             stratum == "S2_within_disease" & null == "A"]
prim[, p_holm_2test := p.adjust(p_emp_signed, method = "holm")]
prim[, F_b_beats_matched_null := observed_rho > null_q95]
prim[, F_a_inside_floor := inside_within_disease_floor_a]
c4 <- CMP[cohort %in% PRIMARY & outcome == "fibrosis_stage" & stratum == "S2_within_disease",
          .(F_c_beats_all_competitors = all(win & boot_p_holm <= 0.05)), by = cohort]
prim <- merge(prim, c4, by = "cohort", all.x = TRUE)
fwrite(prim, file.path(OUT, "PRIMARY_transfer_verdict.tsv"), sep = "\t")

writeLines(c(
  "CFT-v1 control-free transfer of the frozen 139-gene loading",
  "",
  "COUNTING-CONVENTION AUDIT:", capture.output(print(rbindlist(audit_rows))),
  "", "RAW featureCounts -> RELEASED SCORE CHAIN:", capture.output(print(rep_dt)),
  "", "COVERAGE:", capture.output(print(rbindlist(covrows))),
  "", "GENE-UNIVERSE SENSITIVITY:", capture.output(print(rbindlist(sens_rows))),
  "", "ANCHORS:", capture.output(print(rbindlist(anchors))),
  "", "PRIMARY DECISION FAMILY (S2 within-disease, fibrosis, Null A):",
  capture.output(print(prim[, .(cohort, claim, n, observed_rho, null_median, null_q95,
                                margin_over_median, p_emp_signed, p_holm_2test,
                                F_a_inside_floor, F_b_beats_matched_null,
                                F_c_beats_all_competitors, tie_ceiling, obs_over_ceiling)])),
  "", "ALL NULL CELLS:", capture.output(print(NR[, .(cohort, stratum, outcome, null, n,
                                observed_rho, null_median, null_q95, p_emp_signed,
                                margin_over_median, tie_ceiling)])),
  "", "COMPETITORS:", capture.output(print(CMP)),
  "", "PARTIALS:", capture.output(print(rbindlist(partrows)))),
  file.path(OUT, "VERDICT.txt"))

writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
say("DONE -> %s", OUT)
