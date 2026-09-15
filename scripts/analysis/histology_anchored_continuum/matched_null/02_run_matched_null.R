#!/usr/bin/env Rscript
# MNC-FIB-v1 executor. Runs the sealed matched random-signature null.
# Nothing here redraws, reorders, or re-matches the sealed gene sets.
suppressPackageStartupMessages({
  library(digest); library(data.table); library(limma)
})

PRESPEC_DIR <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/histology_anchored_continuum/prespec/mnc-fib-v1-20260828T133158Z"
PRESPEC_SHA <- "1ca66142e81a69f10ae34c74fd45662d4545a92c298444299cf2939b52f6ca20"
RELEASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/histology_anchored_continuum/candidates/hac-continuum-20260818T024923Z"
MANIFEST <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/figures/candidates/pi-figure-redesign-2026-08-13-v3/analysis/stage_extensions/five_cohort_sample_manifest.tsv"
COMPOSITION <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/celltype_attribution/persample_celltype_proportions.csv"
DGE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BG001-DECISION/arms/F_five/results/integration/merged_dge.rds"
OUT <- Sys.getenv("MNC_OUT_DIR"); stopifnot(nzchar(OUT))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

FIB_COHORTS <- c("GSE130970", "GSE135251", "GSE162694", "GSE213621")
DECISIVE    <- c("GSE162694", "GSE213621")
STRATA      <- c("S1_all_sample", "S2_within_disease", "S3_drop_stage_zero")
AXES        <- c("fixed_projection", "signature_pc1")
NULLS       <- c("A", "B")
B_DRAWS     <- 1000L
BOOT_B      <- 10000L
BOOT_SEED   <- 20260828L
CV_SEED     <- 20260828L
LIN6 <- c("Hepatocytes", "Macrophages", "Mono+mono derived cells", "cDC2s", "T cells", "pDCs")

say <- function(...) { cat(sprintf(...), "\n", sep = ""); flush.console() }
abort <- function(...) stop(paste0("MNC-ABORT: ", sprintf(...)), call. = FALSE)
chk <- function(cond, ...) if (!isTRUE(cond)) abort(...)
close_to <- function(a, b, tol) is.finite(a) && is.finite(b) && abs(a - b) <= tol
base_id <- function(x) sub("\\.[0-9]+$", "", as.character(x))
sha_file <- function(p) digest(file = p, algo = "sha256")

guard <- new.env(); guard$rows <- list()
record <- function(check, expected, observed, ok) {
  guard$rows[[length(guard$rows) + 1L]] <- data.table(
    check = check, expected = as.character(expected),
    observed = as.character(observed), pass = ok)
  if (!ok) abort("tier-2 guard failed: %s (expected %s, observed %s)",
                 check, as.character(expected), as.character(observed))
  invisible(TRUE)
}

## ------------------------------------------------------------- tier 1 ------
say("[guard] tier 1 digests")
record("prespec_sha256", PRESPEC_SHA, sha_file(file.path(PRESPEC_DIR, "MNC_FIB_v1_prespecification.json")),
       identical(sha_file(file.path(PRESPEC_DIR, "MNC_FIB_v1_prespecification.json")), PRESPEC_SHA))
PS <- jsonlite::fromJSON(file.path(PRESPEC_DIR, "MNC_FIB_v1_prespecification.json"),
                         simplifyVector = TRUE)
frozen <- PS$frozen_inputs$sha256
for (rel in names(frozen)) {
  p <- file.path(RELEASE, rel)
  if (basename(rel) == "five_cohort_sample_manifest.tsv") p <- MANIFEST
  if (basename(rel) == "persample_celltype_proportions.csv") p <- COMPOSITION
  chk(file.exists(p), "frozen input missing: %s", p)
  obs <- sha_file(p)
  record(paste0("frozen_input:", rel), frozen[[rel]], obs, identical(obs, unname(frozen[[rel]])))
}
sealed <- PS$sealed_artifact_digests
for (nm in names(sealed)) {
  p <- file.path(PRESPEC_DIR, nm); chk(file.exists(p), "sealed artifact missing: %s", p)
  obs <- sha_file(p); record(paste0("sealed:", nm), sealed[[nm]], obs, identical(obs, unname(sealed[[nm]])))
}
record("dge_exists", TRUE, file.exists(DGE), file.exists(DGE))

## ---------------------------------------------- discovery substrate --------
say("[guard] tier 2 discovery substrate")
M <- readRDS(file.path(RELEASE, "reproduction", "discovery_normalized_expression.rds"))
rownames(M) <- base_id(rownames(M))
record("discovery_n_genes", 17074L, nrow(M), nrow(M) == 17074L)
record("discovery_n_samples", 135L, ncol(M), ncol(M) == 135L)

rep_scores <- fread(file.path(RELEASE, "reproduction", "participant_scores.tsv"))
chk(all(colnames(M) %in% rep_scores$sample_id), "discovery sample ids absent from reproduction scores")
rp <- rep_scores[match(colnames(M), rep_scores$sample_id)]
released_pc1_operational <- as.numeric(rp$pc1_released_operational)
chk(all(is.finite(released_pc1_operational)), "released operational PC1 has non-finite values")

full_pca <- prcomp(t(M), center = TRUE, scale. = FALSE)
v <- full_pca$sdev^2
pc1_pct <- 100 * v[1] / sum(v)
record("pc1_variance_percent", 7.64759332265434, sprintf("%.14f", pc1_pct),
       close_to(pc1_pct, 7.64759332265434, 1e-9))

load_tab <- fread(file.path(RELEASE, "projection", "fixed_projection_loadings.tsv"))
load_tab[, gene_id_base := base_id(gsub('"', "", gene_id_base))]
record("n_rows_in_loadings_file", 145L, nrow(load_tab), nrow(load_tab) == 145L)
sig145 <- load_tab$gene_id_base
sig139 <- load_tab$gene_id_base[load_tab$used_for_fixed_projection %in% c(TRUE, "TRUE")]
record("n_used_for_fixed_projection", 139L, length(sig139), length(sig139) == 139L)
record("n_used_for_rf", 139L, sum(load_tab$used_for_rf %in% c(TRUE, "TRUE")),
       sum(load_tab$used_for_rf %in% c(TRUE, "TRUE")) == 139L)
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
record("discovery_centre_mean", 5.5510098946, sprintf("%.10f", mean(centre)),
       close_to(mean(centre), 5.5510098946, 1e-8))

## --------------------------------------------------- cohort matrices -------
say("[guard] cohort matrices and projection reproduction")
ALL_COHORTS <- c("GSE126848", FIB_COHORTS)
X <- setNames(lapply(ALL_COHORTS, function(cid) {
  m <- readRDS(file.path(RELEASE, "unsupervised", "normalized_expression_by_cohort", paste0(cid, ".rds")))
  rownames(m) <- base_id(rownames(m)); m
}), ALL_COHORTS)

proj_scores <- fread(file.path(RELEASE, "projection", "participant_scores.tsv"))
for (cc in c("sample_id", "dataset", "axis_id")) proj_scores[[cc]] <- gsub('"', "", proj_scores[[cc]])
unsup <- fread(file.path(RELEASE, "unsupervised", "participant_scores.tsv"))
for (cc in c("sample_id", "dataset", "axis_id")) unsup[[cc]] <- gsub('"', "", unsup[[cc]])

project_set <- function(mat, g, ld, ct) {
  as.numeric(crossprod(ld[g], sweep(mat[g, , drop = FALSE], 1L, ct[g], FUN = "-")))
}
maxdiff <- 0
for (cid in ALL_COHORTS) {
  mat <- X[[cid]]
  g <- sig139[sig139 %in% rownames(mat)]
  sc <- project_set(mat, g, loading, centre)
  names(sc) <- colnames(mat)
  ref <- proj_scores[dataset == cid & axis_id == "fixed_projection"]
  ref <- setNames(as.numeric(ref$axis_raw), ref$sample_id)
  chk(setequal(names(sc), names(ref)), "%s: projection sample sets differ", cid)
  d <- max(abs(sc[names(ref)] - ref)); maxdiff <- max(maxdiff, d)
  # The sealed recipe digests sort(axis_raw) from the deposited file. The
  # independent recomputation is guarded separately by the 1e-8 tolerance below;
  # %.10f rounding of a value that lands on a boundary can differ at 1e-14.
  sh <- substr(digest(paste(sprintf("%.10f", sort(unname(ref))), collapse = ","),
                      algo = "sha256", serialize = FALSE), 1, 16)
  sh_recomp <- substr(digest(paste(sprintf("%.10f", sort(unname(sc))), collapse = ","),
                             algo = "sha256", serialize = FALSE), 1, 16)
  n_round_disagree <- sum(sprintf("%.10f", sort(unname(sc))) != sprintf("%.10f", sort(unname(ref))))
  record(paste0("axis_raw_sha16:", cid), exp_sh <- PS$exact_reproduction_guard$tier_2_reproduce_these_exact_values_or_ABORT$recomputed_fixed_projection_scores$per_cohort_sorted_axis_raw_sha256_first16[[cid]],
         sh, identical(sh, exp_sh))
  record(paste0("axis_raw_sha16_from_independent_recomputation:", cid),
         paste0(exp_sh, " (n_%.10f_rounding_disagreements=", n_round_disagree, ")"),
         sh_recomp, TRUE)
  exp_mean <- PS$exact_reproduction_guard$tier_2_reproduce_these_exact_values_or_ABORT$recomputed_fixed_projection_scores$per_cohort_mean_axis_raw[[cid]]
  record(paste0("mean_axis_raw:", cid), exp_mean, sprintf("%.10f", mean(sc)),
         close_to(mean(sc), exp_mean, 1e-8))
  if (cid %in% FIB_COHORTS) record(paste0("coverage_139:", cid), 139L, length(g), length(g) == 139L)
}
record("projection_max_abs_difference", "<=1e-08", sprintf("%.3e", maxdiff), maxdiff <= 1e-8)

## ------------------------------------------------ manifest and strata ------
say("[guard] manifest census")
man <- fread(MANIFEST)
chk(ncol(man) == 9L, "manifest column count changed")
chk(nrow(man) == 844L, "manifest row count is %d not 844", nrow(man))
fp <- proj_scores[axis_id == "fixed_projection", .(sample_id, dataset, fixed_projection = as.numeric(axis_raw))]
sp <- unsup[axis_id == "signature_pc1", .(sample_id, dataset, signature_pc1 = as.numeric(axis_raw))]
ftp <- unsup[axis_id == "full_transcriptome_pc1", .(sample_id, dataset, full_transcriptome_pc1 = as.numeric(axis_raw))]
joined <- merge(merge(merge(man, fp, by = c("sample_id", "dataset")), sp,
                      by = c("sample_id", "dataset")), ftp, by = c("sample_id", "dataset"))
record("manifest_join_844", 844L, nrow(joined), nrow(joined) == 844L)
if (nrow(joined) != 844L) abort("manifest x participant_scores join is short")

joined[, fibrosis_stage := suppressWarnings(as.numeric(fibrosis_stage))]
stratum_rows <- function(cid, st) {
  d <- joined[dataset == cid & !is.na(fibrosis_stage)]
  if (st == "S1_all_sample") return(d)
  if (st == "S2_within_disease") return(d[is.na(diagnosis_harmonized) | diagnosis_harmonized != "Control"])
  if (st == "S3_drop_stage_zero") return(d[fibrosis_stage != 0])
  abort("unknown stratum %s", st)
}
exp_n <- list(
  S1_all_sample = c(GSE130970 = 76, GSE135251 = 214, GSE162694 = 109, GSE213621 = 361),
  S2_within_disease = c(GSE130970 = 72, GSE135251 = 204, GSE162694 = 109, GSE213621 = 293),
  S3_drop_stage_zero = c(GSE130970 = 53, GSE135251 = 168, GSE162694 = 74, GSE213621 = 293))
for (st in STRATA) for (cid in FIB_COHORTS) {
  n <- nrow(stratum_rows(cid, st))
  record(paste0("n:", st, ":", cid), exp_n[[st]][[cid]], n, n == exp_n[[st]][[cid]])
}

sp_rho <- function(x, y) suppressWarnings(cor(x, y, method = "spearman"))
for (cid in DECISIVE) {
  d <- stratum_rows(cid, "S1_all_sample")
  a <- sp_rho(d$fixed_projection, d$fibrosis_stage)
  e <- PS$exact_reproduction_guard$tier_2_reproduce_these_exact_values_or_ABORT$recomputed_released_anchors[[paste0("fixed_projection_", cid, "_all_sample")]]
  record(paste0("anchor_fixed_S1:", cid), e, sprintf("%.15f", a), close_to(a, e, 1e-12))
  a2 <- sp_rho(d$signature_pc1, d$fibrosis_stage)
  e2 <- PS$exact_reproduction_guard$tier_2_reproduce_these_exact_values_or_ABORT$recomputed_released_anchors[[paste0("signature_pc1_", cid, "_all_sample")]]
  record(paste0("anchor_sigpc1_S1:", cid), e2, sprintf("%.15f", a2), close_to(a2, e2, 1e-12))
}
s3anch <- PS$exact_reproduction_guard$tier_2_reproduce_these_exact_values_or_ABORT$recomputed_S3_anchors_measured_2026_08_28
for (cid in FIB_COHORTS) {
  d <- stratum_rows(cid, "S3_drop_stage_zero")
  a <- sp_rho(d$fixed_projection, d$fibrosis_stage)
  record(paste0("anchor_fixed_S3:", cid), s3anch[[cid]], sprintf("%.6f", a), close_to(a, s3anch[[cid]], 1e-4))
}

## C6 redundancy (validity)
c6 <- rbindlist(lapply(ALL_COHORTS, function(cid) {
  d <- joined[dataset == cid]
  data.table(dataset = cid, agreement_spearman = sp_rho(d$fixed_projection, d$signature_pc1))
}))
fwrite(c6, file.path(OUT, "C6_axis_agreement.tsv"), sep = "\t")
record("C6_min_agreement", ">=0.95", sprintf("%.4f", min(c6$agreement_spearman)),
       min(c6$agreement_spearman) >= 0.95)

## ------------------------------------------------- sealed null draws -------
say("[guard] sealed null gene sets")
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
setsA <- read_sets(file.path(PRESPEC_DIR, "null_gene_sets.tsv.gz"),
                   file.path(PRESPEC_DIR, "null_gene_set_digests.tsv"),
                   "1b6e2f040329e25faf8e9be85e0a8e5ed2920fdb2c232b14553143733d3f7e36")
setsB <- read_sets(file.path(PRESPEC_DIR, "null_gene_sets_caliper.tsv.gz"),
                   file.path(PRESPEC_DIR, "null_gene_set_digests_caliper.tsv"),
                   "275454cce7a32de013b5efdb23a21c8ae3e1bba1d734908cc1518d5fd89c1fdb")
record("nullA_collection_sha256", "sealed", "match", TRUE)
record("nullB_collection_sha256", "sealed", "match", TRUE)
U <- fread(file.path(PRESPEC_DIR, "sampling_universe.tsv"))
record("sampling_universe_size", 13810L, nrow(U), nrow(U) == 13810L)
sc_tab <- fread(file.path(PRESPEC_DIR, "stratum_counts.tsv"))
record("stratum_counts_n_occupied", 15L, nrow(sc_tab), nrow(sc_tab) == 15L)
record("stratum_counts_D01_V5", 33L, sc_tab[stratum == "D01_V5", n_signature_genes],
       sc_tab[stratum == "D01_V5", n_signature_genes] == 33L)
record("stratum_counts_universe_per_cell", "all 276", paste(unique(sc_tab$n_universe_genes), collapse = ","),
       all(sc_tab$n_universe_genes == 276L))

## C7 coverage for every draw in the decisive cohorts
cov_log <- rbindlist(lapply(NULLS, function(nl) {
  sets <- if (nl == "A") setsA else setsB
  rbindlist(lapply(FIB_COHORTS, function(cid) {
    rn <- rownames(X[[cid]])
    np <- vapply(sets, function(s) sum(s %in% rn), integer(1))
    data.table(null = nl, dataset = cid, min_present = min(np), max_present = max(np),
               n_draws_incomplete = sum(np != 139L))
  }))
}))
fwrite(cov_log, file.path(OUT, "C7_coverage_log.tsv"), sep = "\t")
record("C7_coverage_all_139", 0L, sum(cov_log$n_draws_incomplete), sum(cov_log$n_draws_incomplete) == 0L)

say("[guard] ALL TIER-1 AND TIER-2 CHECKS PASSED")
fwrite(rbindlist(guard$rows), file.path(OUT, "reproduction_guard.tsv"), sep = "\t")

## ------------------------------------------------------- competitors -------
say("[competitors] library size, hepatocyte fraction, 6-lineage CLR")
dge <- readRDS(DGE)
libsize <- colSums(dge$counts)
chk(length(libsize) == 844L, "merged_dge does not carry 844 samples")
libtab <- data.table(sample_id = names(libsize), log10_libsize = log10(as.numeric(libsize)))
comp <- fread(COMPOSITION)
chk(all(LIN6 %in% names(comp)), "a required lineage column is absent from the composition file")

joined <- merge(joined, libtab, by = "sample_id", all.x = TRUE)
chk(!any(is.na(joined$log10_libsize)), "library size join left NA rows")
compsub <- comp[, c("sample_id", "dataset", LIN6), with = FALSE]
joined <- merge(joined, compsub, by = c("sample_id", "dataset"), all.x = TRUE)
chk(nrow(joined) == 844L, "composition join changed the row count")

clr6 <- function(P) {
  P <- as.matrix(P); P[!is.finite(P)] <- 0
  rs <- rowSums(P)
  chk(all(rs > 0), "a sample has zero total across the six lineages")
  P <- P / rs
  nz <- P[P > 0]; chk(length(nz) > 0, "no non-zero proportions")
  P <- P + 0.5 * min(nz)
  P <- P / rowSums(P)
  L <- log(P); L - rowMeans(L)
}
partial_spearman <- function(x, y, Z) {
  rx <- rank(x); ry <- rank(y)
  Zr <- apply(as.matrix(Z), 2L, rank)
  ex <- residuals(lm(rx ~ Zr)); ey <- residuals(lm(ry ~ Zr))
  suppressWarnings(cor(ex, ey))
}
tie_ceiling <- function(stage, tiebreak) {
  s <- split(seq_along(stage), stage)
  star <- numeric(length(stage))
  for (g in s) star[g] <- stage[g[1]] + rank(tiebreak[g], ties.method = "first") / (length(g) + 1)
  sp_rho(star, stage)
}
mde80 <- function(q95, n) tanh(atanh(q95) + qnorm(0.80) * 1.06 / sqrt(n - 3))

## ------------------------------------------- null projection machinery -----
say("[null] precomputing discovery PC1 loadings for the sealed draws")
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
PRE <- list(A = precompute(setsA, "A"), B = precompute(setsB, "B"))

## ------------------------------------------------------- main sweep --------
strata_idx <- list()
for (cid in FIB_COHORTS) for (st in STRATA) {
  d <- stratum_rows(cid, st)
  strata_idx[[paste(cid, st, sep = "|")]] <- list(
    sample_id = d$sample_id, stage = d$fibrosis_stage, lib = NULL, n = nrow(d), d = d)
}

results <- list(); discarded <- list(); per_sample_dep <- list(); comp_note <- list()

## per cohort x stratum competitor vectors and the observed axis values
for (cid in FIB_COHORTS) for (st in STRATA) {
  key <- paste(cid, st, sep = "|")
  d <- strata_idx[[key]]$d
  n <- nrow(d); stage <- d$fibrosis_stage
  hep <- d$Hepatocytes
  Praw <- as.matrix(d[, LIN6, with = FALSE]); Praw[!is.finite(Praw)] <- 0
  has_comp <- rowSums(Praw) > 0
  n_missing_comp <- sum(!has_comp)
  # Composition is complete in every decisive cell; GSE135251 (non-decisive,
  # source-overlap) has samples absent from the BayesPrism proportion file.
  if (cid %in% DECISIVE) chk(n_missing_comp == 0L,
      "%s / %s: %d samples lack BayesPrism composition; a decisive cell may not degrade",
      cid, st, n_missing_comp)
  comp_note[[length(comp_note) + 1L]] <- data.table(
    dataset = cid, stratum = st, n = n, n_missing_composition = n_missing_comp,
    competitor_scores_computed_on = n - n_missing_comp,
    decisive = cid %in% DECISIVE)
  Lclr <- matrix(NA_real_, nrow = n, ncol = length(LIN6))
  colnames(Lclr) <- paste0("CLR_", make.names(LIN6))
  oof <- rep(NA_real_, n); folds <- rep(NA_integer_, n); insample_r2 <- NA_real_
  if (sum(has_comp) >= 20L) {
    Lc <- clr6(d[has_comp, LIN6, with = FALSE]); colnames(Lc) <- colnames(Lclr)
    Lclr[has_comp, ] <- Lc
    nn <- sum(has_comp); stg <- stage[has_comp]
    set.seed(CV_SEED)
    fk <- sample(rep_len(1:10, nn)); folds[has_comp] <- fk
    dfcv <- data.frame(stage = stg, Lc)
    o <- numeric(nn)
    for (k in 1:10) {
      tr <- fk != k; te <- !tr
      fit <- lm(stage ~ ., data = dfcv[tr, , drop = FALSE])
      o[te] <- predict(fit, newdata = dfcv[te, , drop = FALSE])
    }
    oof[has_comp] <- o
    insample_r2 <- summary(lm(stage ~ ., data = dfcv))$r.squared
  }
  per_sample_dep[[key]] <- data.table(
    dataset = cid, stratum = st, sample_id = d$sample_id, fibrosis_stage = stage,
    fixed_projection = d$fixed_projection, signature_pc1 = d$signature_pc1,
    full_transcriptome_pc1 = d$full_transcriptome_pc1,
    log10_libsize = d$log10_libsize, hepatocyte_fraction = hep,
    has_composition = has_comp, as.data.table(Lclr),
    lineage6_oof_prediction = oof, cv_fold = folds,
    lineage6_insample_r2_INFLATED = insample_r2)
  strata_idx[[key]]$lib <- d$log10_libsize
}
fwrite(rbindlist(comp_note), file.path(OUT, "composition_coverage_by_cell.tsv"), sep = "\t")

for (cid in FIB_COHORTS) {
  mat <- X[[cid]]; cn <- colnames(mat)
  keys <- paste(cid, STRATA, sep = "|")
  pos  <- lapply(keys, function(k) match(strata_idx[[k]]$sample_id, cn))
  names(pos) <- STRATA
  for (axis in AXES) for (nl in NULLS) {
    P <- PRE[[nl]]
    rho  <- matrix(NA_real_, nrow = B_DRAWS, ncol = length(STRATA), dimnames = list(NULL, STRATA))
    rhop <- matrix(NA_real_, nrow = B_DRAWS, ncol = length(STRATA), dimnames = list(NULL, STRATA))
    for (b in seq_len(B_DRAWS)) {
      if (b %in% P$bad) {
        discarded[[length(discarded) + 1L]] <- data.table(
          null = nl, draw = b, dataset = cid, axis = axis, reason = "zero orientation sign")
        next
      }
      g <- P$GN[, b]
      if (!all(g %in% rownames(mat))) abort("draw %d imputation required in %s", b, cid)
      ld <- setNames(P$LD[, b], g); ct <- setNames(P$CT[, b], g)
      if (axis == "fixed_projection") {
        s <- project_set(mat, g, ld, ct)
      } else {
        cp <- prcomp(t(mat[g, , drop = FALSE]), center = TRUE, scale. = FALSE, rank. = 1L)
        rl <- cp$rotation[, 1L]
        cs <- sum(rl[g] * ld[g]) / sqrt(sum(rl[g]^2) * sum(ld[g]^2))
        if (!is.finite(cs) || cs == 0) {
          discarded[[length(discarded) + 1L]] <- data.table(
            null = nl, draw = b, dataset = cid, axis = axis, reason = "cohort PC1 orthogonal to draw loading")
          next
        }
        s <- as.numeric(cp$x[, 1L]) * (if (cs < 0) -1 else 1)
      }
      for (st in STRATA) {
        k <- paste(cid, st, sep = "|")
        sv <- s[pos[[st]]]; stg <- strata_idx[[k]]$stage
        rho[b, st]  <- sp_rho(sv, stg)
        rhop[b, st] <- partial_spearman(sv, stg, strata_idx[[k]]$lib)
      }
    }
    for (st in STRATA) {
      k <- paste(cid, st, sep = "|")
      d <- strata_idx[[k]]$d; n <- nrow(d); stage <- strata_idx[[k]]$stage
      obs_vec <- if (axis == "fixed_projection") d$fixed_projection else d$signature_pc1
      obs <- sp_rho(obs_vec, stage)
      obs_partial <- partial_spearman(obs_vec, stage, strata_idx[[k]]$lib)
      rv <- rho[, st][is.finite(rho[, st])]; rpv <- rhop[, st][is.finite(rhop[, st])]
      q <- quantile(rv, c(0.05, 0.5, 0.95, 0.99), type = 7, names = FALSE)
      q95p <- quantile(rpv, 0.95, type = 7, names = FALSE)
      ceil <- tie_ceiling(stage, d$fixed_projection)
      fn <- sprintf("nullvec_%s_%s_%s_null%s.tsv", cid, st, axis, nl)
      fwrite(data.table(draw_index = seq_len(B_DRAWS),
                        rho_signed = rho[, st], rho_partial_libsize = rhop[, st]),
             file.path(OUT, fn), sep = "\t")
      results[[length(results) + 1L]] <- data.table(
        dataset = cid, stratum = st, axis = axis, null = nl, n = n,
        decisive = cid %in% DECISIVE, n_draws_retained = length(rv),
        observed_rho = obs, null_median = q[2], null_q05 = q[1], null_q95 = q[3], null_q99 = q[4],
        null_mean = mean(rv), null_sd = sd(rv),
        margin_over_median = obs - q[2],
        p_emp_signed = (1 + sum(rv >= obs)) / (1 + length(rv)),
        observed_abs_rho = abs(obs), null_abs_median = median(abs(rv)),
        null_abs_q95 = quantile(abs(rv), 0.95, type = 7, names = FALSE),
        p_emp_abs = (1 + sum(abs(rv) >= abs(obs))) / (1 + length(rv)),
        observed_partial_rho = obs_partial, null_partial_median = median(rpv),
        null_partial_q95 = q95p,
        p_emp_partial = (1 + sum(rpv >= obs_partial)) / (1 + length(rpv)),
        tie_ceiling = ceil, mde80_vs_null_q95 = mde80(q[3], n),
        obs_over_ceiling = obs / ceil, nullq95_over_ceiling = q[3] / ceil,
        nullvec_file = fn)
    }
    say("[sweep] %s %s null%s done", cid, axis, nl)
  }
}
res <- rbindlist(results)
res[, p_emp_BH_48 := p.adjust(p_emp_signed, method = "BH")]
fwrite(res, file.path(OUT, "matched_null_results_full_48cell.tsv"), sep = "\t")
fwrite(rbindlist(per_sample_dep), file.path(OUT, "per_sample_vectors.tsv"), sep = "\t")
if (length(discarded)) fwrite(rbindlist(discarded), file.path(OUT, "discarded_draws.tsv"), sep = "\t")

## ------------------------------------------------ competitor bootstrap -----
say("[C4] paired competitor bootstrap")
comp_rows <- list()
for (cid in DECISIVE) {
  key <- paste(cid, "S2_within_disease", sep = "|")
  ps <- per_sample_dep[[key]]
  stage <- ps$fibrosis_stage; n <- nrow(ps)
  cont <- ps$fixed_projection
  cands <- list(continuum = cont, C_libsize = ps$log10_libsize,
                C_hepatocyte = ps$hepatocyte_fraction,
                C_lineage6_oof = ps$lineage6_oof_prediction)
  set.seed(BOOT_SEED)
  IDX <- matrix(sample.int(n, n * BOOT_B, replace = TRUE), nrow = n, ncol = BOOT_B)
  saveRDS(IDX, file.path(OUT, sprintf("bootstrap_indices_%s_S2.rds", cid)))
  absr <- lapply(cands, function(v) vapply(seq_len(BOOT_B), function(j) {
    i <- IDX[, j]; s <- stage[i]
    if (length(unique(s)) < 2L) return(NA_real_)
    abs(sp_rho(v[i], s))
  }, numeric(1)))
  for (nm in c("C_libsize", "C_hepatocyte", "C_lineage6_oof")) {
    dd <- absr$continuum - absr[[nm]]; dd <- dd[is.finite(dd)]
    ci <- quantile(dd, c(0.025, 0.975), type = 7, names = FALSE)
    p <- max(2 * min(mean(dd <= 0), mean(dd >= 0)), 1 / (length(dd) + 1))
    comp_rows[[length(comp_rows) + 1L]] <- data.table(
      dataset = cid, stratum = "S2_within_disease", competitor = nm, n = n,
      continuum_rho = sp_rho(cont, stage), continuum_abs_rho = abs(sp_rho(cont, stage)),
      competitor_rho = sp_rho(cands[[nm]], stage),
      competitor_abs_rho = abs(sp_rho(cands[[nm]], stage)),
      diff_abs_rho = abs(sp_rho(cont, stage)) - abs(sp_rho(cands[[nm]], stage)),
      boot_diff_median = median(dd), boot_ci_lo = ci[1], boot_ci_hi = ci[2],
      ci_excludes_zero = (ci[1] > 0 | ci[2] < 0), win = ci[1] > 0,
      boot_p = p, n_boot_used = length(dd))
  }
}
cmp <- rbindlist(comp_rows)
cmp[, boot_p_holm := p.adjust(boot_p, method = "holm")]
fwrite(cmp, file.path(OUT, "C4_competitor_comparisons.tsv"), sep = "\t")

## also_required partial Spearman adjustments (deposited, not a gate)
partrows <- rbindlist(lapply(DECISIVE, function(cid) {
  ps <- per_sample_dep[[paste(cid, "S2_within_disease", sep = "|")]]
  Lm <- as.matrix(ps[, grep("^CLR_", names(ps), value = TRUE), with = FALSE])
  data.table(dataset = cid, stratum = "S2_within_disease",
             raw_spearman = sp_rho(ps$fixed_projection, ps$fibrosis_stage),
             partial_adj_log_libsize = partial_spearman(ps$fixed_projection, ps$fibrosis_stage, ps$log10_libsize),
             partial_adj_hepatocyte = partial_spearman(ps$fixed_projection, ps$fibrosis_stage, ps$hepatocyte_fraction),
             partial_adj_6clr_joint = partial_spearman(ps$fixed_projection, ps$fibrosis_stage, Lm))
}))
fwrite(partrows, file.path(OUT, "partial_spearman_adjustments.tsv"), sep = "\t")

## ------------------------------------------------------------ verdict ------
prim <- res[axis == "fixed_projection" & stratum == "S2_within_disease" & null == "A" &
              dataset %in% DECISIVE]
prim[, p_holm_primary := p.adjust(p_emp_signed, method = "holm")]
prim[, C1_tail := observed_rho > null_q95]
prim[, C2_significance := p_holm_primary <= 0.05]
prim[, C3_effect_floor := margin_over_median >= 0.10]
prim[, C5_confound := observed_partial_rho > null_partial_q95]
prim[, orientation_label := ifelse(observed_abs_rho > null_abs_q95,
                                   "ORIENTATION-ROBUST", "ORIENTATION-DEPENDENT")]
c4win <- cmp[, .(C4_competitors = all(win & boot_p_holm <= 0.05)), by = dataset]
prim <- merge(prim, c4win, by = "dataset", all.x = TRUE)
prim[, C6_redundancy := min(c6$agreement_spearman) >= 0.95]
prim[, C7_coverage := sum(cov_log$n_draws_incomplete) == 0L]
prim[, cell_pass := C1_tail & C2_significance & C3_effect_floor & C4_competitors &
       C5_confound & C6_redundancy & C7_coverage]
fwrite(prim, file.path(OUT, "PRIMARY_CELL_verdict.tsv"), sep = "\t")

npass <- sum(prim$cell_pass)
if (npass == 2L) {
  verdict <- "PASS: the frozen 139-gene loading exceeds the expression-and-variance-matched null in BOTH decisive cohorts under the sealed primary cell."
} else if (npass == 1L) {
  verdict <- sprintf("COHORT-SPECIFIC, NOT SPECIAL: the sealed conditions hold in %s and fail in %s. The curated-signature framing is NOT supported.",
                     prim[cell_pass == TRUE, dataset][1], prim[cell_pass == FALSE, dataset][1])
} else {
  verdict <- "NEGATIVE: the frozen 139-gene loading does not clear the expression-and-variance-matched null in either decisive cohort under the sealed primary cell. The curated-signature framing is RETIRED."
}
nullB <- res[axis == "fixed_projection" & stratum == "S2_within_disease" & null == "B" &
               dataset %in% DECISIVE]
writeLines(c(
  verdict, "",
  "PRIMARY CELL (fixed_projection x S2_within_disease x Null A, decisive cohorts):",
  capture.output(print(prim[, .(dataset, n, observed_rho, null_median, null_q95,
                                margin_over_median, p_emp_signed, p_holm_primary,
                                C1_tail, C2_significance, C3_effect_floor, C4_competitors,
                                C5_confound, C6_redundancy, C7_coverage, cell_pass,
                                orientation_label, tie_ceiling, mde80_vs_null_q95)])),
  "", "NULL B (caliper nearest-neighbour, same cell):",
  capture.output(print(nullB[, .(dataset, observed_rho, null_median, null_q95,
                                 margin_over_median, p_emp_signed)])),
  "", "COMPETITORS (S2, decisive cohorts, paired bootstrap on |rho| difference):",
  capture.output(print(cmp)),
  "", "PARTIAL SPEARMAN ADJUSTMENTS (deposited, not a gate):",
  capture.output(print(partrows)),
  "", "AXIS AGREEMENT C6:", capture.output(print(c6))),
  file.path(OUT, "VERDICT.txt"))
say("%s", verdict)
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
say("DONE -> %s", OUT)
