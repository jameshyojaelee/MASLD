#!/usr/bin/env Rscript
# Freeze the expression/variance-matched sampling universe and the 1,000
# pre-registered random 139-gene draws for the MNC-FIB-v1 prespecification.
# THIS SCRIPT FITS NOTHING. No PCA, no projection, no correlation with any
# outcome. It emits gene identifiers only.
suppressPackageStartupMessages({library(digest)})

RELEASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/histology_anchored_continuum/candidates/hac-continuum-20260818T024923Z"
OUT     <- Sys.getenv("MNC_OUT_DIR"); stopifnot(nzchar(OUT))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
SEED_BASE   <- 20260828L
N_DRAWS     <- 1000L
FIB_COHORTS <- c("GSE130970", "GSE135251", "GSE162694", "GSE213621")
ALL_COHORTS <- c("GSE126848", FIB_COHORTS)
stopifnot(!file.exists(file.path(OUT, "null_gene_sets.tsv.gz")))

base_id <- function(x) sub("\\.[0-9]+$", "", as.character(x))
say <- function(...) cat(sprintf(...), "\n", sep = "")

M <- readRDS(file.path(RELEASE, "reproduction", "discovery_normalized_expression.rds"))
say("discovery matrix: %d genes x %d samples", nrow(M), ncol(M))
stopifnot(nrow(M) == 17074L, ncol(M) == 135L)
rownames(M) <- base_id(rownames(M))
stopifnot(!anyDuplicated(rownames(M)))

res_ids <- base_id(read.delim(file.path(RELEASE, "inputs", "resource_gene_ids.tsv"),
                              stringsAsFactors = FALSE)[[1]])
say("resource substrate: %d gene ids", length(res_ids))

cohort_rows <- lapply(ALL_COHORTS, function(cid) {
  m <- readRDS(file.path(RELEASE, "unsupervised", "normalized_expression_by_cohort",
                         paste0(cid, ".rds")))
  say("%s matrix: %d genes x %d samples", cid, nrow(m), ncol(m))
  base_id(rownames(m))
}); names(cohort_rows) <- ALL_COHORTS

load_tab <- read.delim(file.path(RELEASE, "projection", "fixed_projection_loadings.tsv"),
                       stringsAsFactors = FALSE)
load_tab$gene_id_base <- base_id(gsub('"', "", load_tab$gene_id_base))
sig145 <- load_tab$gene_id_base
sig139 <- load_tab$gene_id_base[load_tab$used_for_fixed_projection %in% c(TRUE, "TRUE")]
stopifnot(length(sig145) == 145L, length(sig139) == 139L)
sig_sha <- digest(paste0(paste(sort(sig139), collapse = "\n"), "\n"),
                  algo = "sha256", serialize = FALSE)
say("signature sha256 (sorted 139): %s", sig_sha)
stopifnot(identical(sig_sha,
  "165b9168648651b16d6ffe59b4d1f1c18515308e2fd63701cc21d1c6a43ba8ff"))

fib_common <- Reduce(intersect, cohort_rows[FIB_COHORTS])
all_common <- Reduce(intersect, cohort_rows[ALL_COHORTS])
pool_all   <- Reduce(intersect, list(rownames(M), res_ids, fib_common))
if (!all(sig139 %in% pool_all)) stop("signature genes absent from the eligible pool")
U <- setdiff(pool_all, sig145)
say("fib-cohort common: %d | five-cohort common: %d", length(fib_common), length(all_common))
say("eligible pool (incl signature): %d | sampling universe U: %d", length(pool_all), length(U))
if (length(U) < 5000L) stop("sampling universe implausibly small")

mu  <- rowMeans(M)
vr  <- apply(M, 1L, stats::var)
muU <- mu[U]; vrU <- vr[U]

mean_breaks <- unname(stats::quantile(muU, probs = seq(0, 1, by = 0.1), type = 7))
mean_breaks[1] <- -Inf; mean_breaks[length(mean_breaks)] <- Inf
assign_decile <- function(v) as.integer(cut(v, breaks = mean_breaks, labels = FALSE,
                                            include.lowest = TRUE, right = TRUE))
decU <- assign_decile(muU)

var_breaks <- lapply(1:10, function(d) {
  b <- unname(stats::quantile(vrU[decU == d], probs = seq(0, 1, by = 0.2), type = 7))
  b[1] <- -Inf; b[length(b)] <- Inf; b
})
assign_vbin <- function(v, d) as.integer(cut(v, breaks = var_breaks[[d]], labels = FALSE,
                                             include.lowest = TRUE, right = TRUE))
vbinU <- vapply(seq_along(U), function(i) assign_vbin(vrU[i], decU[i]), integer(1))
cellU <- sprintf("D%02d_V%d", decU, vbinU)

muS <- mu[sig139]; vrS <- vr[sig139]
decS <- assign_decile(muS)
vbinS <- vapply(seq_along(sig139), function(i) assign_vbin(vrS[i], decS[i]), integer(1))
cellS <- sprintf("D%02d_V%d", decS, vbinS)
need <- table(cellS)
say("occupied strata: %d of 50 | max n_c = %d", length(need), max(need))

pool_by_cell <- split(U, cellU)
fallbacks <- list()
neighbours <- function(cell) {
  d <- as.integer(substr(cell, 2, 3)); v <- as.integer(substr(cell, 6, 6))
  vs <- c(v, v + 1, v - 1, v + 2, v - 2); vs <- vs[vs >= 1 & vs <= 5]
  ds <- c(d, d + 1, d - 1, d + 2, d - 2); ds <- ds[ds >= 1 & ds <= 10]
  cands <- c(sprintf("D%02d_V%d", d, vs), sprintf("D%02d_V%d", ds, v))
  unique(cands[cands != cell])
}
for (cell in names(need)) {
  have <- length(pool_by_cell[[cell]]); nc <- as.integer(need[[cell]])
  if (is.na(have)) have <- 0L
  if (have < nc) fallbacks[[cell]] <- data.frame(cell = cell, n_required = nc,
                                                 n_available = have,
                                                 stringsAsFactors = FALSE)
}
if (length(fallbacks)) {
  say("WARNING: %d strata are under-populated; the widening ladder will be used",
      length(fallbacks))
  write.table(do.call(rbind, fallbacks), file.path(OUT, "stratum_fallbacks.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

draw_one <- function(b) {
  set.seed(SEED_BASE + b)
  picked <- character(0)
  for (cell in sort(names(need))) {
    nc <- as.integer(need[[cell]])
    avail <- setdiff(pool_by_cell[[cell]], picked); if (is.null(avail)) avail <- character(0)
    if (length(avail) >= nc) {
      picked <- c(picked, sample(avail, nc))
    } else {
      take <- avail; picked <- c(picked, take); short <- nc - length(take)
      for (nb in neighbours(cell)) {
        if (short == 0L) break
        av2 <- setdiff(pool_by_cell[[nb]], picked); if (is.null(av2)) av2 <- character(0)
        k <- min(short, length(av2))
        if (k > 0L) { picked <- c(picked, sample(av2, k)); short <- short - k }
      }
      if (short > 0L) stop(sprintf("draw %d: cell %s could not be filled", b, cell))
    }
  }
  stopifnot(length(picked) == 139L, !anyDuplicated(picked),
            !any(picked %in% sig145))
  picked
}

sets <- lapply(seq_len(N_DRAWS), draw_one)
per_sha <- vapply(sets, function(s) digest(paste(sort(s), collapse = "\n"),
                                           algo = "sha256", serialize = FALSE),
                  character(1))
stopifnot(!anyDuplicated(per_sha))

gz <- gzfile(file.path(OUT, "null_gene_sets.tsv.gz"), "w")
write.table(data.frame(draw_index = rep(seq_len(N_DRAWS), each = 139L),
                       gene_id_base = unlist(sets), stringsAsFactors = FALSE),
            gz, sep = "\t", quote = FALSE, row.names = FALSE)
close(gz)
write.table(data.frame(draw_index = seq_len(N_DRAWS), set_sha256 = per_sha,
                       stringsAsFactors = FALSE),
            file.path(OUT, "null_gene_set_digests.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
write.table(data.frame(gene_id_base = U, discovery_mean = unname(muU),
                       discovery_var = unname(vrU), mean_decile = decU,
                       var_quintile_within_decile = vbinU, stratum = cellU,
                       stringsAsFactors = FALSE),
            file.path(OUT, "sampling_universe.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
write.table(data.frame(gene_id_base = sig139, discovery_mean = unname(muS),
                       discovery_var = unname(vrS), mean_decile = decS,
                       var_quintile_within_decile = vbinS, stratum = cellS,
                       stringsAsFactors = FALSE),
            file.path(OUT, "signature_strata.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
write.table(data.frame(stratum = names(need), n_signature_genes = as.integer(need),
                       n_universe_genes = as.integer(vapply(names(need),
                         function(c) length(pool_by_cell[[c]]), integer(1))),
                       stringsAsFactors = FALSE),
            file.path(OUT, "stratum_counts.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
write.table(data.frame(
    key = c("seed_base", "n_draws", "n_discovery_genes", "n_discovery_samples",
            "n_resource_ids", "n_fib_cohort_common", "n_five_cohort_common",
            "n_eligible_pool", "n_sampling_universe", "n_occupied_strata",
            "signature_sorted_sha256", "gene_sets_collection_sha256"),
    value = c(SEED_BASE, N_DRAWS, nrow(M), ncol(M), length(res_ids),
              length(fib_common), length(all_common), length(pool_all), length(U),
              length(need), sig_sha,
              digest(paste(per_sha, collapse = "\n"), algo = "sha256", serialize = FALSE)),
    stringsAsFactors = FALSE),
  file.path(OUT, "universe_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
say("DONE -> %s", OUT)
