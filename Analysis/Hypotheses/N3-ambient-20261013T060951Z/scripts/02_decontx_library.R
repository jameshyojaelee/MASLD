#!/usr/bin/env Rscript
# N3 ambient check, step 2: decontX on one library with the empty-droplet profile.
# Usage: 02_decontx_library.R <TMPDIR/LIBRARY> <OUT_RDS> <REPORT_GENES_TXT> <GENES_TXT>
suppressMessages({ library(celda); library(Matrix); library(jsonlite) })
argv <- commandArgs(trailingOnly = TRUE)
stopifnot(length(argv) == 4L)
tmp <- argv[[1]]; out_rds <- argv[[2]]
if (file.exists(out_rds)) stop("Refusing to overwrite ", out_rds)
report_genes <- readLines(argv[[3]]); genes <- readLines(argv[[4]])
SEED <- 20261013L
MIN_EMPTY_UMIS <- 10000
MIN_CLUSTER_CELLS <- 10L
set.seed(SEED)

st <- fromJSON(file.path(tmp, "stats.json"))
cells <- read.delim(file.path(tmp, "cells.tsv"), stringsAsFactors = FALSE)
ng <- length(genes); nc <- nrow(cells)
stopifnot(ng == st$n_genes, nc == st$n_cells)
i <- readBin(file.path(tmp, "indices.bin"), integer(), n = st$nnz, size = 4, endian = "little")
p <- readBin(file.path(tmp, "indptr.bin"), integer(), n = nc + 1L, size = 4, endian = "little")
x <- as.double(readBin(file.path(tmp, "data.bin"), integer(), n = st$nnz, size = 4, endian = "little"))
stopifnot(length(i) == st$nnz, length(x) == st$nnz, length(p) == nc + 1L)
raw <- new("dgCMatrix", i = i, p = p, x = x, Dim = c(ng, nc))
rownames(raw) <- genes; colnames(raw) <- cells$barcode
stopifnot(abs(sum(raw@x) - st$atlas_umis) < 0.5)
bg <- readBin(file.path(tmp, "background.bin"), numeric(), n = ng, size = 8, endian = "little")
stopifnot(length(bg) == ng, abs(sum(bg) - st$empty_umis_atlas_genes) < 0.5)

tab <- table(cells$cell_type)
small <- names(tab)[tab < MIN_CLUSTER_CELLS]
zlab <- cells$cell_type; zlab[zlab %in% small] <- "__pooled_small__"
z <- as.integer(factor(zlab))

use_bg <- sum(bg) >= MIN_EMPTY_UMIS
bgm <- NULL
if (use_bg) {
  # decontX uses only Matrix::rowSums(countsBackground); two columns keep it a matrix.
  bgm <- Matrix(cbind(bg, 0), sparse = TRUE); rownames(bgm) <- genes; colnames(bgm) <- c("empty_sum", "zero")
}
run <- function(background) {
  tryCatch(celda::decontX(x = raw, z = z, background = background, seed = SEED, verbose = FALSE),
           error = function(e) { cat("[decontX] error:", conditionMessage(e), "\n"); NULL })
}
t0 <- Sys.time()
fit <- run(bgm)
profile_source <- if (use_bg) "empty_droplets" else "other_clusters_no_background"
if (is.null(fit) && use_bg) { fit <- run(NULL); profile_source <- "other_clusters_after_background_error" }
minutes <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
if (is.null(fit)) {
  dec <- raw; contam <- rep(NA_real_, nc); status <- "uncorrected_passthrough"; profile_source <- "none"
} else {
  dec <- as(fit$decontXcounts, "CsparseMatrix"); contam <- as.numeric(fit$contamination); status <- "corrected"
}
stopifnot(identical(dim(dec), dim(raw)))

ind <- function(group) {
  lv <- sort(unique(group))
  m <- sparseMatrix(i = seq_along(group), j = match(group, lv), x = 1, dims = c(length(group), length(lv)))
  colnames(m) <- lv; m
}
G <- ind(cells$cell_type)
ri <- match(report_genes, genes); stopifnot(!anyNA(ri))
ct_raw <- as.matrix(raw[ri, , drop = FALSE] %*% G); ct_dec <- as.matrix(dec[ri, , drop = FALSE] %*% G)
rownames(ct_raw) <- rownames(ct_dec) <- report_genes
tot_raw <- Matrix::colSums(raw); tot_dec <- Matrix::colSums(dec)
ct_summary <- data.frame(
  cell_type = colnames(G), n_cells = as.integer(Matrix::colSums(G)),
  raw_umis = as.numeric(tot_raw %*% G), dec_umis = as.numeric(tot_dec %*% G),
  mean_cell_contamination = as.numeric((ifelse(is.na(contam), 0, contam) %*% G) / Matrix::colSums(G)),
  stringsAsFactors = FALSE)
if (status != "corrected") ct_summary$mean_cell_contamination <- NA_real_
hep <- cells$cell_type == "Hepatocytes"
res <- list(
  library = st$library, dataset = st$dataset, status = status, profile_source = profile_source,
  seed = SEED, decontx_minutes = minutes, n_clusters = length(unique(z)), pooled_small_types = small,
  export_stats = st, ct_summary = ct_summary, ct_raw = ct_raw, ct_dec = ct_dec,
  hep_raw = Matrix::rowSums(raw[, hep, drop = FALSE]), hep_dec = Matrix::rowSums(dec[, hep, drop = FALSE]),
  background = bg, celda_version = as.character(packageVersion("celda")))
saveRDS(res, out_rds)
cat(sprintf("[decontX] %s status=%s profile=%s cells=%d hep=%d empty_umis=%.0f hep_contam=%.4f min=%.2f\n",
            st$library, status, profile_source, nc, sum(hep), sum(bg),
            1 - sum(res$hep_dec) / max(1, sum(res$hep_raw)), minutes))
