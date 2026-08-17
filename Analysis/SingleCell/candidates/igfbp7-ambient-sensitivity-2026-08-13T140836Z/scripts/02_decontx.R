#!/usr/bin/env Rscript
# Step 2/3 of the IGFBP7 ambient-RNA SENSITIVITY ANALYSIS.
#
# SCOPE (binding): sensitivity analysis only. No Hotspot program is
# rediscovered, refit, reweighted, renamed or re-selected. The frozen
# 117-program registry is never read or written here.
#
# Per SEQUENCING DATASET (never pooled across batches), reconstruct the raw
# genes x cells matrix exported by 01_export_counts.py, verify it against the
# Python per-gene checksum, and run celda::decontX with z = annotated cell type
# so that transcripts characteristic of one cell type which appear in another
# cell type's barcodes are modelled as ambient contamination.
#
# Outputs
#   results/02_ambient_by_gene_celltype.tsv  per gene x cell type x dataset:
#       raw/decontaminated sums, detection, and ambient_fraction
#   results/02_contamination_per_cell.tsv.gz per-cell decontX contamination
#   work/<dataset>/dec_*.bin                 decontaminated counts for the
#       frozen hepatocyte scoring universe on the frozen hepatocyte program
#       gene set, for re-scoring in step 3.

suppressMessages({
  library(celda)
  library(Matrix)
  library(jsonlite)
})

SEED <- 20260813L
set.seed(SEED)

ROOT <- Sys.getenv("MASLD_PROJECT_ROOT")
CAND <- Sys.getenv("CAND_ROOT")
stopifnot(nzchar(ROOT), nzchar(CAND))
WORK <- file.path(CAND, "work")
RES <- file.path(CAND, "results")
dir.create(RES, showWarnings = FALSE, recursive = TRUE)

genes_all <- readLines(file.path(WORK, "genes.txt"))
ngenes <- length(genes_all)
report_genes <- readLines(file.path(WORK, "report_genes.txt"))
score_genes <- readLines(file.path(WORK, "score_genes.txt"))
manifest <- fromJSON(file.path(WORK, "manifest.json"), simplifyVector = FALSE)
datasets <- names(manifest$datasets)

# Optional per-dataset sharding so the 7 datasets can run as parallel SLURM
# jobs. Each shard writes its own result files; 02c_merge.py concatenates them.
only <- Sys.getenv("DECONTX_ONLY_DATASET", unset = "")
SHARD <- nzchar(only)
if (SHARD) {
  stopifnot(only %in% datasets)
  datasets <- only
  cat(sprintf("[decontx] SHARD MODE: %s only\n", only))
}
suffix <- if (SHARD) paste0("__", only) else ""

cat(sprintf("[decontx] %d genes; %d report genes; %d score genes\n",
            ngenes, length(report_genes), length(score_genes)))
cat(sprintf("[decontx] datasets: %s\n", paste(datasets, collapse = ", ")))
cat(sprintf("[decontx] skipped:  %s\n", paste(names(manifest$skipped), collapse = ", ")))

report_idx <- match(report_genes, genes_all)
score_idx <- match(score_genes, genes_all)
stopifnot(!anyNA(report_idx), !anyNA(score_idx))

read_matrix <- function(dsdir, info) {
  mats <- list()
  for (ch in info$chunks) {
    ci <- ch$chunk; nnz <- ch$nnz; nc <- ch$ncells
    i <- readBin(file.path(dsdir, sprintf("chunk%d_indices.bin", ci)),
                 what = integer(), n = nnz, size = 4, endian = "little", signed = TRUE)
    p <- readBin(file.path(dsdir, sprintf("chunk%d_indptr.bin", ci)),
                 what = integer(), n = nc + 1L, size = 4, endian = "little", signed = TRUE)
    x <- as.double(readBin(file.path(dsdir, sprintf("chunk%d_data.bin", ci)),
                           what = integer(), n = nnz, size = 4,
                           endian = "little", signed = TRUE))
    stopifnot(length(i) == nnz, length(x) == nnz, length(p) == nc + 1L)
    mats[[length(mats) + 1L]] <- new("dgCMatrix", i = i, p = p, x = x,
                                     Dim = c(as.integer(ngenes), as.integer(nc)))
  }
  m <- if (length(mats) == 1L) mats[[1]] else do.call(cbind, mats)
  rownames(m) <- genes_all
  m
}

verify <- function(m, info) {
  ok <- TRUE
  for (g in names(info$check_rowsum)) {
    want <- as.numeric(info$check_rowsum[[g]])
    got <- sum(m[g, ])
    if (abs(want - got) > 1e-3 * max(1, abs(want))) {
      cat(sprintf("  [VERIFY-FAIL] %s: py=%.1f R=%.1f\n", g, want, got)); ok <- FALSE
    }
  }
  tw <- as.numeric(info$total_count_sum); tg <- sum(m@x)
  if (abs(tw - tg) > 1e-3 * max(1, abs(tw))) {
    cat(sprintf("  [VERIFY-FAIL] total: py=%.1f R=%.1f\n", tw, tg)); ok <- FALSE
  }
  ok
}

# group-wise column sums of a sparse matrix: returns genes x groups
group_sums <- function(m, grp) {
  lev <- sort(unique(grp))
  ind <- sparseMatrix(i = match(grp, lev), j = seq_along(grp), x = 1,
                      dims = c(length(lev), length(grp)))
  out <- as.matrix(m %*% t(ind))
  colnames(out) <- lev
  out
}

gene_rows <- list()
contam_rows <- list()
run_log <- list()

for (dataset in datasets) {
  dsdir <- file.path(WORK, dataset)
  info <- fromJSON(file.path(dsdir, "dims.json"), simplifyVector = FALSE)
  meta <- read.csv(file.path(dsdir, "meta.csv.gz"), stringsAsFactors = FALSE)
  cat(sprintf("\n[decontx] %s: reconstructing %d cells ...\n", dataset, info$ncells))
  m <- read_matrix(dsdir, info)
  stopifnot(ncol(m) == nrow(meta))
  if (!verify(m, info)) stop(sprintf("reconstruction checksum failed for %s", dataset))

  # decontX needs populated clusters; fold cell types with < MIN_CT cells into a
  # single pooled cluster so the EM does not degenerate. Recorded, not silent.
  MIN_CT <- 10L
  tab <- table(meta$cell_type)
  small <- names(tab)[tab < MIN_CT]
  ct_used <- meta$cell_type
  if (length(small)) ct_used[ct_used %in% small] <- "__pooled_small__"
  z <- as.integer(factor(ct_used))

  # The ambient pool is a property of the SEQUENCING RUN, so decontX estimates a
  # separate ambient profile per sample via `batch`. Samples too small to
  # support their own profile are pooled into one per-dataset batch rather than
  # dropped, so the frozen cell census is preserved exactly.
  MIN_BATCH <- 100L
  stab <- table(meta$sample)
  small_samples <- names(stab)[stab < MIN_BATCH]
  batch_used <- meta$sample
  if (length(small_samples)) {
    batch_used[batch_used %in% small_samples] <- "__pooled_small_samples__"
  }
  n_batches <- length(unique(batch_used))
  cat(sprintf("[decontx] %s: %d z clusters (pooled small types: %s); %d ambient batches (pooled small samples: %d)\n",
              dataset, length(unique(z)),
              if (length(small)) paste(small, collapse = ",") else "none",
              n_batches, length(small_samples)))

  set.seed(SEED)
  t0 <- Sys.time()
  corrected <- TRUE
  fail_reason <- ""
  res <- tryCatch(
    celda::decontX(x = m, z = z, batch = batch_used, seed = SEED, verbose = TRUE),
    error = function(e) {
      cat(sprintf("[decontx] %s: decontX FAILED: %s\n", dataset, conditionMessage(e)))
      NULL
    }
  )
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  if (is.null(res)) {
    # Carry raw counts through unchanged so the frozen cell census and donor
    # census are preserved exactly. Reported, never silent.
    corrected <- FALSE
    fail_reason <- "decontX_error"
    dec <- m
    contam <- rep(NA_real_, ncol(m))
    cat(sprintf("[decontx] %s: UNCORRECTED passthrough (%d cells)\n", dataset, ncol(m)))
  } else {
    dec <- res$decontXcounts
    rownames(dec) <- genes_all
    contam <- as.numeric(res$contamination)
    cat(sprintf("[decontx] %s: decontX done in %.1f min; contamination median=%.4f mean=%.4f\n",
                dataset, elapsed, median(contam), mean(contam)))
  }

  # ---- per gene x cell type diagnostics (raw vs decontaminated) ----------
  mr <- m[report_idx, , drop = FALSE]
  dr <- dec[report_idx, , drop = FALSE]
  grp <- meta$cell_type
  # Binary detection matrices built by rewriting the sparse value slot; this
  # avoids the deprecated logical-to-double sparse coercions in Matrix.
  binarize <- function(x, thr) { y <- x; y@x <- as.numeric(y@x > thr); y }
  s_raw <- group_sums(mr, grp)
  s_dec <- group_sums(dr, grp)
  d_raw <- group_sums(binarize(mr, 0), grp)
  d_dec <- group_sums(binarize(dr, 0), grp)
  d_dec_strict <- group_sums(binarize(dr, 0.5), grp)
  n_by_ct <- as.integer(table(factor(grp, levels = colnames(s_raw))))

  for (j in seq_len(ncol(s_raw))) {
    gene_rows[[length(gene_rows) + 1L]] <- data.frame(
      dataset = dataset,
      corrected = corrected,
      cell_type = colnames(s_raw)[j],
      gene = report_genes,
      n_cells = n_by_ct[j],
      sum_raw = s_raw[, j],
      sum_decont = s_dec[, j],
      n_pos_raw = d_raw[, j],
      n_pos_decont = d_dec[, j],
      n_pos_decont_strict = d_dec_strict[, j],
      stringsAsFactors = FALSE
    )
  }
  rm(mr, dr, s_raw, s_dec, d_raw, d_dec, d_dec_strict); gc()

  # ---- per-cell contamination for the hepatocyte scoring universe -------
  hepsel <- which(meta$in_hep_universe == 1L)
  contam_rows[[length(contam_rows) + 1L]] <- data.frame(
    dataset = dataset,
    cell_id = meta$cell_id[hepsel],
    cell_type = meta$cell_type[hepsel],
    sample = meta$sample[hepsel],
    contamination = contam[hepsel],
    umi_raw = Matrix::colSums(m)[hepsel],
    umi_decont = Matrix::colSums(dec)[hepsel],
    stringsAsFactors = FALSE
  )

  # ---- transport corrected counts for re-scoring ------------------------
  if (length(hepsel)) {
    sub <- dec[score_idx, hepsel, drop = FALSE]
    sub <- as(sub, "dgCMatrix")
    writeBin(as.integer(sub@i), file.path(dsdir, "dec_indices.bin"), size = 4, endian = "little")
    writeBin(as.integer(sub@p), file.path(dsdir, "dec_indptr.bin"), size = 4, endian = "little")
    writeBin(as.double(sub@x), file.path(dsdir, "dec_data.bin"), size = 8, endian = "little")
    raw_sub <- m[score_idx, hepsel, drop = FALSE]
    raw_sub <- as(raw_sub, "dgCMatrix")
    writeBin(as.integer(raw_sub@i), file.path(dsdir, "raw_indices.bin"), size = 4, endian = "little")
    writeBin(as.integer(raw_sub@p), file.path(dsdir, "raw_indptr.bin"), size = 4, endian = "little")
    writeBin(as.double(raw_sub@x), file.path(dsdir, "raw_data.bin"), size = 8, endian = "little")
    write.csv(data.frame(cell_id = meta$cell_id[hepsel]),
              file.path(dsdir, "dec_cells.csv"), row.names = FALSE)
    dims <- list(dataset = dataset, ngenes_score = length(score_idx),
                 ncells = length(hepsel),
                 dec_nnz = length(sub@x), raw_nnz = length(raw_sub@x),
                 dec_sum = sum(sub@x), raw_sum = sum(raw_sub@x))
    write(toJSON(dims, auto_unbox = TRUE, digits = 12),
          file.path(dsdir, "dec_dims.json"))
    cat(sprintf("[decontx] %s: wrote corrected counts %d genes x %d cells (nnz=%d)\n",
                dataset, length(score_idx), length(hepsel), length(sub@x)))
    rm(sub, raw_sub)
  }

  run_log[[length(run_log) + 1L]] <- data.frame(
    dataset = dataset, n_cells = ncol(m), n_clusters = length(unique(z)),
    pooled_small_types = paste(small, collapse = ";"),
    n_ambient_batches = n_batches, n_pooled_small_samples = length(small_samples),
    n_hep_universe = length(hepsel),
    corrected = corrected, fail_reason = fail_reason,
    contamination_median = median(contam, na.rm = TRUE),
    contamination_mean = mean(contam, na.rm = TRUE),
    contamination_q10 = as.numeric(quantile(contam, 0.10, na.rm = TRUE)),
    contamination_q90 = as.numeric(quantile(contam, 0.90, na.rm = TRUE)),
    decontx_minutes = elapsed, seed = SEED, stringsAsFactors = FALSE
  )

  rm(m, dec, res, contam); gc()
}

gene_tab <- do.call(rbind, gene_rows)
gene_tab$ambient_fraction <- ifelse(gene_tab$sum_raw > 0,
                                    1 - gene_tab$sum_decont / gene_tab$sum_raw, NA_real_)
write.table(gene_tab,
            file.path(RES, sprintf("02_ambient_by_gene_celltype%s.tsv", suffix)),
            sep = "\t", quote = FALSE, row.names = FALSE)

contam_tab <- do.call(rbind, contam_rows)
gz <- gzfile(file.path(RES, sprintf("02_contamination_per_cell%s.tsv.gz", suffix)), "w")
write.table(contam_tab, gz, sep = "\t", quote = FALSE, row.names = FALSE)
close(gz)

write.table(do.call(rbind, run_log),
            file.path(RES, sprintf("02_decontx_run_log%s.tsv", suffix)),
            sep = "\t", quote = FALSE, row.names = FALSE)

si <- capture.output(sessionInfo())
writeLines(si, file.path(RES, sprintf("02_sessionInfo%s.txt", suffix)))

cat(sprintf("\n[decontx] wrote %d gene x cell type x dataset rows\n", nrow(gene_tab)))
cat("[decontx] DONE\n")
