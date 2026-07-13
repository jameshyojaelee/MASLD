#!/usr/bin/env Rscript
# GLP-1RA decontX ambient correction (celda::decontX), per dataset, Hep+Endo.
#
# For each dataset: reconstruct the raw genes x cells dgCMatrix from binary CSC
# chunks written by glp1ra_decontx_export.py (VERIFIED against a per-gene checksum),
# run decontX(x, z = cell_type-as-integer, batch = NULL) so hepatocyte ambient
# (ALB, ~85% of the pool) is modeled into the endothelial cells, then aggregate
# raw vs decontaminated detection / magnitude for the incretin axis + controls.
#
# Pre-registered honesty: GLP1R %-expressing may DROP (ambient / off-parenchyma)
# or SURVIVE (real-but-rare). Positive controls fix the interpretation:
#   GCGR, DPP4  -> genuinely expressed, should largely SURVIVE
#   ALB         -> hepatocyte gene, should DROP sharply in ENDOTHELIAL if ambient
#   PECAM1,STAB2-> endothelial markers, should SURVIVE in endothelial
# If GCGR/DPP4 also collapse, decontX is over-stripping -> FLAG (run suspect).

suppressMessages({
  library(celda)
  library(Matrix)
})

ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
WORK <- file.path(ROOT, "RNA-seq/results/glp1ra/scrna_incretin_axis/decontx_work")
OUT  <- file.path(ROOT, "RNA-seq/results/glp1ra/scrna_incretin_axis")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

AXIS <- c("GLP1R", "GIPR", "GCGR", "DPP4", "GCG", "GLP2R")
CTRL <- c("ALB", "PECAM1", "STAB2")
GENES <- c(AXIS, CTRL)
CELLTYPES <- c("Hepatocytes", "Endothelial cells")
SEED <- 42

genes_all <- readLines(file.path(WORK, "genes.txt"))
ngenes <- length(genes_all)
manifest <- jsonlite::fromJSON(file.path(WORK, "manifest.json"), simplifyVector = FALSE)
datasets <- names(manifest$datasets)
cat(sprintf("[run] %d genes; datasets: %s\n", ngenes, paste(datasets, collapse = ", ")))
cat(sprintf("[run] skipped (degenerate): %s\n",
            paste(names(manifest$skipped), collapse = ", ")))

read_matrix <- function(dsdir, info) {
  # reconstruct genes x cells dgCMatrix from CSC chunks, then cbind
   mats <- list()
  for (ch in info$chunks) {
    ci <- ch$chunk; nnz <- ch$nnz; nc <- ch$ncells
    i <- readBin(file.path(dsdir, sprintf("chunk%d_indices.bin", ci)),
                 what = integer(), n = nnz, size = 4, endian = "little", signed = TRUE)
    p <- readBin(file.path(dsdir, sprintf("chunk%d_indptr.bin", ci)),
                 what = integer(), n = nc + 1L, size = 4, endian = "little", signed = TRUE)
    x <- readBin(file.path(dsdir, sprintf("chunk%d_data.bin", ci)),
                 what = double(), n = nnz, size = 8, endian = "little")
    stopifnot(length(i) == nnz, length(x) == nnz, length(p) == nc + 1L)
    m <- new("dgCMatrix", i = i, p = p, x = x, Dim = c(as.integer(ngenes), as.integer(nc)))
    mats[[length(mats) + 1L]] <- m
  }
  m <- if (length(mats) == 1L) mats[[1]] else do.call(cbind, mats)
  rownames(m) <- genes_all
  m
}

verify <- function(m, info) {
  # per-gene rowSum must match Python checksum (reconstruction integrity)
  ok <- TRUE
  for (g in names(info$check_rowsum)) {
    want <- as.numeric(info$check_rowsum[[g]])
    got <- sum(m[g, ])
    if (abs(want - got) > 1e-3 * max(1, abs(want))) {
      cat(sprintf("  [VERIFY-FAIL] %s: py=%.1f R=%.1f\n", g, want, got)); ok <- FALSE
    }
  }
  tot_want <- as.numeric(info$total_count_sum); tot_got <- sum(m@x)
  if (abs(tot_want - tot_got) > 1e-3 * max(1, abs(tot_want))) {
    cat(sprintf("  [VERIFY-FAIL] total: py=%.1f R=%.1f\n", tot_want, tot_got)); ok <- FALSE
  }
  ok
}

# accumulators: (gene, cell_type) -> summed quantities across datasets
key <- function(g, ct) paste(g, ct, sep = "||")
acc <- new.env(parent = emptyenv())
for (g in GENES) for (ct in CELLTYPES) {
  acc[[key(g, ct)]] <- list(n_cells = 0, n_raw_pos = 0, n_dec_pos = 0,
                            n_dec_pos_strict = 0, sum_raw = 0, sum_dec = 0,
                            sum_cp10k_raw = 0, sum_cp10k_dec = 0)
}
contam_by_ct <- list(Hepatocytes = numeric(0), `Endothelial cells` = numeric(0))
per_dataset_rows <- list()

for (dataset in datasets) {
  dsdir <- file.path(WORK, dataset)
  info <- jsonlite::fromJSON(file.path(dsdir, "dims.json"), simplifyVector = FALSE)
  meta <- read.csv(file.path(dsdir, "meta.csv.gz"), stringsAsFactors = FALSE)
  cat(sprintf("\n[run] %s: reconstructing %d cells ...\n", dataset, info$ncells))
  m <- read_matrix(dsdir, info)
  stopifnot(ncol(m) == nrow(meta))
  if (!verify(m, info)) stop(sprintf("reconstruction checksum failed for %s", dataset))
  cat(sprintf("[run] %s: matrix %d x %d, checksum OK. running decontX ...\n",
              dataset, nrow(m), ncol(m)))

  z <- as.integer(meta$z)  # 1 = Hepatocytes, 2 = Endothelial cells
  set.seed(SEED)
  res <- celda::decontX(x = m, z = z, batch = NULL, seed = SEED, verbose = TRUE)
  dec <- res$decontXcounts             # genes x cells, decontaminated
  contam <- res$contamination          # per-cell contamination fraction
  rownames(dec) <- genes_all

  tot_raw <- Matrix::colSums(m)
  tot_dec <- Matrix::colSums(dec)
  tot_raw[tot_raw == 0] <- NA_real_
  tot_dec[tot_dec == 0] <- NA_real_

  for (ct in CELLTYPES) {
    cells <- which(meta$cell_type == ct)
    if (length(cells) == 0) next
    contam_by_ct[[ct]] <- c(contam_by_ct[[ct]], contam[cells])
    for (g in GENES) {
      rg <- as.numeric(m[g, cells])
      dg <- as.numeric(dec[g, cells])
      cp_raw <- rg / tot_raw[cells] * 1e4
      cp_dec <- dg / tot_dec[cells] * 1e4
      a <- acc[[key(g, ct)]]
      a$n_cells          <- a$n_cells + length(cells)
      a$n_raw_pos        <- a$n_raw_pos + sum(rg > 0)
      a$n_dec_pos        <- a$n_dec_pos + sum(dg > 0)
      a$n_dec_pos_strict <- a$n_dec_pos_strict + sum(dg > 0.5)
      a$sum_raw          <- a$sum_raw + sum(rg)
      a$sum_dec          <- a$sum_dec + sum(dg)
      a$sum_cp10k_raw    <- a$sum_cp10k_raw + sum(cp_raw, na.rm = TRUE)
      a$sum_cp10k_dec    <- a$sum_cp10k_dec + sum(cp_dec, na.rm = TRUE)
      acc[[key(g, ct)]]  <- a
      per_dataset_rows[[length(per_dataset_rows) + 1L]] <- data.frame(
        dataset = dataset, gene = g, cell_type = ct, n_cells = length(cells),
        raw_pct = 100 * mean(rg > 0), decont_pct = 100 * mean(dg > 0),
        decont_pct_strict = 100 * mean(dg > 0.5),
        raw_mean_cp10k = mean(cp_raw, na.rm = TRUE),
        decont_mean_cp10k = mean(cp_dec, na.rm = TRUE),
        ambient_fraction = ifelse(sum(rg) > 0, 1 - sum(dg) / sum(rg), NA_real_),
        stringsAsFactors = FALSE)
    }
  }
  cat(sprintf("[run] %s decontX contamination: median=%.3f mean=%.3f (per-cell)\n",
              dataset, median(contam), mean(contam)))
  rm(m, dec, res); gc()
}

# ---- pool across datasets (sum counts) ----
rows <- list()
for (g in GENES) for (ct in CELLTYPES) {
  a <- acc[[key(g, ct)]]
  if (a$n_cells == 0) next
  rows[[length(rows) + 1L]] <- data.frame(
    gene = g, cell_type = ct, n_cells = a$n_cells,
    raw_pct = 100 * a$n_raw_pos / a$n_cells,
    decont_pct = 100 * a$n_dec_pos / a$n_cells,
    decont_pct_strict = 100 * a$n_dec_pos_strict / a$n_cells,
    raw_mean_cp10k = a$sum_cp10k_raw / a$n_cells,
    decont_mean_cp10k = a$sum_cp10k_dec / a$n_cells,
    ambient_fraction = ifelse(a$sum_raw > 0, 1 - a$sum_dec / a$sum_raw, NA_real_),
    stringsAsFactors = FALSE)
}
final <- do.call(rbind, rows)
final <- final[order(match(final$gene, GENES), final$cell_type), ]
csv_path <- file.path(OUT, "decontx_ambient_correction.csv")
write.csv(final, csv_path, row.names = FALSE)
write.csv(do.call(rbind, per_dataset_rows),
          file.path(OUT, "decontx_ambient_correction_per_dataset.csv"), row.names = FALSE)
cat(sprintf("\n[run] wrote %s (%d rows)\n", csv_path, nrow(final)))

# ---- honest verdict summary ----
gv <- function(g, ct, col) {
  r <- final[final$gene == g & final$cell_type == ct, ]
  if (nrow(r) == 0) return(NA_real_) else return(r[[col]])
}
sink(file.path(OUT, "decontx_summary.txt"))
cat("decontX ambient-RNA correction — GLP-1RA incretin axis (human liver atlas)\n")
cat("==========================================================================\n")
cat(sprintf("Datasets pooled: %s\n", paste(datasets, collapse = ", ")))
cat(sprintf("Skipped (degenerate, <50 hep or <50 endo): %s\n",
            paste(names(manifest$skipped), collapse = ", ")))
cat("z = cell_type (Hepatocytes / Endothelial cells); batch = per-dataset loop.\n")
cat("Hep+Endo run together so hepatocyte (ALB) ambient contaminates endothelial.\n\n")

cat("Per-cell decontX contamination fraction (pooled):\n")
for (ct in CELLTYPES) {
  v <- contam_by_ct[[ct]]
  cat(sprintf("  %-18s n=%d  median=%.3f mean=%.3f  q10=%.3f q90=%.3f\n",
              ct, length(v), median(v), mean(v), quantile(v, .1), quantile(v, .9)))
}

fmt <- function(g, ct) sprintf(
  "  %-7s %-18s n=%-7d raw%%=%6.2f -> decont%%=%6.2f (strict%%=%6.2f)  cp10k %6.3f->%6.3f  ambient_frac=%.2f",
  g, ct, gv(g, ct, "n_cells"), gv(g, ct, "raw_pct"), gv(g, ct, "decont_pct"),
  gv(g, ct, "decont_pct_strict"), gv(g, ct, "raw_mean_cp10k"),
  gv(g, ct, "decont_mean_cp10k"), gv(g, ct, "ambient_fraction"))

cat("\nAXIS + CONTROLS (pooled, raw vs decontaminated):\n")
for (g in GENES) { for (ct in CELLTYPES) cat(fmt(g, ct), "\n") }

cat("\n---- VERDICT ----\n")
d_glp_hep  <- gv("GLP1R", "Hepatocytes", "raw_pct") - gv("GLP1R", "Hepatocytes", "decont_pct")
d_glp_endo <- gv("GLP1R", "Endothelial cells", "raw_pct") - gv("GLP1R", "Endothelial cells", "decont_pct")
cat(sprintf("GLP1R hepatocyte %%-expr: %.2f -> %.2f (Δ=%.2f pp; ambient_frac=%.2f)\n",
            gv("GLP1R","Hepatocytes","raw_pct"), gv("GLP1R","Hepatocytes","decont_pct"),
            d_glp_hep, gv("GLP1R","Hepatocytes","ambient_fraction")))
cat(sprintf("GLP1R endothelial %%-expr: %.2f -> %.2f (Δ=%.2f pp; ambient_frac=%.2f)\n",
            gv("GLP1R","Endothelial cells","raw_pct"), gv("GLP1R","Endothelial cells","decont_pct"),
            d_glp_endo, gv("GLP1R","Endothelial cells","ambient_fraction")))
cat(sprintf("POS-CTRL GCGR hep ambient_frac=%.2f (survive?%s); DPP4 hep ambient_frac=%.2f (survive?%s)\n",
            gv("GCGR","Hepatocytes","ambient_fraction"),
            ifelse(gv("GCGR","Hepatocytes","ambient_fraction") < 0.5, "YES", "NO-FLAG"),
            gv("DPP4","Hepatocytes","ambient_fraction"),
            ifelse(gv("DPP4","Hepatocytes","ambient_fraction") < 0.5, "YES", "NO-FLAG")))
cat(sprintf("POS-CTRL ALB endothelial ambient_frac=%.2f (drop if ambient? %s)\n",
            gv("ALB","Endothelial cells","ambient_fraction"),
            ifelse(gv("ALB","Endothelial cells","ambient_fraction") > 0.5, "YES-ambient-confirmed", "NO")))
cat(sprintf("POS-CTRL PECAM1 endothelial ambient_frac=%.2f; STAB2 endothelial ambient_frac=%.2f (should SURVIVE, <0.5)\n",
            gv("PECAM1","Endothelial cells","ambient_fraction"),
            gv("STAB2","Endothelial cells","ambient_fraction")))
sink()
cat("[run] wrote decontx_summary.txt\n")

# ---- tidy up large binaries (keep meta / dims / outputs) ----
bins <- list.files(WORK, pattern = "\\.bin$", recursive = TRUE, full.names = TRUE)
if (length(bins)) { file.remove(bins); cat(sprintf("[run] removed %d binary work files\n", length(bins))) }
cat("[run] DONE\n")
