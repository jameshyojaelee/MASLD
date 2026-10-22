#!/usr/bin/env Rscript
# N3 ambient check, step 3: library -> donor aggregation, fidelity check against N3, model inputs.
# Usage: 03_aggregate.R <RESULT_DIR> <N3_DIR>
suppressPackageStartupMessages({ library(data.table); library(jsonlite) })
argv <- commandArgs(trailingOnly = TRUE); stopifnot(length(argv) == 2L)
res <- argv[[1]]; n3 <- argv[[2]]
outs <- file.path(res, c("hep_donor_counts.tsv.gz", "library_ambient_qc.tsv", "donor_ambient.tsv",
                         "donor_gene_celltype_counts.tsv.gz", "donor_matrices.rds", "donor_meta.tsv"))
if (any(file.exists(outs))) stop("Refusing to overwrite aggregate outputs in ", res)

genes <- readLines(file.path(res, "genes.txt"))
report_genes <- readLines(file.path(res, "report_genes.txt"))
dm <- fread(file.path(n3, "donor_meta.tsv"))
lm <- fread(file.path(n3, "library_meta.tsv"))
inc <- dm[included == TRUE]
libs <- lm[donor %in% inc$donor]
stopifnot(nrow(inc) == 70L, nrow(libs) == 208L)
files <- file.path(res, "per_library", paste0(libs$library, ".rds"))
if (!all(file.exists(files))) stop("missing per-library results: ", paste(libs$library[!file.exists(files)], collapse = ","))

donors <- sort(inc$donor)
mk <- function() matrix(0, length(genes), length(donors), dimnames = list(genes, donors))
hep_raw <- mk(); hep_dec <- mk(); bg <- mk()
qc <- list(); long <- list()
for (k in seq_len(nrow(libs))) {
  r <- readRDS(files[k]); d <- libs$donor[k]
  stopifnot(r$library == libs$library[k], length(r$hep_raw) == length(genes))
  hep_raw[, d] <- hep_raw[, d] + r$hep_raw
  hep_dec[, d] <- hep_dec[, d] + r$hep_dec
  if (r$profile_source == "empty_droplets") bg[, d] <- bg[, d] + r$background
  st <- r$export_stats
  cs <- as.data.table(r$ct_summary)
  h <- cs[cell_type == "Hepatocytes"]
  qc[[k]] <- data.table(
    library = r$library, dataset = r$dataset, donor = d, status = r$status, profile_source = r$profile_source,
    n_cells = st$n_cells, n_hep_cells = h$n_cells, n_clusters = r$n_clusters,
    pooled_small_types = paste(r$pooled_small_types, collapse = ";"),
    n_raw_barcodes = st$n_raw_barcodes, n_cellranger_called = st$n_cellranger_called,
    atlas_cells_in_raw_matrix = st$atlas_cells_in_raw_matrix,
    atlas_cells_total_equal_raw_matrix = st$atlas_cells_total_equal_raw_matrix,
    n_empty_droplets = st$n_empty_droplets, empty_umis_atlas_genes = st$empty_umis_atlas_genes,
    raw_matrix_umis = st$raw_matrix_umis, called_cell_umis = st$called_cell_umis,
    frac_umis_outside_called_cells = 1 - st$called_cell_umis / st$raw_matrix_umis,
    hep_raw_umis = h$raw_umis, hep_dec_umis = h$dec_umis,
    hep_contamination = 1 - h$dec_umis / h$raw_umis,
    all_cells_contamination = 1 - sum(cs$dec_umis) / sum(cs$raw_umis),
    decontx_minutes = r$decontx_minutes, seed = r$seed, celda = r$celda_version)
  a <- as.data.table(as.table(r$ct_raw)); setnames(a, c("symbol", "cell_type", "raw"))
  b <- as.data.table(as.table(r$ct_dec)); setnames(b, c("symbol", "cell_type", "dec"))
  ab <- merge(a, b, by = c("symbol", "cell_type"))
  ab <- merge(ab, cs[, .(cell_type, n_cells, ct_raw_umis = raw_umis, ct_dec_umis = dec_umis)], by = "cell_type")
  ab[, donor := d]
  long[[k]] <- ab
}
qc <- rbindlist(qc)
fwrite(qc, file.path(res, "library_ambient_qc.tsv"), sep = "\t")
long <- rbindlist(long)
dl <- long[, .(raw = sum(raw), dec = sum(dec)), by = .(donor, cell_type, symbol)]
# cell-type totals per donor: sum over libraries (one row per library x cell type x gene; take one gene)
tot <- long[symbol == report_genes[1], .(n_cells = sum(n_cells), ct_raw_umis = sum(ct_raw_umis), ct_dec_umis = sum(ct_dec_umis)),
            by = .(donor, cell_type)]
dl <- merge(dl, tot, by = c("donor", "cell_type"))
fwrite(dl, file.path(res, "donor_gene_celltype_counts.tsv.gz"), sep = "\t")

# ------------------------------------------------ fidelity against N3 pseudobulk
old <- fread(file.path(n3, "hep_donor_counts.tsv.gz"))
stopifnot(identical(old$symbol, genes))
old_m <- as.matrix(old[, donors, with = FALSE]); storage.mode(old_m) <- "double"
max_abs_diff <- max(abs(old_m - hep_raw))
fid <- data.table(n_genes = length(genes), n_donors = length(donors), max_abs_diff_raw_vs_N3 = max_abs_diff,
                  n_cells_differ = sum(old_m != hep_raw), total_umis_N3 = sum(old_m), total_umis_here = sum(hep_raw))
fwrite(fid, file.path(res, "hep_donor_counts_raw_check.tsv"), sep = "\t")
print(fid)
if (max_abs_diff != 0) stop("Uncorrected hepatocyte donor sums differ from N3; stopping (PRESPEC fidelity check).")

out <- data.table(symbol = genes, ensembl = old$ensembl, round(hep_dec, 3))
fwrite(out, file.path(res, "hep_donor_counts.tsv.gz"), sep = "\t")

dq <- qc[, .(n_libraries = .N, n_passthrough = sum(status != "corrected"),
             n_no_empty_profile = sum(profile_source != "empty_droplets"),
             hep_raw_umis = sum(hep_raw_umis), hep_dec_umis = sum(hep_dec_umis),
             raw_matrix_umis = sum(raw_matrix_umis), called_cell_umis = sum(called_cell_umis),
             n_empty_droplets = sum(n_empty_droplets), empty_umis = sum(empty_umis_atlas_genes)), by = donor]
dq[, hep_contamination := 1 - hep_dec_umis / hep_raw_umis]
dq[, frac_umis_outside_called_cells := 1 - called_cell_umis / raw_matrix_umis]
dq <- merge(inc[, .(donor, dataset, disease_stage_coarse, stage_ord, sex_call, n_hep_cells)], dq, by = "donor")
stopifnot(all(abs(dq$hep_raw_umis - colSums(hep_raw)[dq$donor]) < 0.5))
fwrite(dq, file.path(res, "donor_ambient.tsv"), sep = "\t")
saveRDS(list(hep_raw = hep_raw, hep_dec = hep_dec, background = bg), file.path(res, "donor_matrices.rds"))

# inputs that 03_donor_model.R reads from its RESULT_DIR (unchanged copies of the N3 files)
for (f in c("donor_meta.tsv", "target_list.tsv", "target_celltype_summary.tsv")) {
  stopifnot(file.copy(file.path(n3, f), file.path(res, f), overwrite = FALSE))
}
cat("libraries:", nrow(qc), " corrected:", sum(qc$status == "corrected"),
    " empty-droplet profile:", sum(qc$profile_source == "empty_droplets"),
    " donors with a passthrough library:", sum(dq$n_passthrough > 0), "\n")
writeLines(capture.output(sessionInfo()), file.path(res, "sessionInfo_aggregate.txt"))
