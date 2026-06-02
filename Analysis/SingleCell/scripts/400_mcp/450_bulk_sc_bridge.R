#!/usr/bin/env Rscript
# 450_bulk_sc_bridge.R
# Bridge bulk k=6 NMF programs to scRNA cNMF gene spectra.
# Correlate bulk W (genes x 6) with each scRNA cNMF program gene-spectra vector.
# Identify "scRNA-unique" programs where max |r| < 0.3 to any bulk program.

suppressPackageStartupMessages({
  library(data.table)
  library(NMF)
})

root <- Sys.getenv("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
# NMF_CACHE_PATH env var lets Bravo remediation variants drive this bridge
# against alternate bulk W matrices (e.g. protonly / nonprotonly) without
# clobbering canonical outputs.
bulk_cache <- Sys.getenv("NMF_CACHE_PATH",
                          file.path(root, "RNA-seq/results/subtypes/nmf_results_cache_clean.rds"))
# Output suffix derived from cache basename (empty for canonical clean).
.cache_base <- sub("\\.rds$", "", basename(bulk_cache))
SUFFIX <- sub("^nmf_results_cache_clean", "", .cache_base)  # "" or "_protonly" / "_nonprotonly"
BULK_K_TAG_LBL <- if (as.integer(Sys.getenv("NMF_BULK_K", "6")) == 6L) "" else
                  sprintf("_k%s", Sys.getenv("NMF_BULK_K", "6"))
bulk_labels_f <- file.path(root,
  sprintf("RNA-seq/results/subtypes/program_labels%s%s.csv", SUFFIX, BULK_K_TAG_LBL))
cnmf_spectra_dir <- file.path(root, "Analysis/SingleCell/results_gpu_v2/mcp/cnmf_runs")
out_dir <- file.path(root, "Analysis/SingleCell/results_gpu_v2/mcp/integration")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Bulk k=6 W
BULK_K <- as.integer(Sys.getenv("NMF_BULK_K", "6"))
bulk <- readRDS(bulk_cache)
if (!is.null(bulk$nmf_results) && !is.null(bulk$nmf_results[[as.character(BULK_K)]])) {
  W_bulk <- basis(bulk$nmf_results[[as.character(BULK_K)]])
} else {
  stop(sprintf("k=%d not in bulk cache %s", BULK_K, bulk_cache))
}
rownames(W_bulk) <- sub("\\..*$", "", rownames(W_bulk))  # ensembl_base

# Labels file may have leading '#' comment lines (written by 95c_relabel_programs.R).
labels <- NULL
if (file.exists(bulk_labels_f)) {
  lines <- readLines(bulk_labels_f)
  lines <- lines[!grepl("^#", lines)]
  labels <- fread(text = paste(lines, collapse = "\n"))
}
# biological_label (95c) takes precedence; fall back to legacy column name if present.
label_col <- NULL
if (!is.null(labels)) {
  if ("biological_label" %in% names(labels)) label_col <- "biological_label"
  else if ("program_label" %in% names(labels)) label_col <- "program_label"
}
if (!is.null(label_col) && nrow(labels) == ncol(W_bulk)) {
  # Prepend letter code so we can still reference bulk programs by A/B/C (historical)
  # while carrying the biology label in the same string.
  colnames(W_bulk) <- paste0(LETTERS[seq_len(ncol(W_bulk))], "=",
                             labels[[label_col]][seq_len(ncol(W_bulk))])
} else {
  colnames(W_bulk) <- LETTERS[seq_len(ncol(W_bulk))]
}
cat(sprintf("[450] bulk cache: %s (k=%d, suffix='%s')\n", bulk_cache, BULK_K, SUFFIX))
cat("[450] bulk program labels: ", paste(colnames(W_bulk), collapse = " | "), "\n")

# Load gencode metadata to map symbol <-> ensembl
gene_meta <- fread(file.path(root, "data/gencode_v49_gene_metadata.tsv.gz"))
sym_to_ens <- setNames(gene_meta$ensembl_base, gene_meta$gene_name)

# Iterate cNMF runs + k values available
args <- commandArgs(trailingOnly = TRUE)
name_k <- if (length(args) >= 2) args[1:2] else c("global", "20")
name <- name_k[1]; k <- as.integer(name_k[2])
spectra_f <- file.path(cnmf_spectra_dir, name, name,
                       sprintf("%s.gene_spectra_score.k_%d.dt_0_03.txt", name, k))
if (!file.exists(spectra_f)) {
  spectra_f <- file.path(cnmf_spectra_dir, name,
                         sprintf("%s.gene_spectra_score.k_%d.dt_0_03.txt", name, k))
}
if (!file.exists(spectra_f)) stop(sprintf("spectra file not found: %s", spectra_f))
spectra <- as.matrix(fread(spectra_f), rownames = 1)  # programs x genes (symbols)

# Map cNMF genes (symbols) to ensembl
cnmf_gene_names <- colnames(spectra)
cnmf_ens <- sym_to_ens[cnmf_gene_names]
keep <- !is.na(cnmf_ens)
spectra <- spectra[, keep]
colnames(spectra) <- cnmf_ens[keep]

# Intersect with bulk W rownames
common <- intersect(colnames(spectra), rownames(W_bulk))
cat(sprintf("[450] common ensembl: %d\n", length(common)))
W_bulk <- W_bulk[common, ]
spectra <- spectra[, common]

# Pearson r matrix: cNMF programs x bulk programs
R <- cor(t(spectra), W_bulk, method = "pearson", use = "pairwise.complete.obs")
out <- as.data.frame(as.table(R))
setDT(out)
setnames(out, c("cnmf_program", "bulk_program", "pearson_r"))
# Output filename carries BOTH the scRNA-cNMF k (positional arg `k`) AND the
# bulk NMF k (BULK_K), so protonly k=6 and protonly k=8 bridges don't collide.
# When BULK_K is the default 6 we keep the legacy unsuffixed filename for
# backward compatibility.
bulk_k_tag <- if (BULK_K == 6L) "" else sprintf("_bulkk%d", BULK_K)
long_out <- file.path(out_dir,
  sprintf("bulk_sc_bridge_%s_k%d%s%s.tsv", name, k, SUFFIX, bulk_k_tag))
fwrite(out, long_out, sep = "\t")

# Per-cnmf-program max absolute correlation
summary <- out[, .(max_abs_r = max(abs(pearson_r)),
                   max_bulk_program = bulk_program[which.max(abs(pearson_r))]),
               by = cnmf_program]
summary[, scRNA_unique := max_abs_r < 0.3]
sum_out <- file.path(out_dir,
  sprintf("bulk_sc_bridge_%s_k%d%s%s_summary.tsv", name, k, SUFFIX, bulk_k_tag))
fwrite(summary, sum_out, sep = "\t")

cat(sprintf("[450] wrote %s\n", long_out))
cat(sprintf("[450] wrote %s\n", sum_out))
cat(sprintf("[450] scRNA-unique programs: %d / %d\n",
            sum(summary$scRNA_unique), nrow(summary)))
