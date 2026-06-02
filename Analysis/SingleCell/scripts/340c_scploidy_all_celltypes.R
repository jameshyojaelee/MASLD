#!/usr/bin/env Rscript
# 340c_scploidy_all_celltypes.R
# Re-run scPloidy on ALL cell types per donor (not just hepatocytes).
# Rationale: scPloidy's mixture EM needs a diploid reference. Restricting to hepatocytes
# (highly polyploid) starves the model of 2N anchor cells, biasing calls upward.
# Including non-hepatocyte cells (~2N in normal liver biology) lets the model calibrate.
#
# Inputs: full scATAC barcode table (all cell types passing label transfer), fragments
# Outputs: scploidy_percell_allCT_<donor>.csv.gz, scploidy_celltype_summary.csv

suppressPackageStartupMessages({
  library(scPloidy)
  library(GenomicRanges)
  library(IRanges)
  library(data.table)
  library(readr)
})

HG38_AUTO_LEN <- c(
  chr1=248956422,  chr2=242193529,  chr3=198295559,  chr4=190214555,
  chr5=181538259,  chr6=170805979,  chr7=159345973,  chr8=145138636,
  chr9=138394717,  chr10=133797422, chr11=135086622, chr12=133275309,
  chr13=114364328, chr14=107043718, chr15=101991189, chr16=90338345,
  chr17=83257441,  chr18=80373285,  chr19=58617616,  chr20=64444167,
  chr21=46709983,  chr22=50818468
)

PROJECT_ROOT <- Sys.getenv("MASLD_PROJECT_ROOT",
                           "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
setwd(PROJECT_ROOT)

FRAG_DIR    <- "Analysis/ATAC/Human_Multiome/results/fragments"
ALL_BC_TSV  <- "Analysis/SingleCell/results_gpu_v2/ploidy/all_celltype_barcodes.tsv.gz"
SIMPLE_REPEAT <- "data/ucsc_hg38/simpleRepeat.hg38.txt.gz"
OUT_DIR     <- "Analysis/SingleCell/results_gpu_v2/ploidy"
dir.create(file.path(OUT_DIR, "allCT"), showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "qc_allCT"), showWarnings = FALSE)

DONOR_ARG <- Sys.getenv("DONOR_ID", "")
PLOIDY_LEVELS <- c(2, 4, 8)
TN5_OFFSET    <- c(1, 0)
MIN_CELLS     <- 200

autosomes <- names(HG38_AUTO_LEN)
target_regions <- GRanges(seqnames = autosomes,
                          ranges   = IRanges(start = 1L, end = HG38_AUTO_LEN))

cat("[", format(Sys.time()), "] Loading simpleRepeat\n", sep = "")
sr <- fread(cmd = paste("zcat", shQuote(SIMPLE_REPEAT)),
            header = FALSE, sep = "\t", select = c(2, 3, 4),
            col.names = c("chrom", "chromStart", "chromEnd"))
sr[, chromStart := chromStart + 1L]
sr <- sr[chrom %in% autosomes]
simpleRepeat <- GenomicRanges::makeGRangesFromDataFrame(
  sr, seqnames.field = "chrom", start.field = "chromStart", end.field = "chromEnd"
)
cat("  simpleRepeat regions:", length(simpleRepeat), "\n")

bc <- fread(ALL_BC_TSV)
stopifnot(all(c("barcode", "donor_id", "condition", "cell_type") %in% names(bc)))
donor_list <- if (nzchar(DONOR_ARG)) DONOR_ARG else sort(unique(bc$donor_id))

run_one_donor <- function(donor) {
  cat("\n=== Donor", donor, "===\n")
  bc_d <- bc[donor_id == donor]
  cat("  barcodes:", nrow(bc_d), "cells across",
      length(unique(bc_d$cell_type)), "cell types\n")
  if (nrow(bc_d) < MIN_CELLS) {
    return(data.frame(donor_id = donor, status = "too_few_cells",
                      n_cells = nrow(bc_d), stringsAsFactors = FALSE))
  }
  frag_in <- file.path(FRAG_DIR, paste0(donor, "_fragments.tsv.gz"))
  if (!file.exists(frag_in)) {
    return(data.frame(donor_id = donor, status = "missing_fragments",
                      stringsAsFactors = FALSE))
  }
  cat("[", format(Sys.time()), "] fragmentoverlapcount() on", nrow(bc_d), "cells\n", sep = "")
  fo <- tryCatch(
    fragmentoverlapcount(
      file = frag_in, targetregions = target_regions,
      excluderegions = simpleRepeat,
      targetbarcodes = bc_d$barcode, Tn5offset = TN5_OFFSET
    ),
    error = function(e) { cat("  ERROR fo:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(fo)) return(data.frame(donor_id = donor, status = "fail_fo", stringsAsFactors = FALSE))
  cat("  fo done, n=", nrow(fo), "\n")
  cat("[", format(Sys.time()), "] ploidy() ...\n", sep = "")
  pl <- tryCatch(ploidy(fo, PLOIDY_LEVELS),
                 error = function(e) { cat("  ERROR ploidy:", conditionMessage(e), "\n"); NULL })
  if (is.null(pl)) return(data.frame(donor_id = donor, status = "fail_em",
                                     stringsAsFactors = FALSE))

  # Annotate with cell_type + nfrags from fragmentoverlap
  pl <- as.data.table(pl)
  fo_dt <- as.data.table(fo)[, .(barcode, nfrags = as.numeric(nfrags))]
  pl <- merge(pl, bc_d[, .(barcode, cell_type, condition)], by = "barcode", all.x = TRUE)
  pl <- merge(pl, fo_dt, by = "barcode", all.x = TRUE)
  pl[, donor_id := donor]
  out_csv <- file.path(OUT_DIR, "allCT", paste0("scploidy_percell_allCT_", donor, ".csv.gz"))
  fwrite(pl, out_csv, compress = "gzip")
  cat("  wrote", out_csv, "\n")

  # Per-celltype summary
  ct_tab <- pl[, .(n = .N,
                   frac_2N = mean(ploidy.moment == 2),
                   frac_4N = mean(ploidy.moment == 4),
                   frac_8N = mean(ploidy.moment == 8),
                   mean_frac_ploidy = mean(ploidy.momentfractional, na.rm = TRUE),
                   median_frac_ploidy = median(ploidy.momentfractional, na.rm = TRUE),
                   median_nfrags = median(nfrags, na.rm = TRUE)),
               by = cell_type]
  ct_tab[, donor_id := donor]
  ct_tab[, condition := unique(pl$condition)[1]]
  ct_tab[, status := "ok"]

  # Quick QC plot: per-celltype ploidy distribution
  pdf(file.path(OUT_DIR, "qc_allCT", paste0(donor, "_allCT_diagnostic.pdf")),
      width = 10, height = 6)
  par(mfrow = c(2, 1), mar = c(8, 4, 2, 1))
  cts <- sort(unique(pl$cell_type))
  prop_mat <- sapply(cts, function(ct) {
    t <- table(factor(pl[cell_type == ct, ploidy.moment], levels = PLOIDY_LEVELS))
    if (sum(t) == 0) return(c(0,0,0)); as.numeric(t / sum(t))
  })
  rownames(prop_mat) <- paste0(PLOIDY_LEVELS, "N")
  barplot(prop_mat, las = 2, col = c("#4d9221","#e7298a","#7570b3"),
          main = paste(donor, "ploidy by cell type"), cex.names = 0.7,
          ylim = c(0, 1), ylab = "fraction", legend.text = TRUE,
          args.legend = list(x = "topright", cex = 0.7))
  # Sample sizes
  n_per_ct <- pl[, .N, by = cell_type]
  barplot(setNames(n_per_ct$N, n_per_ct$cell_type), las = 2, cex.names = 0.7,
          main = "cells per cell type", ylab = "n", log = "y")
  dev.off()

  ct_tab[]
}

results <- lapply(donor_list, run_one_donor)
ct_results <- do.call(rbind, lapply(results, function(x) {
  if (!"cell_type" %in% names(x)) {
    stat <- if (!is.null(x$status)) as.character(x$status) else "unknown"
    data.frame(donor_id = x$donor_id, cell_type = NA, n = NA,
               frac_2N = NA, frac_4N = NA, frac_8N = NA,
               mean_frac_ploidy = NA, median_frac_ploidy = NA, median_nfrags = NA,
               condition = NA, status = stat, stringsAsFactors = FALSE)
  } else as.data.frame(x)
}))
sum_path <- file.path(OUT_DIR, "scploidy_celltype_summary.csv")
if (file.exists(sum_path) && nzchar(DONOR_ARG)) {
  prev <- fread(sum_path)
  prev <- prev[!(donor_id %in% ct_results$donor_id)]
  ct_results <- rbind(prev, ct_results, fill = TRUE)
}
fwrite(ct_results, sum_path)
cat("\n[", format(Sys.time()), "] DONE\n", sep = "")
print(head(ct_results, 30))
