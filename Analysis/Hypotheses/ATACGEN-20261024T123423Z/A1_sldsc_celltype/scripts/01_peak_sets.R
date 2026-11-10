#!/usr/bin/env Rscript
# A1 step 1: ATAC v3 lineage peak sets for S-LDSC (adapted from H22 01_peak_sets.R).
#
# Inputs: v3 consensus peak manifest (hg38, 500 bp, five lineages), corrected v3
# DA table (MASH minus normal, per cohort), v3 donor-by-peak counts (for
# accessibility), GENCODE v49 GTF (TSS), BSgenome hg38 (GC), UCSC hg38ToHg19.
#
# Output: <OUT>/peaks/peaks.tsv, one row per kept lineage peak, with hg19
# coordinates (1-based closed), set flags (only, prom, tested, gain/lost
# tertiles, 10% tails) and per-peak values for tested peaks (z_comb,
# accessibility, GC, stratum, standardized continuous annotations).
# Plus peaks/liftover_report.tsv and peaks/set_sizes_hg38.tsv.
#
# Usage: 01_peak_sets.R <OUT> <V3_DIR> <V3_CORR_DIR> <GTF_GZ> <CHAIN>
suppressPackageStartupMessages({
  library(data.table); library(GenomicRanges); library(rtracklayer)
  library(BSgenome.Hsapiens.UCSC.hg38); library(Biostrings); library(Matrix)
})
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 5)
OUT <- args[1]; V3 <- args[2]; CORR <- args[3]; GTF <- args[4]; CHAIN <- args[5]
LINEAGES <- c("cholangiocyte", "hepatocyte", "macrophage", "stellate", "t_nk")
DA_LINEAGES <- c("hepatocyte", "stellate")
COHORTS <- c("GSE244832", "GSE281367")
TSS_FLANK <- 1000L
dir.create(file.path(OUT, "peaks"), recursive = TRUE, showWarnings = FALSE)

# ---- consensus peaks ---------------------------------------------------------
man <- fread(file.path(V3, "consensus_peak_manifest.tsv"),
             select = c("lineage", "peak_id", "chrom", "start0", "end", "blacklist_overlap"))
stopifnot(setequal(unique(man$lineage), LINEAGES))
auto <- paste0("chr", 1:22)
cat("manifest peaks:", nrow(man), "\n")
man[, autosome := chrom %in% auto]
man[, blacklist_overlap := as.logical(blacklist_overlap)]
pk <- man[autosome == TRUE & blacklist_overlap == FALSE]
pk[, coord := paste0(chrom, ":", start0, "-", end)]
stopifnot(!anyDuplicated(pk[, paste(lineage, coord)]))
cat("kept (chr1-22, not blacklisted):", nrow(pk), "\n")
gr <- GRanges(pk$chrom, IRanges(pk$start0 + 1L, pk$end))

# ---- lineage-only: no overlap with any other lineage's kept peak ------------
pk[, only := FALSE]
for (L in LINEAGES) {
  i <- which(pk$lineage == L)
  j <- which(pk$lineage != L)
  hit <- overlapsAny(gr[i], gr[j])
  pk[i, only := !hit]
}

# ---- promoter: overlaps TSS +/- 1 kb of a protein_coding or lncRNA transcript --
awk <- paste0("zcat ", shQuote(GTF), " | awk -F'\\t' '$3==\"transcript\" && $1 ~ /^chr([1-9]|1[0-9]|2[0-2])$/'",
              " | grep -E 'gene_type \"(protein_coding|lncRNA)\"' | cut -f1,4,5,7")
tx <- fread(cmd = awk, header = FALSE, col.names = c("chr", "s", "e", "strand"))
tx[, tss := fifelse(strand == "+", s, e)]
tx <- unique(tx[, .(chr, tss)])
cat("unique protein_coding/lncRNA TSS on chr1-22:", nrow(tx), "\n")
tss_gr <- GRanges(tx$chr, IRanges(pmax(1L, tx$tss - TSS_FLANK), tx$tss + TSS_FLANK))
pk[, prom := overlapsAny(gr, tss_gr)]

# ---- DA-tested peaks: combined MASH-minus-normal z -----------------------------
da <- fread(file.path(CORR, "da/da_peak_results.tsv.gz"),
            select = c("lineage", "peak_coordinate", "logFC_gse244832", "pvalue_gse244832",
                       "logFC_gse281367", "pvalue_gse281367", "standard_error_gse244832",
                       "standard_error_gse281367"))
stopifnot(setequal(unique(da$lineage), DA_LINEAGES))
zf <- function(lfc, p) sign(lfc) * qnorm(p / 2, lower.tail = FALSE)
da[, z244832 := zf(logFC_gse244832, pvalue_gse244832)]
da[, z281367 := zf(logFC_gse281367, pvalue_gse281367)]
da[, z_comb := (z244832 + z281367) / sqrt(2)]
# check: z from p agrees with logFC / corrected SE (moderated t; df differ, so only rank agreement)
chk <- da[, .(spearman_z_vs_t_244832 = cor(z244832, logFC_gse244832 / standard_error_gse244832, method = "spearman"),
              spearman_z_vs_t_281367 = cor(z281367, logFC_gse281367 / standard_error_gse281367, method = "spearman"),
              n = .N), by = lineage]
print(chk)
fwrite(chk, file.path(OUT, "peaks/da_z_check.tsv"), sep = "\t")
stopifnot(all(is.finite(da$z_comb)))
setnames(da, "peak_coordinate", "coord")
pk <- merge(pk, da[, .(lineage, coord, z244832, z281367, z_comb)], by = c("lineage", "coord"), all.x = TRUE, sort = FALSE)
pk[, tested := !is.na(z_comb)]
n_da_unmatched <- da[!pk[tested == TRUE], on = c("lineage", "coord"), .N]
cat("DA peaks not among kept consensus peaks (blacklist or non-autosomal):", n_da_unmatched, "\n")

# ---- accessibility of tested peaks: mean log2(CPM + 0.5) over DA donors ------
qc <- fread(file.path(V3, "da/da_model_qc.tsv"))
qc <- qc[model == "primary_condition_only"]
pk[, access := NA_real_]
for (L in DA_LINEAGES) {
  per_cohort <- list()
  for (C in COHORTS) {
    pre <- file.path(V3, "counts", C, L)
    don <- fread(paste0(pre, ".donors.tsv"))
    pks <- fread(paste0(pre, ".peaks.tsv"))
    m <- readMM(gzcon(file(paste0(pre, ".mtx.gz"), "rb")))
    stopifnot(nrow(m) == nrow(don), ncol(m) == nrow(pks))
    q <- qc[cohort == C & lineage == L]
    ri <- match(q$donor_id, don$donor_id)
    stopifnot(!anyNA(ri))
    ci <- match(pk[lineage == L & tested == TRUE, coord], pks$peak_coordinate)
    stopifnot(!anyNA(ci))
    x <- as.matrix(m[ri, ci, drop = FALSE])
    lcpm <- log2(sweep(x, 1, q$raw_library_size, "/") * 1e6 + 0.5)
    per_cohort[[C]] <- colMeans(lcpm)
  }
  pk[lineage == L & tested == TRUE, access := (per_cohort[[1]] + per_cohort[[2]]) / 2]
}

# ---- GC of tested peaks ----------------------------------------------------------
pk[, gc := NA_real_]
ti <- which(pk$tested)
gr_t <- GRanges(pk$chrom[ti], IRanges(pk$start0[ti] + 1L, pk$end[ti]))
sq <- getSeq(BSgenome.Hsapiens.UCSC.hg38, gr_t)
pk[ti, gc := as.numeric(letterFrequency(sq, "GC", as.prob = TRUE))]

# ---- strata and rank-based sets within tested peaks --------------------------------
ntile5 <- function(x) as.integer(ceiling(5 * frank(x, ties.method = "first") / length(x)))
setorder(pk, lineage, chrom, start0)
pk[, `:=`(stratum = NA_integer_, gain = FALSE, lost = FALSE, gain10 = FALSE, lost10 = FALSE)]
for (L in DA_LINEAGES) {
  i <- which(pk$lineage == L & pk$tested)
  s <- (ntile5(pk$access[i]) - 1L) * 5L + ntile5(pk$gc[i])
  pk[i, stratum := s]
}
pk[tested == TRUE, r := frank(z_comb, ties.method = "first"), by = .(lineage, stratum)]
pk[tested == TRUE, n_s := .N, by = .(lineage, stratum)]
pk[tested == TRUE, `:=`(lost = r <= floor(n_s / 3), gain = r > n_s - floor(n_s / 3),
                        lost10 = r <= floor(0.1 * n_s), gain10 = r > n_s - floor(0.1 * n_s))]
std <- function(x) (x - mean(x)) / sd(x)
pk[, `:=`(resp_std = NA_real_, sign_std = NA_real_, access_std = NA_real_, gc_std = NA_real_)]
pk[tested == TRUE, `:=`(resp_std = std(abs(z_comb)), sign_std = std(z_comb),
                        access_std = std(access), gc_std = std(gc)), by = lineage]

# ---- liftover hg38 -> hg19 (H22 rule) ----------------------------------------------
gr38 <- GRanges(pk$chrom, IRanges(pk$start0 + 1L, pk$end))
names(gr38) <- paste(pk$lineage, pk$coord, sep = "|")
ch <- import.chain(CHAIN)
lo <- unlist(liftOver(gr38, ch))
lo <- data.table(pkey = names(lo), chr19 = as.character(seqnames(lo)), s19 = start(lo), e19 = end(lo))
lo <- lo[, .(n_chr = uniqueN(chr19), chr19 = chr19[1], s19 = min(s19), e19 = max(e19)), by = pkey]
pk[, pkey := paste(lineage, coord, sep = "|")]
pk <- merge(pk, lo, by = "pkey", all.x = TRUE, sort = FALSE)
pk[, lifted := !is.na(chr19) & n_chr == 1L & chr19 == chrom & (e19 - s19 + 1L) <= 2L * (end - start0)]
rep <- pk[, .(total = .N, mapped_any = sum(!is.na(chr19)), lifted_kept = sum(lifted)), by = lineage]
fwrite(rep, file.path(OUT, "peaks/liftover_report.tsv"), sep = "\t")
print(rep)

setorder(pk, lineage, chrom, start0)
pk[, chr19n := as.integer(sub("chr", "", chr19))]
out <- pk[, .(lineage, coord, chrom, start0, end, lifted, chr19 = chr19n, s19, e19,
              only, prom, tested, gain, lost, gain10, lost10, stratum,
              z244832, z281367, z_comb, access, gc, resp_std, sign_std, access_std, gc_std)]
fwrite(out, file.path(OUT, "peaks/peaks.tsv"), sep = "\t", na = "NA")

sz <- rbind(
  pk[, .(set = "ALL", n_hg38 = .N, n_lifted = sum(lifted)), by = lineage],
  pk[only == TRUE, .(set = "ONLY", n_hg38 = .N, n_lifted = sum(lifted)), by = lineage],
  pk[only == TRUE & prom == TRUE, .(set = "PROM", n_hg38 = .N, n_lifted = sum(lifted)), by = lineage],
  pk[only == TRUE & prom == FALSE, .(set = "DIST", n_hg38 = .N, n_lifted = sum(lifted)), by = lineage],
  pk[tested == TRUE, .(set = "TESTED", n_hg38 = .N, n_lifted = sum(lifted)), by = lineage],
  pk[gain == TRUE, .(set = "GAIN", n_hg38 = .N, n_lifted = sum(lifted)), by = lineage],
  pk[lost == TRUE, .(set = "LOST", n_hg38 = .N, n_lifted = sum(lifted)), by = lineage],
  pk[gain10 == TRUE, .(set = "GAIN10", n_hg38 = .N, n_lifted = sum(lifted)), by = lineage],
  pk[lost10 == TRUE, .(set = "LOST10", n_hg38 = .N, n_lifted = sum(lifted)), by = lineage])
fwrite(sz, file.path(OUT, "peaks/set_sizes_hg38.tsv"), sep = "\t")
print(sz)
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo_01.txt"))
cat("A1 step 1 done\n")
