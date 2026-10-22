#!/usr/bin/env Rscript
# H22 step 1: donor-level peak sets for S-LDSC from GSE296875 multiome ATAC.
#
# Inputs: donor-by-lineage pseudobulk fragment counts in the authors' 306,706
# peaks (built for H24 from the authors' peaks assay), donor phenotypes, UCSC
# hg38ToHg19 chain. Units with >= 20 nuclei only; donors are the unit.
#
# Outputs under <OUT>/peaks and <OUT>/da:
#   da/<contrast>_<lineage>.tsv.gz   limma-voom tables (one row per tested peak)
#   da/specificity_<lineage>.tsv.gz  lineage-vs-rest moderated t
#   peaks/peaks_hg19.tsv             lifted autosomal peaks (row order of membership)
#   peaks/membership.u8              uint8, annotation-major (n_annot x n_peak)
#   peaks/annotations.tsv            annotation names, type, lineage, sizes
#   peaks/liftover_report.tsv
#
# Usage: 01_peak_sets.R <OUT_DIR> <H24_GSE296875_DIR> <DONOR_ENDPOINTS_TSV> <CHAIN>
suppressPackageStartupMessages({
  library(data.table); library(edgeR); library(limma)
  library(GenomicRanges); library(rtracklayer)
  library(BSgenome.Hsapiens.UCSC.hg38); library(Biostrings)
})
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 4)
OUT <- args[1]; PB <- args[2]; ENDPOINTS <- args[3]; CHAIN <- args[4]
SEED <- 20260923L; K_NULL <- 50L; TOP_FRAC <- 0.10; MIN_NUCLEI <- 20L
MIN_PER_GROUP <- 8L; MIN_PEAKS <- 1000L
dir.create(file.path(OUT, "peaks"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "da"), recursive = TRUE, showWarnings = FALSE)

# ---- counts --------------------------------------------------------------
dims <- jsonlite::fromJSON(file.path(PB, "dims.json"))$atac
peaks <- fread(file.path(PB, "peaks.tsv"))$peak_id
units <- fread(file.path(PB, "units.tsv"), colClasses = list(character = "donor_id"))
stopifnot(length(peaks) == dims[1], nrow(units) == dims[2])
con <- file(file.path(PB, "atac_peak_by_unit.f32"), "rb")
cnt <- readBin(con, "numeric", n = dims[1] * dims[2], size = 4, endian = "little")
close(con)
cnt <- matrix(cnt, nrow = dims[1], ncol = dims[2])
stopifnot(all(abs(colSums(cnt) - units$atac_in_peak_fragments) < 1))
rownames(cnt) <- peaks
colnames(cnt) <- paste(units$donor_id, units$lineage_id, sep = "|")

ep <- fread(ENDPOINTS, colClasses = list(character = "donor_id"))
ep[, steatotic := fifelse(steatosis_observed, steatosis_numeric >= 5, NA)]
ep[, fibrosis := fifelse(fibrosis_observed, as.logical(fibrosis_any), NA)]
units <- merge(units, ep[, .(donor_id, steatotic, fibrosis, reported_sex, age_in_yr)],
               by = "donor_id", all.x = TRUE, sort = FALSE)
setorder(units, unit_index)
stopifnot(identical(units$unit_index, 0:(ncol(cnt) - 1)))
units <- units[nuclei >= MIN_NUCLEI]
cnt <- cnt[, units$unit_index + 1L]
lineages <- sort(unique(units$lineage_id))
cat("units kept:", nrow(units), " donors:", uniqueN(units$donor_id), "\n")

# ---- liftover and GC -------------------------------------------------------
pk <- data.table(peak_id = peaks)
pk[, c("chr", "s", "e") := tstrsplit(peak_id, "-", fixed = TRUE)]
pk[, `:=`(s = as.integer(s), e = as.integer(e))]
pk[, width38 := e - s + 1L]
auto <- paste0("chr", 1:22)
pk[, autosome := chr %in% auto]
gr38 <- GRanges(pk[autosome == TRUE, chr], IRanges(pk[autosome == TRUE, s], pk[autosome == TRUE, e]))
names(gr38) <- pk[autosome == TRUE, peak_id]
seq38 <- getSeq(BSgenome.Hsapiens.UCSC.hg38, gr38)
gc <- as.numeric(letterFrequency(seq38, "GC", as.prob = TRUE))
pk[autosome == TRUE, gc := gc]
ch <- import.chain(CHAIN)
lo <- unlist(liftOver(gr38, ch))
lo <- data.table(peak_id = names(lo), chr19 = as.character(seqnames(lo)),
                 s19 = start(lo), e19 = end(lo))
lo <- lo[, .(n_chr = uniqueN(chr19), chr19 = chr19[1], s19 = min(s19), e19 = max(e19)), by = peak_id]
pk <- merge(pk, lo, by = "peak_id", all.x = TRUE, sort = FALSE)
pk[, lifted := autosome & !is.na(chr19) & n_chr == 1L & chr19 == chr &
                (e19 - s19 + 1L) <= 2L * width38]
rep <- data.table(total = nrow(pk), autosomal = sum(pk$autosome),
                  mapped_any = sum(!is.na(pk$chr19)), lifted_kept = sum(pk$lifted))
fwrite(rep, file.path(OUT, "peaks/liftover_report.tsv"), sep = "\t")
print(rep)
pk <- pk[match(peaks, peak_id)]  # original order

# ---- helpers ----------------------------------------------------------------
top_sets <- function(t, ids) {
  n_top <- ceiling(TOP_FRAC * length(t))
  i_up <- order(-t)[seq_len(n_top)]; i_up <- i_up[t[i_up] > 0]
  i_dn <- order(t)[seq_len(n_top)];  i_dn <- i_dn[t[i_dn] < 0]
  list(gain = ids[i_up], lost = ids[i_dn])
}
fit_contrast <- function(lin, var, acc) {
  u <- units[lineage_id == lin & !is.na(get(var))]
  grp <- u[[var]]
  if (sum(grp) < MIN_PER_GROUP || sum(!grp) < MIN_PER_GROUP) return(NULL)
  y <- DGEList(cnt[acc, paste(u$donor_id, u$lineage_id, sep = "|"), drop = FALSE])
  y <- calcNormFactors(y, method = "TMM")
  df <- data.frame(g = factor(grp, levels = c(FALSE, TRUE)), sex = factor(u$reported_sex),
                   age = u$age_in_yr, well = droplevels(factor(u$well_id)))
  design <- model.matrix(~ g + sex + age + well, data = df)
  stopifnot(qr(design)$rank == ncol(design))
  v <- voomWithQualityWeights(y, design)
  fit <- eBayes(lmFit(v, design))
  tt <- topTable(fit, coef = "gTRUE", number = Inf, sort.by = "none")
  data.table(peak_id = rownames(tt), logFC = tt$logFC, AveExpr = tt$AveExpr, t = tt$t,
             P.Value = tt$P.Value, adj.P.Val = tt$adj.P.Val,
             n_true = sum(grp), n_false = sum(!grp), df_resid = fit$df.residual[1])
}

ann <- list(); meta <- list()
add_ann <- function(name, type, lin, ids, contrast = NA, direction = NA, draw = NA) {
  ann[[name]] <<- ids
  meta[[name]] <<- data.table(annotation = name, type = type, lineage = lin, contrast = contrast,
                              direction = direction, draw = draw, n_peaks = length(ids),
                              n_lifted = sum(pk$lifted[match(ids, pk$peak_id)]))
}
add_ann("UNION", "union", NA, pk$peak_id)

# ---- accessible universe per lineage (Amendment 1: H24 peak rule) ------------
# A peak is accessible in lineage L if raw count >= 5 in >= 50% of L's donor units.
acc_list <- list()
for (L in lineages) {
  cl <- cnt[, units$lineage_id == L, drop = FALSE]
  acc_list[[L]] <- rownames(cl)[rowSums(cl >= 5) >= 0.5 * ncol(cl)]
  cat(L, "accessible peaks:", length(acc_list[[L]]), "\n")
}
spec_universe <- unique(unlist(acc_list))

# ---- specificity --------------------------------------------------------------
lin_f <- factor(units$lineage_id, levels = lineages)
y <- DGEList(cnt[spec_universe, , drop = FALSE])
y <- calcNormFactors(y, method = "TMM")
design <- model.matrix(~ 0 + lin_f + factor(units$donor_id))
colnames(design) <- make.names(colnames(design))
v <- voomWithQualityWeights(y, design)
fit <- lmFit(v, design)
lcols <- paste0("lin_f", lineages)
for (L in lineages) {
  cvec <- setNames(rep(0, ncol(design)), colnames(design))
  cvec[paste0("lin_f", L)] <- 1
  cvec[setdiff(lcols, paste0("lin_f", L))] <- -1 / (length(lineages) - 1)
  f2 <- eBayes(contrasts.fit(fit, matrix(cvec, ncol = 1, dimnames = list(names(cvec), L))))
  tt <- topTable(f2, coef = 1, number = Inf, sort.by = "none")
  fwrite(data.table(peak_id = rownames(tt), logFC = tt$logFC, AveExpr = tt$AveExpr, t = tt$t,
                    P.Value = tt$P.Value),
         file.path(OUT, sprintf("da/specificity_%s.tsv.gz", L)), sep = "\t")
  s <- top_sets(tt$t, rownames(tt))
  add_ann(paste0("SPEC_", L), "spec", L, s$gain)
}
rm(y, v, fit); gc()

# ---- per-lineage accessible universe and DA ------------------------------------
null_inputs <- list()
for (L in lineages) {
  acc <- acc_list[[L]]
  add_ann(paste0("ACC_", L), "acc", L, acc)
  for (ct in c("steatotic", "fibrosis")) {
    tab <- fit_contrast(L, ct, acc)
    if (is.null(tab)) { cat(L, ct, "not DA-testable\n"); next }
    fwrite(tab, file.path(OUT, sprintf("da/%s_%s.tsv.gz", ct, L)), sep = "\t")
    tag <- if (ct == "steatotic") "STEA" else "FIB"
    s <- top_sets(tab$t, tab$peak_id)
    add_ann(sprintf("%s_GAIN_%s", tag, L), "gain", L, s$gain, ct, "gain")
    add_ann(sprintf("%s_LOST_%s", tag, L), "lost", L, s$lost, ct, "lost")
    if (ct == "steatotic") {
      fg <- tab[adj.P.Val < 0.05 & t > 0, peak_id]; fl <- tab[adj.P.Val < 0.05 & t < 0, peak_id]
      add_ann(sprintf("FDR_GAIN_%s", L), "fdr_gain", L, fg, ct, "gain")
      add_ann(sprintf("FDR_LOST_%s", L), "fdr_lost", L, fl, ct, "lost")
      null_inputs[[L]] <- list(tab = tab, gain = s$gain, lost = s$lost)
    }
  }
}

# ---- seeded permutation null ------------------------------------------------
set.seed(SEED, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
qbin <- function(x) as.integer(ceiling(5 * frank(x, ties.method = "first") / length(x)))
strata_out <- list()
for (L in sort(names(null_inputs))) {
  ni <- null_inputs[[L]]
  d <- merge(ni$tab[, .(peak_id, AveExpr)], pk[lifted == TRUE, .(peak_id, width38, gc)], by = "peak_id")
  setorder(d, peak_id)
  d[, stratum := 100L * qbin(width38) + 10L * qbin(gc) + qbin(AveExpr)]
  d[, label := fifelse(peak_id %chin% ni$gain, "gain", fifelse(peak_id %chin% ni$lost, "lost", "other"))]
  strata_out[[L]] <- d[, .(n = .N, n_gain = sum(label == "gain"), n_lost = sum(label == "lost")),
                       by = stratum][, lineage := L]
  setorder(d, stratum, peak_id)
  idx_by_s <- split(seq_len(nrow(d)), d$stratum)
  for (k in seq_len(K_NULL)) {
    lab <- d$label
    for (sname in names(idx_by_s)) {
      ii <- idx_by_s[[sname]]
      lab[ii] <- lab[ii][sample.int(length(ii))]
    }
    add_ann(sprintf("NULL_GAIN_%s_%d", L, k), "null_gain", L, d$peak_id[lab == "gain"], "steatotic", "gain", k)
    add_ann(sprintf("NULL_LOST_%s_%d", L, k), "null_lost", L, d$peak_id[lab == "lost"], "steatotic", "lost", k)
  }
}
fwrite(rbindlist(strata_out), file.path(OUT, "peaks/null_strata.tsv"), sep = "\t")

# ---- write lifted peaks and membership ----------------------------------------
meta <- rbindlist(meta)
meta[, testable := n_lifted >= MIN_PEAKS]
fwrite(meta, file.path(OUT, "peaks/annotations.tsv"), sep = "\t")
pl <- pk[lifted == TRUE]
pl[, chrn := as.integer(sub("chr", "", chr19))]
setorder(pl, chrn, s19, e19)
fwrite(pl[, .(peak_id, chr = chrn, start = s19, end = e19, width38, gc)],
       file.path(OUT, "peaks/peaks_hg19.tsv"), sep = "\t")
con <- file(file.path(OUT, "peaks/membership.u8"), "wb")
for (a in meta$annotation) {
  writeBin(as.raw(as.integer(pl$peak_id %chin% ann[[a]])), con)
}
close(con)
cat("annotations:", nrow(meta), " lifted peaks:", nrow(pl), "\n")
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo_01.txt"))
cat("H22 step 1 done\n")
