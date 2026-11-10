#!/usr/bin/env Rscript
# A3 step 2: liftover, hg38 allele orientation, lineage peak overlap, GC, TSS distance,
# scoring sequences; JASPAR2024 motif export (same call as ATAC v3 06_run_chromvar.R);
# motifmatchr reference calls on a validation sample.

suppressPackageStartupMessages({
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(data.table)
  library(GenomicRanges)
  library(rtracklayer)
  library(JASPAR2024)
  library(TFBSTools)
  library(motifmatchr)
  library(SummarizedExperiment)
  library(Matrix)
  library(digest)
})
set.seed(42)

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
A3 <- Sys.getenv("A3_OUT", file.path(BASE, "Analysis/Hypotheses/ATACGEN-20261024T123423Z/A3_motif_disruption"))
ATAC <- file.path(BASE, "Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-candidate-2026-08-11-r1")
CHAIN <- file.path(BASE, "data/broadaway_eqtl/hg19ToHg38.over.chain")
GENES <- file.path(BASE, "RNA-seq/results/noncoding_resource/candidates/noncoding-bulk-2026-08-11/annotation/GENCODE_v49/run_19749483/gene_identity.tsv")
LINEAGES <- c("hepatocyte", "stellate", "macrophage", "cholangiocyte", "t_nk")
WORK <- file.path(A3, "work")
RES <- file.path(A3, "results")

log_msg <- function(...) message(format(Sys.time(), "%H:%M:%S "), ...)

# ---- motifs: identical retrieval to the chromVAR run ----
jaspar_sqlite <- db(JASPAR2024())
motifs <- getMatrixSet(jaspar_sqlite,
  opts = list(collection = "CORE", tax_group = "vertebrates", all_versions = FALSE))
motif_ids <- names(motifs)
cv <- fread(cmd = paste("zcat", file.path(ATAC, "chromvar/chromvar_results.tsv.gz")))
stopifnot(setequal(motif_ids, unique(cv$motif_id)), length(motif_ids) == 879L)
lens <- vapply(motifs, function(m) ncol(TFBSTools::Matrix(m)), integer(1))
LMAX <- max(lens)
log_msg("motifs: ", length(motifs), "; max length ", LMAX)

meta <- data.table(
  motif_id = motif_ids,
  name = vapply(motifs, name, character(1)),
  length = lens,
  family = vapply(motifs, function(m) paste(unlist(tags(m)$family), collapse = ";"), character(1)),
  class = vapply(motifs, function(m) paste(unlist(m@matrixClass), collapse = ";"), character(1))
)
fwrite(meta, file.path(WORK, "motif_meta.tsv"), sep = "\t")
pfm_long <- rbindlist(lapply(motif_ids, function(id) {
  x <- TFBSTools::Matrix(motifs[[id]])
  data.table(motif_id = id, pos = seq_len(ncol(x)), A = x["A", ], C = x["C", ], G = x["G", ], T = x["T", ])
}))
fwrite(pfm_long, file.path(WORK, "motif_pfm.tsv"), sep = "\t")

# ---- background: union of the five lineage consensus peak sets ----
peaks <- lapply(LINEAGES, function(l) {
  b <- fread(file.path(ATAC, "peaks/consensus", paste0(l, ".bed")), header = FALSE)
  GRanges(b$V1, IRanges(b$V2 + 1L, b$V3))
})
names(peaks) <- LINEAGES
union_peaks <- reduce(do.call(c, unname(peaks)))
af <- colSums(alphabetFrequency(getSeq(BSgenome.Hsapiens.UCSC.hg38, union_peaks))[, c("A", "C", "G", "T")])
af <- af / sum(af)
bg <- c(A = (af[["A"]] + af[["T"]]) / 2, C = (af[["C"]] + af[["G"]]) / 2,
        G = (af[["C"]] + af[["G"]]) / 2, T = (af[["A"]] + af[["T"]]) / 2)
fwrite(data.table(base = names(bg), freq = as.numeric(bg)), file.path(WORK, "bg.tsv"), sep = "\t")
log_msg("bg: ", paste(names(bg), round(bg, 4), collapse = " "))

pwms <- do.call(PWMatrixList, lapply(motifs, toPWM))
pwm_long <- rbindlist(lapply(motif_ids, function(id) {
  x <- motifmatchr:::convert_pwm(pwms[[id]], as.numeric(bg))
  data.table(motif_id = id, pos = seq_len(ncol(x)), A = x[1, ], C = x[2, ], G = x[3, ], T = x[4, ])
}))
fwrite(pwm_long, file.path(WORK, "motif_pwm_log2.tsv"), sep = "\t")

# ---- variants ----
u <- fread(cmd = paste("zcat", shQuote(file.path(WORK, "universe_hg19.tsv.gz"))))
u <- u[vclass %in% c("high", "member_low", "null")]
log_msg("universe rows kept (high/member_low/null): ", nrow(u))
gr19 <- GRanges(paste0("chr", u$chr), IRanges(u$pos19, width = 1L))
chain <- import.chain(CHAIN)
lo <- liftOver(gr19, chain)
nmap <- lengths(lo)
one <- nmap == 1L
lo1 <- unlist(lo[one])
ok <- width(lo1) == 1L & as.character(seqnames(lo1)) == paste0("chr", u$chr[one])
idx <- which(one)[ok]
u <- u[idx]
u[, pos38 := start(lo1)[ok]]
log_msg("1:1 liftover kept: ", nrow(u), " of ", length(gr19))
lift_qc <- data.table(n_input = length(gr19), n_unmapped = sum(nmap == 0L), n_multi = sum(nmap > 1L),
                      n_one_to_one_same_chr = nrow(u))

u <- u[!(chr == 6L & pos38 >= 25e6 & pos38 <= 34e6)]
gr <- GRanges(paste0("chr", u$chr), IRanges(u$pos38, width = 1L))

for (l in LINEAGES) {
  hits <- findOverlaps(gr, peaks[[l]], select = "first")
  u[, (paste0("peak_", l)) := hits]
}
in_any <- Reduce(`|`, lapply(LINEAGES, function(l) !is.na(u[[paste0("peak_", l)]])))
keep <- in_any | u$vclass %in% c("high", "member_low")
u <- u[keep]
gr <- gr[keep]
log_msg("in any lineage peak or gated member: ", nrow(u))

ref <- as.character(getSeq(BSgenome.Hsapiens.UCSC.hg38, gr))
comp <- c(A = "T", C = "G", G = "C", T = "A")
pal <- (u$a_lo == "A" & u$a_hi == "T") | (u$a_lo == "C" & u$a_hi == "G")
alt <- rep(NA_character_, nrow(u))
flip <- rep(FALSE, nrow(u))
m1 <- u$a_lo == ref
m2 <- u$a_hi == ref
alt[m1] <- u$a_hi[m1]
alt[m2] <- u$a_lo[m2]
cm <- !(m1 | m2) & !pal
c_lo <- comp[u$a_lo]
c_hi <- comp[u$a_hi]
f1 <- cm & c_lo == ref
f2 <- cm & c_hi == ref
alt[f1] <- c_hi[f1]
alt[f2] <- c_lo[f2]
flip[f1 | f2] <- TRUE
u[, `:=`(ref = ref, alt = alt, strand_flipped = flip, palindromic = pal)]
orient_qc <- data.table(n = nrow(u), ref_match = sum(m1 | m2), strand_flip = sum(flip),
                        no_match_dropped = sum(is.na(alt)), palindromic = sum(pal))
u <- u[!is.na(alt)]
gr <- GRanges(paste0("chr", u$chr), IRanges(u$pos38, width = 1L))
log_msg("oriented: ", nrow(u))

# GC of 101 bp
gcseq <- getSeq(BSgenome.Hsapiens.UCSC.hg38, resize(gr, 101L, fix = "center"))
fr <- letterFrequency(gcseq, c("G", "C", "A", "T"))
u[, gc101 := (fr[, "G"] + fr[, "C"]) / pmax(rowSums(fr), 1)]

# TSS distance (GENCODE v49 protein_coding + lncRNA, canonical chromosomes)
g <- fread(GENES)
g <- g[toupper(as.character(is_canonical_chromosome)) == "TRUE" & gene_type %in% c("protein_coding", "lncRNA")]
tss <- GRanges(g$chromosome, IRanges(ifelse(g$strand == "+", g$start_1based, g$end_1based), width = 1L),
               gene_name = g$gene_name, gene_id = g$gene_id_base, gene_type = g$gene_type)
d <- distanceToNearest(gr, tss, ignore.strand = TRUE)
u[, tss_dist := NA_integer_]
u[queryHits(d), tss_dist := mcols(d)$distance]
u[queryHits(d), nearest_tss_gene := mcols(tss)$gene_name[subjectHits(d)]]

# peak signal: max native pooled_qvalue over retained native peaks (either cohort)
# overlapping the consensus peak
SUPPORT <- file.path(ATAC, "peaks/native_peak_support.tsv")
support <- fread(SUPPORT, select = c("cohort", "lineage", "chrom", "start0", "end", "pooled_qvalue", "retained"))
support <- support[toupper(as.character(retained)) == "TRUE"]
for (l in LINEAGES) {
  pcol <- paste0("peak_", l)
  sel <- which(!is.na(u[[pcol]]))
  cons <- peaks[[l]][unique(u[[pcol]][sel])]
  nat <- support[lineage == l]
  ngr <- GRanges(nat$chrom, IRanges(nat$start0 + 1L, nat$end))
  ov <- findOverlaps(cons, ngr)
  sig <- tapply(nat$pooled_qvalue[subjectHits(ov)], queryHits(ov), max)
  sigv <- rep(NA_real_, length(cons))
  sigv[as.integer(names(sig))] <- as.numeric(sig)
  map <- match(u[[pcol]][sel], unique(u[[pcol]][sel]))
  u[sel, (paste0("signal_", l)) := sigv[map]]
  u[sel, (paste0("peakid_", l)) := paste0(seqnames(peaks[[l]])[u[[pcol]][sel]], ":",
                                          start(peaks[[l]])[u[[pcol]][sel]] - 1L, "-",
                                          end(peaks[[l]])[u[[pcol]][sel]])]
}

# scoring sequence: ref allele, +-(LMAX-1)
in_any <- Reduce(`|`, lapply(LINEAGES, function(l) !is.na(u[[paste0("peak_", l)]])))
u[, in_any_peak := in_any]
sgr <- resize(gr[in_any], 2L * LMAX - 1L, fix = "center")
u[in_any_peak == TRUE, seq := as.character(getSeq(BSgenome.Hsapiens.UCSC.hg38, sgr))]
stopifnot(all(substr(u$seq[in_any], LMAX, LMAX) == u$ref[in_any]))

for (l in LINEAGES) u[, (paste0("peak_", l)) := !is.na(get(paste0("peak_", l)))]
fwrite(u, file.path(WORK, "variants_hg38.tsv.gz"), sep = "\t")
fwrite(cbind(lift_qc, orient_qc), file.path(RES, "liftover_orientation_qc.tsv"), sep = "\t")
cnt <- rbindlist(lapply(LINEAGES, function(l) {
  x <- u[get(paste0("peak_", l)) == TRUE]
  data.table(lineage = l, high = sum(x$vclass == "high"), member_low = sum(x$vclass == "member_low"),
             null = sum(x$vclass == "null"))
}))
fwrite(cnt, file.path(RES, "in_peak_class_counts.tsv"), sep = "\t")
print(cnt)

# validation sample for the Python scorer: motifmatchr calls on reference sequences
vi <- sample(which(u$in_any_peak), 2000L)
vseq <- DNAStringSet(u$seq[vi])
mm <- matchMotifs(pwms, vseq, bg = as.numeric(bg), out = "scores")
sc <- as(assays(mm)$motifScores, "TsparseMatrix")
mt <- as(assays(mm)$motifMatches, "TsparseMatrix")
fwrite(data.table(row = seq_along(vi), vkey = u$key[vi], seq = u$seq[vi]), file.path(WORK, "validation_seqs.tsv"), sep = "\t")
fwrite(data.table(row = mt@i + 1L, motif_id = colnames(mm)[mt@j + 1L]), file.path(WORK, "validation_motifmatchr_matches.tsv"), sep = "\t")
fwrite(data.table(row = sc@i + 1L, motif_id = colnames(mm)[sc@j + 1L], score = sc@x), file.path(WORK, "validation_motifmatchr_scores.tsv"), sep = "\t")

ins <- data.table(path = c(CHAIN, GENES, file.path(ATAC, "chromvar/chromvar_results.tsv.gz"), jaspar_sqlite, SUPPORT,
                           file.path(ATAC, "peaks/consensus", paste0(LINEAGES, ".bed"))))
ins[, sha256 := vapply(path, function(p) digest(p, algo = "sha256", file = TRUE, serialize = FALSE), character(1))]
fwrite(ins, file.path(RES, "input_sha256.tsv"), sep = "\t", append = TRUE, col.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(RES, "sessionInfo_chromvar_env.txt"))
log_msg("done")
