#!/usr/bin/env Rscript
# Slice the isoform-selective windows (from 12c vs-canonical + 12d vs-expressed-set) out of each
# dominant transcript's spliced cDNA and emit the mature-mRNA SENSE sequence + reverse-complement
# (= RfxCas13d crRNA spacer orientation, since the spacer is antisense to target mRNA). The cDNA is
# produced by gffread from the GRCm39 genome (run separately just before this). Exon windows = full
# tileable region; junction windows = 30 nt centered on the spliced seam.
suppressPackageStartupMessages({ library(data.table); library(Biostrings) })
PROJ <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
MRES <- file.path(PROJ, "RNA-seq/results/isoform_diversity/mouse"); FLANK <- 15

cdna <- readDNAStringSet(file.path(MRES, "actionable_dominant_cdna.fa"))
names(cdna) <- sub("\\..*$", "", sub("\\s.*", "", names(cdna)))

w1f <- file.path(MRES, "actionable_unique_windows_detail.tsv")     # vs canonical (12c)
w2f <- file.path(MRES, "actionable_expressedset_windows.tsv")      # vs expressed set (12d)
w  <- rbindlist(list(if(file.exists(w1f)) fread(w1f), if(file.exists(w2f)) fread(w2f)), fill=TRUE)
cat(sprintf("[12e] %d windows over %d genes\n", nrow(w), uniqueN(w$gene)))

res <- rbindlist(lapply(seq_len(nrow(w)), function(i) {
  tx <- sub("\\..*$","", w$transcript[i]); s <- cdna[[tx]]; if (is.null(s)) return(NULL); L <- length(s)
  if (w$window_type[i] == "exon") { a <- w$cdna_start[i]; b <- min(w$cdna_end[i], L) }
  else { p <- w$cdna_start[i]; a <- max(1, p-FLANK+1); b <- min(L, p+FLANK) }
  sq <- subseq(s, a, b)
  data.table(gene=w$gene[i], mouse=w$mouse[i], transcript=tx, comparison=w$comparison[i],
             window_type=w$window_type[i], cdna_from=a, cdna_to=b, win_len=b-a+1,
             g_start=w$g_start[i], g_end=w$g_end[i],
             mrna_sense=as.character(sq), crRNA_spacer_revcomp=as.character(reverseComplement(sq)))
}))
setorder(res, gene, comparison, -win_len)
fwrite(res, file.path(MRES, "actionable_window_sequences.tsv"), sep="\t")
fa <- unlist(lapply(seq_len(nrow(res)), function(i)
  c(sprintf(">%s|%s|%s|%s|cdna%d-%d|len%d", res$gene[i], res$mouse[i], res$transcript[i],
            res$window_type[i], res$cdna_from[i], res$cdna_to[i], res$win_len[i]), res$mrna_sense[i])))
writeLines(fa, file.path(MRES, "actionable_window_sequences.fa"))
cat("\n==== windows per gene (vs_expressed_set) ====\n")
print(res[comparison=="vs_expressed_set", .N, by=.(gene, window_type)][order(gene)], row.names=FALSE)
cat("\n-> actionable_window_sequences.tsv / .fa\n")
