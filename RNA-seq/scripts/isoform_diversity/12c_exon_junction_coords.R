#!/usr/bin/env Rscript
# Exon + splice-junction coordinates for each target's DOMINANT (D) vs CANONICAL (C) transcript,
# from GENCODE vM38. Computes isoform-UNIQUE exons + UNIQUE junctions (D-vs-C set difference) = the
# windows a Cas13 guide can tile to discriminate the dominant isoform from the canonical one.
# Targets + dominant transcript read DYNAMICALLY from the targetability table; canonical from
# tx2gene (is_ensembl_canonical). Cas13 binds mature mRNA: a "junction" = the seam where two exons
# adjacent only in this isoform meet. Library-design only (not manuscript).
suppressPackageStartupMessages({ library(data.table) })
GTF <- commandArgs(trailingOnly=TRUE)[1]
PROJ <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
MRES <- file.path(PROJ, "RNA-seq/results/isoform_diversity/mouse")
OUT  <- MRES; strip <- function(x) sub("\\..*$","",x)

tg <- fread(file.path(MRES, "gse213621_targets_mouse_isoform_structure.tsv"))[isoform_selective_possible == TRUE]
t2g <- fread(file.path(PROJ, "RNA-seq/results/isoform_diversity/_index/tx2gene_vM38_primary.tsv.gz"))
t2g[, mensg := strip(geneid)]
canon_tx <- function(mensg) { c <- t2g[mensg==strip(mensg) & is_ensembl_canonical==1, txname]; if(length(c)) strip(c[1]) else NA }
pairs <- data.table(human=tg$symbol, mouse=tg$mouse_symbol,
                    D=strip(tg$mouse_dom_iso),
                    C=vapply(tg$mouse_ensg, canon_tx, character(1)))
want <- unique(c(pairs$D, pairs$C)); want <- want[!is.na(want)]

cat("[12c] reading", GTF, "for", nrow(pairs), "targets\n")
g <- fread(cmd=paste("zcat", GTF, "| awk -F'\\t' '$3==\"exon\"'"), header=FALSE, sep="\t", quote="")
setnames(g, c("chr","src","feat","start","end","score","strand","frame","attr"))
g[, txid := strip(sub('.*transcript_id "([^"]+)".*','\\1', attr))]
ex <- g[txid %in% want, .(chr, start, end, strand, txid)]

tx_exons <- function(t){ if(is.na(t)) return(NULL); e<-ex[txid==t]; if(!nrow(e)) return(NULL)
  e <- if(e$strand[1]=="+") e[order(start)] else e[order(-start)]
  e[, exon:=.I]; e[, len:=end-start+1]; e[, cdna_end:=cumsum(len)]; e[, cdna_start:=cdna_end-len+1]; e[] }
tx_junctions <- function(e){ if(is.null(e)||nrow(e)<2) return(data.table()); plus<-e$strand[1]=="+"
  rbindlist(lapply(seq_len(nrow(e)-1), function(i){
    d<-if(plus)e$end[i] else e$start[i]; a<-if(plus)e$start[i+1] else e$end[i+1]
    data.table(up_exon=e$exon[i],dn_exon=e$exon[i+1],donor_g=d,acceptor_g=a,cdna_pos=e$cdna_end[i],jkey=paste0(d,"_",a))})) }
exkey <- function(e) paste0(e$start,"_",e$end)

allrows <- list(); summ <- list(); uniqwin <- list()
for (i in seq_len(nrow(pairs))) {
  hd<-pairs$human[i]; ms<-pairs$mouse[i]; Dt<-pairs$D[i]; Ct<-pairs$C[i]
  eD<-tx_exons(Dt); eC<-tx_exons(Ct); jD<-tx_junctions(eD); jC<-tx_junctions(eC)
  same <- !is.na(Ct) && Dt==Ct; noC <- is.null(eC)
  uD_ex <- if(!same && !noC) eD[!(exkey(eD) %in% exkey(eC))] else eD[0]
  uD_j  <- if(!same && !noC && nrow(jD)) jD[!(jkey %in% jC$jkey)] else jD[0]
  summ[[hd]] <- data.table(gene=hd, mouse=ms, chr=eD$chr[1], strand=eD$strand[1],
    D_tx=Dt, D_exons=nrow(eD), D_len=max(eD$cdna_end),
    C_tx=Ct, C_exons=if(noC)0L else nrow(eC), C_len=if(noC)NA_integer_ else max(eC$cdna_end),
    uniqD_exons=nrow(uD_ex), uniqD_junctions=nrow(uD_j), DequalsC=same, C_missing=noC)
  if (!same && !noC && nrow(uD_ex)) uniqwin[[paste0(hd,"E")]] <- cbind(
    gene=hd, mouse=ms, transcript=Dt, comparison="vs_canonical", window_type="exon",
    uD_ex[, .(chr, strand, cdna_start, cdna_end, g_start=start, g_end=end, len)])
  if (!same && !noC && nrow(uD_j)) uniqwin[[paste0(hd,"J")]] <- cbind(
    gene=hd, mouse=ms, transcript=Dt, comparison="vs_canonical", window_type="junction",
    uD_j[, .(chr=eD$chr[1], strand=eD$strand[1], cdna_start=cdna_pos, cdna_end=cdna_pos,
             g_start=donor_g, g_end=acceptor_g, len=NA_integer_)])
  for (tag in c("D","C")) { e<-get(paste0("e",tag)); if(is.null(e)) next
    allrows[[paste0(hd,tag)]] <- cbind(gene=hd,mouse=ms,role=tag,txid=get(paste0(tag,"t")),
      e[, .(exon,chr,start,end,strand,len,cdna_start,cdna_end)]) }
}
sumdt <- rbindlist(summ)
cat("\n================ DvsC SUMMARY ================\n"); print(sumdt, row.names=FALSE)
fwrite(rbindlist(allrows), file.path(OUT,"actionable_exon_coords.tsv"), sep="\t")
fwrite(sumdt, file.path(OUT,"actionable_DvsC_summary.tsv"), sep="\t")
if (length(uniqwin)) fwrite(rbindlist(uniqwin, fill=TRUE), file.path(OUT,"actionable_unique_windows_detail.tsv"), sep="\t")
# dominant transcript IDs (versioned) for the gffread step in 12e
fwrite(data.table(transcript=strip(tg$mouse_dom_iso), dom_iso=tg$mouse_dom_iso), file.path(OUT,"actionable_dominant_ids.tsv"), sep="\t")
cat("\n-> actionable_exon_coords.tsv / actionable_DvsC_summary.tsv / actionable_unique_windows_detail.tsv\n")
