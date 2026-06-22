#!/usr/bin/env Rscript
# Exon + splice-junction windows that make each target's DOMINANT isoform selective vs the UNION of
# all OTHER EXPRESSED isoforms of the same gene (the comparison that actually matters for a Cas13
# guide that must hit the dominant disease isoform while sparing everything else present in mouse
# liver). Generalizes the earlier CERKL/RCAN1-only version to ALL targets. Expressed isoform set is
# read from 12b output (actionable_isoform_detail.tsv); exon coords from GENCODE vM38.
suppressPackageStartupMessages({ library(data.table) })
GTF <- commandArgs(trailingOnly=TRUE)[1]
PROJ <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
MRES <- file.path(PROJ, "RNA-seq/results/isoform_diversity/mouse"); strip <- function(x) sub("\\..*$","",x)

det <- fread(file.path(MRES, "actionable_isoform_detail.tsv"))     # per-transcript expressed flag
tg  <- fread(file.path(MRES, "gse213621_targets_mouse_isoform_structure.tsv"))[isoform_selective_possible==TRUE]
det[, txid := strip(transcript)]
det <- merge(det, tg[, .(human=symbol, dom=strip(mouse_dom_iso))], by="human")
want <- unique(det[expressed==TRUE, txid])

g <- fread(cmd=paste("zcat", GTF, "| awk -F'\\t' '$3==\"exon\"'"), header=FALSE, sep="\t", quote="")
setnames(g, c("chr","src","feat","start","end","score","strand","frame","attr"))
g[, txid := strip(sub('.*transcript_id "([^"]+)".*','\\1', attr))]
ex <- g[txid %in% want, .(chr, start, end, strand, txid)]

tx_exons <- function(t){ e<-ex[txid==t]; if(!nrow(e)) return(NULL)
  e <- if(e$strand[1]=="+") e[order(start)] else e[order(-start)]
  e[, exon:=.I]; e[, len:=end-start+1]; e[, cdna_end:=cumsum(len)]; e[, cdna_start:=cdna_end-len+1]; e[] }
tx_junc <- function(e){ if(is.null(e)||nrow(e)<2) return(data.table()); plus<-e$strand[1]=="+"
  rbindlist(lapply(seq_len(nrow(e)-1), function(i){ d<-if(plus)e$end[i] else e$start[i]; a<-if(plus)e$start[i+1] else e$end[i+1]
    data.table(up_exon=e$exon[i],dn_exon=e$exon[i+1],donor_g=d,acceptor_g=a,cdna_pos=e$cdna_end[i],jkey=paste0(d,"_",a))})) }
exkey <- function(e) paste0(e$start,"_",e$end)

allwin <- list()
for (hd in unique(det$human)) {
  d <- det[human==hd]; Dt <- d$dom[1]
  others <- setdiff(d[expressed==TRUE, txid], Dt)
  eD <- tx_exons(Dt); jD <- tx_junc(eD); if(is.null(eD)) next
  oth_ex <- unique(unlist(lapply(others, function(t) exkey(tx_exons(t)))))
  oth_j  <- unique(unlist(lapply(others, function(t){ej<-tx_junc(tx_exons(t)); if(nrow(ej))ej$jkey else character()})))
  uE <- eD[!(exkey(eD) %in% oth_ex)]
  uJ <- if(nrow(jD)) jD[!(jkey %in% oth_j)] else jD
  cat(sprintf("%-9s %s: %d other expressed isoforms | unique-to-dominant: %d exons, %d junctions\n",
              hd, d$mouse[1], length(others), nrow(uE), nrow(uJ)))
  if(nrow(uE)) allwin[[paste0(hd,"E")]] <- cbind(gene=hd,mouse=d$mouse[1],transcript=Dt,
    comparison="vs_expressed_set", window_type="exon",
    uE[,.(chr,strand,cdna_start,cdna_end,g_start=start,g_end=end,len)])
  if(nrow(uJ)) allwin[[paste0(hd,"J")]] <- cbind(gene=hd,mouse=d$mouse[1],transcript=Dt,
    comparison="vs_expressed_set", window_type="junction",
    uJ[,.(chr=eD$chr[1],strand=eD$strand[1],cdna_start=cdna_pos,cdna_end=cdna_pos,
          g_start=donor_g,g_end=acceptor_g,len=NA_integer_)])
}
if(length(allwin)) fwrite(rbindlist(allwin,fill=TRUE), file.path(MRES,"actionable_expressedset_windows.tsv"), sep="\t")
cat("\n-> actionable_expressedset_windows.tsv\n")
