#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(rtracklayer)})
run <- 'Analysis/Hypotheses/COLOC-IDENTITY-20261107T150017Z'
out <- file.path(run,'gram_inputs');stopifnot(!dir.exists(out));dir.create(out)
d <- fread(file.path(run,'mprip_refit_inputs/eqtl_aligned.tsv'))
pos <- sort(unique(d$position))
gr <- GRanges('chr17',IRanges(pos,width=1),strand='+')
forward <- liftOver(gr,import.chain('data/broadaway_eqtl/hg19ToHg38.over.chain'))
reverse_chain <- import.chain('data/broadaway_eqtl/hg38ToHg19.over.chain')
rows <- lapply(seq_along(pos),function(i) {
  f <- forward[[i]]
  ok <- length(f)==1L && as.character(seqnames(f))=='chr17' && as.character(strand(f))=='+'
  back_ok <- FALSE
  if(ok) {
    b <- liftOver(f,reverse_chain)[[1]]
    back_ok <- length(b)==1L && as.character(seqnames(b))=='chr17' && start(b)==pos[i] && as.character(strand(b))=='+'
  }
  data.table(position=pos[i],position_hg38=if(ok) start(f) else NA_integer_,reciprocal_unique=ok&&back_ok)
})
fwrite(merge(d,rbindlist(rows),by='position',sort=FALSE),file.path(out,'mapped_reference_events.tsv'),sep='\t')
capture.output(sessionInfo(),file=file.path(out,'liftover_sessionInfo.txt'))
