#!/usr/bin/env Rscript
# Exact GRCh38 native peak scan. Motif compatibility does not imply TF binding.
# Args: native peaks.tsv, genes.tsv, comma-separated Ensembl targets, NEW output.
a <- commandArgs(TRUE)
stopifnot(length(a)==4L,nzchar(Sys.getenv("SLURM_JOB_ID")),!dir.exists(a[4]))
suppressPackageStartupMessages({
  library(TFBSTools); library(JASPAR2024); library(motifmatchr)
  library(BSgenome.Hsapiens.UCSC.hg38); library(GenomicRanges)
  library(Matrix); library(digest)
})
peaks <- read.delim(a[1],stringsAsFactors=FALSE)
genes <- read.delim(a[2],stringsAsFactors=FALSE)
targets <- strsplit(a[3],",",fixed=TRUE)[[1]]
stopifnot(length(targets)<=5L,all(targets %in% genes$gene_id),
          !anyDuplicated(peaks$peak_id),!anyDuplicated(genes$gene_id))
sel <- integer()
for (id in targets) {
  g <- genes[genes$gene_id==id,]
  cis <- which(peaks$chrom==g$chrom & peaks$start0<g$tss0+250001 &
    peaks$end>max(0,g$tss0-250000) &
    !(peaks$start0<g$tss0+2001 & peaks$end>max(0,g$tss0-2000)))
  distance <- abs((peaks$start0[cis]+peaks$end[cis])/2-g$tss0)
  sel <- union(sel,head(cis[order(distance,peaks$peak_id[cis])],25L))
}
stopifnot(length(sel)>0L,length(sel)<=125L)
p <- peaks[sel,]
gr <- GRanges(p$chrom,IRanges(p$start0+1,p$end))
database <- db(JASPAR2024())
motifs <- getMatrixSet(database,opts=list(collection="CORE",species=9606,all_versions=FALSE))
stopifnot(length(motifs)>0L)
started <- proc.time()[3]
matches <- motifMatches(matchMotifs(motifs,gr,genome=BSgenome.Hsapiens.UCSC.hg38,
  out="matches",p.cutoff=5e-5,bg="genome"))
elapsed <- proc.time()[3]-started
ij <- summary(matches)
mapped <- list(); annotations <- list()
for (j in seq_along(motifs)) {
  name <- TFBSTools::name(motifs[[j]])
  components <- strsplit(name,"::",fixed=TRUE)[[1]]
  ids <- list()
  for(symbol in components) {
    m <- unique(genes$gene_id[genes$gene_name==symbol])
    # Ambiguous gene symbols are not silently resolved.
    if(length(m)==1L) ids[[length(ids)+1L]] <- data.frame(tf_gene_id=m,
      mapping_state=if(length(components)>1L) "complex_component" else "exact_symbol")
  }
  annotation <- TFBSTools::tags(motifs[[j]])
  family <- paste(annotation$family,collapse=";")
  class <- paste(annotation$class,collapse=";")
  annotations[[j]] <- data.frame(motif_id=names(motifs)[j],motif_name=name,
    motif_family=family,motif_class=class,
    n_components=length(components),n_mapped_components=length(ids),
    binding_inferred=FALSE,unique_TF_inferred=FALSE)
  hits <- ij$i[ij$j==j]
  if(length(hits) && length(ids)) {
    tf <- do.call(rbind,ids)
    z <- merge(data.frame(peak_id=p$peak_id[hits]),tf,by=NULL)
    z$motif_id <- names(motifs)[j]; z$motif_name <- name
    z$motif_family <- family; z$motif_class <- class
    z$unique_TF_inferred <- FALSE; z$binding_inferred <- FALSE
    mapped[[length(mapped)+1L]] <- z
  }
}
dir.create(a[4],recursive=TRUE)
writeMM(matches,file.path(a[4],"matches.mtx"))
write.table(p,file.path(a[4],"peaks.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write.table(do.call(rbind,annotations),file.path(a[4],"motifs.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
if(length(mapped)) out <- do.call(rbind,mapped) else out <- data.frame(
  peak_id=character(),tf_gene_id=character(),mapping_state=character(),motif_id=character(),
  motif_name=character(),motif_family=character(),motif_class=character(),
  unique_TF_inferred=logical(),binding_inferred=logical())
write.table(out,file.path(a[4],"motif_tf_map.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
# Hash sequences actually scanned; genome package/version also recorded below.
seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38,gr)
source <- data.frame(database_path=database,
  database_sha256=digest(database,algo="sha256",file=TRUE),
  native_peak_file_sha256=digest(a[1],algo="sha256",file=TRUE),
  scanned_sequences_sha256=digest(paste(as.character(seqs),collapse="\n"),algo="sha256",serialize=FALSE),
  genome_package="BSgenome.Hsapiens.UCSC.hg38",
  genome_package_version=as.character(packageVersion("BSgenome.Hsapiens.UCSC.hg38")),
  genome_2bit_path=system.file("extdata/single_sequences.2bit",package="BSgenome.Hsapiens.UCSC.hg38"),
  genome_2bit_sha256=digest(system.file("extdata/single_sequences.2bit",package="BSgenome.Hsapiens.UCSC.hg38"),algo="sha256",file=TRUE),
  genome="GRCh38/UCSC_hg38",coordinate_input="BED0_halfopen",
  p_cutoff=5e-5,background="genome",n_peaks=nrow(p),n_motifs=length(motifs),
  match_seconds=elapsed,motif_tf_mapping="exact_unique_annotation_symbol_or_complex_component",
  inference="motif_compatibility_only")
write.table(source,file.path(a[4],"sources.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
capture.output(sessionInfo(),file=file.path(a[4],"sessionInfo.txt"))
