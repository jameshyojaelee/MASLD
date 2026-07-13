#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table); library(metafor)})

stop("RETIRED NONCANONICAL: bespoke per-intron Welch-PSI/metafor analysis is not assay-native. Use Scripts 108-110 joint LeafCutter 2.0.3; cohort runs are replication only.")

root <- Sys.getenv("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
prod <- file.path(root, "GWAS/finemapping/results/seqfunc/disease_splicing/production")
out <- file.path(prod, "meta_analysis"); dir.create(out, recursive=TRUE, showWarnings=FALSE)
cohorts <- c("GSE126848","GSE130970","GSE135251","GSE162694","GSE213621")

parse_counts <- function(dataset) {
  f <- file.path(prod,dataset,paste0(dataset,"_perind.counts.gz"))
  sig <- file.path(prod,dataset,paste0(dataset,"_ds_cluster_significance.txt"))
  stopifnot(file.exists(f), file.exists(sig), file.info(f)$size > 0, file.info(sig)$size > 0)
  x <- fread(f); setnames(x, 1L, "intron")
  samples <- names(x)[-1L]
  grp <- fread(file.path(prod,dataset,"groups.tsv"), header=FALSE,
               col.names=c("sample","group"))
  stopifnot(!anyDuplicated(grp$sample), setequal(samples, grp$sample))
  grp <- grp[match(samples, sample)]
  vals <- lapply(x[, -1L], function(z) {
    p <- tstrsplit(z, "/", fixed=TRUE); as.numeric(p[[1]]) / pmax(as.numeric(p[[2]]), 1)
  })
  psi <- as.matrix(as.data.table(vals)); storage.mode(psi) <- "double"
  den <- as.matrix(x[, lapply(.SD, function(z) as.numeric(tstrsplit(z,"/",fixed=TRUE)[[2]])), .SDcols=-1L])
  usable <- den >= 10
  ic <- which(grp$group=="Control"); id <- which(grp$group=="Disease")
  calc <- function(ii) {
    a <- psi[ii,ic][usable[ii,ic]]; b <- psi[ii,id][usable[ii,id]]
    nc <- length(a); nd <- length(b)
    if(nc < 5L || nd < 5L) return(c(nc,nd,NA,NA,NA,NA))
    va <- var(a); vb <- var(b); se <- sqrt(va/nc + vb/nd)
    c(nc,nd,mean(a),mean(b),mean(b)-mean(a),se)
  }
  st <- t(vapply(seq_len(nrow(x)), calc, numeric(6)))
  z <- st[,5]/st[,6]
  ans <- data.table(dataset=dataset, intron=x$intron, n_control=st[,1], n_disease=st[,2],
                    psi_control=st[,3], psi_disease=st[,4], delta_psi=st[,5], se=st[,6],
                    p=2*pnorm(-abs(z)))
  ans[, cluster := sub("^[^:]+:[0-9]+:[0-9]+:", "", intron)]
  ans[, c("chrom","start","end") := tstrsplit(sub(":clu_.*$","",intron), ":", fixed=TRUE)]
  ans[, `:=`(start=as.integer(start),end=as.integer(end))]
  ans[, cohort_fdr := p.adjust(p, "BH")]
  ans
}

cohort_effects_file <- file.path(out,"cohort_intron_effects.tsv")
if (file.exists(cohort_effects_file) && file.info(cohort_effects_file)$size > 0) {
  message("Reusing completed cohort effects: ", cohort_effects_file)
  per <- fread(cohort_effects_file)
} else {
  per <- rbindlist(lapply(cohorts, parse_counts), fill=TRUE)
  fwrite(per, cohort_effects_file, sep="\t")
}
# LeafCutter cluster numbers are assigned independently within each cohort and
# therefore must never be part of the cross-cohort key. Its displayed junction
# coordinates are the flanking exon bases (STAR intron start-1/end+1).
per[, strand := sub(".*_([+-])$", "\\1", intron)]
stopifnot(all(per$strand %chin% c("+", "-")))
per[, junction_id := paste(chrom, start, end, strand, sep=":")]

meta_one <- function(d) {
  d <- d[is.finite(delta_psi) & is.finite(se) & se>0]
  if(nrow(d)<2L) return(NULL)
  fit <- tryCatch(rma.uni(yi=delta_psi, sei=se, method="REML", data=d), error=function(e) NULL)
  if(is.null(fit)) return(NULL)
  data.table(k=nrow(d), meta_delta_psi=as.numeric(fit$b),
             meta_se=fit$se, ci_low=fit$ci.lb, ci_high=fit$ci.ub, p=fit$pval,
             tau2=fit$tau2, I2=fit$I2,
             n_direction_positive=sum(d$delta_psi>0), n_direction_negative=sum(d$delta_psi<0),
             datasets=paste(d$dataset,collapse=";"))
}
meta <- per[, meta_one(.SD), by=junction_id]
meta[, meta_fdr := p.adjust(p,"BH")]

# Exact splice-boundary annotation from GENCODE v49 transcript exon adjacency.
gtf <- "/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"
g <- fread(cmd=paste("zgrep -v '^#'",shQuote(gtf)), sep="\t", header=FALSE, quote="", fill=TRUE,
           select=c(1,3,4,5,7,9), col.names=c("chrom","feature","start","end","strand","attr"))
ex <- g[feature=="exon"]
ex[, transcript_id := sub('.*transcript_id "([^"]+)".*','\\1',attr)]
ex[, gene_id := sub('.*gene_id "([^"]+)".*','\\1',attr)]
ex[, gene_name := fifelse(grepl('gene_name "',attr),sub('.*gene_name "([^"]+)".*','\\1',attr),gene_id)]
setorder(ex, transcript_id, start)
# LeafCutter IDs encode the two splice-site exon boundaries, not the interior
# intron interval: upstream exon end and downstream exon start.
ann <- ex[, .(chrom=chrom[1], strand=strand[1], start=end[-.N], end=start[-1L],
              gene_id=gene_id[1], gene_name=gene_name[1]), by=transcript_id]
ann <- ann[start<=end, .(gene_id=paste(unique(gene_id),collapse=";"),
                         gene_name=paste(unique(gene_name),collapse=";"),
                         transcript_id=paste(unique(transcript_id),collapse=";")),
           by=.(chrom,start,end,strand)]
meta[, c("chrom","leaf_start","leaf_end","strand") := tstrsplit(junction_id, ":", fixed=TRUE)]
meta[, `:=`(leaf_start=as.integer(leaf_start), leaf_end=as.integer(leaf_end))]
meta[, `:=`(start=leaf_start, end=leaf_end)]
meta <- merge(meta, ann, by=c("chrom","start","end","strand"), all.x=TRUE)
meta[, annotated_exact := !is.na(gene_id)]
if (nrow(meta) < 10000L)
  stop("Cross-cohort key audit failed: fewer than 10,000 shared testable introns")
if (mean(meta$annotated_exact) < 0.20)
  stop("Coordinate/annotation audit failed: fewer than 20% exact GENCODE v49 introns")
meta <- meta[order(meta_fdr, -abs(meta_delta_psi))]
fwrite(meta, file.path(out,"random_effects_intron_meta.tsv"), sep="\t")

qc <- per[, .(n_introns=.N, n_testable=sum(is.finite(se)&se>0), n_nominal=sum(p<.05,na.rm=TRUE),
              n_fdr=sum(cohort_fdr<.05,na.rm=TRUE), median_abs_delta=median(abs(delta_psi),na.rm=TRUE)), by=dataset]
qc <- rbind(qc, data.table(dataset="META",n_introns=nrow(meta),n_testable=nrow(meta),
                           n_nominal=sum(meta$p<.05),n_fdr=sum(meta$meta_fdr<.05),
                           median_abs_delta=median(abs(meta$meta_delta_psi))))
fwrite(qc,file.path(out,"qc_summary.tsv"),sep="\t")
writeLines(c("# Disease splicing random-effects meta-analysis","",
  "Cohorts are analyzed independently; no raw junction counts are pooled across studies.",
  "Per-intron PSI uses LeafCutter numerator/cluster-denominator counts and requires denominator >=10 in >=5 samples per group.",
  "Disease-minus-control PSI and Welch standard errors are combined by REML random-effects meta-analysis for exact coordinate-matched introns present in >=2 cohorts.",
  "GENCODE v49 annotation requires exact upstream-exon-end/downstream-exon-start splice-boundary and strand agreement. Uniform-prior seqfunc firewall remains unchanged."), file.path(out,"README.md"))
message("Wrote ", nrow(per), " cohort effects and ",nrow(meta)," meta-analyzed introns")
