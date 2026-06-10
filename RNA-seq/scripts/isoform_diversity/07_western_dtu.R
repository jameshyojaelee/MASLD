#!/usr/bin/env Rscript
# Cas13 'Western' diet-group DTU (DIAMOND + WD + WD_Fructose), reusing existing
# vM38 kallisto quant (same index as mouse arm; within-group contrast is
# strand-cancelling so unstranded is valid). satuRn + stageR, same as 04_dtu.R.
suppressPackageStartupMessages({
  library(data.table); library(tximport); library(DRIMSeq); library(satuRn)
  library(SummarizedExperiment); library(stageR)
})
PROJ <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
ID   <- file.path(PROJ, "RNA-seq/scripts/isoform_diversity")
RES  <- file.path(PROJ, "RNA-seq/results/isoform_diversity/mouse")
GROUP <- Sys.getenv("CAS13GROUP", "Western")

drv  <- fread(file.path(ID, "western_driver.tsv"))
samp <- fread(file.path(ID, "western_samples.tsv"))
samp <- merge(samp, drv[, .(sample_id, abundance_path)], by = "sample_id")
samp <- samp[cas13_group == GROUP]
tx2gene <- fread(file.path(PROJ, "RNA-seq/results/isoform_diversity/_index/tx2gene_vM38_primary.tsv.gz"))
cat(sprintf("[western] %s: %d samples (%dC/%dD), %d datasets\n", GROUP, nrow(samp),
            sum(samp$group_binary=="Control"), sum(samp$group_binary=="Disease"),
            uniqueN(samp$dataset)))

files <- samp$abundance_path; names(files) <- samp$sample_id
stopifnot(all(file.exists(files)))
txi <- tximport(files, type = "kallisto", txOut = TRUE,
                countsFromAbundance = "dtuScaledTPM",
                tx2gene = tx2gene[, .(txname, geneid)], ignoreAfterBar = TRUE,
                ignoreTxVersion = FALSE, dropInfReps = TRUE)
cts <- txi$counts; rownames(cts) <- sub("\\|.*$", "", rownames(cts))
tx2gene <- tx2gene[match(rownames(cts), tx2gene$txname)]

# --- DRIMSeq filter ---
n <- nrow(samp); n.small <- min(table(samp$group_binary))
cdf <- data.frame(gene_id = tx2gene$geneid, feature_id = tx2gene$txname,
                  as.data.frame(cts), check.names = FALSE)
sdat <- data.frame(sample_id = samp$sample_id,
                   group = factor(samp$group_binary, levels = c("Control","Disease")),
                   cohort = factor(samp$dataset), sex = factor(samp$sex))
d <- dmFilter(dmDSdata(counts = cdf, samples = sdat),
              min_samps_feature_expr = n.small, min_feature_expr = 10,
              min_samps_feature_prop = n.small, min_feature_prop = 0.10,
              min_samps_gene_expr = n.small, min_gene_expr = 10)
cat(sprintf("[western] after dmFilter: %d tx / %d genes\n",
            nrow(counts(d)), uniqueN(counts(d)$gene_id)))

fc <- as.matrix(counts(d)[, -c(1,2)]); rownames(fc) <- counts(d)$feature_id
se <- SummarizedExperiment(assays = list(counts = fc), colData = DataFrame(sdat),
        rowData = data.frame(isoform_id = counts(d)$feature_id, gene_id = counts(d)$gene_id,
                             row.names = counts(d)$feature_id))
adj <- c(if (nlevels(droplevels(sdat$cohort)) > 1) "cohort",
         if (nlevels(droplevels(sdat$sex))    > 1) "sex")
form <- as.formula(paste(c("~ 0 + group", adj), collapse = " + "))
design <- model.matrix(form, as.data.frame(colData(se)))
L <- matrix(0, ncol(design), 1, dimnames = list(colnames(design), "Dis_vs_Ctrl"))
L["groupDisease",1] <- 1; L["groupControl",1] <- -1
se <- satuRn::fitDTU(object = se, formula = form, parallel = TRUE,
                     BPPARAM = BiocParallel::MulticoreParam(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK","8"))))
se <- satuRn::testDTU(object = se, contrasts = L, diagplot1 = FALSE, diagplot2 = FALSE, sort = FALSE)
res <- rowData(se)[["fitDTUResult_Dis_vs_Ctrl"]]
res$isoform_id <- rownames(res)
res$gene_id <- counts(d)$gene_id[match(rownames(res), counts(d)$feature_id)]
setDT(res)

# --- Delta-proportion (Disease - Control) on tested counts ---
gid_f <- counts(d)$gene_id[match(rownames(fc), counts(d)$feature_id)]
gsum <- rowsum(fc, gid_f); prop <- fc / gsum[match(gid_f, rownames(gsum)), ]
gv <- as.character(sdat$group)
dprop <- setNames(rowMeans(prop[,gv=="Disease",drop=FALSE],na.rm=TRUE)-rowMeans(prop[,gv=="Control",drop=FALSE],na.rm=TRUE), rownames(fc))
res[, dprop := dprop[isoform_id]]
# --- stageR on empirical + theoretical ---
simes <- function(p){p<-sort(p);if(!length(p))return(NA_real_);min(length(p)*p/seq_along(p))}  # sort() drops NA; guard all-NA gene
run_stage <- function(pv){
  names(pv)<-res$isoform_id; pscr<-tapply(pv,res$gene_id,simes); pscr<-setNames(as.numeric(pscr),names(pscr)); pscr<-pscr[!is.na(pscr)]
  so<-stageWiseAdjustment(stageRTx(pScreen=pscr,pConfirmation=matrix(pv,ncol=1,dimnames=list(res$isoform_id,"tx")),
        pScreenAdjusted=FALSE,tx2gene=data.frame(transcript=res$isoform_id,gene=res$gene_id)),method="dtu",alpha=0.05,allowNA=TRUE)
  p<-as.data.table(getAdjustedPValues(so,onlySignificantGenes=FALSE,order=FALSE))
  setnames(p,c("txID","gene","transcript"),c("isoform_id","gscreen","tconf")); p[,.(isoform_id,gscreen,tconf)] }
emp<-run_stage(res$empirical_pval); setnames(emp,c("gscreen","tconf"),c("gene_screen_padj","tx_confirm_padj"))
theo<-run_stage(res$pval); setnames(theo,c("gscreen","tconf"),c("gene_screen_padj_theo","tx_confirm_padj_theo"))
res <- merge(merge(res,emp,by="isoform_id",all.x=TRUE),theo,by="isoform_id",all.x=TRUE)
res <- merge(res, tx2gene[, .(isoform_id = txname, tx_biotype)], by = "isoform_id", all.x = TRUE)
res[, dtu_confirmed_effect := !is.na(tx_confirm_padj) & tx_confirm_padj<0.05 & abs(dprop)>=0.10]
tt<-res$t[is.finite(res$t)]
nullp<-tryCatch({l<-locfdr::locfdr(tt,plot=0);setNames(as.numeric(l$fp0["mlest",c("delta","sigma","p0")]),c("delta","sigma","p0"))},
                error=function(e)c(delta=median(tt),sigma=mad(tt),p0=NA_real_))
null_status<-if(!is.na(nullp["sigma"])&&(nullp["sigma"]>1.5||(!is.na(nullp["p0"])&&nullp["p0"]>1)))"COLLAPSED-use-theoretical" else if(abs(nullp["delta"])>0.5)"SHIFTED" else "ok"
fwrite(res, file.path(RES, sprintf("dtu_mouse_group_cas13%s.tsv", GROUP)), sep = "\t")
ndtu<-uniqueN(res[gene_screen_padj<0.05,gene_id]); ndtuT<-uniqueN(res[gene_screen_padj_theo<0.05,gene_id]); neff<-uniqueN(res[dtu_confirmed_effect==TRUE,gene_id])
cat(sprintf("[western] %s | empirical genes=%d (effect %d) | theoretical=%d | null delta=%.2f sigma=%.2f p0=%.2f -> %s\n",
            GROUP, ndtu, neff, ndtuT, nullp["delta"], nullp["sigma"], nullp["p0"], null_status))
diag_f<-file.path(RES,"dtu_null_diagnostics.tsv")
fwrite(data.table(stratum=sprintf("mouse_group_cas13%s",GROUP),n_tx=nrow(res),emp_dtu_genes=ndtu,emp_effect_genes=neff,theo_dtu_genes=ndtuT,
                  null_delta=round(nullp["delta"],3),null_sigma=round(nullp["sigma"],3),null_p0=round(nullp["p0"],3),
                  emp_frac_lt05=round(mean(res$empirical_pval<0.05,na.rm=TRUE),3),null_status=null_status,
                  primary=ifelse(null_status=="COLLAPSED-use-theoretical","theoretical","empirical")),
       diag_f,sep="\t",append=file.exists(diag_f))
