set.seed(20260923)
setwd('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
D <- 'RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls'
m <- readRDS(file.path(D,'meta_matched.rds'))
cat('meta rows', nrow(m), '\n'); print(names(m))
a <- read.delim('RNA-seq/Human/Patient_Cohorts/results/rna_genotypes_restricted/obs-20260926T155651Z/ancestry/full-20260929T142807Z/library_ancestry_pcs.tsv')
x <- merge(m, a[,c('run','ancestry_call','axis_EUR_AFR','PC1','PC2','PC3','PC4','n_sites_used')], by.x='sample_id', by.y='run', all.x=TRUE)
cat('joined', nrow(x), 'missing call', sum(is.na(x$ancestry_call)), '\n')
print(table(x$ancestry_call, x$group_binary, useNA='ifany'))
print(table(x$dataset, x$group_binary))
print(table(x$dataset, x$ancestry_call=='AFR', x$group_binary))
y <- x[x$dataset!='GSE126848',]
t <- table(factor(y$ancestry_call=='AFR',levels=c(TRUE,FALSE)), factor(y$group_binary,levels=c('Control','Disease')), y$dataset)
print(mantelhaen.test(t, exact=TRUE))
# PC difference within cohort
for (p in c('PC1','PC2','PC3','PC4')) { f <- lm(as.formula(paste(p,'~ dataset + group_binary')), data=x); print(p); print(coef(summary(f))['group_binaryDisease',]) }
summary(x$n_sites_used)
# canonical calls
d <- read.csv(file.path(D,'deg_results.csv'))
cat('genes', nrow(d), 'calls', sum(d$padj<0.05 & abs(d$logFC)>0.5, na.rm=TRUE), '\n')
print(d[d$symbol %in% c('GSTM1','GSTT1','UGT2B17','HLA-DRB5','GSTT2B','CYP2D6','APOL1','DARC','ACKR1'), c('symbol','logFC','padj')])
# residual SD median for MDE
w <- read.delim(file.path(D,'lvqw_sample_weights.tsv')); print(head(w,2))
# PNPLA3 positive control
g <- read.delim('Analysis/Hypotheses/H3-20261009T180331Z/results/geno-20261009T180507Z/person_genotypes.tsv')
lib <- strsplit(g$libraries, '[,;|]')
gl <- data.frame(lib=unlist(lib), call=rep(g$rs738409_call, lengths(lib)))
z <- merge(x[,c('sample_id','ancestry_call')], gl, by.x='sample_id', by.y='lib')
z <- z[!is.na(z$call),]
z$G <- c(RR=0,RA=1,AA=2)[z$call]
cat('PNPLA3 callable joined', nrow(z), '\n')
print(aggregate(G~ancestry_call, z, function(v) c(n=length(v), Gfreq=mean(v)/2)))
af <- tapply(z$G, z$ancestry_call, function(v) mean(v)/2)
# AFR vs EUR allele-count test
s <- z[z$ancestry_call %in% c('AFR','EUR','AMR'),]
for (k in c('AFR','AMR')) { q <- s[s$ancestry_call %in% c(k,'EUR'),]; tb <- rbind(tapply(q$G, q$ancestry_call, sum), tapply(2-q$G, q$ancestry_call, sum)); print(k); print(tb); print(fisher.test(tb)$p.value) }
sessionInfo()
