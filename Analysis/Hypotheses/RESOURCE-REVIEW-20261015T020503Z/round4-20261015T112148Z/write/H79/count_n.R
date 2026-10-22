# H79 prespec counts: design-side only. Reads no expression and no fit output beyond deg_results.csv summary columns.
setwd('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
D <- 'RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls'
A <- 'RNA-seq/Human/Patient_Cohorts/results/rna_genotypes_restricted/obs-20260926T155651Z/ancestry/full-20260929T142807Z/library_ancestry_pcs.tsv'
m <- readRDS(file.path(D,'meta_matched.rds'))
a <- read.delim(A, stringsAsFactors=FALSE)
cat('pcs rows', nrow(a), 'unique run', length(unique(a$run)), '\n')
x <- merge(m, a, by.x='sample_id', by.y='run', all.x=TRUE, suffixes=c('','.anc'))
cat('joined', nrow(x), 'NA PC1', sum(is.na(x$PC1)), 'NA call', sum(is.na(x$ancestry_call)), 'NA inferred_sex', sum(is.na(x$inferred_sex)), '\n')
cat('dataset == cohort', all(x$dataset==x$cohort), '\n')
for (v in c('run_role','exclusion_reason','passes_min_sites','possibly_mixed','cohort_contaminated','far_from_centroid','intermediate')) { cat('--', v, '\n'); print(table(x[[v]], x$group_binary, useNA='ifany')) }
cat('-- ancestry_call x dataset x group\n'); print(ftable(table(x$dataset, x$group_binary, x$ancestry_call)))
cat('-- mean PCs by call\n'); print(aggregate(cbind(PC1,PC2,PC3,PC4,axis_EUR_AFR)~ancestry_call, x, function(v) round(mean(v),4)))
cat('-- cor PC1..4 with axis_EUR_AFR\n'); print(round(cor(x[,c('PC1','PC2','PC3','PC4')], x$axis_EUR_AFR),3))
cat('-- PC SD\n'); print(round(sapply(x[,c('PC1','PC2','PC3','PC4')], sd),4))
# design collinearity (no expression)
x$group_binary <- factor(x$group_binary, levels=c('Control','Disease'))
X0 <- model.matrix(~ dataset + inferred_sex + group_binary, x)
X1 <- model.matrix(~ dataset + inferred_sex + PC1+PC2+PC3+PC4 + group_binary, x)
cat('rank X0', qr(X0)$rank, ncol(X0), 'rank X1', qr(X1)$rank, ncol(X1), '\n')
r2 <- summary(lm(as.numeric(group_binary=='Disease') ~ dataset + inferred_sex + PC1+PC2+PC3+PC4, x))$r.squared
r20 <- summary(lm(as.numeric(group_binary=='Disease') ~ dataset + inferred_sex, x))$r.squared
cat('R2 group on base', round(r20,4), 'with PCs', round(r2,4), 'unweighted variance inflation of group coef', round((1-r20)/(1-r2),4), '\n')
# EUR-only arm
e <- x[x$ancestry_call=='EUR',]
cat('-- EUR arm dataset x group\n'); print(table(e$dataset, e$group_binary)); cat('EUR total', nrow(e), '\n')
cat('EUR sex x group\n'); print(table(e$inferred_sex, e$group_binary))
XE <- model.matrix(~ dataset + inferred_sex + group_binary, e); cat('rank XE', qr(XE)$rank, ncol(XE), '\n')
cat('-- inferred_sex x group (all)\n'); print(table(x$inferred_sex, x$group_binary))
# canonical calls near threshold
d <- read.csv(file.path(D,'deg_results.csv'))
call <- d$padj<0.05 & abs(d$logFC)>0.5
cat('genes', nrow(d), 'calls', sum(call), 'up', sum(call & d$logFC>0), 'down', sum(call & d$logFC<0), '\n')
cat('calls |lfc|<0.6', sum(call & abs(d$logFC)<0.6), '<0.75', sum(call & abs(d$logFC)<0.75), ' calls padj>0.01', sum(call & d$padj>0.01), '\n')
cat('median SE of calls', round(median(d$SE[call]),4), '\n')
cat('non-calls padj<0.05 & 0.4<|lfc|<=0.5', sum(d$padj<0.05 & abs(d$logFC)>0.4 & abs(d$logFC)<=0.5, na.rm=TRUE), '\n')
print(d[d$symbol %in% c('GSTM1','UGT2B17','ACKR1','PNPLA3','GSTT1'), c('gene','symbol','logFC','SE','padj','AveExpr')])
writeLines(capture.output(sessionInfo()), '/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/write/H79/count_n.sessionInfo.txt')
