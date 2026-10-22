# H79 audit: independent recount of every design-side number in PRESPEC_H79.md.
# Reads metadata, ancestry PCs, PNPLA3 calls and the summary columns of deg_results.csv.
# Reads no expression values and fits no expression model.
suppressPackageStartupMessages(library(data.table))
setwd('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
D <- 'RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls'
A <- 'RNA-seq/Human/Patient_Cohorts/results/rna_genotypes_restricted/obs-20260926T155651Z/ancestry/full-20260929T142807Z/library_ancestry_pcs.tsv'
G <- 'Analysis/Hypotheses/H3-20261009T180331Z/results/geno-20261009T180507Z/person_genotypes.tsv'

m <- as.data.table(readRDS(file.path(D, 'meta_matched.rds')))
a <- fread(A)
cat('meta rows', nrow(m), 'unique sample_id', uniqueN(m$sample_id), '\n')
cat('pcs rows', nrow(a), 'unique run', uniqueN(a$run), 'cohorts in pcs file:', paste(sort(unique(a$cohort)), collapse = ','), '\n')
cat('GSE289173 rows in pcs file', sum(a$cohort == 'GSE289173'), '\n')
i <- match(m$sample_id, a$run)
x <- cbind(m, a[i, .(cohort, run_role, possibly_mixed, far_from_centroid, ancestry_call, PC1, PC2, PC3, PC4, axis_EUR_AFR, passes_min_sites, n_sites_used)])
cat('joined', nrow(x), 'unmatched', sum(is.na(i)), 'NA PC1-4', sum(!complete.cases(x[, .(PC1, PC2, PC3, PC4)])), '\n')
cat('dataset == cohort', all(x$dataset == x$cohort), '\n')
cat('run_role x group\n'); print(table(x$run_role, x$group_binary))
cat('possibly_mixed x group\n'); print(table(x$possibly_mixed, x$group_binary))
cat('far_from_centroid x group\n'); print(table(x$far_from_centroid, x$group_binary))
x[, anc3 := fifelse(ancestry_call %in% c('AFR', 'EUR'), ancestry_call, 'other')]
cat('ancestry (AFR/EUR/other) x dataset x group\n'); print(dcast(x[, .N, by = .(dataset, group_binary, anc3)], dataset + group_binary ~ anc3, value.var = 'N', fill = 0))
print(dcast(x[, .N, by = .(group_binary, anc3)], group_binary ~ anc3, value.var = 'N', fill = 0))
cat('AFR fraction control', round(mean(x[group_binary == 'Control']$ancestry_call == 'AFR'), 4),
    'disease', round(mean(x[group_binary == 'Disease']$ancestry_call == 'AFR'), 4), '\n')
cat('n AFR', sum(x$ancestry_call == 'AFR'), 'n EUR', sum(x$ancestry_call == 'EUR'),
    'SE factor sqrt(1/nAFR+1/nEUR)', round(sqrt(1/sum(x$ancestry_call == 'AFR') + 1/sum(x$ancestry_call == 'EUR')), 4), '\n')
cat('inferred_sex x group\n'); print(table(x$inferred_sex, x$group_binary))
e <- x[ancestry_call == 'EUR']
cat('EUR arm dataset x group\n'); print(table(e$dataset, e$group_binary)); cat('EUR n', nrow(e), '\n')
print(table(e$inferred_sex, e$group_binary))

# CMH exact over the 4 cohorts with AFR units
y <- x[dataset != 'GSE126848']
t3 <- table(factor(y$ancestry_call == 'AFR', levels = c(TRUE, FALSE)), factor(y$group_binary, levels = c('Control', 'Disease')), y$dataset)
mh <- mantelhaen.test(t3, exact = TRUE)
cat('CMH exact OR', round(mh$estimate, 3), 'CI', round(mh$conf.int, 3), 'p', signif(mh$p.value, 4), '\n')
for (p in c('PC1', 'PC2', 'PC3', 'PC4')) {
  f <- lm(as.formula(paste(p, '~ dataset + group_binary')), data = x)
  cat(p, 'disease minus control shift', signif(coef(f)['group_binaryDisease'], 4), 't', round(summary(f)$coef['group_binaryDisease', 't value'], 2), '\n')
}
cat('mean PCs by call\n'); print(x[, lapply(.SD, function(v) round(mean(v), 4)), by = ancestry_call, .SDcols = c('PC1', 'PC2', 'PC3', 'PC4')][order(ancestry_call)])
cat('cor PC1 with axis_EUR_AFR', round(cor(x$PC1, x$axis_EUR_AFR), 3), ' SDs', round(sapply(x[, .(PC1, PC2, PC3, PC4)], sd), 4), '\n')

# design
x[, group_binary := factor(group_binary, levels = c('Control', 'Disease'))]
X0 <- model.matrix(~ dataset + inferred_sex + group_binary, x)
X1 <- model.matrix(~ dataset + inferred_sex + PC1 + PC2 + PC3 + PC4 + group_binary, x)
XE <- model.matrix(~ dataset + inferred_sex + group_binary, e)
cat('rank X0', qr(X0)$rank, '/', ncol(X0), ' X1', qr(X1)$rank, '/', ncol(X1), ' XE', qr(XE)$rank, '/', ncol(XE), '\n')
r20 <- summary(lm(as.numeric(group_binary == 'Disease') ~ dataset + inferred_sex, x))$r.squared
r21 <- summary(lm(as.numeric(group_binary == 'Disease') ~ dataset + inferred_sex + PC1 + PC2 + PC3 + PC4, x))$r.squared
cat('R2', round(r20, 4), round(r21, 4), 'VIF', round((1 - r20)/(1 - r21), 4), '\n')

# How much of each PC, after dataset + sex + group, is cohort-by-group cell structure (kept by the
# within-cell permutation) and how much is within-cell variation (destroyed by it).
x[, cell := interaction(dataset, group_binary, drop = TRUE)]
cat('cells', nlevels(x$cell), '\n')
for (p in c('PC1', 'PC2', 'PC3', 'PC4')) {
  r <- resid(lm(as.formula(paste(p, '~ dataset + inferred_sex + group_binary')), x))
  r2c <- summary(lm(r ~ x$cell))$r.squared
  cat(p, 'share of residual variance (after dataset+sex+group) that is cell-level interaction', round(r2c, 4), '\n')
}

# PNPLA3 rs738409 by ancestry call, overall and within disease group
g <- fread(G)
gl <- g[, .(lib = unlist(strsplit(libraries, '[,;|]'))), by = .(person_id, rs738409_call)]
z <- merge(x[, .(sample_id, ancestry_call, group_binary, dataset)], gl, by.x = 'sample_id', by.y = 'lib')
cat('PNPLA3 rows joined', nrow(z), 'unique units', uniqueN(z$sample_id), 'unique persons', uniqueN(z$person_id), '\n')
z <- z[!is.na(rs738409_call)]
z[, Gc := c(RR = 0, RA = 1, AA = 2)[rs738409_call]]
cat('callable', nrow(z), '\n')
print(z[, .(n = .N, Gfreq = round(mean(Gc)/2, 4)), by = .(ancestry_call)][order(ancestry_call)])
print(z[ancestry_call %in% c('AFR', 'EUR'), .(n = .N, Gfreq = round(mean(Gc)/2, 4)), by = .(group_binary, ancestry_call)][order(group_binary, ancestry_call)])
q <- z[ancestry_call %in% c('AFR', 'EUR')]
tb <- rbind(G = tapply(q$Gc, q$ancestry_call, sum), C = tapply(2 - q$Gc, q$ancestry_call, sum))
cat('allele Fisher AFR vs EUR, all units: p', signif(fisher.test(tb)$p.value, 4), '\n'); print(tb)
for (gr in c('Control', 'Disease')) {
  qq <- q[group_binary == gr]
  tb <- rbind(G = tapply(qq$Gc, qq$ancestry_call, sum), C = tapply(2 - qq$Gc, qq$ancestry_call, sum))
  ft <- fisher.test(tb)
  cat('allele Fisher AFR vs EUR within', gr, ': OR', round(ft$estimate, 3), 'p', signif(ft$p.value, 4), '\n'); print(tb)
}
al <-rbind(q[, .(allele = 'G', w = Gc, anc = ancestry_call, grp = as.character(group_binary))],
            q[, .(allele = 'C', w = 2 - Gc, anc = ancestry_call, grp = as.character(group_binary))])
t4 <- xtabs(w ~ factor(allele, levels = c('G', 'C')) + factor(anc, levels = c('AFR', 'EUR')) + grp, al)
mh2 <- mantelhaen.test(t4, exact = TRUE)
cat('allele CMH AFR vs EUR stratified by group: OR', round(mh2$estimate, 3), 'CI', round(mh2$conf.int, 3), 'p', signif(mh2$p.value, 4), '\n')
# person-level (not allele-level) check: carriers of G, AFR vs EUR, stratified by group
q[, carrier := Gc > 0]
t5 <- table(factor(q$carrier, levels = c(TRUE, FALSE)), factor(q$ancestry_call, levels = c('AFR', 'EUR')), q$group_binary)
mh3 <- mantelhaen.test(t5, exact = TRUE)
cat('person-level G carrier CMH AFR vs EUR stratified by group: OR', round(mh3$estimate, 3), 'p', signif(mh3$p.value, 4), '\n')

# Calls and the named genes
d <- fread(file.path(D, 'deg_results.csv'))
call <- d$padj < 0.05 & abs(d$logFC) > 0.5
cat('genes', nrow(d), 'calls', sum(call), 'up', sum(call & d$logFC > 0), 'down', sum(call & d$logFC < 0),
    'calls |lfc|<0.6', sum(call & abs(d$logFC) < 0.6), '\n')
print(d[symbol %in% c('GSTM1', 'ACKR1'), .(symbol, logFC, padj)])

# T2 subset feasibility and empirical rejection thresholds
cat('EUR cell counts (control, disease) and pool sizes\n')
print(merge(x[, .(pool = .N), by = .(dataset, group_binary)], e[, .(eur = .N), by = .(dataset, group_binary)]))
for (B in c(200, 100)) for (al in c(0.025, 0.05)) {
  k <- floor(al * (B + 1) - 1)
  cat(sprintf('B=%d alpha=%.3f: reject iff #{ref >= obs} <= %d, i.e. obs greater than the %d-th largest reference value\n', B, al, k, k + 1))
}
writeLines(capture.output(sessionInfo()), 'Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/audit/H79/audit_recount.sessionInfo.txt')
