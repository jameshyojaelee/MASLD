#!/usr/bin/env Rscript
# Correct the working EISA ratio model; preserve all November 2 results.
# This is a new analysis definition, not an equivalent reproduction of its stage-only model.
source('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/scripts/s3a_common.R')
library(jsonlite)
a <- parse_args(commandArgs(trailingOnly = TRUE), 'out', defaults = list(permutations = '200'))
out <- new_outdir(a$out); nperm <- as.integer(a$permutations)
SX <- file.path(ROOT, 'Analysis/Hypotheses/SSX-20261102T133212Z')
oldroot <- file.path(SX, 'out_b-20261102T144813Z')
old <- fread(file.path(oldroot, 'eisa_all_filtered_genes.tsv.gz'))
U <- fread(file.path(SS, 'out_s2_qc/units.tsv'))[toupper(as.character(include)) == 'TRUE']
U[, sample_alias := as.character(sample_alias)]
U <- merge(U, read_stage(file.path(ROOT, 'data/external/SteatoSite/fastq/PRJEB58625_fibrosis_stage_20261025T133450Z.tsv')), by = 'sample_alias')[!is.na(stage_code)]
setorder(U, sample_alias); stopifnot(nrow(U) == 493L)
U <- add_labels(U, U$stage_code)
info <- U[, .(sample_alias, batch, inferred_sex, bin)]
fwrite(U[, .(sample_alias, batch, inferred_sex, stage_code, bin)], file.path(out, 'donor_roster.tsv'), sep = '\t')
d <- readRDS(file.path(SS, 'out_s2_qc/dge_included.rds'))
E <- (if (inherits(d, 'DGEList')) d$counts else as.matrix(d))[, U$sample_alias]; rm(d)
I <- do.call(cbind, lapply(U$sample_alias, function(s) {
  z <- fread(file.path(SX, 'out_b_counts', paste0(s, '.intron.counts.gz')), header = FALSE, col.names = c('gene', 'n'))
  setNames(z$n, z$gene)
}))
colnames(I) <- U$sample_alias
libE <- colSums(E); libI <- colSums(I)
nfE <- calcNormFactors(DGEList(E))$samples$norm.factors
nfI <- calcNormFactors(DGEList(I))$samples$norm.factors
keep <- rowMeans(E[rownames(I), ] >= 10) >= .5 & rowMeans(I >= 10) >= .5
K <- rownames(I)[keep]
stopifnot(length(K) == 12010L, setequal(K, old$gene))
logE <- log2(sweep(E[K, ] + .5, 2, libE * nfE, '/') * 1e6)
RR <- logE - log2(sweep(I[K, ] + .5, 2, libI * nfI, '/') * 1e6)
AE <- rowMeans(logE); rm(E, I, logE); invisible(gc())
ratio_fit <- function(bin) {
  bin_values <- factor(bin, levels = BINS)
  z <- copy(info); z[, bin := bin_values]
  stopifnot(identical(z$bin, bin_values))
  X <- design_of(z, 'bin')
  stopifnot(qr(X)$rank == ncol(X))
  f <- lmFit(RR, X); f$Amean <- AE; f <- eBayes(f, trend = TRUE)
  tb <- topTable(f, coef = 'binF3F4', number = Inf, sort.by = 'none')
  data.table(gene = K, dEmI = tb$logFC, t = tb$t, P = tb$P.Value, FDR = p.adjust(tb$P.Value, 'BH'))
}
fit <- ratio_fit(U$bin)
identity_fit <- ratio_fit(as.character(U$bin))
stopifnot(identical(fit, identity_fit))
cmp <- merge(fit, old[, .(gene, saved_dEmI = dEmI_batch_sex, saved_q = eisa_FDR_batch_sex)], by = 'gene')
err <- c(effect = max(abs(cmp$dEmI - cmp$saved_dEmI)), q = max(abs(cmp$FDR - cmp$saved_q)))
stopifnot(all(is.finite(err)), all(err < 1e-12))
fwrite(data.table(quantity = names(err), max_absolute_difference = err), file.path(out, 'saved_adjusted_reproduction.tsv'), sep = '\t')
R <- copy(old)
setnames(R, c('dEmI', 'eisa_FDR'), c('dEmI_stage_only', 'eisa_FDR_stage_only'))
R[, c('dEmI_batch_sex', 'eisa_FDR_batch_sex', 'class') := NULL]
R <- merge(R, fit[, .(gene, dEmI, eisa_FDR = FDR, eisa_P = P)], by = 'gene')
R[, class := fifelse(eisa_FDR < .05, 'exon_intron_ratio_change_supported',
  fifelse(sign(dI) == sign(dE) & dI_FDR < .05, 'intron_change_supported_ratio_unresolved', 'unresolved'))]
fwrite(R, file.path(out, 'eisa_all_filtered_genes.tsv.gz'), sep = '\t')
SC <- fread(file.path(SS, 'out_s3b-20261028T035515Z/steatosite_catalog_sidecar.tsv.gz'))
SC[, gene := gene_id_versioned]
REP <- merge(SC[steatosite_stage_state == 'replicated', .(gene, symbol, primary_evidence_class, S3a_F3F4 = steatosite_F3F4_log2FC)], R[, !'symbol'], by = 'gene', all.x = TRUE)
ann <- fread(file.path(oldroot, 'replicated_gene_classes.tsv.gz'), select = c('gene', 'class'))
REP[is.na(class), class := ann$class[match(gene, ann$gene)]]
stopifnot(nrow(REP) == 7214L, !anyNA(REP$class))
fwrite(REP, file.path(out, 'replicated_gene_classes.tsv.gz'), sep = '\t')
fwrite(REP[, .N, by = class][order(-N)], file.path(out, 'class_counts.tsv'), sep = '\t')
NM <- R[symbol %in% c('GNMT', 'MAT1A', 'COL1A1', 'IGFBP7', 'TNFSF10', 'LUM')]
fwrite(NM, file.path(out, 'named_genes.tsv'), sep = '\t')
set.seed(20261103L)
strat <- paste(info$batch, info$inferred_sex)
null_hits <- integer(nperm)
draw_labels <- matrix(NA_character_, nperm, nrow(U))
first_null_effect_difference <- NA_real_
for (i in seq_len(nperm)) {
  b <- as.character(U$bin)
  for (s in unique(strat)) { j <- which(strat == s); b[j] <- b[j][sample.int(length(j))] }
  draw_labels[i, ] <- b
  stopifnot(all(table(strat, factor(b, levels = BINS)) == table(strat, factor(U$bin, levels = BINS))))
  nf <- ratio_fit(b)
  if (i == 1L) first_null_effect_difference <- max(abs(nf$dEmI - fit$dEmI))
  null_hits[i] <- sum(nf$FDR < .05)
  if (i %% 20L == 0L) cat('Completed adjusted permutations:', i, '/', nperm, '\n')
}
if (nperm) {
  stopifnot(any(draw_labels != rep(as.character(U$bin), each = nperm)),
            first_null_effect_difference > 1e-8,
            !all(null_hits == sum(R$eisa_FDR < .05)))
}
fwrite(data.table(permutation = seq_len(nperm), calls = null_hits), file.path(out, 'perm_null_hits.tsv'), sep = '\t')
fwrite(data.table(permutation = seq_len(nperm), as.data.table(draw_labels)), file.path(out, 'perm_labels.tsv.gz'), sep = '\t')
summ <- list(model = '~ batch + inferred_sex + bin', n_donors = nrow(U), n_genes = nrow(R), ratio_calls = sum(R$eisa_FDR < .05),
  ratio_calls_abs_log2_ge_0_5 = sum(R$eisa_FDR < .05 & abs(R$dEmI) >= .5),
  median_abs_ratio_change_all_calls = median(abs(R[eisa_FDR < .05, dEmI])),
  median_abs_ratio_change_replicated_calls = median(abs(REP[eisa_FDR < .05, dEmI])),
  n_permutations = nperm, permutation_seed = 20261103L,
  p_global_call_count = if (nperm) (1 + sum(null_hits >= sum(R$eisa_FDR < .05))) / (1 + nperm) else NA_real_,
  null_median = if (nperm) median(null_hits) else NA_real_, null_max = if (nperm) max(null_hits) else NA_real_,
  first_null_max_effect_difference = first_null_effect_difference,
  historical_stage_only_calls = sum(R$eisa_FDR_stage_only < .05), shared_calls = sum(R$eisa_FDR < .05 & R$eisa_FDR_stage_only < .05),
  effects_equal_saved_adjusted_within = max(err), classification_note = 'No equivalence or RNA-stability absence is inferred from nonsignificance.')
write_json(summ, file.path(out, 'summary.json'), auto_unbox = TRUE, pretty = TRUE, digits = NA)
writeLines(capture.output(sessionInfo()), file.path(out, 'sessionInfo.txt'))
print(summ)
