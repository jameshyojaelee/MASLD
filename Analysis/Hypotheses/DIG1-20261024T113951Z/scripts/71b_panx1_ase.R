# DIG1 job 71b (written before job 71a counts were read). PANX1 allele-specific expression in bulk MASLD biopsies.
# Background: the PANX1 coloc haplotype (lead rs4593973, hg38 chr11:94144210 A>G; EUR AF 0.10) colocalizes with ALT/AST
# in 6 studies (§43). The alt haplotype lowers liver PANX1 (GTEx liver slope -0.39 at 94150836 G; Broadaway z +9.5 for
# the ref allele). 3'UTR SNPs 94181312 G>C and 94181824 T>G tag it (EUR r +0.74; r 1.00 with each other).
# Calls (per SNP): depth >= 10; alt fraction <= 0.1 -> 0, >= 0.9 -> 2, else 1.
# 3'UTR haplotype per sample: both SNPs het -> counts summed (in phase, r 1.00); one het -> that SNP.
# P1 (primary): in 3'UTR heterozygotes with >= 20 reads, the alt allele is under-represented. Statistic: median per-sample
#   alt fraction, two-sided Wilcoxon against 0.5. Mapping-bias control: exon-1 SNP 94128854 C>T, whose alt allele is on
#   the opposite, higher-expression side (GTEx slope +0.21; r -0.26 with the lead), so generic reference bias would push
#   both below 0.5. P1 passes if the 3'UTR median < 0.5 (p < 0.05) and lies below the exon-1 median (Wilcoxon rank-sum
#   p < 0.05).
# S (Holm over 4, two-sided): is the allelic imbalance in 3'UTR hets larger with the inflamed state? Quasibinomial GLM of
#   alt/ref counts on INFLAMED_HEP, then disease group, NAS, fibrosis stage separately; each + dataset.
# D (descriptive): PANX1 log2 CPM ~ 3'UTR alt dosage + dataset + sex (HC3); INFLAMED_HEP, NAS, fibrosis ~ dosage +
#   dataset + sex. Intronic SNPs are used as genotypes only if call rate >= 0.8 (expected low in poly(A) RNA).
# Unit: one biopsy per participant. No ancestry covariates. Cross-sectional.
# DEVIATION (after the first run crashed; only site_summary had printed): the exon-1 control SNP has no reads
#   (call rate 0), so P1 has no mapping-bias control and is reported as descriptive. The S tests compare samples at the
#   same SNP, so a constant mapping bias does not affect them.
# Usage: Rscript 71b_panx1_ase.R <job71a_out> <new_out_dir>
suppressPackageStartupMessages({library(data.table); library(edgeR)})
args <- commandArgs(TRUE); stopifnot(length(args) == 2L, !dir.exists(args[2])); out <- args[2]; dir.create(out, recursive = TRUE)
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"; RUN <- file.path(ROOT, "Analysis/Hypotheses/DIG1-20261024T113951Z")
A <- fread(file.path(args[1], "allele_counts.tsv"), header = FALSE, col.names = c("sample_id", "pos", "ref", "alt", "ad"))
EXP <- data.table(pos = c(94128854, 94132981, 94137172, 94144210, 94150836, 94181312, 94181824),
                  eref = c("C", "T", "C", "A", "C", "G", "T"), ealt = c("T", "G", "T", "G", "G", "C", "G"))
A <- merge(A, EXP, by = "pos")
cnt <- function(ad, alt, ealt) { v <- as.numeric(strsplit(ad, ",")[[1]]); al <- strsplit(alt, ",")[[1]]
  k <- match(ealt, al); c(v[1], if (is.na(k)) 0 else v[k + 1]) }
cc <- t(mapply(cnt, A$ad, A$alt, A$ealt)); A[, `:=`(nref = cc[, 1], nalt = cc[, 2])]
stopifnot(all(A$ref == A$eref))
D <- fread(file.path(RUN, "out_genotype_state/sample_scores_genotypes.tsv"))
A <- A[sample_id %in% D$sample_id]
A[, dp := nref + nalt][, f := nalt / dp]
A[, call := fifelse(dp < 10, NA_real_, fifelse(f <= 0.1, 0, fifelse(f >= 0.9, 2, 1)))]
calls <- dcast(A, sample_id ~ pos, value.var = "call"); setnames(calls, -1, paste0("g", names(calls)[-1]))
cr <- A[, .(n = .N, call_rate = sum(!is.na(call)) / nrow(D), af = mean(call, na.rm = TRUE) / 2,
            n_het = sum(call == 1, na.rm = TRUE), med_dp = median(dp)), by = pos]
fwrite(cr, file.path(out, "site_summary.tsv"), sep = "\t"); print(cr)
# 3'UTR haplotype
U <- A[pos %in% c(94181312, 94181824) & call == 1]
U <- U[, .(nref = sum(nref), nalt = sum(nalt), n_snps = .N), by = sample_id][nref + nalt >= 20]
U[, f := nalt / (nref + nalt)]
X <- A[pos == 94128854 & call == 1 & dp >= 20, .(sample_id, nref, nalt, f)]
w1 <- wilcox.test(U$f, mu = 0.5); NAp <- list(p.value = NA_real_)
w2 <- if (nrow(X) > 1) wilcox.test(X$f, mu = 0.5) else NAp; w3 <- if (nrow(X) > 1) wilcox.test(U$f, X$f) else NAp
P1 <- data.table(test = c("3UTR hets vs 0.5", "exon1 hets vs 0.5", "3UTR vs exon1"), n = c(nrow(U), nrow(X), NA),
                 median_f = c(median(U$f), if (nrow(X)) median(X$f) else NA, NA), pooled_f = c(sum(U$nalt) / sum(U$nref + U$nalt), sum(X$nalt) / sum(X$nref + X$nalt), NA),
                 p = c(w1$p.value, w2$p.value, w3$p.value))
P1[, pass := NA]
fwrite(P1, file.path(out, "P1_ase.tsv"), sep = "\t"); print(P1)
# S: imbalance vs state
U <- merge(U, D, by = "sample_id")
S <- rbindlist(lapply(c("INFLAMED_HEP", "group_binary", "nas_score", "fibrosis_stage"), function(v) {
  d <- U[!is.na(get(v))]; d[, x := if (is.character(get(v))) as.numeric(get(v) != sort(unique(get(v)))[1]) else as.numeric(get(v))]
  if (length(unique(d$dataset)) > 1) fit <- glm(cbind(nalt, nref) ~ x + dataset, family = quasibinomial, data = d)
  else fit <- glm(cbind(nalt, nref) ~ x, family = quasibinomial, data = d)
  co <- summary(fit)$coefficients["x", ]
  data.table(term = v, n = nrow(d), beta_logit = co[1], se = co[2], p = co[4], level1 = if (is.character(U[[v]])) sort(unique(na.omit(U[[v]])))[2] else "")
}))
S[, holm := p.adjust(p, "holm")]; fwrite(S, file.path(out, "S_imbalance_vs_state.tsv"), sep = "\t"); print(S)
# D: dosage models
REL <- file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls")
d <- readRDS(file.path(REL, "merged_dge.rds")); stopifnot(identical(colnames(d$counts), as.character(D$sample_id)))
d <- calcNormFactors(d); lc <- cpm(d, log = TRUE, prior.count = 2)
D$PANX1 <- lc[grep("^ENSG00000110218", rownames(lc)), ]
U2 <- A[pos %in% c(94181312, 94181824), .(dos = if (all(is.na(call))) NA_real_ else max(call, na.rm = TRUE)), by = sample_id]
D <- merge(D, U2, by = "sample_id", all.x = TRUE)
hc3 <- function(fit) { X <- model.matrix(fit); h <- hatvalues(fit); e <- residuals(fit); B <- solve(crossprod(X))
  sqrt(diag(B %*% crossprod(X * (e / (1 - h))) %*% B)) }
DM <- rbindlist(lapply(c("PANX1", "INFLAMED_HEP", "nas_score", "fibrosis_stage"), function(y) {
  dd <- D[!is.na(get(y)) & !is.na(dos) & !is.na(sex_final)]
  fit <- lm(as.formula(paste(y, "~ dos + dataset + sex_final")), data = dd); se <- hc3(fit)
  b <- coef(fit)["dos"]; s <- se[names(coef(fit)) == "dos"]
  data.table(outcome = y, n = nrow(dd), n_carriers = sum(dd$dos > 0), beta = b, se = s, p = 2 * pnorm(-abs(b / s)))
}))
fwrite(DM, file.path(out, "D_dosage_models.tsv"), sep = "\t"); print(DM)
fwrite(U, file.path(out, "utr_het_counts.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
