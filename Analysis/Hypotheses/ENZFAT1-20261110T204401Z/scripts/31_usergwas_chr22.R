#!/usr/bin/env Rscript
# ENZFAT1 S3 step 2: Genomic SEM userGWAS on chromosome 22 HapMap3 SNPs, compared with the closed-form z.
# Covariance structure: the alt_ukb and pdff block of the S1 9-trait ldsc() output (same S, V, I).
# Reference MAF: eur_w_ld_chr chromosome 22 MAF with hm3 alleles (no 1000G reference file on disk).
# Pass rule (prespecified): correlation of z above 0.99.
# Usage: 31_usergwas_chr22.R <S1 out dir> <S3 out dir> <cores>
suppressMessages({library(GenomicSEM); library(data.table)})
args <- commandArgs(trailingOnly = TRUE)
s1 <- args[1]; out <- args[2]; cores <- as.integer(args[3])
set.seed(20261110)
setwd(out)
cov9 <- readRDS(file.path(s1, "ldsc_all.rds"))
tr <- colnames(cov9$S); t2 <- match(c("pdff", "alt_ukb"), tr)
idx <- matrix(NA_integer_, length(tr), length(tr)); idx[lower.tri(idx, diag = TRUE)] <- seq_len(sum(lower.tri(idx, diag = TRUE)))
idx[upper.tri(idx)] <- t(idx)[upper.tri(idx)]
sel <- c(idx[t2[1], t2[1]], idx[t2[2], t2[1]], idx[t2[2], t2[2]])
# userGWAS reads covstruc by position: [[1]] V, [[2]] S, [[3]] I (as ldsc() returns it)
cs <- list(V = cov9$V[sel, sel], S = cov9$S[t2, t2], I = cov9$I[t2, t2], N = cov9$N[sel], m = cov9$m)
dimnames(cs$S) <- list(c("pdff", "alt_ukb"), c("pdff", "alt_ukb")); colnames(cs$S) <- c("pdff", "alt_ukb")
print(cs$S); print(cs$I)
snps <- sumstats(files = c(file.path(out, "chr22_pdff.tsv"), file.path(out, "chr22_alt_ukb.tsv")),
                 ref = file.path(out, "chr22_ref_eur_w_ld_hm3.tsv"), trait.names = c("pdff", "alt_ukb"),
                 se.logit = c(FALSE, FALSE), OLS = c(TRUE, TRUE), linprob = c(FALSE, FALSE))
cat("SNPs after sumstats():", nrow(snps), "\n")
model <- "F1 =~ NA*pdff + alt_ukb
F2 =~ NA*alt_ukb
F1 ~~ 1*F1
F2 ~~ 1*F2
F1 ~~ 0*F2
pdff ~~ 0*pdff
alt_ukb ~~ 0*alt_ukb
pdff ~~ 0*alt_ukb
F1 ~ SNP
F2 ~ SNP"
g <- userGWAS(covstruc = cs, SNPs = snps, estimation = "DWLS", model = model, sub = c("F2~SNP"),
              cores = cores, parallel = TRUE)
g <- as.data.table(g[[1]])
fwrite(g, file.path(out, "usergwas_chr22_F2.tsv"), sep = "\t")
cf <- fread(file.path(out, "alt_not_fat_component_GRCh37.tsv.gz"))[CHR == 22]
m <- merge(g[, .(SNP, A1u = A1, A2u = A2, z_user = Z_Estimate, est_user = est)], cf[, .(SNP, A1, A2, Z)], by = "SNP")
m[, z_cf := ifelse(A1u == A1, Z, ifelse(A1u == A2, -Z, NA_real_))]
m <- m[is.finite(z_user) & is.finite(z_cf)]
r <- cor(m$z_user, m$z_cf)
chk <- data.table(snps_compared = nrow(m), cor_z = r, slope_z = coef(lm(z_user ~ z_cf, m))[2],
                  max_abs_diff = max(abs(m$z_user - m$z_cf)), pass_cor_gt_0.99 = r > 0.99)
fwrite(chk, file.path(out, "usergwas_chr22_check.tsv"), sep = "\t"); print(chk)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
