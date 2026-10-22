#!/usr/bin/env Rscript
# H3 step 5: DerSimonian-Laird pooling re-derived with metafor::rma on the cohort estimates.
# Usage: 05_verify.R TESTS_DIR
suppressPackageStartupMessages({library(data.table); library(metafor)})
tests <- commandArgs(trailingOnly = TRUE)[[1]]
est <- fread(file.path(tests, "h3_cohort_estimates.tsv"))
pool <- fread(file.path(tests, "h3_pooled.tsv"))
chk <- rbindlist(lapply(seq_len(nrow(pool)), function(i) {
  e <- est[rsid == pool$rsid[i] & program == pool$program[i]]
  m <- rma(yi = e$beta, sei = e$se, method = "DL")
  data.table(gene = pool$gene[i], program = pool$program[i], beta_script = pool$beta[i], beta_metafor = as.numeric(m$beta),
             se_script = pool$se[i], se_metafor = m$se, p_script = pool$p[i], p_metafor = m$pval, tau2_script = pool$tau2[i], tau2_metafor = m$tau2)
}))
print(chk)
stopifnot(max(abs(chk$beta_script - chk$beta_metafor), abs(chk$se_script - chk$se_metafor), abs(chk$p_script - chk$p_metafor), abs(chk$tau2_script - chk$tau2_metafor)) < 1e-8)
h <- fread(file.path(tests, "h31_cohort_estimates.tsv")); hp <- fread(file.path(tests, "h31_pooled.tsv"))
m <- rma(yi = h$beta, sei = h$se, method = "DL")
cat("H31 script beta", hp$beta, "metafor", as.numeric(m$beta), " p", hp$p, m$pval, "\n")
stopifnot(abs(hp$beta - as.numeric(m$beta)) < 1e-8, abs(hp$p - m$pval) < 1e-8)
writeLines(capture.output(sessionInfo()), file.path(tests, "sessionInfo_verify.txt"))
cat("verification passed\n")
