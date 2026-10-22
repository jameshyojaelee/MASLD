# H55 planning: semi-synthetic power for the prespecified gene family (canonical & filterByExpr, 1,002 genes)
# on TMM logCPM of GTEx v8 liver counts. The exposure is ALWAYS permuted within strata (center x Hardy group)
# before use; the unpermuted exposure is never related to expression, so no observed statistic is computed.
# Signal injected: H8 pooled steatosis-grade logFC (0 where H8 did not test), times attenuation a.
# Exposure proxy: note ordinal (no_mention 0, unspecified 2), as in gate-L02-01/simpower.py.
suppressMessages({library(edgeR); library(data.table)})
set.seed(20260923)
B <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
R <- paste0(B, "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/")
S <- paste0(R, "scratch/plan-H55/")
meta <- fread(paste0(S, "design_meta_counts_only.tsv")); meta <- meta[!is.na(SMCENTER) & SMCENTER != ""]
ng <- fread(paste0(R, "scratch/gate-L02-01/note_grades.tsv"))
meta <- merge(meta, ng[, .(donor, ord)], by = "donor", sort = FALSE)
meta[, x := ifelse(grade == "no_mention", 0, ifelse(grade == "unspecified", 2, ord))]
meta[, ageg := fcase(AGE %in% c("20-29","30-39"), "20-39", AGE == "40-49", "40-49", AGE == "50-59", "50-59", default = "60-79")]
meta[, hg := fcase(H == 0, "H0", H %in% 1:2, "H12", default = "H34")]
meta[, stratum := paste(SMCENTER, hg)]
print(table(meta$stratum))
cnt <- fread(paste0(S, "inputs/gene_reads_2017-06-05_v8_liver.gct.gz"), skip = 2)
cnt <- cnt[!grepl("_PAR_Y$", Name)]
M <- as.matrix(cnt[, meta$SAMPID, with = FALSE]); rownames(M) <- sub("\\..*$", "", cnt$Name)
y <- DGEList(M); keep <- filterByExpr(y); y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE])
L <- cpm(y, log = TRUE, prior.count = 2)
fam <- fread(paste0(S, "primary_family_genes.tsv")); pc <- fread(paste0(S, "pc_composite_genes.tsv"))
h8 <- fread(paste0(B, "Analysis/Hypotheses/H8-20261009T015529Z/run-20261009T020920Z/meta_steatosis.tsv.gz"))
bmap <- setNames(h8$mu, h8$gene_id)
design_for <- function(d) {
  f <- ~ ageg + factor(SEX) + factor(H) + SMTSISCH + SMRIN + SMCENTER
  if (length(unique(d$H)) == 1) f <- ~ ageg + factor(SEX) + SMTSISCH + SMRIN + SMCENTER
  model.matrix(f, data = d)
}
perm_x <- function(d) { xp <- d$x; for (s in unique(d$stratum)) { i <- which(d$stratum == s); xp[i] <- d$x[i][sample.int(length(i))] }; xp }
run <- function(label, idx, genes, lfc, nsim = 400, atts = c(1, 0.5, 0.25)) {
  d <- meta[idx]; Z <- design_for(d); Q <- qr.Q(qr(Z))
  Y <- t(L[genes, idx, drop = FALSE]); beta <- bmap[genes]; beta[is.na(beta)] <- 0
  Yr0 <- Y - Q %*% crossprod(Q, Y)
  stat <- function(Yr, xp) { xr <- xp - Q %*% crossprod(Q, xp); b <- drop(crossprod(xr, Yr)) / sum(xr^2); cor(lfc, b, method = "spearman") }
  null <- numeric(nsim); alt <- matrix(NA, nsim, length(atts))
  for (k in seq_len(nsim)) {
    xp <- perm_x(d); null[k] <- stat(Yr0, xp)
    for (j in seq_along(atts)) { Ya <- Y + outer(xp, atts[j] * beta); Yr <- Ya - Q %*% crossprod(Q, Ya); alt[k, j] <- stat(Yr, xp) }
  }
  lo <- quantile(null, 0.025); hi <- quantile(null, 0.975)
  cat(sprintf("%s: donors %d, design cols %d, genes %d (H8-tested %d); null sd %.3f [%.3f, %.3f]; MDE~2.8sd %.2f", label, nrow(d), ncol(Z) + 1, length(genes), sum(genes %in% h8$gene_id), sd(null), lo, hi, 2.8 * sd(null)))
  for (j in seq_along(atts)) cat(sprintf("; a=%.2f mean rho %.3f power %.2f", atts[j], mean(alt[, j]), mean(alt[, j] > hi | alt[, j] < lo)))
  cat("\n")
}
all <- seq_len(nrow(meta))
tpm <- fread(paste0(B, "data/external/allelic_refs/gtex_v8/gene_tpm_2017-06-05_v8_liver.gct.gz"), skip = 2)
tpm <- tpm[!grepl("_PAR_Y$", Name)]; Tm <- as.matrix(tpm[, meta$SAMPID, with = FALSE]); rownames(Tm) <- sub("\\..*$", "", tpm$Name)
med <- apply(Tm[fam$ens, ], 1, median)
f2 <- fam[med >= 1]
run("SENS family median TPM>=1 within filterByExpr family (logCPM)", all, f2$ens, f2$bulk_logFC)
ex <- fread(paste0(S, "exclusion_set_tnfa_acutephase.tsv"))
f3 <- fam[!(ens %in% ex$ensembl_gene)]
run("SECONDARY family minus TNFA/acute-phase set", all, f3$ens, f3$bulk_logFC)
