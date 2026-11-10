# SS S3h (F-SS8, PRESPEC_SS section 5 and PRESPEC_SS_addendum_S3h.md; exploratory; fork of DIG1 jobs 71b/72b and
# scripts/analysis/in_biopsy_allelic/a5_allelic). Stage-dependent allelic expression in the 426 genetically anchored or
# concordant r8 Catalog genes, from allele counts at the S0d sites (s3_allele_counts.py).
# Het (jobs 71b/72b): reads for the expected REF + ALT >= 10 and 0.1 < ALT fraction < 0.9. C>T and G>A sites are not
#   in the S0d list.
# One site per gene, chosen without stage: most het libraries among all included libraries; ties: higher median het
#   depth, then lower position. Testable: >= 30 staged het individuals with inferred_sex at that site.
# Stage test (family; BH over testable genes): quasibinomial GLM cbind(ALT, REF) ~ stage_num + batch + inferred_sex
#   (single-level terms dropped), two-sided t on residual df. Mapping-bias null: the S0d null sites (non-catalog
#   protein-coding non-eGenes, one SNV per gene) run through the same het, testability and model rules;
#   lambda_null = median(t^2) / qchisq(0.5, 1). If lambda_null > 1 every catalog t is divided by sqrt(lambda_null)
#   before the normal p (p_gc); p_gc is the family p. Beside: empirical p against null sites in the same median-het-depth
#   quintile x EUR MAF bin (0.05-0.15, 0.15-0.30, 0.30-0.50; depth quintile only if the cell has < 20 null sites).
# Descriptive: intercept model (~ 1): pooled ALT fraction against the matched null-site pooled fractions (job 72b
#   design); folded imbalance |log((ALT + 0.5) / (REF + 0.5))| ~ stage_num + batch + sex (OLS, HC3), which does not
#   depend on the unknown phase between tag and regulatory variant.
# Unit: genetic individual (S2 Q3). One site per gene, so genes are tested once. Cross-sectional.
# Usage: Rscript s3h_stage_ase.R --input <S2 dir> --stage <stage tsv> --counts <s3_allele_counts dir>
#        --targets <out_s0d_targets dir> --out <new dir>
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(trailingOnly = TRUE)
opt <- function(k) { i <- match(k, a); if (is.na(i)) NA_character_ else a[i + 1] }
inp <- opt("--input"); stf <- opt("--stage"); cdir <- opt("--counts"); tdir <- opt("--targets"); out <- opt("--out")
stopifnot(!anyNA(c(inp, stf, cdir, tdir, out)), !dir.exists(out)); dir.create(out, recursive = TRUE)
SS <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z"
stopifnot(system(paste("cd", file.path(SS, "frozen"), "&& sha256sum -c --quiet SHA256SUMS")) == 0)
for (l in readLines(file.path(tdir, "SHA256SUMS"))) { h <- sub(" .*", "", l); f <- sub("^[0-9a-f]+  ", "", l)
  stopifnot(sub(" .*", "", system2("sha256sum", f, stdout = TRUE)) == h) }
U <- fread(file.path(inp, "units.tsv"))[toupper(as.character(include)) == "TRUE"]
U[, sample_alias := as.character(sample_alias)]
st <- fread(cmd = paste("tr -d '\\r' <", shQuote(stf)), colClasses = "character")
U[, stage_num := suppressWarnings(as.numeric(st$nash_crn_kleiner_fibrosis_stage[match(sample_alias, st$sample_alias)]))]
TG <- fread(file.path(tdir, "s1b_targets_annotated.tsv"))[role %in% c("s3h", "s3h_null")]
TG <- TG[, .(chrom, pos, role, gene, id, af_eur, maf = pmin(af_eur, 1 - af_eur))]
A <- fread(file.path(cdir, "ase_counts.tsv.gz"), colClasses = c(sample_alias = "character"))
A <- merge(A[sample_alias %in% U$sample_alias], TG, by = c("chrom", "pos"))
A[, dp := nref + nalt][, f := nalt / dp]
H <- A[dp >= 10 & f > 0.1 & f < 0.9]
H <- merge(H, U[, .(sample_alias, batch, inferred_sex, stage_num)], by = "sample_alias")
H[, unit := fifelse(role == "s3h", id, paste0("null:", chrom, ":", pos))]
site <- H[, .(n_het_all = .N, med_dp = as.numeric(median(dp)), n_het_staged = sum(!is.na(stage_num) & !is.na(inferred_sex))),
          by = .(unit, role, gene, chrom, pos, maf)]
setorder(site, unit, -n_het_all, -med_dp, pos)
chosen <- site[, .SD[1], by = unit]
chosen[, testable := n_het_staged >= 30]
fwrite(site, file.path(out, "site_het_counts.tsv"), sep = "\t")
cat_ids <- unique(TG[role == "s3h", id])
miss <- data.table(unit = setdiff(cat_ids, chosen[role == "s3h", unit]))
used_cov <- function(z, vars) vars[vapply(vars, function(v) length(unique(z[[v]])) > 1, logical(1))]
one <- function(k) {
  z <- H[unit == chosen$unit[k] & pos == chosen$pos[k] & !is.na(stage_num) & !is.na(inferred_sex)]
  cv <- used_cov(z, c("batch", "inferred_sex"))
  res <- data.table(unit = chosen$unit[k], n = nrow(z), beta = NA_real_, se = NA_real_, t = NA_real_, df = NA_real_,
                    pooled_f = sum(z$nalt) / sum(z$dp), int_logit = NA_real_, int_p = NA_real_, fold_beta = NA_real_, fold_p = NA_real_)
  if (length(unique(z$stage_num)) < 2) return(res)
  g <- glm(reformulate(c("stage_num", cv), "cbind(nalt, nref)"), family = quasibinomial, data = z)
  co <- summary(g)$coefficients["stage_num", ]
  res[, `:=`(beta = co[1], se = co[2], t = co[3], df = g$df.residual)]
  g0 <- summary(glm(cbind(nalt, nref) ~ 1, family = quasibinomial, data = z))$coefficients[1, ]
  res[, `:=`(int_logit = g0[1], int_p = g0[4])]
  z[, yfold := abs(log((nalt + 0.5) / (nref + 0.5)))]
  m <- lm(reformulate(c("stage_num", cv), "yfold"), data = z); X <- model.matrix(m); e <- residuals(m); h <- hatvalues(m)
  XtXi <- solve(crossprod(X)); V <- XtXi %*% crossprod(X * (e / (1 - h))) %*% XtXi
  res[, `:=`(fold_beta = coef(m)[["stage_num"]], fold_p = 2 * pt(-abs(coef(m)[["stage_num"]] / sqrt(V["stage_num", "stage_num"])), m$df.residual))]
  res
}
R <- rbindlist(lapply(which(chosen$testable), one))
R <- merge(chosen, R, by = "unit", all.x = TRUE)
R[, p := 2 * pt(-abs(t), df)]
NL <- R[role == "s3h_null" & is.finite(t)]
lambda <- median(NL$t^2) / qchisq(0.5, 1)
R[, p_gc := if (lambda > 1) 2 * pnorm(-abs(t) / sqrt(lambda)) else 2 * pnorm(-abs(t))]
# matched null (depth quintile x MAF bin)
TT <- R[is.finite(t)]
brk <- unique(quantile(TT$med_dp, seq(0, 1, 0.2))); R[, dbin := cut(med_dp, brk, include.lowest = TRUE, labels = FALSE)]
R[, mbin := cut(maf, c(0.05, 0.15, 0.30, 0.50), include.lowest = TRUE, right = FALSE, labels = FALSE)]
nullp <- function(x, db, mb, col, centre) {
  pool <- R[role == "s3h_null" & is.finite(get(col)) & dbin == db & mbin == mb][[col]]
  if (length(pool) < 20) pool <- R[role == "s3h_null" & is.finite(get(col)) & dbin == db][[col]]
  if (!is.finite(x) || !length(pool)) return(c(NA_real_, length(pool)))
  m <- if (centre) median(pool) else 0
  c((1 + sum(abs(pool - m) >= abs(x - m))) / (1 + length(pool)), length(pool)) }
R[, plogit := qlogis(pmin(pmax(pooled_f, 1e-4), 1 - 1e-4))]
R[, c("p_emp_stage", "n_null_stage", "p_emp_intercept", "n_null_intercept") := NA_real_]
for (i in which(R$role == "s3h" & is.finite(R$t))) {
  e1 <- nullp(R$t[i], R$dbin[i], R$mbin[i], "t", FALSE); e2 <- nullp(R$plogit[i], R$dbin[i], R$mbin[i], "plogit", TRUE)
  set(R, i, c("p_emp_stage", "n_null_stage", "p_emp_intercept", "n_null_intercept"), as.list(c(e1, e2)))
}
F <- R[role == "s3h" & testable == TRUE & is.finite(t)]
F[, bh := p.adjust(p_gc, "BH")][, bh_fold := p.adjust(fold_p, "BH")]
fwrite(F[order(p_gc)], file.path(out, "f_ss8_catalog_genes.tsv"), sep = "\t")
fwrite(R[role == "s3h_null"], file.path(out, "null_sites.tsv"), sep = "\t")
summ <- c(paste("staged units with sex:", sum(!is.na(U$stage_num) & !is.na(U$inferred_sex))),
          paste("catalog genes in S0d universe:", length(cat_ids), "; with a het site:", sum(chosen$role == "s3h"),
                "; no het site:", nrow(miss)),
          paste("testable catalog genes:", nrow(F), "; null sites testable:", nrow(NL), "of", sum(chosen$role == "s3h_null")),
          paste("lambda_null:", signif(lambda, 3), "; null pooled ALT fraction median:", signif(median(NL$pooled_f), 3)),
          paste("catalog genes BH < 0.05 (p_gc):", sum(F$bh < 0.05), "; folded BH < 0.05:", sum(F$bh_fold < 0.05, na.rm = TRUE)))
writeLines(summ, file.path(out, "summary.txt"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
options(width = 220); cat(summ, sep = "\n"); print(head(F[order(p_gc)], 30), digits = 3)
