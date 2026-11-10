# SS S3f step 2 (F-SS6, PRESPEC_SS section 5; fork of DIG1 jobs 63, 71b, 72b). Genotype family in SteatoSITE.
# Genotypes: s3f_genotype_table.py (job 62b calls and QC, HWE once; HSD17B13 junction readout of job 65).
# Score INFLAMED_HEP = FCAMR, TRIM31, CXCL10, IL32, CCL20 (UBD left out, as job 63): z of log2 TMM CPM (prior 2) over
#   all included libraries (one cohort), mean over genes. PANX1 = log2 TMM CPM (prior 2) of ENSG00000110218.
# Units: staged included libraries with inferred_sex (stage codes 1-4; unstaged enter no test). OLS, HC3, t on
#   residual df; terms with one level among complete cases are dropped (job 63 rule).
# PC-geno (read first; gate): stage_num ~ PNPLA3_148M + batch + inferred_sex, estimate > 0 and p < 0.05.
# Family (Holm over 3, two-sided): T1 INFLAMED_HEP ~ PNPLA3_148M; T2 INFLAMED_HEP ~ HSD17B13_junction; T3 PANX1 ~
#   PANX1 3'UTR alt dosage (94181312 G>C / 94181824 T>G, max of the two calls, job 71b rule); each + batch + inferred_sex.
#   HSD17B13 junction genotype (job 63 rule): >= 30 junction reads; skip ratio <= 0.02 -> 0, 0.04-0.32 -> 1, >= 0.45 -> 2;
#   used if call rate >= 0.90 and HWE chi-square p >= 1e-4.
# Beside: + stage_num (mediation check); per-batch estimates; T3 allelic imbalance (job 71b P1 and S) and the job 72b
#   reference-bias null on the S0d null sites; PANX1 dosage -> stage_num (expected null).
# --dataset-term (774-unit reproduction of DIG1 job 63 only): z within dataset, + dataset in every model, HWE per
#   dataset (genotype step), and no staged-only restriction (job 63 used all units with a call).
# Unit: genetic individual (S2 Q3). No ancestry covariates. Cross-sectional.
# Usage: Rscript s3f_genotype_models.R --input <S2 dir> --stage <stage tsv> --geno <s3f_genotype_table dir>
#        --counts <s3_allele_counts dir> --targets <out_s0d_targets dir> --out <new dir> [--dataset-term]
suppressPackageStartupMessages({library(edgeR); library(data.table)})
a <- commandArgs(trailingOnly = TRUE)
opt <- function(k) { i <- match(k, a); if (is.na(i)) NA_character_ else a[i + 1] }
DSET <- "--dataset-term" %in% a
inp <- opt("--input"); stf <- opt("--stage"); gdir <- opt("--geno"); cdir <- opt("--counts"); tdir <- opt("--targets"); out <- opt("--out")
stopifnot(!anyNA(c(inp, stf, gdir, cdir, tdir, out)), !dir.exists(out)); dir.create(out, recursive = TRUE)
SS <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z"
stopifnot(system(paste("cd", file.path(SS, "frozen"), "&& sha256sum -c --quiet SHA256SUMS")) == 0)
U <- fread(file.path(inp, "units.tsv"))[toupper(as.character(include)) == "TRUE"]
U[, sample_alias := as.character(sample_alias)]
if (!DSET) U[, dataset := "one_cohort"]
st <- fread(cmd = paste("tr -d '\\r' <", shQuote(stf)), colClasses = "character")
U[, stage_num := suppressWarnings(as.numeric(st$nash_crn_kleiner_fibrosis_stage[match(sample_alias, st$sample_alias)]))]
d <- readRDS(file.path(inp, "dge_included.rds")); stopifnot(identical(colnames(d$counts), U$sample_alias))
deg <- fread(file.path(SS, "frozen/release/deg_results.csv")); sym <- setNames(deg$symbol, deg$gene)
rn <- rownames(d$counts); g_sym <- ifelse(rn %in% names(sym), sym[rn], rn)
INF <- c("FCAMR", "TRIM31", "CXCL10", "IL32", "CCL20")
SETS <- list(INFLAMED_HEP = INF, IMMUNE = c("PTPRC", "CD3E", "CD2", "CD68", "CD163"), FIBRO = c("COL1A1", "COL1A2", "COL3A1", "LUM"),
             DUCTULAR = c("KRT7", "KRT19", "EPCAM", "SOX9"), HEPMASS = c("ALB", "APOA1", "CYP2E1", "HP"))
d <- calcNormFactors(d); lc <- cpm(d, log = TRUE, prior.count = 2)
for (s in names(SETS)) {
  idx <- which(g_sym %in% SETS[[s]]); v <- rep(NA_real_, ncol(lc))
  for (co in unique(U$dataset)) { j <- which(U$dataset == co); v[j] <- colMeans(t(scale(t(lc[idx, j, drop = FALSE]))), na.rm = TRUE) }
  set(U, j = s, value = v)
}
U[, PANX1 := lc[grep("^ENSG00000110218", rn), ]]
writeLines(sapply(names(SETS), function(s) paste(s, paste(intersect(SETS[[s]], g_sym), collapse = ","))), file.path(out, "set_genes_found.txt"))
# genotypes
W <- fread(file.path(gdir, "genotype_calls_wide.tsv")); QC <- fread(file.path(gdir, "site_qc.tsv"))
pass <- QC$rsid[toupper(as.character(QC$pass_qc)) == "TRUE"]
U[, PNPLA3_148M := if ("rs738409" %in% pass) W$rs738409[match(sample_alias, W$sample_alias)] else NA_real_]
J <- fread(file.path(gdir, "hsd17b13_skip6.tsv"), colClasses = c(sample_alias = "character"))
J[, tot := skip6 + canon67]
J[, jg := fifelse(tot < 30, NA_real_, fifelse(skip_ratio <= 0.02, 0, fifelse(skip_ratio >= 0.04 & skip_ratio <= 0.32, 1, fifelse(skip_ratio >= 0.45, 2, NA_real_))))]
U[, HSD17B13_junction := J$jg[match(sample_alias, J$sample_alias)]]
hwe_chi <- function(g) { g <- g[!is.na(g)]; n <- c(sum(g == 0), sum(g == 1), sum(g == 2)); N <- sum(n); p <- (2 * n[1] + n[2]) / (2 * N)
  e <- N * c(p^2, 2 * p * (1 - p), (1 - p)^2); if (any(e == 0)) return(NA_real_); pchisq(sum((n - e)^2 / e), 1, lower.tail = FALSE) }
jqc <- U[, .(n = .N, called = sum(!is.na(HSD17B13_junction)), n0 = sum(HSD17B13_junction == 0, na.rm = TRUE),
             n1 = sum(HSD17B13_junction == 1, na.rm = TRUE), n2 = sum(HSD17B13_junction == 2, na.rm = TRUE), hwe_p = hwe_chi(HSD17B13_junction)), by = dataset]
jpass <- mean(!is.na(U$HSD17B13_junction)) >= 0.90 && all(jqc$hwe_p >= 1e-4, na.rm = TRUE)
fwrite(jqc, file.path(out, "hsd17b13_junction_qc.tsv"), sep = "\t")
if (!jpass) U[, HSD17B13_junction := NA_real_]
# PANX1 3'UTR dosage and het counts (job 71b)
A <- fread(file.path(cdir, "ase_counts.tsv.gz"), colClasses = c(sample_alias = "character"))
A <- A[sample_alias %in% U$sample_alias]
A[, dp := nref + nalt][, f := nalt / dp]
A[, call := fifelse(dp < 10, NA_real_, fifelse(f <= 0.1, 0, fifelse(f >= 0.9, 2, 1)))]
PX <- A[chrom == "chr11" & pos %in% c(94181312, 94181824)]
U2 <- PX[, .(dos = if (all(is.na(call))) NA_real_ else max(call, na.rm = TRUE)), by = sample_alias]
U[, PANX1_3UTR := U2$dos[match(sample_alias, U2$sample_alias)]]
px_sites <- A[chrom == "chr11" & pos >= 94128854 & pos <= 94181824, .(n = .N, call_rate = sum(!is.na(call)) / nrow(U),
               af = mean(call, na.rm = TRUE) / 2, n_het = sum(call == 1, na.rm = TRUE), med_dp = as.numeric(median(dp))), by = pos]
fwrite(px_sites, file.path(out, "panx1_site_summary.tsv"), sep = "\t")
fwrite(U, file.path(out, "sample_scores_genotypes.tsv"), sep = "\t")

fit <- function(f, dat) {
  v <- all.vars(f); cc <- dat[complete.cases(dat[, v, with = FALSE])]
  rhs <- attr(terms(f), "term.labels")
  rhs <- rhs[!vapply(rhs, function(t) t %in% c("dataset", "batch", "inferred_sex") && length(unique(cc[[t]])) < 2, logical(1))]
  if (nrow(cc) < 10 || length(unique(cc[[rhs[1]]])) < 2)
    return(data.frame(term = rhs[1], estimate = NA_real_, se = NA_real_, p = NA_real_, n = nrow(cc)))
  mod <- lm(reformulate(rhs, response = v[1]), data = cc); X <- model.matrix(mod); e <- residuals(mod); h <- hatvalues(mod)
  XtXi <- solve(crossprod(X)); V <- XtXi %*% crossprod(X * (e / (1 - h))) %*% XtXi
  b <- coef(mod); se <- sqrt(diag(V))
  data.frame(term = names(b), estimate = b, se = se, p = 2 * pt(-abs(b / se), df.residual(mod)), n = nobs(mod), row.names = NULL)
}
COV <- c(if (DSET) "dataset", "batch", "inferred_sex")
D <- if (DSET) U[!is.na(inferred_sex)] else U[!is.na(stage_num) & !is.na(inferred_sex)]
one <- function(y, x, dat = D, extra = NULL) { r <- fit(reformulate(c(x, extra, COV), y), dat); r[r$term == x, , drop = FALSE] }
pcg <- cbind(test = "PC-geno", outcome = "stage_num", one("stage_num", "PNPLA3_148M", U[!is.na(stage_num) & !is.na(inferred_sex)]))
pcg$pass <- isTRUE(pcg$estimate > 0 && pcg$p < 0.05)
fwrite(pcg, file.path(out, "pc_geno.tsv"), sep = "\t")
FAM <- rbind(cbind(test = "T1", outcome = "INFLAMED_HEP", one("INFLAMED_HEP", "PNPLA3_148M")),
             cbind(test = "T2", outcome = "INFLAMED_HEP", one("INFLAMED_HEP", "HSD17B13_junction")),
             cbind(test = "T3", outcome = "PANX1", one("PANX1", "PANX1_3UTR")))
FAM$p_holm <- p.adjust(FAM$p, "holm")
FAM$dig1_estimate <- c(0.176, -0.137, -0.124); FAM$dig1_se <- c(0.041, 0.045, 0.0375); FAM$dig1_n <- c(757, 764, 726)
FAM$read_as_result <- pcg$pass
fwrite(FAM, file.path(out, "f_ss6_family.tsv"), sep = "\t")
# beside
B <- list()
for (k in seq_len(nrow(FAM))) {
  x <- FAM$term[k]; y <- FAM$outcome[k]
  B[[length(B) + 1]] <- cbind(test = FAM$test[k], analysis = "plus_stage_num", outcome = y,
                               one(y, x, D[!is.na(stage_num)], extra = "stage_num"))
  for (b in sort(unique(D$batch))) {
    r <- fit(reformulate(c(x, "inferred_sex"), y), D[batch == b]); r <- r[r$term == x, , drop = FALSE]
    B[[length(B) + 1]] <- cbind(test = FAM$test[k], analysis = paste0("batch_", b), outcome = y, r)
  }
  for (s in setdiff(names(SETS), "INFLAMED_HEP")) if (y == "INFLAMED_HEP")
    B[[length(B) + 1]] <- cbind(test = FAM$test[k], analysis = "comparator_score", outcome = s, one(s, x))
}
B[[length(B) + 1]] <- cbind(test = "T3", analysis = "dosage_to_stage_expected_null", outcome = "stage_num",
                            one("stage_num", "PANX1_3UTR", U[!is.na(stage_num) & !is.na(inferred_sex)]))
B[[length(B) + 1]] <- cbind(test = "T2", analysis = "dosage_to_stage", outcome = "stage_num",
                            one("stage_num", "HSD17B13_junction", U[!is.na(stage_num) & !is.na(inferred_sex)]))
BES <- rbindlist(B, fill = TRUE); fwrite(BES, file.path(out, "beside_models.tsv"), sep = "\t")
# T3 allelic imbalance (job 71b P1 and S) and reference-bias null (job 72b)
H <- PX[call == 1][, .(nref = sum(nref), nalt = sum(nalt), n_snps = .N), by = sample_alias][nref + nalt >= 20]
H[, f := nalt / (nref + nalt)]
P1 <- data.table(test = "3UTR hets vs 0.5", n = nrow(H), median_f = if (nrow(H)) median(H$f) else NA_real_,
                 pooled_f = if (nrow(H)) sum(H$nalt) / sum(H$nref + H$nalt) else NA_real_,
                 p = if (nrow(H) > 1) wilcox.test(H$f, mu = 0.5)$p.value else NA_real_)
H <- merge(H, U, by = "sample_alias")
S <- rbindlist(lapply(c("INFLAMED_HEP", "stage_num"), function(v) {
  dd <- H[!is.na(get(v))]
  if (nrow(dd) < 10) return(data.table(term = v, n = nrow(dd), beta_logit = NA_real_, se = NA_real_, p = NA_real_))
  rhs <- c(v, if (length(unique(dd$batch)) > 1) "batch", if (DSET && length(unique(dd$dataset)) > 1) "dataset")
  co <- summary(glm(reformulate(rhs, "cbind(nalt, nref)"), family = quasibinomial, data = dd))$coefficients[v, ]
  data.table(term = v, n = nrow(dd), beta_logit = co[1], se = co[2], p = co[4]) }))
TG <- fread(file.path(tdir, "s1b_targets_annotated.tsv"))[role == "s3h_null"]
NS <- merge(A, TG[, .(chrom, pos, eref = ref, ealt = alt)], by = c("chrom", "pos"))[dp >= 20 & f > 0.1 & f < 0.9]
N <- NS[, .(n_het = .N, pooled_f = sum(nalt) / sum(dp), cls = paste0(eref[1], ">", ealt[1])), by = .(chrom, pos)][n_het >= 20]
CL <- c("G>C", "C>G", "T>G", "A>C"); px <- P1$pooled_f
NUL <- data.table(set = c("all null SNVs", "same base-change classes"), n_sites = c(nrow(N), sum(N$cls %in% CL)),
                  null_median_of_pooled = c(median(N$pooled_f), median(N[cls %in% CL]$pooled_f)), panx1_pooled = px,
                  p_emp = c((1 + sum(N$pooled_f <= px)) / (1 + nrow(N)), (1 + sum(N[cls %in% CL]$pooled_f <= px)) / (1 + sum(N$cls %in% CL))))
fwrite(P1, file.path(out, "panx1_P1_ase.tsv"), sep = "\t"); fwrite(S, file.path(out, "panx1_S_imbalance_vs_state.tsv"), sep = "\t")
fwrite(NUL, file.path(out, "panx1_vs_null.tsv"), sep = "\t"); fwrite(N, file.path(out, "null_site_summary.tsv"), sep = "\t")
writeLines(c(paste("mode:", if (DSET) "dataset-term reproduction" else "SteatoSITE one cohort"),
             paste("sites passing QC:", paste(pass, collapse = ",")), paste("HSD17B13 junction passes:", jpass),
             paste("units in analysis set:", nrow(D))), file.path(out, "run_notes.txt"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
options(width = 220)
print(readLines(file.path(out, "run_notes.txt"))); print(jqc); print(pcg, digits = 3); print(FAM, digits = 3)
print(BES, digits = 3); print(P1); print(S); print(NUL)
