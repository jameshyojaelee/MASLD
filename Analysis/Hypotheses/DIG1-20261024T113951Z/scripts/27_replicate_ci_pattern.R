# DIG1 job 27 (exploratory; replication of jobs 22, 22b, 26). Is the Complex I protein loss with
# Complex II spared seen in two other human liver proteomes, and is it specific to MASLD?
#   A  PXD027722: 20 NASH, 15 controls (bariatric and non-bariatric), 10 cirrhosis. Deposited
#      log2 Spectronaut matrix (no missing values in the deposit); the single "Gene name" column is used.
#   B  PXD065892: 15 NAFLD, 9 healthy, 19 autoimmune hepatitis (AIH). Deposited minprob-imputed
#      matrix (liver_sn_minprob_processed.tsv), columns GENE~ACCESSION.
# Processing as job 22: single-gene protein groups, quantile normalised, median per gene; scores are
# z within protein across the cohort's samples, mean over members; sets fixed as in job 22.
# Primary family (Holm over 2): CI - CII, disease minus control, permutation p (20,000 label
# shuffles, seed 20261024), A: NASH vs control; B: NAFLD vs healthy. Covariates are not used
# (age, BMI missing for some controls in A; none deposited in B).
# Secondary (BH over the table): every set and contrast, adding A: cirrhosis vs control and
# B: AIH vs healthy (is the pattern specific to metabolic disease?). Partner check within the
# disease samples of each cohort: Spearman of job 26 top partners and sets with CI - CII.
# Usage: Rscript 27_replicate_ci_pattern.R <new_out_dir>
suppressPackageStartupMessages({library(data.table); library(limma)})
args <- commandArgs(TRUE); stopifnot(length(args) == 1L, !dir.exists(args[1]))
out <- args[1]; dir.create(out, recursive = TRUE)
set.seed(20261024)
EXT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/data/external"
to_gene <- function(x, gsym) {
  ok <- !is.na(gsym) & gsym != "" & !grepl(";", gsym, fixed = TRUE)
  x <- x[ok, , drop = FALSE]; gsym <- gsym[ok]
  x[!is.finite(x)] <- NA_real_
  x <- normalizeBetweenArrays(x, method = "quantile")
  idx <- split(seq_len(nrow(x)), gsym)
  t(vapply(idx, function(i) if (length(i) == 1) x[i, ] else apply(x[i, , drop = FALSE], 2, median, na.rm = TRUE), numeric(ncol(x)))) }
# cohort A
pa <- fread(file.path(EXT, "PXD027722/download_20261017T012300Z/PatientSamples_proteome.csv"))
aa <- fread(file.path(EXT, "PXD027722/download_20261017T012300Z/PatientSamples_annotation_file.csv"), select = 1:6)
setnames(aa, c("patient", "sample", "group", "bmi", "age", "sex"))
sa <- intersect(aa$sample, names(pa)); aa <- aa[match(sa, sample)]
XA <- as.matrix(pa[, ..sa]); storage.mode(XA) <- "double"
GA <- to_gene(XA, pa[["Gene name"]]); colnames(GA) <- sa
aa[, grp := fcase(group == "NASH", "disease", group == "control", "control", group == "Cirrhosis", "cirrhosis")]
# cohort B
pb <- fread(file.path(EXT, "PXD065892/liver_20261017T003606Z/liver_sn_minprob_processed.tsv"))
XB <- t(as.matrix(pb[, -(1:3)])); storage.mode(XB) <- "double"
GB <- to_gene(XB, sub("~.*$", "", rownames(XB))); colnames(GB) <- pb$sample
ab <- data.table(sample = pb$sample, subject = pb$subject, group = pb$group)
ab[, grp := fcase(group == "NAFLD", "disease", group == "Healthy", "control", group == "AIH", "AIH")]
stopifnot(!anyDuplicated(ab$subject), !anyDuplicated(aa$patient))
sets_for <- function(g) {
  CI_all <- setdiff(g[grepl("^NDUF[ABSCV][0-9]", g)], "NDUFA4")
  S <- list(CI = CI_all, CII = intersect(c("SDHA", "SDHB", "SDHC", "SDHD"), g), CIII = g[grepl("^UQCR|^CYC1$", g)],
    CIV = union(g[grepl("^COX[4-8]", g)], intersect("NDUFA4", g)), CV = g[grepl("^ATP5", g)], mtDNA = g[grepl("^MT-", g)],
    mt_ribosome = g[grepl("^MRP[LS][0-9]", g)],
    mito_mass = intersect(c("VDAC1", "VDAC2", "VDAC3", "SLC25A4", "SLC25A5", "SLC25A6", "TOMM20", "TOMM40"), g),
    hep_mito = intersect(c("CPS1", "OTC", "HMGCS2", "GLS2"), g),
    FAO = intersect(c("CPT2", "ACADM", "ACADL", "ACADVL", "ACADS", "HADHA", "HADHB", "ECHS1", "HADH", "ETFA", "ETFB", "ETFDH"), g),
    proteasome_19S = intersect(c(paste0("PSMC", 1:6), paste0("PSMD", 1:14)), g),
    mt_proteases = intersect(c("LONP1", "CLPP", "CLPX", "YME1L1", "AFG3L2", "SPG7", "OMA1", "PARL", "HTRA2"), g),
    MDB_keratins = intersect(c("KRT8", "KRT18", "SQSTM1"), g), LD_coat = intersect(c("PLIN2", "PLIN4", "ABHD5"), g))
  S[lengths(S) > 0] }
PARTNERS <- c("OPA1", "LONP1", "PSMD2", "PSMC1", "ACTR1A", "DCTN1", "DYNLT3", "PPP1CA", "PPP2R1A", "TPD52", "PGRMC2", "ERP29", "MT1G")
score_tab <- function(G, ann, cohort) {
  S <- sets_for(rownames(G)); Z <- t(scale(t(G)))
  sc <- as.data.table(sapply(S, function(s) colMeans(Z[s, , drop = FALSE], na.rm = TRUE)))
  sc[, CI_minus_CII := CI - CII][, CI_minus_mass := CI - mito_mass]
  for (p in intersect(PARTNERS, rownames(G))) sc[, (p) := Z[p, ]]
  list(sets = data.table(cohort = cohort, set = names(S), n = lengths(S), genes = sapply(S, paste, collapse = ",")),
       scores = cbind(ann, sc)) }
A <- score_tab(GA, aa, "PXD027722"); B <- score_tab(GB, ab, "PXD065892")
fwrite(rbind(A$sets, B$sets), file.path(out, "sets_used.tsv"), sep = "\t")
fwrite(A$scores, file.path(out, "scores_PXD027722.tsv"), sep = "\t"); fwrite(B$scores, file.path(out, "scores_PXD065892.tsv"), sep = "\t")
perm <- function(y, lab, nperm = 20000) {
  ok <- !is.na(y); y <- y[ok]; lab <- lab[ok]; d <- mean(y[lab]) - mean(y[!lab])
  nd <- replicate(nperm, { l <- sample(lab); mean(y[l]) - mean(y[!l]) })
  c(diff = d, sd_pooled = sd(y), p = (1 + sum(abs(nd) >= abs(d) - 1e-12)) / (nperm + 1), n1 = sum(lab), n0 = sum(!lab)) }
contr <- list(list("PXD027722", A$scores, "disease", "control"), list("PXD065892", B$scores, "disease", "control"),
              list("PXD027722", A$scores, "cirrhosis", "control"), list("PXD065892", B$scores, "AIH", "control"))
prim <- rbindlist(lapply(contr[1:2], function(k) { d <- k[[2]][grp %in% c(k[[3]], k[[4]])]
  r <- perm(d$CI_minus_CII, d$grp == k[[3]]); data.table(cohort = k[[1]], contrast = paste(k[[3]], "vs", k[[4]]), t(r)) }))
prim[, p_holm := p.adjust(p, "holm")]
fwrite(prim, file.path(out, "primary_holm.tsv"), sep = "\t")
Y <- function(sc) setdiff(names(sc), c(names(aa), names(ab)))
sec <- rbindlist(lapply(contr, function(k) { d <- k[[2]][grp %in% c(k[[3]], k[[4]])]
  rbindlist(lapply(Y(d), function(y) { r <- perm(d[[y]], d$grp == k[[3]], 5000)
    data.table(cohort = k[[1]], contrast = paste(k[[3]], "vs", k[[4]]), score = y, t(r)) })) }))
sec[, bh := p.adjust(p, "BH")]
fwrite(sec, file.path(out, "all_contrasts.tsv"), sep = "\t")
within <- rbindlist(lapply(list(list("PXD027722", A$scores, c("disease", "cirrhosis")), list("PXD065892", B$scores, "disease"),
                                list("PXD065892", B$scores, "AIH")), function(k) {
  d <- k[[2]][grp %in% k[[3]]]
  rbindlist(lapply(setdiff(Y(d), c("CI_minus_CII", "CI", "CII")), function(y) { ct <- suppressWarnings(cor.test(d[[y]], d$CI_minus_CII, method = "spearman"))
    data.table(cohort = k[[1]], samples = paste(k[[3]], collapse = "+"), n = nrow(d), score = y, rho = unname(ct$estimate), p = ct$p.value) })) }))
within[, bh := p.adjust(p, "BH")]
fwrite(within, file.path(out, "within_disease_partners.tsv"), sep = "\t")
sub <- rbindlist(lapply(list(list("PXD027722", GA, aa), list("PXD065892", GB, ab)), function(k) {
  G <- k[[2]]; an <- k[[3]]; S <- sets_for(rownames(G)); gg <- c(S$CI, S$CII, S$CIV, S$CV, S$mt_ribosome)
  rbindlist(lapply(gg, function(p) { lab <- an$grp == "disease"; keep <- an$grp %in% c("disease", "control")
    data.table(cohort = k[[1]], gene = p, set = names(S)[sapply(S, function(s) p %in% s)][1],
               diff = mean(G[p, keep & lab]) - mean(G[p, keep & !lab])) })) }))
fwrite(sub, file.path(out, "per_subunit_disease_minus_control_log2.tsv"), sep = "\t")
fwrite(sub[, .(n = .N, mean_diff = mean(diff), frac_down = mean(diff < 0)), by = .(cohort, set)], file.path(out, "per_subunit_summary.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
print(prim, digits = 3); print(sec[score %in% c("CI", "CII", "CIII", "CIV", "CV", "mt_ribosome", "mito_mass", "hep_mito", "CI_minus_CII", "CI_minus_mass", "OPA1", "LONP1", "proteasome_19S", "LD_coat", "MDB_keratins")], digits = 2)
print(within[order(p)][1:30], digits = 2); print(fread(file.path(out, "per_subunit_summary.tsv")), digits = 2)
