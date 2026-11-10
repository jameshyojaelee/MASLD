# DIG1 job 63 (exploratory; written before any genotype from jobs 62/62b was read). Do liver-disease alleles,
# called from the bulk RNA reads, shift the inflamed-hepatocyte bulk score?
# Genotypes: count-based calls of job 62b. A site is used only if it passes the 62b QC (call rate >= 0.90, exact
# HWE p >= 1e-4 in every cohort with >= 20 calls, agreement with bcftools call >= 0.95).
# Allele coded (count of): PNPLA3 rs738409 G (148M); PNPLA3 rs2294918 G (434E); TM6SF2 rs58542926 T (167K);
#   GCKR rs1260326 T (446L; GRCh38 reference, so 2 - alt); MTARC1 rs2642438 A (165T; reference, 2 - alt);
#   MBOAT7 rs641738 T (reference, 2 - alt); APOE e4 = rs429358 C; APOE e2 = rs7412 T; HSD17B13 rs72613567 AA;
#   SERPINA1 rs28929474 T (Z), carrier vs non-carrier (descriptive only).
# Score: INFLAMED_HEP exactly as job 58 (FCAMR, TRIM31, CXCL10, IL32, CCL20; UBD absent from the matrix):
#   within-cohort z of log2 TMM CPM (prior 2), mean over genes. Comparators IMMUNE, FIBRO, DUCTULAR, HEPMASS.
# Model G: score ~ allele count + cohort + sex, all 774 units with a call; OLS, HC3 errors, t on residual df.
# Primary (one test): PNPLA3 148M, two-sided alpha 0.05; expected positive (148M raises NASH activity).
# Secondary family (Holm): the other QC-passing variants with pooled allele frequency 0.05-0.95 (APOE e4, e2 separate).
# Call-validity control, read before the primary: NAS ~ 148M + cohort + sex (GSE130970, GSE135251, GSE162694),
#   steatosis grade (GSE130970) and fibrosis (four cohorts) must not be negative for 148M (literature: positive).
# Descriptive: model G + NAS + fibrosis (cohorts with both; these may be mediators, so this is not a causal
#   decomposition); disease units only; per-cohort estimates with I2; comparator scores; each INFLAMED_HEP gene;
#   SERPINA1 Z carriers; PNPLA3 148M x 434E.
# Limits fixed now: no genome-wide genotypes, so ancestry within a cohort is not modelled; biopsy cohorts are
#   selected on disease; one sample per participant (release units). Cross-sectional.
# Deviation added before any genotype or association was read (after job 65 skip-ratio histogram only):
#   HSD17B13_junction = rs72613567 genotype read from exon 6 skipping (Ma 2019, doi 10.1002/hep.30350: the
#   minor allele produces exon-6-skipped transcripts). Job 65 skip ratio = exon5->7 / (exon5->7 + exon6->7) STAR
#   unique junction reads; >= 30 reads; <= 0.02 -> 0, 0.04-0.32 -> 1, >= 0.45 -> 2, otherwise not called.
#   Used if call rate >= 0.90 and HWE chi-square p >= 1e-4 in every cohort; it then replaces the read-level
#   insertion call (HSD17B13_ins, descriptive only) as the HSD17B13 test in the secondary family.
#   Expected direction negative (the allele protects against ballooning and lobular inflammation: Pirola 2018,
#   doi 10.1194/jlr.P089953). Agreement with the insertion calls is reported.
# Usage: Rscript 63_genotype_state.R <job62b_out_dir> <new_out_dir>
suppressPackageStartupMessages(library(edgeR))
args <- commandArgs(trailingOnly = TRUE); stopifnot(length(args) == 2)
gdir <- args[1]; out <- args[2]; stopifnot(!dir.exists(out)); dir.create(out, recursive = TRUE)
rel <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls"
d <- readRDS(file.path(rel, "merged_dge.rds")); m <- readRDS(file.path(rel, "meta_matched.rds"))
m <- m[match(colnames(d$counts), as.character(m$sample_id)), ]
stopifnot(identical(colnames(d$counts), as.character(m$sample_id)), identical(as.character(d$samples$dataset), as.character(m$dataset)))
deg <- read.csv(file.path(rel, "deg_results.csv")); sym <- setNames(deg$symbol, deg$gene)
rn <- rownames(d$counts); g_sym <- ifelse(rn %in% names(sym), sym[rn], rn)
INF <- c("FCAMR", "TRIM31", "UBD", "CXCL10", "IL32", "CCL20")
SETS <- c(list(INFLAMED_HEP = INF, IMMUNE = c("PTPRC", "CD3E", "CD2", "CD68", "CD163"), FIBRO = c("COL1A1", "COL1A2", "COL3A1", "LUM"),
               DUCTULAR = c("KRT7", "KRT19", "EPCAM", "SOX9"), HEPMASS = c("ALB", "APOA1", "CYP2E1", "HP")),
          setNames(as.list(setdiff(INF, "UBD")), paste0("g_", setdiff(INF, "UBD"))))
d <- calcNormFactors(d); lc <- cpm(d, log = TRUE, prior.count = 2)
sc <- sapply(SETS, function(gs) {
  idx <- which(g_sym %in% gs); v <- rep(NA_real_, ncol(lc))
  for (co in unique(m$dataset)) { j <- which(m$dataset == co); v[j] <- colMeans(t(scale(t(lc[idx, j, drop = FALSE]))), na.rm = TRUE) }
  v
})
D <- cbind(m[, c("sample_id", "dataset", "group_binary", "sex_final", "fibrosis_stage", "nas_score", "steatosis_grade")], as.data.frame(sc))
W <- read.delim(file.path(gdir, "genotype_calls_wide.tsv"), check.names = FALSE)
QC <- read.delim(file.path(gdir, "site_qc.tsv"))
pass <- QC$rsid[QC$pass_qc %in% c(TRUE, "True")]
CODE <- data.frame(rsid = c("rs738409", "rs2294918", "rs58542926", "rs1260326", "rs2642438", "rs641738", "rs429358", "rs7412", "rs72613567", "rs28929474"),
                   label = c("PNPLA3_148M", "PNPLA3_434E", "TM6SF2_167K", "GCKR_446L", "MTARC1_165T", "MBOAT7_rs641738T", "APOE_e4", "APOE_e2",
                             "HSD17B13_ins", "SERPINA1_Z"),
                   flip = c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE))
for (i in seq_len(nrow(CODE))) {
  r <- CODE$rsid[i]; v <- if (r %in% names(W)) W[[r]][match(D$sample_id, W$sample_id)] else rep(NA_real_, nrow(D))
  if (!(r %in% pass)) v[] <- NA_real_
  D[[CODE$label[i]]] <- if (CODE$flip[i]) 2 - v else v
}
D$SERPINA1_Z <- as.numeric(D$SERPINA1_Z > 0)
J <- read.delim(file.path(dirname(gdir), "out_hsd17b13_junctions/hsd17b13_skip6.tsv"))
J$tot <- J$skip6 + J$canon67
jg <- ifelse(J$tot < 30, NA, ifelse(J$skip_ratio <= 0.02, 0, ifelse(J$skip_ratio >= 0.04 & J$skip_ratio <= 0.32, 1,
             ifelse(J$skip_ratio >= 0.45, 2, NA))))
D$HSD17B13_junction <- jg[match(D$sample_id, J$sample_id)]
hwe_chi <- function(g) { g <- g[!is.na(g)]; n <- c(sum(g == 0), sum(g == 1), sum(g == 2)); N <- sum(n); p <- (2 * n[1] + n[2]) / (2 * N)
  e <- N * c(p^2, 2 * p * (1 - p), (1 - p)^2); if (any(e == 0)) return(NA_real_); pchisq(sum((n - e)^2 / e), 1, lower.tail = FALSE) }
jqc <- do.call(rbind, lapply(split(D, D$dataset), function(x) data.frame(dataset = x$dataset[1], n = nrow(x), called = sum(!is.na(x$HSD17B13_junction)),
  n0 = sum(x$HSD17B13_junction == 0, na.rm = TRUE), n1 = sum(x$HSD17B13_junction == 1, na.rm = TRUE), n2 = sum(x$HSD17B13_junction == 2, na.rm = TRUE),
  hwe_p = hwe_chi(x$HSD17B13_junction))))
jpass <- mean(!is.na(D$HSD17B13_junction)) >= 0.90 && all(jqc$hwe_p >= 1e-4, na.rm = TRUE)
ins_raw <- if ("rs72613567" %in% names(W)) W[["rs72613567"]][match(D$sample_id, W$sample_id)] else rep(NA_real_, nrow(D))
agree <- table(junction = D$HSD17B13_junction, insertion_call_before_qc = ins_raw, useNA = "ifany")
write.table(jqc, file.path(out, "hsd17b13_junction_qc.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(paste("junction genotype passes:", jpass), capture.output(print(agree))), file.path(out, "hsd17b13_junction_vs_insertion.txt"))
if (!jpass) D$HSD17B13_junction <- NA_real_
write.table(D, file.path(out, "sample_scores_genotypes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
fit <- function(f, dat) {
  # drop dataset / sex terms that have a single level among complete cases
  v <- all.vars(f); cc <- dat[complete.cases(dat[, v, drop = FALSE]), , drop = FALSE]
  rhs <- attr(terms(f), "term.labels"); rhs <- rhs[!(rhs %in% c("dataset", "sex_final") & sapply(rhs, function(t) t %in% names(cc) && length(unique(cc[[t]])) < 2))]
  f <- reformulate(rhs, response = v[1])
  mod <- lm(f, data = cc); X <- model.matrix(mod); e <- residuals(mod); h <- hatvalues(mod)
  XtXi <- solve(crossprod(X)); V <- XtXi %*% crossprod(X * (e / (1 - h))) %*% XtXi
  b <- coef(mod); se <- sqrt(diag(V))
  data.frame(term = names(b), estimate = b, se = se, p = 2 * pt(-abs(b / se), df.residual(mod)), n = nobs(mod), row.names = NULL)
}
freq <- sapply(CODE$label, function(l) mean(D[[l]], na.rm = TRUE) / 2)
usable <- CODE$label[CODE$rsid %in% pass & CODE$label != "SERPINA1_Z" & freq >= 0.05 & freq <= 0.95]
desc_only <- character(0)
if (jpass) { if ("HSD17B13_ins" %in% usable) desc_only <- "HSD17B13_ins"; usable <- c(setdiff(usable, "HSD17B13_ins"), "HSD17B13_junction") }
rows <- list(); add <- function(r, ...) { a <- list(...); for (k in names(a)) r[[k]] <- a[[k]]; rows[[length(rows) + 1]] <<- r }
# call-validity control
if ("PNPLA3_148M" %in% usable) {
  for (y in c("nas_score", "steatosis_grade", "fibrosis_stage")) {
    dd <- D[!is.na(D[[y]]) & !is.na(D$PNPLA3_148M), ]
    f <- if (length(unique(dd$dataset)) > 1) paste(y, "~ PNPLA3_148M + dataset + sex_final") else paste(y, "~ PNPLA3_148M + sex_final")
    r <- fit(as.formula(f), dd); add(r[r$term == "PNPLA3_148M", ], analysis = "control_histology", outcome = y, variant = "PNPLA3_148M")
  }
}
# model G, all scores, all usable variants
for (v in c(usable, desc_only, if ("rs28929474" %in% pass) "SERPINA1_Z")) for (s in names(SETS)) {
  dd <- D[!is.na(D[[v]]), ]
  r <- fit(as.formula(paste(s, "~", v, "+ dataset + sex_final")), dd); add(r[r$term == v, ], analysis = "G", outcome = s, variant = v)
  if (s == "INFLAMED_HEP") {
    d2 <- dd[!is.na(dd$nas_score) & !is.na(dd$fibrosis_stage), ]
    r <- fit(as.formula(paste(s, "~", v, "+ dataset + sex_final + nas_score + fibrosis_stage")), d2)
    add(r[r$term == v, ], analysis = "G_plus_histology", outcome = s, variant = v)
    r <- fit(as.formula(paste(s, "~", v, "+ dataset + sex_final")), d2); add(r[r$term == v, ], analysis = "G_histology_subset", outcome = s, variant = v)
    d3 <- dd[dd$group_binary == "Disease", ]
    r <- fit(as.formula(paste(s, "~", v, "+ dataset + sex_final")), d3); add(r[r$term == v, ], analysis = "G_disease_only", outcome = s, variant = v)
    for (co in unique(dd$dataset)) {
      d4 <- dd[dd$dataset == co, ]
      if (length(unique(d4[[v]])) < 2) next
      r <- fit(as.formula(paste(s, "~", v, "+ sex_final")), d4); add(r[r$term == v, ], analysis = "per_cohort", outcome = s, variant = v, cohort = co)
    }
  }
}
R <- do.call(rbind, lapply(rows, function(r) { for (k in c("cohort")) if (is.null(r[[k]])) r[[k]] <- NA; r }))
P <- subset(R, analysis == "G" & outcome == "INFLAMED_HEP" & variant == "PNPLA3_148M")
S <- subset(R, analysis == "G" & outcome == "INFLAMED_HEP" & variant %in% setdiff(usable, "PNPLA3_148M"))
S$p_holm <- p.adjust(S$p, "holm")
het <- do.call(rbind, lapply(unique(R$variant[R$analysis == "per_cohort"]), function(v) {
  x <- subset(R, analysis == "per_cohort" & variant == v & is.finite(se) & se > 0); w <- 1 / x$se^2; est <- sum(w * x$estimate) / sum(w)
  Q <- sum(w * (x$estimate - est)^2); data.frame(variant = v, k = nrow(x), fe = est, fe_se = sqrt(1 / sum(w)), Q = Q,
                                                 I2 = if (Q > 0) max(0, (Q - (nrow(x) - 1)) / Q) else 0)
}))
write.table(R, file.path(out, "genotype_models.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(P, file.path(out, "primary_pnpla3.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(S, file.path(out, "secondary_holm.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(het, file.path(out, "per_cohort_heterogeneity.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(paste("usable:", paste(usable, collapse = ",")), paste("descriptive only:", paste(desc_only, collapse = ",")), paste("freq:", paste(names(freq), round(freq, 3), collapse = "; "))),
           file.path(out, "variants_used.txt"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
options(width = 220)
print(readLines(file.path(out, "variants_used.txt"))); print(jqc); print(readLines(file.path(out, "hsd17b13_junction_vs_insertion.txt")))
print(subset(R, analysis == "control_histology"), digits = 3)
print(P, digits = 3); print(S, digits = 3)
print(subset(R, analysis != "per_cohort" & outcome == "INFLAMED_HEP"), digits = 3)
print(subset(R, analysis == "G" & outcome != "INFLAMED_HEP"), digits = 3)
print(subset(R, analysis == "per_cohort"), digits = 3); print(het, digits = 3)
