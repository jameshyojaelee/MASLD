# REL1 (PRESPEC_REL1.md): hepatic ALT abundance by stage (A1 snRNA, A2 liver protein) and
# hepatocyte proteins in plasma (B1 baseline ballooning, B2 change after bariatric surgery).
# Part C is untestable by the prespecified donor-count rule and is not computed.
# Usage: Rscript rel1.R <new_out_dir>
suppressPackageStartupMessages({
  library(data.table); library(limma); library(edgeR); library(metafor); library(survival); library(ordinal)
})
args <- commandArgs(trailingOnly = TRUE)
OUT <- args[1]
if (is.na(OUT) || dir.exists(OUT)) stop("output directory missing or already exists")
dir.create(OUT)
SEED <- 20261110
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
HYP <- file.path(ROOT, "Analysis/Hypotheses")
V5 <- file.path(HYP, "V5-20261021T122723Z/outputs")
V3 <- file.path(HYP, "V3-20261008T015239Z/pseudobulk")
NF <- file.path(HYP, "NUCLEUS-FIX-20261107/unique_counts")
C5 <- file.path(HYP, "CONF1-20261109T221846Z/out_c5_extract_c1b")
PX <- file.path(ROOT, "data/PXD051911")
IN <- c(
  atlas_counts = file.path(NF, "pseudobulk/Hepatocytes_pseudobulk.csv"),
  atlas_donors = file.path(NF, "metadata/donor_metadata.tsv"),
  atlas_libraries = file.path(NF, "metadata/library_metadata.tsv"),
  atlas_sex = file.path(HYP, "RV1-20261031T155620Z/out-20261102T130905Z/donors.tsv"),
  v3_counts = file.path(V3, "pseudobulk_hepatocyte.tsv.gz"),
  v3_donors = file.path(V3, "donor_metadata.tsv"),
  v3_cells = file.path(V3, "nuclei_per_donor_lineage.tsv"),
  c5_counts = file.path(C5, "pseudobulk_hep.tsv.gz"),
  c5_donors = file.path(C5, "donors.tsv"),
  v5_counts = file.path(V5, "job3/pseudobulk_hep.tsv.gz"),
  v5_donors = file.path(V5, "donor_qc.tsv"),
  v5_cells = file.path(V5, "job3/cells_per_donor_set.tsv"),
  annotation = file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BULK-F-FIVE/frozen_model_inputs/gencode_v49_gene_metadata.tsv.gz"),
  px_meta = file.path(PX, "meta_data.txt"),
  px_liver = file.path(PX, "liver_protein_quant.txt"),
  px_plasma_raw = file.path(PX, "plasma_protein_quant.txt"),
  px_plasma_bc = file.path(PX, "plasma_protein_quant_filtered_and_batch_corrected.txt"))
stopifnot(all(file.exists(IN)))
wr <- function(x, f) fwrite(as.data.table(x), file.path(OUT, f), sep = "\t", na = "NA")
sha <- function(f) sub(" .*", "", system2("sha256sum", shQuote(f), stdout = TRUE))
wr(data.table(input = names(IN), path = IN, sha256 = vapply(IN, sha, "")), "inputs_sha256.tsv")

TARGETS <- c("GPT", "GOT1", "GOT2", "GGT1")
families <- list(); controls <- list()

# =============================== A1: snRNA per hepatocyte ===============================
MIN_HEP <- 30; DETECT <- 0.80
annot <- fread(IN["annotation"], select = c("ensembl_base", "gene_name"))
atlas_pb <- fread(IN["atlas_counts"])
load_atlas <- function(ds) {
  dm <- fread(IN["atlas_donors"])[dataset == ds & toupper(as.character(exclude)) == "FALSE"]
  lib <- fread(IN["atlas_libraries"])[donor %in% dm$donor]
  sex <- fread(IN["atlas_sex"], select = c("donor", "sex"))
  stopifnot(all(dm$donor %in% sex$donor), all(lib$sample %in% names(atlas_pb)))
  counts <- vapply(dm$donor, function(d) rowSums(as.matrix(atlas_pb[, lib[donor == d, sample], with = FALSE])),
                   numeric(nrow(atlas_pb)))
  meta <- data.table(donor = dm$donor, female = as.integer(sex$sex[match(dm$donor, sex$donor)] == "Female"),
                     stage_ord = as.numeric(dm$stage_ordinal),
                     n_hep = vapply(dm$donor, function(d) sum(lib[donor == d][["n_Hepatocytes"]]), numeric(1)))
  list(meta = meta, counts = counts, symbol = atlas_pb$gene)
}
load_v3 <- function() {
  m <- fread(IN["v3_donors"]); n3 <- fread(IN["v3_cells"]); pb <- fread(IN["v3_counts"])
  meta <- data.table(donor = m$donor, female = as.integer(m$sex == "female"), stage_ord = as.numeric(m$stage_ordinal),
                     n_hep = n3$hepatocyte[match(m$donor, n3$donor)])
  list(meta = meta, counts = as.matrix(pb[, meta$donor, with = FALSE]), symbol = pb$symbol)
}
load_c5 <- function() {
  pb <- fread(IN["c5_counts"]); m <- fread(IN["c5_donors"])
  meta <- data.table(donor = m$donor, female = as.integer(m$female), stage_ord = as.numeric(m$stage_ord), n_hep = m$n_hep)
  list(meta = meta, counts = as.matrix(pb[, meta$donor, with = FALSE]), symbol = pb$symbol)
}
load_v5 <- function() {
  dq <- fread(IN["v5_donors"])[included == TRUE]
  cells <- fread(IN["v5_cells"])[set == "hep"]
  pb <- fread(IN["v5_counts"])
  stopifnot(all(dq$disease == 0 | dq$fibrosis >= 2))
  meta <- data.table(donor = dq$gsm, female = as.integer(dq$sex_f),
                     stage_ord = ifelse(dq$disease == 0, 0, ifelse(dq$fibrosis >= 4, 3, 2)),
                     n_hep = cells$n_cells[match(dq$gsm, cells$gsm)])
  list(meta = meta, counts = as.matrix(pb[, meta$donor, with = FALSE]),
       symbol = annot$gene_name[match(sub("\\..*$", "", pb$gene_id), annot$ensembl_base)])
}
CO <- list(GSE202379 = load_atlas("GSE202379"), GSE244832 = load_atlas("GSE244832"),
           GSE256398 = load_v3(), GSE289173 = load_c5(), GSE298719 = load_v5())
rm(atlas_pb); invisible(gc())
POOLED <- c("GSE202379", "GSE244832", "GSE256398", "GSE289173")
A1_GENES <- c(TARGETS, "CYP2C19", "VCP", "GPT2")

voom_fit <- function(x, xvar, drop_stage1 = FALSE) {
  m <- x$meta[n_hep >= MIN_HEP]
  stopifnot(!anyNA(m[, .(female, stage_ord, n_hep)]))
  if (drop_stage1) m <- m[stage_ord != 1]
  m[, log10_n := log10(n_hep)]
  m[, xv := if (xvar == "stage") stage_ord else as.numeric(stage_ord >= 2)]
  uniq <- !is.na(x$symbol) & !(x$symbol %in% x$symbol[duplicated(x$symbol)])
  cnt <- x$counts[uniq, m$donor, drop = FALSE]; rownames(cnt) <- x$symbol[uniq]
  covs <- c("xv", if (sd(m$female) > 0) "female", "log10_n")
  design <- model.matrix(as.formula(paste("~", paste(covs, collapse = " + "))), data = m)
  dge <- DGEList(cnt)
  keep <- filterByExpr(dge, design)
  detect <- rowMeans(cnt > 0)
  dge <- calcNormFactors(dge[keep, , keep.lib.sizes = FALSE])
  v <- voom(dge, design)
  fit <- eBayes(lmFit(v, design))
  tt <- as.data.table(topTable(fit, coef = "xv", number = Inf, sort.by = "none", confint = TRUE), keep.rownames = "gene")
  tt[, `:=`(se = logFC / t, detect_frac = detect[gene], df = fit$df.total[1])]
  tt[, testable := detect_frac >= DETECT]
  list(tt = tt, n = nrow(m), stage_counts = paste(table(factor(m$stage_ord, 0:3)), collapse = "/"),
       covariates = paste(covs, collapse = "+"))
}
A1 <- lapply(names(CO), function(k) { cat("A1 voom", k, "\n"); voom_fit(CO[[k]], "stage") })
names(A1) <- names(CO)
per_cohort <- rbindlist(lapply(names(A1), function(k) {
  a <- A1[[k]]
  rbindlist(lapply(A1_GENES, function(g) {
    r <- a$tt[gene == g]
    if (nrow(r) == 0) return(data.table(cohort = k, gene = g, n_donors = a$n, stage_counts_0123 = a$stage_counts,
                                         in_filter = FALSE, testable = FALSE))
    r[, .(cohort = k, gene, n_donors = a$n, stage_counts_0123 = a$stage_counts, in_filter = TRUE, testable,
          detect_frac, AveExpr, slope = logFC, ci_low = CI.L, ci_high = CI.R, se, t, p = P.Value, covariates = a$covariates)]
  }), fill = TRUE)
}), fill = TRUE)
per_cohort[, assay := fifelse(cohort == "GSE298719", "Flex probe (secondary, not pooled)", "snRNA 3prime")]
wr(per_cohort, "A1_per_cohort.tsv")

pool <- function(d) {
  d <- d[testable == TRUE]
  if (nrow(d) < 3) return(data.table(k = nrow(d), pooled = FALSE))
  r <- rma.uni(yi = d$slope, sei = d$se, method = "DL", test = "knha")
  data.table(k = nrow(d), pooled = TRUE, n_donors = sum(d$n_donors), slope = as.numeric(r$b), ci_low = r$ci.lb,
             ci_high = r$ci.ub, p = r$pval, tau2 = r$tau2, Q = r$QE, Q_p = r$QEp, I2 = r$I2,
             cohorts = paste(d$cohort, collapse = ";"), n_negative = sum(d$slope < 0))
}
pooled <- rbindlist(lapply(A1_GENES, function(g) cbind(gene = g, analysis = "primary_4_snRNA",
                                                       pool(per_cohort[gene == g & cohort %in% POOLED]))), fill = TRUE)
no289 <- rbindlist(lapply(A1_GENES, function(g) cbind(gene = g, analysis = "without_GSE289173",
                                                      pool(per_cohort[gene == g & cohort %in% setdiff(POOLED, "GSE289173")]))), fill = TRUE)
wr(rbindlist(list(pooled, no289), fill = TRUE), "A1_pooled.tsv")
fa1 <- pooled[gene %in% TARGETS, .(family = "A1 snRNA pooled stage slope", item = gene, estimate = slope, ci_low, ci_high,
                                    p = fifelse(pooled, p, 1), tested = pooled)]
fa1[is.na(p), p := 1]; fa1[, holm := p.adjust(p, "holm")]
families$A1 <- fa1
controls$A1_pos <- pooled[gene == "CYP2C19", .(control = "A1 positive CYP2C19 pooled slope < 0, CI excludes 0",
  estimate = slope, ci_low, ci_high, pass = isTRUE(pooled) && slope < 0 && ci_high < 0)]
controls$A1_neg <- pooled[gene == "VCP", .(control = "A1 negative VCP |pooled slope| < 0.15",
  estimate = slope, ci_low, ci_high, pass = isTRUE(pooled) && abs(slope) < 0.15)]

# dilution check: percentile of slope among testable genes in the same AveExpr decile, per cohort
dil <- rbindlist(lapply(names(A1), function(k) {
  tt <- A1[[k]]$tt[testable == TRUE]
  tt[, decile := cut(AveExpr, quantile(AveExpr, 0:10 / 10), include.lowest = TRUE, labels = FALSE)]
  rbindlist(lapply(intersect(A1_GENES, tt$gene), function(g) { r <- tt[gene == g]; same <- tt[decile == r$decile]
    data.table(cohort = k, gene = g, decile = r$decile, slope = r$logFC, pct_genes_lower = mean(same$logFC < r$logFC),
               decile_median_slope = median(same$logFC)) }))
}))
wr(dil, "A1_dilution.tsv")
# disease (stage >= 2) vs healthy, stage 1 donors left out
A1d <- lapply(names(CO), function(k) voom_fit(CO[[k]], "disease", drop_stage1 = TRUE)); names(A1d) <- names(CO)
dis <- rbindlist(lapply(names(A1d), function(k) A1d[[k]]$tt[gene %in% A1_GENES,
  .(cohort = k, gene, n_donors = A1d[[k]]$n, testable, contrast = "stage>=2 vs healthy", logFC, ci_low = CI.L, ci_high = CI.R, se, p = P.Value)]))
dis_pool <- rbindlist(lapply(A1_GENES, function(g) { d <- dis[gene == g & cohort %in% POOLED]; setnames(d, "logFC", "slope")
  cbind(gene = g, analysis = "pooled_disease_vs_healthy", pool(d)) }), fill = TRUE)
wr(rbindlist(list(dis, dis_pool), fill = TRUE), "A1_secondary.tsv")

# =============================== A2: liver protein ======================================
meta <- fread(IN["px_meta"])
liv <- fread(IN["px_liver"])
lm_meta <- meta[sample_group == "initial_sample" & !is.na(liver_proteomics_filename) & liver_proteomics_filename != ""]
stopifnot(nrow(lm_meta) == 58, all(lm_meta$liver_proteomics_filename %in% names(liv)))
lm_meta[, `:=`(saf_ord = c(No_MASLD = 0, MASL = 1, MASH = 2)[saf_diagnosis],
               batch = ifelse(grepl("^2019", liver_proteomics_filename), "b2019", "b2020"),
               female = as.integer(gender == "Female"), fib = as.numeric(sub("F", "", kleiner_fibrosis_grade)),
               masld = as.integer(saf_diagnosis != "No_MASLD"))]
stopifnot(!anyNA(lm_meta[, .(saf_ord, alder, bmi, female, fib)]))
LX <- log2(as.matrix(liv[, lm_meta$liver_proteomics_filename, with = FALSE]))
LX[!is.finite(LX)] <- NA
rownames(LX) <- paste(liv$ProteinAccessions, liv$Genes, sep = "|")
LX <- LX[rowMeans(is.na(LX)) <= 0.5, ]
liv_fit <- function(xvar) {
  d <- model.matrix(as.formula(paste("~", xvar, "+ batch + alder + bmi + female")), data = lm_meta)
  f <- eBayes(lmFit(LX, d), robust = TRUE)
  tt <- as.data.table(topTable(f, coef = xvar, number = Inf, sort.by = "none", confint = TRUE), keep.rownames = "row")
  tt[, gene := sub(".*\\|", "", row)][, n_proteins := nrow(LX)][, predictor := xvar]
}
A2_GENES <- c(TARGETS, "MAT1A", "VCP", "GPT2")
A2 <- rbindlist(lapply(c("saf_ord", "fib", "masld"), function(xv) liv_fit(xv)[gene %in% A2_GENES,
  .(predictor, gene, protein_group = sub("\\|.*", "", row), n = 58L, n_proteins, log2_effect = logFC, ci_low = CI.L,
    ci_high = CI.R, p = P.Value, bh_all_proteins = adj.P.Val)]))
wr(A2, "A2_liver_protein.tsv")
stopifnot(all(TARGETS %in% A2[predictor == "saf_ord", gene]))
fa2 <- A2[predictor == "saf_ord" & gene %in% TARGETS, .(family = "A2 liver protein SAF slope", item = gene,
           estimate = log2_effect, ci_low, ci_high, p, tested = TRUE)]
fa2[, holm := p.adjust(p, "holm")]; families$A2 <- fa2
controls$A2_pos <- A2[predictor == "saf_ord" & gene == "MAT1A", .(control = "A2 positive MAT1A SAF slope < 0, CI excludes 0",
  estimate = log2_effect, ci_low, ci_high, pass = log2_effect < 0 & ci_high < 0)]
controls$A2_neg <- A2[predictor == "saf_ord" & gene == "VCP", .(control = "A2 negative VCP |SAF slope| < 0.15",
  estimate = log2_effect, ci_low, ci_high, pass = abs(log2_effect) < 0.15)]

# =============================== B: plasma ==============================================
raw <- fread(IN["px_plasma_raw"]); bc <- fread(IN["px_plasma_bc"])
v1 <- meta[sample_group == "initial_sample" & !is.na(plasma_proteomics_filename) & plasma_proteomics_filename != ""]
v3 <- meta[sample_group == "follow_up_sample"]
stopifnot(nrow(v1) == 136, nrow(v3) == 41, all(v3$patient_name %in% v1$patient_name),
          !anyNA(v1[, .(alder, bmi, gender, hepatocellular_ballooning_score, plasma_batch_effects)]))
v1[, `:=`(balloon = as.integer(hepatocellular_ballooning_score >= 1), female = as.integer(gender == "Female"),
          batch = factor(plasma_batch_effects))]
ACC <- c(ALDOB = "P05062", ADH1B = "P00325;P07327", ADH4 = "P08319", GOT1 = "P17174",
         SHBG = "P04278", ALB = "P02768", LPA = "P08519")
RAW <- log2(as.matrix(raw[, -(1:3)])); RAW[!is.finite(RAW)] <- NA; rownames(RAW) <- raw$ProteinAccessions
BC <- as.matrix(bc[, -1]); rownames(BC) <- bc$ProteinAccessions
stopifnot(all(ACC %in% rownames(RAW)))
miss_base <- sapply(ACC, function(a) mean(is.na(RAW[a, v1$plasma_proteomics_filename])))
use_bc <- sapply(ACC, function(a) a %in% rownames(BC)) & miss_base <= 0.20
rule <- data.table(protein = names(ACC), accession = ACC, raw_missing_baseline = miss_base,
                   in_bc = ACC %in% rownames(BC), model = ifelse(use_bc, "limma_bc", "tobit_raw"))
wr(rule, "B_matrix_rule.tsv")

# B1 limma on the batch-corrected matrix (all proteins), baseline samples
d1 <- model.matrix(~ balloon + batch + alder + female + bmi, data = v1)
fb <- eBayes(lmFit(BC[, v1$unique_identifier], d1), robust = TRUE)
tb <- as.data.table(topTable(fb, coef = "balloon", number = Inf, sort.by = "none", confint = TRUE), keep.rownames = "acc")
tf <- as.data.table(topTable(fb, coef = "female", number = Inf, sort.by = "none", confint = TRUE), keep.rownames = "acc")
tobit <- function(a, X) {
  y <- RAW[a, X$plasma_proteomics_filename]; det <- !is.na(y)
  okb <- names(which(tapply(det, X$batch, sum) > 0)); keep <- as.character(X$batch) %in% okb
  X <- droplevels(X[keep]); y <- y[keep]; det <- det[keep]
  fl <- tapply(y, X$batch, min, na.rm = TRUE)
  yc <- ifelse(det, y, fl[as.character(X$batch)])
  res <- data.table(n = length(y), n_detected = sum(det), n_dropped_batch = sum(!keep))
  if (sum(det) < 20) return(cbind(res, tested = FALSE, note = "fewer than 20 detected"))
  f <- tryCatch(survreg(Surv(yc, det, type = "left") ~ balloon + batch + alder + female + bmi, data = X, dist = "gaussian"),
                error = function(e) NULL)
  if (is.null(f)) return(cbind(res, tested = FALSE, note = "Tobit failed"))
  b <- coef(f)["balloon"]; s <- sqrt(vcov(f)["balloon", "balloon"])
  lg <- glm(det ~ balloon + batch + alder + female + bmi, data = X, family = binomial)
  cbind(res, tested = TRUE, note = "", estimate = b, ci_low = b - 1.96 * s, ci_high = b + 1.96 * s, p = 2 * pnorm(-abs(b / s)),
        sigma = f$scale, detect_logOR = coef(lg)["balloon"], detect_logOR_se = sqrt(vcov(lg)["balloon", "balloon"]),
        detect_p = summary(lg)$coefficients["balloon", 4])
}
B1 <- rbindlist(lapply(c("ALDOB", "ADH1B", "ADH4", "GOT1", "ALB"), function(g) {
  a <- ACC[g]
  if (use_bc[g]) {
    r <- tb[acc == a]
    data.table(protein = g, model = "limma_bc", n = sum(!is.na(BC[a, v1$unique_identifier])), n_detected = sum(!is.na(BC[a, v1$unique_identifier])),
               n_dropped_batch = 0L, tested = TRUE, note = "", estimate = r$logFC, ci_low = r$CI.L, ci_high = r$CI.R, p = r$P.Value)
  } else cbind(protein = g, model = "tobit_raw", tobit(a, v1))
}), fill = TRUE)
wr(B1, "B1_baseline.tsv")
fb1 <- B1[protein %in% c("ALDOB", "ADH1B", "ADH4", "GOT1"), .(family = "B1 baseline ballooning any vs none", item = protein,
          estimate, ci_low, ci_high, p = fifelse(tested, p, 1), tested)]
fb1[, holm := p.adjust(p, "holm")]; families$B1 <- fb1
# B1 secondary: ALDOB raw log2 OLS; ALDOB ordinal model of ballooning 0/1/2
v1[, aldob_raw := RAW[ACC["ALDOB"], plasma_proteomics_filename]][, aldob_bc := BC[ACC["ALDOB"], unique_identifier]]
ols <- lm(aldob_raw ~ balloon + batch + alder + female + bmi, data = v1)
cm <- clm(factor(hepatocellular_ballooning_score, ordered = TRUE) ~ aldob_bc + batch + alder + female + bmi, data = v1[!is.na(aldob_bc)])
sc <- summary(cm)$coefficients
wr(data.table(analysis = c("ALDOB raw log2 OLS, ballooning any vs none", "ALDOB bc, proportional-odds ballooning 0/1/2, log-odds per log2 unit"),
              n = c(nobs(ols), nobs(cm)), estimate = c(coef(ols)["balloon"], sc["aldob_bc", 1]),
              ci_low = c(confint(ols)["balloon", 1], sc["aldob_bc", 1] - 1.96 * sc["aldob_bc", 2]),
              ci_high = c(confint(ols)["balloon", 2], sc["aldob_bc", 1] + 1.96 * sc["aldob_bc", 2]),
              p = c(summary(ols)$coefficients["balloon", 4], sc["aldob_bc", 4])), "B1_secondary.tsv")
controls$B_sex <- tf[acc == ACC["SHBG"], .(control = "B1 SHBG female minus male > 0, CI excludes 0", estimate = logFC,
                                           ci_low = CI.L, ci_high = CI.R, pass = logFC > 0 & CI.L > 0)]
controls$B_neg <- B1[protein == "ALB", .(control = "B1 negative ALB |ballooning difference| < 0.15", estimate, ci_low, ci_high,
                                         pass = abs(estimate) < 0.15)]

# B2 pairs
P <- data.table(patient = v3$patient_name)
b1 <- v1[match(P$patient, patient_name)]
P[, `:=`(id1 = b1$unique_identifier, id3 = v3$unique_identifier, f1 = b1$plasma_proteomics_filename,
         f3 = v3$plasma_follow_up_proteomics_filename, batch1 = droplevels(b1$batch), age = b1$alder, female = b1$female,
         d_steatosis = v3$steatosis_score - b1$steatosis_score, d_nas = v3$nafld_activity_score - b1$nafld_activity_score,
         d_bmi = v3$bmi - b1$bmi)]
pair_vals <- function(a, matrix) if (matrix == "bc") cbind(BC[a, P$id1], BC[a, P$id3]) else cbind(RAW[a, P$f1], RAW[a, P$f3])
hc3_slope <- function(y, x, Z) {
  X <- cbind(x = x, Z); XtXi <- solve(crossprod(X)); C <- XtXi %*% t(X); b <- C %*% y; e <- y - X %*% b
  h <- rowSums((X %*% XtXi) * X); se <- sqrt(sum(C[1, ]^2 * e^2 / (1 - h)^2)); df <- nrow(X) - ncol(X)
  c(estimate = b[1], se = se, ci_low = b[1] - qt(0.975, df) * se, ci_high = b[1] + qt(0.975, df) * se,
    p_hc3 = 2 * pt(-abs(b[1] / se), df), df = df)
}
NPERM <- 20000L
B2_one <- function(g, hist, offset) {
  a <- ACC[g]; mtx <- if (use_bc[g]) "bc" else "raw"
  V <- pair_vals(a, mtx); ok <- !is.na(V[, 1]) & !is.na(V[, 2])
  base <- data.table(protein = g, histology = hist, matrix = mtx, n_pairs_both = sum(ok))
  if (sum(ok) < 20) return(cbind(base, tested = FALSE, note = "fewer than 20 pairs detected at both visits"))
  d <- P[ok]; y <- V[ok, 2] - V[ok, 1]; x <- d[[hist]]
  Z <- model.matrix(~ batch1 + age + female, data = droplevels(d))
  st <- hc3_slope(y, x, Z)
  MZ <- diag(nrow(Z)) - Z %*% solve(crossprod(Z), t(Z))
  slope_of <- function(xx) { r <- MZ %*% xx; sum(r * y) / sum(r^2) }
  ry <- MZ %*% rank(y)
  rho_of <- function(xx) { r <- MZ %*% rank(xx); sum(r * ry) / sqrt(sum(r^2) * sum(ry^2)) }
  set.seed(SEED + offset)
  strata <- split(seq_along(x), as.character(d$batch1))
  perm <- replicate(NPERM, { xx <- x; for (s in strata) xx[s] <- x[s][sample.int(length(s))]; c(slope_of(xx), rho_of(xx)) })
  b0 <- slope_of(x); r0 <- rho_of(x)
  cbind(base, tested = TRUE, note = "", estimate = st[["estimate"]], ci_low = st[["ci_low"]], ci_high = st[["ci_high"]],
        p = (1 + sum(abs(perm[1, ]) >= abs(b0) - 1e-12)) / (NPERM + 1), p_hc3 = st[["p_hc3"]],
        partial_spearman = r0, p_spearman = (1 + sum(abs(perm[2, ]) >= abs(r0) - 1e-12)) / (NPERM + 1),
        sd_delta_hist = sd(x), mean_delta_protein = mean(y))
}
B2 <- rbindlist(lapply(seq_along(c("d_steatosis", "d_nas")), function(j) {
  h <- c("d_steatosis", "d_nas")[j]
  rbindlist(lapply(seq_along(c("ALDOB", "ADH1B", "ADH4", "GOT1")), function(i)
    B2_one(c("ALDOB", "ADH1B", "ADH4", "GOT1")[i], h, 10 * j + i)), fill = TRUE)
}), fill = TRUE)
wr(B2, "B2_pairs.tsv")
for (h in c("d_steatosis", "d_nas")) {
  f <- B2[histology == h, .(family = paste("B2 change vs", h), item = protein, estimate, ci_low, ci_high,
                            p = fifelse(tested, p, 1), tested)]
  f[, holm := p.adjust(p, "holm")]; families[[h]] <- f
}
# B2 secondary: + change in BMI; raw log2 matrix
sec <- rbindlist(lapply(c("d_steatosis", "d_nas"), function(h) {
  a <- ACC["ALDOB"]; out <- list()
  for (mtx in c("bc", "raw")) {
    V <- pair_vals(a, mtx); ok <- !is.na(V[, 1]) & !is.na(V[, 2]); d <- P[ok]; y <- V[ok, 2] - V[ok, 1]
    for (adj in c("base", "plus_delta_BMI")) {
      Z <- model.matrix(as.formula(if (adj == "base") "~ batch1 + age + female" else "~ batch1 + age + female + d_bmi"), data = droplevels(d))
      out[[length(out) + 1]] <- data.table(protein = "ALDOB", histology = h, matrix = mtx, adjustment = adj, n_pairs = sum(ok),
                                           as.list(hc3_slope(y, d[[h]], Z)))
    }
  }
  rbindlist(out)
}))
wr(sec, "B2_secondary.tsv")
# pairing control: LPA baseline vs follow-up
L <- pair_vals(ACC["LPA"], "bc"); okL <- !is.na(L[, 1]) & !is.na(L[, 2])
rho <- cor(L[okL, 1], L[okL, 2], method = "spearman")
set.seed(SEED + 99)
null <- replicate(10000, cor(L[okL, 1], sample(L[okL, 2]), method = "spearman"))
controls$B_pair <- data.table(control = "B2 pairing LPA Spearman > 0.7 and > 99th pct of 10,000 re-pairings",
                              estimate = rho, ci_low = NA_real_, ci_high = quantile(null, 0.99), pass = rho > 0.7 & rho > quantile(null, 0.99))

# =============================== write families and controls ===========================
wr(rbindlist(families, fill = TRUE), "families.tsv")
wr(rbindlist(controls, fill = TRUE), "controls.tsv")
wr(data.table(part = "C", result = "could not be tested",
              reason = "GSE256398 groups Healthy 6, MASLD F0 3, MASH F2-F3 4, MASH cirrhosis 4; rule needs >= 5 donors in each group"),
   "C_status.tsv")
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
print(rbindlist(families, fill = TRUE)); print(rbindlist(controls, fill = TRUE))
cat("REL1 DONE\n")
