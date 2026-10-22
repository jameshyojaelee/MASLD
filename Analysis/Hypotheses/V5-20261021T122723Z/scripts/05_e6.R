# V5 job 3c: E6 (117 programs, descriptive) and the E7 availability note. 04_outcome.R stopped at E6
# on a data.table scoping error (local 'set' masked by the column 'set'); this reruns only E6 with
# the same functions, permutation seed and design. Generated from 04_outcome.R lines 1-55, 69-105, 239-261.
suppressPackageStartupMessages({library(data.table); library(edgeR)})
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
RUN <- file.path(ROOT, "Analysis/Hypotheses/V5-20261021T122723Z")
OUT <- file.path(RUN, "outputs")
J3 <- file.path(OUT, "job3")
SEED <- 20260923
NPERM <- 19999
source(file.path(ROOT, "scripts/manuscript/program_observability_map/analysis_lib.R"))
MEMBERSHIP <- file.path(ROOT, "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot/program_membership_v2.tsv")
REGISTRY <- file.path(ROOT, "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot/program_registry_v2.tsv")
ANNOT <- file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BULK-F-FIVE/frozen_model_inputs/gencode_v49_gene_metadata.tsv.gz")
PRIMARY <- "hotspot_hepatocytes_48f39dd4d817a10e"
SECONDARY <- "hotspot_hepatocytes_f05c535ae5bbc0b9"
NC1 <- "hotspot_hepatocytes_bece6cd6d8fc9c67"
LINEAGE_SET <- c(hepatocytes = "hep", fibroblasts = "fibroblasts", cholangiocytes = "cholangiocytes",
                 macrophages = "macrophages", tcells = "tcells")

chk <- readLines(file.path(RUN, "CHECKED_BEFORE_OUTCOME.txt"))
if (!any(grepl("PROCEED: yes", chk, fixed = TRUE))) stop("STOP: CHECKED_BEFORE_OUTCOME.txt not PROCEED: yes")
if (sha256_file(file.path(RUN, "PRESPEC_V5.md")) != strsplit(readLines(file.path(RUN, "PRESPEC_V5.md.sha256")), " ")[[1]][1])
  stop("STOP: prespec hash differs")
if (sha256_file(MEMBERSHIP) != "feec3fc9ccaaffaa9abc1ac5d593cc06845ed8140460c8873a51bc8652d7336b") stop("membership hash")
if (sha256_file(REGISTRY) != "b134b5f8e46a716e857af22271aa54e1294432b9f9ab735559d7f603cdd29fe7") stop("registry hash")
if (file.exists(file.path(J3, "E6_117_programs_descriptive.tsv"))) stop("refusing to overwrite E6")

wr <- function(x, f) fwrite(as.data.table(x), file.path(J3, f), sep = "\t", na = "NA")
membership <- fread(MEMBERSHIP)
registry <- fread(REGISTRY)
annotation <- fread(ANNOT, select = c("gene_id", "gene_name"))
dq <- fread(file.path(OUT, "donor_qc.tsv"), colClasses = list(character = "id"))
gate <- setnames(fread(file.path(OUT, "gate_G0_G3.tsv"), header = TRUE), c("V1", "value"))
primary_model <- gate[V1 == "primary_model", value]
cells <- fread(file.path(J3, "cells_per_donor_set.tsv"))
det <- fread(file.path(J3, "detection_by_set.tsv.gz"))
cv <- dq[included == TRUE][order(gsm)]
cv[, log10_n_hep := log10(n_hep)]
n <- nrow(cv)

# ---- scoring (V3 route, detection-restricted rows) ---------------------------------
score_set <- function(set, donors, mem = membership, reg = registry) {
  pb <- fread(file.path(J3, sprintf("pseudobulk_%s.tsv.gz", set)))
  counts <- as.matrix(pb[, ..donors]); rownames(counts) <- pb$gene_id
  lc <- cpm(calcNormFactors(DGEList(counts)), log = TRUE, prior.count = 1)
  keep <- det$gene_id[det[[set]] >= 0.01]
  lc <- lc[rownames(lc) %in% keep, , drop = FALSE]
  S <- collapse_symbols(lc, rownames(lc), annotation)
  sm <- data.table(sample_id = donors, dataset = "GSE298719")
  list(S = S, sc = score_programs(S, sm, mem, reg, 0.80))
}
# ---- fitting helpers (copied from H90 05_fit_h90.R) ----------------------------------
hc3_fit <- function(Y, X, j) {
  Y <- as.matrix(Y)
  XtXi <- solve(crossprod(X))
  Cm <- XtXi %*% t(X)
  B <- Cm %*% Y
  E <- Y - X %*% B
  h <- rowSums((X %*% XtXi) * X)
  W <- E^2 / (1 - h)^2
  se_hc3 <- sqrt(as.vector((Cm[j, ]^2) %*% W))
  df <- nrow(X) - ncol(X)
  sigma <- sqrt(colSums(E^2) / df)
  list(beta = B[j, ], se = se_hc3, t = B[j, ] / se_hc3, df = df, sigma = sigma,
       se_ols = sigma * sqrt(XtXi[j, j]), h = h)
}
fl_null <- function(y, X, j, P) {
  X0 <- X[, -j, drop = FALSE]
  f0 <- X0 %*% solve(crossprod(X0), crossprod(X0, y))
  r0 <- as.vector(y - f0)
  hc3_fit(matrix(r0[t(P)], nrow = length(y)) + as.vector(f0), X, j)$t
}
perm_p <- function(t, tnull) (1 + sum(abs(tnull) >= abs(t))) / (length(tnull) + 1)
full_test <- function(y, X, j, P, label) {
  f <- hc3_fit(y, X, j)
  tn <- fl_null(y, X, j, P)
  q <- qt(0.975, f$df)
  list(row = data.table(test = label, n = length(y), beta = f$beta, hc3_se = f$se, hc3_t = f$t,
                        ci_low = f$beta - q * f$se, ci_high = f$beta + q * f$se, resid_sd = f$sigma,
                        perm_p = perm_p(f$t, tn), hc3_p_printed_not_used = 2 * pt(-abs(f$t), f$df),
                        df = f$df, max_leverage = max(f$h)),
       tnull = tn, fit = f)
}
draw_perm <- function(n) { set.seed(SEED); t(replicate(NPERM, sample.int(n))) }
design <- function(d, nhep = d$log10_n_hep, model = primary_model) {
  if (model == "score ~ disease") cbind(1, disease = d$disease) else
    cbind(1, disease = d$disease, sex_f = d$sex_f, log10_n_hep = nhep)
}
P <- draw_perm(n)
ce7 <- cells[set == "hep_h90umi"][match(cv$gsm, gsm), n_cells]
s7 <- score_set("hep_h90umi", cv$gsm[ce7 >= 30])$sc
wr(data.table(arm = "E7", donors_ge30 = sum(ce7 >= 30), primary_scored = PRIMARY %in% rownames(s7$primary) && all(is.finite(as.numeric(s7$primary[PRIMARY, ]))), coverage = s7$cohort_coverage[feature_id == PRIMARY]$cohort_weight_coverage), "E7_availability.tsv")
# ---- E6: 117 programs in their own lineage ---------------------------------------------
e6 <- rbindlist(lapply(names(LINEAGE_SET), function(ct) {
  pb_set <- LINEAGE_SET[[ct]]
  nc <- cells[set == pb_set][match(cv$gsm, gsm), n_cells]
  ok <- nc >= 30
  reg <- registry[cell_type == ct]
  base <- data.table(program_uid = reg$program_uid, cell_type = ct, module_name = reg$module_name, n = sum(ok),
                     testable = FALSE, beta = NA_real_, hc3_t = NA_real_, perm_p = NA_real_)
  if (sum(ok & cv$disease == 0) < 3 || sum(ok & cv$disease == 1) < 3) return(base)
  scm <- score_set(pb_set, cv$gsm[ok], membership[program_uid %in% reg$program_uid], reg)$sc$primary
  d <- cv[ok]; Xe <- cbind(1, disease = d$disease, sex_f = d$sex_f, log10_n = log10(nc[ok]))
  Pe <- if (all(ok)) P else draw_perm(sum(ok))
  for (i in seq_len(nrow(base))) {
    ye <- as.numeric(scm[base$program_uid[i], d$gsm])
    if (!all(is.finite(ye))) next
    r <- full_test(ye, Xe, 2, Pe, base$program_uid[i])$row
    base[i, `:=`(testable = TRUE, beta = r$beta, hc3_t = r$hc3_t, perm_p = r$perm_p)]
  }
  base
}))
e6[, bh_q_testable := NA_real_]
e6[testable == TRUE, bh_q_testable := p.adjust(perm_p, "BH")]
wr(e6, "E6_117_programs_descriptive.tsv")
writeLines(capture.output(sessionInfo()), file.path(J3, "sessionInfo_job3c.txt"))
print(e6[testable == TRUE][order(perm_p)][1:15])
