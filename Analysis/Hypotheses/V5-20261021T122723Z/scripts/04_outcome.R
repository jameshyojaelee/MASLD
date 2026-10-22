# V5 job 3 (b): scores, primary test, Q1, G5, PC4, NC1-NC4, fixed-sequence secondary,
# descriptive arms E1-E10. Follows PRESPEC_V5.md.
# Unit: donor. Score: V3 route (TMM log-CPM pseudobulk, rows restricted to genes detected in
# >= 1% of the set's cells, collapse_symbols, score_programs). Model: HC3 closed form,
# Freedman-Lane permutation p (19,999 draws, seed 20260923). Fitting helpers copied from
# H90 05_fit_h90.R.
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
if (file.exists(file.path(J3, "results_primary.tsv"))) stop("refusing to overwrite results_primary.tsv")

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
hep <- score_set("hep", cv$gsm)
cov_hep <- hep$sc$cohort_coverage
wr(cov_hep, "program_coverage_hep.tsv")
sc_hep <- hep$sc$primary
getscore <- function(scm, uid) as.numeric(scm[uid, cv$gsm])

# realised Ductular genes vs job 1 G3 list
d3 <- fread(file.path(OUT, "D3_confirmatory_gene_detection.tsv"))
g3_genes <- sort(d3[program_uid == PRIMARY & scored == TRUE, gene])
real <- sort(intersect(membership[program_uid == PRIMARY & !is.na(mapped_symbol), mapped_symbol], rownames(hep$S)))
wr(data.table(job1_g3_genes = paste(g3_genes, collapse = ";"), realised_genes = paste(real, collapse = ";"),
              identical = identical(g3_genes, real)), "realised_gene_set_check.tsv")

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

# ---- primary, G5 ------------------------------------------------------------------
y <- getscore(sc_hep, PRIMARY)
if (any(!is.finite(y))) stop("primary score missing (coverage below 0.80 after restriction?)")
X <- design(cv)
J <- 2
P <- draw_perm(n)
gz <- gzfile(file.path(J3, "permutation_matrix.tsv.gz"), "w"); write.table(P, gz, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE); close(gz)
prim <- full_test(y, X, J, P, sprintf("primary: Ductular injury, %s", primary_model))
X0 <- X[, -J, drop = FALSE]; f0 <- X0 %*% solve(crossprod(X0), crossprod(X0, y)); r0 <- as.vector(y - f0)
t_id <- hc3_fit(as.vector(f0) + r0[seq_len(n)], X, J)$t
g5 <- abs(t_id - prim$fit$t) < 1e-10
wr(data.table(check = "identity permutation reproduces observed HC3 t", diff = abs(t_id - prim$fit$t), pass = g5), "G5_check.tsv")
if (!g5) { writeLines("untestable: G5 failed", file.path(J3, "STOP_G5.txt")); stop("STOP G5") }
gz <- gzfile(file.path(J3, "null_t_primary.tsv.gz"), "w"); writeLines(format(prim$tnull, digits = 15), gz); close(gz)

# ---- PC4 --------------------------------------------------------------------------
q95 <- quantile(abs(prim$tnull), 0.95, names = FALSE)
pc4 <- rbindlist(lapply(c(design_MDE = 1.75, V3_C2_beta = 1.08, discovery_SH_mean = 0.83), function(delta) {
  hit <- vapply(1:1000, function(k) {
    dp <- cv$disease[P[k, ]]
    Xk <- X; Xk[, "disease"] <- dp
    # a relabelling can equal the male indicator (7 MASLD, 7 male), making the design singular; skipped
    if (qr(Xk)$rank < ncol(Xk)) return(NA)
    f <- hc3_fit(y + delta * dp * prim$fit$sigma, Xk, J)
    f$beta > 0 && abs(f$t) > q95
  }, logical(1))
  data.table(delta_sd = delta, detection = mean(hit, na.rm = TRUE), n_draws_used = sum(!is.na(hit)), null_abs_t_q95 = q95)
}))
wr(pc4, "PC4_detection.tsv")

# ---- NC1, NC2, NC3, NC4, Q1 ----------------------------------------------------------
y_nc1 <- getscore(sc_hep, NC1)
nc1 <- if (all(is.finite(y_nc1))) full_test(y_nc1, X, J, P, "NC1: Sinusoidal endothelial (module 2)")$row else data.table(test = "NC1", note = "not available (coverage)")
nc1[, call := if ("perm_p" %in% names(nc1)) perm_p < 0.05 else NA]
wr(nc1, "NC1_null_program.tsv")
nc2 <- mean(2 * pt(-abs(prim$tnull), prim$fit$df) < 0.05)
wr(data.table(fraction_null_hc3_p_below_0.05 = nc2, label = ifelse(nc2 > 0.08, "anticonservative", "not anticonservative")), "NC2_calibration.tsv")

S <- hep$S
Z <- zscore_rows(S)
variable <- rownames(Z)[rowSums(is.finite(Z)) == ncol(Z)]
mean_expr <- rowMeans(S[variable, , drop = FALSE])
dsym <- merge(det[, .(gene_id = sub("\\..*$", "", gene_id), d = hep)],
              unique(annotation[, .(gene_id = sub("\\..*$", "", gene_id), symbol = toupper(gene_name))]), by = "gene_id")
dsym <- dsym[, .(d = max(d)), by = symbol]
det_v <- dsym$d[match(variable, dsym$symbol)]
prog_genes <- unique(toupper(membership[program_uid %in% c(PRIMARY, SECONDARY, NC1) & !is.na(mapped_symbol), mapped_symbol]))
pool <- setdiff(variable, prog_genes)
dec <- function(v) cut(rank(v, ties.method = "first"), 10, labels = FALSE)
cell <- setNames(paste(dec(mean_expr), dec(det_v)), variable)
mm <- membership[program_uid == PRIMARY & mapped_symbol %in% variable, .(original_l1_weight = sum(original_l1_weight)), by = mapped_symbol]
set.seed(20261007)
rmem <- rbindlist(lapply(seq_len(1000), function(i) {
  used <- character()
  pick <- vapply(mm$mapped_symbol, function(g) {
    cand <- setdiff(pool[cell[pool] == cell[[g]]], used)
    if (!length(cand)) cand <- setdiff(pool, used)  # empty cell: any unused pool gene (counted below)
    p <- cand[sample.int(length(cand), 1L)]; used <<- c(used, p); p
  }, character(1))
  data.table(program_uid = sprintf("random_%04d", i), mapped_symbol = pick, original_l1_weight = mm$original_l1_weight)
}))
n_empty <- sum(vapply(mm$mapped_symbol, function(g) !any(cell[pool] == cell[[g]]), logical(1)))
rs <- score_programs(S, data.table(sample_id = cv$gsm, dataset = "GSE298719"), rmem,
                     data.table(program_uid = unique(rmem$program_uid)), 0.80)$primary
t_sets <- hc3_fit(t(rs[, cv$gsm]), X, J)$t
nc3_pct <- 100 * mean(t_sets < prim$fit$t)
wr(data.table(n_sets = length(t_sets), observed_t = prim$fit$t, percentile = nc3_pct, set_t_median = median(t_sets),
              set_t_q95 = quantile(t_sets, 0.95, names = FALSE), program_genes_with_empty_cell = n_empty), "NC3_matched_sets.tsv")

hz <- cv[disease == 0]
yh <- y[cv$disease == 0]
nc4 <- anova(lm(yh ~ factor(hz$plate)))
wr(data.table(F = nc4$`F value`[1], df1 = nc4$Df[1], df2 = nc4$Df[2], p_descriptive = nc4$`Pr(>F)`[1],
              plates = paste(unique(hz$plate), collapse = ";")), "NC4_plate_spread_healthy.tsv")
mixed <- cv$plate %in% c("E35250325_1", "E36240625_1")
Xq <- cbind(1, disease = cv$disease[mixed], plate = as.numeric(cv$plate[mixed] == "E36240625_1"))
q1f <- hc3_fit(y[mixed], Xq, 2)
q1_ok <- sign(q1f$beta) == sign(prim$fit$beta) && abs(q1f$beta) >= 0.5 * abs(prim$fit$beta)
wr(data.table(n = sum(mixed), beta_within_plate = q1f$beta, hc3_se = q1f$se, primary_beta = prim$fit$beta,
              robust_to_plate = q1_ok), "Q1_plate_qualifier.tsv")

# ---- primary result and secondary ---------------------------------------------------
pr <- prim$row
verdict <- if (pr$beta > 0 && pr$perm_p < 0.05) "pass" else if (pr$beta < 0 && pr$perm_p < 0.05) "opposite direction" else "fail"
y_sec <- getscore(sc_hep, SECONDARY)
g3s <- gate[V1 == "G3_secondary_pass", value] == "True" && all(is.finite(y_sec))
sec <- if (g3s) full_test(y_sec, X, J, P, "secondary: Stromal ECM (module 8)") else NULL
sec_status <- if (is.null(sec)) "untestable (G3 coverage)" else if (verdict != "pass") "not tested (primary did not pass); descriptive only" else
  if (sec$row$beta > 0 && sec$row$perm_p < 0.05) "pass" else "fail"
res <- copy(pr)[, verdict := verdict]
if (!is.null(sec)) res <- rbind(res, copy(sec$row)[, verdict := sec_status])
wr(data.table(gsm = cv$gsm, id = cv$id, group = cv$group, plate = cv$plate, ductular = y,
              stromal_ecm = y_sec, nc1 = y_nc1, n_hep = cv$n_hep), "donor_scores.tsv")
wr(res, "results_primary.tsv")

# ---- descriptive arms ---------------------------------------------------------------
E <- list()
sub_test <- function(keep, yv, label, Xfun = function(d) design(d)) {
  d <- cv[keep]; Pk <- if (all(keep)) P else draw_perm(sum(keep))
  full_test(yv[keep], Xfun(d), 2, Pk, label)$row
}
E$E1 <- sub_test(!(cv$id %in% c("12", "13")), y, "E1 without F4 donors 12, 13")
E$E2 <- full_test(y, cbind(1, disease = cv$disease), 2, P, "E2 no covariates")$row
e3 <- fread(file.path(J3, "E3_projection_scores.tsv"))
y3 <- e3[program_uid == PRIMARY][match(cv$gsm, gsm), projection]
E$E3 <- full_test(y3, X, J, P, "E3 H90 projection score")$row
lodo <- t(vapply(seq_len(n), function(i) { f <- hc3_fit(y[-i], X[-i, ], J); c(f$beta, f$t) }, numeric(2)))
wr(data.table(dropped = cv$gsm, beta = lodo[, 1], hc3_t = lodo[, 2]), "E4_lodo.tsv")
E$E4 <- data.table(test = "E4 leave-one-donor-out", n = n - 1, beta = min(lodo[, 1]), ci_low = NA, ci_high = max(lodo[, 1]),
                   note = "beta = min over deletions, ci_high column holds max")
cho <- cells[set == "cholangiocytes"]
ok5 <- cho[match(cv$gsm, gsm), n_cells] >= 20
if (sum(ok5 & cv$disease == 0) >= 3 && sum(ok5 & cv$disease == 1) >= 3) {
  s5 <- score_set("cholangiocytes", cv$gsm[ok5], membership[program_uid == PRIMARY], registry[program_uid == PRIMARY])$sc$primary
  y5 <- rep(NA_real_, n); y5[ok5] <- as.numeric(s5[PRIMARY, cv$gsm[ok5]])
  if (all(is.finite(y5[ok5]))) E$E5 <- sub_test(ok5, y5, "E5 cholangiocyte Ductular score",
    function(d) cbind(1, disease = d$disease, sex_f = d$sex_f, log10_n = log10(cho[match(d$gsm, gsm), n_cells])))
}
for (arm in list(list("E7", "hep_h90umi", "E7 H90 200-10,000 UMI hepatocytes"), list("E8", "hep_marker", "E8 marker-label hepatocytes"))) {
  ce <- cells[set == arm[[2]]][match(cv$gsm, gsm), n_cells]
  ok <- ce >= 30
  s <- score_set(arm[[2]], cv$gsm[ok])$sc$primary
  ya <- rep(NA_real_, n); ya[ok] <- as.numeric(s[PRIMARY, cv$gsm[ok]])
  if (all(is.finite(ya[ok])) && sum(ok & cv$disease == 0) >= 3 && sum(ok & cv$disease == 1) >= 3)
    E[[arm[[1]]]] <- sub_test(ok, ya, arm[[3]], function(d) design(d, log10(ce[match(d$gsm, cv$gsm)])))
}
E$E9 <- full_test(y, cbind(X, age = cv$age), J, P, "E9 primary + age")$row
ok10 <- cv$n_hep >= 50
if (sum(ok10 & cv$disease == 0) >= 3 && sum(ok10 & cv$disease == 1) >= 3)
  E$E10 <- sub_test(ok10, y, "E10 donors with >= 50 hepatocytes (score from 30-cell set)")
wr(rbindlist(E, fill = TRUE), "E_descriptive_arms.tsv")

# ---- E6: 117 programs in their own lineage ---------------------------------------------
e6 <- rbindlist(lapply(names(LINEAGE_SET), function(ct) {
  set <- LINEAGE_SET[[ct]]
  nc <- cells[set == (if (set == "hep") "hep" else set)][match(cv$gsm, gsm), n_cells]
  ok <- nc >= 30
  reg <- registry[cell_type == ct]
  base <- data.table(program_uid = reg$program_uid, cell_type = ct, module_name = reg$module_name, n = sum(ok),
                     testable = FALSE, beta = NA_real_, hc3_t = NA_real_, perm_p = NA_real_)
  if (sum(ok & cv$disease == 0) < 3 || sum(ok & cv$disease == 1) < 3) return(base)
  scm <- score_set(set, cv$gsm[ok], membership[program_uid %in% reg$program_uid], reg)$sc$primary
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

writeLines(capture.output(sessionInfo()), file.path(J3, "sessionInfo.txt"))
cat("verdict:", verdict, "| perm p", pr$perm_p, "| beta", pr$beta, "| Q1 robust", q1_ok, "| NC3 pct", nc3_pct, "\n")
