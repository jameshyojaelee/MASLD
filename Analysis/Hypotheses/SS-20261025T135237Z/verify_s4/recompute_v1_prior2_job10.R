# Independent recomputation of SS S3d / S3e H38 / S3g headline statistics.
# Written from PRESPEC_SS.md and addenda S3d, S3e, S3g only; scripts/ not read.
suppressPackageStartupMessages({
  library(data.table); library(edgeR); library(limma)
})
set.seed(20261027)
SS  <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z"
OUT <- file.path(SS, "verify_s4")
STAGE <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/data/external/SteatoSite/fastq/PRJEB58625_fibrosis_stage_20261025T133450Z.tsv"
GENEMETA <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BULK-F-FIVE/frozen_model_inputs/gencode_v49_gene_metadata.tsv.gz"

rows <- list()
add <- function(stat, mine, theirs, count = FALSE) {
  mine <- as.numeric(mine); theirs <- as.numeric(theirs)
  if (count) ok <- isTRUE(mine == theirs) else {
    rel <- abs(mine - theirs) / max(abs(theirs), 1e-300)
    ok <- isTRUE(rel < 1e-6)
  }
  rows[[length(rows) + 1]] <<- data.table(statistic = stat, mine = mine, theirs = theirs,
                                          rel_diff = if (count) abs(mine - theirs) else abs(mine - theirs) / max(abs(theirs), 1e-300),
                                          match = ok)
}

## ---------- shared unit table ----------
units <- fread(file.path(SS, "out_s2_qc/units.tsv"))
units <- units[include == TRUE]
st <- fread(STAGE, colClasses = "character")
st[] <- lapply(st, function(x) gsub("\r", "", x))
setnames(st, gsub("\r", "", names(st)))
st[, stage_num := suppressWarnings(as.integer(nash_crn_kleiner_fibrosis_stage))]
st <- st[stage_num %in% 1:4, .(sample_alias, stage_num)]
u <- merge(units, st, by = "sample_alias")
u <- u[!is.na(inferred_sex) & inferred_sex %in% c("F", "M")]
u[, batch := factor(batch)]; u[, inferred_sex := factor(inferred_sex)]
cat("included", nrow(units), "staged with sex", nrow(u), "\n")
add("staged units n", nrow(u), 493, TRUE)

ols_hc3 <- function(y, X) {
  # HC3 by hand: (X'X)^-1 X' diag(e^2/(1-h)^2) X (X'X)^-1
  colnames(X) <- paste0("X", colnames(X))
  q <- qr(X); stopifnot(q$rank == ncol(X))
  b <- qr.coef(q, y); e <- as.vector(y - X %*% b)
  XtXi <- chol2inv(qr.R(q))
  h <- rowSums((X %*% XtXi) * X)
  meat <- crossprod(X * (e / (1 - h)))
  V <- XtXi %*% meat %*% XtXi
  se <- sqrt(diag(V)); names(se) <- names(b) <- colnames(X)
  list(b = b, se = se, df = nrow(X) - ncol(X))
}

## ---------- 1. S3d ----------
ps <- fread(file.path(SS, "out_s3d-20261026T233706Z/program_scores.tsv.gz"))
d <- merge(u[, .(sample_alias, batch, inferred_sex, stage_num)], ps, by = "sample_alias")
X <- model.matrix(~ stage_num + batch + inferred_sex, d)
progs <- setdiff(names(ps), "sample_alias")
res <- rbindlist(lapply(progs, function(p) {
  y <- d[[p]]
  if (all(is.na(y))) return(data.table(feature_id = p, testable = FALSE))
  stopifnot(!anyNA(y))
  r <- ols_hc3(y, X)
  b <- r$b["Xstage_num"]; se <- r$se["Xstage_num"]; t <- b / se
  data.table(feature_id = p, testable = TRUE, n = nrow(X), df = r$df, beta = b, se = se, t = t,
             p2 = 2 * pt(-abs(t), r$df), p1 = pt(t, r$df, lower.tail = FALSE))
}), fill = TRUE)
res[testable == TRUE, q := p.adjust(p2, "BH")]
fwrite(res, file.path(OUT, "s3d_recomputed_programs.tsv"), sep = "\t")
f4 <- fread(file.path(SS, "out_s3d-20261026T233706Z/FSS4_confirmatory_programs.tsv"))[model == "primary"]
for (i in seq_len(nrow(f4))) {
  m <- res[feature_id == f4$feature_id[i]]
  lab <- f4$program[i]
  add(paste0("S3d F-SS4 ", lab, " beta"), m$beta, f4$beta[i])
  add(paste0("S3d F-SS4 ", lab, " SE HC3"), m$se, f4$se_hc3[i])
  add(paste0("S3d F-SS4 ", lab, " t"), m$t, f4$t[i])
  add(paste0("S3d F-SS4 ", lab, " n"), m$n, f4$n[i], TRUE)
  add(paste0("S3d F-SS4 ", lab, " df"), m$df, f4$df[i], TRUE)
}
holm <- p.adjust(res[feature_id %in% f4$feature_id]$p1, "holm")
cat("Holm one-sided F-SS4:", holm, "\n")
s5 <- fread(file.path(SS, "out_s3d-20261026T233706Z/FSS5_summary.tsv"))
add("S3d F-SS5 testable programs", sum(res$testable), s5$n_testable, TRUE)
add("S3d F-SS5 BH q<0.05", sum(res$q < 0.05, na.rm = TRUE), s5$n_bh_q05, TRUE)
add("S3d F-SS5 BH q<0.05 up", sum(res$q < 0.05 & res$beta > 0, na.rm = TRUE), s5$n_bh_q05_up, TRUE)
add("S3d F-SS5 BH q<0.05 down", sum(res$q < 0.05 & res$beta < 0, na.rm = TRUE), s5$n_bh_q05_down, TRUE)
allp <- fread(file.path(SS, "out_s3d-20261026T233706Z/FSS5_all_programs.tsv"))
cmp <- merge(res[testable == TRUE], allp[, .(feature_id, beta_o = beta, se_o = se_hc3, q_o = q_bh)], by = "feature_id")
add("S3d F-SS5 max rel diff beta (112 programs) < 1e-6 -> 0/1", as.integer(max(abs(cmp$beta - cmp$beta_o) / abs(cmp$beta_o)) < 1e-6), 1, TRUE)
cat("F-SS5 max rel diff beta", max(abs(cmp$beta - cmp$beta_o) / abs(cmp$beta_o)),
    "se", max(abs(cmp$se - cmp$se_o) / cmp$se_o), "q", max(abs(cmp$q - cmp$q_o) / cmp$q_o), "\n")

## ---------- expression (needed for H38 B score and S3g) ----------
dge <- readRDS(file.path(SS, "out_s2_qc/dge_included.rds"))
cnt <- dge$counts
gm <- fread(GENEMETA)
ens <- sub("\\..*$", "", rownames(cnt))
sym <- gm$gene_name[match(ens, gm$ensembl_base)]
cat("rows without symbol:", sum(is.na(sym)), "\n")
lib <- colSums(cnt)
nf_tmm <- calcNormFactors(DGEList(cnt, lib.size = lib), method = "TMM")$samples$norm.factors
logcpm <- function(prior, nf) cpm(DGEList(cnt, lib.size = lib, norm.factors = nf), log = TRUE, prior.count = prior)
L1 <- logcpm(1, nf_tmm); L2 <- logcpm(2, nf_tmm)
gene_row <- function(L, g) {
  idx <- which(sym == g)
  if (length(idx) == 0) return(NULL)
  if (length(idx) > 1) cat("note: symbol", g, "has", length(idx), "rows; using first\n")
  L[idx[1], ]
}
zs <- function(v) (v - mean(v)) / sd(v)
setscore <- function(L, genes) {
  M <- do.call(rbind, lapply(genes, function(g) { v <- gene_row(L, g); if (is.null(v)) NULL else zs(v) }))
  colMeans(M)
}

## ---------- 2. H38 ----------
h <- fread(file.path(SS, "out_s3e-20261026T233706Z/h38_person_table.tsv"))
Bmine <- setscore(L1, c("CD79A", "MS4A1"))
h[, Bmine := Bmine[match(sample_alias, colnames(cnt))]]
cat("H38 Bscore max abs diff vs table:", max(abs(h$Bmine - h$Bscore)), "\n")
add("H38 B-cell score max abs diff (0 = identical)", max(abs(h$Bmine - h$Bscore)), 0)
hs <- h[!is.na(stage_num) & n >= 10]
add("H38 persons", nrow(hs), 493, TRUE)
Xh <- model.matrix(~ stage_num + sexM + Bscore + factor(batch), hs)
bb_nll <- function(par, k, n, X) {
  eta <- X %*% par[1:ncol(X)]; mu <- plogis(eta)
  rho <- plogis(par[ncol(X) + 1]); th <- (1 - rho) / rho
  a <- mu * th; b <- (1 - mu) * th
  -sum(lbeta(k + a, n - k + b) - lbeta(a, b) + lchoose(n, k))
}
start <- c(coef(glm(cbind(hs$k, hs$n - hs$k) ~ Xh - 1, family = binomial)), qlogis(0.05))
o <- optim(start, bb_nll, k = hs$k, n = hs$n, X = Xh, method = "BFGS",
           control = list(maxit = 10000, reltol = 1e-16))
o <- optim(o$par, bb_nll, k = hs$k, n = hs$n, X = Xh, method = "Nelder-Mead",
           control = list(maxit = 20000, reltol = 1e-16))
o <- optim(o$par, bb_nll, k = hs$k, n = hs$n, X = Xh, method = "BFGS",
           control = list(maxit = 10000, reltol = 1e-16))
H <- optimHess(o$par, bb_nll, k = hs$k, n = hs$n, X = Xh)
se <- sqrt(diag(solve(H)))
bi <- which(colnames(Xh) == "stage_num")
hb <- o$par[bi]; hse <- se[bi]; hz <- hb / hse
cat(sprintf("H38 b %.10f se %.10f z %.8f rho %.8f conv %d nll %.8f\n", hb, hse, hz, plogis(o$par[ncol(Xh) + 1]), o$convergence, o$value))
add("H38 T2 b", hb, 0.1478355); add("H38 T2 se", hse, 0.02331166); add("H38 T2 z", hz, 6.3416977)
add("H38 T2 rho", plogis(o$par[ncol(Xh) + 1]), 0.06504335)

## ---------- 3. S3g ----------
ss <- u[match(colnames(cnt), sample_alias)]
keep <- !is.na(ss$sample_alias)
su <- ss[keep]
cov <- model.matrix(~ batch + inferred_sex, su)
pspear <- function(x, y, Z = NULL, rank_Z = FALSE) {
  rx <- rank(x); ry <- rank(y)
  C <- cov
  if (!is.null(Z)) { Z <- as.matrix(Z); if (rank_Z) Z <- apply(Z, 2, rank); C <- cbind(C, Z) }
  ex <- lm.fit(C, rx)$residuals; ey <- lm.fit(C, ry)$residuals
  r <- cor(ex, ey); df <- length(x) - 2 - (ncol(C) - 1)
  list(r = r, t = r * sqrt(df / (1 - r^2)), df = df)
}
pdir <- function(t, df, dir) if (dir == "<0") pt(t, df) else if (dir == ">0") pt(t, df, lower.tail = FALSE) else 2 * pt(-abs(t), df)
g1 <- function(g) gene_row(L1, g)[keep]
stage <- su$stage_num
tests <- fread(file.path(SS, "out_s3fgh-20261026T233706Z/s3g/f_ss7_tests.tsv"))
mine <- list()
put <- function(id, r, dir) mine[[id]] <<- data.table(id = id, est = r$r, stat = r$t, df = r$df, p = pdir(r$t, r$df, dir))
put("s45a", pspear(g1("TNFSF10"), stage), "<0")
put("s45b", pspear(g1("TNFRSF10B"), stage), ">0")
put("s45c", pspear(g1("TNFRSF10A"), stage), ">0")
put("s30a", pspear(g1("THRB"), stage), "<0")
put("s30b", pspear(g1("SLC16A10"), stage), "<0")
# job 10 sets, prior 2, z over all included units
bil <- setscore(L2, c("SOX9", "EPCAM", "KRT7", "KRT19", "CD24", "SPP1", "KRT23"))[keep]
hep <- setscore(L2, c("ALB", "APOA1", "APOC3", "SERPINA1", "TF", "HP", "AHSG", "FGB"))[keep]
p21 <- zs(gene_row(L2, "CDKN1A"))[keep]; gnmt <- zs(gene_row(L2, "GNMT"))[keep]
put("s11a", pspear(bil, stage), ">0")
put("s11b", pspear(p21, stage), ">0")
put("s11c", pspear(gnmt, stage), "<0")
put("s9a", pspear(gnmt, bil, hep), "<0")
put("s9b", pspear(gnmt, p21, hep), "<0")
r9 <- pspear(gnmt, bil, hep, rank_Z = TRUE); cat("s9a with ranked hep proxy r", r9$r, "\n")
# s38 job 58 sets, prior 2
infl <- setscore(L2, c("FCAMR", "TRIM31", "CXCL10", "IL32", "CCL20"))[keep]
imm  <- setscore(L2, c("PTPRC", "CD3E", "CD2", "CD163"))[keep]
duct <- setscore(L2, c("KRT7", "KRT19", "EPCAM", "SOX9"))[keep]
hepm <- setscore(L2, c("ALB", "APOA1", "CYP2E1", "HP"))[keep]
X38 <- cbind(model.matrix(~ stage_num + batch + inferred_sex, su), IMMUNE = imm, HEPMASS = hepm, DUCTULAR = duct)
r38 <- ols_hc3(infl, X38)
t38 <- r38$b["Xstage_num"] / r38$se["Xstage_num"]
mine[["s38"]] <- data.table(id = "s38", est = r38$b["Xstage_num"], stat = t38, df = r38$df, p = pt(t38, r38$df, lower.tail = FALSE))
# score columns vs saved unit_scores
us <- fread(file.path(SS, "out_s3fgh-20261026T233706Z/s3g/unit_scores.tsv"))
us <- us[match(su$sample_alias, sample_alias)]
for (cc in list(c("INFLAMED_HEP", "infl"), c("IMMUNE", "imm"), c("DUCTULAR", "duct"), c("HEPMASS", "hepm"),
                c("j10_biliary", "bil"), c("j10_hep_proxy", "hep"), c("j10_GNMT", "gnmt"), c("j10_p21", "p21")))
  cat("unit score", cc[1], "max abs diff", max(abs(us[[cc[1]]] - get(cc[2]))), "\n")

# loci: limma voomWithQualityWeights ~ batch + sex + stage_num, filterByExpr(group = bin)
cs <- cnt[, keep]
bin <- factor(ifelse(stage == 1, "F0F1", ifelse(stage == 2, "F2", "F3F4")))
Xl <- model.matrix(~ batch + inferred_sex + stage_num, su)
loci_t <- function(method) {
  y <- DGEList(cs, lib.size = colSums(cs))
  kf <- filterByExpr(y, group = bin); y <- y[kf, , keep.lib.sizes = TRUE]
  y <- calcNormFactors(y, method = method)
  v <- voomWithQualityWeights(y, Xl)
  fit <- eBayes(lmFit(v, Xl))
  data.table(gene = rownames(y), ensembl = sub("\\..*$", "", rownames(y)), t = fit$t[, "stage_num"])
}
lt_tmm <- loci_t("TMM"); lt_rle <- loci_t("RLE")
lo <- zcat <- fread(file.path(SS, "out_s3fgh-20261026T233706Z/s3g/limma_stage_num_t.tsv.gz"))
for (nm in c("tmm", "rle")) {
  lt <- get(paste0("lt_", nm)); mm <- merge(lt, lo, by = "ensembl")
  cat("limma", nm, "genes mine", nrow(lt), "theirs", nrow(lo), "overlap", nrow(mm), "max abs t diff", max(abs(mm$t.x - mm$t.y)), "\n")
}
cg <- fread(file.path(SS, "frozen/dig1/coloc_genes_lead_z.tsv"))
lc <- fread(file.path(SS, "frozen/dig1/loci.tsv"))
cg <- merge(unique(cg[, .(ensembl, gene, locus)]), lc[, .(locus, cls)], by = "locus")
cg <- cg[cls %in% c("release", "disease")]
loci_test <- function(lt) {
  g <- unique(merge(cg, lt[, .(ensembl, t)], by = "ensembl")[, .(ensembl, locus, cls, t)])
  w <- wilcox.test(abs(g[cls == "release"]$t), abs(g[cls == "disease"]$t), conf.int = TRUE, exact = FALSE)
  list(n = nrow(g), nrel = sum(g$cls == "release"), ndis = sum(g$cls == "disease"), W = unname(w$statistic), p = w$p.value, est = unname(w$estimate), g = g)
}
for (nm in c("tmm", "rle")) { z <- loci_test(get(paste0("lt_", nm))); cat("loci", nm, "n", z$n, z$nrel, z$ndis, "W", z$W, "p", z$p, "est", z$est, "\n") }
LT <- loci_test(lt_rle)
mine[["loci"]] <- data.table(id = "loci", est = LT$est, stat = LT$W, df = NA, p = LT$p)

M <- rbindlist(mine, fill = TRUE)
fwrite(M, file.path(OUT, "s3g_recomputed_tests.tsv"), sep = "\t")
for (i in seq_len(nrow(M))) {
  tt <- tests[id == M$id[i]]
  add(paste0("S3g ", M$id[i], " est"), M$est[i], tt$est)
  add(paste0("S3g ", M$id[i], " stat"), M$stat[i], tt$stat)
  add(paste0("S3g ", M$id[i], " p"), M$p[i], tt$p)
}
# BH over all 13 listed tests: their p's, with mine substituted where recomputed (s33 uses theirs)
pp <- tests$p; names(pp) <- tests$id
add("S3g family size", length(pp), 13, TRUE)
bh_theirs_p <- p.adjust(pp, "BH")
pp2 <- pp; pp2[M$id] <- M$p; bh_mine_p <- p.adjust(pp2, "BH")
for (i in seq_along(pp)) {
  add(paste0("S3g BH ", names(pp)[i], " (BH of their p)"), bh_theirs_p[i], tests$bh[i])
  add(paste0("S3g BH ", names(pp)[i], " (BH of my p)"), bh_mine_p[i], tests$bh[i])
}

R <- rbindlist(rows)
fwrite(R, file.path(OUT, "comparison.tsv"), sep = "\t")
print(R, nrows = 500)
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
