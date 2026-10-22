# Gate L02-02, login-node follow-up (light): null spread of candidate and alternative
# global statistics under a sex-restricted label permutation. Observed statistics are
# NOT computed. Proxy model: TMM logCPM, pooled-variance t per region.
suppressMessages({library(edgeR); library(data.table)})
set.seed(20261016)
root <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
sc   <- file.path(root, "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/gate-L02-02")
h3   <- file.path(root, "Analysis/MASLD_Model_Benchmark/executions/gse267145-matrix-audit-21064333/raw/GSE267119_H3K27ac.txt.gz")
join <- fread(file.path(root, "Analysis/MASLD_Model_Benchmark/executions/gse267145-authoritative-join-21064930/participant_join.tsv"))
gm   <- fread(file.path(sc, "gene_promoter_map.tsv"))
hdr <- scan(h3, what = "", nlines = 1, quiet = TRUE)
dt  <- fread(cmd = paste("zcat", h3), skip = 1, header = FALSE)
cnt <- as.matrix(dt[, -1, with = FALSE]); rownames(cnt) <- gsub('"', "", dt[[1]])
pid <- sub("_.*", "", hdr); colnames(cnt) <- pid
st  <- join$stage[match(pid, join$participant_id)]; sex <- join$sex[match(pid, join$participant_id)]
y <- calcNormFactors(DGEList(cnt), method = "TMM"); lc <- cpm(y, log = TRUE, prior.count = 2)
rm(dt, cnt); invisible(gc())

hit <- gm[canonical == TRUE & n_reg_gene_tss > 0]
# negative control: non-call genes matched 1:1 on bulk AveExpr and promoter mean count (nearest, no replacement)
pool <- gm[canonical == FALSE & n_reg_gene_tss > 0 & !(best_region %in% hit$best_region)]
zs <- function(a, ref) (a - mean(ref)) / sd(ref)
hx <- cbind(zs(hit$AveExpr, pool$AveExpr), zs(log2(hit$best_region_mean_count), log2(pool$best_region_mean_count)))
px <- cbind(zs(pool$AveExpr, pool$AveExpr), zs(log2(pool$best_region_mean_count), log2(pool$best_region_mean_count)))
used <- logical(nrow(pool)); mi <- integer(nrow(hit))
for (i in sample(nrow(hit))) {
  d <- (px[, 1] - hx[i, 1])^2 + (px[, 2] - hx[i, 2])^2; d[used] <- Inf
  mi[i] <- which.min(d); used[mi[i]] <- TRUE
}
nc <- pool[mi]
cat("matched negative-control genes", nrow(nc), "; |bulk logFC| median call", round(median(abs(hit$logFC)), 3),
    "non-call", round(median(abs(nc$logFC)), 3), "\n")

Xh <- lc[hit$best_region, ]; Xn <- lc[nc$best_region, ]
tstat <- function(M, g) {
  n1 <- sum(g); n0 <- sum(!g)
  m1 <- rowMeans(M[, g]); m0 <- rowMeans(M[, !g])
  v <- (rowSums((M[, g] - m1)^2) + rowSums((M[, !g] - m0)^2)) / (n1 + n0 - 2)
  list(d = m1 - m0, t = (m1 - m0) / sqrt(v * (1 / n1 + 1 / n0)), df = n1 + n0 - 2)
}
stat_set <- function(g, idx) {
  a <- tstat(Xh[, idx], g); b <- tstat(Xn[, idx], g)
  p <- 2 * pt(-abs(a$t), a$df); q <- p.adjust(p, "BH"); sg <- q < 0.05
  conc <- sum(sg & sign(a$d) == sign(hit$logFC)); disc <- sum(sg & sign(a$d) != sign(hit$logFC))
  c(rho_call = cor(hit$logFC, a$d, method = "spearman"),
    rho_nc = cor(nc$logFC, b$d, method = "spearman"),
    n_sig = sum(sg), conc_minus_disc = conc - disc,
    sign_agree_frac = mean(sign(a$d) == sign(hit$logFC)))
}
B <- 2000
fem <- sex == "F"
# (A) candidate's scheme: all 99, labels shuffled, 75/24 kept
dis <- st != "NOR"
nullA <- t(replicate(B, stat_set(sample(dis), rep(TRUE, 99))))
# (B) sex-restricted: men are all disease in the data, so permute only among the 85 women (61 disease / 24 NOR)
idxF <- which(fem); disF <- dis[fem]
nullB <- t(replicate(B, stat_set(sample(disF), idxF)))
summ <- function(M, lab) {
  for (s in colnames(M)) cat(sprintf("%s %-16s mean %8.3f sd %7.3f q2.5 %8.3f q97.5 %8.3f\n", lab, s,
      mean(M[, s]), sd(M[, s]), quantile(M[, s], .025), quantile(M[, s], .975)))
}
summ(nullA, "A_all99 "); summ(nullB, "B_women85")
cat(sprintf("B: MDE for rho_call ~ %.3f ; for rho_call - rho_nc ~ %.3f (sd of difference %.3f)\n",
    2.8 * sd(nullB[, "rho_call"]), 2.8 * sd(nullB[, "rho_call"] - nullB[, "rho_nc"]), sd(nullB[, "rho_call"] - nullB[, "rho_nc"])))
cat("B: fraction of draws with >=1 BH call:", mean(nullB[, "n_sig"] > 0), "; with >=20:", mean(nullB[, "n_sig"] >= 20),
    "; 95th pct n_sig", quantile(nullB[, "n_sig"], .95), "; 95th pct conc_minus_disc", quantile(nullB[, "conc_minus_disc"], .95), "\n")
# label-blind check of why the null is wide: correlation of bulk logFC with region loading on PC1 of promoter logCPM
pc <- prcomp(t(Xh), center = TRUE, scale. = FALSE)
cat("PC1 of 807 promoter regions explains", round(100 * pc$sdev[1]^2 / sum(pc$sdev^2), 1), "% of variance; Spearman(bulk logFC, PC1 loading) =",
    round(cor(hit$logFC, pc$rotation[hit$best_region, 1], method = "spearman"), 3), "\n")
fwrite(data.table(nullB), file.path(sc, "null_sexrestricted_proxy.tsv"), sep = "\t")
fwrite(nc[, .(gene, symbol, logFC, AveExpr, best_region, best_region_mean_count)], file.path(sc, "matched_negative_controls.tsv"), sep = "\t")
